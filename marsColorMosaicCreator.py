
import os
import sys
import re
import subprocess
import numpy
import copy
import multiprocessing
import threading
import logging
import datetime

import IrgGeoFunctions
#import copyGeoTiffInfo
import mosaicTileManager # TODO: Normalize caps!
import MosaicUtilities
import hrscImageManager
import hrscFileCacher

import stackImagePyramid

"""

Existing tools:
- RegisterHrsc.cpp
    - Input  = Basemap, HRSC
    - Output = spatialTransform
- writeHrscColorPairs.cpp
    - Input  = Basemap, HRSC, spatialTransform
    - Output = File containing pixel color pairs
- transformHrscImageColor.cpp
    - Input  = HRSC, colorTransform
    - Output = Color transformed HRSC image


Total Mars map width meters = 21338954.2
height = 10669477.2

Input basemap is 128x256 tiles, each tile 45x45, mpp = ~1852
- Total 32768 tiles.
- Between +/-60 degrees, 86x256 tiles = 22016 tiles
With 32x  increase, each tile is 1440x1440,   ~6MB,  190GB, mpp = ~58
With 64x  increase, each tile is 2880x2880,  ~24MB,  760GB, mpp = ~29
With 128x increase, each tile is 5760x5760,  ~95MB,    3TB, mpp = ~14.5 <-- Probably fine
With 160x increase, each tile is 7200x7200, ~150MB,  4.6TB, mpp = ~11.6 <--- PLENTY of resolution!

Currently using 64x !!
TODO: Go down to 128!


If there are about 3600 HRSC images (more if we fetch updates from the last few months)
  at a batch size of 20, that is 180 batches!  To finish in two months (60 days) this 
  means 3 completed batches per day.

Batch procedure:
- Keep a count of how many images we have downloaded
- Keep processing until we have downloaded COUNT images
  - TODO: Run through the DB and add all images missing URL's to the bad image list.
- For each HRSC image, update all the tiles it overlaps.
- Before each tile is touched for the first time, make a backup of it!
- After a batch finishes, make the kml pyramid and send out an email.
- MANUAL INTERVENTION
  - Look at the kml pyramid and make sure things are ok.
  - Is there a way to highlight the modified tiles? 
- When a batch is confirmed, clean up the data and remove the backed up tiles.
- Start the next batch!

"""

#----------------------------------------------------------------------------
# Constants

# TODO: If utilization is poor we can improve the parallel processing system!
NUM_DOWNLOAD_THREADS = 5 # There are five files we download per data set
NUM_PROCESS_THREADS  = 8


IMAGE_BATCH_SIZE = 2 # This should be set equal to the HRSC cache size

# TODO: Need to manage the processed HRSC folders, not just the download folders!
#       - The batch management can take care of this



# Set up the log path here.
# - Log tiles are timestamped as is each line in the log file
LOG_FORMAT_STR = '%(asctime)s %(name)s %(message)s'
currentTime = datetime.datetime.now()
logPath = ('/byss/smcmich1/data/hrscMosaicLogs/hrscMosaicLog_%s.txt' % currentTime.isoformat() )
logging.basicConfig(filename=logPath,
                    format=LOG_FORMAT_STR,
                    level=logging.DEBUG)

BAD_HRSC_FILE_PATH = '/byss/smcmich1/repo/MassUpload/badHrscSets.csv'

# Currently used to control the area we operate over
#HRSC_FETCH_ROI = None # Fetch ALL hrsc images
#HRSC_FETCH_ROI = MosaicUtilities.Rectangle(-180.0, 180.0, -60.0, 60.0) # No Poles
#HRSC_FETCH_ROI = MosaicUtilities.Rectangle(-116.0, -110.0, -2.0, 3.5) # Restrict to a mountain region
HRSC_FETCH_ROI = MosaicUtilities.Rectangle(133.0, 142.0, 46, 50.0) # Viking 2 lander region

#-----------------------------------------------------------------------------------------
# Functions



def cacheManagerThreadFunction(databasePath, outputFolder, inputQueue, outputQueue):
    '''Thread to allow downloading of HRSC data in parallel with image processing.
       The input queue recieves three types of commands:
           "STOP" --> Finish current tasks, then exit.
           "KILL" --> Immediately kill all the threads and exit.
           "FETCH data_set_name" --> Fetch the specified data set
           "FINISHED data_set_name" --> Signals that this data set can be safely deleted
'''

    logger = logging.getLogger('DownloadThread')

    # Initialize a process pool to be managed by this thread
    downloadPool = None
    if NUM_DOWNLOAD_THREADS > 1:
        downloadPool = multiprocessing.Pool(processes=NUM_DOWNLOAD_THREADS)

    # Set up the HRSC file manager object
    logger.info('Initializing HRSC file caching object')
    hrscFileFetcher = hrscFileCacher.HrscFileCacher(databasePath, outputFolder, BAD_HRSC_FILE_PATH, downloadPool)

    while True:

        # Fetch the next requested data set from the input queue
        request = inputQueue.get() 

        # Handle stop request
        if request == 'STOP':
            logger.info('Download thread manager received stop request, stopping download threads...')
            # Gracefully wait for current work to finish
            if downloadPool:
                downloadPool.close()
                downloadPool.join()  
            break
            
        if request == 'KILL':
            logger.info('Download thread manager received kill request, killing download threads...')
            # Immediately stop all work
            if downloadPool:
                downloadPool.terminate()
                downloadPool.join()  
            break
      
        if 'FETCH' in request:
            dataSet = request[len('FETCH'):].strip()
            logger.info('Got request to fetch data set ' + dataSet)
            # Download this HRSC image using the thread pool
            # - TODO: Allow download overlap of multiple data sets at once!
            try:
                hrscInfoDict = hrscFileFetcher.fetchHrscDataSet(dataSet)
            except Exception, e:
                # When we fail to fetch a data set, send out a failure message and keep going.
                logger.error('Caught exception fetching data set ' + dataSet + '\n' + str(e))
                hrscInfoDict = {'setName': dataSet, 'error': True}
                outputQueue.put(hrscInfoDict)
                continue    
            logger.info('Finished fetching data set ' + dataSet)
            # Put the output information on the output queue
            outputQueue.put(hrscInfoDict)

   
        
        #--> Need to make sure we never delete an image until we are finished using it

        # TODO: Implement 'finished' message?

    # We only get here when we break out of the main loop
    outputQueue.put('STOPPED')
    logger.info('Download manager thread stopped.')


# For debugging only, the hrscFileCacher class does the actual call from the database.
def getHrscImageList():
    '''For just returns a fixed list of HRSC images for testing'''   
    return ['h0022_0000',
            'h0506_0000',
            'h2411_0000']#,
            #'h6419_0000']



#def getBasemapTileSet():
#    '''Returns a limited basemap ROI for debugging'''
#    return MosaicUtilities.Rectangle(100, 102, 64, 66)




def getCoveredOutputTiles(basemapInstance, hrscInstance):
    '''Return a bounding box containing all the output tiles covered by the HRSC image'''
    
    hrscBoundingBoxDegrees = hrscInstance.getBoundingBoxDegrees()
    
    # DEBUG!  Restrict to a selected area.
    hrscBoundingBoxDegrees = HRSC_FETCH_ROI.getIntersection(hrscBoundingBoxDegrees)
    
    intersectRect = basemapInstance.getIntersectingTiles(hrscBoundingBoxDegrees)
    return intersectRect
    
    #return intersectRect.getIntersection(debugRect)
    #return MosaicUtilities.Rectangle(196, 197, 92, 93) # DEBUG


def getHrscTileUpdateDict(basemapInstance, tileIndex, hrscInstance):
    '''Gets the dictionary of HRSC tiles that need to update the given basemap tile index'''

    thisTileBounds = basemapInstance.getTileRectDegree(tileIndex)

    # Get the tile information from the HRSC image
    tileDict = hrscInstance.getTileInfo(thisTileBounds, tileIndex.getPostfix())
    #print 'Found these tile intersections:'
    #for hrscTile in tileDict.itervalues():
    #    print hrscTile['prefix']

    return tileDict
    
    

def updateTileWithHrscImage(hrscTileInfoDict, outputTilePath, tileLogPath):
    '''Update a single output tile with the given HRSC image'''

    # TODO: This C++ program can do multiple tiles in one call.

    # For each tile...
    hrscTiles = ''
    for hrscTile in hrscTileInfoDict.itervalues():    
        #try:
        cmd = ('./hrscMosaic ' + outputTilePath +' '+ outputTilePath +' '+ hrscTile['newColorPath'] +' '+
                                  hrscTile['tileMaskPath'] +' '+ hrscTile['tileToTileTransformPath'])
        MosaicUtilities.cmdRunner(cmd, outputTilePath, True)
        #raise Exception('DEBUG')
        hrscTiles += hrscTile['prefix'] + ', '

    # Return the path to log the success to
    return (tileLogPath, hrscTiles)
    
    


def updateTilesContainingHrscImage(basemapInstance, hrscInstance, pool=None):
    '''Updates all output tiles containing this HRSC image'''

    logger = logging.getLogger('MainProgram')

    # Find all the output tiles that intersect with this
    outputTilesRect = getCoveredOutputTiles(basemapInstance, hrscInstance)

    hrscSetName = hrscInstance.getSetName()
    mainLogPath = basemapInstance.getMainLogPath()
    
    # Skip this function if we have completed adding this HRSC image
    if basemapInstance.checkLog(mainLogPath, hrscSetName):
        logger.info('Have already completed adding HRSC image ' + hrscSetName + ',  skipping it.')
        return
    
    logger.info('Started updating tiles for HRSC image ' + hrscSetName)

    logger.info('Found overlapping output tiles:  ' + str(outputTilesRect))
    if pool:
        logger.info('Initializing tile output tasks...')
    
    # Loop through all the tiles
    tileResults = []
    for row in range(outputTilesRect.minY, outputTilesRect.maxY):
        for col in range(outputTilesRect.minX, outputTilesRect.maxX):
    
            # Set up the til information
            tileIndex  = MosaicUtilities.TileIndex(row, col) #basemapInstance.getTileIndex(98.5, -27.5)
            tileBounds = basemapInstance.getTileRectDegree(tileIndex)
            
            logger.info('Using HRSC image ' + hrscSetName + ' to update tile: ' + str(tileIndex))
            logger.info('--> Tile bounds = ' + str(tileBounds))

            #logger.info('\nMaking sure basemap info is present...')
            
            # Now that we have selected a tile, generate all of the tile images for it.
            # - The first time this is called for a tile it generates the backup image for the tile.
            (smallTilePath, largeTilePath, grayTilePath, outputTilePath, tileLogPath) =  \
                        basemapInstance.generateTileImages(tileIndex, False)
        
            #print '\nPasting on HRSC tiles...'

            # Have we already written this HRSC image to this tile?
            comboAlreadyWritten = basemapInstance.checkLog(tileLogPath, hrscSetName)
            if comboAlreadyWritten:
                logger.info('-- Skipping already written tile!') #Don't want to double-write the same image.
                continue
        
            # Get information about which HRSC tiles to paste on to the basemap
            hrscTileInfoDict = getHrscTileUpdateDict(basemapInstance, tileIndex, hrscInstance)
            if not hrscTileInfoDict: # If there are no tiles to use, move on to the next output tile!
                continue
        
            # Update the selected tile with the HRSC image
            if pool:
                # Send the function and arguments to the thread pool
                dictCopy = copy.copy(hrscTileInfoDict)
                tileResults.append(pool.apply_async(updateTileWithHrscImage,
                                                    args=(dictCopy, outputTilePath, tileLogPath)))
            else: # Just run the function
                updateTileWithHrscImage(hrscTileInfoDict, outputTilePath, tileLogPath)
            
            # DEBUG breaks
            #break
        #break


    if pool: # Wait for all the tasks to complete
        logger.info('Finished initializing tile output tasks.')
        logger.info('Waiting for tile processes to complete...')
        for result in tileResults:
            # Each task finishes by returning the log path for that tile.
            # - Record that we have used this HRSC/tile combination.
            # - This requires that tiles with no HRSC tiles do not get assigned a task.
            (tileLogPath, hrscTilePrefixList) = result.get()
            basemapInstance.updateLog(tileLogPath, hrscSetName, hrscTilePrefixList)
            
            
        logger.info('All tile writing processes have completed')

    #raise Exception('DEBUG')
        
    # Log the fact that we have finished adding this HRSC image    
    basemapInstance.updateLog(mainLogPath, hrscSetName)

    logger.info('Finished updating tiles for HRSC image ' + hrscSetName)

#-----------------------------------------------------------------------------------------

# Laptop
#testDirectory    = '/home/smcmich1/data/hrscMapTest/'
#fullBasemapPath  = testDirectory + 'projection_space_basemap.tif'
#sourceHrscFolder = testDirectory + 'external_data'
#hrscOutputFolder = testDirectory + 'hrscFiles'
#outputTileFolder = testDirectory + 'outputTiles'
#databasePath     = 'FAIL'

# Lunokhod 2
fullBasemapPath  = '/byss/smcmich1/data/hrscBasemap/projection_space_basemap.tif'
sourceHrscFolder = '/home/smcmich1/data/hrscDownloadCache'
hrscOutputFolder = '/home/smcmich1/data/hrscProcessedFiles'
outputTileFolder = '/byss/smcmich1/data/hrscBasemap/outputTiles_64'
databasePath     = '/byss/smcmich1/data/google/googlePlanetary.db'
kmlPyramidFolder = '/byss/docroot/smcmich1/hrscMosaicKml'
kmlPyramidWebAddress = 'http://byss.arc.nasa.gov/smcmich1/hrscMosaicKml/0.kml'

print 'Starting basemap enhancement script...'

logger = logging.getLogger('MainProgram')

# Echo logging to stdout
echo = logging.StreamHandler(sys.stdout)
echo.setLevel(logging.DEBUG)
echo.setFormatter(logging.Formatter(LOG_FORMAT_STR))
logger.addHandler(echo)


# Initialize the multi-threading worker pool
processPool  = None
if NUM_PROCESS_THREADS > 1:
    processPool = multiprocessing.Pool(processes=NUM_PROCESS_THREADS)

logger.info('==== Initializing the base map object ====')
basemapInstance = mosaicTileManager.MarsBasemap(fullBasemapPath, outputTileFolder)
mainLogPath = basemapInstance.getMainLogPath()
logger.info('--- Finished initializing the base map object ---\n')


# Get a list of the HRSC images we are testing with
#fullImageList = getHrscImageList()
tempFileFinder = hrscFileCacher.HrscFileCacher(databasePath, sourceHrscFolder, BAD_HRSC_FILE_PATH)
fullImageList = tempFileFinder.getHrscSetList(HRSC_FETCH_ROI)
tempFileFinder = None # Delete this temporary object

logger.info('Identified ' + str(len(fullImageList)) + ' HRSC images in the requested region:\n'+
            str(fullImageList))

# DEBUG --> Test this image!
#fullImageList = fullImageList[6:8]
#print fullImageList

## Prune out all the HRSC images that we have already added to the mosaic.
#hrscImageList = []
#for hrscSetName in fullImageList:
#    if False:#basemapInstance.checkLog(mainLogPath, hrscSetName):
#        logger.info('Have already completed adding HRSC image ' + hrscSetName + ',  skipping it.')
#    else:
#        hrscImageList.append(hrscSetName)

hrscImageList = fullImageList

# Restrict the image list to the batch size
# - It would be more accurate to only count valid images but this is good enough
hrscImageList = hrscImageList[0:IMAGE_BATCH_SIZE]   #['h3276_0000']
batchName     = hrscImageList[0] # TODO: Assign the batches a number.
print 'Image list for this batch: ' + str(hrscImageList)

## Set up the HRSC file manager object
logger.info('Starting communication queues')
#hrscFileFetcher = hrscFileCacher.HrscFileCacher(databasePath, sourceHrscFolder, downloadPool)
downloadCommandQueue  = multiprocessing.Queue()
downloadResponseQueue = multiprocessing.Queue()
logger.info('Initializing HRSC file caching thread')
downloadThread = threading.Thread(target=cacheManagerThreadFunction,
                                  args  =(databasePath, sourceHrscFolder,            
                                          downloadCommandQueue, downloadResponseQueue)
                                 )
downloadThread.daemon = True # Needed for ctrl-c to work
logger.info('Running thread...')
downloadThread.start()


# Go ahead and send a request to fetch the first HRSC image
logger.info('Sending FETCH command: ' + hrscImageList[0])
downloadCommandQueue.put('FETCH ' + hrscImageList[0])


# Loop through input HRSC images
numHrscDataSets = len(hrscImageList) 
for i in range(0,numHrscDataSets): 
    
    # Get the name of this and the next data set
    hrscSetName = hrscImageList[i]
    nextSetName = None
    if i < numHrscDataSets-1:
        nextSetName = hrscImageList[i+1]
        # Go ahead and submit the fetch request for the next set name.
        logger.info('Sending FETCH command: ' + nextSetName)
        downloadCommandQueue.put('FETCH ' + nextSetName)

    # Notes on downloading:
    # - Each iteration of this loop commands one download, and waits for one download.
    # - The queues keep things in order, and the download thread handles one data set at a time.
    # - This means that we can have one data set downloading while one data set is being processed.
    # - The next improvement to be made would be to download multiple data sets at the same time.

   
    ## Pick a location to store the data for this HRSC image
    thisHrscFolder = os.path.join(hrscOutputFolder, hrscSetName)

    #try:

    logger.info('=== Fetching HRSC image ' + hrscSetName + ' ===')

    # Fetch the HRSC data from the web
    #hrscFileInfoDict = hrscFileFetcher.fetchHrscDataSet(hrscSetName)
    hrscFileInfoDict = downloadResponseQueue.get() # Wait for the parallel thread to provide the data
    if not 'setName' in hrscFileInfoDict:
        raise Exception('Ran out of HRSC files, processing stopped!!!')
    if hrscFileInfoDict['setName'] != hrscSetName:
        raise Exception('Set fetch mismatch!  Expected %s, got %s instead!' % 
                         (hrscSetName, hrscFileInfoDict['setName']))
    logger.info('Received fetch information for ' + hrscSetName)

    if 'error' in hrscFileInfoDict:
        logger.info('Skipping data set ' + hrscSetName + ' which could not be fetched.')
        continue

    #print 'SKIPPING IMAGE PROCESSING!!!'
    #continue

    logger.info('\n=== Initializing HRSC image ' + hrscSetName + ' ===')

    # Preprocess the HRSC image
    hrscInstance = hrscImageManager.HrscImage(hrscFileInfoDict, thisHrscFolder, basemapInstance, False, processPool)

    # TODO: Need to make the HRSC manager clean up the processed folder too!

    logger.info('--- Now initializing high res HRSC content ---')

    # Complete the high resolution components
    hrscInstance.prepHighResolutionProducts()
    
    logger.info('--- Finished initializing HRSC image ---\n')

    #raise Exception('DEBUG')

    #continue # DEBUG - Just update the registration

    # Call the function to update all the output images for this HRSC image
    updateTilesContainingHrscImage(basemapInstance, hrscInstance, processPool)

    logger.info('<<<<< Finished writing all tiles for this HRSC image! >>>>>')

    #raise Exception('DEBUG')

PROCESS_POOL_KILL_TIMEOUT = 180 # The pool should not be doing any work at this point!
if processPool:
    logger.info('Cleaning up the processing thread pool...')
    processPool.close()
    processPool.join(PROCESS_POOL_KILL_TIMEOUT)

downloadCommandQueue.put('STOP') # Stop the download thread
downloadThread.join()


# TODO: Batch cleanup stuff!

# Generate a KML pyramid of the tiles for diagnostics
stackImagePyramid.main(outputTileFolder, kmlPyramidFolder)

# Send a message notifiying that the output needs to be reviewed!
msgText = '''
KML pyramid link:
'''+kmlPyramidWebAddress+'''
To undo the tile changes:
cp -r '''+basemapInstance.getBackupFolder()+' '+outputTileFolder+''' 

TODO: Also need to update the input log files in the output folder!

To accept the tile changes:
rm * '''+basemapInstance.getBackupFolder()+''''

To start the next batch:
source /byss/smcmich1/run_hrsc_basemap_script.sh
'''
MosaicUtilities.sendEmail('scott.t.mcmichael@nasa.gov', 
                          'HRSC map batch '+batchName+' completed',
                          msgText)

logger.info('Basemap generation script completed!')

