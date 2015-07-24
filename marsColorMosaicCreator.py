
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
import time
import traceback

import IrgGeoFunctions
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

NUM_DOWNLOAD_THREADS = 5 # There are five files we download per data set
NUM_PROCESS_THREADS  = 20


IMAGE_BATCH_SIZE = 1 # This should be set equal to the HRSC cache size



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
#HRSC_FETCH_ROI = MosaicUtilities.Rectangle(133.0, 142.0, 46, 50.0) # Viking 2 lander region
#HRSC_FETCH_ROI = MosaicUtilities.Rectangle(-78.0, -63.0, -13.0, -2.5) # Candor Chasma region
HRSC_FETCH_ROI = MosaicUtilities.Rectangle(-161.0, -154.0, -60.0, -50.0) # Region near -60 lat

#-----------------------------------------------------------------------------------------
# Functions

def getDiskUsage():
    '''Return simpe disk space usage information'''
    cmd = ['df', '-h']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()
    return textOutput


def cacheManagerThreadFunction(databasePath, hrscDownloadFolder, hrscProcessedFolder, inputQueue, outputQueue):
    '''Thread to allow downloading of HRSC data in parallel with image processing.
       The input queue recieves three types of commands:
           "STOP" --> Finish current tasks, then exit.
           "KILL" --> Immediately kill all the threads and exit.
           "FETCH data_set_name" --> Fetch the specified data set
           "FINISHED data_set_name" --> Signals that this data set can be safely deleted
'''

    logger = logging.getLogger('DownloadThread')

    # Echo logging to stdout
    echo = logging.StreamHandler(sys.stdout)
    echo.setLevel(logging.DEBUG)
    echo.setFormatter(logging.Formatter(LOG_FORMAT_STR))
    logger.addHandler(echo)



    # Initialize a process pool to be managed by this thread
    downloadPool = None
    if NUM_DOWNLOAD_THREADS > 1:
        downloadPool = multiprocessing.Pool(processes=NUM_DOWNLOAD_THREADS)

    # Set up the HRSC file manager object
    logger.info('Initializing HRSC file caching object')
    hrscFileFetcher = hrscFileCacher.HrscFileCacher(databasePath, hrscDownloadFolder, hrscProcessedFolder,
                                                    BAD_HRSC_FILE_PATH, downloadPool)

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
                logger.error('Caught exception fetching data set ' + dataSet + '\n' + 
                             str(e) + '\n' + str(sys.exc_info()[0]) + '\n')
                logger.error(sys.exc_info()[0])
                print(traceback.format_exc())
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



def getCoveredOutputTiles(basemapInstance, hrscInstance):
    '''Return a bounding box containing all the output tiles covered by the HRSC image'''
    
    hrscBoundingBoxDegrees = hrscInstance.getBoundingBoxDegrees()
    #print 'HRSC BB = ' + str(hrscBoundingBoxDegrees)
    
    # Expand the computed bounding box a little bit to insure that
    #  we don't miss any tiles around the edges.
    BUFFER_SIZE = 0.1 # BB buffer size in degrees
    hrscBoundingBoxDegrees.expand(BUFFER_SIZE, BUFFER_SIZE)
    
    
    # DEBUG!  Restrict to a selected area.
    hrscBoundingBoxDegrees = HRSC_FETCH_ROI.getIntersection(hrscBoundingBoxDegrees)
    #print 'overlap BB = ' + str(hrscBoundingBoxDegrees)
    
    intersectRect = basemapInstance.getIntersectingTiles(hrscBoundingBoxDegrees)
    return intersectRect
    

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

    # Append all the tiles into one big command line call
    hrscTiles = ''
    cmd = './hrscMosaic ' + outputTilePath +' '+ outputTilePath
    for hrscTile in hrscTileInfoDict.itervalues():    
        #try:
        # This pastes the HRSC tile on top of the current output tile.  Another function
        #  will have made sure the correct output tile is in place.
        cmd += (' '+ hrscTile['newColorPath'] +' '+
                  hrscTile['tileMaskPath'] +' '+ hrscTile['tileToTileTransformPath'])
        hrscTiles += hrscTile['prefix'] + ', '

    # Execute the command line call
    MosaicUtilities.cmdRunner(cmd, outputTilePath, True)

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
    
    # Do all the basemap calls first using the pool before doing the HRSC work
    # - This will make sure that the proper file exists in the output directory to 
    #   paste incoming HRSC tiles on top of.
    logger.info('Making sure we have required basemap tiles for HRSC image ' + hrscSetName)
    basemapInstance.generateMultipleTileImages(outputTilesRect, pool, force=False)
    
    if pool:
        logger.info('Initializing tile output tasks...')
    
    # Loop through all the tiles
    tileResults = []
    for row in range(outputTilesRect.minY, outputTilesRect.maxY):
        for col in range(outputTilesRect.minX, outputTilesRect.maxX):
    
            # Set up the til information
            tileIndex  = MosaicUtilities.TileIndex(row, col)
            tileBounds = basemapInstance.getTileRectDegree(tileIndex)
            
            logger.info('Using HRSC image ' + hrscSetName + ' to update tile: ' + str(tileIndex))
            logger.info('--> Tile bounds = ' + str(tileBounds))

            #logger.info('\nMaking sure basemap info is present...')
            
            # Retrieve the needed paths for this tile
            (smallTilePath, largeTilePath, grayTilePath, outputTilePath, tileLogPath, tileBackupPath) = \
                basemapInstance.getPathsForTile(tileIndex) 
        
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

#================================================================================

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
outputTileFolder = '/byss/smcmich1/data/hrscBasemap/outputTiles_128'
backupFolder     = '/byss/smcmich1/data/hrscBasemap/output_tile_backups'
databasePath     = '/byss/smcmich1/data/google/googlePlanetary.db'
kmlPyramidFolder = '/byss/docroot/smcmich1/hrscMosaicKml'


# --- Folder notes ---
# - sourceHrscFolder holds the downloaded and preprocessed HRSC data
# - hrscOutputFolder holds the fully processed HRSC files
# - The current crop of tiles is written to outputTileFolder
# - The persistent set of final output tiles is kept in backupFolder

#================================================================================


print 'Starting basemap enhancement script...'

startTime = time.time()

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
basemapInstance = mosaicTileManager.MarsBasemap(fullBasemapPath, outputTileFolder, backupFolder)
basemapInstance.copySupportFilesFromBackupDir() # Copies the main log from the backup dir to output dir
basemapInputsUsedLog = basemapInstance.getMainLogPath()
logger.info('--- Finished initializing the base map object ---\n')


# Get a list of the HRSC images we are testing with
#fullImageList = getHrscImageList()
tempFileFinder = hrscFileCacher.HrscFileCacher(databasePath, sourceHrscFolder, 
                                               hrscOutputFolder, BAD_HRSC_FILE_PATH)

# Run-once code to find all the incomplete data sets in one pass
#fullImageList = tempFileFinder.getHrscSetList()
#tempFileFinder.findIncompleteSets(fullImageList)
#raise Exception('DONE FINDING BAD SETS')


fullImageList = tempFileFinder.getHrscSetList(HRSC_FETCH_ROI)
tempgFileFinder = None # Delete this temporary object

logger.info('Identified ' + str(len(fullImageList)) + ' HRSC images in the requested region:\n'+
            str(fullImageList))


# Prune out all the HRSC images that we have already added to the mosaic.
hrscImageList = []
for hrscSetName in fullImageList:
    if basemapInstance.checkLog(basemapInputsUsedLog, hrscSetName):
        logger.info('Have already completed adding HRSC image ' + hrscSetName + ',  skipping it.')
    else:
        hrscImageList.append(hrscSetName)
hrscImageList = ['h2410_0000'] # DEBUG

# Restrict the image list to the batch size
# - It would be more accurate to only count valid images but this is good enough
hrscImageList = hrscImageList[0:IMAGE_BATCH_SIZE]   #['h3276_0000']
batchName     = hrscImageList[0] # TODO: Assign the batches a number.
print 'Image list for this batch: ' + str(hrscImageList)

# Set up the HRSC file manager thread
logger.info('Starting communication queues')
downloadCommandQueue  = multiprocessing.Queue()
downloadResponseQueue = multiprocessing.Queue()
logger.info('Initializing HRSC file caching thread')
downloadThread = threading.Thread(target=cacheManagerThreadFunction,
                                  args  =(databasePath, sourceHrscFolder, hrscOutputFolder,       
                                          downloadCommandQueue, downloadResponseQueue)
                                 )
downloadThread.daemon = True # Needed for ctrl-c to work
logger.info('Running thread...')
downloadThread.start()


# Go ahead and send a request to fetch the first HRSC image
logger.info('Sending FETCH command: ' + hrscImageList[0])
downloadCommandQueue.put('FETCH ' + hrscImageList[0])


# Loop through input HRSC images
numHrscDataSets   = len(hrscImageList) 
processedDataSets = []
failedDataSets    = []
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

    try:

        logger.info('=== Fetching HRSC image ' + hrscSetName + ' ===')

        # Fetch the HRSC data from the web
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

        logger.info('\n=== Initializing HRSC image ' + hrscSetName + ' ===')

        # Preprocess the HRSC image
        hrscInstance = hrscImageManager.HrscImage(hrscFileInfoDict, thisHrscFolder, basemapInstance, False, processPool)

        logger.info('--- Now initializing high res HRSC content ---')

        # Complete the high resolution components
        hrscInstance.prepHighResolutionProducts()
        
        logger.info('--- Finished initializing HRSC image ---\n')


        # Call the function to update all the output images for this HRSC image
        updateTilesContainingHrscImage(basemapInstance, hrscInstance, processPool)

        logger.info('<<<<< Finished writing all tiles for this HRSC image! >>>>>')
        
        # Record that we finished processing this HRSC image
        processedDataSets.append( (hrscSetName, hrscInstance.getBoundingBoxDegrees()) )

    except Exception, e:
        # When we fail to fetch a data set, log a failure message and keep going.
        failedDataSets.append(hrscSetName)
        logger.error('Caught exception processing data set ' + hrscSetName + '\n' + 
                     str(e) + '\n' + str(sys.exc_info()[0]) + '\n')
        logger.error(traceback.format_exc())
            

numHrscImagesProcessed = len(processedDataSets)

PROCESS_POOL_KILL_TIMEOUT = 180 # The pool should not be doing any work at this point!
if processPool:
    logger.info('Cleaning up the processing thread pool...')
    processPool.close()
    processPool.join()

downloadCommandQueue.put('STOP') # Stop the download thread
downloadThread.join()

# Compute the run time for the output message
SECONDS_TO_HOURS = 1.0 / (60.0*60.0)
stopTime = time.time()
runTime  = (stopTime - startTime) * SECONDS_TO_HOURS

if numHrscImagesProcessed > 0:
    # Generate a KML pyramid of the tiles for diagnostics
    kmlPyramidLocalPath  = stackImagePyramid.main(outputTileFolder, kmlPyramidFolder, processedDataSets)
    pos                  = kmlPyramidLocalPath.find('/smcmich1')
    kmlPyramidWebAddress = 'http://byss.arc.nasa.gov' + kmlPyramidLocalPath[pos:]

    # Send a message notifiying that the output needs to be reviewed!
    msgText = '''
Finished processing ''' +str(numHrscImagesProcessed) + ''' HRSC images!

elapsed time = ''' + str(runTime) + ''' hours.

KML pyramid link:
'''+kmlPyramidWebAddress+'''
To undo the tile changes:
rm  '''+outputTileFolder+'''/*

To accept the tile changes:
rsync  -urv '''+outputTileFolder +' '+ backupFolder+'''
rm  '''+outputTileFolder+'''/*

To start the next batch, run:
/byss/smcmich1/run_hrsc_basemap_script.sh

Disk usage info:
''' + getDiskUsage()+'''
Processed image list:
'''
    for i in processedDataSets:
        msgText += i[0] + '\n'
    msgText += '\n Failed image list:\n'
    for i in failedDataSets:
        msgText += i + '\n'
else:
    msgText = '''ERROR: No HRSC images in the batch could be processed!\n''' + str(failedDataSets)
    
    
MosaicUtilities.sendEmail('scott.t.mcmichael@nasa.gov', 
                          'HRSC map batch '+batchName+' completed',
                          msgText)

logger.info('Basemap generation script completed!')

