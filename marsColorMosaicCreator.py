
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
import shutil
import optparse

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
- Between +/-60 degrees, 86x256 tiles = 22016 tiles (actually between +/-60.46875)
With 32x  increase, each tile is 1440x1440,   ~6MB,  190GB, mpp = ~58
With 64x  increase, each tile is 2880x2880,  ~24MB,  760GB, mpp = ~29
With 128x increase, each tile is 5760x5760,  ~95MB,    3TB, mpp = ~14.5 <-- Using this!
With 160x increase, each tile is 7200x7200, ~150MB,  4.6TB, mpp = ~11.6



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

# Control how many threads are in the thread pools.
NUM_DOWNLOAD_THREADS = 5 # There are five files we download per data set

# This limits the program to parsing this many HRSC files before stopping.
IMAGE_BATCH_SIZE = 1 # This should be set equal to the HRSC cache size




#-----------------------------------------------------------------------------------------
# Functions

def getDiskUsage():
    '''Return simple disk space usage information'''
    cmd = ['df', '-h']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()
    return textOutput


def cacheManagerThreadFunction(DATABASE_PATH, hrscDownloadFolder, hrscProcessedFolder, inputQueue, outputQueue):
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
    echo.setFormatter(logging.Formatter(MosaicUtilities.LOG_FORMAT_STR))
    logger.addHandler(echo)



    # Initialize a process pool to be managed by this thread
    downloadPool = None
    if NUM_DOWNLOAD_THREADS > 1:
        downloadPool = multiprocessing.Pool(processes=NUM_DOWNLOAD_THREADS)

    # Set up the HRSC file manager object
    logger.info('Initializing HRSC file caching object')
    hrscFileFetcher = hrscFileCacher.HrscFileCacher(DATABASE_PATH, hrscDownloadFolder, hrscProcessedFolder,
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

    # We only get here when we break out of the main loop
    outputQueue.put('STOPPED')
    logger.info('Download manager thread stopped.')



def getCoveredOutputTiles(basemapInstance, hrscInstance):
    '''Return a bounding box containing all the output tiles covered by the HRSC image'''
    
    ## This bounding box can be in either the +/-180 range or the 0-360 range
    hrscBoundingBoxProjected = hrscInstance.getBoundingBoxProjected()
    print 'HRSC BB = ' + str(hrscBoundingBoxProjected)
    
    # Expand the computed bounding box a little bit to insure that
    #  we don't miss any tiles around the edges.
    BUFFER_SIZE = 5000 # BB buffer size in meters
    hrscBoundingBoxProjected.expand(BUFFER_SIZE, BUFFER_SIZE)

    # DEBUG!  Restrict to a selected area.
    hrscBoundingBoxProjected = HRSC_FETCH_PROJ_ROI.getIntersection(hrscBoundingBoxProjected)

    intersectTileList = basemapInstance.getIntersectingTiles(hrscBoundingBoxProjected)
    return intersectTileList
    

def updateTileWithHrscImage(hrscTileInfoDict, outputTilePath, tileLogPath):
    '''Update a single output tile with the given HRSC image'''

    # Write the output to a new temporary file in case we wreck it!
    tempFilePath = outputTilePath + '_temp.tif'

    # Append all the tiles into one big command line call
    hrscTiles = ''
    cmd = './hrscMosaic ' + outputTilePath +' '+ tempFilePath
    for hrscTile in hrscTileInfoDict.itervalues():    

        # This pastes the HRSC tile on top of the current output tile.  Another function
        #  will have made sure the correct output tile is in place.
        cmd += (' '+ hrscTile['newColorPath'] +' '+
                  hrscTile['tileMaskPath'] +' '+ hrscTile['tileToTileTransformPath'])
        hrscTiles += hrscTile['prefix'] + ', '

    # Execute the command line call
    MosaicUtilities.cmdRunner(cmd, tempFilePath, True)

    # Make sure that the command did not ruin the output file!
    if MosaicUtilities.isImageFileValid(tempFilePath):
        # If we didn't ruin the output file, copy it to the proper location.
        shutil.copyfile(tempFilePath, outputTilePath)
        os.remove(tempFilePath)
    else:
        # On a failure, just log the error so that we can keep going.
        logger = logging.getLogger('MainProgram')
        logger.error('Running command:  ' + cmd + '\n' +
                     'has broken file: ' + tempFilePath)
        return (tileLogPath, 'FAILED')

    # Return the path to log the success to
    return (tileLogPath, hrscTiles)
    
    


def updateTilesContainingHrscImage(basemapInstance, hrscInstance, pool=None):
    '''Updates all output tiles containing this HRSC image'''

    logger = logging.getLogger('MainProgram')

    # Find all the output tiles that intersect with this
    outputTilesList = getCoveredOutputTiles(basemapInstance, hrscInstance)

    #for t in outputTilesList:
    #    print t
    #raise Exception('DEBUG')

    hrscSetName = hrscInstance.getSetName()
    mainLogPath = basemapInstance.getMainLogPath()
    
    # Skip this function if we have completed adding this HRSC image
    if basemapInstance.checkLog(mainLogPath, hrscSetName):
        logger.info('Have already completed adding HRSC image ' + hrscSetName + ',  skipping it.')
        basemapInstance.updateLog(mainLogPath, hrscSetName) # This should have already been logged!
        return

    logger.info('Started updating tiles for HRSC image ' + hrscSetName)
    #logger.info('Found overlapping output tiles:  ' + str(outputTilesList))
    
    # Do all the basemap calls first using the pool before doing the HRSC work
    # - This will make sure that the proper file exists in the output directory to 
    #   paste incoming HRSC tiles on top of.
    logger.info('Making sure we have required basemap tiles for HRSC image ' + hrscSetName)
    basemapInstance.generateMultipleTileImages(outputTilesList, pool, force=False)
    #basemapInstance.generateMultipleTileImages(outputTilesList, None, force=False) # Single thread debug

    if pool:
        logger.info('Initializing tile output tasks...')
    
    # Loop through all the tiles
    tileResults = []
    for tileIndex in outputTilesList:

        tileBounds = basemapInstance.getTileRectProjected(tileIndex)
                
        # Retrieve the needed paths for this tile
        (smallTilePath, largeTilePath, grayTilePath, outputTilePath, tileLogPath, tileBackupPath) = \
            basemapInstance.getPathsForTile(tileIndex) 
    
        if not os.path.exists(tileBackupPath):
            print 'Skipping non-existant output tile ' + str(tileIndex)
            continue  # Skip non-existant basemap tiles
    
        # Have we already written this HRSC image to this tile?
        comboAlreadyWritten = basemapInstance.checkLog(tileLogPath, hrscSetName)
        if comboAlreadyWritten: #Don't want to double-write the same image.
            logger.info('-- Skipping already written tile: ' + str(tileIndex)) 
            continue

        logger.info('Using HRSC image ' + hrscSetName + ' to update tile: ' + str(tileIndex))
        logger.info('--> Tile bounds = ' + str(tileBounds))
        #logger.info('\nMaking sure basemap info is present...')
    
        # Get information about which HRSC tiles to paste on to the basemap
        hrscTileInfoDict = hrscInstance.getTileInfo(basemapInstance, tileBounds, tileIndex.getPostfix())
        if not hrscTileInfoDict: # If there are no HRSC tiles to use, move on to the next output tile!
            continue
    
        # Update the selected tile with the HRSC image
        if pool:
            # Send the function and arguments to the thread pool
            dictCopy = copy.copy(hrscTileInfoDict)
            tileResults.append(pool.apply_async(updateTileWithHrscImage,
                                                args=(dictCopy, outputTilePath, tileLogPath)))
        else: # Just run the function
            updateTileWithHrscImage(hrscTileInfoDict, outputTilePath, tileLogPath)
        
        #print 'DEBUG - only updating one tile!'
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
        
    # Log the fact that we have finished adding this HRSC image
    print 'Log path: #' + mainLogPath+ '#'
    print 'Set name: #' + hrscSetName+ '#'
    basemapInstance.updateLog(mainLogPath, hrscSetName)

    logger.info('Finished updating tiles for HRSC image ' + hrscSetName)


# TODO: Need a polar version of this!

def generateAllUpsampledBasemapTiles(basemapInstance, pool):
    '''Generate all the basemap tiles from the input low-res image.
       There are 128*256 = 32,768 tiles in the full image, about 20,000 tiles
       in the +/-60 version.'''

    print 'GENERATING ALL BASEMAP TILES'

    # Get all the tiles we are interested in.
    # Should have tiles from -180 to 180 in the +/-60 range. (on /byss)
    allTileList = basemapInstance.getIntersectingTiles(HRSC_FETCH_PROJ_ROI)
    
    # Make all the outputs.  This will take a while!
    basemapInstance.generateMultipleTileImages(allTileList, pool, force=False)
    

#================================================================================


def mainProcessingFunction(options):

    print 'Starting basemap enhancement script...'

    startTime = time.time()

    logger = logging.getLogger('MainProgram')

    # Echo logging to stdout
    echo = logging.StreamHandler(sys.stdout)
    echo.setLevel(logging.DEBUG)
    echo.setFormatter(logging.Formatter(MosaicUtilities.LOG_FORMAT_STR))
    logger.addHandler(echo)


    # Initialize the multi-threading worker pool
    processPool  = None
    if options.numThreads > 1:
        processPool = multiprocessing.Pool(processes=options.numThreads)

    logger.info('==== Initializing the base map object ====')
    basemapInstance = mosaicTileManager.MarsBasemap(FULL_BASEMAP_PATH, NEW_OUTPUT_TILE_FOLDER, BACKUP_FOLDER,
                                                    MosaicUtilities.PROJ_TYPE_NORMAL)
    basemapInstance.copySupportFilesFromBackupDir() # Copies the main log from the backup dir to output dir
    basemapInputsUsedLog = basemapInstance.getMainLogPath()
    print 'CHECKING LOG PATH: ' + basemapInputsUsedLog

    # Open a seperate text file to record the names of failed input data sets
    failedSetsLogPath = os.path.join(NEW_OUTPUT_TILE_FOLDER, 'failed_input_list.txt')
    failedSetsLogFile = open(failedSetsLogPath, 'a')

    # Create additional basemap instances.
    # - The 180 degree centered map is never for creating output tiles, only for
    #   registration and preprocessing of wraparound HRSC images.
    # - The other two basemaps have their own set of tiles and are used in place of the normal basemap.
    dummyFolder = '/dev/null'
    basemapInstance180   = mosaicTileManager.MarsBasemap(FULL_BASEMAP_PATH180,    
                                                         dummyFolder, dummyFolder, 
                                                         MosaicUtilities.PROJ_TYPE_360)
    basemapInstanceNorth = mosaicTileManager.MarsBasemap(FULL_BASEMAP_PATH_NORTH, 
                                                         NEW_OUTPUT_TILE_FOLDER, BACKUP_FOLDER_NORTH,
                                                         MosaicUtilities.PROJ_TYPE_NORTH_POLE)
    basemapInstanceSouth = mosaicTileManager.MarsBasemap(FULL_BASEMAP_PATH_SOUTH, 
                                                         NEW_OUTPUT_TILE_FOLDER, BACKUP_FOLDER_SOUTH,
                                                         MosaicUtilities.PROJ_TYPE_SOUTH_POLE)
    basemapForOutputTiles = basemapInstance
    if options.mapType == MosaicUtilities.PROJ_TYPE_NORTH_POLE:
        basemapForOutputTiles = basemapInstanceNorth
    elif options.mapType == MosaicUtilities.PROJ_TYPE_SOUTH_POLE:
        basemapForOutputTiles = basemapInstanceSouth


    logger.info('--- Finished initializing the base map object ---\n')

    # Run once code to generate all of the starting basemap tiles!
    if options.genBasemapTiles:
        generateAllUpsampledBasemapTiles(basemapForOutputTiles, processPool)
        print 'DONE GENERATING ALL INPUT TILES'
        return False # Stop with no email message sent

    # Get a list of the HRSC images we are testing with
    tempFileFinder = hrscFileCacher.HrscFileCacher(DATABASE_PATH, HRSC_DOWNLOAD_FOLDER, 
                                                   HRSC_PROCESSING_FOLDER, BAD_HRSC_FILE_PATH)

    # Run-once code to find all the incomplete data sets in one pass
    #fullImageList = tempFileFinder.getHrscSetList()
    #tempFileFinder.findIncompleteSets(fullImageList)
    #raise Exception('DONE FINDING BAD SETS')


    fullImageList = tempFileFinder.getHrscSetList(HRSC_FETCH_DEG_ROI)
    tempgFileFinder = None # Delete this temporary object

    logger.info('Identified ' + str(len(fullImageList)) + ' HRSC images in the requested region.')


    # Prune out all the HRSC images that we have already added to the mosaic.
    hrscImageList = []
    for hrscSetName in fullImageList:
        if basemapForOutputTiles.checkLog(basemapInputsUsedLog, hrscSetName):
            logger.info('Have already completed adding HRSC image ' + hrscSetName + ',  skipping it.')
        else:
            hrscImageList.append(hrscSetName)
    #hrscImageList = ['h1167_0000'] # DEBUG

    numDataSetsRemainingToProcess = len(hrscImageList)
    logger.info('Num data sets remaining to process = ' + str(numDataSetsRemainingToProcess))

    # Restrict the image list to the batch size
    # - It would be more accurate to only count valid images but this is good enough
    hrscImageList = hrscImageList[0:IMAGE_BATCH_SIZE]
    try:
      batchName = hrscImageList[0]
    except:
      batchName = 'Default Name'
    logger.info('Image list for this batch: ' + str(hrscImageList))

    if len(hrscImageList) > 0:
      # Set up the HRSC file manager thread
      logger.info('Starting communication queues')
      downloadCommandQueue  = multiprocessing.Queue()
      downloadResponseQueue = multiprocessing.Queue()
      logger.info('Initializing HRSC file caching thread')
      downloadThread = threading.Thread(target=cacheManagerThreadFunction,
                                        args  =(DATABASE_PATH, HRSC_DOWNLOAD_FOLDER, HRSC_PROCESSING_FOLDER,       
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
    setProcessTimes   = []
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

       
        # Pick a location to store the data for this HRSC image
        thisHrscFolder = os.path.join(HRSC_PROCESSING_FOLDER, hrscSetName)

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

            #raise Exception('STOPPED AFTER DOWNLOAD')

            setStartTime = time.time()
            logger.info('\n=== Initializing HRSC image ' + hrscSetName + ' ===')

            # Preprocess the HRSC image
            hrscInstance = hrscImageManager.HrscImage(hrscFileInfoDict, thisHrscFolder,
                                                      basemapInstance, basemapInstance180,
                                                      basemapInstanceNorth, basemapInstanceSouth,
                                                      False, processPool, options.mapType)

            logger.info('--- Now initializing high res HRSC content ---')

            # Complete the high resolution components
            hrscInstance.prepHighResolutionProducts()
            
            logger.info('--- Finished initializing HRSC image ---\n')

            # Call the function to update all the output images for this HRSC image
            updateTilesContainingHrscImage(basemapForOutputTiles, hrscInstance, processPool)

            logger.info('<<<<< Finished writing all tiles for this HRSC image! >>>>>')
            
            # Record that we finished processing this HRSC image
            dataSetResults = (hrscSetName, hrscInstance.getBoundingBoxDegrees())
            processedDataSets.append(dataSetResults)

            # Back up debug thumbnails from this data set to a permanent location
            recordThumbnails(hrscSetName)

            # Record how long this data set took to process
            setStopTime = time.time()
            setProcessTimes.append(setStopTime - setStartTime)

        except Exception, e:
            # When we fail to fetch a data set, log a failure message and keep going.
            failedDataSets.append(hrscSetName)
            logger.error('Caught exception processing data set ' + hrscSetName + '\n' + 
                         str(e) + '\n' + str(sys.exc_info()[0]) + '\n')
            logger.error(traceback.format_exc())
            failedSetsLogFile.write(hrscSetName)
                

    numHrscImagesProcessed = len(processedDataSets)

    failedSetsLogFile.close() # Close the failure list log file

    PROCESS_POOL_KILL_TIMEOUT = 5 # The pool should not be doing any work at this point!
    if processPool:
        logger.info('Cleaning up the processing thread pool...')
        # Give the pool processes a little time to stop, them kill them.
        processPool.close()
        time.sleep(PROCESS_POOL_KILL_TIMEOUT)
        processPool.terminate()
        processPool.join()

    try:
        downloadCommandQueue.put('STOP') # Stop the download thread
        downloadThread.join()
    except:
        print 'Exception thrown shutting down the downloader!'

    #raise Exception('DEBUG - SKIP HTML GENERATION!!')

    # Call output reporting/logging functions

    sentEmail = generateTreeAndEmail(startTime, numHrscImagesProcessed, setProcessTimes, 
                                     processedDataSets, failedDataSets, numDataSetsRemainingToProcess, options)

    logger.info('Basemap generation script completed!')
    return sentEmail


def recordThumbnails(setName):
    '''Copy some thumbnail files to a centralized location for records and debugging'''

    # TODO: The HRSC manager should take care of this!
    try:
        # Copy the low res nadir image overlaid on a section of the low res mosaic
        debugImageInputPath = os.path.join(HRSC_PROCESSING_FOLDER, setName+'/'+setName+'_registration_debug_mosaic.tif')
        debugImageCopyPath  = os.path.join(OUTPUT_REGISTRATION_FOLDER, setName+'_registration_image.tif')
        shutil.copy(debugImageInputPath, debugImageCopyPath)

        # Copy the low res nadir image to a thumbnail folder
        debugImageInputPath = os.path.join(HRSC_PROCESSING_FOLDER, setName+'/'+setName+'_nd3_basemap_res.tif')
        debugImageCopyPath  = os.path.join(OUTPUT_THUMBNAIL_FOLDER, setName+'_nadir_thumbnail.tif')
        shutil.copy(debugImageInputPath, debugImageCopyPath)
    except:
        logger.error('Error copying a debug file for set ' + setName)


def generateTreeAndEmail(startTime, numHrscImagesProcessed, setProcessTimes, 
               processedDataSets, failedDataSets, numDataSetsRemainingToProcess, options):
    '''Generate an HRSC image pyramid and then send out a status email'''

    # Compute the run time for the output message
    SECONDS_TO_HOURS = 1.0 / (60.0*60.0)
    stopTime = time.time()
    runTime  = (stopTime - startTime) * SECONDS_TO_HOURS

    if (numHrscImagesProcessed <= 0) and (IMAGE_BATCH_SIZE != 0):
        # Brief message when no data was processed
        msgText = '''ERROR: No HRSC images in the batch could be processed!\n''' + str(failedDataSets)
    else: # Full data handling
        if not options.skipKmlPyramid:
            try:
                # Generate a KML pyramid of the tiles for diagnostics
                kmlPyramidLocalPath  = stackImagePyramid.main(NEW_OUTPUT_TILE_FOLDER, KML_PYRAMID_FOLDER, processedDataSets)
                pos                  = kmlPyramidLocalPath.find('/smcmich1')
                kmlPyramidWebAddress = 'http://byss.arc.nasa.gov' + kmlPyramidLocalPath[pos:]

                if options.uploadBucket:
                    # --> GSutil needs to be provided or on the path!
                    # -- gsutil also needs to be configured to upload to the correct bucket
                    # --gsutil-path /byss/smcmich1/programs/gsutil_install/gsutil
                    cmd = ('python ../sendToGoogleBucket.py sync-parallel  --dir '+KML_PYRAMID_FOLDER+
                              ' -p 1 --chunk-size 200 --gsutil-path /home/scott_mcmichael/MOUNT_POINT/programs/gsutil/gsutil')
                    if options.bucketPrefix:
                        cmd += (' --prepend-path '+options.bucketPrefix)
                    print cmd
                    os.system(cmd)

            except:
                logger.error('Error generating image pyramid!')
                kmlPyramidWebAddress = 'FAILED!'
        else:
            kmlPyramidWebAddress = 'Skipped production.'
            
        # Send a message notifiying that the output needs to be reviewed!
        msgText = '''
    Finished processing ''' +str(numHrscImagesProcessed) + ''' HRSC images!

    elapsed time = ''' + str(runTime) + ''' hours.

    Number remaining data sets = ''' + str(numDataSetsRemainingToProcess) + '''

    KML pyramid link:
    '''+kmlPyramidWebAddress+'''

    Registration debug images are here:
    '''+OUTPUT_REGISTRATION_FOLDER+'''
    --> To clear:
    rm '''+OUTPUT_REGISTRATION_FOLDER+'''/*.tif

    Image thumbnails are stored here:
    '''+OUTPUT_THUMBNAIL_FOLDER+'''
    Don't clear these!
    -------

    To undo the tile changes:
    rm  '''+NEW_OUTPUT_TILE_FOLDER+'''/*

    To accept the tile changes:
    rsync  --update --existing -avz '''+NEW_OUTPUT_TILE_FOLDER +'/ '+ BACKUP_FOLDER+'''/
    rm  '''+NEW_OUTPUT_TILE_FOLDER+'''/*

    To start the next batch, run:
    /byss/smcmich1/run_hrsc_basemap_script.sh

    Disk usage info:
    ''' + getDiskUsage()+'''
    Processed image list:
    '''
        index = 0
        for i in processedDataSets:
            msgText += i[0] + ' in ' + str(setProcessTimes[index]/60.0)+ ' minutes\n'
            index += 1
        if failedDataSets:
          msgText += '\n Failed image list:\n'
          for i in failedDataSets:
              msgText += i + '\n'
        
        
    MosaicUtilities.sendEmail('scott.t.mcmichael@nasa.gov', 
                              'HRSC map batch completed',
                              msgText)

    return True


    # Commands for generating the 180 centered image
    # gdal_translate projection_space_basemap.tif left.tif -srcwin 0 0 5760 5760
    # gdal_translate projection_space_basemap.tif right.tif -srcwin 5760 0 5760 5760
    # montage -mode Concatenate -tile 2x1 -background black  -depth 8  right.tif left.tif center180.tif
    #gdal_translate center180.tif projection_space_basemap180.tif -a_srs "+proj=eqc +lon_0=180 +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m" -a_ullr -10669477.100 5334738.600 10669477.100 -5334738.600


def setGlobalConfigs(argsIn):
    '''Parse the input arguments and set global variables'''
    
    # Define global variables we will use
    global FULL_BASEMAP_PATH
    global FULL_BASEMAP_PATH180
    global FULL_BASEMAP_PATH_NORTH
    global FULL_BASEMAP_PATH_SOUTH
    global BACKUP_FOLDER
    global BACKUP_FOLDER_NORTH
    global BACKUP_FOLDER_SOUTH
    global DATABASE_PATH
    global RUN_LOG_FOLDER
    
    global BAD_HRSC_FILE_PATH
    
    global NEW_OUTPUT_TILE_FOLDER
    global HRSC_DOWNLOAD_FOLDER
    global HRSC_PROCESSING_FOLDER
    
    global KML_PYRAMID_FOLDER
    global OUTPUT_THUMBNAIL_FOLDER
    global OUTPUT_REGISTRATION_FOLDER
    
    global HRSC_FETCH_DEG_ROI
    global HRSC_FETCH_PROJ_ROI
    
    
    # Input argument parsing
    usage = "usage: marsColorMosaicCreator.py [--help]\n"
    parser = optparse.OptionParser(usage=usage)
    
    parser.add_option('--generate-basemap-tiles', action='store_true', 
                      dest='genBasemapTiles', default=False,
                      help='Generate the basemap tiles instead of normal processing..')
    parser.add_option('--skip-kml-pyramid', action='store_true', 
                      dest='skipKmlPyramid', default=False,
                      help='Generate the basemap tiles instead of normal processing..')
    parser.add_option("--upload-bucket", dest="uploadBucket", default='',
                      help="If provided, upload the KML tree to a Google Cloud Storage Bucket.")
    parser.add_option("--threads", type="int", dest="numThreads", default=16,
                      help="Number of threads to use for processing.")
    parser.add_option("--safe-folder", dest="safeFolder", default='/byss/smcmich1/data',
                      help="Folder to store output files in.")
    parser.add_option("--volatile-folder", dest="volatileFolder", default='/home/smcmich1/data',
                      help="Folder to store temporary files in.")
    parser.add_option("--repo-folder", dest="repoFolder", default='/byss/smcmich1/repo',
                      help="Folder where the repository is installed.")
    parser.add_option("--map-type", dest="mapRaw", default='normal',
                      help="Type of map to make (normal, north, south).")
    parser.add_option("--bucket-prefix", dest="bucketPrefix", default='',
                      help="Optional upload bucket prefix to keep results seperate.")
    
    parser.add_option("--node-index", type="int", dest="nodeIndex", default=0,
                      help="Designated Google Compute processing node, sets region and bucket prefix.")
    
    (options, args) = parser.parse_args()
    
    # Parse out the desired map type
    if options.mapRaw == 'normal':
        options.mapType = MosaicUtilities.PROJ_TYPE_NORMAL
        raise Exception('DEBUG ONLY DO POLES!')
    elif options.mapRaw == 'north':
        options.mapType = MosaicUtilities.PROJ_TYPE_NORTH_POLE
    elif options.mapRaw == 'south':
        options.mapType = MosaicUtilities.PROJ_TYPE_SOUTH_POLE
    else:
        raise Exception('Map type not recognized!')
    
    SAFE_FOLDER     = options.safeFolder     # Permanent files go here
    VOLATILE_FOLDER = options.volatileFolder # Temporary files go here
    REPO_FOLDER     = options.repoFolder     # Source code folder
    
    FULL_BASEMAP_PATH       = os.path.join(SAFE_FOLDER, 'hrscBasemap/projection_space_basemap.tif')
    FULL_BASEMAP_PATH180    = os.path.join(SAFE_FOLDER, 'hrscBasemap180/projection_space_basemap180.tif')
    FULL_BASEMAP_PATH_NORTH = os.path.join(SAFE_FOLDER, 'hrscBasemapNorth/map_north_pole_low.tif')
    FULL_BASEMAP_PATH_SOUTH = os.path.join(SAFE_FOLDER, 'hrscBasemapSouth/map_south_pole_low.tif')
    
    BACKUP_FOLDER       = os.path.join(SAFE_FOLDER, 'hrscBasemap/output_tile_backups')
    BACKUP_FOLDER_NORTH = os.path.join(SAFE_FOLDER, 'hrscBasemapNorth/output_tile_backups')
    BACKUP_FOLDER_SOUTH = os.path.join(SAFE_FOLDER, 'hrscBasemapSouth/output_tile_backups')
    
    DATABASE_PATH  = os.path.join(SAFE_FOLDER, 'google/googlePlanetary.db')
    RUN_LOG_FOLDER = os.path.join(SAFE_FOLDER, 'hrscMosaicLogs')
    
    BAD_HRSC_FILE_PATH = os.path.join(REPO_FOLDER, 'MassUpload/badHrscSets.csv')
    
    NEW_OUTPUT_TILE_FOLDER = os.path.join(VOLATILE_FOLDER, 'hrscNewOutputTiles')
    HRSC_DOWNLOAD_FOLDER   = os.path.join(VOLATILE_FOLDER, 'hrscDownloadCache')
    HRSC_PROCESSING_FOLDER = os.path.join(VOLATILE_FOLDER, 'hrscProcessedFiles')
    
    if 'byss' in SAFE_FOLDER: # On byss, write where we can post online.
        KML_PYRAMID_FOLDER = '/byss/docroot/smcmich1/hrscMosaicKml'
    else: # Just create somewhere
        KML_PYRAMID_FOLDER = os.path.join(SAFE_FOLDER, 'hrscMosaicKml')
    
    OUTPUT_THUMBNAIL_FOLDER    = os.path.join(SAFE_FOLDER, 'hrscThumbnails')
    OUTPUT_REGISTRATION_FOLDER = os.path.join(SAFE_FOLDER, 'hrscRegistration')
    
    
    # --- Folder notes ---
    # - HRSC_DOWNLOAD_FOLDER holds the downloaded and preprocessed HRSC data
    # - HRSC_PROCESSING_FOLDER holds the fully processed HRSC files
    # - The current crop of tiles is written to NEW_OUTPUT_TILE_FOLDER
    # - The persistent set of final output tiles is kept in BACKUP_FOLDER
    
    
    # TODO: Split this ROI based on an input parameter!
    
    
    # Used to control the area we operate over
    #HRSC_FETCH_ROI = None # Fetch ALL hrsc images
    
    # These are the full ROI's for the north and south pole images
    if options.mapType == MosaicUtilities.PROJ_TYPE_NORTH_POLE:
        HRSC_FETCH_PROJ_ROI = MosaicUtilities.Rectangle(-1959288.087, 1958410.536, -1958437.302, 1959261.357)
        HRSC_FETCH_DEG_ROI  = MosaicUtilities.Rectangle(-180, 180, 60, 90)
    if options.mapType == MosaicUtilities.PROJ_TYPE_SOUTH_POLE:
        HRSC_FETCH_PROJ_ROI = MosaicUtilities.Rectangle(-1959439.740, 1958258.882, -1958280.029, 1959418.630)
        HRSC_FETCH_DEG_ROI  = MosaicUtilities.Rectangle(-180, 180, -90, -60)
    
    if options.nodeIndex > 0:
        # Then this is a specified Google Compute node, and has a designated processing 
        #  region and options.
        # - The polar nodes are as follows: 17-20 (North), 21-24(South) [broken up into quadrants]
        
        def grabQuadrant(rect, quadrant):
            '''Simple function to grab a certain quadrant out of a rectangle.'''
            if quadrant==0: return MosaicUtilities.Rectangle(rect.minX, 0.0, 0.0, rect.maxY) # TL
            if quadrant==1: return MosaicUtilities.Rectangle(0.0, rect.maxX, 0.0, rect.maxY) # TR
            if quadrant==2: return MosaicUtilities.Rectangle(0.0, rect.maxX, rect.minY, 0.0) # BR
            if quadrant==3: return MosaicUtilities.Rectangle(rect.minX, 0.0, rect.minY, 0.0) # BL
            raise Exception('Illegal quadrant!')
        
        def grabLonSlice(rect, quadrant):
            '''Similar function specialized for longitude.'''
            center = rect.getCenterCoord()
            if quadrant==0: return MosaicUtilities.Rectangle(-180, -90, rect.minY, rect.maxY) # TL
            if quadrant==1: return MosaicUtilities.Rectangle(  90, 180, rect.minY, rect.maxY) # TR
            if quadrant==2: return MosaicUtilities.Rectangle(   0,  90, rect.minY, rect.maxY) # BR
            if quadrant==3: return MosaicUtilities.Rectangle( -90,   0, rect.minY, rect.maxY) # BL
            raise Exception('Illegal quadrant!')
        
        if (options.nodeIndex == 17) or (options.nodeIndex == 21):
            HRSC_FETCH_PROJ_ROI = grabQuadrant(HRSC_FETCH_PROJ_ROI, 0)
            HRSC_FETCH_DEG_ROI  = grabLonSlice(HRSC_FETCH_DEG_ROI,  0)
        if (options.nodeIndex == 18) or (options.nodeIndex == 22):
            HRSC_FETCH_PROJ_ROI = grabQuadrant(HRSC_FETCH_PROJ_ROI, 1)
            HRSC_FETCH_DEG_ROI  = grabLonSlice(HRSC_FETCH_DEG_ROI,  1)
        if (options.nodeIndex == 19) or (options.nodeIndex == 23):
            HRSC_FETCH_PROJ_ROI = grabQuadrant(HRSC_FETCH_PROJ_ROI, 2)
            HRSC_FETCH_DEG_ROI  = grabLonSlice(HRSC_FETCH_DEG_ROI,  2)
        if (options.nodeIndex == 20) or (options.nodeIndex == 24):
            HRSC_FETCH_PROJ_ROI = grabQuadrant(HRSC_FETCH_PROJ_ROI, 3)
            HRSC_FETCH_DEG_ROI  = grabLonSlice(HRSC_FETCH_DEG_ROI,  3)
    
        options.bucketPrefix = 'node_' + str(options.nodeIndex)
        options.uploadBucket = 'hrsc_map_storage'
    
    print 'Using fetch ROI ' + str(HRSC_FETCH_DEG_ROI)
    print HRSC_FETCH_PROJ_ROI
    
    # Also set up logging here.
    # - Log tiles are timestamped as is each line in the log file
    currentTime = datetime.datetime.now()
    logPath = os.path.join(RUN_LOG_FOLDER, ('hrscMosaicLog_%s.txt' % currentTime.isoformat()) )
    logging.basicConfig(filename=logPath,
                        format=MosaicUtilities.LOG_FORMAT_STR,
                        level=logging.DEBUG)
    
    
    return options

def main(argsIn):
    '''Outer try-catch handling to make sure the program does not silently die.
       Also parses input options.'''
    sentMessage = False
    try:

        # Parse the input arguments and set global params
        # - Using global params to avoid passing a giant list of paths around.
        options = setGlobalConfigs(argsIn)

        sentMessage = mainProcessingFunction(options)
        text = 'Outer processing function exited without exception.'
    except Exception, e:
        # When we fail to fetch a data set, send out a failure message and keep going.
        text = ('Caught fatal exception in outer HRSC mosaic handler!\n' + 
                     str(e) + '\n' + str(sys.exc_info()[0]) + '\n')
        print(traceback.format_exc())

    # Send out an email if we did not already do so.
    if not sentMessage:
        MosaicUtilities.sendEmail('scott.t.mcmichael@nasa.gov', 
                                  'HRSC map program stopped.',
                                  text)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))




