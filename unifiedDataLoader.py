#!/usr/bin/env python
# __BEGIN_LICENSE__
#  Copyright (c) 2009-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NGT platform is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance with the
#  License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__

import sys

#from BeautifulSoup import BeautifulSoup
from bs4 import BeautifulSoup

import os, glob, optparse, re, shutil, subprocess, string, time, urllib, urllib2

import multiprocessing, traceback

#import sqlite3
from pysqlite2 import dbapi2 as sqlite3

import mapsEngineUpload, IrgStringFunctions, IrgGeoFunctions

#------------------------------------------------------------------
# Global definitions

SENSOR_TYPE_HiRISE = 0
SENSOR_TYPE_HRSC   = 1
SENSOR_TYPE_CTX    = 2
SENSOR_TYPE_THEMIS = 3

# Main status codes describing upload status
STATUS_NONE      = 0
STATUS_UPLOADED  = 1
STATUS_CONFIRMED = 2
STATUS_ERROR     = -1

# Extra status codes, equal to STATUS_NONE plus assigned to a specific machine
STATUS_ASSIGNED_LUNOKHOD1 = 1000
STATUS_ASSIGNED_LUNOKHOD2 = 1001
STATUS_ASSIGNED_M         = 1002
STATUS_ASSIGNED_BYSS      = 1003
STATUS_ASSIGNED_ALDERAAN  = 1004
STATUS_ASSIGNED_CHIP      = 1005
STATUS_ASSIGNED_DALE      = 1006

SENSOR_CODES = {'hirise' : SENSOR_TYPE_HiRISE,
                'hrsc'  : SENSOR_TYPE_HRSC,
                'ctx'   : SENSOR_TYPE_CTX,
                'THEMIS': SENSOR_TYPE_THEMIS}

# Set this code equal to the machine's STATUS code to only process those files!
#THIS_MACHINE_CODE = STATUS_ASSIGNED_LUNOKHOD2

#----------------------------------------------------------------

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, ''' Script for grabbing and uploading data files'''

    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class TableRecord:
    '''Helper class for parsing table records'''
    
    # Variables
    data = None
    
    def __init__(self, row):
        self.data = row

    def tableId(self): # Unique ID
        return self.data[0]
    def sensor(self): # Sensor code
        return int(self.data[1])
    def subtype(self): # String identifying the sub-sensor
        return self.data[2]
    def setName(self): # Data set name
        return self.data[3]
    def acqTime(self): # Time image taken
        return self.data[4]
    def status(self): # Upload status
        return self.data[5]
    def version(self): # Version tag
        return self.data[6]
    def remoteURL(self): # Original file source
        return self.data[7]
    def assetID(self): # Asset ID assigned by Maps Engine
        return self.data[8]
    def uploadTime(self): # Time file uploaded to Maps Engine
        return self.data[9]
    #def minLat(self): # TODO: Nice Bounding Box wrapper
    #    return self.data[10]

# Function that is passed in to findAllDataSets below.
def addDataRecord(db, sensor, subType, setName, remoteURL):
    '''Adds a sensor data entry in the database'''

    # Do nothing if this dataset is already in the database
    cursor = db.cursor()
    cursor.execute("SELECT * FROM Files WHERE sensor=? AND subtype=? AND setname=?",
                      (str(sensor), subType, setName))
    if cursor.fetchone() != None:
        return True

    print 'Adding ' + setName + ',  ' + subType
    
    cursor.execute("INSERT INTO Files VALUES(null, ?, ?, ?, null, 0, null, ?, null, null, null, null, null, null)",
               (str(sensor), subType, setName, remoteURL))
    db.commit()
    return True
    
#--------------------------------------------------------------------------------
# List of functions that need to be provided for each of the sensors!


#def getCreationTime(fileList):
#    """Extract the file creation time and return in YYYY-MM-DDTHH:MM:SSZ format"""
#    Takes the list of files returned by fetchAndPrepFile
#    return null

#def getBoundingBox(fileList):
#    """Return the bounding box for this data set in the format (minLon, maxLon, minLat, maxLat)"""
#    Takes the list of files returned by fetchAndPrepFile
#    return null

#def findAllDataSets(db, dataAddFunctionCall, sensorCode):
#    '''Add all known data sets to the SQL database'''
#    return False

#def fetchAndPrepFile(setName, subtype, remoteURL, workDir):
#    '''Retrieves a remote file and prepares it for upload'''
#    returns a list of created files with the first one being the one to upload

#def getUploadList(fileList):
#    '''Returns the subset of the fileList that needs to be uploaded to GME'''
#    return fileList[0]

#--------------------------------------------------------------------------------
# Other functions

# TODO: Move to a helper file!
def getCurrentTimeString():
    """Return the current time in YYYY-MM-DDTHH:MM:SSZ format"""    
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

def logWriter(logQueue, logPath):
    '''Listens for messages on the queue and writes them to a file.'''

    print 'Starting log writer'
    f = open(logPath, 'a') 
    while 1: # Run until the process is killed
        message = logQueue.get() # Wait for a new message
        if message == 'stop_queue': # Check for the quit signal
            break
        f.write(str(message))
        f.flush()
    f.close()
    print 'Writer stopped'


def uploadFile(dbPath, fileInfo, logQueue, workDir):
    """Uploads a remote file to maps engine"""
    
    print 'Uploading file ' + fileInfo.remoteURL()

    try:

        db = sqlite3.connect(dbPath)
        cursor = db.cursor()
        
        # Choose the "library" to use based on the sensor type
        # - The current version of each sensor's data files is set here.
        if   fileInfo.sensor() == SENSOR_TYPE_HiRISE:
            from hiriseDataLoader import fetchAndPrepFile, getCreationTime, getBoundingBox, getUploadList
            version = 1
        elif fileInfo.sensor() == SENSOR_TYPE_HRSC:
            from hrscDataLoader import fetchAndPrepFile, getCreationTime, getBoundingBox, getUploadList
            version = 1
        elif fileInfo.sensor() == SENSOR_TYPE_CTX:
            from ctxDataLoader import fetchAndPrepFile, getCreationTime, getBoundingBox, getUploadList
            version = 1
        else:
            raise Exception('Sensor type ' + fileInfo.sensor() + ' is not supported!')
        
        # Call sensor-specific function to fetch and prepare the file.

        localFileList = fetchAndPrepFile(fileInfo.setName(), fileInfo.subtype(), fileInfo.remoteURL(), workDir)
        if len(localFileList) == 0:
            raise Exception('Failed to retrieve any local files!')

        # In order for this to work, all modules must set the subtype consistently.
        isDem = (fileInfo.subtype() == 'DEM')

        preppedFilePath = localFileList[0]
        print preppedFilePath
        
        if not os.path.exists(preppedFilePath):
            raise Exception('Prepped file does not exist: ' + preppedFilePath)

        # Check the file size
        MAX_UPLOAD_SIZE = 1024*1024*1024*10 # 10 GB upload limit
        statinfo = os.stat(preppedFilePath)
        fileSizeBytes = statinfo.st_size
        if (fileSizeBytes > MAX_UPLOAD_SIZE):
            raise Exception('Processed file is too large to upload!') # TODO: Automated split of file?
        

        # Need the time string before uploading
        timeString = getCreationTime(localFileList)

        # Go ahead and get the bounding box       
        fileBbox = getBoundingBox(localFileList)

        # Upload all of the files
        uploadList = getUploadList(localFileList)
        cmdArgs = ['--sensor', str(fileInfo.sensor()), '--acqTime', timeString]
        if isDem: # For DEM files, add an additional DEM tag.
            cmdArgs.append('--tag') 
            cmdArgs.append('DEM')
        for f in uploadList:
            cmdArgs.append(f)
        #print cmdArgs
        assetId = mapsEngineUpload.main(cmdArgs)
        #assetId = 12345 # DEBUG!

        # The file won't make it up every time so the --checkUploads call will be
        #   needed to catch stragglers.

        # Extract a timestamp from the file
        cursor.execute("UPDATE Files SET acqTime=? WHERE idx=?", (timeString, str(fileInfo.tableId())))
        
        # Log the bounding box
        bboxString = ('Bbox: ' + str(fileBbox[0]) +' '+ str(fileBbox[1]) +' '+ str(fileBbox[2]) +' '+ str(fileBbox[3]))
        cursor.execute("UPDATE Files SET minLon=?, maxLon=?, minLat=?, maxLat=? WHERE idx=?",
                       (str(fileBbox[0]), str(fileBbox[1]), str(fileBbox[2]), str(fileBbox[3]), str(fileInfo.tableId())))

        # Log the time string
        currentTimeString = getCurrentTimeString()
        cursor.execute("UPDATE Files SET status=?, uploadTime=?, version=?, assetID=? WHERE idx=?",
                        (str(STATUS_UPLOADED), currentTimeString, str(version), str(assetId), str(fileInfo.tableId())))

        db.commit()
        db.close() 

        # Record that we uploaded the file
        logString = fileInfo.setName() +', '+ str(assetId) +', '+ bboxString + '\n' # Log path and the Maps Engine asset ID
        #print logString
        logQueue.put(logString)

            
        # Delete all the local files left by the prep function
        for f in localFileList:
            os.remove(f)
    
    except Exception, e:# Make sure an error does not stop other uploads
    
        print traceback.format_exc()
        print str(e)

        # Record the error in the database
        cursor.execute("UPDATE Files SET status=? WHERE idx=?", 
                       (str(STATUS_ERROR), str(fileInfo.tableId())))
        db.commit()
        db.close()

        try: # Delete localFileList if it exists
            for f in localFileList:
                os.remove(f)
                #pass
        except:
            pass

        return -1

    print 'Finished uploading data file!'
    return assetId

# TODO: Write a wrapper script to keep this going
def uploadNextFile(dbPath, sensorCode, outputFolder, numFiles=1, numThreads=1):
    """Determines the next file to upload, uploads it, and logs it"""
    
    print 'Searching for next file to upload...'

    db = sqlite3.connect(dbPath)
    cursor = db.cursor()
        
    # Query the SQL database for one or more entries for this sensor which have not been uploaded yet.
    # - Using THIS_MACHINE_CODE means that only files allocated to this machine will be handled!
    cursor.execute('SELECT * FROM Files WHERE sensor=? AND status=? LIMIT ?',
                   #(str(sensorCode), str(THIS_MACHINE_CODE), str(numFiles)))
                   (str(sensorCode), str(STATUS_NONE), str(numFiles)))
    rows = cursor.fetchall()
    db.close()
    if rows == []: # Make sure we found the next lines
        raise Exception('Could not find any data files left to upload!')
    
    # Create processing pool
    # Limit number of threads to numFiles+1
    if numThreads > numFiles+1:
        numThreads = numFiles+1
    print 'Spawning ' + str(numThreads) + ' worker threads'
    pool = multiprocessing.Pool(processes=numThreads+1) # One extra thread for the logWriter

    # Create multiprocessing manager and a queue
    manager = multiprocessing.Manager()
    queue   = manager.Queue()    

    # Start up the log writing thread
    logPath = os.path.join(outputFolder, 'activitylog.txt')
    print 'Writing output to file ' + logPath
    logResult = pool.apply_async(logWriter, args=(queue, logPath))
    
    # For each item we fetched, spawn a process to upload it.
    jobResults = []
    for line in rows:
        fileInfo = TableRecord(line) # Wrap the data line
        jobResults.append(pool.apply_async(uploadFile, args=(dbPath, fileInfo, queue, outputFolder)))
    
    
    # Wait until all threads have finished
    print 'Waiting for all threads to complete...'
    for r in jobResults:
        r.get()
    
    # Stop the queue and all the threads
    print 'Cleaning up...'
    queue.put('stop_queue')
    pool.close()
    pool.join()
    
    print 'All threads finished!'
   
    return True
    
    
def checkUploads(db, sensorType):

    print 'Checking the status of uploaded files...'    

    cursor = db.cursor()

    # Get server authorization and hold on to the token
    bearerToken = mapsEngineUpload.authorize()

    MAX_NUM_RETRIES     = 2  # Max number of times to retry (in case server is busy)
    SLEEP_TIME          = 1.5 # Time to wait between retries (Google handles only one operation/second)
    FAILURE_RESET_LIMIT = 20

    # Query the data base for all sensor data files which have been uploaded but not confirmed 
    cursor.execute('SELECT * FROM Files WHERE sensor=? AND status=?', (str(sensorType), str(STATUS_UPLOADED)))
    rows = cursor.fetchall()
    print 'Found ' + str(len(rows)) + ' entries'
    skip = 0
    failureCount = 0
    for row in rows:
        if skip > 0:
            skip = skip - 1
            continue

        time.sleep(1.1) # Minimum sleep time

        # Wrap the row to make it easier to get information
        o = TableRecord(row)
        
        # Check if this asset was uploaded
        # TODO: Option to force checking for confirmed files
        print 'Checking asset ID = ' + o.assetID()
        for i in range(1,MAX_NUM_RETRIES):
            status, responseCode = mapsEngineUpload.checkIfFileIsLoaded(bearerToken, o.assetID())
            if (responseCode == 403) or (responseCode == 503) or (responseCode == 401):
                print 'Server is busy, sleeping ' + str(SLEEP_TIME) + ' seconds...'
                failureCount += 1 # Keep track of successive failures
                time.sleep(SLEEP_TIME) 
            else: # Got a valid response from Google
                failureCount = 0 # Reset the failure count
                break
        gotValidResponse = (failureCount == 0) # Record if we got a real answer
        
        if failureCount >= FAILURE_RESET_LIMIT: # Too many failures in a row, reset our connection!
            print 'Hit successive failure count limit, resetting our connection!!!!!!!!!!!!!!!!!!!!'
            time.sleep(10)
            bearerToken = mapsEngineUpload.authorize()
            time.sleep(10)
            failureCount = 0
            
        if status:
            # Mark the file upload as confirmed so we don't check it again
            cursor.execute("UPDATE Files SET status=? WHERE idx=?",
                           (str(STATUS_CONFIRMED), str(o.tableId())))
            db.commit()
        else: # Did not get confirmation of successful upload
            if gotValidResponse: # Did get confirmation of failed upload
                print 'Data set ' + o.setName() + ' was not uploaded correctly!'
                # Update the file info in the database to show it was not updated correctly.
                # - TODO: Is there a way to make sure we re-upload to the same asset ID?
                cursor.execute("UPDATE Files SET status=? WHERE idx=?",
                               (str(STATUS_NONE), str(o.tableId())))
                db.commit()
            else:
                print 'Never got a valid response, failure count = ' + str(failureCount)

    print 'Finished checking uploaded files.'
    

def updateDbFromWeb(db, sensorType):
    '''Update the database from files already uploaded'''
    # This is required if the database file is lost after data has been uploaded!

    print 'Querying Maps Engine for uploaded file list...'

    cursor = db.cursor()

    MAX_NUM_RETRIES = 4  # Max number of times to retry (in case server is busy)
    SLEEP_TIME      = 1.1 # Time to wait between retries (Google handles only one operation/second)

    # Get server authorization and hold on to the token.
    bearerToken = mapsEngineUpload.authorize()

    # Get the list of assets
    assetList = mapsEngineUpload.getRasterList(bearerToken)
    if not assetList:
        raise Exception('Failed to detect any uploaded assets!')
    print 'Found ' + str(len(assetList)) + ' existing raster assets in Maps Engine'

    # Add each asset to the database if it is not already there.
    for asset in assetList:

        uploadName = asset['name']
        if uploadName[-4:] == '.jp2': # These are based on ASU's images
            print 'Skipping uploaded jp2 file'
            continue

        # Fetch the the acquisition time
        for i in range(1,MAX_NUM_RETRIES):
            status, detailedInfo = mapsEngineUpload.queryUploadedFile(bearerToken, asset['assetID'])
            if (not status) and (detailedInfo == 403) or (detailedInfo == 503):
                print 'Server is busy, sleeping ' + str(SLEEP_TIME) + ' seconds...'
                time.sleep(SLEEP_TIME)
            else:
                break

        # CTX
        sensor     = SENSOR_TYPE_CTX
        version    = 1 
        subType    = '' # This is the volume which we can't get from online
        setName    = uploadName[4:-4] # Extract this from the uploaded asset name
        remoteURL  = '' # Not stored online!
        acqTime    = detailedInfo['acquisitionTime']['start'] #TODO: Check syntax!
        uploadTime = asset['uploadTime']

        print 'Adding ' + setName + ',  ' + asset['assetID']

        # Check if the asset is already loaded in to the database.
        cursor.execute('SELECT * FROM Files WHERE assetId=?', [asset['assetID']])
        rows = cursor.fetchall()
        if len(rows) >= 1:            
            record = TableRecord(rows[0]) # Wrap the row to make it easier to get information

            # Check if this table entry has already been updated
            if record.status() == STATUS_UPLOADED:
                print 'This entry has already been updated!'
                #continue

            # Update the table row with the upload information
            cursor.execute("UPDATE Files SET acqTime=?, status=?, version=?, assetId=?, uploadTime=?, minLon=?, maxLon=?, minLat=?, maxLat=? WHERE idx=?",
                           (acqTime, str(STATUS_UPLOADED), version, asset['assetID'], uploadTime, 
                            asset['minLat'], asset['minLon'], asset['maxLat'], asset['maxLon'], str(record.tableId())))
        else:
            # We haven't seen the asset ID before, but what about the set name?
            cursor.execute('SELECT * FROM Files WHERE setname=?', [setName])
            rows = cursor.fetchall()

            if len(rows) > 1: # Check for this error case
                raise Exception('Multiple rows found with set name ' + setname)

            if len(rows) < 1: # Set name not in the database
                # Otherwise we need to add the asset to the database.
                # - Maybe this should never happen?
                raise Exception('WARNING: This data set not found in the database!')
                cursor.execute("INSERT INTO Files VALUES(null, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                               (str(sensor), subType, setName, acqTime, str(STATUS_UPLOADED), version, remoteURL, asset['assetID'], 
                                uploadTime, asset['minLat'], asset['minLon'], asset['maxLat'], asset['maxLon']))

            else: # We found the set in the database, update it.
                record = TableRecord(rows[0]) # Wrap the row to make it easier to get information

                cursor.execute("UPDATE Files SET acqTime=?, status=?, version=?, assetId=?, uploadTime=?, minLon=?, maxLon=?, minLat=?, maxLat=? WHERE idx=?",
                               (acqTime, str(STATUS_UPLOADED), version, asset['assetID'], uploadTime, 
                                asset['minLat'], asset['minLon'], asset['maxLat'], asset['maxLon'], str(record.tableId())))

        # Execute whatever action we decidede on.
        db.commit()

        #raise Exception('DEBUG FORMATS!')

        # Note that the upload status is not confirmed.
        # - The upload checker script should be run to take care of this.


    print 'Finished checking uploaded files.'


def getDataList(db, sensorCode):
    '''Update the list of available data sets from the given sensor'''

    # Choose the "library" to use based on the sensor type
    # - The current version of each sensor's data files is set here.
    if   sensorCode == SENSOR_TYPE_HiRISE:
        from hiriseDataLoader import findAllDataSets
    elif sensorCode == SENSOR_TYPE_HRSC:
        from hrscDataLoader import findAllDataSets
    elif sensorCode == SENSOR_TYPE_CTX:
        from ctxDataLoader import findAllDataSets
    else:
        raise Exception('Sensor type ' + sensorCode + ' is not supported!')
    
    return findAllDataSets(db, addDataRecord, sensorCode)


def checkForBadUploads(sensorCode, db = None):
    '''Search the local SQL database to find bad uploads that the DB has missed.'''

    # TODO: Parse user inputs to get these!
    tag         = 'HiRISE'
    cacheFolder = '/home/pirl/smcmich1/tempCache2' # Web API query results are backed up in to this folder as JSON files.
    
    print 'Retrieving asset list for sensor ' + tag + ' in folder ' + cacheFolder
    if not os.path.exists(cacheFolder):
        os.mkdir(cacheFolder)

    # Get server authorization and hold on to the token.
    bearerToken = mapsEngineUpload.authorize()

    # First get a list of all assets.
    assetList = mapsEngineUpload.findAllRasterUploads(bearerToken, cacheFolder, tag)
    if not assetList:
        raise Exception('Failed to detect any uploaded assets!')
    print 'Found ' + str(len(assetList)) + ' existing raster assets in Maps Engine'


    if not db: # If the DB is not present just look up the files.
        return

    print 'Looking for uploads not contained in the database...'

    cursor = db.cursor()

    # Now make a list of all assets which are not found in the database
    extraUploads = []
    for upload in assetList:
        thisAssetId = str(upload['assetID']).strip()
        
        # Query the SQL database to find this upload
        cursor.execute('SELECT * FROM Files WHERE sensor=? and assetID=?', (str(sensorCode), thisAssetId))
        rows = cursor.fetchall()
        if rows == []: # An orphan upload, what we are looking for!
            #print 'WARNING: Zero DB rows found for asset ID ' + thisAssetId
            extraUploads.append(upload)
        if len(rows) > 1: # Duplicate DB entries!  This is unexpected.
            print 'WARNING: More than one DB row found for asset ID ' + thisAssetId
            for row in rows:
                print row

    db.close() # Done with the database

    print 'Found ' + str(len(extraUploads)) + ' unmatched raster assets!\n\n'

    # Print out each extra upload name and asset id, one per line.
    for u in extraUploads:
        print u['name'] +', '+ u['assetID']


    # TODO: Do something about the bad uploads!

    #if assetInfo['processingStatus'] != 'complete':
    #    print '\n===========================\n'
    #    print f
    #    #badAssetList.append(assetInfo)

    #print 'Finished fetching results, now searching for duplicates...'
    #   
    ## Loop through all the dictionaries and find duplicate entries!
    #numAssets = len(assetList)
    #for i in range(0,numAssets): # For each asset
    #    entry = assetList[i]
    #    matchCount = 0
    #    for j in range(i+1,numAssets): # See if any remaining assets match
    #        other = assetList[j]
    #        if (entry['name'] == other['name']) and (entry['assetID'] != other['assetID']):
    #            matchCount += 1
    #    badAssetList.append(entry) # This returns the first duplicate of each pair
    #   
    #print 'Finished searching for duplicates!'
    #   
    #return badAssetList

    #print assetList

#--------------------------------------------------------------------------------

def main():

    print "Started unifiedDataLoader.py"

    # ----- Parse input arguments -----
    usage = "usage: unifiedDataLoader.py <sensor name> [--help][--manual]\n  "
    parser = optparse.OptionParser(usage=usage)
    
    parser.add_option("-u", "--upload", dest="upload", type=int,
                      help="Upload this many files instead of fetching the list.")

    parser.add_option("--checkUploads", action="store_true", default=False,
                                dest="checkUploads",  help="Check that all uploaded files actually made it up.")

    parser.add_option('--verifyWebUploads', action='store_true', dest='verifyWebUploads', default=False,
                      help='Check the web API for bad uploads.')

    parser.add_option("--threads", type="int", dest="numThreads", default=1,
                      help="Number of threads to use.")
                      
    parser.add_option("--manual", action="callback", callback=man,
                      help="Read the manual.")
    (options, args) = parser.parse_args()

    # Make sure the user passed in the name of the sensor
    try:
        options.sensorType = SENSOR_CODES[args[0].lower()]
    except:
        raise Exception('Did not recognize sensor name: ' + args[0])

    ## Now check for the output folder (working directory)
    #if len(args) != 2:
    #    raise Exception('Missing output folder!')
    #    return 1;
    #options.outputFolder = args[1]
    # The output path is hardcoded for now with subfolders for each sensor
    options.outputFolder = os.path.join('/home/pirl/smcmich1/Data/google/', args[0].lower())
    # -- Done parsing input arguments --

    # Check the database connection
    # - Default should be to db = a thread-safe connection
    # - TODO: Find this database without hard coding it!
    dbPath = '/home/pirl/smcmich1/Data/google/googlePlanetary.db'
    db = sqlite3.connect(dbPath)
    print 'Connected to database'
    
    print "Beginning processing....."
    
    startTime = time.time()
    
    # Make sure the working directory exists
    if not os.path.exists(options.outputFolder):
        os.mkdir(options.outputFolder)

    if options.verifyWebUploads:
        # Rarely used option to search for online only problem files!
        # - This does not touch the database, only Maps Engine online.
        checkForBadUploads(options.sensorType, db)
        return 0

    ## Rarely used option to update after a local database loss!    
    #updateDbFromWeb(db, options.sensorType)
    #return 0
    
    if options.checkUploads: # Check to see if uploaded files made it up ok
        checkUploads(db, options.sensorType)
        
    elif options.upload: # Upload one or more files
        db.close()
        uploadNextFile(dbPath, options.sensorType, options.outputFolder, options.upload, options.numThreads)
        
    else: # Update the database of data files
        getDataList(db, options.sensorType)
    
    
    endTime = time.time()
    
    print "Finished in " + str(endTime - startTime) + " seconds."
    return 0


if __name__ == "__main__":
    sys.exit(main())
