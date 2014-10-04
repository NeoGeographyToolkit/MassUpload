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

from BeautifulSoup import BeautifulSoup

import os, glob, optparse, re, shutil, subprocess, string, time, urllib, urllib2

import multiprocessing

import mapsEngineUpload, IrgStringFunctions, IrgGeoFunctions


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, ''' Script for grabbing HRSC data files'''

    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


#--------------------------------------------------------------------------------


def getCreationTime(fileList):
    """Extract the file creation time and return in YYYY-MM-DDTHH:MM:SSZ format"""
    
    if len(fileList) < 2:
        raise Exception('Error, missing original file path!')
    filePath = fileList[1]
    
    # Use subprocess to parse the command output
    cmd = ['gdalinfo', filePath]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    cmdOut, err = p.communicate()

    # Find the time string in the text
    timeString = IrgStringFunctions.getLineAfterText(cmdOut, 'PRODUCT_CREATION_TIME=')
    
    # Get to the correct format
    timeString = timeString.strip()
    timeString = timeString[:-5] + 'Z'
    
    return timeString
    
def getBoundingBox(fileList):
    """Return the bounding box for this data set in the format (minLon, maxLon, minLat, maxLat)"""
    return IrgGeoFunctions.getImageBoundingBox(fileList[0])

def findAllDataSets(db, dataAddFunctionCall, sensorCode):
    '''Add all known data sets to the SQL database'''
    
    print 'Updating HRSC PDS data list...'
          
    baseUrl = "http://pds-geosciences.wustl.edu/mex/mex-m-hrsc-5-refdr-mapprojected-v2/mexhrsc_1001/data/"

    # Parse the top PDS level
    print baseUrl
    parsedIndexPage = BeautifulSoup(urllib2.urlopen((baseUrl)).read())
    
    # Loop through outermost directory
    for line in parsedIndexPage.findAll('a'):
        
        # Skip links we are not interested in
        if len(line.string.strip()) != 4:
            continue
        
        dataPrefix = 'h' + line.string
        
        subFolderUrl = baseUrl + line.string + '/'
        print subFolderUrl
        parsedDataPage = BeautifulSoup(urllib2.urlopen((subFolderUrl)).read())
        
        # Loop through the data files
        # - There is a core set of files on each page but there can be
        #   more with incremented image numbers.
        for d in parsedDataPage.findAll('a'):
            if '[' in d.string: # Skip some bad links
                continue

            print d.string

            dataFileName = d.string[:-4]     # Lop off the .img portion
            subtype      = dataFileName[-3:] # Extract the type

            # Easier to just call this function to get the URL
            url = generatePdsPath(dataFileName)

            dataAddFunctionCall(db, sensorCode, subtype,   dataFileName, url)

    
def fetchAndPrepFile(setName, subtype, remoteURL, workDir):
    '''Retrieves a remote file and prepares it for upload'''
    
    #print 'Uploading file ' + remoteURL
        
    localFileName = os.path.splitext(os.path.basename(remoteURL))[0]+'.tif'
    localFilePath = os.path.join(workDir, localFileName)
    downloadPath  = os.path.join(workDir, os.path.basename(remoteURL))
    
    if not os.path.exists(localFilePath):
        # Download the file
        cmd = 'wget ' + remoteURL + ' -O ' + downloadPath
        print cmd
        os.system(cmd)
    
    # Convert to GTiff format
    cmd = 'gdal_translate -of GTiff ' + downloadPath + ' ' + localFilePath
    print cmd
    os.system(cmd)
    
    # Insert a delay here to make sure that the tiny files do not get uploaded too fast for Google.
    # - For the larger files this will barely be noticable.
    time.sleep(5)

    return [localFilePath, downloadPath]
    

#--------------------------------------------------------------------------------





# fileType is the file name after the prefix
def generatePdsPath(filePrefix):
    """Generate the full PDS path for a given HRSC data file"""
       
    # File prefix looks like this: hHXXX_DDDD_SSS
    fileType = '.img'
    
    # Extract the run number --> HXXX
    runNum = filePrefix[1:5]
       
    filename = filePrefix + fileType
    baseUrl  = "http://pds-geosciences.wustl.edu/mex/mex-m-hrsc-5-refdr-mapprojected-v2/mexhrsc_1001/data/"
    fullUrl  = baseUrl + runNum +"/"+ filename

    #print filePrefix + fileType + ' -> ' + fullUrl
    return fullUrl


fileBbox = IrgGeoFunctions.getImageBoundingBox(preppedFilePath)

def getDataList(outputFilePath):
    """Populates a text file list of all available HRSC RDR data"""
       
    print 'Updating HRSC PDS data list...'
          
    baseUrl = "http://pds-geosciences.wustl.edu/mex/mex-m-hrsc-5-refdr-mapprojected-v2/mexhrsc_1001/data/"

    # Parse the top PDS level
    parsedIndexPage = BeautifulSoup(urllib2.urlopen((baseUrl)).read())
    
    outputFile = open(outputFilePath, 'w')
    
    # Loop through outermost directory
    for line in parsedIndexPage.findAll('a'):
        
        dataPrefix = 'h' + line.string
        
        subFolderUrl = baseUrl + line.string + '/'
        parsedDataPage = BeautifulSoup(urllib2.urlopen((subFolderUrl)).read())
        
        # Loop through the data files
        # - There is a core set of files on each page but there can be
        #   more with incremented image numbers.
        for d in parsedDataPage.findAll('a'):
            dataFileName = d.string[:-4] # Lop off the .img portion
        
            outputFile.write(dataFileName + '\n')

    outputFile.close()
    
    print 'Wrote updated data list to ' + outputFilePath


def uploadFile(filePrefix, remoteFilePath, logQueue, tempDir):
    """Uploads a remote file to maps engine"""
    
    print 'Uploading file ' + remoteFilePath
        
    localFileName = os.path.splitext(os.path.basename(remoteFilePath))[0]+'.tif'
    localFilePath = os.path.join(tempDir, localFileName)
    downloadPath  = os.path.join(tempDir, os.path.basename(remoteFilePath))
    
    if not os.path.exists(localFilePath):
        # Download the file
        cmd = 'wget ' + remoteFilePath + ' -O ' + downloadPath
        print cmd
        os.system(cmd)
    
    # Extract the image creation time
    timeString = getCreationTime(downloadPath)
    
    # Convert to GTiff format
    cmd = 'gdal_translate -of GTiff ' + downloadPath + ' ' + localFilePath
    print cmd
    os.system(cmd)
    
    # Upload the file
    cmdArgs = [localFilePath, '--sensor', '1', '--acqTime', timeString]
    #print cmdArgs
    assetId = mapsEngineUpload.main(cmdArgs)
    #assetId = 12345
    
    #TODO: Check to make sure the file made it up!
    
    # Find out the bounding box of the file and generate a log string
    fileBbox = IrgGeoFunctions.getImageBoundingBox(localFilePath)
    bboxString = ('Bbox: ' + str(fileBbox[0]) +' '+ str(fileBbox[1]) +' '+ str(fileBbox[2]) +' '+ str(fileBbox[3]))
    
    # Record that we uploaded the file
    logString = filePrefix +', '+ str(assetId) +', '+ bboxString + '\n' # Log path and the Maps Engine asset ID
    #print logString
    logQueue.put(logString)

        
    # Delete the file
    print 'rm ' + localFilePath
    #os.remove(localFilePath)
    
    print 'Finished uploading HRSC data file'
    return assetId


def logWriter(logQueue, logPath):
    '''listens for messages on the q, writes to file. '''

    print 'Starting writer'
    f = open(logPath, 'a') 
    while 1: # Run until the process is killed
        message = logQueue.get() # Wait for a new message
        if message == 'stop_queue': # Check for the quit signal
            break
        f.write(str(message))
        f.flush()
    f.close()
    print 'Writer stopped'


# TODO: Select from different file types?
def uploadNextFile(dataListPath, outputFolder, numFiles=1, numThreads=1):
    """Determines the next file to upload, uploads it, and logs it"""
    
    print 'Searching for next file to upload...'
    
    # Set up the output paths    
    logPath = os.path.join(outputFolder, 'uploadedPatchFilesTest.txt')
    
    inFile = open(dataListPath,    'r')

    # Get the last line read (could be more efficient)
    lastUploadedLine = ''
    if os.path.exists(logPath):
        outFile = open(logPath, 'r')
        for line in outFile:
            lastUploadedLine = line
        outFile.close()
    lastUploadedLine = lastUploadedLine.split(',')[0].strip() # Remove the asset ID from the string
    #print '#' + lastUploadedLine + '#'
    
    # Now find that line in the input file list
    linesToProcess = []
    breakNext = -1
    for line in inFile:
        line = line.strip() # Remove \n from the line
        #print '#' + line + '#'
        
        if breakNext > 0:
            linesToProcess.append(line) # Record this line and move on to the next
            breakNext = breakNext - 1
            continue
        if breakNext == 0: # No more lines to get!
            break

        if (lastUploadedLine == ''): # This is the first file!
            linesToProcess.append(line)
            breakNext = numFiles - 1 # Still get the next N-1 files
            continue
        if (line == lastUploadedLine): # Found the last one downloaded
            breakNext = numFiles # Get the next N-1 files
            continue
    
    if linesToProcess == []: # Make sure we found the next lines
        raise Exception('Failed to find line that comes after: ' + lastUploadedLine)
    
    # Create processing pool
    print 'Spawning ' + str(numThreads) + ' worker threads'
    pool = multiprocessing.Pool(processes=numThreads+1) # One extra thread for the logWriter

    # Create multiprocessing manager and a queue
    manager = multiprocessing.Manager()
    queue   = manager.Queue()    

    # Start up the log writing thread
    print 'Writing output to file ' + logPath
    logResult = pool.apply_async(logWriter, args=(queue, logPath))
    
    # For each line generate the full download path
    jobResults = []
    for line in linesToProcess:
        fullPath = generatePdsPath(line.lower())
        jobResults.append(pool.apply_async(uploadFile, args=(line, fullPath, queue, outputFolder)))
    
    
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
    
def checkUploads(logPath):

    print 'Checking the status of uploaded files...'    

    # Get server authorization and hold on to the token
    bearerToken = mapsEngineUpload.authorize()

    if not os.path.exists(logPath):
        raise Exception('Input log file ' + logPath + ' does not exist!')
        
    outFile = open(logPath, 'r')
    for line in outFile:
        # Extract the asset ID
        prefix, assetId, bbox = line.split(',')
        
        # Check if this asset was uploaded
        print 'Checking asset ID = ' + assetId.strip()
        status, responseCode = mapsEngineUpload.checkIfFileIsLoaded(bearerToken, assetId.strip())
        
        if not status:
            print 'Prefix ' + prefix + ' was not uploaded correctly!'
            raise Exception('oaeunteosu')
            # TODO: Do something about it!
        
    outFile.close()
    


#--------------------------------------------------------------------------------

def main():

    print "Started hrscDataLoader.py"


    try:
        try:
            usage = "usage: hrscDataLoader.py <output folder> [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            
            parser.add_option("-u", "--upload", dest="upload", type=int,
                              help="Upload this many files instead of fetching the list.")

            parser.add_option("--checkUploads", action="store_true", default=False,
                                        dest="checkUploads",  help="Verify that all uploaded files actually made it up.")

            parser.add_option("--threads", type="int", dest="numThreads", default=1,
                              help="Number of threads to use.")
                              
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            (options, args) = parser.parse_args()

            if len(args) < 1:
                raise Exception('Missing output folder!')
                return 1;
            options.outputFolder = args[0]


        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        if not os.path.exists(options.outputFolder):
            os.mkdir(options.outputFolder)

        # These are input lists that just show the available images
        #pspListPath = os.path.join(options.outputFolder, 'pspList.csv')

        # TODO: Set this
        inputListPath = 'hrscImageList_smallPatch.csv'

        if options.checkUploads:
            checkUploads(os.path.join(options.outputFolder, 'uploadedPatchFilesTest.txt'))
        elif options.upload:
            uploadNextFile(inputListPath, options.outputFolder, options.upload, options.numThreads)
        else: # Update the data list
            # TODO: Clean up file paths!
            getDataList(inputListPath)


        endTime = time.time()

        print "Finished in " + str(endTime - startTime) + " seconds."
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

	# To more easily debug this program, comment out this catch block.
    # except Exception, err:
    #     sys.stderr.write( str(err) + '\n' )
    #     return 1


if __name__ == "__main__":
    sys.exit(main())
