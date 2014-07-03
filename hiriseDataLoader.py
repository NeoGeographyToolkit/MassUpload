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

import GoogleMapsEngine


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, ''' Script for grabbing HiRISE data files'''

    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


#--------------------------------------------------------------------------------

# fileType is the file name after the prefix
def generatePdsPath(filePrefix, fileType):
    """Generate the full PDS path for a given HiRISE data file"""
    
    # File prefix looks like this: PSP_009716_1755 or ESP_011984_1755
    
    # Extract the mission code (only a few possibilities)
    missionCode = filePrefix[0:3]
    
    # Determine which ORB folder this will be in (each contains 100 files)
    # - Looks like this: ORB_001300_001399
    frameNumber = filePrefix[4:9]
    orbDir = 'ORB_'+ frameNumber[0:4] +'00_'+ frameNumber[0:4] +'99'
    
    filename = filePrefix + fileType
    baseUrl  = "http://hirise-pds.lpl.arizona.edu/PDS/RDR/"
    fullUrl  = baseUrl + missionCode +"/"+ orbDir +"/"+ filePrefix +"/"+ filename

    #print filePrefix + fileType + ' -> ' + fullUrl
    return fullUrl

# missionCode should be ESP or PSP
def getDataList(outputFilePath, missionCode):
    """Populates a text file list of all available HiRISE RDR data"""
       
    print 'Updating HiRISE PDS data list...'
       
    #http://hirise-pds.lpl.arizona.edu/PDS/RDR/<PSP or ESP>/<ORB_PATH>/<PSP or ESP PATH>/<FILE_NAME.JP2>
    #http://hirise-pds.lpl.arizona.edu/PDS/RDR/ESP/ORB_011900_011999/ESP_011984_1755/
    #http://hirise-pds.lpl.arizona.edu/download/PDS/RDR/PSP/ORB_001400_001499/PSP_001430_1815/PSP_001430_1815_RED.JP2
    
   
    # TODO: Do another loop for ESP!
    baseUrl = "http://hirise-pds.lpl.arizona.edu/PDS/RDR/"+missionCode+"/"

    # Parse the top PDS level
    parsedIndexPage = BeautifulSoup(urllib2.urlopen((baseUrl)).read())

    outputFile = open(outputFilePath, 'w')

    # Loop through outermost directory
    for line in parsedIndexPage.findAll('a'):
        orbName = line.string
        if not 'ORB' in orbName: # Skip link up a directory
            continue
        orbPath = baseUrl + orbName
        
        #print 'Scanning directory ' + orbPath
        
        # Parse next directory level
        orbPage = BeautifulSoup(urllib2.urlopen((orbPath)).read())
        
        #print orbPage.prettify()
        
        # Loop through inner directory
        for line in orbPage.findAll('a'):
            pspName = line.string
            if not 'PSP' in pspName: # Skip link up a directory
                continue
            
            ## No need to loop through the next level; the file names are fixed.
            #pspPath = orbPath + pspName + pspName[:-1]
            
            # TODO: Don't store paths, just store the name.
            #outputFile.write(pspPath + '\n')
            outputFile.write(pspName[:-1] + '\n')
            
    outputFile.close()
    
    print 'Wrote updated data list to ' + outputFilePath


def uploadFile(filePrefix, remoteFilePath, logQueue, tempDir):
    """Uploads a remote file to maps engine"""
    
    print 'Uploading file ' + remoteFilePath
    
    
    localFilePath = os.path.join(tempDir, os.path.basename(remoteFilePath))
    if not os.path.exists(localFilePath):
        # Download the file
        cmd = 'wget ' + remoteFilePath + ' -O ' + localFilePath
        print cmd
        os.system(cmd)
    
    # Upload the file
    cmdArgs = [localFilePath, '--sensor', 0]
    #print cmdArgs
    assetId = GoogleMapsEngine.main(cmdArgs)
    #assetId = 12345
    
    #TODO: Check to make sure the file made it up!
    
    # Record that we uploaded the file
    logString = filePrefix + ', ' + str(assetId) + '\n' # Log path and the Maps Engine asset ID
    #print logString
    logQueue.put(logString)

        
    # Delete the file
    print 'rm ' + localFilePath
    #os.remove(localFilePath)
    
    print 'Finished uploading HiRISE data file'
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


# TODO: Select from different file types
def uploadNextFile(dataListPath, outputFolder, getColor=False, numFiles=1, numThreads=1):
    """Determines the next file to upload, uploads it, and logs it"""
    
    print 'Searching for next file to upload...'
    
    # Set up the output paths    
    #uploadedColorPath   = os.path.join(outputFolder, 'uploadedPatchFilesTest.txt')
    uploadedRedPath   = os.path.join(outputFolder, 'uploadedRed.csv')
    uploadedColorPath = os.path.join(outputFolder, 'uploadedColor.csv')
    
    logPath   = uploadedRedPath
    targetEnd = '_RED.JP2'
    if getColor:
        logPath   = uploadedColorPath
        targetEnd = '_COLOR.JP2'
    
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
        fullPath = generatePdsPath(line, targetEnd)
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
    





#--------------------------------------------------------------------------------

def main():

    print "Started hiriseDataLoader.py"

    try:
        try:
            usage = "usage: hiriseDataLoader.py <output folder> [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("-u", "--upload", dest="upload", type=int,
                              help="Upload this many files instead of fetching the list.")
            #parser.add_option("-o", "--output-folder", dest="outputFolder",
            #                  help="Specifies the folder to copy the data to.",
            #                  default='./')
            #parser.add_option("-n", "--name", dest="name",
            #                  help="Only get the data for the DTM with this name.",
            #                  default='')
            
            parser.add_option("--color", action="store_true", dest="getColor",
                              help="Retrieve COLOR image instead of RED images.")

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
        pspListPath = os.path.join(options.outputFolder, 'pspList.csv')
        espListPath = os.path.join(options.outputFolder, 'espList.csv')
        fullListFile = os.path.join(options.outputFolder, 'fullList.csv')

        # If we are not uploading data, update the data list
        if not options.upload:
            getDataList(pspListPath, 'PSP')
            getDataList(espListPath, 'ESP')
            # TODO: Concatenate files!!!
            #cmd = 'cat'
            #os.system(cmd)
        else:
            uploadNextFile('hiriseImageList_smallPatch.txt', options.outputFolder, options.getColor, options.upload, options.numThreads)
            #uploadNextFile(fullListFile, options.getColor, options.upload, options.numThreads)


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
