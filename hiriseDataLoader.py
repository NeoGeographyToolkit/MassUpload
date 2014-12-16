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

from bs4 import BeautifulSoup

import os, glob, optparse, re, shutil, subprocess, string, time, urllib, urllib2

import multiprocessing

import mapsEngineUpload, IrgStringFunctions, IrgGeoFunctions, IrgFileFunctions

THIS_FOLDER  = os.path.dirname(os.path.abspath(__file__))
FIX_JP2_TOOL = os.path.join(THIS_FOLDER, 'fix_jp2')


#--------------------------------------------------------------------------------
# Functions needed for the unified data loader

def getUploadList(fileList):
    '''Returns the subset of the fileList that needs to be uploaded to GME'''
    return [fileList[0]] # No need to upload the lable file

def getCreationTime(fileList):
    """Extract the file creation time and return in YYYY-MM-DDTHH:MM:SSZ format"""
    
    if 'DTE' in fileList[0]:  # DEM files do not have a creation time!
        #raise Exception('Error, missing label file path!')
        return '1900-01-01T00:00:00Z' # Return a flag date
    filePath = fileList[1]
    
    timeString = ''
    f = open(filePath, 'r')
    for line in f:
        if 'STOP_TIME' in line:
            timeString = IrgStringFunctions.getLineAfterText(line, '=')
            break
    f.close()
  
    if not timeString:
        raise Exception('Unable to find time string in file ' + filePath)
  
    # Get to the correct format
    timeString = timeString.strip()
    timeString = timeString[:-4] + 'Z'
  
    return timeString

def getBoundingBox(fileList):
    """Return the bounding box for this data set in the format (minLon, maxLon, minLat, maxLat)"""
    if len(fileList) == 2: # Read BB from the label file
        return IrgGeoFunctions.getBoundingBoxFromIsisLabel(fileList[1])
    else: # No label file, read it from the the main file
        return IrgGeoFunctions.getImageBoundingBox(fileList[0]) # This information is also available in the IMG file header


def findAllDataSets(db, dataAddFunctionCall, sensorCode):
    '''Add all known data sets to the SQL database'''
    
    print 'Updating HiRISE PDS data list...'
       
    #http://hirise-pds.lpl.arizona.edu/PDS/RDR/<PSP or ESP>/<ORB_PATH>/<PSP or ESP PATH>/<FILE_NAME.JP2>
    #http://hirise-pds.lpl.arizona.edu/PDS/RDR/ESP/ORB_011900_011999/ESP_011984_1755/
    #http://hirise-pds.lpl.arizona.edu/download/PDS/RDR/PSP/ORB_001400_001499/PSP_001430_1815/PSP_001430_1815_RED.JP2
    
    # We are interested in the primary and extended mission phases
    missionCodeList = ['PSP', 'ESP']
    for missionCode in missionCodeList:
    
        print 'Finding ' + missionCode + ' image files...'
    
        baseUrl = "http://hirise-pds.lpl.arizona.edu/PDS/RDR/"+missionCode+"/"
    
        # Parse the top PDS level
        parsedIndexPage = BeautifulSoup(urllib2.urlopen((baseUrl)).read())
       
        # Loop through outermost directory
        for line in parsedIndexPage.findAll('a'):

            orbName = line.string
            if not 'ORB' in orbName: # Skip link up a directory
                continue
            orbPath = baseUrl + orbName
            
            # Parse next directory level
            orbPage = BeautifulSoup(urllib2.urlopen((orbPath)).read())
            
            # Loop through inner directory
            for line in orbPage.findAll('a'):
                xspName = line.string
                if not missionCode in xspName: # Skip link up a directory
                    continue
                dataPrefix = xspName[:-1] # Strip off the trailing '/'
                
                ## No need to loop through the next level; the file names are fixed.
                redUrl   = orbPath + dataPrefix + '/' + dataPrefix + '_RED.JP2'
                colorUrl = orbPath + dataPrefix + '/' + dataPrefix + '_COLOR.JP2'
                
                # Add both the RED and COLOR files to the database.
                dataAddFunctionCall(db, sensorCode, 'RED',   dataPrefix, redUrl)
                dataAddFunctionCall(db, sensorCode, 'COLOR', dataPrefix, colorUrl)


        print 'Finding ' + missionCode + ' DEM files...'

        # Now go through and add the DTMs for this mission code
        baseDemUrl = 'http://www.uahirise.org/PDS/DTM/' + missionCode + '/'
    
        # Parse the top PDS level
        parsedIndexPage = BeautifulSoup(urllib2.urlopen((baseDemUrl)).read())
    
        # Loop through outermost directory
        for line in parsedIndexPage.findAll('a'):
            orbName = line.string
            if not 'ORB' in orbName: # Skip link up a directory
                continue
            orbPath = baseDemUrl + orbName
            
            #print 'Scanning directory ' + orbPath
            
            # Parse next directory level
            orbPage = BeautifulSoup(urllib2.urlopen((orbPath)).read())
            
            #print orbPage.prettify()
            
            # Loop through inner directory
            for line in orbPage.findAll('a'):
                xspName = line.string
                if not missionCode in xspName: # Skip link up a directory
                    continue
                dataPrefix = xspName[:-1] # Strip off the trailing '/'
                
                stereoPageUrl = orbPath + dataPrefix + '/'
                stereoPage    = BeautifulSoup(urllib2.urlopen((stereoPageUrl)).read())
                
                # Pick out the link with the DEM
                for line in stereoPage.findAll('a'):
                    fileName = line.string
                    if '.IMG' in fileName: # Only the DEM has a .IMG extension
                        demUrl = stereoPageUrl + fileName
                        dataAddFunctionCall(db, sensorCode, 'DEM', dataPrefix, demUrl)


def fetchAndPrepFile(setName, subtype, remoteURL, workDir):
    '''Retrieves a remote file and prepares it for upload'''
    
    #print 'Uploading file ' + setName
  
    if subtype != 'DEM': # Handles RED and COLOR images
        # The label file URL is the same as the image but with a different extension
        remoteLabelURL = getLabelPathFromImagePath(remoteURL)
    
        localFilePath  = os.path.join(workDir, os.path.basename(remoteURL))
        localLabelPath = os.path.join(workDir, os.path.basename(remoteLabelURL))

        if not os.path.exists(localLabelPath):
            # Try to get the label locally!
            pdsStart     = remoteLabelURL.find('PDS')
            localPdsPath = os.path.join('/HiRISE/Data/', remoteLabelURL[pdsStart:])
            print localPdsPath
            if os.path.exists(localPdsPath): # File available locally, just copy it!
                cmd = 'cp ' + localPdsPath +' '+ localLabelPath
            else:  # Download the image
                cmd = 'wget ' + remoteLabelURL + ' -O ' + localLabelPath
            print cmd
            os.system(cmd)

        if not IrgFileFunctions.fileIsNonZero(localLabelPath):
            raise Exception('Unable to download from URL: ' + remoteLabelURL)

        # Check the projection type
        projType = IrgGeoFunctions.getProjectionFromIsisLabel(localLabelPath)
        if projType == 'POLAR STEREOGRAPHIC':
            os.remove(localLabelPath)
            raise Exception('POLAR STEREOGRAPHIC images on hold until Google fixes a bug!')
    
        if not os.path.exists(localFilePath): # Need to get it from somewhere
            # Try to get the image locally!
            pdsStart     = remoteURL.find('PDS')
            localPdsPath = os.path.join('/HiRISE/Data/', remoteURL[pdsStart:])
            if os.path.exists(localPdsPath): # File available locally, just copy it!
                cmd = 'cp ' + localPdsPath +' '+ localFilePath
            else:  # Download the image
                cmd = 'wget ' + remoteURL + ' -O ' + localFilePath
            print cmd
            os.system(cmd)


        if not IrgFileFunctions.fileIsNonZero(localFilePath):
            raise Exception('Unable to download from URL: ' + remoteURL)
    
        # Call code to fix the header information in the JP2 file!
        cmd = FIX_JP2_TOOL +' '+ localFilePath
        print cmd
        os.system(cmd)

        # First file is for upload, second contains the timestamp.
        return [localFilePath, localLabelPath]

    else: # Handle DEMs
        
        # For DEMs there is no label file
        localFilePath = os.path.join(workDir, os.path.basename(remoteURL))
        if not os.path.exists(localFilePath):
            # Download the image
            cmd = 'wget ' + remoteURL + ' -O ' + localFilePath
            print cmd
            os.system(cmd)
        if not IrgFileFunctions.fileIsNonZero(localFilePath): # Make sure we got the file
            raise Exception('Unable to download from URL: ' + remoteURL)

        # Generate a header file from the IMG file
        localLabelPath = localFilePath[:-4] + '.LBL'
        cmd = 'head -n 90 ' + localFilePath +' >  '+ localLabelPath
        print cmd
        os.system(cmd)

        # Check if this is a polar stereographic image
        f = open(localLabelPath)
        for line in f:
            if ("MAP_PROJECTION_TYPE" in line) and ("POLAR STEREOGRAPHIC" in line):
                raise Exception('POLAR STEREOGRAPHIC DEMs on hold until Google fixes a bug!')
        f.close()
            
        # Convert from IMG to TIF
        tiffFilePath = localFilePath[:-4] + '.TIF'
        if not os.path.exists(tiffFilePath):
            cmd = 'gdal_translate -of GTiff ' + localFilePath +' '+ tiffFilePath
            print cmd
            os.system(cmd)
        
            # Correct projected coordinates problems
            cmd = 'python /home/pirl/smcmich1/repo/Tools/geoTiffTool.py --normalize-eqc-lon ' + tiffFilePath
            print cmd
            os.system(cmd)

        os.remove(localFilePath) # Clean up source image
        return [tiffFilePath, localLabelPath] 


#--------------------------------------------------------------------------------

def getLabelPathFromImagePath(imagePath):
    '''Given the image path, return the corresponding label path'''
    # Just replace the extension!
    return (imagePath[:-4] + '.LBL')


# fileType is the file name after the prefix
def generatePdsPaths(filePrefix, fileType):
    """Generate the full PDS path for a given HiRISE data file"""
    
    # File prefix looks like this: PSP_009716_1755 or ESP_011984_1755
    
    # Extract the mission code (only a few possibilities)
    missionCode = filePrefix[0:3]
    
    # Determine which ORB folder this will be in (each contains 100 files)
    # - Looks like this: ORB_001300_001399
    frameNumber = filePrefix[4:9]
    orbDir = 'ORB_'+ frameNumber[0:4] +'00_'+ frameNumber[0:4] +'99'
    
    fileName  = filePrefix + fileType
    labelName = os.path.splitext(fileName)[0] + '.LBL'
    baseUrl   = "http://hirise-pds.lpl.arizona.edu/PDS/RDR/"
    imageUrl  = baseUrl + missionCode +"/"+ orbDir +"/"+ filePrefix +"/"+ fileName
    labelUrl  = baseUrl + missionCode +"/"+ orbDir +"/"+ filePrefix +"/"+ labelName
    
    return (imageUrl, labelUrl)
#
## missionCode should be ESP or PSP
#def getDataList(outputFilePath, missionCode):
#    """Populates a text file list of all available HiRISE RDR data"""
#       
#    print 'Updating HiRISE PDS data list...'
#       
#    #http://hirise-pds.lpl.arizona.edu/PDS/RDR/<PSP or ESP>/<ORB_PATH>/<PSP or ESP PATH>/<FILE_NAME.JP2>
#    #http://hirise-pds.lpl.arizona.edu/PDS/RDR/ESP/ORB_011900_011999/ESP_011984_1755/
#    #http://hirise-pds.lpl.arizona.edu/download/PDS/RDR/PSP/ORB_001400_001499/PSP_001430_1815/PSP_001430_1815_RED.JP2
#    
#   
#    # TODO: Do another loop for ESP!
#    baseUrl = "http://hirise-pds.lpl.arizona.edu/PDS/RDR/"+missionCode+"/"
#
#    # Parse the top PDS level
#    parsedIndexPage = BeautifulSoup(urllib2.urlopen((baseUrl)).read())
#
#    outputFile = open(outputFilePath, 'w')
#
#    # Loop through outermost directory
#    for line in parsedIndexPage.findAll('a'):
#        orbName = line.string
#        if not 'ORB' in orbName: # Skip link up a directory
#            continue
#        orbPath = baseUrl + orbName
#        
#        #print 'Scanning directory ' + orbPath
#        
#        # Parse next directory level
#        orbPage = BeautifulSoup(urllib2.urlopen((orbPath)).read())
#        
#        #print orbPage.prettify()
#        
#        # Loop through inner directory
#        for line in orbPage.findAll('a'):
#            pspName = line.string
#            if not 'PSP' in pspName: # Skip link up a directory
#                continue
#            
#            ## No need to loop through the next level; the file names are fixed.
#            #pspPath = orbPath + pspName + pspName[:-1]
#            
#            # TODO: Don't store paths, just store the name.
#            #outputFile.write(pspPath + '\n')
#            outputFile.write(pspName[:-1] + '\n')
#            
#    outputFile.close()
#    
#    print 'Wrote updated data list to ' + outputFilePath
#
#
#def uploadFile(filePrefix, remoteFilePath, remoteLabelPath, logQueue, tempDir):
#    """Uploads a remote file to maps engine"""
#    
#    print 'Uploading file ' + remoteFilePath
#    
#    localFilePath  = os.path.join(tempDir, os.path.basename(remoteFilePath))
#    localLabelPath = os.path.join(tempDir, os.path.basename(remoteLabelPath))
#    
#    
#    if not os.path.exists(localFilePath):
#        # Download the file
#        cmd = 'wget ' + remoteFilePath + ' -O ' + localFilePath
#        print cmd
#        os.system(cmd)
#
#    if not os.path.exists(localLabelPath):
#        # Download the file
#        cmd = 'wget ' + remoteLabelPath + ' -O ' + localLabelPath
#        print cmd
#        os.system(cmd)
#
#
#    timeString = getCreationTime(localLabelPath)
#    
#    # Upload the file
#    cmdArgs = [localFilePath, '--sensor', '0', '--acqTime', timeString]
#    print cmdArgs
#    assetId = mapsEngineUpload.main(cmdArgs)
#    #assetId = 12345
#    
#    #TODO: Check to make sure the file made it up!
#
#    # Find out the bounding box of the file and generate a log string
#    fileBbox = IrgGeoFunctions.getImageBoundingBox(localFilePath)
#    bboxString = ('Bbox: ' + str(fileBbox[0]) +' '+ str(fileBbox[1]) +' '+ str(fileBbox[2]) +' '+ str(fileBbox[3]))
#    
#    # Record that we uploaded the file
#    logString = filePrefix +', '+ str(assetId) +', '+ bboxString + '\n' # Log path and the Maps Engine asset ID
#    print logString
#    logQueue.put(logString)
#
#        
#    # Delete the file
#    print 'rm ' + localFilePath
#    #os.remove(localFilePath)
#    
#    print 'Finished uploading HiRISE data file'
#    return assetId
#
#
#def logWriter(logQueue, logPath):
#    '''listens for messages on the q, writes to file. '''
#
#    print 'Starting writer'
#    f = open(logPath, 'a') 
#    while 1: # Run until the process is killed
#        message = logQueue.get() # Wait for a new message
#        if message == 'stop_queue': # Check for the quit signal
#            break
#        f.write(str(message))
#        f.flush()
#    f.close()
#    print 'Writer stopped'
#
#
## TODO: Select from different file types
#def uploadNextFile(dataListPath, outputFolder, getColor=False, numFiles=1, numThreads=1):
#    """Determines the next file to upload, uploads it, and logs it"""
#    
#    print 'Searching for next file to upload...'
#    
#    # Set up the output paths    
#    #uploadedColorPath   = os.path.join(outputFolder, 'uploadedPatchFilesTest.txt')
#    uploadedRedPath   = os.path.join(outputFolder, 'uploadedRed.csv')
#    uploadedColorPath = os.path.join(outputFolder, 'uploadedColor.csv')
#    
#    logPath   = uploadedRedPath
#    targetEnd = '_RED.JP2'
#    if getColor:
#        logPath   = uploadedColorPath
#        targetEnd = '_COLOR.JP2'
#    
#    inFile = open(dataListPath,    'r')
#
#    # Get the last line read (could be more efficient)
#    lastUploadedLine = ''
#    if os.path.exists(logPath):
#        outFile = open(logPath, 'r')
#        for line in outFile:
#            lastUploadedLine = line
#        outFile.close()
#    lastUploadedLine = lastUploadedLine.split(',')[0].strip() # Remove the asset ID from the string
#    #print '#' + lastUploadedLine + '#'
#    
#    # Now find that line in the input file list
#    linesToProcess = []
#    breakNext = -1
#    for line in inFile:
#        line = line.strip() # Remove \n from the line
#        #print '#' + line + '#'
#        
#        if breakNext > 0:
#            linesToProcess.append(line) # Record this line and move on to the next
#            breakNext = breakNext - 1
#            continue
#        if breakNext == 0: # No more lines to get!
#            break
#
#        if (lastUploadedLine == ''): # This is the first file!
#            linesToProcess.append(line)
#            breakNext = numFiles - 1 # Still get the next N-1 files
#            continue
#        if (line == lastUploadedLine): # Found the last one downloaded
#            breakNext = numFiles # Get the next N-1 files
#            continue
#    
#    if linesToProcess == []: # Make sure we found the next lines
#        raise Exception('Failed to find line that comes after: ' + lastUploadedLine)
#    
#    # Create processing pool
#    print 'Spawning ' + str(numThreads) + ' worker threads'
#    pool = multiprocessing.Pool(processes=numThreads+1) # One extra thread for the logWriter
#
#    # Create multiprocessing manager and a queue
#    manager = multiprocessing.Manager()
#    queue   = manager.Queue()    
#
#    # Start up the log writing thread
#    print 'Writing output to file ' + logPath
#    logResult = pool.apply_async(logWriter, args=(queue, logPath))
#    
#    # For each line generate the full download path
#    jobResults = []
#    for line in linesToProcess:
#        remoteImagePath, remoteLabelPath = generatePdsPaths(line, targetEnd)
#        jobResults.append(pool.apply_async(uploadFile, args=(line, remoteImagePath, remoteLabelPath, queue, outputFolder)))
#    
#    
#    # Wait until all threads have finished
#    print 'Waiting for all threads to complete...'
#    for r in jobResults:
#        r.get()
#    
#    # Stop the queue and all the threads
#    print 'Cleaning up...'
#    queue.put('stop_queue')
#    pool.close()
#    pool.join()
#    
#    print 'All threads finished!'
#   
#    
#    return True
#    
#
#
#def checkUploads(logPath):
#
#    print 'Checking the status of uploaded files...'    
#
#    # Get server authorization and hold on to the token
#    bearerToken = mapsEngineUpload.authorize()
#
#    if not os.path.exists(logPath):
#        raise Exception('Input log file ' + logPath + ' does not exist!')
#        
#    outFile = open(logPath, 'r')
#    for line in outFile:
#        # Extract the asset ID
#        prefix, assetId, bbox = line.split(',')
#        
#        # Check if this asset was uploaded
#        status = mapsEngineUpload.checkIfFileIsLoaded(bearerToken, assetId)
#        
#        if not status:
#            print 'Prefix ' + prefix + ' was not uploaded correctly!'
#            # TODO: Do something about it!
#        
#    outFile.close()
#
#
##--------------------------------------------------------------------------------
#
#def main():
#
#    print "Started hiriseDataLoader.py"
#
#    try:
#        try:
#            usage = "usage: hiriseDataLoader.py <output folder> [--help][--manual]\n  "
#            parser = optparse.OptionParser(usage=usage)
#            parser.add_option("-u", "--upload", dest="upload", type=int,
#                              help="Upload this many files instead of fetching the list.")
#            #parser.add_option("-o", "--output-folder", dest="outputFolder",
#            #                  help="Specifies the folder to copy the data to.",
#            #                  default='./')
#            #parser.add_option("-n", "--name", dest="name",
#            #                  help="Only get the data for the DTM with this name.",
#            #                  default='')
#            
#            parser.add_option("--color", action="store_true", dest="getColor",
#                              help="Retrieve COLOR image instead of RED images.")
#
#            parser.add_option("--checkUploads", action="store_true", default=False,
#                                        dest="checkUploads",  help="Verify that all uploaded files actually made it up.")
#
#            parser.add_option("--threads", type="int", dest="numThreads", default=1,
#                              help="Number of threads to use.")
#                              
#            parser.add_option("--manual", action="callback", callback=man,
#                              help="Read the manual.")
#            (options, args) = parser.parse_args()
#
#            if len(args) < 1:
#                raise Exception('Missing output folder!')
#                return 1;
#            options.outputFolder = args[0]
#
#
#        except optparse.OptionError, msg:
#            raise Usage(msg)
#
#        print "Beginning processing....."
#
#        startTime = time.time()
#
#        if not os.path.exists(options.outputFolder):
#            os.mkdir(options.outputFolder)
#
#        # These are input lists that just show the available images
#        pspListPath = os.path.join(options.outputFolder, 'pspList.csv')
#        espListPath = os.path.join(options.outputFolder, 'espList.csv')
#        fullListFile = os.path.join(options.outputFolder, 'fullList.csv')
#
#        if options.checkUploads:
#            checkUploads(os.path.join(options.outputFolder, 'uploadedPatchFilesTest.txt'))
#        elif options.upload:
#            uploadNextFile('hiriseImageList_smallPatch.txt', options.outputFolder, options.getColor, options.upload, options.numThreads)
#            #uploadNextFile(fullListFile, options.getColor, options.upload, options.numThreads)
#        else: # Update the data list
#            # TODO: Clean up file paths!
#            getDataList(pspListPath, 'PSP')
#            getDataList(espListPath, 'ESP')
#            # TODO: Concatenate files!!!
#            #cmd = 'cat'
#            #os.system(cmd)
#
#
#        endTime = time.time()
#
#        print "Finished in " + str(endTime - startTime) + " seconds."
#        return 0
#
#    except Usage, err:
#        print >>sys.stderr, err.msg
#        return 2
#
#	# To more easily debug this program, comment out this catch block.
#    # except Exception, err:
#    #     sys.stderr.write( str(err) + '\n' )
#    #     return 1
#
#
#if __name__ == "__main__":
#    sys.exit(main())
