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

import mapsEngineUpload, IrgStringFunctions, IrgGeoFunctions, IrgIsisFunctions


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, ''' Script for grabbing CTX data files'''

    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


#--------------------------------------------------------------------------------

def getCreationTime(fileList):
    """Extract the file creation time and return in YYYY-MM-DDTHH:MM:SSZ format"""
    
    if len(fileList) < 2:
        raise Exception('Error, missing label file path!')
    filePath = fileList[1]
    
    timeString = ''
    f = open(filePath, 'r')
    for line in f:
        if 'StartTime' in line:
            print line
            timeString = IrgStringFunctions.getLineAfterText(line, '=')
            break
    f.close()
  
    if not timeString:
        raise Exception('Unable to find time string in file ' + filePath)
  
    # Get to the correct format
    timeString = timeString.strip()
    timeString = timeString[:-4] + 'Z'
  
    # The time string is almost in the correct format
    return timeString


def findAllDataSets(db, dataAddFunctionCall, sensorCode):
    '''Add all known data sets to the SQL database'''

    print 'Updating CTX PDS data list...'
   
    baseUrl = "http://pds-imaging.jpl.nasa.gov/data/mro/mars_reconnaissance_orbiter/ctx/"

    # Parse the top PDS level
    parsedIndexPage = BeautifulSoup(urllib2.urlopen((baseUrl)).read())  
    #print parsedIndexPage.prettify()

    # Loop through outermost directory
    for line in parsedIndexPage.findAll('a'):
        volumeName = line.string
        if (not 'mrox_' in volumeName) or ('txt' in volumeName): # Skip other links
            continue
        volumePath = baseUrl + volumeName + 'data/'
        
        print 'Scanning directory ' + volumePath
        
        # Parse next directory level
        volumePage = BeautifulSoup(urllib2.urlopen((volumePath)).read())
        print volumePage.prettify()
        raise Exception('Checking VOLUME PAGE DEBUG')
        
        # Loop through inner directory
        for line in volumePage.findAll('a'):
            dataName = line.string
            if not '.IMG' in dataName: # Skip other links
                continue
            
            ## Now store the file prefix and the volume we found it in
            #outputFile.write(dataName[:-4] +', '+ volumeName[:-1] +'\n')
            
            # Load volume name as subtype, data prefix as data set name.
            dataAddFunctionCall(db, sensorCode, volumeName[:-1], dataName[:-4], line)

    print 'Added CTX data files to database!'
    
    

def fetchAndPrepFile(setName, subtype, remoteURL, workDir):
    '''Retrieves a remote file and prepares it for upload'''
    
    print 'Uploading file ' + setName
    
    asuImagePath  = os.path.join(workDir, setName + '_noGeo.jp2') # Map projected image from ASU
    asuLabelPath  = os.path.join(workDir, setName + '_noGeo.lbl') # Label from ASU
    edrPath       = os.path.join(workDir, setName + '.IMG')       # Raw image from PDS
    #cubPath     = os.path.join(workDir, setName + '.cub')        # Output of mroctx2isis
    #calPath     = os.path.join(workDir, setName + '.cal.cub')    # Output of ctxcal
    mapPath       = os.path.join(workDir, setName + '.map.cub')   # Output of cam2map
    mapLabelPath  = os.path.join(workDir, setName + '.map.pvl')   # Specify projection to cam2map
    localFilePath = os.path.join(workDir, setName + '.tif')       # The output file we will upload
    
    # Generate the remote URLs from the data prefix and volume stored in these parameters
    asuImageUrl, asuLabelPath, edrUrl = generatePdsPath(setName, subtype)
    
    # We are using the label path in both projection cases
    if not os.path.exists(asuLabelPath):
        # Download the label file
        cmd = 'wget ' + labelUrl + ' -O ' + asuLabelPath
        print cmd
        os.system(cmd)

    if True: # Map project the EDR ourselves <-- Going with this approach!
        print 'Projecting the EDR image using ISIS...'
        
        if not os.path.exists(mapLabelPath):
            # Generate the map label file
            numLinesCopied = extractPvlSection(asuLabelPath, mapLabelPath, "Mapping")
            if numLinesCopied < 10:
                raise Exception('Failed to copy map data from file ' + asuLabelPath)
        
        if not os.path.exists(edrPath):
            # Download the EDR file
            cmd = 'wget ' + edrUrl + ' -O ' + edrPath
            print cmd
            os.system(cmd)
      
        # Convert and apply calibration to the CTX file
        calPath = IrgIsisFunctions.prepareCtxImage(edrPath, workDir, True)
        
        if not os.path.exists(mapPath):
            # Generate the map projected file
            cmd = 'cam2map matchmap=True from=' + calPath + ' to=' + mapPath + ' map='+mapLabelPath
            print cmd
            os.system(cmd)
        
        if not os.path.exists(localFilePath):
            # Generate the final image to upload
            cmd = 'gdal_translate -of GTiff ' + mapPath + ' ' + localFilePath
            print cmd
            os.system(cmd)
        
        # Clean up intermediate files    
        #os.remove(mapLabelPath)
        #os.remove(edrPath)
        #os.remove(calPath)
        #os.remove(mapPath)
        
    else: # Use the map projected image from the ASU web site
        print 'Using ASU projected image...'
        
        if not os.path.exists(asuImagePath):
            # Download the image file
            cmd = 'wget ' + asuImageUrl + ' -O ' + asuImagePath
            print cmd
            os.system(cmd)

        if not os.path.exists(localFilePath):
            # Correct the file - The JP2 file from ASU needs the geo data from the label file!
            cmd = 'addGeoToAsuCtxJp2.py --label '+ asuLabelPath +' '+ asuImagePath +' '+ localFilePath
            print cmd
            os.system(cmd)
            
            if not os.path.exists(localFilePath):
                raise Exception('Script to add geo data to JP2 file failed!')
            
        # Clean up
        #os.remove(asuImagePath)
    
    # Two local files are left around, the first should be uploaded.
    return [localFilePath, asuLabelPath]
        


#--------------------------------------------------------------------------------


def extractPvlSection(inputPath, outputPath, sectionName):
    """Copies a section of a PVL file to another file"""
    
    startLine = "Group = " + sectionName

    inputFile  = open(inputPath, 'r')
    outputFile = open(outputPath, 'w')
    copyingLines = False
    numLines = 0
    for line in inputFile:
        if startLine in line:
            copyingLines = True # Start copying this section
            
        if copyingLines: # Copy this line to the output file
            
            # Specific checks for the mapping section!
            
            # Force this value to zero
            if "CenterLongitude" in line: 
                line = "CenterLongitude    = 0.0\n"
            # Force this value to -180 to 180 range
            if "LongitudeDomain" in line: 
                line = "LongitudeDomain    = 180\n"
            
            outputFile.write(line)
            numLines = numLines + 1
            
        if copyingLines and ("End_Group" in line):
            break # Quit parsing the file
    inputFile.close()
    outputFile.close()
    
    return numLines # Return the number of lines copied

def generatePdsPath(filePrefix, volume):
    """Generate the full PDS path for a given CTX data file"""
    
    # File prefix looks something like this: B08_012841_1751_XN_04S222W
    # Volume ID looks like this: mrox_0738
    
    # Grab the file directly from the ASU map projected database.
    imageUrl = ("http://image.mars.asu.edu/stream/"+filePrefix+
               ".jp2?image=/mars/images/ctx/"+volume+"/prj_full/"+filePrefix+".jp2")

    labelUrl = ("http://image.mars.asu.edu/stream/"+filePrefix+
                ".scyl.isis.hdr?image=/mars/images/ctx/"+volume+"/stage/"
                +filePrefix+".scyl.isis.hdr")

    edrUrl   = ("http://pds-imaging.jpl.nasa.gov/data/mro/mars_reconnaissance_orbiter/ctx/"
                +volume+"/data/"+filePrefix+".IMG")

    return (imageUrl, labelUrl, edrUrl)

def getDataList(outputFilePath):
    """Populates a text file list of all available CTX data and its volume"""
       
    print 'Updating CTX PDS data list...'
       
   
    baseUrl = "http://pds-imaging.jpl.nasa.gov/data/mro/mars_reconnaissance_orbiter/ctx/"

    # Parse the top PDS level
    parsedIndexPage = BeautifulSoup(urllib2.urlopen((baseUrl)).read())  
    #print parsedIndexPage.prettify()

    outputFile = open(outputFilePath, 'w')

    # Loop through outermost directory
    for line in parsedIndexPage.findAll('a'):
        volumeName = line.string
        if (not 'mrox_' in volumeName) or ('txt' in volumeName): # Skip other links
            continue
        volumePath = baseUrl + volumeName + 'data/'
        
        print 'Scanning directory ' + volumePath
        
        # Parse next directory level
        volumePage = BeautifulSoup(urllib2.urlopen((volumePath)).read())
        #print volumePage.prettify()
        
        # Loop through inner directory
        for line in volumePage.findAll('a'):
            dataName = line.string
            if not '.IMG' in dataName: # Skip other links
                continue
            
            # Now store the file prefix and the volume we found it in
            outputFile.write(dataName[:-4] +', '+ volumeName[:-1] +'\n')

    outputFile.close()
    
    print 'Wrote updated data list to ' + outputFilePath


def uploadFile(filePrefix, imageUrl, labelUrl, edrUrl, reproject, logQueue, tempDir):
    """Uploads a remote file to maps engine"""
    
    print 'Uploading file ' + filePrefix
    
    asuImagePath  = os.path.join(tempDir, filePrefix + '_noGeo.jp2') # Map projected image from ASU
    asuLabelPath  = os.path.join(tempDir, filePrefix + '_noGeo.lbl') # Label from ASU
    edrPath       = os.path.join(tempDir, filePrefix + '.IMG')       # Raw image from PDS
    #cubPath     = os.path.join(tempDir, filePrefix + '.cub')        # Output of mroctx2isis
    #calPath     = os.path.join(tempDir, filePrefix + '.cal.cub')    # Output of ctxcal
    mapPath       = os.path.join(tempDir, filePrefix + '.map.cub')   # Output of cam2map
    mapLabelPath  = os.path.join(tempDir, filePrefix + '.map.pvl')   # Specify projection to cam2map
        
    # We are using the label path in both projection cases
    if not os.path.exists(asuLabelPath):
        # Download the label file
        cmd = 'wget ' + labelUrl + ' -O ' + asuLabelPath
        print cmd
        os.system(cmd)

    timeString = getCreationTime(asuLabelPath)

    if reproject: # Map project the EDR ourselves
        print 'Projecting the EDR image using ISIS...'
        
        if not os.path.exists(mapLabelPath):
            # Generate the map label file
            numLinesCopied = extractPvlSection(asuLabelPath, mapLabelPath, "Mapping")
            if numLinesCopied < 10:
                raise Exception('Failed to copy map data from file ' + asuLabelPath)
        
        if not os.path.exists(edrPath):
            # Download the EDR file
            cmd = 'wget ' + edrUrl + ' -O ' + edrPath
            print cmd
            os.system(cmd)
      
        # Convert and apply calibration to the CTX file
        calPath = IrgIsisFunctions.prepareCtxImage(edrPath, tempDir, True)
        
        if not os.path.exists(mapPath):
            # Generate the map projected file
            cmd = 'cam2map matchmap=True from=' + calPath + ' to=' + mapPath + ' map='+mapLabelPath
            print cmd
            os.system(cmd)
        
        localFilePath = os.path.join(tempDir, filePrefix + '.tif') # The output file we will upload
        if not os.path.exists(localFilePath):
            # Generate the final image to upload
            cmd = 'gdal_translate -of GTiff ' + mapPath + ' ' + localFilePath
            print cmd
            os.system(cmd)
        
        # Clean up intermediate files    
        #os.remove(mapLabelPath)
        #os.remove(edrPath)
        #os.remove(calPath)
        #os.remove(mapPath)
        
    else: # Use the map projected image from the ASU web site
        print 'Using ASU projected image...'
        
        if not os.path.exists(asuImagePath):
            # Download the image file
            cmd = 'wget ' + imageUrl + ' -O ' + asuImagePath
            print cmd
            os.system(cmd)

        localFilePath = os.path.join(tempDir, filePrefix + '.jp2') # The output file we will upload
        if not os.path.exists(localFilePath):
            # Correct the file - The JP2 file from ASU needs the geo data from the label file!
            cmd = 'addGeoToAsuCtxJp2.py --label '+ asuLabelPath +' '+ asuImagePath +' '+ localFilePath
            print cmd
            os.system(cmd)
            
            if not os.path.exists(localFilePath):
                raise Exception('Script to add geo data to JP2 file failed!')
            
        # Clean up
        #os.remove(asuImagePath)
    
    # Upload the file
    cmdArgs = [localFilePath, '--sensor', '2', '--acqTime', timeString]
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
    #os.remove(asuLabelPath)
    
    print 'Finished uploading CTX data file'
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
def uploadNextFile(dataListPath, outputFolder, reproject=False, numFiles=1, numThreads=1):
    """Determines the next file to upload, uploads it, and logs it"""
    
    print 'Searching for next file to upload...'
    
    # Set up the output paths    
    logPath = os.path.join(outputFolder, 'uploadedPatchFiles.txt')
    
    inFile = open(dataListPath,    'r')

    # Get the last line read (could be more efficient)
    lastUploadedLine = ''
    if os.path.exists(logPath):
        outFile = open(logPath, 'r')
        for line in outFile:
            lastUploadedLine = line
        outFile.close()
    lastUploadedLine = lastUploadedLine.split(',')[0].strip() # Remove the asset ID from the string
    print '#' + lastUploadedLine + '#'
    
    # Now find that line in the input file list
    linesToProcess = []
    breakNext = -1
    for line in inFile:
        # Remove /n
        line = line.strip()
        prefix = line.split(',')[0].strip() # Strip off the volume label
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
        if (prefix == lastUploadedLine): # Found the last one downloaded
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
        # Each line contains a prefix and a volume label seperated by a comma
        prefix, volume = line.split(',')
        prefix = prefix.strip()
        volume = volume.strip()
        imageUrl, labelUrl, edrUrl = generatePdsPath(prefix, volume)
        jobResults.append(pool.apply_async(uploadFile, args=(prefix, imageUrl, labelUrl, edrUrl, reproject, queue, outputFolder)))
    
    
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
        status = mapsEngineUpload.checkIfFileIsLoaded(bearerToken, assetId)
        
        if not status:
            print 'Prefix ' + prefix + ' was not uploaded correctly!'
            # TODO: Do something about it!
        
    outFile.close()


#--------------------------------------------------------------------------------

def main():

    print "Started ctxDataLoader.py"

    try:
        try:
            usage = "usage: ctxDataLoader.py <output folder> [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("-u", "--upload", dest="upload", type=int,
                              help="Upload this many files instead of fetching the list.")

            parser.add_option("--reproject", action="store_true", dest="reproject", default=False,
                              help="Project the images ourselves instead of using the ASU projections.")

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
        linkListPath = os.path.join(options.outputFolder, 'ctxDataList.csv')

        if options.checkUploads:
            checkUploads(os.path.join(options.outputFolder, 'uploadedPatchFilesTest.txt'))
        elif options.upload:
            uploadNextFile('ctxImageList_smallPatch.csv', options.outputFolder, options.reproject, options.upload, options.numThreads)
            #uploadNextFile(fullListFile, options.getColor, options.upload, options.numThreads)
        else: # Update the data list
            getDataList(linkListPath)


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
