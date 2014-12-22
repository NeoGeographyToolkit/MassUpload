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

import common

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


def findAllDataSets(db, sensorCode):
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
                common.addDataRecord(db, sensorCode, 'RED',   dataPrefix, redUrl)
                common.addDataRecord(db, sensorCode, 'COLOR', dataPrefix, colorUrl)


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
                        common.addDataRecord(db, sensorCode, 'DEM', dataPrefix, demUrl)


MAX_POLAR_UPLOAD_WIDTH  = 20000 # TODO: What is this really?
POLAR_WIDTH_CHUNK_SIZE  = 20000
CHUNK_SETNAME_SEPERATOR = '___' # Appended to the setname field to indicate which chuck this is

def getChunkNum(setName):
    '''Determine which chunk number a DB entry is for'''
    if CHUNK_SETNAME_SEPERATOR in setName:
        # The chunk number is an integer after the seperator at the end of setName.
        chunkNumPos = setName.rfind(CHUNK_SETNAME_SEPERATOR) + len(CHUNK_SETNAME_SEPERATOR)
        chunkNum    = int(setName[chunkNumPos:])
    else: # The main DB entry is associated with chunk index zero
        chunkNum = 0
    return chunkNum

def makeChunkSetName(setName, chunkNum):
    '''Returns the set name used for a given chunk'''
    if chunkNum == 0:
        return setName
    else:
        return setname + CHUNK_SETNAME_SEPERATOR + str(chunkNum)

def getChunkAreaString(width, height, chunkNum):
    # Each chunk is the full height of the image
    y      = 0
    yRange = height
    x      = chunkNum*POLAR_WIDTH_CHUNK_SIZE
    xRange = POLAR_WIDTH_CHUNK_SIZE # Don't need to get this exactly right since the crop tool clamps to the border
    return  '%d,%d,%d,%d' % (x, y, xRange, yRange)


def jp2ToImgAndBack(jp2Path, labelPath, outputJp2Path, areaString=None):
    '''Convert a Jp2 to an IMG and back.  Can crop and restore missing metadata.'''
    
    cmd = ('/home/pirl/smcmich1/programs/PDS_JP2-3.17_Linux-x86_64/bin/JP2_to_PDS -Force -Input '  + jp2Path +
                                                                                       ' -Output ' + tempImgPath +
                                                                                       ' -LAbel '  + labelPath)
    if areaString:
        cmd += ' -Area ' + areaString
    print cmd
    os.system(cmd)
    if not IrgFileFunctions.fileIsNonZero(tempImgPath):
        raise Exception('Failed to convert chunk to IMG file: ' + jp2Path)
    
    cmd = ('/home/pirl/smcmich1/programs/PDS_JP2-3.17_Linux-x86_64/bin/PDS_to_JP2 -Force -Geotiff -Input '  + tempImgPath +
                                                                                                ' -Output ' + outputJp2Path)
    print cmd
    os.system(cmd)
    if not IrgFileFunctions.fileIsNonZero(localChunkPath):
        raise Exception('Failed to convert chunk to JP2 file: ' + outputJp2Path)
    
    os.remove(tempImgPath)
    
    return True

def writeFakeLabelFromJp2(jp2Path):
    '''For JP2 files with no associated LBL file, create an artificial one.'''

    # Read information about the input JP2 file
    infoDict = getImageGeoInfo(jp2Path, False)

    # Create the new label file
    labelPath = os.path.splitext(jp2Path)[0] + '.LBL'
    f = open(labelPath, 'w')
    
    # Add only the fields that we require
    f.write('/*Fake label file for input file: ' + jp2Path + '*/')
    
    f.write('MAP_PROJECTION_TYPE          = "' + infoDict['Projection'] + '"')
    
    # We could get these but it is not worth the effort =)
    f.write('MAXIMUM_LATITUDE             = 1.000000000000 <DEG>')
    f.write('MINIMUM_LATITUDE             = 0.000000000000 <DEG>')
    f.write('EASTERNMOST_LONGITUDE        = 1.000000000000 <DEG>')
    f.write('WESTERNMOST_LONGITUDE        = 0.000000000000 <DEG>')
    
    # This one we can't get
    f.write('STOP_TIME                    = 2008-08-24T08:08:08.000')
    
    f.close()
    
    return labelPath



def fetchAndPrepFile(setName, subtype, remoteURL, workDir):
    '''Retrieves a remote file and prepares it for upload'''
    
    #print 'Uploading file ' + setName
  
    if subtype != 'DEM': # Handles RED and COLOR images
        # The label file URL is the same as the image but with a different extension
        remoteLabelURL = getLabelPathFromImagePath(remoteURL)
    
        localFilePath  = os.path.join(workDir, os.path.basename(remoteURL))
        localLabelPath = os.path.join(workDir, os.path.basename(remoteLabelURL))

        # Retrieve the header file
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

        # Retrieve the image file
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

        # Check if there is geo data in the JP2 file
        jp2HasGeoData = IrgGeoFunctions.doesImageHaveGeoData(localFilePath)

        # If there is no label file, try to generate an artificial one.
        fakeLabelFile = False
        if not IrgFileFunctions.fileIsNonZero(localLabelPath):
            #raise Exception('Unable to download from URL: ' + remoteLabelURL)
            print 'WARNING: Unable to download from URL: ' + remoteLabelURL
            if not jp2HasGeoData:
                raise Exception('No geo data in JP2 and no Label file!  Cannot handle this!')
            print 'Generating a fake LBL file to proceed...'
            localLabelPath = writeFakeLabelFromJp2(localFilePath)
            fakeLabelFile  = True
        # At this point we always have a label file but it may be fake

        # Some images are missing geo information but we can add it from the label file
        if not jp2HasGeoData:
            # Correct the local file, then remove the old (bad) file
            tempPath = localFilePath + '.corrected'
            jp2ToImgAndBack(localFilePath, localLabelPath, tempPath)
            print 'Deleting JP2 file without metadata!'
            os.remove(localFilePath)
            os.mv(tempPath, localFilePath)
    
        # Call code to fix the header information in the JP2 file!
        cmd = FIX_JP2_TOOL +' '+ localFilePath
        print cmd
        os.system(cmd)


        # Check the projection type
        projType        = IrgGeoFunctions.getProjectionFromIsisLabel(localLabelPath)
        (width, height) = IrgIsisFunctions.getImageSize(localFilePath)
        if (projType == 'POLAR STEREOGRAPHIC') and (width < MAX_POLAR_UPLOAD_WIDTH):
            # Google has trouble digesting these files so handle them differently.

            #os.remove(localLabelPath)
            #raise Exception('POLAR STEREOGRAPHIC images on hold until Google fixes a bug!')
            print 'Special handling for POLAR STEROGRAPHIC image!'

            if fakeLabelFile:
                print 'Cannot reprocess polar image without a label file!'
                print 'All we can do is upload the file and hope for the best.'
                # First file is for upload, second contains the timestamp.
                return [localFilePath, localLabelPath]
            
            # Compute how many chunks are needed for this image
            numChunks = ceil(width / POLAR_WIDTH_CHUNK_SIZE)
            
            # Determine which chunk this DB entry is for
            chunkNum = getChunkNum(setName)
            print 'This is chunk number ' + str(chunkNum)
            
            if chunkNum >= numChunks: # Check for chunk number error
                raise Exception('Illegal chunk number: ' + setName)
                
            # If this is the main DB entry, we need to make sure the other DB entries exist!
            if chunkNum == 0:
                # Go ahead and try to add each chunk, the call will only go through if it does not already exist.
                for i in range(1,numChunks):
                    chunkSetName = makeChunkSetName(setName, i)
                    print 'Add chunk set name to DB: ' + chunkSetName
                    #common.addDataRecord(db, common.SENSOR_TYPE_HiRISE, subtype, chunkSetName, remoteURL)
                    
            raise Exception('DEBUG')
            
            # Now actually generate the desired chunk
            # - Need to use PIRL tools to extract a chunk to an IMG format, then convert that back to JP2 so that Google can read it.
            fileBasePath    = os.path.splitext(localFilePath)[0]
            localImgPath    = fileBasePath + '.IMG'
            localChunkPath  = fileBasePath + '_' + str(chunkNum) + '.JP2'
            chunkAreaString = getChunkAreaString(width, height, chunkNum)
            jp2ToImgAndBack(localFilePath, localLabelPath, localChunkPath, chunkAreaString=None)           
            
            # Just use the same label file, we don't care if the DB has per-chunk boundaries.
            return [localChunkPath, localLabelPath]
            
        else: # A normal, non-polar file.
            # First file is for upload, second contains the timestamp.
            return [localFilePath, localLabelPath]

    
    # TODO: Handle POLAR DEMS

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
