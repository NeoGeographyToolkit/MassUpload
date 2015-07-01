
import os
import sys
import time
import common
import logging
import sqlite3
#from pysqlite2 import dbapi2 as sqlite3
    
import MosaicUtilities

# Set the size limit of our data set cache
MAX_STORED_DATA_SETS = 20
    


# This function has to be outside the class so it can work with multiprocessing
def downloadHrscFile(remoteURL, localFilePath):
    '''Download a single HRSC file and convert to TIFF format'''
    
    force = False

    # Download the file
    downloadPath = localFilePath + '_download.IMG'
    if force or not os.path.exists(localFilePath):
      cmd = 'wget ' + remoteURL + ' -O ' + downloadPath
      MosaicUtilities.cmdRunner(cmd, downloadPath, force)

    # Convert to GTiff format
    cmd = 'gdal_translate -of GTiff ' + downloadPath + ' ' + localFilePath
    MosaicUtilities.cmdRunner(cmd, localFilePath, force)
    
    # Clean up the download file
    if os.path.exists(downloadPath):
        os.remove(downloadPath)


class BadHrscFileChecker():
    '''Simple class to monitor a list of known bad HRSC files'''

    def __init__(self, csvPath, logger):
        '''Init with list of bad files'''
        # Load the sets into an internal list
        self._badList = []
        with open(csvPath, 'r') as handle:
            for line in handle:
                name = line.strip()
                if '_' not in line: # Lines may not include the 'default' second part
                    name += '_0000' 
                if 'h' not in line: # Lines may not include the leading 'h'
                    name = 'h' + name
                self._badList.append(name)

        logger.info('Loaded bad HRSC list: \n' + str(self._badList))

    def isSetBad(self, setName):
        '''Returns True if this is a known bad data set'''
        return (setName in self._badList)



class HrscFileCacher():
    '''Class to manage which HRSC images are stored locally on disk'''
    
    
    def __init__(self, dbPath, outputFolder, badHrscFilePath, pool=None):
        '''Load the DB from disk'''

        # The input dictionary contains some channels we don't need for the basemap
        self._NEEDED_CHANNELS = ['nd3', 're3', 'bl3', 'gr3', 'ir3']
        
        self._logger = logging.getLogger('hrscFileCacher')
        self._logger.info('Initializing hrscFileCacher')

        if not os.path.exists(dbPath):
            raise Exception('ERROR: SQL database does not exist!')
        self._db = sqlite3.connect(dbPath)
        
        self._pool = pool

        self._badChecker = BadHrscFileChecker(badHrscFilePath, self._logger)

        # Create the output folder
        self._outputFolder = outputFolder
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder)
        
        # Create a working directory
        self._workDir = os.path.join(outputFolder, 'working_directory')
        if not os.path.exists(self._workDir):
            os.mkdir(self._workDir)
        
        self._cachedDataSets     = [] # Fully downloaded data sets
        self._incompleteDataSets = [] # List of partially downloaded data sets
        
        # Load list of cached sets from disk.
        self._findCachedFiles()
        

    def __del__(self):
        '''Close the open DB'''
        self._logger.info('Cleaning up hrscFileCacher')
        self._db.close()
    
    def _findCachedFiles(self):
        '''Identify all the cached files currently on disk'''
    
        # Just search the output directory for folders and each one is a data set.
        # - This will fail if any other junk gets in the output folder
        itemsInOutputDir = os.listdir(self._outputFolder)
        currentTime = time.time()
        for f in itemsInOutputDir:
            if (len(f) > 3) and (f != 'working_directory'): # Skip junk and work directory
                setFolder = os.path.join(self._outputFolder, f)
                if self._checkIfSetIsComplete(f):
                    self._cachedDataSets.append( (f, currentTime) )
                    self._logger.info('hrscFileCacher: Found existing cached file ' + f)
                else:
                    self._incompleteDataSets.append(setFolder)
                    self._logger.warning('Incomplete data set found: ' + setFolder)
    
    def _checkIfSetIsComplete(self, setName):
        setFolder = self._getStorageFolder(setName)
        for c in self._NEEDED_CHANNELS:
            fileName = self._makeFileName(setName, c)
            filePath = os.path.join(setFolder, fileName)
            if not os.path.exists(filePath):
                return False # A file is missing!
            #else:
            #    print 'Confirmed file ' + filePath
        return True # All files found

    def getHrscSetList(self, lonlatRect=None):
        '''Retrieve a list of HRSC data sets from the local database'''
        
        hrscSetList = []
        cursor = self._db.cursor()
    
        # Only need to retrieve the nadir files here
        query = ('SELECT * FROM Files WHERE sensor=%d AND subtype="nd3" AND status=%d' %
                      (common.SENSOR_TYPE_HRSC, common.STATUS_CONFIRMED))
        if lonlatRect != None: # Incorporate a bounding box
            query += (' AND minLat<%d AND maxLat>%d AND minLon<%d AND maxLon>%d' %
                      (lonlatRect.maxY, lonlatRect.minY, lonlatRect.maxX, lonlatRect.minX))
        query += ' AND (maxLon - minLon) < 10.0' # FOR NOW: Skip wraparound images.  Later we need to fix!
        #query += ' LIMIT 5'  # DEBUG!!!
        cursor.execute(query)
        rows = cursor.fetchall()
        
        if rows == []: # Make sure we found the next lines
            raise Exception('Could not find any HRSC files!')
    
        # Just extract the set name
        for line in rows:
            fileInfo = common.TableRecord(line) # Wrap the data line
            setName = fileInfo.setName()[:-4] # Strip off '_nd3' from the end

            if not self._badChecker.isSetBad(setName):
                hrscSetList.append(setName)
                #print fileInfo.setName() + ' --> ' + fileInfo.bbString()
            else:
                self._logger.info('Skipping known bad set ' + setName)
            
        return hrscSetList
    
    def fetchHrscDataSet(self, setName):
        '''Downloads all the files for a data set (if needed) and returns paths to them.'''
        
        # Look up all the URLs
        urlDict  = self._getUrlDictForSet(setName)
        # Download the files if needed
        pathDict = self._retrieveDataSetForHrscMap(urlDict)
        return pathDict
    
    
    def _getUrlDictForSet(self, setName):
        '''Return the URLs for all HRSC files in a single set'''

        cursor = self._db.cursor()
    
        # Only need to retrieve the nadir files here
        setMatchString = '"' + setName + '%"'
        searchString = ('SELECT * FROM Files WHERE sensor=%d AND setname LIKE %s AND status=%d' %
                            (common.SENSOR_TYPE_HRSC, setMatchString, common.STATUS_CONFIRMED))
        cursor.execute(searchString)
        rows = cursor.fetchall()
      
        if rows == []: # Make sure we found the next lines
            self._logger.info(searchString + '\n' + str(rows))
            raise Exception('Could not find any data files for data set: ' + setName)
   
        # Build a dict containing the URL for each channel
        dataDict = {}
        dataDict['setName'] = setName
        for line in rows:
            fileInfo = common.TableRecord(line) # Wrap the data line
            dataDict[fileInfo.subtype()] = fileInfo.remoteURL()
        
        return dataDict
    
    def _getStorageFolder(self, setName):
        '''Get the location where this data set will be stored'''
        return os.path.join(self._outputFolder, setName)

    def _deleteDataSet(self, folder):
        self._logger.info('Deleting data folder: ' + folder)
        cmd = 'rm -rf ' + folder
        os.system(cmd)

    def _makeRoomForNewDataSet(self):
        '''Delete a cached data set to make room for a new one'''
        
        if len(self._cachedDataSets) < MAX_STORED_DATA_SETS:
            return # We still have room, don't delete anything.
        
        # Clear out any incomplete data sets
        for setFolder in self._incompleteDataSets:
            self._deleteDataSet(setFolder)
        self._incompleteDataSets = []

        # Find the last accessed data set
        oldestIndex = -1
        oldestTime  = 0
        currentTime = time.time()
        for i in range(0,len(self._cachedDataSets)):
            pair = self._cachedDataSets[i]
            timeDiff = currentTime - pair[1]
            if timeDiff > oldestTime:
                oldestTime  = timeDiff
                oldestIndex = i
        # Remove the data set from the internal record
        oldestPair = self._cachedDataSets.pop(i)
        oldestSet  = oldestPair[1]
        
        # Delete it
        setFolder = self._getStorageFolder(oldestSet)
        self._deleteDataSet(setFolder)
        
    
    def _downloadHrscFile(self, remoteURL, localFilePath):
        '''Download a single HRSC file and convert to TIFF format'''
        
        # Download the file
        downloadPath = os.path.join(self._workDir, 'download.IMG')
        cmd = 'wget ' + remoteURL + ' -O ' + downloadPath
        MosaicUtilities.cmdRunner(cmd, downloadPath)
    
        # Convert to GTiff format
        cmd = 'gdal_translate -of GTiff ' + downloadPath + ' ' + localFilePath
        MosaicUtilities.cmdRunner(cmd, localFilePathPath)
        
        # Clean up the download file
        if os.path.exists(downloadPath):
            os.remove(downloadPath)
    

    def _makeFileName(self, setName, channel):
        return setName + '_' + channel + '.tif'

    def _retrieveDataSetForHrscMap(self, hrscDataDict):
        '''Retrieves the HRSC files needed to build the HRSC basemap'''
    
        setName = hrscDataDict['setName']
    
        if setName not in self._cachedDataSets: # If we don't currently have the data available

            # Make sure all the required channels exist in the database
            for c in self._NEEDED_CHANNELS:
                if not c in hrscDataDict:
                    raise Exception('ERROR: Missing channel ' + c + ' in data set ' + setName)
                
            # Make sure we have enough room to store the new set
            self._makeRoomForNewDataSet()
                
            # Create the output folder
            setOutputFolder = self._getStorageFolder(setName)
            if not os.path.exists(setOutputFolder):
                os.mkdir(setOutputFolder)
            
            # Internally record that we have the data set
            self._cachedDataSets.append( (setName, time.time()) )
            
            # Download each of the files
            poolResults = []
            for c in self._NEEDED_CHANNELS:
                url        = hrscDataDict[c]
                fileName   = self._makeFileName(setName, c)
                outputPath = os.path.join(setOutputFolder, fileName)
                
                if self._pool:
                    poolResults.append(self._pool.apply_async(downloadHrscFile, 
                                                              args=(url, outputPath)))
                else: # Just call the function
                    downloadHrscFile(url, outputPath)

            if self._pool: # Wait for all the tasks to complete
                print 'Waiting for HRSC download processes to complete...'
                for result in poolResults:
                    result.get()
                print 'All HRSC download processes complete for set ' + setName

        # Now we have this HRSC image locally.

        # Build a list of the files in the output folder in an output dictionary
        outputDict = {'setName': setName}
        allChannelList = []
        for c in self._NEEDED_CHANNELS:
            fileName   = self._makeFileName(setName, c)
            outputPath = os.path.join(setOutputFolder, fileName)
            outputDict[c] = outputPath
            allChannelList.append(outputPath)
        outputDict['allChannelPaths'] = allChannelList # Include a list containing all the channels
    
        return outputDict

    
    
    #def queryHrscFilesInLocation():
        
    #    select setname from Files where sensor=1 and minLon<102.64 and maxLon>98.2023 and minLat<-18.6077 and maxLat>-50.1955 and subtype="nd3";    

