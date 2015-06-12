
import os
import sys

import common
    
# Set the size limit of our data set cache
MAX_STORED_DATA_SETS = 10
    

class HrscFileCacher():
    '''Class to manage which HRSC images are stored locally on disk'''
    
    
    def __init__(self, dbPath, outputFolder):
        '''Load the DB from disk'''
        
        if not os.path.exists(dbPath):
            raise Exception('ERROR: SQL database does not exist!')
        self._db = sqlite3.connect(dbPath)
        
        # Create the output folder
        self._outputFolder = outputFolder
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder)
        
        # Create a working directory
        self._workDir = os.path.join(outputFolder, 'working_directory')
        
        self._cachedDataSets = []
        
        # Load list of cached sets from disk.
        self._findCachedFiles()
        

    def __del__(self):
        '''Close the open DB'''
        db.close()
    
    def _findCachedFiles(self):
        '''Identify all the cached files currently on disk'''
    
        # Just search the output directory for folders and each one is a data set.
        # - This will fail if any other junk gets in the output folder
        itemsInOutputDir = os.listdir(self._outputFolder)
        currentTime = time.time()
        for f in itemsInOutputDir:
            if (len(f) > 3) and (f != 'working_directory'): # Skip junk and work directory
                self._cachedDataSets.append( (f, currentTime) )
    
    # TODO: Try this out on Lunokhod2!!
    def getAllHrscSetList():
        '''Retrieve a list of ALL the HRSC data sets from the local database'''
        
        hrscSetList = []
        cursor = self._db.cursor()
    
        # Only need to retrieve the nadir files here
        cursor.execute('SELECT * FROM Files WHERE sensor=? AND subtype=nd3 AND status=? LIMIT 5',
                       (str(common.SENSOR_TYPE_HRSC), str(common.STATUS_CONFIRMED)))
        rows = cursor.fetchall()
        
        if rows == []: # Make sure we found the next lines
            raise Exception('Could not find any HRSC files!')
    
        # Just extract the set name
        for line in rows:
            fileInfo = common.TableRecord(line) # Wrap the data line
            hrscSetList.append(setName)
            
        print hrscSetList
        
        raise Exception('DEBUG')
    
    def fetchHrscDataSet(self, setName):
        '''Downloads all the files for a data set (if needed) and returns paths to them.'''
        
        # Look up all the URLs
        urlDict  = self._getUrlDictForSet(setName)
        # Download the files if needed
        pathDict = self._retrieveDataSetForHrscMap(urlDict)
    
    
    def _getUrlDictForSet(self, setName):
        '''Return the URLs for all HRSC files in a single set'''

        cursor = self._db.cursor()
    
        # Only need to retrieve the nadir files here
        cursor.execute('SELECT * FROM Files WHERE sensor=? AND setname=? AND status=?',
                       (str(common.SENSOR_TYPE_HRSC), setName, str(common.STATUS_CONFIRMED)))
        rows = cursor.fetchall()
        
        if rows == []: # Make sure we found the next lines
            raise Exception('Could not find any data files for data set: ' + setName)
    
        # Build a dict containing the URL for each channel
        dataDict = {}
        dataDict['setName'] = setName
        for line in rows:
            fileInfo = common.TableRecord(line) # Wrap the data line
            dataDict[fileInfo.subtype] = fileInfo.remoteURL
            
        return dataDict
    
    def _getStorageFolder(self, setName):
        '''Get the location where this data set will be stored'''
        return os.path.join(self._outputFolder, setName)

    def _makeRoomForNewDataSet(self):
        '''Delete a cached data set to make room for a new one'''
        
        if len(self._cachedDataSets) < MAX_STORED_DATA_SETS:
            return # We still have room, don't delete anything.
        
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
        print 'Deleting cached data set ' + oldestSet
        setFolder = self._getStorageFolder(oldestSet)
        cmd = 'rm -rf ' + setFolder
        print cmd
        os.system(cmd)
        
    
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
    

    def _retrieveDataSetForHrscMap(self, hrscDataDict):
        '''Retrieves the HRSC files needed to build the HRSC basemap'''
    
        setName = hrscDataDict['setName']

        NEEDED_CHANNELS = ['nd3', 're3', 'bl3', 'gr3', 'nir'] # TODO: Check these!
    
        if setName not in self._cachedDataSets: # If we don't currently have the data available
            
            # Make sure all the required channels exist in the database
            for c in NEEDED_CHANNELS:
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
            for c in NEEDED_CHANNELS:
                url        = hrscDataDict[c]
                fileName   = setName + '_' + c
                outputPath = os.path.join(setOutputFolder, fileName)
                self._downloadHrscFile(url, outputPath)
        # Now we have this HRSC image locally.

        # Build a list of the files in the output folder in an output dictionary
        outputDict = {'setName': setName}
        allChannelList = []
        for c in NEEDED_CHANNELS:
            fileName   = setName + '_' + c
            outputPath = os.path.join(setOutputFolder, fileName)
            outputDict[c] = outputPath
            allChannelList.append(outputPath)
        outputDict['allChannelPaths'] = allChannelList # Include a list containing all the channels
    
        return outputDict

    
    
    #def queryHrscFilesInLocation():
        
    #    select setname from Files where sensor=1 and minLon<102.64 and maxLon>98.2023 and minLat<-18.6077 and maxLat>-50.1955 and subtype="nd3";    

