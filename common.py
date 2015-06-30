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

import os, glob, optparse, re, shutil, subprocess, string, time, urllib, urllib2
import datetime

'''
    Common definitions and utilities
'''


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
    def minLon(self):
        return self.data[10]
    def maxLon(self):
        return self.data[11]
    def minLat(self):
        return self.data[12]
    def maxLat(self):
        return self.data[13]

    def bbString(self):
        s = (str(self.minLon()) + ', ' + str(self.maxLon()) + ', ' + 
             str(self.minLat()) + ', ' + str(self.maxLat()))
        return s

    def __str__(self):
        s = ''
        for d in self.data:
            s += str(d)

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

#def fetchAndPrepFile(db, setName, subtype, remoteURL, workDir):
#    '''Retrieves a remote file and prepares it for upload'''
#    returns a list of created files with the first one being the one to upload

#def getUploadList(fileList):
#    '''Returns the subset of the fileList that needs to be uploaded to GME'''
#    return fileList[0]


#--------------------------------------------------------------------------------
# Function to call to add an entry to the database

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


def removeDataRecord(db, sensor, subType, setName):
    '''Removes a sensor data entry from the database'''

    print 'DELETING RECORD ' + setName + ',  ' + subType
    cursor = db.cursor()
    cursor.execute("DELETE FROM Files WHERE sensor=? AND subtype=? AND setname=?",
                      (str(sensor), subType, setName))
    db.commit()
    return True


#--------------------------------------------------------------------------------
# Miscellaneous functions

def isFileThisOld(filePath, numDays=0, numHours=0):
    '''Return true if a file is old'''
    fileTime = datetime.datetime.fromtimestamp(os.stat(filePath).st_ctime)
    today    = datetime.datetime.now()
    return (today > (fileTime + datetime.timedelta(days=numDays, hours=numHours)) )






