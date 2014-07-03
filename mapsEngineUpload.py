#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

import sys, os, glob, optparse, re, shutil, subprocess, string, time

import json, urllib2, requests, argparse

import apiclient
from apiclient import discovery
import httplib2
from oauth2client import client
from oauth2client import file as oauth2client_file
from oauth2client import tools

# Authorization codes
# TODO: Load these!
API_KEY       = 'AIzaSyAM1ytSqkzubDMzjVWBjM19uawCkIBVvLY'
CLIENT_ID     = '298099604529-69gprkqj67qkm5ncfik32uenug8qgagn.apps.googleusercontent.com'
CLIENT_SECRET = 'kNDmMQi_BH2ttN3XRIY2GA-7'
PROJECT_ID    = '04070367133797133737'

SENSOR_TYPE_HiRISE = 0
SENSOR_TYPE_HRSC   = 1
SENSOR_TYPE_CTX    = 2


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
Tool for uploading raster images to Google Maps Engine
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

# TODO: Is there a way to check for a file without knowing the asset ID?
def checkIfFileIsLoaded(bearerToken, assetId = '04070367133797133737-15079892155897256865'):
    """Determine if a file has already been uploaded into Maps Engine"""

    # Send request for information on this asset
    url = 'https://www.googleapis.com/mapsengine/v1/rasters/' + assetId
    tokenString = 'Bearer '+bearerToken
    headers = {'Authorization': tokenString}
    response = requests.get(url, headers=headers)
    
    # See if we got a response
    if response.status_code != 200:
        return (False, response.status_code)
    
    # TODO: More accurate check?
    # Check if all files are loaded and processing is finished    
    jsonDict = json.loads(response.text)
    status = True
    for f in jsonDict['files']:
        if (f['uploadStatus'] != 'complete'):
            status = False
    if jsonDict['processingStatus'] != 'complete':
        status = False

    print response.text
    return (status, response.status_code)
    
    #url = 'https://www.googleapis.com/mapsengine/v1/assets?projectId='+PROJECT_ID
    #tokenString = 'Bearer '+bearerToken
    #headers = {'Authorization': tokenString
    #response = requests.get(url, headers=headers)
    #print response.text
    # ?
    # raster collection ID = GET /rasterCollections?projectId=<id>
    #   raster list = GET /rasterCollections/#rasterCollectionId$/rasters
    
#--------------------------------------------------------------------------------

# TODO: Need to update these to use the refresh token for "offline" access!

def getCredentials(redo=False):

    # Parse OAuth command line arguments
    parser = argparse.ArgumentParser(parents=[tools.argparser])
    #flags = parser.parse_args()
    flags, unknown = parser.parse_known_args()
    
    print 'Looking for cached credentials'
    # Check if we already have a file with OAuth credentials
    storage = oauth2client_file.Storage('mapsengine.dat')
    credentials = storage.get()
    
    
    if credentials is None or credentials.invalid or redo:
        print 'Getting credentials'
        # Start local server, redirect user to authentication page, receive OAuth
        # credentials on the local server, and store credentials in file
        flow = client.OAuth2WebServerFlow(
                            client_id=CLIENT_ID,
                            client_secret=CLIENT_SECRET,
                            response_type='code',
                            scope='https://www.googleapis.com/auth/mapsengine',
                            user_agent='Google-MapsEngineApiSample/1.0',
                            redirect_uri='urn:ietf:wg:oauth:2.0:oob',
                            access_type='offline',
                            approval_prompt='auto') # Change to 'auto'
        credentials = tools.run_flow(flow, storage, flags)
      
    return credentials

def authorize(redo=False):

    credentials = getCredentials(redo)    
    
    #jsonCredentials = credentials.to_json()
    #jsonDict = json.loads(jsonCredentials)
    #for key, value in jsonDict.iteritems():
    #    print key + ' ---> ' + str(value)
    
    # Set up discovery with authorized credentials
    http = credentials.authorize(httplib2.Http())
    
    service = apiclient.discovery.build('mapsengine', 'v1', http=http, developerKey=API_KEY)

    # It is not clear why but adding this pointless code fixed the authorization problems!
    #print '=================='    
    # Read the location of every Feature in a Table.
    features = service.tables().features()
    request = features.list(id='12421761926155747447-06672618218968397709',
                            maxResults=10, version='published')
    response = request.execute()
    #for feature in response['features']:
    #  print feature['geometry']['coordinates']
    #print '---------------------'

    jsonCredentials = credentials.to_json()
    jsonDict = json.loads(jsonCredentials)
    #for key, value in jsonDict.iteritems():
    #    print key + ' ---> ' + str(value)
    
    bearerToken = jsonDict['access_token']
    #print bearerToken
    return bearerToken
    

## Read the location of every Feature in draft version of Table.
#features = service.tables().features()
#request = features.list(id=TABLE_ID)
#while request is not None:
#  resource = request.execute()
#  for feature in resource['features']:
#    print feature['geometry']['coordinates']
#
#  # Is there an additional page of features to load?
#  request = features.list_next(request, resource)

def createRasterAsset(bearerToken, inputFile, sensorType):
    
    url = 'https://www.googleapis.com/mapsengine/v1/rasters/upload'
    
    justFilename = os.path.basename(inputFile)
    
    if sensorType == SENSOR_TYPE_HiRISE:
        data = ( 
        {
          "projectId": PROJECT_ID,  # REQUIRED, taken from Maps Engine URL
          "name": 'HiRISE_'+justFilename,  # REQUIRED
          "description": "HiRISE map projected RDR data",
          "files": [ # REQUIRED
            { "filename": justFilename }
          ],
          #"acquisitionTime": {
          #  "start": "2010-01-01T12:00:00Z",
          #  "end": "2010-12-01T12:00:00Z",
          #  "precision": "second"
          #},
          "draftAccessList": "Map Editors", # REQUIRED
          "attribution": "NASA Public Domain", # REQUIRED
          "tags": ["Mars", "MRO", "HiRISE"],
          "maskType": "autoMask"
        } )
    elif sensorType == SENSOR_TYPE_HRSC:
        data = ( 
        {
          "projectId": PROJECT_ID,  # REQUIRED, taken from Maps Engine URL
          "name": 'HRSC_'+justFilename,  # REQUIRED
          "description": "HRSC map projected RDR data",
          "files": [ # REQUIRED
            { "filename": justFilename }
          ],
          #"acquisitionTime": {
          #  "start": "2010-01-01T12:00:00Z",
          #  "end": "2010-12-01T12:00:00Z",
          #  "precision": "second"
          #},
          "draftAccessList": "Map Editors", # REQUIRED
          "attribution": "NASA Public Domain", # REQUIRED
          "tags": ["Mars", "MEX", "HRSC"],
          "maskType": "autoMask"
        } )
    elif sensorType == SENSOR_TYPE_CTX:
        data = ( 
        {
          "projectId": PROJECT_ID,  # REQUIRED, taken from Maps Engine URL
          "name": 'CTX_'+justFilename,  # REQUIRED
          "description": "CTX map projected RDR data",
          "files": [ # REQUIRED
            { "filename": justFilename }
          ],
          #"acquisitionTime": {
          #  "start": "2010-01-01T12:00:00Z",
          #  "end": "2010-12-01T12:00:00Z",
          #  "precision": "second"
          #},
          "draftAccessList": "Map Editors", # REQUIRED
          "attribution": "NASA Public Domain", # REQUIRED
          "tags": ["Mars", "MRO", "CTX"],
          "maskType": "autoMask"
        } )
    else:
        raise Exception('Unrecognized sensor type!')
        
    print(data)
    tokenString = 'Bearer '+bearerToken
    headers = {'Authorization': tokenString,
               'Content-Type': 'application/json'}
    print(headers)
    
    response = requests.post(url, data=json.dumps(data), headers=headers)
    
    print response.text
    
    print('Received status code ' + str(response.status_code))
    if response.status_code == 401:
        print('Error: Unauthorized access!')
    if response.status_code == 403:
        print('Error: Forbidden operation!')
    if response.status_code != 200:
        return (False, response.status_code)
    
    jsonDict = json.loads(response.text)
    #for key, value in jsonDict.iteritems():
    #    print key + ' ---> ' + str(value)
       
    # Return the asset ID from the response
    return (True, jsonDict['id'])
    


def uploadFile(bearerToken, assetId, filename):

    # Check input image
    if not os.path.exists(filename):
        raise Exception('Input image file ' + filename + ' is missing!')
    imageSizeBytes = os.path.getsize(filename)
    
    # Get the file extension tag to use
    ext = os.path.splitext(filename)[1]
    if (ext.lower == '.jp2'):
        contentString = 'image/jpg2'
    else: # Default to geotiff type
        contentString = 'image/tiff'
    
    # Set up POST request
    justFilename = os.path.basename(filename)
    url = 'https://www.googleapis.com/upload/mapsengine/v1/rasters/'+str(assetId)+'/files?filename='+justFilename
    tokenString = 'Bearer '+bearerToken
    headers = {'Authorization':  tokenString,
               'Content-Type':   contentString,
               'Content-Length': str(imageSizeBytes)}

    print headers
    print url
    print filename

    # Submit the post request with file data
#    fileList = {'file': open(filename, 'rb')}
#    response = requests.post(url, headers=headers, files=fileList)
 
    with open(filename, 'rb') as f:
        response = requests.post(url, headers=headers, data=f) 
    
    print response.text
    
    # Check response status code
    print('Received status code ' + str(response.status_code))
    if response.status_code != 204:
        print 'Failed to upload file!'
        return False
    print 'File upload started successfully!'
    return True
 

    # To check upload progress:
    # GET https://www.googleapis.com/mapsengine/v1/rasters/{raster_ID}
    # Authorization: Bearer {token}


def main(argsIn):

    print ('#################################################################################')
    print ("Running mapsEngineUpload.py")

    #try:
    #try:
    usage = "usage: mapsEngineUpload.py <input image> [--manual]\n  "
    parser = optparse.OptionParser(usage=usage)

    parser.add_option("--sensor", type="int", dest="sensor", default=0,
                              help="Which sensor? (HiRISE=0, HRSC=1, CTX=2).")

    parser.add_option("--manual", action="callback", callback=man,
                      help="Read the manual.")
    (options, args) = parser.parse_args(argsIn)

    if len(args) < 1: # DEBUG
        #options.inputPath = 'means.png'
        options.inputPath = '/home/smcmich1/data/production/NAC_DTM_M151318807_M181974094/results/output-DEM.tif'
        #options.inputPath = '/home/smcmich1/data/google/mapsengine-cmd-line-sample/PSP_001427_1820_RED.JP2'
    else:
        options.inputPath = args[0]
    
    #except(optparse.OptionError, msg):
    #    raise Usage(msg)
    
    
    if not os.path.exists(options.inputPath):
        raise Exception('Input file does not exist!')

    startTime = time.time()

    # Get server authorization
    bearerToken = authorize()
    
    
    #checkIfFileIsLoaded(bearerToken)
    #raise Exception('----TEST----')
    
    
    # Create empty raster asset request
    (success, assetId) = createRasterAsset(bearerToken, options.inputPath, options.sensor)
    #if not success:
    #    print 'Refreshing access token...'
    #    bearerToken = authorize(True)
    #    (success, assetId) = createRasterAsset(bearerToken, options.inputPath, options.sensor)
    if not success:
        raise Exception('Could not get access token!')
    
    if success:
        print 'Created asset ID ' + str(assetId)
    
        # Load a file associated with the asset
        uploadFile(bearerToken, assetId, options.inputPath)

    
    endTime = time.time()

    print("Finished in " + str(endTime - startTime) + " seconds.")
    print('#################################################################################')
    return assetId

    #except(Usage, err):
    #    print(err)
    #    #print(>>sys.stderr, err.msg)
    #    return 2

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
