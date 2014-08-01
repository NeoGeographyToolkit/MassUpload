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

import os, glob, optparse, re, shutil, subprocess, string, time


import IrgStringFunctions, IrgGeoFunctions


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, ''' Script for grabbing HiRISE data files'''

    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


#--------------------------------------------------------------------------------

    


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

    print "Started prepThemisMosaic.py"

    try:
        try:
            usage = "usage: prepThemisMosaic.py  [--help][--manual]\n  "
            parser = optparse.OptionParser(usage=usage)
            #parser.add_option("-u", "--upload", dest="upload", type=int,
            #                  help="Upload this many files instead of fetching the list.")
            #parser.add_option("-o", "--output-folder", dest="outputFolder",
            #                  help="Specifies the folder to copy the data to.",
            #                  default='./')
            #parser.add_option("-n", "--name", dest="name",
            #                  help="Only get the data for the DTM with this name.",
            #                  default='')
            
            #parser.add_option("--color", action="store_true", dest="getColor",
            #                  help="Retrieve COLOR image instead of RED images.")

            #parser.add_option("--threads", type="int", dest="numThreads", default=1,
            #                  help="Number of threads to use.")
                              
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            (options, args) = parser.parse_args()


        except optparse.OptionError, msg:
            raise Usage(msg)

        print "Beginning processing....."

        startTime = time.time()

        # For now everything is hard coded!
        sourceFolder = '/byss/mars/themis/nightir_100m'
        destFolder   = '/home/smcmich1/data/themisNight100m'
        
        # Force all images to have the same pixel resolution copied from lat00_lon000.tif
        horizontalPixelRes = ' 0.001687051876845 '
        verticalPixelRes   = ' -0.001687004442445 '
        
        
        # Loop through all files in the source folder
        for f in os.listdir(sourceFolder):
            # Skip everything but the .tif files
            if (not f[-4:] == '.tif'):
                continue
                
            
            inputPath   = os.path.join(sourceFolder, f)
            vrtPath     = os.path.join(destFolder, f+'.vrt')
            vrtEditPath = os.path.join(destFolder, f+'.edit.vrt')
            outputPath  = os.path.join(destFolder, f)
            
            # Skip files we already created
            if os.path.exists(outputPath):
                print 'Output file ' + outputPath + ' already exists.'
                cmd = '/home/smcmich1/repot/MassUpload/mapsEngineUpload.py --sensor 3 ' + outputPath
                print cmd
                os.system(cmd)
                continue
            
            # Build a vrt file in the output directory
            cmd = 'gdalbuildvrt ' + vrtPath +' '+ inputPath
            print cmd
            os.system(cmd)
                       
            # Replace the resolution fields in the vrt file
            oldVrtFile = open(vrtPath,     'r')
            newVrtFile = open(vrtEditPath, 'w')
            for line in oldVrtFile:

              if not ('GeoTransform' in line): # Just copy the other lines
                  newVrtFile.write(line)
                  continue
                  
              # Otherwise make the edit before writing the line
              parts   = line.split(',')
              newLine = (parts[0] +','+ horizontalPixelRes +','+ 
                         parts[2] +','+ parts[3] +','+ parts[4] +','+
                         verticalPixelRes + '</GeoTransform>\n')
              newVrtFile.write(newLine)
              #print line
              #print newLine

            # Done writing the new vrt file
            oldVrtFile.close()
            newVrtFile.close()

            # Convert edited VRT file to final output file
            cmd = 'gdal_translate -co compress=lzw ' + vrtEditPath + ' ' + outputPath;
            print cmd
            os.system(cmd)

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
