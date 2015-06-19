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

import os, glob, optparse, re, shutil, subprocess, sys, string, simplekml

import MosaicUtilities

''' This tool takes a grid of image tiles that form the base of an image pyramid and
    constructs the lower resolution image tiles on top of them.
'''
# TODO: Move this up to the Tools directory    
    
    
# Convention: Level 0 is the base level where the tiles start.
    
def makeTileName(tileIndex, isImage=True):
    '''Method of generating the tile name'''
    if isImage:
        return ("tile_%d_%d.tif" % (tileIndex.row, tileIndex.col))
    else:
        return ("tile_%d_%d.kml" % (tileIndex.row, tileIndex.col))
    
def getLevelFolder(outputFolder, level):
    '''Get the folder where a level is stored'''
    return os.path.join(outputFolder, str(level))


class TileLevel:
    
    def __init__(self, outputFolder, level, levelZeroHeight, levelZeroWidth):
        self.outputFolder = outputFolder
        self.level        = level
    
        # TODO: Accept flexible boundaries
        self._bounds = MosaicUtilities.Rectangle(0, 360, -90, 90)
    
        # Compute the number of tiles we need on this level
        
    
    def getTilePath(self, level, tileIndex, isImage=True):
        '''Returns the full path to a tile'''
        folder   = getLevelFolder(self._outputFolder, level)
        filename = makeTileName(tileIndex, isImage)
        return os.path.join(folder, filename)
    
    def getInputTiles(self, tileIndex):
        '''Given a tile index, return the indices of the source tiles'''
        indexList = []
        indexList.append( MosaicUtilities.TileIndex(tileIndex.row*2+0, tileIndex.col*2+0) )
        indexList.append( MosaicUtilities.TileIndex(tileIndex.row*2+0, tileIndex.col*2+1) )
        indexList.append( MosaicUtilities.TileIndex(tileIndex.row*2+1, tileIndex.col*2+0) )
        indexList.append( MosaicUtilities.TileIndex(tileIndex.row*2+1, tileIndex.col*2+1) )
        return indexList
    
    def getInputTileString(self, tileIndex):
        '''Make an input string listing all the input tiles to a new tile'''
        if self.level == 0: # Special handling for the base level
            return None
        
        # Get the indices of all the input tiles
        indexList  = self.getInputTiles(tileIndex)
        inputLevel = self.level - 1
        
        # Build a string containing the paths to each of the input tiles
        tileString = ''
        for index in indexList:
            tilePath = getTilePath(inputLevel, index, True) # Only need the image paths here
            tileString += tilePath + ' '
        
    def makeTile(self, tileIndex):
        '''Create an input tile'''
    
        # Create the image file
        inputTileTiffString = self.getInputTileString(tileIndex)
        outputTilePath  = self.getTilePath(self.level, tileIndex, True)
        cmd = 'montage -resize 50% -background black -mode Concatenate -tile 2x2 ' + inputTileString +' '+ outputTilePath
        print cmd
        os.system(cmd)
        
        # Create the KML file
        inputTileKmlString = inputTileTiffString.replace('.tif', '.kml')
        
    def getLatLonAltBox(self, level, tileIndex):
        '''Gets the ROI for one tile'''
        return simplekml.LatLonAltBox(north=0, south=0, east=0, west=0)
    
    def makeKmlRegion(self, level, tileIndex):
        '''Sets up the KML region object for a single tile'''
        simplekml.Region(latlonaltbox=self.getLatLonAltBox(level, tile),
                         lod=simplekml.Lod(minlodpixels=128, maxlodpixels=-1)) # TODO: Set this?
        
    def makeKmlFile(self, tileIndex):
        '''Creates the KML file for a given image tile'''
        kml     = simplekml.Kml()
        kml.doc = simplekml.Document()
        # Region for this file
        kml.doc.region = self.makeKmlRegion(self.level, tileIndex)
        
        # Network links for lower level files
        if level > 0:
            levelDown = level - 1
            indexList = self.getInputTiles(tileIndex)
            for index in indexList:
                netlink = kml.newnetworklink(name=index.getPostFix(),
                                             region=self.makeKmlRegion(levelDown, index))
                netLink.link.href = self.getTilePath(levelDown, index, False)
                netLink.link.viewrefreshmode = simplekml.ViewRefreshMode.onrequest
        
        # Ground overlay
        tileImagePath = self.getTilePath(self.level, tileIndex, True)
        tileKmlPath   = self.getTilePath(self.level, tileIndex, False)
        kml.doc.groundOverlay = simplekml.GroundOverlay(icon=simplekml.Icon(href=tileImagePath),
                                                        latlonaltbox=self.getLatLonAltBox(level, tile))
        # Save the completed file
        kml.save(tileKmlPath)

    # TODO: Multiprocess function to generate all of the tiles in one level!

def main():

    outputPath = ''


    try:
        usage = "usage: geoTiffTool.py [--help][--manual] geotiffPath\n"# [lat] [lon] [lat2] [lon2]\n  "
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("--manual", action="callback", callback=man,
                          help="Read the manual.")

        ## Many GeoTiff images are stored in projected coordinates which makes this tricky
        #parser.add_option("--crop", action="store_true", dest="doCrop", default=False,
        #                            help="Crops the image between the four lat/lon arguments.")

        parser.add_option("--normalize-eqc-lon", action="store_true", dest="normEqcLon", default=False,
                                    help="Normalizes an EQC image to lie between +/- 180 degrees.")

        # Add other features as the need arises

        #parser.add_option("-o", "--output-path", dest="outputPath",
        #                        help="Output path (default replace extension with .kml")
        
        (options, args) = parser.parse_args()

        # Parse positional arguments
        if not args: parser.error("need input cube file")           
        tiffPath = args[0]

    except optparse.OptionError, msg:
        raise Usage(msg)


    # First make a copy of each input tile in the output folder to form the base layer,
    #  but resize each of the tiles to 512 pixels to make the pyramid creation easier.



    return 0



if __name__ == "__main__":
    sys.exit(main())
