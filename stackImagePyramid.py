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

import os, glob, optparse, re, shutil, subprocess, sys, string
import math
import traceback
import simplekml
import MosaicUtilities

''' This tool takes a grid of image tiles that form the base of an image pyramid and
    constructs the lower resolution image tiles on top of them.
'''
# TODO: Move this up to the Tools directory    
    
    
# Convention: Level 0 is the base level where the tiles start.
    
def makeTileName(tileIndex, isImage=True):
    '''Method of generating the tile name'''
    if isImage:
        return ("tile_%04d_%04d.tif" % (tileIndex.row, tileIndex.col))
    else:
        return ("tile_%04d_%04d.kml" % (tileIndex.row, tileIndex.col))
    
def getLevelFolder(outputFolder, level):
    '''Get the folder where a level is stored'''
    return os.path.join(outputFolder, str(level))


class KmlTreeMaker:
    '''Class to build a Google Earth KML tree for visualizing the output'''
    
    def __init__(self, sourceFolder, outputFolder):
        
        self.sourceFolder = sourceFolder
        self.outputFolder = outputFolder
        self._tileSize    = 512 # Should be 256 or 512
        
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder)
        
        # TODO: Accept flexible boundaries
        self._bounds = MosaicUtilities.Rectangle(-180, 180, -90, 90)
        self._layerTilings = []
        self._addTileLayer(0) # Go ahead and set up the first tile layer
        # Other layers are added as they are requested
        
        
    def makeLevel(self, level):
        '''Make all the tiles for the given level.
           All previous levels must already exist.'''
           
        self._addTileLayer(level) # Make sure we are ready for this layer

        levelFolder = getLevelFolder(self.outputFolder, level)
        if not os.path.exists(levelFolder):
            os.mkdir(levelFolder)

        tileBounds = self._layerTilings[level].getTileIndexRect()
        numTiles   = tileBounds.area()
        print 'Making KML pyramid tiles for level ' + str(level)
        print ' - Found ' + str(numTiles) + ' tiles.'

        # TODO: Run this in parallel!
        force = False
        for tileIndex in tileBounds.indexGenerator():
            try:
                gotTile = self.makeTile(level, tileIndex, force)
            except Exception, e:
                print 'Caught exception processing tile ' + str(tileIndex)
                print str(e)
                traceback.print_exc(file=sys.stdout)
                
            #if gotTile:
            #    raise Exception('DEBUG')

    def _hasTileLayer(self, level):
        return len(self._layerTilings) > level
        
    def _addTileLayer(self, level):
        '''Sets up an internal object to handle on level's tile computations'''
        if self._hasTileLayer(level): # Already have it!
            return
        if len(self._layerTilings) != level: # Missing lower levels, recursively create them!
            self._addTileLayer(level - 1)
            
        levelScale = math.pow(2,level) # Each level halves the tiles in each direction
        
        # TODO: Improve this!
        FULL_BASEMAP_HEIGHT = 5760  # In pixels, low resolution.
        FULL_BASEMAP_WIDTH  = 11520       
        BASEMAP_TILE_HEIGHT = 45 # In pixels, chosen to divide evenly.
        BASEMAP_TILE_WIDTH  = 45
        numTileRows = (FULL_BASEMAP_HEIGHT / BASEMAP_TILE_HEIGHT) / levelScale
        numTileCols = (FULL_BASEMAP_WIDTH  / BASEMAP_TILE_WIDTH ) / levelScale
        tileWidth   = self._bounds.width()  / numTileCols
        tileHeight  = self._bounds.height() / numTileRows
        self._layerTilings.append(MosaicUtilities.Tiling(self._bounds, tileWidth, tileHeight, True))
        print self._layerTilings[level]
    
    
    def getTilePath(self, level, tileIndex, isImage=True):
        '''Returns the full path to a tile'''
        folder   = getLevelFolder(self.outputFolder, level)
        filename = makeTileName(tileIndex, isImage)
        return os.path.join(folder, filename)
    
    def _getInputTiles(self, tileIndex):
        '''Given a tile index, return the indices of the source tiles'''
        indexList = []
        indexList.append( MosaicUtilities.TileIndex(tileIndex.row*2+0, tileIndex.col*2+0) )
        indexList.append( MosaicUtilities.TileIndex(tileIndex.row*2+0, tileIndex.col*2+1) )
        indexList.append( MosaicUtilities.TileIndex(tileIndex.row*2+1, tileIndex.col*2+0) )
        indexList.append( MosaicUtilities.TileIndex(tileIndex.row*2+1, tileIndex.col*2+1) )
        return indexList
    
    def _getInputTileString(self, level, tileIndex):
        '''Make an input string listing all the input tiles to a new tile'''
        
        # Get the indices of all the input tiles
        indexList  = self._getInputTiles(tileIndex)
        inputLevel = level - 1
        
        # Build a string containing the paths to each of the input tiles
        tileString = ''
        for index in indexList:
            tilePath = getTilePath(inputLevel, index, True) # Only need the image paths here
            tileString += tilePath + ' '
        
    def makeTile(self, level, tileIndex, force=False):
        '''Create an input tile.
           Returns True if the tile exists upon exit.'''

        #print 'Making tile L:' + str(level) +' '+ str(tileIndex)

        outputTilePath = self.getTilePath(level, tileIndex, True)

        # On level zero we don't have to make the image
        if level == 0:
            # Resize the existing tile
            
            sourceFile = 'output_' + makeTileName(tileIndex)
            sourcePath = os.path.join(self.sourceFolder, sourceFile)
            if not os.path.exists(sourcePath): 
                return False # Quit here if the input file does not exist
            cmd = ('gdal_translate '+ sourcePath +' '+ outputTilePath +
                         ' -outsize '+ str(self._tileSize) +' '+ str(self._tileSize))
        else:
            # Need to assemble the tile from sources
            # TODO: Verify that all the source images exist!
        
            # Create the image file
            inputTileTiffString = self._getInputTileString(tileIndex)
            cmd = 'montage -resize 50% -background black -mode Concatenate -tile 2x2 ' + inputTileString +' '+ outputTilePath
            
        MosaicUtilities.cmdRunner(cmd, outputTilePath, force)

        # Create the KML file
        return self.makeKmlFile(level, tileIndex)
        
    def _getLatLonAltBox(self, level, tileIndex):
        '''Gets the ROI for one tile'''
        bbox = self._layerTilings[level].getTileBounds(tileIndex)
        return simplekml.LatLonAltBox(north=bbox.maxY, south=bbox.minY, west=bbox.minX, east=bbox.maxX)
    
    def _makeKmlRegion(self, level, tileIndex):
        '''Sets up the KML region object for a single tile'''
        return simplekml.Region(latlonaltbox=self._getLatLonAltBox(level, tileIndex),
                                lod=simplekml.Lod(minlodpixels=128, maxlodpixels=-1)) # TODO: Set this?
        
    def makeKmlFile(self, level, tileIndex):
        '''Creates the KML file for a given image tile.
           Returns True if the file exists upon exit.'''
        
        KML = False # Helper names for getTilePath
        TIF = True

        # Paths for the kml and associated tif file
        tileKmlPath   = self.getTilePath(level, tileIndex, KML)
        tileImagePath = self.getTilePath(level, tileIndex, TIF)
           
        ## There should be no need to recreate the KML files!
        #if os.path.exists(tileKmlPath):
        #    return True
        
        kml     = simplekml.Kml()
        # Region for this file
        kml.document.region = self._makeKmlRegion(level, tileIndex)
        
        # Network links for lower level files
        if level > 0:
            levelDown = level - 1
            indexList = self._getInputTiles(tileIndex)
            for index in indexList:
                sourceTilePath = self.getTilePath(levelDown, index, KML)
                if not os.path.exists(sourceTilePath):
                    continue # Don't link to non-existant tiles
                netlink = doc.newnetworklink(name=index.getPostFix(),
                                             region=self._makeKmlRegion(levelDown, index))
                netLink.link.href = sourceTilePath
                netLink.link.viewrefreshmode = simplekml.ViewRefreshMode.onrequest
        
        # Ground overlay
        kml.document.groundoverlay = kml.newgroundoverlay(icon=simplekml.Icon(href=tileImagePath),
                                                             latlonbox=self._getLatLonAltBox(level, tileIndex))

        # Save the completed file
        kml.save(tileKmlPath)
        return os.path.exists(tileKmlPath)


def main():


    #outputPath = ''
    #
    #
    #try:
    #    usage = "usage: geoTiffTool.py [--help][--manual] geotiffPath\n"# [lat] [lon] [lat2] [lon2]\n  "
    #    parser = optparse.OptionParser(usage=usage)
    #    parser.add_option("--manual", action="callback", callback=man,
    #                      help="Read the manual.")
    #
    #    ## Many GeoTiff images are stored in projected coordinates which makes this tricky
    #    #parser.add_option("--crop", action="store_true", dest="doCrop", default=False,
    #    #                            help="Crops the image between the four lat/lon arguments.")
    #
    #    parser.add_option("--normalize-eqc-lon", action="store_true", dest="normEqcLon", default=False,
    #                                help="Normalizes an EQC image to lie between +/- 180 degrees.")
    #
    #    # Add other features as the need arises
    #
    #    #parser.add_option("-o", "--output-path", dest="outputPath",
    #    #                        help="Output path (default replace extension with .kml")
    #    
    #    (options, args) = parser.parse_args()
    #
    #    # Parse positional arguments
    #    if not args: parser.error("need input cube file")           
    #    tiffPath = args[0]
    #
    #except optparse.OptionError, msg:
    #    raise Usage(msg)


    # First make a copy of each input tile in the output folder to form the base layer,
    #  but resize each of the tiles to 512 pixels to make the pyramid creation easier.

    sourceFolder = '/home/smcmich1/data/hrscMapTest/outputTiles/'
    outputFolder = '/home/smcmich1/data/hrscMapTest/kmlTree/'
    numLevels    = 1
    
    treeMaker = KmlTreeMaker(sourceFolder, outputFolder)

    for level in range(0, numLevels):
        treeMaker.makeLevel(level)

    return 0



if __name__ == "__main__":
    sys.exit(main())