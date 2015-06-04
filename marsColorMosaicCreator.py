
import os
import sys
import re
import subprocess
import numpy
import IrgGeoFunctions
import copyGeoTiffInfo
import mosaicTileManager # TODO: Normalize caps!
import MosaicUtilities
import hrscImageManager

"""
TODO:

- Algorithm to smooth out transitions.
- Generate a good result image.
- Generate a list of features Earth Engine would need to replicate all steps.
- Switch the code to a real tiling scheme.
- Make sure the same set HRSC images line up properly
    - This will become apparent at higher resolutions!
    - We can't just treat them as perfectly aligned!
    --> Treat Nadir image as the standard and align other images to it.
- Better mask handling

Existing tools:
- RegisterHrsc.cpp
    - Input  = Basemap, HRSC
    - Output = spatialTransform
- writeHrscColorPairs.cpp
    - Input  = Basemap, HRSC, spatialTransform
    - Output = File containing pixel color pairs
- transformHrscImageColor.cpp
    - Input  = HRSC, colorTransform
    - Output = Color transformed HRSC image
    - TODO   = Add cleanup/pansharp


"""

#----------------------------------------------------------------------------
# Constants


#-----------------------------------------------------------------------------------------
# Functions


def getHrscImageList():
    '''For just returns a fixed list of HRSC images for testing'''
    
    # TODO: Search our database to perform this function!
    return ['h0022_0000']#,
            #'h0506_0000']#,
            #'h2411_0000',
            #'h6419_0000')

def getCoveredOutputTiles(basemapInstance, hrscInstance):
    '''Return a bounding box containing all the output tiles covered by the HRSC image'''
    
    hrscBoundingBoxDegrees = hrscInstance.getBoundingBoxDegrees()
    return basemapInstance.getIntersectingTiles(hrscBoundingBoxDegrees)

    #return MosaicUtilities.Rectangle(196, 197, 92, 93) # DEBUG

def updateTileWithHrscImage(basemapInstance, tileIndex, outputTilePath, hrscInstance):
    '''Update a single output tile with the given HRSC image'''

    thisTileBounds = basemapInstance.getTileRectDegree(tileIndex)
    #hrscSetName    = hrscInstance.getSetName()


    # Get the tile information from the HRSC image
    tileDict = hrscInstance.getTileInfo(thisTileBounds)
    print 'Found these tile intersections:'
    for hrscTile in tileDict.itervalues():
        print hrscTile['prefix']


    # TODO: This C++ program can do multiple tiles in one call.

    # For each tile...
    i = 0
    for hrscTile in tileDict.itervalues(): 
    
        #try:
        cmd = ('./hrscMosaic ' + outputTilePath +' '+ outputTilePath +' '+ hrscTile['newColorPath'] +' '+
                                  hrscTile['tileMaskPath'] +' '+ hrscTile['tileToTileTransformPath'])
        MosaicUtilities.cmdRunner(cmd, outputTilePath, True)

        #i += 1
        #if i == 3:
        #    raise Exception('DEBUG')
    
        # TODO: Make sure the output tiles still have their geo info





def updateTilesContainingHrscImage(basemapInstance, hrscInstance):
    '''Updates all output tiles containing this HRSC image'''

    # Find all the output tiles that intersect with this
    outputTilesRect = getCoveredOutputTiles(basemapInstance, hrscInstance)

    print 'Found overlapping output tiles:  ' + str(outputTilesRect)
    
    # Loop through all the tiles
    for row in range(outputTilesRect.minY, outputTilesRect.maxY):
        for col in range(outputTilesRect.minX, outputTilesRect.maxX):
    
            # Set up the til information
            tileIndex  = MosaicUtilities.TileIndex(row, col) #basemapInstance.getTileIndex(98.5, -27.5)
            tileBounds = basemapInstance.getTileRectDegree(tileIndex)
            
            print 'Using HRSC image ' + hrscSetName + ' to update tile: ' + str(tileIndex)
            print '--> Tile bounds = ' + str(tileBounds)

            print '\nMaking sure basemap info is present...'
            
            # Now that we have selected a tile, generate all of the tile images for it.
            (smallTilePath, largeTilePath, grayTilePath, outputTilePath) = basemapInstance.generateTileImages(tileIndex, False)
        
            print '\nPasting on HRSC tiles...'

            #raise Exception('DEBUG')
        
            # Update the selected tile with the HRSC image
            updateTileWithHrscImage(basemapInstance, tileIndex, outputTilePath, hrscInstance)
            
            #raise Exception('DEBUG')
        
    print '\n---> Finished updating tiles for HRSC image ' + hrscSetName

#-----------------------------------------------------------------------------------------


fullBasemapPath  = '/home/smcmich1/data/hrscMapTest/projection_space_basemap.tif'

sourceHrscFolder = '/home/smcmich1/data/hrscMapTest/external_data'

hrscOutputFolder = '/home/smcmich1/data/hrscMapTest/hrscFiles'

outputTileFolder = '/home/smcmich1/data/hrscMapTest/outputTiles'



print 'Starting basemap enhancement script...'


# TODO: Go down to 20 to 10 meters!
OUTPUT_RESOLUTION_METERS_PER_PIXEL = 100

print '\n==== Initializing the base map object ===='
basemapInstance = mosaicTileManager.MarsBasemap(fullBasemapPath, outputTileFolder, OUTPUT_RESOLUTION_METERS_PER_PIXEL)

print '--- Finished initializing the base map object ---\n'

# Get a list of all the HRSC images we are testing with
hrscImageList = getHrscImageList()

# Loop through input HRSC images
for hrscSetName in hrscImageList: 
    
    # Pick a location to store the data for this HRSC image
    thisHrscFolder = os.path.join(hrscOutputFolder, hrscSetName)

    #try:

    print '\n=== Initializing HRSC image ' + hrscSetName + ' ==='

    # Fetch and preprocess the HRSC image
    # - TODO: Use the manager class to handle this!
    hrscInstance = hrscImageManager.HrscImage(hrscSetName, sourceHrscFolder, thisHrscFolder, basemapInstance)

    print '--- Now initializing high res HRSC content ---'

    # Complete the high resolution components
    hrscInstance.prepHighResolutionProducts()
    
    print '--- Finished initializing HRSC image ---\n'

    # Call the function to update all the output images for this HRSC image
    updateTilesContainingHrscImage(basemapInstance, hrscInstance)  

    print '<<<<< Finished writing all tiles for this HRSC image! >>>>>'

    # TODO: Clean up if necessary

    raise Exception('DEBUG')



print 'Basemap enhancement script completed!'

"""
== pansharp prototype ==

Generate RED version of Noel's map (ONCE)
Equitorial circumference = 21338954.25548 / 2 = 10669477.1 <--- x span
Polar      circumference = 21338954.25548 / 4 =  5334738.6 <--- y span
gdal_translate mars_full_albedo.tif projection_space_basemap.tif -a_srs "+proj=eqc +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m" -a_ullr -10669477.1 5334738.6 10669477.1 -5334738.6
gdal_translate noel_basemap.tif noel_basemap_red.tif -b 1
DONE

Generate a reduced-resolution of the file to match Noel's map using GDAL_translate.
- Noel's map resolution is (circumference 21,344 km/ 11520 pixels) = 1852.4 meters per pixel.

Get a bounding box surrounding the HRSC image and extract that region from Noel's map.

Manually checked spatial transforms (200% resolution):
h0022 = 282, 1190
h0506 = 80.6, 1323
h2411 = 243, 2389
h6419 = 295, 1829
---> Need to be able to compute these automatically!

Basemap resolution transforms:
h0022 = 41, 49
h0506 = -61, 115
h2411 = 122, 583
h6419 = 

DB command to find a list of overlapping HRSC files:
select setname from Files where sensor=1 and minLon<102.64 and maxLon>98.2023 and minLat<-18.6077 and maxLat>-50.1955 and subtype="nd3";


A professionally made partial mosaic: http://maps.planet.fu-berlin.de/
    - This is a hard problem, our only advantages are Earth Engine and the color reference map!

Query sqlite3 database to get a list of all the HRSC data sets after a certain date
    3500 sets -> 18000 files!


List of overlapping images for h0022_0000_nd3:

Use these:
h0506_0000_nd3 <-- Feels like there should be overlap =)
h2411_0000_nd3 <-- Lots of overlap!
h6419_0000_nd3 <-- Higher resolution shot of middle of image


h0248_0000_nd3
h0300_0000_nd3
h0440_0000_nd3 <-- Big image, but no overlap =(
h0451_0000_nd3 <-- Messed up top!
h0462_0000_nd3
h0506_0000_nd3 <-- Feels like there should be overlap =)
h0637_0000_nd3
h0648_0000_nd3
h1592_0000_nd3 <-- Weird polar image
h1774_0000_nd3
h1786_0001_nd3
h1942_0000_nd3 <-- Small image in the center
h2345_0000_nd3
h2389_0000_nd3
h2400_0001_nd3 <-- Lots of overlap, but not all components =()
h2411_0000_nd3 <-- Lots of overlap!
h2466_0000_nd3 
h2510_0001_nd3
h2619_0000_nd3
h2652_0000_nd3
h2663_0001_nd3
h2726_0000_nd3
h2729_0000_nd3
h4261_0000_nd3 <-- Maybe aligned?
h4272_0000_nd3
h4294_0000_nd3
h4469_0001_nd3
h4642_0000_nd3
h4817_0000_nd3
h6411_0000_nd3
h6419_0000_nd3 <-- Higher resolution shot of middle of image
h6437_0000_nd3
h8429_0000_nd3
h8474_0000_nd3
h8562_0000_nd3
h8990_0000_nd3
ha526_0000_nd3
ha688_0000_nd3
ha776_0000_nd3
hc598_0017_nd3

"""




