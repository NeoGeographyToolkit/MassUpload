
import os
import sys
import re
import subprocess
import IrgGeoFunctions

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

# The proj4 string defining the space the base map is projected in
PROJ4_STRING = "+proj=eqc +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m"

# The order that HRSC channel images are stored and their names in a list
NUM_HRSC_CHANNELS = 5
HRSC_RED   = 0
HRSC_GREEN = 1
HRSC_BLUE  = 2
HRSC_NIR   = 3
HRSC_NADIR = 4
CHANNEL_STRINGS = ['red', 'green', 'blue', 'nir', 'nadir']

#----------------------------------------------------------------------------
# Functions

def getHrscChannelPaths(hrscBasePath):
    '''Get the list of all HRSC channel paths from the base path'''
    return [hrscBasePath+'_re3.tif', # TODO: Take out these space things later
            hrscBasePath+'_gr3.tif',
            hrscBasePath+'_bl3.tif',
            hrscBasePath+'_ir3.tif',
            hrscBasePath+'_nd3.tif']

class CmdRunException(Exception):
    '''Exception type indicating an error with a cmd call'''
    pass

def cmdRunner(cmd, outputPath=None, force=False):
    '''Executes a command if the output file does not exist and
       throws if the output file is not created'''

    if not os.path.exists(outputPath) or force:
        print cmd
        os.system(cmd)
    if outputPath and (not os.path.exists(outputPath)):
        raise CmdRunException('Failed to create output file: ' + outputPath)
    return True


def getImageSize(imagePath):
    """Returns the size [rows, columns] in an image"""
    # Make sure the input file exists
    if not os.path.exists(imagePath):
        raise Exception('Image file ' + imagePath + ' not found!')
       
    # Use subprocess to suppress the command output
    cmd = ['gdalinfo', imagePath]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()

    # Extract the size from the text
    sizePos    = textOutput.find('Size is')
    endPos     = textOutput.find('\n', sizePos+7)
    sizeStr    = textOutput[sizePos+7:endPos]
    sizeStrs   = sizeStr.strip().split(',')
    cols = int(sizeStrs[0])
    rows = int(sizeStrs[1])
    
    return (rows, cols)

def countBlackPixels(imagePath, isGray=True):
    '''Returns the number of black pixels in an image'''
    
    # Call imageMagick to print out the results
    if isGray:
        cmd = ['convert', imagePath, '-fill', 'white', '+opaque', 'gray(0)', '-format', '%c', 'histogram:info:']
    else: # RGB
        cmd = ['convert', imagePath, '-fill', 'white', '+opaque', 'rgb(0,0,0)', '-format', '%c', 'histogram:info:']
    #print cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    text, err = p.communicate()
    
    # Output looks like this:
    #     257066: (  0,  0,  0) #000000 black
    #     317182: (255,255,255) #FFFFFF white

    lines = text.split('\n')
    blackLine = lines[0].strip()
    otherLine = lines[1].strip()
    blackCount = int(blackLine[:blackLine.find(':')])
    #otherCount = int(otherLine[:otherLine.find(':')])
    return blackCount

# TODO: Stuff up here should come from common files
#======================================================================



def warpHrscFile(sourcePath, outputFolder, metersPerPixel):
    '''Warps an HRSC file into the same projection space as the base map'''
    fileName   = sourcePath[sourcePath.rfind('/')+1:]
    warpedPath = os.path.join(outputFolder, fileName)[:-4] + '_resample.tif'
    cmd = (GDAL_DIR+'gdalwarp ' + sourcePath +' '+ warpedPath + ' -r cubicspline '
             ' -t_srs "'+PROJ4_STRING+'" -tr '
             + str(metersPerPixel)+' '+str(metersPerPixel)+' -overwrite')
    cmdRunner(cmd, warpedPath)
    return warpedPath


def getTilePrefix(tileRow, tileCol):
    '''Return a standard string representation of a tile index'''
    return (str(tileRow)+'_'+str(tileCol))


def splitImage(imagePath, outputFolder, tileSize=512):
    '''Splits up an image into a grid of tiles and returns all the tile paths'''
    
    filename     = os.path.basename(imagePath)[:-4] # Strip extension
    outputPrefix = os.path.join(outputFolder, filename + '_tile_')

    # Skip tile creation if the first tile is present
    # - May need to make the decision smarter later on
    firstTilePath = outputPrefix + '0_0.tif'
    if not os.path.exists(firstTilePath):
    
        # This tile size is in the warped image (HRSC) resolution
        # --> May be able to 
        cmd = ('convert %s -crop %dx%d -set filename:tile "%%[fx:page.y/%d]_%%[fx:page.x/%d]" +repage +adjoin "%s%%[filename:tile].tif"'
                  % (imagePath, tileSize, tileSize, tileSize, tileSize, outputPrefix))
        print cmd
        os.system(cmd)
    
    # Build the list of output files
    outputTileInfoList = []
    for f in os.listdir(outputFolder):
        if '_tile_' not in f: # Skip any junk
            continue
        thisPath = os.path.join(outputFolder, f)
        numbers  =  re.findall(r"[\d']+", f) # Extract all numbers from the file name
        
        # Figure out the position of the tile
        tileRow  = int(numbers[3]) # In the tile grid
        tileCol  = int(numbers[4])
        pixelRow = tileRow * tileSize # In pixel coordinates relative to the original image
        pixelCol = tileCol * tileSize
        
        # Get other tile information
        height, width   = getImageSize(thisPath)
        totalNumPixels  = height*width
        blackPixelCount = countBlackPixels(thisPath)
        validPercentage = 1.0 - (blackPixelCount / totalNumPixels)
        
        thisTileInfo = {'path'        : thisPath,
                        'tileRow'     : tileRow,
                        'tileCol'     : tileCol,
                        'pixelRow'    : pixelRow,
                        'pixelCol'    : pixelCol,
                        'height'      : height,
                        'width'       : width,
                        'percentValid': validPercentage,
                        'prefix'      : getTilePrefix(tileRow, tileCol)
                       }
        outputTileInfoList.append(thisTileInfo)

    return outputTileInfoList


def writeSubSpatialTransform(fullPath, outputPath, pixelRow, pixelCol):
    '''Generates a spatial transform file for a single tile'''
    
    if os.path.exists(outputPath):
        return
    
    fIn  = open(fullPath,   'r')
    fOut = open(outputPath, 'w')
    
    fOut.write(fIn.readline()) # Copy the header
    
    xParts = fIn.readline().strip().split(',') # Read in the existing offset
    yParts = fIn.readline().strip().split(',')
    
    # Compute and write the new offsets
    # - We are using translation only affine transforms so this is simple
    newOffsetX = float(xParts[2]) + pixelCol
    newOffsetY = float(yParts[2]) + pixelRow
    fOut.write('%s, %s, %lf\n' % (xParts[0], xParts[1], newOffsetX))
    fOut.write('%s, %s, %lf\n' % (yParts[0], yParts[1], newOffsetY))
    
    fOut.write(fIn.readline()) # Copy the last line
    fIn.close()
    fOut.close()



def writeSubBrightnessGains(fullPath, outputPath, pixelRow, tileHeight):
    '''Generates a brightness gains file for a single tile'''
    
    if os.path.exists(outputPath):
        return

    fIn  = open(fullPath,   'r')
    fOut = open(outputPath, 'w')
    numInputLines = int(fIn.readline().strip())
    fOut.write(str(tileHeight)+'\n') # Write the number of rows in the tile
    for i in range(pixelRow): # Skip to the start of the pixels in the tile
        fIn.readline()
    
    for i in range(tileHeight): # Copy all the lines for this tile
        fOut.write(fIn.readline())
    
    fIn.close()
    fOut.close()
    

def getAdjancentTiles(tile, tileDict):
    '''Gets a list containing all the (still valid) tiles which are adjacent to the provided tile'''
    
    tileRow = tile['tileRow'] # The location of the input tile
    tileCol = tile['tileCol']
    
    # Build a list of the adjacent tiles
    adjacentTileList = []
    for r in range(-1,2):
        for c in range(-1,2):
            if (r==0) and (c==0): # Skip the main tile
                continue
            try:
                prefix  = getTilePrefix(tileRow+r, tileCol+c)
                adjTile = tileDict[prefix]
                if adjTile['stillValid']:
                    adjTile['rowOffset'] = r
                    adjTile['colOffset'] = c
                    adjacentTileList.append(adjTile)
            except:
                pass # This means this tile does not actually exist
    
    return adjacentTileList


def generateNewHrscColor(tile, tileDict, force=False):
    '''Generate a new color image from a single HRSC tile'''
    
    if not tile['stillValid']: # Skip tiles which have already failed
        return False
        
    try:
        # Get a list af adjacent tiles to pass in to the color transformer
        adjacentTiles = getAdjancentTiles(tile, tileDict)
        
        # Compute a weighting for each tile based on the pixel count
        totalWeight = 1.0 # The main image is the reference so it has weight 1 initially
        for adjTile in adjacentTiles:
            totalWeight += (adjTile['percentValid'] / tile['percentValid'])
            
        # Generate the parameter sequence for the next program call
        adjacentTileString = ''
        for adjTile in adjacentTiles:
            tileWeight     = (adjTile['percentValid'] / tile['percentValid']) / totalWeight
            thisTileString = adjTile['colorTransformPath'] +' '+ str(tileWeight) +' '+ str(adjTile['colOffset']) +' '+ str(adjTile['rowOffset'])
            adjacentTileString += (thisTileString + ' ')
        mainWeight = 1.0/totalWeight # Normalize the main weight too
        
        # Transform the HRSC image color
        # - Brightness correction is applied before color transform
        cmd = ('./transformHrscImageColor ' + tile['allTilesStringAndMask'] +' '+ tile['brightnessGainsPath'] +' '+ tile['newColorPath'] +' '+
                                              tile['colorTransformPath'] +' '+ str(mainWeight) +' '+ adjacentTileString
                                              )
        cmdRunner(cmd, tile['newColorPath'], force)
    except CmdRunException:
        tile['stillValid'] = False
    
    return True


def generateHrscColorImage(basemapCropPath, basemapGrayPath, hrscBasePathIn, outputFolder, metersPerPixel):
    '''Convert from HRSC color channels to an RGB image that matches the basemap colors'''

    # Set up all of the paths for this HRSC data set
    setName              = hrscBasePathIn[hrscBasePathIn.rfind('/')+1:]
    hrscBasePathOut      = os.path.join(outputFolder, setName)
    tileFolder           = hrscBasePathOut + '_tiles'
    spatialTransformPath = hrscBasePathOut+'_spatial_transform.csv'
    brightnessGainsPath  = hrscBasePathOut+'_brightness_gains.csv'
    hrscInputPaths       = getHrscChannelPaths(hrscBasePathIn)

    # Make sure the output directories exist
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
    if not os.path.exists(tileFolder):
        os.mkdir(tileFolder)
   
    # Transform the HRSC image to the same projection/resolution as the upsampled base map crop
    hrscWarpedPaths = [warpHrscFile(path, os.path.dirname(hrscBasePathOut), metersPerPixel) for path in hrscInputPaths]
       
    # For convenience generate a string containing all the input channel files
    hrscPathString = ''
    for path in hrscWarpedPaths:
        hrscPathString += path + ' '
    
    # Generate the spatial transform
    # - TODO: This can be done at a lower resolution than all the following steps!
    #         There is no point to running this at a higher resolution than the base map.
    #         Will have to scale the transform up to the full resolution.
    cmd = './RegisterHrsc ' + basemapGrayPath +' '+ hrscWarpedPaths[HRSC_NADIR] +' '+ spatialTransformPath
    cmdRunner(cmd, spatialTransformPath, False)

    
    # Compute the brightness scaling gains
    # - TODO: Do this in a tile-friendly way!
    #         Can do this by smoothing out the corrections after they have been computed per-tile.
    cmd = './computeBrightnessCorrection ' + basemapGrayPath +' '+ hrscPathString +' '+ spatialTransformPath +' '+ brightnessGainsPath
    cmdRunner(cmd, brightnessGainsPath, False)


    # Break the registered image up into tiles
    # - There is one list of tiles per HRSC channel
    # - Each channel gets its own subfolder
    TILE_SIZE = 512
    tileInfoLists = [[], [], [], [], []] # One list per channel
    for c in range(NUM_HRSC_CHANNELS):
        # Get the info for this channel
        warpedPath    = hrscWarpedPaths[c]
        channelString = CHANNEL_STRINGS[c]
        
        # Write all the tiles to a folder and get the info list for those tiles
        channelOutputFolder = os.path.join(tileFolder, channelString)
        if not os.path.exists(channelOutputFolder):
            os.mkdir(channelOutputFolder)
        tileInfoLists[c] = splitImage(warpedPath, channelOutputFolder, TILE_SIZE)
        
    # Verify that each channel generated the same number of tiles
    numTiles = len(tileInfoLists[0])
    for pathList in tileInfoLists:
        assert(len(pathList) == numTiles)
    print 'Generated ' +str(numTiles)+ ' tiles.'


    # Loop through each of the tiles we created and consolidate information across channels
    tileDict = {}
    for i in range(numTiles):

        # Start with the tile info from the first channel
        thisTileInfo   = tileInfoLists[0][i]
        allTilesString = ''
        # Add in paths to the tile for each channel, including a joint string for convenience.
        for c in range(NUM_HRSC_CHANNELS):
            pathKey = CHANNEL_STRINGS[c] + '_path'
            thisTileInfo[pathKey] = tileInfoLists[c][i]['path'],
            allTilesString       += tileInfoLists[c][i]['path'] + ' '
        thisTileInfo['allTilesString'] = allTilesString
        
        # Set up paths for the files we will generate for this tile
        filePrefix = 'tile_' + thisTileInfo['prefix']
        thisTileInfo['colorPairPath'        ] = os.path.join(tileFolder, filePrefix+'_color_pairs.csv')
        thisTileInfo['colorTransformPath'   ] = os.path.join(tileFolder, filePrefix+'_color_transform.csv')
        thisTileInfo['newColorPath'         ] = os.path.join(tileFolder, filePrefix+'_new_color.tif')
        thisTileInfo['spatialTransformPath' ] = os.path.join(tileFolder, filePrefix+'_spatial_transform.csv')
        thisTileInfo['brightnessGainsPath'  ] = os.path.join(tileFolder, filePrefix+'_brightness_gains.csv')
        thisTileInfo['tileMaskPath'         ] = os.path.join(tileFolder, filePrefix+'_tile_mask.tif')
        thisTileInfo['allTilesStringAndMask'] = allTilesString + ' ' + thisTileInfo['tileMaskPath']
        
        thisTileInfo['stillValid'] = True # Set this to false if there is an error processing this tile
        
        # Generate individual tile versions of the spatial transform and brightness gains
        writeSubSpatialTransform(spatialTransformPath, thisTileInfo['spatialTransformPath'],
                                 thisTileInfo['pixelRow'], thisTileInfo['pixelCol'])
        writeSubBrightnessGains(brightnessGainsPath, thisTileInfo['brightnessGainsPath'],
                                thisTileInfo['pixelRow'], thisTileInfo['height'])
    
        # Make a mask of valid pixels for this tile
        # - It would be great if the input images had a mask.
        # - Currently all-black pixels anywhere in the image get masked out!
        cmd = './makeSimpleImageMask ' + thisTileInfo['tileMaskPath'] +' '+ thisTileInfo['allTilesString']
        cmdRunner(cmd, thisTileInfo['tileMaskPath'], False)
        
        key = thisTileInfo['prefix']
        tileDict[key] = thisTileInfo
        
    
    # Now that we have all the per-tile info, compute the spatial transform for each tile.
    for tile in tileDict.itervalues():    
        
        try:
        
            # Generate the color pairs
            # - HRSC colors are written with the brightness correction already applied
            cmd = './writeHrscColorPairs ' + basemapCropPath +' '+ tile['allTilesStringAndMask'] +' '+ tile['spatialTransformPath'] +' '+ tile['brightnessGainsPath'] +' '+ tile['colorPairPath']
            cmdRunner(cmd, tile['colorPairPath'], False)
            
            # Compute the color transform
            cmd = 'python /home/smcmich1/repo/MassUpload/solveHrscColor.py ' + tile['colorTransformPath'] +' '+ tile['colorPairPath']
            cmdRunner(cmd, tile['colorTransformPath'], False)
            
        except CmdRunException: # Flag errors with this tile
            tile['stillValid'] = False

# TODO: Be robust to tile failures -> Some will definately lie outside the image!

    # Now that we have all the spatial transforms, generate the new color image for each tile.
    for tile in tileDict.itervalues():
        
        # New color generation is handled in this function
        generateNewHrscColor(tile, tileDict, True)
        
        #raise Exception('DEBUG')

    #raise Exception('DEBUG')


    # Return all of the generated tile information!
    return tileDict






#-------------------------------------------------------------------


fullBasemapPath      = '/home/smcmich1/data/hrscMapTest/noel_basemap.tif'
cropBasemapSmallPath = '/home/smcmich1/data/hrscMapTest/basemap_crop_small.tif'
cropBasemapPath      = '/home/smcmich1/data/hrscMapTest/basemap_crop.tif'
cropBasemapGrayPath  = '/home/smcmich1/data/hrscMapTest/basemap_crop_red.tif'


hrscBasePathInList = ['/home/smcmich1/data/hrscMapTest/external_data/h0022_0000',
                      '/home/smcmich1/data/hrscMapTest/external_data/h0506_0000',
                      '/home/smcmich1/data/hrscMapTest/external_data/h2411_0000',
                      '/home/smcmich1/data/hrscMapTest/external_data/h6419_0000']

outputFolder    = '/home/smcmich1/data/hrscMapTest/'

GDAL_DIR = '/home/smcmich1/programs/gdal-1.11.0-install/bin/'

print 'Starting basemap enhancement script...'


#---------------------------
# Prep the base map
# - For now we extract an arbitrary chunk, later this will be tiled.
#
## Get the HRSC bounding box and expand it
#HRSC_BB_EXPAND_DEGREES = 1.5
#(minLon, maxLon, minLat, maxLat) = IrgGeoFunctions.getGeoTiffBoundingBox(hrscBasePathInList[0]+'_nd3.tif')
#minLon -= HRSC_BB_EXPAND_DEGREES
#maxLon += HRSC_BB_EXPAND_DEGREES
#minLat -= HRSC_BB_EXPAND_DEGREES
#maxLat += HRSC_BB_EXPAND_DEGREES
#
#print 'Region bounds:' + str((minLon, maxLon, minLat, maxLat))
#
## Convert the bounding box from degrees to the projected coordinate system (meters)
DEGREES_TO_PROJECTION_METERS = 59274.9
NOEL_MAP_METERS_PER_PIXEL = 1852.4
#minX = minLon*DEGREES_TO_PROJECTION_METERS
#maxX = maxLon*DEGREES_TO_PROJECTION_METERS
#minY = minLat*DEGREES_TO_PROJECTION_METERS
#maxY = maxLat*DEGREES_TO_PROJECTION_METERS
#
## Crop out the correct section of the base map
#projCoordString = '%f %f %f %f' % (minX, maxLat, maxX, minY)
#cmd = (GDAL_DIR+'gdal_translate ' + fullBasemapPath +' '+ cropBasemapSmallPath
#                         +' -projwin '+ projCoordString)
#cmdRunner(cmd, cropBasemapSmallPath)
#
## Increase the resolution of the cropped image
## TODO: Can this be done in one step?
RESOLUTION_INCREASE = 200 # In percent
#cmd = (GDAL_DIR+'gdal_translate ' + cropBasemapSmallPath +' '+ cropBasemapPath
#       +' -outsize '+str(RESOLUTION_INCREASE)+'% '+str(RESOLUTION_INCREASE)+'% ')
#cmdRunner(cmd, cropBasemapPath)
#
## Generate the grayscale version of the cropped basemap
#cmd = (GDAL_DIR+'gdal_translate -b 1 ' + cropBasemapPath +' '+ cropBasemapGrayPath)
#cmdRunner(cmd, cropBasemapGrayPath)

# In the future:
#   For each chunk of the basemap, find all the HRSC images that overlap it.
#   Throw out images which are flagged as bad.

#------------------------------------
# Process the individual HRSC images and add them to the mosaic

# TODO: We need HRSC masks to handle the borders properly!!!

for hrscPath in hrscBasePathInList: # Loop through input HRSC images

    # Transform the HRSC image to the same projection/resolution as the upsampled base map crop
    metersPerPixel = NOEL_MAP_METERS_PER_PIXEL / (RESOLUTION_INCREASE/100.0)
    tileDict = generateHrscColorImage(cropBasemapPath, cropBasemapGrayPath, hrscPath, outputFolder, metersPerPixel)

    # TODO: Need to handle tiling of the output mosaic also!
    # TODO: This C++ program can do multiple tiles in one call.
    mosaicPath = '/home/smcmich1/data/hrscMapTest/outputMosaic.tif'
    for tile in tileDict.itervalues():
    
        if not tile['stillValid']: # Skip tiles which have already failed
            continue
    
        try:
            tileParamsStr = mosaicPath +' '+ tile['newColorPath'] +' '+ tile['tileMaskPath'] +' '+ tile['spatialTransformPath']
            if os.path.exists(mosaicPath):
                cmd = './hrscMosaic ' + mosaicPath +' '+ tileParamsStr
            else: # Only the first call
                cmd = './hrscMosaic ' + cropBasemapPath +' '+ tileParamsStr
            cmdRunner(cmd, mosaicPath, True)
            
        except CmdRunException:
            tile['stillValid'] = False

    #raise Exception('DEBUG')

print 'Basemap enhancement script completed!'

"""
== pansharp prototype ==

Generate RED version of Noel's map (ONCE)
Equitorial circumference = 21338954.25548 / 2 = 10669477.1 <--- x span
Polar      circumference = 21338954.25548 / 4 =  5334738.6 <--- y span
gdal_translate RA_Albedo.tif noel_basemap.tif
               -a_srs "+proj=eqc +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m"
               -a_ullr -10669477.1 5334738.6 10669477.1 -5334738.6
gdal_translate noel_basemap.tif noel_basemap_red.tif -b 1
DONE

FETCH HRSC IMAGE

Pick out HRSC image of interest and download ND3 file
DONE

PREP INPUT IMAGES

--> TODO: Try these steps out a increased resolution from the base map!
            2x or 4x should be good.

Generate a reduced-resolution of the file to match Noel's map using GDAL_translate.
- Noel's map resolution is (circumference 21,344 km/ 11520 pixels) = 1852.4 meters per pixel.

Get a bounding box surrounding the HRSC image and extract that region from Noel's map.

Min latitude  = -50.1955
Max latitude  = -18.6077
Min longitude = 98.2023
Max longitude = 102.64

Convert from lat/lon to projected coordinates by multiplying degrees with these constants:
meters per degree lon = 59274.9
meters per degree lat = 59274.9

gdal_translate noel_basemap_red.tif red_crop.tif -projwin 97.2023 -17.6077 103.64 -51.1955
 5761656.61227 -1043694.65673  6143250.636 -3034608.14295


Transform the HRSC image to the same projection/resolution as the base map
~/programs/gdal-1.11.0-install/bin/gdalwarp h0022_0000_nd3.tif h0022_0000_nd3_resample.tif -t_srs "+proj=eqc +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m"  -tr 1852  1852 -overwrite


ALIGN THE IMAGES

Compute a pixel to pixel match between those two images.
- HRSC to NOEL => 24.384001, 32.455997 (basemap res)


GENERATE SHARPENED BASE TILE

Generate an unsharp image from the HRSC pan image.

Iterate through the base tile and for each pixel use the transform to find the matching point
 in the unsharp image, then add the unsharp image value to the base tile.  Also create a mask
 showing which pixels were modified (many won't overlap with the HRSC image).



Manually checked spatial transforms:
h0022 = 282, 1190
h0506 = 80.6, 1323
h2411 = 243, 2389
h6419 = 295, 1829
---> Need to be able to compute these automatically!


"""



"""
DB command to find a list of overlapping HRSC files:
select setname from Files where sensor=1 and minLon<102.64 and maxLon>98.2023 and minLat<-18.6077 and maxLat>-50.1955 and subtype="nd3";


Alignment target = /byss/mars/themis/dayir_100m/
- Noel's map goes down to +/-90 degrees lat, but matching will gets very hard outside 60 degrees!
- Alignment between Noel and Themis is close but not exact, definately possible with ASP tools.



A professionally made partial mosaic: http://maps.planet.fu-berlin.de/
    - This is a hard problem, our only advantages are Earth Engine and the color reference map!

Query sqlite3 database to get a list of all the HRSC data sets after a certain date
    3500 sets -> 18000 files!

From each of these, we will need the following:
    ir3
    nd3
    red
    green
    blue

Align the nd3 image to the THEMIS Day Time global map.
    Vision workbench?
    
Generate a grid of pixel tie points between the nd3 and THEMIS

Generate a grid of pixel tie points between the mars color reference and THEMIS
    Only do this once!

Make a list of pixel pairs: HRSC - COLOR REF

Compute 5x3 transform using pixel pairs

Transfrom HRSC bands to a new RGB representation using the transform

Apply any image enhancement algorithms

Transform the RGB image to align with the THEMIS base map

While we are at it, try applying PAN unsharp mask to the upsampled color base layer!
    - This might be the only way to easily make a good mosaic


OUTSTANDING ISSUES
- Dark side color handling?
    - Ensure fixed brightness across image?
        - What about dark spots and poles?
        - Maybe brightness profile changed to match the color reference map?
    - Transform changes gradually down image?
        - Think of color transform as a rotation that smoothly changes
        - Can this be done between images as well?
    HRSC_h0344_0000_nd3
    HRSC_h0543_0000_nd3
    HRSC_h2011_0001_nd3

- Streak handling and other random crap
    HRSC_h0334_0000_nd3
    HRSC_h1970_0000_nd3

- Alignment verification

- Seam blending --> Making sure the colors match!
    - Images are captured at all sorts of light levels and angles so even a grayscale
      mosaic is hard to get right (see http://maps.planet.fu-berlin.de/)


- Can we do everything in OpenCV?
    - The largest file is 13,500x150,000
    - Image alignment is done at a lower resolution so no problem.
    - Everything else can be done in dumb tiles.


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




