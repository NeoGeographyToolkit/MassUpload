
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

# The proj4 string defining the space the base map is projected in
#PROJ4_STRING = "+proj=eqc +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m"






## TODO: Move to a general file
#def projCoordToPixelCoord(x, y, geoInfo):
#    '''Converts from projected coordinates into pixel coordinates using
#       the results from getImageGeoInfo()'''
#
#    xOffset = x - geoInfo['projection_bounds'][0]
#    yOffset = y - geoInfo['projection_bounds'][3]
#    
#    column = xOffset / geoInfo['pixel_size'][0]
#    row    = yOffset / geoInfo['pixel_size'][1]
#    
#    return (column, row)


#def estimateRegistration(baseImage, otherImage, outputPath):
#    '''Writes an estimated registration transform to a file based on geo metadata.'''
#    # This function assumes the images are in the same projection system!
#    
#    # Get the projection bounds and size in both images
#    baseGeoInfo     = IrgGeoFunctions.getImageGeoInfo(baseImage,  False)
#    otherGeoInfo    = IrgGeoFunctions.getImageGeoInfo(otherImage, False)
#    baseProjBounds  = baseGeoInfo[ 'projection_bounds']
#    otherProjBounds = otherGeoInfo['projection_bounds']
#    baseImageSize   = baseGeoInfo[ 'image_size']
#    otherImageSize  = otherGeoInfo['image_size']
#    
#    #print baseGeoInfo
#    #print '----'
#    #print otherGeoInfo
#    
#    # Now estimate the bounding box of the other image in the base image
#    topLeftCoord = projCoordToPixelCoord(otherProjBounds[0], otherProjBounds[3], baseGeoInfo)
#    
#    # Write the output file
#    with open(outputPath, 'w') as f:
#        f.write('3, 3\n')
#        f.write('1, 0, %lf\n' % (topLeftCoord[0]))
#        f.write('0, 1, %lf\n' % (topLeftCoord[1]))
#        f.write('0, 0, 1\n')
#    if not os.path.exists(outputPath):
#        raise Exception('Failed to create transform file ' + outputPath)
#
#    return topLeftCoord


#def splitImage(imagePath, outputFolder, tileSize=512):
#    '''Splits up an image into a grid of tiles and returns all the tile paths'''
#    
#    filename     = os.path.basename(imagePath)[:-4] # Strip extension
#    outputPrefix = os.path.join(outputFolder, filename + '_tile_')
#
#    # Skip tile creation if the first tile is present
#    # - May need to make the decision smarter later on
#    firstTilePath = outputPrefix + '0_0.tif'
#    if not os.path.exists(firstTilePath):
#    
#        # This tile size is in the warped image (HRSC) resolution
#        # --> May be able to 
#        cmd = ('convert %s -crop %dx%d -set filename:tile "%%[fx:page.y/%d]_%%[fx:page.x/%d]" +repage +adjoin "%s%%[filename:tile].tif"'
#                  % (imagePath, tileSize, tileSize, tileSize, tileSize, outputPrefix))
#        print cmd
#        os.system(cmd)
#    
#    # Build the list of output files
#    outputTileInfoList = []
#    for f in os.listdir(outputFolder):
#        if '_tile_' not in f: # Skip any junk
#            continue
#        thisPath = os.path.join(outputFolder, f)
#        numbers  =  re.findall(r"[\d']+", f) # Extract all numbers from the file name
#        
#        # Figure out the position of the tile
#        tileRow  = int(numbers[3]) # In the tile grid
#        tileCol  = int(numbers[4])
#        pixelRow = tileRow * tileSize # In pixel coordinates relative to the original image
#        pixelCol = tileCol * tileSize
#        
#        # Get other tile information
#        height, width   = getImageSize(thisPath)
#        totalNumPixels  = height*width
#        blackPixelCount = countBlackPixels(thisPath)
#        validPercentage = 1.0 - (blackPixelCount / totalNumPixels)
#        
#        thisTileInfo = {'path'        : thisPath,
#                        'tileRow'     : tileRow,
#                        'tileCol'     : tileCol,
#                        'pixelRow'    : pixelRow,
#                        'pixelCol'    : pixelCol,
#                        'height'      : height,
#                        'width'       : width,
#                        'percentValid': validPercentage,
#                        'prefix'      : getTilePrefix(tileRow, tileCol)
#                       }
#        outputTileInfoList.append(thisTileInfo)
#
#    return outputTileInfoList
#
#
#def writeSubSpatialTransform(fullPath, outputPath, pixelRow, pixelCol, force):
#    '''Generates a spatial transform file for a single tile'''
#    
#    if os.path.exists(outputPath) and (not force):
#        return
#    
#    fIn  = open(fullPath,   'r')
#    fOut = open(outputPath, 'w')
#    
#    fOut.write(fIn.readline()) # Copy the header
#    
#    xParts = fIn.readline().strip().split(',') # Read in the existing offset
#    yParts = fIn.readline().strip().split(',')
#    
#    # Compute and write the new offsets
#    # - We are using translation only affine transforms so this is simple
#    newOffsetX = float(xParts[2]) + pixelCol
#    newOffsetY = float(yParts[2]) + pixelRow
#    fOut.write('%s, %s, %lf\n' % (xParts[0], xParts[1], newOffsetX))
#    fOut.write('%s, %s, %lf\n' % (yParts[0], yParts[1], newOffsetY))
#    
#    fOut.write(fIn.readline()) # Copy the last line
#    fIn.close()
#    fOut.close()
#    
#
#def scaleSpatialTransform(inputPath, outputPath, scale, offsetX, offsetY, force=False):
#    '''Scale up a spatial transform for a higher resolution'''
#    
#    if os.path.exists(outputPath) and (not force):
#        return
#    
#    fIn  = open(inputPath,  'r')
#    fOut = open(outputPath, 'w')
#    
#    fOut.write(fIn.readline()) # Copy the header
#    
#    xParts = fIn.readline().strip().split(',') # Read in the existing offset
#    yParts = fIn.readline().strip().split(',')
#    
#    # Compute and write the new offsets
#    # - We are using translation only affine transforms so this is simple
#    #print xParts
#    #print yParts
#    #print offsetX
#    #print offsetY
#    newOffsetX = (float(xParts[2]) + offsetX)* scale
#    newOffsetY = (float(yParts[2]) + offsetY) * scale
#    fOut.write('%s, %s, %f\n' % (xParts[0], xParts[1], newOffsetX))
#    fOut.write('%s, %s, %f\n' % (yParts[0], yParts[1], newOffsetY))
#    
#    fOut.write(fIn.readline()) # Copy the last line
#    fIn.close()
#    fOut.close()
#
#
#def computeSpatialRegistration(basemapInstance, tileIndex, hrscPath, outputPath, outputPathLowRes, tempFolder='', force=False):
#    '''Compute the spatial registration from the HRSC image to the base image'''
#
#    # Estimate the spatial transform using the image metadata
#    # - This is always done against the full image basemap
#    estimatedTransformPath     = os.path.join(tempFolder, 'spatial_transform_basemap_estimated.csv')
#    basemapLowResTransformPath = os.path.join(tempFolder, 'spatial_transform_basemap_low.csv')
#    estX, estY = estimateRegistration(basemapInstance.getRegistrationFile(), hrscPath, estimatedTransformPath)
#    
#    #raise Exception('DEBUG')
#
#    # TODO: Check the number of inliers!    
#    # Refine the spatial transform using image data
#    # - This is computed at the base map resolution and then scaled up to the output resolution
#    cmd = ('./RegisterHrsc ' + basemapInstance.getRegistrationFile() +' '+ hrscPath
#           +' '+ basemapLowResTransformPath +' '+ str(1) +' '+ estimatedTransformPath)
#    cmdRunner(cmd, basemapLowResTransformPath, force)
#    
#    #raise Exception('DEBUG')
#
#    # Convert the transform to be relative to the output tile
#    # - Do this in both low (basemap) and high (output) resolution
#    tileBounds = basemapInstance.getLowResPixelBounds(tileIndex)
#    scaleSpatialTransform(basemapLowResTransformPath, outputPath, basemapInstance.getResolutionIncrease(),
#                          -tileBounds.minX, -tileBounds.minY, force)
#    scaleSpatialTransform(basemapLowResTransformPath, outputPathLowRes, 1.0,
#                          -tileBounds.minX, -tileBounds.minY, force)
#    
#    # Clear temporary files
#    #os.remove(estimatedTransformPath)
#    #os.remove(basemapLowResTransformPath)
#    
#    #raise Exception('DEBUG')


#def writeSubBrightnessGains(fullPath, outputPath, pixelRow, tileHeight, scale):
#    '''Generates a brightness gains file for a single tile'''
#    
#    if os.path.exists(outputPath):
#        return
#    
#    fIn = open(fullPath,   'r')
#    fOut = open(outputPath, 'w')
#    numInputLines = int(fIn.readline().strip())
#    fOut.write(str(tileHeight)+'\n') # Write the number of rows in the tile
#    for i in range(pixelRow): # Skip to the start of the pixels in the tile
#        fIn.readline()
#    
#    for i in range(tileHeight): # Copy all the lines for this tile
#        fOut.write(fIn.readline())
#    
#    fIn.close()
#    fOut.close()



#def splitScaleBrightnessGains(fullPath, tileDict, scale):
#    '''Generates a brightness gains file for a single tile'''
#
#    print 'Generating split brightness gains...'
#    
#    # Read in the entire input file
#    lowResVals = numpy.loadtxt(fullPath, skiprows=1, delimiter=',', usecols=(0,))   
#    lowResRows = range(0,len(lowResVals))
#    
#    for tile in tileDict.itervalues():
#        
#        outputPath = tile['brightnessGainsPath']
#        if os.path.exists(outputPath):
#            continue
#        
#        # Compute the row range in the input tile
#        tileHeight       = tile['height'  ]
#        pixelRow         = tile['pixelRow']
#        fullSizeStartRow = pixelRow
#        fullSizeStopRow  = (pixelRow+tileHeight)
#        
#        # For each desired ouput value, compute the location in the input values
#        thisTileRowsInInput = numpy.empty([tileHeight])
#        index = 0
#        for r in range(fullSizeStartRow, fullSizeStopRow):
#            thisTileRowsInInput[index] = r / scale
#            index += 1
#        # Use numpy to interpolate values 
#        thisTileVals = numpy.interp(thisTileRowsInInput, lowResRows, lowResVals)
#        zeroCol      = [0 for i in thisTileVals]
#            
#        # Write out the interpolated values
#        numpy.savetxt(outputPath, thisTileVals, header=str(tileHeight),  fmt='%1.6f, 0.0', comments='')  
#    
#
#
#def getAdjacentTiles(tile, tileDict):
#    '''Gets a list containing all the (still valid) tiles which are adjacent to the provided tile'''
#    
#    tileRow = tile['tileRow'] # The location of the input tile
#    tileCol = tile['tileCol']
#    
#    # Build a list of the adjacent tiles
#    adjacentTileList = []
#    for r in range(-1,2):
#        for c in range(-1,2):
#            if (r==0) and (c==0): # Skip the main tile
#                continue
#            try:
#                prefix  = getTilePrefix(tileRow+r, tileCol+c)
#                adjTile = tileDict[prefix]
#                if adjTile['stillValid']:
#                    adjTile['rowOffset'] = r
#                    adjTile['colOffset'] = c
#                    adjacentTileList.append(adjTile)
#            except:
#                pass # This means this tile does not actually exist
#    
#    return adjacentTileList


def generateNewHrscColor(tile, tileDict, force=False):
    '''Generate a new color image from a single HRSC tile'''
    
    if not tile['stillValid']: # Skip tiles which have already failed
        return False
        
    try:
        # Get a list af adjacent tiles to pass in to the color transformer
        adjacentTiles = getAdjacentTiles(tile, tileDict)
        
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


def generateHrscColorImage(basemapInstance, thisTileIndex, hrscPrefixIn, outputFolder):
    '''Convert from HRSC color channels to an RGB image that matches the basemap colors'''

    # Set up all of the paths for this HRSC data set
    setName                    = hrscPrefixIn[hrscPrefixIn.rfind('/')+1:]
    hrscBasePathOut            = os.path.join(outputFolder, setName)
    tileFolder                 = hrscBasePathOut + '_tiles'
    brightnessGainsPath        = hrscBasePathOut+'_brightness_gains.csv'
    lowResMaskPath             = hrscBasePathOut+'_low_res_mask.tif'
    hrscInputPaths             = getHrscChannelPaths(hrscPrefixIn)

    tileFolder                 = basemapInstance.getTileFolder(thisTileIndex)
    spatialTransformPathLowRes = os.path.join(tileFolder, setName+'_spatial_transform_lowRes.csv')
    spatialTransformPath       = os.path.join(tileFolder, setName+'_spatial_transform.csv')

    forceFromHere = False # Force recomputation

    # Make sure the output directories exist
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
    if not os.path.exists(tileFolder):
        os.mkdir(tileFolder)
           
    #forceFromHere = True
           
    # Transform the HRSC image to the same projection/resolution as the upsampled base map crop
    hrscWarpedPaths = [warpHrscFile(path, os.path.dirname(hrscBasePathOut), '_output_res', basemapInstance.getHighResMpp(), False)
                           for path in hrscInputPaths]
    # Also do a set at the original base map resolution
    lowResWarpedPaths = [warpHrscFile(path, os.path.dirname(hrscBasePathOut), '_basemap_res', basemapInstance.getLowResMpp(), False)
                           for path in hrscInputPaths]
    
    # For convenience generate strings containing all the input channel files
    hrscPathString = ''
    for path in hrscWarpedPaths:
        hrscPathString += path + ' '
    lowResHrscPathString = ''
    for path in lowResWarpedPaths:
        lowResHrscPathString += path + ' '

    # Make a mask at the low resolution
    cmd = './makeSimpleImageMask ' + lowResMaskPath +' '+ lowResHrscPathString
    cmdRunner(cmd, lowResMaskPath, forceFromHere)

    # Compute the spatial registration from the HRSC image to the basemap
    computeSpatialRegistration(basemapInstance, thisTileIndex, lowResWarpedPaths[HRSC_NADIR],
                               spatialTransformPath, spatialTransformPathLowRes, outputFolder, forceFromHere)
    
    # Compute the brightness scaling gains
    # - This is done at low resolution
    # - The low resolution output is smoothed out later to avoid jagged edges.
    cmd = './computeBrightnessCorrection ' + lowResGrayBaseTile +' '+ lowResHrscPathString +' '+ spatialTransformPathLowRes +' '+ brightnessGainsPath
    cmdRunner(cmd, brightnessGainsPath, forceFromHere)


    raise Exception('DEBUG')

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

        # Check the percent valid to see if there is any content in this tile
        # - Percent valid is the only field that varies by channel
        # - All the other information can just be read from the first channel
        thisTileInfo = tileInfoLists[0][i]
        percentValid = 1.0
        for c in range(0,NUM_HRSC_CHANNELS):
            chanInfo = tileInfoLists[c][i]
            if chanInfo['percentValid'] < percentValid:
                percentValid = chanInfo['percentValid']
        thisTileInfo['percentValid'] = percentValid
        if percentValid < 0.01:
            print 'Dropping empty tile' + thisTileInfo['prefix']
            continue

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
                                 thisTileInfo['pixelRow'], thisTileInfo['pixelCol'], forceFromHere)
    
        # Make a mask of valid pixels for this tile
        # - It would be great if the input images had a mask.
        # - Currently all-black pixels anywhere in the image get masked out!
        cmd = './makeSimpleImageMask ' + thisTileInfo['tileMaskPath'] +' '+ thisTileInfo['allTilesString']
        cmdRunner(cmd, thisTileInfo['tileMaskPath'], forceFromHere)
        
        key = thisTileInfo['prefix']
        tileDict[key] = thisTileInfo

    # Generate a "personalized" brightness file for each tile
    splitScaleBrightnessGains(brightnessGainsPath, tileDict, scaleRatio)

    
    # Now that we have all the per-tile info, compute the color transform for each tile.
    for tile in tileDict.itervalues():    
        
        try:
        
            # Generate the color pairs
            # - HRSC colors are written with the brightness correction already applied
            cmd = ('./writeHrscColorPairs ' + highResColorBaseTile +' '+ tile['allTilesStringAndMask']
                   +' '+ tile['spatialTransformPath'] +' '+ tile['brightnessGainsPath'] +' '+ tile['colorPairPath'])
            cmdRunner(cmd, tile['colorPairPath'], forceFromHere)
            
            #raise Exception('DEBUG')
            
            # Compute the color transform
            cmd = ('python /home/smcmich1/repo/MassUpload/solveHrscColor.py ' + tile['colorTransformPath']
                                                                              +' '+ tile['colorPairPath'])
            cmdRunner(cmd, tile['colorTransformPath'], forceFromHere)
            
        except CmdRunException: # Flag errors with this tile
            tile['stillValid'] = False

# TODO: Be robust to tile failures -> Some will definately lie outside the image!

    # Now that we have all the spatial transforms, generate the new color image for each tile.
    for tile in tileDict.itervalues():
        
        # New color generation is handled in this function
        generateNewHrscColor(tile, tileDict, forceFromHere)
        
        #raise Exception('DEBUG')

    #raise Exception('DEBUG')


    # Return all of the generated tile information!
    return tileDict

#-----------------------------------------------------------------------------------------
# Helper classes


   

def findOverlappingHrscImages(tileBounds):
    '''Return a list of all the HRSC image prefixes that overlap a given region'''
    
    # TODO: Search our database to perform this function!
    return ['h0022_0000']#,
            #'h0506_0000',
            #'h2411_0000',
            #'h6419_0000')


#-----------------------------------------------------------------------------------------


fullBasemapPath      = '/home/smcmich1/data/hrscMapTest/projection_space_basemap.tif'

sourceHrscFolder = '/home/smcmich1/data/hrscMapTest/external_data'

hrscOutputFolder = '/home/smcmich1/data/hrscMapTest/hrscFiles'



print 'Starting basemap enhancement script...'



#---------------------------
# Prep the base map
# - For now we extract an arbitrary chunk, later this will be tiled.
#
# Get the HRSC bounding box and expand it
#HRSC_BB_EXPAND_DEGREES = 1.5
#(minLon, maxLon, minLat, maxLat) = IrgGeoFunctions.getGeoTiffBoundingBox(hrscBasePathInList[0]+'_nd3.tif')
#minLon -= HRSC_BB_EXPAND_DEGREES
#maxLon += HRSC_BB_EXPAND_DEGREES
#minLat -= HRSC_BB_EXPAND_DEGREES
#maxLat += HRSC_BB_EXPAND_DEGREES


# Break up the base map into tiles

# DEBUG: Set a single high res fixed tile for testing!
# TODO:  Iterate through a set of tiles! Eventually all of them!

# TODO: Go down to 20 to 10 meters!
OUTPUT_RESOLUTION_METERS_PER_PIXEL = 100

basemapInstance = mosaicTileManager.MarsBasemap(fullBasemapPath, OUTPUT_RESOLUTION_METERS_PER_PIXEL)

# TODO: Move the tileIndex class?
thisTileIndex      = MosaicUtilities.TileIndex(88, 199) #basemapInstance.getTileIndex(98.5, -27.5)
thisTileBounds     = basemapInstance.getTileRectDegree(thisTileIndex)

print 'Tile index  = ' + str(thisTileIndex)
print 'Tile bounds = ' + str(thisTileBounds)

# Now that we have selected a tile, generate all of the tile images for it.
(smallTilePath, largeTilePath, grayTilePath, outputTilePath) = basemapInstance.generateTileImages(thisTileIndex)

## DEBUG Clear the existing output tile
#os.remove(outputTilePath)

# Now find all of the HRSC images that may overlap the tile.
hrscPrefixList = findOverlappingHrscImages(thisTileBounds)

#------------------------------------
# Process the individual HRSC images and add them to the mosaic

# TODO: This C++ program can do multiple tiles in one call.



for hrscPrefix in hrscPrefixList: # Loop through input HRSC images

    thisHrscPrefix = os.path.join(sourceHrscFolder, hrscPrefix)
    thisHrscFolder = os.path.join(hrscOutputFolder, hrscPrefix)

    #try:

    # Transform the HRSC image to the same projection/resolution as the upsampled base map crop
    #metersPerPixel = NOEL_MAP_METERS_PER_PIXEL / (RESOLUTION_INCREASE/100.0)
    #tileDict = generateHrscColorImage(basemapInstance, thisTileIndex, thisHrscPrefix, thisHrscFolder)

    # Fetch and preprocess the HRSC image
    # - TODO: Use the manager class to handle this!
    hrscObject = hrscImageManager.HrscImage(hrscPrefix, sourceHrscFolder, thisHrscFolder, basemapInstance)

    # Complete the high resolution components
    hrscObject.prepHighResolutionProducts()

    # Get the tile information from the HRSC image
    tileDict = hrscObject.getTileInfo(thisTileBounds)
    print 'Found these tile intersections:'
    for hrscTile in tileDict.itervalues():
        print hrscTile['prefix']

    #raise Exception('DEBUG')

    #except: # Testing registration
    #    continue

    # TODO: Choose tile order based on edge content?

    # For each tile...
    i = 0
    for hrscTile in tileDict.itervalues(): 
    
        #try:
        cmd = ('./hrscMosaic ' + outputTilePath +' '+ outputTilePath +' '+ hrscTile['newColorPath'] +' '+
                                  hrscTile['tileMaskPath'] +' '+ hrscTile['tileToTileTransformPath'])
        MosaicUtilities.cmdRunner(cmd, outputTilePath, True)
            
        #except CmdRunException:
        #    hrscTile['stillValid'] = False

        #i += 1
        #if i == 3:
        #    raise Exception('DEBUG')
    
        # TODO: Make sure the output tiles still have their geo info

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




