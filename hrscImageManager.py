



import os
import sys

import IrgGeoFunctions

import mosaicTileManager


# TODO: Move to a general file
def projCoordToPixelCoord(x, y, geoInfo):
    '''Converts from projected coordinates into pixel coordinates using
       the results from getImageGeoInfo()'''

    xOffset = x - geoInfo['projection_bounds'][0]
    yOffset = y - geoInfo['projection_bounds'][3]
    
    column = xOffset / geoInfo['pixel_size'][0]
    row    = yOffset / geoInfo['pixel_size'][1]
    
    return (column, row)



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



#-----------------------------------------------------------------------


# The order that HRSC channel images are stored and their names in a list
NUM_HRSC_CHANNELS = 5
HRSC_RED   = 0
HRSC_GREEN = 1
HRSC_BLUE  = 2
HRSC_NIR   = 3
HRSC_NADIR = 4
CHANNEL_STRINGS = ['red', 'green', 'blue', 'nir', 'nadir']

HRSC_HIGH_RES_TILE_SIZE = 1024

class HrscImage():
    '''
       Class to manage an input HRSC image.
       
       In addition to loading and tiling the HRSC image components,
       this class handles some low resolution adjustment to match the
       base map.
       
       This class does not do any operations that depend on specific
       high resolution basemap data.
    '''
    
    def __init__(self, setName, inputFolder, outputFolder, basemapInstance, force=False):
        '''Set up all the low resolution HRSC products.'''
        
        print 'Initializing HRSC image: ' + setName
        
        
        # Set up some paths
        self._hrscBasePathIn      = os.path.join(inputFolder,  setName)
        self._hrscBasePathOut     = os.path.join(outputFolder, setName)
        self._tileFolder          = self._hrscBasePathOut + '_tiles'
        self._lowResMaskPath      = self._hrscBasePathOut + '_low_res_mask.tif'
        self._brightnessGainsPath = self._hrscBasePathOut + '_brightness_gains.csv'
        self._basemapCropPath     = self._hrscBasePathOut + '_low_res_spatial_transform.csv' # A crop of the basemap used in several places
        self._colorPairPath       = self._hrscBasePathOut + '_low_res_color_pairs.csv'
        self._lowResSpatialRegistrationPath  = self._hrscBasePathOut + '_low_res_spatial_transform_basemap.csv'
        self._highResSpatialRegistrationPath = self._hrscBasePathOut + '_high_res_spatial_transform_basemap.csv'
        self._lowResSpatialCroppedRegistrationPath = self._hrscBasePathOut + '_low_res_cropped_spatial_transform.csv'
        
        
        # Get full list of input paths
        self._inputHrscPaths = _getHrscChannelPaths(self._hrscBasePathIn)
        
        # Compute the HRSC bounding box
        lowResNadirPath = self._lowResWarpedPaths[HRSC_NADIR]
        (minLon, maxLon, minLat, maxLat) = IrgGeoFunctions.getImageBoundingBox(lowResNadirPath)
        hrscBoundingBoxDegrees = Rectangle(minLon, maxLon, minLat, maxLat)
        
        
        # Record input parameters
        self._basemapMpp  = basemapInstance.getLowResMpp()
        self._outputMpp   = basemapInstance.getHighResMpp()
        self._proj4String = basemapInstance.getProj4String()
        self._basemapRegistrationPath = basemapInstance.getRegistrationFile() # Path to the grayscale low res entire base map
        
        # Cut out a region from the basemap around the location of the HRSC image
        CROP_BUFFER_LAT = 0.5
        CROP_BUFFER_LON = 0.5
        expandedBoundingBox = hrscBoundingBoxDegrees
        expandedBoundingBox.expand(CROP_BUFFER_LON, CROP_BUFFER_LAT)
        basemapInstance.getCroppedRegionDegrees(expandedBoundingBox, self._basemapCropPath)

        print 'Generating low res image copies...'
        
        # Generate a copy of each input HRSC channel at the low basemap resolution
        self._lowResWarpedPaths = [self._warpToProjection(path, outputFolder, '_basemap_res',
                                                          self._basemapMpp, force)
                                   for path in self._inputHrscPaths]
        
        # Build up a string containing all the low res paths for convenience
        self._lowResPathString = ''
        for path in self._lowResWarpedPaths:
            self._lowResPathString += path + ' '
            
        print 'Generating low resolution mask...'
            
        # Make a mask at the low resolution
        cmd = './makeSimpleImageMask ' + self._lowResMaskPath +' '+ self._lowResPathString
        cmdRunner(cmd, self._lowResMaskPath, force)            
        self._lowResPathStringAndMask = self._lowResPathString +' '+ self._lowResMaskPath
        self._lowResMaskImageSize = IrgGeoFunctions.getImageSize(self._lowResMaskPath)
        
        
        # Compute the spatial registration from the HRSC image to the base map
        self._computeSpatialRegistration(basemapInstance, lowResNadirPath, force)

        # Compute the spatial transform to the cropped region
        isHighRes = False
        basemapInstance.updateTransformToBoundsDegrees(self._lowResSpatialRegistrationPath, self._lowResSpatialCroppedRegistrationPath, expandedBoundingBox, isHighRes)
        
        # Compute the brightness scaling gains relative to the cropped base map
        # - This is done at low resolution
        # - The low resolution output is smoothed out later to avoid jagged edges.
        cmd = './computeBrightnessCorrection ' + self._basemapCropPath +' '+ lowResHrscPathString +' '+ self._lowResSpatialCroppedRegistrationPath +' '+ self._brightnessGainsPath
        cmdRunner(cmd, brightnessGainsPath, forceFromHere)        
        
        
        
        # Now we have done everything we plan to with the low resolution maps
        
        
    def prepHighResolutionProducts(force=False):
        '''Generates all of the high resolution HRSC products'''
        
        
        # Convert each high res input image into the output format
        
        # Generate a copy of each input HRSC channel at the high output resolution
        print 'Generating high resolution warped channel images...'
        self._highResWarpedPaths = [self._warpToProjection(path, outputFolder, '_output_res',
                                                          self._outputMpp, force)
                                   for path in self._inputHrscPaths]
        
        # Build up a string containing all the low res paths for convenience
        self._highResPathString = ''
        for path in self._highResWarpedPaths:
            self._highResPathString += path + ' '
    
        # Make a mask at the outputresolution
        print 'Generating high resolution mask...'
        cmd = './makeSimpleImageMask ' + self._highResMaskPath +' '+ self._highResPathString
        cmdRunner(cmd, self._highResMaskPath, force)            
        self._highResPathStringAndMask = self._highResPathString +' '+ self._highResMaskPath        
        
        
        
        
        # Split the image up into tiles at the full output resolution       
        # - There is one list of tiles per HRSC channel
        # - Each channel gets its own subfolder
        this._tileFolder = os.path.join(os.path.dirname(self._hrscBasePathOut), 'tiles') # TODO: Do we need no divide up the folders more?        
        tileInfoLists = [[], [], [], [], []] # One list per channel
        for c in range(NUM_HRSC_CHANNELS):
            # Get the info for this channel
            warpedPath    = self._highResWarpedPaths[c]
            channelString = CHANNEL_STRINGS[c]
            
            # Write all the tiles to a folder and get the info list for those tiles
            channelOutputFolder = os.path.join(this._tileFolder, channelString)
            if not os.path.exists(channelOutputFolder):
                os.mkdir(channelOutputFolder)
            tileInfoLists[c] = splitImage(imagePath, channelOutputFolder, HRSC_HIGH_RES_TILE_SIZE)
            
        # Verify that each channel generated the same number of tiles
        numTiles = len(tileInfoLists[0])
        for pathList in tileInfoLists:
            assert(len(pathList) == numTiles)
        print 'Generated ' +str(numTiles)+ ' tiles.'
        
        # Loop through each of the tiles we created and consolidate information across channels
        self._tileDict = {}
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
            filePrefix       = 'tile_' + thisTileInfo['prefix']
            tileInfoBasePath = os.path.join(this._tileFolder, filePrefix)
            thisTileInfo['colorPairPath'        ] = tileInfoBasePath+'_color_pairs.csv'
            thisTileInfo['colorTransformPath'   ] = tileInfoBasePath+'_color_transform.csv'
            thisTileInfo['newColorPath'         ] = tileInfoBasePath+'_new_color.tif'
            thisTileInfo['spatialTransformPath' ] = tileInfoBasePath+'_spatial_transform.csv'
            thisTileInfo['brightnessGainsPath'  ] = tileInfoBasePath+'_brightness_gains.csv'
            thisTileInfo['tileMaskPath'         ] = tileInfoBasePath, filePrefix+'_tile_mask.tif'
            thisTileInfo['allTilesStringAndMask'] = allTilesString + ' ' + thisTileInfo['tileMaskPath']
            
            thisTileInfo['stillValid'] = True # Set this to false if there is an error processing this tile
            
            # Generate individual tile versions of the spatial transform and brightness gains
            self._getHighResSpatialTransformForTile(thisTileInfo['pixelRow'], thisTileInfo['pixelCol'], thisTileInfo['spatialTransformPath'], force)
        
            # Make a mask of valid pixels for this tile
            # - It would be great if the input images had a mask.
            # - Currently all-black pixels anywhere in the image get masked out!
            cmd = './makeSimpleImageMask ' + thisTileInfo['tileMaskPath'] +' '+ thisTileInfo['allTilesString']
            cmdRunner(cmd, thisTileInfo['tileMaskPath'], force)
            
            key = thisTileInfo['prefix']
            self._tileDict[key] = thisTileInfo
    
                
            
            
            
        # Generate per-tile spatial transforms.
        
        # Generate per-tile brightness transforms.
        
        # Fill out tile information dictionary.
        
        
        # Generate the color pairs
        # - HRSC colors are written with the brightness correction already applied
        # - Color pairs are computed at the input basemap (low) resolution.
        cmd = ('./writeHrscColorPairs ' + self._basemapCropPath +' '+ self._lowResPathStringAndMask
               +' '+ self._lowResSpatialCroppedRegistrationPath +' '+ self._brightnessGainsPath +' '+ self._colorPairPath)
        cmdRunner(cmd, self._colorPairPath, forceFromHere)        
        
        
        
        TODO!!!
        
          
        
        
    def getTransformToTile():
        '''Computes the transform to a specific output tile'''
        pass
        ## Convert the transform to be relative to the output tile
        ## - Do this in both low (basemap) and high (output) resolution
        #tileBounds = basemapInstance.getLowResPixelBounds(tileIndex)
        #scaleSpatialTransform(basemapLowResTransformPath, outputPath, basemapInstance.getResolutionIncrease(),
        #                      -tileBounds.minX, -tileBounds.minY, force)
        #scaleSpatialTransform(basemapLowResTransformPath, outputPathLowRes, 1.0,
        #                      -tileBounds.minX, -tileBounds.minY, force)
        #        
        
        
        
    def _getHrscChannelPaths(hrscBasePath):
        '''Get the list of all HRSC channel paths from the base path'''
        # This skips the channels we are not interested in
        return (hrscBasePath+'_re3.tif',
                hrscBasePath+'_gr3.tif',
                hrscBasePath+'_bl3.tif',
                hrscBasePath+'_ir3.tif',
                hrscBasePath+'_nd3.tif')
    
    
    def _warpToProjection(sourcePath, outputFolder, postfix, proj4String, metersPerPixel, force=False):
        '''Warps an HRSC file into the specified projection space'''
        fileName   = sourcePath[sourcePath.rfind('/')+1:]
        warpedPath = os.path.join(outputFolder, fileName)[:-4] + postfix +'.tif'
        cmd = ('gdalwarp ' + sourcePath +' '+ warpedPath + ' -r cubicspline '
                 ' -t_srs "'+proj4String+'" -tr '
                 + str(metersPerPixel)+' '+str(metersPerPixel)+' -overwrite')
        cmdRunner(cmd, warpedPath, force)
        return warpedPath
            
        
    def _estimateRegistration(baseImage, otherImage, outputPath):
        '''Writes an estimated registration transform to a file based on geo metadata.'''
        # This function assumes the images are in the same projection system!
        
        # Get the projection bounds and size in both images
        baseGeoInfo     = IrgGeoFunctions.getImageGeoInfo(baseImage,  False)
        otherGeoInfo    = IrgGeoFunctions.getImageGeoInfo(otherImage, False)
        baseProjBounds  = baseGeoInfo[ 'projection_bounds']
        otherProjBounds = otherGeoInfo['projection_bounds']
        baseImageSize   = baseGeoInfo[ 'image_size']
        otherImageSize  = otherGeoInfo['image_size']
        
        #print baseGeoInfo
        #print '----'
        #print otherGeoInfo
        
        # Now estimate the bounding box of the other image in the base image
        topLeftCoord = projCoordToPixelCoord(otherProjBounds[0], otherProjBounds[3], baseGeoInfo)
        
        transform = SpatialTransform()
        transform.setShift(topLeftCoord[0], topLeftCoord[1])
        tranform.write(outputPath)
    
        return topLeftCoord    
    
    
    def _computeSpatialRegistration(basemapInstance, hrscPath, force=False):
        '''Compute the spatial registration from the HRSC image to the base image'''
    
        # Estimate the spatial transform using the image metadata
        # - This is always done against the full image basemap
        estimatedTransformPath = self._hrscBasePathOut + 'spatial_transform_basemap_estimated.csv'
        estX, estY = estimateRegistration(self._basemapRegistrationPath, hrscPath, estimatedTransformPath)
        
        #raise Exception('DEBUG')
    
        # TODO: Check the number of inliers!    
        # Refine the spatial transform using image data
        # - This is computed at the base map resolution and then scaled up to the output resolution
        cmd = ('./RegisterHrsc ' + self._basemapRegistrationPath +' '+ hrscPath
               +' '+ self._lowResSpatialRegistrationPath +' '+ str(1) +' '+ estimatedTransformPath)
        cmdRunner(cmd, self._lowResSpatialRegistrationPath, force)
        
        # We have the low res spatial transform, now compute the high res spatial transform to the base map.
        lowResSpatialTransform = SpatialTransform(self._lowResSpatialRegistrationPath)
        tx, ty = lowResSpatialTransform.getShift()
        self._boundingBoxLowResPixels = Rectangle(tx, tx+self._lowResMaskImageSize[0], ty, ty+self._lowResMaskImageSize[1])
        self._boundingBoxProjMeters   = basemapInstance.lowResPixelBoxToProjectedBox(self._boundingBoxLowResPixels)
        self._boundingBoxDegrees      = basemapInstance.projectedBoxToDegreeBox(self._boundingBoxProjMeters)
        
        
        self._highResSpatialRegistrationPath
        
        #raise Exception('DEBUG')

    #def _scaleSpatialTransform(inputPath, outputPath, scale, offsetX, offsetY, force=False):
    #    '''Scale up a spatial transform for a higher resolution'''
    #    
    #    if os.path.exists(outputPath) and (not force):
    #        return
    #    
    #    transform = SpatialTransform(inputPath)
    #    dx, dy = transform.getShift()
    #    dx = (dx + offsetX) * scale
    #    dy = (dy + offsetY) * scale
    #    transform.setShift(dx, dy)
    #    
    #    transform.write(outputPath)
        

    
    def _getHighResSpatialTransformForTile(self, row, col, outputPath, force):
        '''We already have the '''
        pass
    



def splitScaleBrightnessGains(fullPath, tileDict, scale):
    '''Generates a brightness gains file for a single tile'''

    print 'Generating split brightness gains...'
    
    # Read in the entire input file
    lowResVals = numpy.loadtxt(fullPath, skiprows=1, delimiter=',', usecols=(0,))   
    lowResRows = range(0,len(lowResVals))
    
    for tile in tileDict.itervalues():
        
        outputPath = tile['brightnessGainsPath']
        if os.path.exists(outputPath):
            continue
        
        # Compute the row range in the input tile
        tileHeight       = tile['height'  ]
        pixelRow         = tile['pixelRow']
        fullSizeStartRow = pixelRow
        fullSizeStopRow  = (pixelRow+tileHeight)
        
        # For each desired ouput value, compute the location in the input values
        thisTileRowsInInput = numpy.empty([tileHeight])
        index = 0
        for r in range(fullSizeStartRow, fullSizeStopRow):
            thisTileRowsInInput[index] = r / scale
            index += 1
        # Use numpy to interpolate values 
        thisTileVals = numpy.interp(thisTileRowsInInput, lowResRows, lowResVals)
        zeroCol      = [0 for i in thisTileVals]
            
        # Write out the interpolated values
        numpy.savetxt(outputPath, thisTileVals, header=str(tileHeight),  fmt='%1.6f, 0.0', comments='')  
    


def getAdjacentTiles(tile, tileDict):
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
    
    
    
    
    
class HrscImageManager():
    '''Class to manage which HRSC images are stored locally on disk'''
    pass