



import os
import sys
import shutil
import logging

import copyGeoTiffInfo
import MosaicUtilities


class MarsBasemap:
    '''
       Class to manage the tiles of the input basemap and output image.
       
       Since the output image is extremely large it is broken up into a large number
       of reasonably sized tiles.
       
       The top left tile is index 0,0.
    '''
    
    def __init__(self, fullBasemapPath, outputTileFolder, backupFolder, center180=False):
        
        # Info about Noel's base map
        DEGREES_TO_PROJECTION_METERS = 59274.9
        FULL_BASEMAP_HEIGHT = 5760  # In pixels, low resolution.
        FULL_BASEMAP_WIDTH  = 11520
        
        BASEMAP_TILE_HEIGHT = 45 # In pixels, chosen to divide evenly.
        BASEMAP_TILE_WIDTH  = 45
        # --> Basemap size is 256 x 128 tiles
        
        #self.NOEL_MAP_METERS_PER_PIXEL = 1852.340625 # TODO: Make sure this is accurate before reprojecting everything        

        # Derived output parameters
        self.resolutionIncrease = 128
        outputHeight = FULL_BASEMAP_HEIGHT*self.resolutionIncrease
        outputWidth  = FULL_BASEMAP_WIDTH *self.resolutionIncrease
        numTileRows  = FULL_BASEMAP_HEIGHT / BASEMAP_TILE_HEIGHT
        numTileCols  = FULL_BASEMAP_WIDTH  / BASEMAP_TILE_WIDTH

        self._logger = logging.getLogger('mosaicTileManager')

        self._logger.info('Splitting basemap into %dx%d tiles for %d total tiles.'
                          % (numTileCols, numTileRows, numTileRows*numTileCols) )

        # Initialize two class instance to manage the coordinate systems
        # - The pixel operations are different but GDC/projection calls will return the same results.
        self._isCenter180 = center180
        self._lowResImage  = MosaicUtilities.TiledGeoRefImage(DEGREES_TO_PROJECTION_METERS, FULL_BASEMAP_WIDTH, FULL_BASEMAP_HEIGHT, numTileCols, numTileRows, center180)
        self._highResImage = MosaicUtilities.TiledGeoRefImage(DEGREES_TO_PROJECTION_METERS, outputWidth,        outputHeight,        numTileCols, numTileRows, center180)

        # Set up image products
        self.fullBasemapPath     = fullBasemapPath
        self.fullBasemapGrayPath = fullBasemapPath[:-4] + '_gray.tif'
        self._makeGrayBasemap()
        
        self._baseTileFolder = os.path.join(os.path.dirname(fullBasemapPath), 'basemap_tiles')
        if not os.path.exists(self._baseTileFolder):
            os.mkdir(self._baseTileFolder)

        self._backupFolder = backupFolder
        if not os.path.exists(self._backupFolder):
            os.mkdir(self._backupFolder)
            
        self._outputTileFolder = outputTileFolder
        if not os.path.exists(outputTileFolder):
            os.mkdir(outputTileFolder)
            
        # Create the main log file
        mainLogPath = self.getMainLogPath()
        cmd = 'touch ' + mainLogPath
        os.system(cmd)
        
    def getMainLogPath(self):
        return os.path.join(self._outputTileFolder, 'main_log.txt')

    def copySupportFilesFromBackupDir(self):
        '''Moves any needed files from the backup folder to the output folder'''
        backupMainLogPath = os.path.join(self._backupFolder, 'main_log.txt')
        cmd = 'rsync --update ' + backupMainLogPath +' '+ self.getMainLogPath()
        MosaicUtilities.cmdRunner(cmd, self.getMainLogPath(), True)

    def _makeGrayBasemap(self):
        '''Creates a grayscale version of the basemap if it does not already exist'''
        if not os.path.exists(self.fullBasemapGrayPath):
            cmd = MosaicUtilities.GDAL_TRANSLATE_PATH+' -b 1 ' + self.fullBasemapPath +' '+ self.fullBasemapGrayPath
            MosaicUtilities.cmdRunner(cmd, self.fullBasemapGrayPath, False)

    #------------------------------------------------------------
    # Helper functions

    def getBackupFolder(self):
        return self._backupFolder

    def getColorBasemapPath(self):
        '''Get the path to the original full basemap image'''
        return self.fullBasemapPath

    #def getGrayBasemapPath(self):
    #    '''Get the path to the grayscale full basemap image'''
    #    return self.fullBasemapGrayPath

    def getLowResMpp(self):
        return self._lowResImage.getMetersPerPixelX() # Same in both dimensions
    
    def getHighResMpp(self):
        return self._highResImage.getMetersPerPixelX() # Same in both dimensions
    
    def getResolutionIncrease(self):
        return self.resolutionIncrease
    
    def getProj4String(self):
        '''This is the projection system used for the global Mars map'''
        if self._isCenter180:
            return "+proj=eqc +lon_0=180 +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m"
        else: # Center 0
            return "+proj=eqc +lon_0=0   +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m"

    def getTileRectDegree(self, tileIndex):
        '''Get the bounding box of a tile in degrees'''
        return self._lowResImage.getTileRectDegree(tileIndex)

    def degreeRoiToPixelRoi(self, roi, isHighRes=False):
        if isHighRes:
            return self._highResImage.degreeRectToPixelRect(roi)
        else:
            return self._lowResImage.degreeRectToPixelRect(roi)
        
    def pixelRoiToDegreeRoi(self, roi, isHighRes=False):
        if isHighRes:
            return self._highResImage.pixelRectToDegreeRect(roi)
        else:
            return self._lowResImage.pixelRectToDegreeRect(roi)

    #------------------------------------------------------------
    # Tile creation functions
    
    def makeCroppedRegionProjMeters(self, boundingBoxProj, outputPath, force=False):
        '''Crops out a region of the original basemap image.'''
        
        (minX, maxX, minY, maxY) = boundingBoxProj.getBounds()
        if self._isCenter180: # Hack because the 180 version uses shifted projection coordinates!
            minX -= 10669477.100
            maxX -= 10669477.100
        projCoordString = '%f %f %f %f' % (minX, maxY, maxX, minY)
        cmd = (MosaicUtilities.GDAL_TRANSLATE_PATH+' '+ self.fullBasemapPath +' '+ outputPath
                                 +' -projwin '+ projCoordString)
        MosaicUtilities.cmdRunner(cmd, outputPath, force)
        
    def makeCroppedRegionDegrees(self, boundingBoxDegrees, outputPath, force=False):
        '''Crops out a region of the original basemap image.'''
        boundingBoxProj = self._highResImage.degreeRectToProjectedRect(boundingBoxDegrees)
        self.makeCroppedRegionProjMeters(boundingBoxProj, outputPath, force)
    

   
    #------------------------------------------------------------
    # Tile utility functions
    
    def getTileFolder(self, tileIndex):
        '''Get the folder for storing a given tile'''
        tilename = 'tile_' + tileIndex.getPostfix()
        return os.path.join(self._baseTileFolder, tilename)
    
    def getOutputTilePath(self, tileIndex):
        return os.path.join(self.getTileFolder(tileIndex), 'output_tile.tif')
    
    def lowResTransformToHighRes(self, baseTransformPath):
        '''Loads a low res transform and converts it to high res'''
        transform = SpatialTransform(baseTransformPath)
        dx, dy = transform.getShift()
        dx *= self._resolutionIncrease
        dy *= self._resolutionIncrease
        transform.setShift(dx, dy)
        return transform
    
    def convertPixelRoiResolution(self, pixelRectIn, inputIsHighRes=False):
        '''Converts a pixel ROI from the low to the high resolution image or vice versa'''
        
        if inputIsHighRes:
            degreeRect = self._highResImage.pixelRectToDegreeRect(pixelRectIn)
            return self._lowResImage.degreeRectToPixelRect(degreeRect)
        else:
            degreeRect = self._lowResImage.pixelRectToDegreeRect(pixelRectIn)
            return self._highResImage.degreeRectToPixelRect(degreeRect)
    
    def getIntersectingTiles(self, rectDegrees):
        '''Returns an iterator containing all the tiles which intersect the input rectangle'''
        return self._highResImage.getIntersectingTiles(rectDegrees)
    

    def generateMultipleTileImages(self, tileList, pool=None, force=False):
        '''Generate all the tile images in a range'''
        
        # Loop through all the tiles
        cmdList = []
        for tileIndex in tileList:
            
            tileBounds = self.getTileRectDegree(tileIndex)
            
            # Now that we have selected a tile, generate all of the tile images for it.
            # - The first time this is called for a tile it generates the backup image for the tile.
            (smallTilePath, largeTilePath, grayTilePath, outputTilePath, 
              tileLogPath, tileBackupPath, cmd1 , cmd2) =  \
                        self._generateImagesForTile(tileIndex, force=False)
            
            if not pool: # Go ahead and run the last two commands
                MosaicUtilities.cmdRunner(cmd1, grayTilePath, force)
                MosaicUtilities.cmdRunner(cmd2, outputTilePath, force)
            else: # Add the commands to a list
                # The seconds command seems prone to failure, so try it more than once.
                NUM_CONVERT_RETRIES = 5
                cmdList.append( (cmd1, grayTilePath,   force) )
                cmdList.append( (cmd2, outputTilePath, force, NUM_CONVERT_RETRIES) ) # Normal operation
                #cmdList.append( (cmd2, tileBackupPath, force, NUM_CONVERT_RETRIES) ) # BATCH GENERATE ONLY!
                
    
        if not pool: # No pool, we are finished.
            return True
        
        # Otherwise send all the commands to the processor.
        print 'Sending ' + str(len(cmdList)) + ' mosaic tile commands to the pool.'
        pool.map(MosaicUtilities.cmdRunnerWrapper, cmdList)
        return True
        
    def getPathsForTile(self, tileIndex):
        '''For a given tile index, returns some relevant file paths.'''
    
        tileFolder = self.getTileFolder(tileIndex)
        if not os.path.exists(tileFolder):
            os.mkdir(tileFolder)

        smallTilePath  = os.path.join(tileFolder, 'basemap_orig_res.tif')
        grayTilePath   = os.path.join(tileFolder, 'basemap_orig_res_gray.tif')
        largeTilePath  = os.path.join(tileFolder, 'basemap_output_res.tif') #DEFUNCT
        
        outputTileName = 'output_tile_'+tileIndex.getPostfix()+'.tif'
        outputTilePath = os.path.join(self._outputTileFolder, outputTileName)
        tileBackupPath = os.path.join(self._backupFolder,     outputTileName)
        tileLogPath    = os.path.join(self._outputTileFolder, 'output_tile_'+tileIndex.getPostfix()+'_log.txt')
    
        return (smallTilePath, largeTilePath, grayTilePath, outputTilePath, tileLogPath, tileBackupPath)
    
    def _generateImagesForTile(self, tileIndex, force=False):
        '''Generate all the basemap sourced images for a tile and return paths.
           Also generate a backup copy of the output tile if it does not already exist.
           In order to facilitate processing pool usage, this returns the last two
           commands that need to be executed.'''

        # Retrieve the needed paths
        (smallTilePath, largeTilePath, grayTilePath, outputTilePath, tileLogPath, tileBackupPath) = \
            self.getPathsForTile(tileIndex)        

        degreeRoi = self.getTileRectDegree(tileIndex)
        self._logger.info('MosaicTileManager: Generating tile images for region: ' + str(degreeRoi))

        # Crop out the section of the original base map for this tile
        self.makeCroppedRegionDegrees(degreeRoi, smallTilePath)

        # Generate a grayscale version of the small copy of this tile
        cmd1 = (MosaicUtilities.GDAL_TRANSLATE_PATH+' -b 1 ' + smallTilePath +' '+ grayTilePath)

        # Generate a copy of this tile at the full output resolution
        # - The image is blurred as it is upsampled so it does not look pixelated
        # - This operation is expensive and only needs to happen once per tile.
        # - All future tile updates will be pasted on top of this tile.
        if not os.path.exists(tileBackupPath):
            self._logger.warning('MosaicTileManager: Generating missing base tile: ' + tileBackupPath)
            raise Exception('There should not be any tiles missing from the backup folder!\n'
                            +tileBackupPath) # Turn on to prevent unexpected tile generation!
            if os.path.exists(outputTilePath):
                raise Exception('Output tile should never exist without backup file!')
            cmd2 = ('/byss/smcmich1/programs/ImageMagick-6.9.1_install/bin/convert -monitor -define filter:blur=0.88 -filter quadratic -resize ' 
                   + str(self.resolutionIncrease*100)+'% ' + smallTilePath +' '+ tileBackupPath)
        else: # Just make this a dummy command
            cmd2 = ':'
        # If this tile does not yet exist in the output folder, copy the latest backup there
        if not os.path.exists(outputTilePath):
            cmd2 += ' && cp ' + tileBackupPath +' '+ outputTilePath

        # Create the empty tile log file
        cmd = 'touch ' + tileLogPath
        os.system(cmd)
        
        return (smallTilePath, largeTilePath, grayTilePath, outputTilePath, 
                    tileLogPath, tileBackupPath, cmd1, cmd2)
        
    
    def checkLog(self, logPath, name):
        '''Return True if the name exists in the log file'''
        with open(logPath, 'r') as f:
            for line in f:
                if name in line: # Additional information is ignored
                    return True
        return False
    
    def updateLog(self, logPath, name, extra=None):
        '''Update a log file with a name'''
        with open(logPath, 'a+') as f:
            # Check if we have already written this HRSC file here
            for line in f:
                if name in line:
                    print 'WARNING: ' + name + ' added to tile but is already present!'
            # Append the hrsc set name to the tile log
            if extra: # Append additional information
                f.write(name + ' - ' + extra + '\n')
            else:
                f.write(name + '\n')
    
    
    
    
