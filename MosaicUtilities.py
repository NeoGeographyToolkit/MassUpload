



import os
import sys


# TODO: Organize further!



#----------------------------------------------------------------------------
# Functions



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
#=======================================================================================
#=======================================================================================



class TileIndex:
    '''Class used for indexing tiles'''
    def __init__(self, row=0, col=0):
        self.row = row
        self.col = col
    
    def __str__(self):
        return ('(R:%d, C:%d)' % (self.row, self.col))
    
    def getPostfix(self):
        '''Get a string usable as a postfix'''
        return ('%04d_%04d' % (self.row, self.col))
    
    
class Rectangle:
    '''Simple rectangle class for ROIs. Max values are NON-INCLUSIVE'''
    def __init__(self, minX=0, maxX=0, minY=0, maxY=0):
        self.minX = minX
        self.maxX = maxX
        self.minY = minY
        self.maxY = maxY
    
    def __str__(self):
        return ('minX: %f, maxX: %f, minY: %f, maxY: %f' % (self.minX, self.maxX, self.minY, self.maxY))
    
    def getBounds(self):
        '''Returns (minX, maxX, minY, maxY)'''
        return (self.minX, self.maxX, self.minY, self.maxY)
    
    def width(self):
        return maxX - minX
    def height(self):
        return maxY - minY
    
    def getMinCoord(self):
        return (self.minX, self.minY)
    def getMaxCoord(self):
        return (self.maxX, self.maxY)
    
    def scaleByConstant(self, xScale, yScale=None):
        '''Scale the units by a constant'''
        if yScale == None:
            yScale = xScale
        self.minX *= xScale
        self.maxX *= xScale
        self.minY *= yScale
        self.maxY *= yScale
        
    def expand(self, left, down, right=None, up=None):
        '''Expand the box by an amount in each direction'''
        self.minX -= left
        self.minY -= down
        if right == None: # If right and up are not passed in, use left and down for both sides.
            right = left
        if up == None:
            up = down
        self.maxX += right
        self.maxY += up


class Tiling:
    '''Sets up a tiling scheme'''
    
    def __init__(self, numTileCols, numTileRows, width, height):
        
        self._numTileCols = numTileCols
        self._numTileRows = numTileRows
        self._tileWidth   = width  / numTileCols
        self._tileHeight  = height / numTileRows
        
    def getTileWidth(self):
        return self._tileWidth
    
    def getTileHeight(self):
        return self._tileHeight
        
    def getTile(self, x, y):
        '''Computes the tile that contains the input location'''
        tileRow = x / self._tileWidth
        tileCol = y / self._tileHeight
        return TileIndex(tileRow, tileCol)
     
    def getTileBounds(self, tile):
        '''Returns the boundaries of a given tile'''
        xStart = tile.col * self._tileWidth
        yStart = tile.row * self._tileHeight
        return Rectangle(xStart, xStart+self._tileWidth,
                         yStart, yStart+self._tileHeight)

class SpatialTransform:
    '''Class to represent a spatial transform in image coordinates.
       Currently it is just a translation.'''
    
    def __init__(self, path=None):
        self.values = [1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0]
        
        if path: # Load from file if one was passed in
            self.load(path)
        
    #def shift(self, dx, dy):
    #    '''Add a translation to the transform'''
    #    values[2] += dx
    #    values[5] += dy
    
    def getShift(self):
        return (values[2], values[5])
    
    def setShift(self, dx, dy):
        values[2] = dx
        values[5] = dy
    
    def addShift(self, dx, dy):
        values[2] += dx
        values[5] += dy
    
    def load(self, path):
        '''Read the transform from a file'''

        with open(path, 'r') as f:
            f.readline() # Skip the header
            values[0], values[1], values[2] = fIn.readline().strip().split(',')
            values[3], values[4], values[5] = fIn.readline().strip().split(',')
            values[6], values[7], values[8] = fIn.readline().strip().split(',')

    def write(self, path):
        '''Save the transform to a file'''
        
        # Write the output file
        with open(path, 'w') as f:
            f.write('3, 3\n')
            f.write('%lf, %lf, %lf\n' % (values[0], values[1], values[2]))
            f.write('%lf, %lf, %lf\n' % (values[3], values[4], values[5]))
            f.write('%lf, %lf, %lf\n' % (values[6], values[7], values[8]))
        if not os.path.exists(path):
            raise Exception('Failed to create transform file ' + path)


class GeoReference:
    '''Handles GDC / projected space transforms.
       Currently only works for a simple, global Mercator transform.'''
    
    def __init__(self, degreesToMeters):
        self._degreesToMeters = degreesToMeters
        self._lonLatBounds    = Rectangle(-180, 180, -90, 90)
        self._projectionBounds = self._lonLatBounds # TODO: Force copies!
        self._projectionBounds.scaleByConstant(degreesToMeters)
    
    def degreesToProjected(self, lon, lat):
        '''Given a (lon, lat) coordinate, convert to the projected coordinate system.'''
        return (lon * self._degreesToMeters,
                lat * self._degreesToMeters)
   
    def projectedToDegrees(self, projX, projY):
        '''Given a projected coordinate, returns (lon, lat) in degrees'''
        return (proxX / self._degreesToMeters,
                projY / self._degreesToMeters)
    
    def degreeRectToProjectedRect(self, degreeRect):
        '''Convert a bounding box in degrees to one in projected coordinates'''
        projRect = degreeRect
        projRect.scaleByConstant(self._degreesToMeters)
        return projRect
    
    def projectedRectToDegreeRect(self, degreeRect):
        '''Convert a bounding box in degrees to one in projected coordinates'''
        degreeRect = projRect
        degreeRect.scaleByConstant(1.0/self._degreesToMeters)
        return degreeRect
    
    def getLonLatBounds(self):
        return self._lotLatBounds
 
    def getProjectionBounds(self):
        return self._projectionBounds


class ImageCoverage:
    '''Handles XY space to image space transforms'''
    
    def __init(self, numCols, numRows, geoBounds):
        self._geoBounds = geoBounds
        self._numRows   = numRows
        self._numCols   = numCols
        self._metersPerPixelX = self._geoBounds.width()  / numCols
        self._metersPerPixelY = self._geoBounds.height() / numRows
        
        
    def numRows(self):
        return self._numRows
    def numCols(self):
        return self._numCols
    
    def getMetersPerPixelX(self):
        return self._metersPerPixelX
    def getMetersPerPixelY(self):
        return self._metersPerPixelY
    
    # Point conversion functions ----------------------------------------
    def projectedToPixel(self, projX, projY):
        '''Convert projected coordinates to pixels, low res'''
        return ( (projX - self._geoBounds.minX) / self._metersPerPixelX,
                 (self._geoBounds.maxY - projY) / self._metersPerPixelY)
    
    def pixelToProjected(self, col, row):
        '''Convert projected coordinates to pixels, low res'''
        return ( col*self._metersPerPixelX + self._geoBounds.minX,
                 self._geoBounds.maxY - row*self._metersPerPixelY)
    
    # ROI conversion functions ---------------------------------------------
    def pixelRectToProjectedRect(self, pixelRect):
        '''From a pixel ROI, computes the projected ROI.'''
        (projMinX, projMaxY) = pixelToProjected(pixelRect.minX, pixelRect.minY)
        (projMaxX, projMinY) = pixelToProjected(pixelRect.maxX, pixelRect.maxY)
        return Rectangle(projMinX, projMaxX, projMinY, projMaxY)

    def projectedRectToPixelRect(self, projectedRect):
        '''From a projected ROI, compute the pixel bounding box.'''
        (pixelMinX, pixelMinY) = projectedToPixel(projectedRect.minX, projectedRect.maxY)
        (pixelMaxX, pixelMaxY) = projectedToPixel(projectedRect.maxX, projectedRect.minY)
        return Rectangle(pixelMinX, pixelMaxX, pixelMinY, pixelMaxY)
    


class ImageWithGeoRef(GeoReference, ImageCoverage):
    '''Adds an image to a GeoReference'''

    def __init__(self, degreesToMeters, numRows, numCols):
        
        # Iniatialize with the GeoReference and the image size
        # - GeoReference happens first so we can call getProjectionBounds()
        GeoReference.__init__(self, degreesToMeters)
        ImageCoverage.__init(self, numRows, numCols, self.getProjectionBounds())

    # Condensed conversion functions
    
    def pixelToDegrees(self, col, row):
        projX, projY = pixelToProjected(col, row)
        return projectedToDegrees(projX, projY)

    def degreesToPixel(self, lon, lat):
        projX, projY = degreesToProjected(lon, lat)
        return projectedToPixel(projX, projY)
    
    def pixelRectToDegreeRect(self, pixelRect):
        projRect = self.pixelRectToProjectedRect(pixelRect)
        return self.projectedRectToDegreeRect(projRect)

    def degreeRectToPixelRect(self, degreeRect):
        projRect = self.degreeRectToProjectedRect(degreeRect)
        return self.projectedRectToPixelRect(projRect)




class TiledGeoRefImage(ImageWithGeoRef):
    '''ImageWithGeoRef with tiles added'''
    
    def __init__(self, degreesToMeters, numCols, numRows, numTileCols, numTileRows):
        GeoReference.__init__(degreesToMeters, numCols, numRows) # Init this first!
        self._tiling = Tiling(numTileCols, numTileRows, self.numCols(), self.numRows())


    def getTileSizePixels(self):
        return (self._tiling.getTileWidth(), self._tiling.getTileHeight())

        
    def getTileAtLonLat(self, lon, lat):
        '''Return the tile containing a given lon/lat location'''
        col, row = self.degreesToPixel(lon, lat)
        return self._tiling.getTile(col, row)

    def getTileRectPixel(self, tile):
        '''Returns the boundaries of a given tile in pixels'''
        return self._tiling.getTileBounds(tile)
     
    def getTileRectDegree(self, tile):
        '''Returns the boundaries of a given tile in degrees'''
        pixelRect = self.getTileRectPixels(tile)
        return self.pixelRectToDegreeRect(pixelRect)
        


    
    