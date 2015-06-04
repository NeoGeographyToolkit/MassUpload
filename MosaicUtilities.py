



import os
import sys
import copy
import math
import subprocess

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


#def getImageSize(imagePath):
#    """Returns the size [rows, columns] in an image"""
#    # Make sure the input file exists
#    if not os.path.exists(imagePath):
#        raise Exception('Image file ' + imagePath + ' not found!')
#       
#    # Use subprocess to suppress the command output
#    cmd = ['gdalinfo', imagePath]
#    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
#    textOutput, err = p.communicate()
#
#    # Extract the size from the text
#    sizePos    = textOutput.find('Size is')
#    endPos     = textOutput.find('\n', sizePos+7)
#    sizeStr    = textOutput[sizePos+7:endPos]
#    sizeStrs   = sizeStr.strip().split(',')
#    cols = int(sizeStrs[0])
#    rows = int(sizeStrs[1])
#    
#    return (rows, cols)

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
    #print text
    
    # Output looks like this, but lines are missing if there are no pixels in them.
    #     257066: (  0,  0,  0) #000000 black
    #     317182: (255,255,255) #FFFFFF white

    lines = text.split('\n')
    if not ('black' in text): # No black pixels!
        return 0
    else:   
        blackLine = lines[0].strip()
        blackCount = int(blackLine[:blackLine.find(':')])
        return blackCount
        
    # Currently we don't care about the other count
    #otherLine = lines[1].strip()
    #otherCount = int(otherLine[:otherLine.find(':')])
    


# TODO: Stuff up here should come from common files
#=======================================================================================
#=======================================================================================



class TileIndex:
    '''Class used for indexing tiles'''
    def __init__(self, row=0, col=0):
        self.row = int(row)
        self.col = int(col)
    
    def __str__(self):
        return ('(R:%d, C:%d)' % (self.row, self.col))
    
    def getPostfix(self):
        '''Get a string usable as a postfix'''
        return ('%04d_%04d' % (self.row, self.col))
    
    
# TODO: Make changes so this can work with integers?
class Rectangle:
    '''Simple rectangle class for ROIs. Max values are NON-INCLUSIVE'''
    def __init__(self, minX=0, maxX=0, minY=0, maxY=0):
        self.minX = minX
        self.maxX = maxX
        self.minY = minY
        self.maxY = maxY
        #
        #if not self.hasArea(): # Debug helper
        #    print 'RECTANGLE WARNING: ' + str(self)
    
    def __str__(self):
        return ('minX: %f, maxX: %f, minY: %f, maxY: %f' % (self.minX, self.maxX, self.minY, self.maxY))
    
    def getBounds(self):
        '''Returns (minX, maxX, minY, maxY)'''
        return (self.minX, self.maxX, self.minY, self.maxY)
    
    def width(self):
        return self.maxX - self.minX
    def height(self):
        return self.maxY - self.minY
    
    def hasArea(self):
        '''Returns true if the rectangle contains any area.'''
        return (self.width() > 0) and (self.height() > 0)
    
    def perimiter(self):
        return 2*self.width() + 2*self.height()
    
    def area(self):
        '''Returns the valid area'''
        if not self.hasArea():
            return 0
        else:
            return self.height() * self.width()

    def getMinCoord(self):
        return (self.minX, self.minY)
    def getMaxCoord(self):
        return (self.maxX, self.maxY)
    
    def shift(self, dx, dy):
        '''Shifts the entire box'''
        self.minX += dx
        self.maxX += dx
        self.minY += dy
        self.maxY += dy
    
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

    def expandToContain(self, x, y):
        '''Expands the rect to contain the given point'''
        if x < self.minX: self.minX = x
        if x > self.maxX: self.maxX = x
        if y < self.minY: self.minY = y
        if y > self.maxY: self.maxY = y
        
    def getIntersection(self, otherRect):
        '''Returns the overlapping region of two rectangles'''
        overlap = Rectangle(max(self.minX, otherRect.minX),
                            min(self.maxX, otherRect.maxX),
                            max(self.minY, otherRect.minY),
                            min(self.maxY, otherRect.maxY))
        return overlap
        
    def overlaps(self, otherRect):
        '''Returns true if there is any overlap between this and another rectangle'''
        overlapArea = self.getIntersection(otherRect)
        return overlapArea.hasArea()
    

class Tiling:
    '''Sets up a tiling scheme'''
    
    def __init__(self, boundsRect, tileWidth, tileHeight):
        '''Init with the region to cover and the tile size.
           All units are the same as is used in boundsRect.'''
        self._bounds      = boundsRect
        self._numTileCols = math.ceil(boundsRect.width()  / tileWidth)
        self._numTileRows = math.ceil(boundsRect.height() / tileHeight)
        self._tileWidth   = tileWidth  # Nominal size, some tiles will be smaller.
        self._tileHeight  = tileHeight
        #print 'Init tiling = ' + self.__str__()
        
    def getNominalTileSize(self):
        '''Get the nominal tile size'''
        return (self._tileWidth, self.TileWidth)
        
    def getTileSize(self, tileIndex):
        '''Returns the actual size of a selected tile'''
        bb = self.getTileBounds(tileIndex)
        return (bb.width(), bb.height())
        
    def getTile(self, x, y):
        '''Computes the tile that contains the input location'''
        tileCol = (x - self._bounds.minX) / self._tileWidth
        tileRow = (y - self._bounds.minY) / self._tileHeight
        return TileIndex(tileRow, tileCol)
    
    def getIntersectingTiles(self, rect):
        '''Returns a rectangle containing all the intersecting tiles'''
        ti = self.getTile(rect.minX, rect.minY)
        tileRect = Rectangle(ti.col, ti.col, ti.row, ti.row)
        ti = self.getTile(rect.maxX, rect.minY);  tileRect.expandToContain(ti.col, ti.row);
        ti = self.getTile(rect.maxX, rect.maxY);  tileRect.expandToContain(ti.col, ti.row);
        ti = self.getTile(rect.minX, rect.maxY);  tileRect.expandToContain(ti.col, ti.row);
        tileRect.maxX += 1 # Until the rect is update for integers we need to do this!
        tileRect.maxY += 1
        return tileRect
    
    
    def getTileBounds(self, tileIndex):
        '''Returns the boundaries of a given tile'''
        # Compute the nominal bounds
        xStart = self._bounds.minX + tileIndex.col * self._tileWidth
        yStart = self._bounds.minY + tileIndex.row * self._tileHeight
        bb     = Rectangle(xStart, xStart+self._tileWidth,
                           yStart, yStart+self._tileHeight)
        # Restrict the bounds to the initialized boundary
        return bb.getIntersection(self._bounds)
    
    def __str__(self):
        return ('numTilesCols: %d, numTileRows: %d, bounds: %s' %
                (self._numTileCols, self._numTileRows, str(self._bounds)))

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
        return (self.values[2], self.values[5])
    
    def setShift(self, dx, dy):
        self.values[2] = dx
        self.values[5] = dy
    
    def addShift(self, dx, dy):
        self.values[2] += dx
        self.values[5] += dy
        
    def setScaling(self, scaling):
        '''The scaling from input units to output units'''
        self.values[0] = scaling
        self.values[4] = scaling
    
    def transform(self, x, y):
        '''Apply the transform to an input point'''
        xOut = self.values[0]*x + self.values[1]*y + self.values[2]
        yOut = self.values[3]*x + self.values[4]*y + self.values[5]
        return (xOut, yOut)
    
    def load(self, path):
        '''Read the transform from a file'''

        with open(path, 'r') as f:
            f.readline() # Skip the header
            self.values[0], self.values[1], self.values[2] = f.readline().strip().split(',')
            self.values[3], self.values[4], self.values[5] = f.readline().strip().split(',')
            self.values[6], self.values[7], self.values[8] = f.readline().strip().split(',')
        for i in range(0,len(self.values)): # Convert the strings to floats
            self.values[i] = float(self.values[i])

    def write(self, path):
        '''Save the transform to a file'''
        
        # Write the output file
        with open(path, 'w') as f:
            f.write('3, 3\n')
            f.write('%lf, %lf, %lf\n' % (self.values[0], self.values[1], self.values[2]))
            f.write('%lf, %lf, %lf\n' % (self.values[3], self.values[4], self.values[5]))
            f.write('%lf, %lf, %lf\n' % (self.values[6], self.values[7], self.values[8]))
        if not os.path.exists(path):
            raise Exception('Failed to create transform file ' + path)


def getTransformedBoundingBox(transform, rectangle):
        '''Compute a bounding box for a transformed rectangle'''
        x,y = transform(rectangle.minX,rectangle.minY);  bb = Rectangle(x,x,y,y)   # Init a rect at the first point
        x,y = transform(rectangle.maxX,rectangle.minY);  bb.expandToContain(x,y);  # Expand to contain all of the other points
        x,y = transform(rectangle.maxX,rectangle.maxY);  bb.expandToContain(x,y);
        x,y = transform(rectangle.minX,rectangle.maxY);  bb.expandToContain(x,y);
        return bb

class GeoReference:
    '''Handles GDC / projected space transforms.
       Currently only works for a simple, global Mercator transform.'''
    
    def __init__(self, degreesToMeters):
        self._degreesToMeters = degreesToMeters
        self._lonLatBounds    = Rectangle(-180, 180, -90, 90)
        self._projectionBounds = copy.copy(self._lonLatBounds)
        self._projectionBounds.scaleByConstant(degreesToMeters)
        print 'GeoReference init: ' + str(self)
    
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
        projRect = copy.copy(degreeRect)
        projRect.scaleByConstant(self._degreesToMeters)
        return projRect
    
    def projectedRectToDegreeRect(self, projRect):
        '''Convert a bounding box in degrees to one in projected coordinates'''
        degreeRect = copy.copy(projRect)
        degreeRect.scaleByConstant(1.0/self._degreesToMeters)
        return degreeRect
    
    def getLonLatBounds(self):
        return self._lotLatBounds
 
    def getProjectionBounds(self):
        return self._projectionBounds

    def __str__(self):
        return ('lonlatBounds: %s, projectionBounds: %s') % (str(self._lonLatBounds), str(self._projectionBounds))
    

class ImageCoverage:
    '''Handles XY space to image space transforms'''
    
    def __init__(self, numCols, numRows, geoBounds):
        self._geoBounds = geoBounds
        self._numCols   = numCols
        self._numRows   = numRows
        self._metersPerPixelX = self._geoBounds.width()  / numCols
        self._metersPerPixelY = self._geoBounds.height() / numRows
        print 'ImageCoverage init: ' + str(self)
        
    def __str__(self):
        return ('geoBounds: %s, numCols: %d, numRows: %d, mppX: %lf, mppy: %lf'
                   % (str(self._geoBounds), self._numCols, self._numRows,
                      self._metersPerPixelX, self._metersPerPixelY))
    
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
        (projMinX, projMaxY) = self.pixelToProjected(pixelRect.minX, pixelRect.minY)
        (projMaxX, projMinY) = self.pixelToProjected(pixelRect.maxX, pixelRect.maxY)
        return Rectangle(projMinX, projMaxX, projMinY, projMaxY)

    # TODO: Handle rounding!
    def projectedRectToPixelRect(self, projectedRect):
        '''From a projected ROI, compute the pixel bounding box.'''
        (pixelMinX, pixelMinY) = self.projectedToPixel(projectedRect.minX, projectedRect.maxY)
        (pixelMaxX, pixelMaxY) = self.projectedToPixel(projectedRect.maxX, projectedRect.minY)
        return Rectangle(pixelMinX, pixelMaxX, pixelMinY, pixelMaxY)


class ImageWithGeoRef(GeoReference, ImageCoverage):
    '''Adds an image to a GeoReference'''

    def __init__(self, degreesToMeters, numCols, numRows):
        
        # Iniatialize with the GeoReference and the image size
        # - GeoReference happens first so we can call getProjectionBounds()
        GeoReference.__init__(self, degreesToMeters)
        ImageCoverage.__init__(self, numCols, numRows, self.getProjectionBounds())

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
        '''These are the total image height/width, not per tile.'''
        # This class requires that the pixel dimensions work out exactly!
        ImageWithGeoRef.__init__(self, degreesToMeters, numCols, numRows) # Init this first!
        pixelBounds  = Rectangle(0, numCols, 0, numRows)
        tileWidth    = numCols / numTileCols
        tileHeight   = numRows / numTileRows
        self._tiling = Tiling(pixelBounds, tileWidth, tileHeight) # These are total pixel sizes

    def getTileSizePixels(self):
        return self._tiling.getNominalTileSize()

        
    def getTileAtLonLat(self, lon, lat):
        '''Return the tile containing a given lon/lat location'''
        col, row = self.degreesToPixel(lon, lat)
        return self._tiling.getTile(col, row)

    def getTileRectPixel(self, tile):
        '''Returns the boundaries of a given tile in pixels'''
        return self._tiling.getTileBounds(tile)
     
    def getTileRectDegree(self, tile):
        '''Returns the boundaries of a given tile in degrees'''
        pixelRect = self.getTileRectPixel(tile)
        return self.pixelRectToDegreeRect(pixelRect)
        
    def getIntersectingTiles(self, rectDegrees):
        '''Returns a Rectangle containing all the tiles intersecting the input ROI'''
        # Convert to pixels, then just use the tiling function.
        rectPixels = self.degreeRectToPixelRect(rectDegrees)
        return self._tiling.getIntersectingTiles(rectPixels)
    


    
    