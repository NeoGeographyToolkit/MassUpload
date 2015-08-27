

import os
import sys
import numpy

'''
This file computes an HRSC -> RGB transform given pixel pairs from
the HRSC image and the input RGB basemap.

Each tile has its own computed transform.

Two transforms are computed: The primary transform, and a pure
intensity scale transform.  When applied, the pixel first has
the HRSC -> RGB transform applied, then the intensity is substituted
with the scaled HRSC NADIR intensity using a YCbCr transform.

This transform method gets a close match to the basemap color while
preserving the high quality of the NADIR channel.  The other channels
are more prone to artifacts which become visible if they influence the
intensity channel.

'''

def rgb2ycbcr(rgb):
    '''Converts a single RGB pixel to YCbCr'''    
    ycbcr = [0, 0, 0]
    ycbcr[0] =       0.299   *rgb[0] + 0.587   *rgb[1] + 0.114   *rgb[2]
    ycbcr[1] = 128 - 0.168736*rgb[0] - 0.331264*rgb[1] + 0.5     *rgb[2]
    ycbcr[2] = 128 + 0.5     *rgb[0] - 0.418688*rgb[1] - 0.081312*rgb[2]

    for i in range(0,3):
        if ycbcr[i] < 0:
            ycbcr[i] = 0
        if ycbcr[i] > 255:
            ycbcr[i] = 255

    return ycbcr
    
def ycbcr2rgb(ycbcr):
    '''Converts a single YCbCr pixel to RGB'''
    rgb = [0, 0, 0]
    rgb[0] = ycbcr[0]                              + 1.402   * (ycbcr[2] - 128)
    rgb[1] = ycbcr[0] - 0.34414 * (ycbcr[1] - 128) - 0.71414 * (ycbcr[2] - 128)
    rgb[2] = ycbcr[0] + 1.772   * (ycbcr[1] - 128)

    for i in range(0,3):
        if rgb[i] < 0:
            rgb[i] = 0
        if rgb[i] > 255:
            rgb[i] = 255

    return rgb


def solveTransform(inputPathList, outputPath):
    '''Solve for a 5x3 transform to convert HRSC to RGB'''

    # Load the input data
    # --> Format is "baseR, baseG, baseB, R, G, B, NIR, NADIR" 
    print 'Loading input files'
    basePixelList = []
    hrscPixelList = []
    
    basePixelListY = []
    hrscPixelListY = []
    
    for path in inputPathList: # Loop through input files

        # Loop through the lines in the file
        fileHandle = open(path, 'r')
        for line in fileHandle:
            # Break the line up into two seperate strings
            parts = line.strip().split(',')
            basePixel = [int(parts[0]), int(parts[1]), int(parts[2])]
            hrscPixel = [int(parts[3]), int(parts[4]), int(parts[5]), int(parts[6]), int(parts[7])]
            
            # Convert the base pixel from RGB to YCbCr
            basePixelYCC = rgb2ycbcr(basePixel)
              
            # Build a list of lists
            basePixelList.append(basePixel) # RGB
            hrscPixelList.append(hrscPixel)

            # Keep a seperate list of just intensity values
            basePixelListY.append([basePixelYCC[0]]) # Y
            hrscPixelListY.append([hrscPixel[4]]) # Nadir
            
        fileHandle.close()

    targets  = numpy.matrix(basePixelList) # Strip last semicolons
    inputs   = numpy.matrix(hrscPixelList)
    targetsY = numpy.matrix(basePixelListY) # Strip last semicolons
    inputsY  = numpy.matrix(hrscPixelListY)

    if (inputs.shape[1] == 0) or (inputs.shape[0] != targets.shape[0]):
        print 'ERROR--->'
        print targets.shape
        print inputs.shape
        print inputPathList
        print outputPath
        raise Exception('Error: Invalid input pixel pair files!')


    print 'Computing transform...'

    #print 'inputs'
    #print inputs
    #print 'targets'
    #print targets

    # Solve ax = b (5x3 matrix)
    transform, residuals, rank, s = numpy.linalg.lstsq(inputs, targets)
    print '--- Transform ---'
    print transform
    #print '--- Residuals ---'
    #print residuals

    # Compute Y scaling (just one number)
    print 'Computing Y Scaling...'
    transformY, residualsY, rank, s = numpy.linalg.lstsq(inputsY, targetsY)
    #print 'YYYYY'
    #print numpy.mean(inputsY)
    #print numpy.mean(targetsY)

    numRows = transform.shape[0] + 1
    numCols = transform.shape[1]

    #print 'Show results:'
    #print inputs*transform

    #print 'Show diff:'
    #diff = inputs*transform - targets
    #print diff  
    #print 'Mean abs(diff)'
    #print numpy.mean(abs(diff), axis=0)

    # Write the output transform to file
    f = open(outputPath, 'w')
    f.write(str(numRows) +', '+ str(numCols) + '\n')
    for r in range(numRows-1):
        #print str(transform[r])
        for c in range(numCols-1):
            f.write(str(transform[r,c]) + ', ')
        f.write(str(transform[r,numCols-1]) + '\n')
    f.write(str(float(transformY)) + ', 0.0, 0.0') # Extra row to store Y scaling
    f.close()

    print 'Finished writing out the transform file.'

def main():

    # Check the input
    if len(sys.argv) < 3:
      raise Exception('Must pass in the output path and the path to at least one pixel pair file!')

    # Read all the input paths
    inputPaths = []
    numInputFiles = len(sys.argv) - 2
    for i in range(numInputFiles):
      inputPaths.append(sys.argv[i+2])

    outputPath = sys.argv[1]

    solveTransform(inputPaths, outputPath)

if __name__ == "__main__":
    sys.exit(main())


