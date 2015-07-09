

import os
import sys
import numpy


def solveTransform(inputPathList, outputPath):
    '''Solve for a 5x3 transform to convert HRSC to RGB'''

    # Load the input data
    # --> Format is "baseR, baseG, baseB, R, G, B, NIR, NADIR" 
    print 'Loading input files'
    basePixelList = []
    hrscPixelList = []
    for path in inputPathList: # Loop through input files

        # Loop through the lines in the file
        fileHandle = open(path, 'r')
        for line in fileHandle:
            # Break the line up into two seperate strings
            parts = line.strip().split(',')
            basePixel = [int(parts[0]), int(parts[1]), int(parts[2])]
            hrscPixel = [int(parts[3]), int(parts[4]), int(parts[5]), int(parts[6]), int(parts[7])]
              
            # Build a list of lists
            basePixelList.append(basePixel)
            hrscPixelList.append(hrscPixel)
        fileHandle.close()

    targets  = numpy.matrix(basePixelList) # Strip last semicolons
    inputs = numpy.matrix(hrscPixelList)

    if (inputs.shape[1] == 0) or (inputs.shape[0] != targets.shape[0]):
        print targets.shape
        print inputs.shape
        raise Exception('Error: Invalid input pixel pair files!')


    print 'Computing transform...'

    # Solve ax = b
    transform, residuals, rank, s = numpy.linalg.lstsq(inputs, targets)

    #transform = transform.transpose

    print transform

    print 'Residuals:'
    print residuals

    numRows = transform.shape[0]
    numCols = transform.shape[1]


    # Write the output transform to file
    f = open(outputPath, 'w')
    f.write(str(numRows) +', '+ str(numCols) + '\n')
    for r in range(numRows):
        print str(transform[r])
        for c in range(numCols-1):
            f.write(str(transform[r,c]) + ', ')
        f.write(str(transform[r,numCols-1]) + '\n')
    f.close()


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


