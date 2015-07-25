

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <HrscCommon.h>


//=============================================================


bool loadInputImages(int argc, char** argv, cv::Mat &basemapImage, std::vector<cv::Mat> &hrscChannels, 
                     cv::Mat &hrscMask, cv::Mat &transform, BrightnessCorrector &corrector,
                     std::string &outputPath)
{
  std::vector<std::string> hrscPaths(NUM_HRSC_CHANNELS);
  std::string baseImagePath = argv[1];
  hrscPaths[0] = argv[2]; // R
  hrscPaths[1] = argv[3]; // G
  hrscPaths[2] = argv[4]; // B
  hrscPaths[3] = argv[5]; // NIR
  hrscPaths[4] = argv[6]; // NADIR
  std::string hrscMaskPath         = argv[7];
  std::string spatialTransformPath = argv[8];
  std::string brightnessPath       = argv[9];
  outputPath = argv[10];
  

  const int LOAD_GRAY = 0;
  const int LOAD_RGB  = 1;
  
  // Load the base map
  if (!readOpenCvImage(baseImagePath, basemapImage, LOAD_RGB))
    return false;

  // Load all of the HRSC images
  for (size_t i=0; i<NUM_HRSC_CHANNELS; ++i)
  {
    if (!readOpenCvImage(hrscPaths[i], hrscChannels[i], LOAD_GRAY))
      return false;
  }

  // Load the HRSC mask
  if (!readOpenCvImage(hrscMaskPath, hrscMask, LOAD_GRAY))
      return false;
  
  // Load the spatial transform
  if (!readTransform(spatialTransformPath, transform))
    return false;

  // Load brightness correction data
  if (!corrector.readProfileCorrection(brightnessPath))
    return false;

  //std::cout << "DONE LOADING" << std::endl;
  return true;
}

/// Generate a list of matched pixels for the base map and HRSC images
/// --> Format is "baseR, baseG, baseB, R, G, B, NIR, NADIR" 
int writeColorPairs(const cv::Mat &basemapImage, const cv::Mat &spatialTransform,
                     const BrightnessCorrector &corrector,
                     const std::vector<cv::Mat> hrscChannels, const cv::Mat &hrscMask,
                     const std::string outputPath,
                     const int sampleDist = 25)
{
  // Control what percentage of the pixel pairs we use
  const int numRows = hrscMask.rows;
  const int numCols = hrscMask.cols;

  std::ofstream outputFile;

  // Iterate over the pixels of the HRSC image
  bool gotValue;
  cv::Vec3b baseValues;
  size_t numPairs = 0;
  for (int r=0; r<numRows; r+=sampleDist)
  {
    for (int c=0; c<numCols; c+=sampleDist)
    {         
      // Skip masked out HRSC pixels
      if (hrscMask.at<MASK_DATA_TYPE>(r,c) == 0)
        continue;
    
      // Compute the equivalent location in the basemap image
      float baseX, baseY;
      affineTransform(spatialTransform, c, r, baseX, baseY);
      
      // Extract all of the basemap values at that location
      baseValues = interpPixelRgb(basemapImage, baseX, baseY, gotValue);
      if (!gotValue) 
        continue; // Failed to interpolate a value here
      
      // Write all the values to a line in the file
      if (numPairs == 0)
        outputFile.open(outputPath.c_str(), std::ofstream::out);
      
      for (size_t i=0; i<NUM_BASE_CHANNELS; ++i)
      {
        outputFile << static_cast<int>(baseValues[i]) <<", "; // Cast to int so we don't print ASCII
      }
      // The HRSC values have a brightness correction applied
      for (size_t i=0; i<NUM_HRSC_CHANNELS-1; ++i)
      {
        outputFile << static_cast<int>(corrector.correctPixel(hrscChannels[i].at<unsigned char>(r,c), r))  <<", ";
      }
      outputFile << static_cast<int>(corrector.correctPixel(hrscChannels[NUM_HRSC_CHANNELS-1].at<unsigned char>(r,c), r)) << std::endl;
      numPairs++;
  
    } // End col loop
  } // End row loop
  
  if (numPairs > 0)
    outputFile.close();

  //std::cout << "numPairs = " << numPairs << std::endl;
  return numPairs;
}

//============================================================================


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 11)
  {
    printf("usage: WriteColorPairs <Base Image Path> <HRSC Red> <HRSC Green> <HRSC Blue> <HRSC NIR> <HRSC Nadir> <HRSC Mask> <Transform File Path> <Brightness File Path> <Output Path>\n");
    return -1;
  }
  
  //printf("Loading input images...\n");
  cv::Mat hrscMask, basemapImage, spatialTransform, gain;
  std::string outputPath;
  std::vector<cv::Mat> hrscChannels(NUM_HRSC_CHANNELS);
  BrightnessCorrector corrector;
  if (!loadInputImages(argc, argv, basemapImage, hrscChannels, hrscMask, spatialTransform, corrector, outputPath))
    return -1;

  // TODO: The spatial transform should be from HRSC to BASEMAP

  // Generate the list of color pairs
  //printf("Writing pixel pairs...\n");
  const int MIN_PIXEL_PAIRS = 80;
  int sampleDist = 64;
  int numSamples = 0;
  while (sampleDist > 0)
  {
    numSamples = writeColorPairs(basemapImage, spatialTransform, corrector, hrscChannels, hrscMask, outputPath, sampleDist);
        
    if (numSamples >= MIN_PIXEL_PAIRS) // We got enough samples, we are finished!
        return 0; 
    // Otherwise reduce the sample distance 
    sampleDist /= 4;
    printf("Trying again with sample distance = %d\n", sampleDist);
  }
  if (numSamples > 0) // At least we got something, for now we just use this.
  {
    printf("Warning: Only found %d color pairs!\n", numSamples);
    return 0;
  }
  
  // Did not get any points even with a sample distance of 1!
  printf("Failed to detect any color pairs!!!\n");
  return -1;
}


