

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <HrscCommon.h>


//=============================================================


bool loadInputImages(int argc, char** argv, cv::Mat &basemapImage, std::vector<cv::Mat> &hrscChannels, 
                     cv::Mat &transform, BrightnessCorrector &corrector, std::string &outputPath)
{
  std::vector<std::string> hrscPaths(NUM_HRSC_CHANNELS);
  std::string baseImagePath = argv[1];
  hrscPaths[0] = argv[2]; // R
  hrscPaths[1] = argv[3]; // G
  hrscPaths[2] = argv[4]; // B
  hrscPaths[3] = argv[5]; // NIR
  hrscPaths[4] = argv[6]; // NADIR
  std::string spatialTransformPath = argv[7];
  std::string brightnessPath = argv[8];
  outputPath  = argv[9];
  

  const int LOAD_GRAY = 0;
  const int LOAD_RGB  = 1;
  
  // Load the base map
  basemapImage = cv::imread(baseImagePath, LOAD_RGB);
  if (!basemapImage.data)
  {
    printf("Failed to load base map image!\n");
    return false;
  }

  // Load all of the HRSC images
  for (size_t i=0; i<NUM_HRSC_CHANNELS; ++i)
  {
    hrscChannels[i] = cv::imread(hrscPaths[i], LOAD_GRAY);
    if (!hrscChannels[i].data)
    {
      printf("Failed to load HRSC image: %s\n", hrscPaths[i].c_str());
      return false;
    }
  }

  // Load the spatial transform
  if (!readTransform(spatialTransformPath, transform))
    return false;

  // Load brightness correction data
  if (!corrector.readProfileCorrection(brightnessPath))
    return false;

  return true;
}

/// Generate a list of matched pixels for the base map and HRSC images
/// --> Format is "baseR, baseG, baseB, R, G, B, NIR, NADIR" 
bool writeColorPairs(const cv::Mat &basemapImage, const cv::Mat &spatialTransform, const BrightnessCorrector &corrector,
                     const std::vector<cv::Mat> hrscChannels, const std::string outputPath)
{
  // Control what percentage of the pixel pairs we use
  const int SAMPLE_DIST = 25;
 

  std::ofstream outputFile(outputPath.c_str());

  // Iterate over the pixels of the HRSC image
  bool gotValue;
  cv::Vec3b baseValues;
  for (int r=0; r<hrscChannels[0].rows; r+=SAMPLE_DIST)
  {
  
    for (int c=0; c<hrscChannels[0].cols; c+=SAMPLE_DIST)
    {     
      // TODO: Do we need to check all channels?
      // Skip empty pixels in the HRSC images
      if (hrscChannels[0].at<unsigned char>(r,c) == 0)
        continue;
    
      // Compute the equivalent location in the basemap image
      float baseX, baseY;
      affineTransform(spatialTransform, c, r, baseX, baseY);
      
      // Extract all of the basemap values at that location
      baseValues = interpPixelRgb(basemapImage, baseX, baseY, gotValue);
           
      if (gotValue) // Write all the values to a line in the file
      {
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
      }
    } // End col loop
  } // End row loop
  
  outputFile.close();

  return true;
}

//============================================================================


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 10)
  {
    printf("usage: WriteColorPairs <Base Image Path> <HRSC Red> <HRSC Green> <HRSC Blue> <HRSC NIR> <HRSC Nadir> <Transform File Path> <Brightness File Path> <Output Path>\n");
    return -1;
  }
  
  printf("Loading input images...\n");
  cv::Mat basemapImage, spatialTransform, gain;
  std::string outputPath;
  std::vector<cv::Mat> hrscChannels(NUM_HRSC_CHANNELS);
  BrightnessCorrector corrector;
  if (!loadInputImages(argc, argv, basemapImage, hrscChannels, spatialTransform, corrector, outputPath))
    return -1;

  // TODO: The spatial transform should be from HRSC to BASEMAP

  // Generate the list of color pairs
  printf("Writing pixel pairs...\n");
  writeColorPairs(basemapImage, spatialTransform, corrector, hrscChannels, outputPath);

  return 0;
}


