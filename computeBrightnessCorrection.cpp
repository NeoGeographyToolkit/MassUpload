

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <HrscCommon.h>


//=============================================================


bool loadInputImages(int argc, char** argv, cv::Mat &basemapImage, std::vector<cv::Mat> &hrscChannels, 
                     cv::Mat &hrscMask, cv::Mat &transform, std::string &outputPath)
{
  std::vector<std::string> hrscPaths(NUM_HRSC_CHANNELS);
  std::string baseImagePath = argv[1];
  hrscPaths[0] = argv[2]; // R
  hrscPaths[1] = argv[3]; // G
  hrscPaths[2] = argv[4]; // B
  hrscPaths[3] = argv[5]; // NIR
  hrscPaths[4] = argv[6]; // NADIR
  std::string hrscMaskPath = argv[7];
  std::string spatialTransformPath = argv[8];
  outputPath  = argv[9];
  

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

  return true;
}

/// Takes the mean of an RGB image across its horizontal axis, leaving only a vertical line of values.
bool rgbVertProfile(const std::vector<cv::Mat> hrscChannels, const cv::Mat &spatialTransform, 
                    const cv::Mat &hrscMask,
                    const cv::Mat &inputImage, cv::Mat &outputProfile)
{ 
  // Set the ROI in the basemap according to the transformed HRSC footprint
  cv::Rect_<int> baseBounds = getboundsInOtherImage(inputImage, hrscChannels[0], spatialTransform);

  // Allocate the output storage
  outputProfile.create(baseBounds.height, 1, CV_32FC1);

  // Used to compute the matching location in the HRSC image
  const int colOffset = static_cast<int>(spatialTransform.at<float>(0, 2));
  const int rowOffset = static_cast<int>(spatialTransform.at<float>(1, 2));

  // Loop through the input image
  for (int r=baseBounds.y; r<baseBounds.y+baseBounds.height; r++)
  {
    int hrscRow = r - rowOffset;
    // Compute the mean grayscale value of each row
    float  thisRowMean = 0.0;
    size_t numUsedCols = 0;
    for (int c=baseBounds.x; c<baseBounds.x+baseBounds.width; c++)
    {      
      int hrscCol = c - colOffset;
      
      // Skip masked out HRSC pixels
      if (hrscMask.at<unsigned char>(hrscRow, hrscCol) == 0)
        continue;
    
      cv::Vec3b inputPixel = inputImage.at<cv::Vec3b>(r,c);
      thisRowMean += (inputPixel[0] + inputPixel[1] + inputPixel[2])/3.0;
      ++numUsedCols;
    
    } // End col loop
    if (numUsedCols == 0) // Avoid divide by zero
      thisRowMean = 128.0;
    else
      thisRowMean /= static_cast<float>(numUsedCols);
    
    // Store the mean value for this row in the output image
    int i = r - baseBounds.y;
    outputProfile.at<float>(i,0) = static_cast<float>(thisRowMean);
  } // End row loop
  
  return true;
}

/// Takes the mean of an HRSC image across its horizontal axis, leaving only a vertical line of values.
/// - The mean is taken across all HRSC channels.
bool hrscVertProfile(const std::vector<cv::Mat> hrscChannels, const cv::Mat &hrscMask,
                     cv::Mat &outputProfile)
{
  // Allocate the output storage
  const int numRows = hrscMask.rows;
  const int numCols = hrscMask.cols;
  outputProfile.create(numRows, 1, CV_32FC1);
  
  // Loop through the input image
  for (int r=0; r<numRows; r++)
  {
    // Compute the mean grayscale value of each row
    float thisRowMean = 0.0;
    size_t numColsUsed = 0;
    for (int c=0; c<numCols; c++)
    {     
      
      // Skip masked out HRSC pixels
      if (hrscMask.at<unsigned char>(r,c) == 0)
        continue;

      // For each pixel the brightness is the mean value across all channels
      float thisPixelSum = 0.0;
      for (int channel=0; channel<NUM_HRSC_CHANNELS; channel++)
      {         
         thisPixelSum += static_cast<float>(hrscChannels[channel].at<unsigned char>(r,c));
      }
      thisRowMean += thisPixelSum/static_cast<float>(NUM_HRSC_CHANNELS);
      ++numColsUsed;
    
    } // End col loop
    if (numColsUsed == 0) // Avoid divide by zero
      thisRowMean = 128.0;
    else
      thisRowMean /= static_cast<float>(numColsUsed);
    
    // Store the mean value for this row in the output image
    outputProfile.at<float>(r,0) = static_cast<float>(thisRowMean);
  } // End row loop
  
  return true; 
}

/// Compute a gain/offset to make the HRSC brightness match the basemap brightness
void computeGainOffsets(const cv::Mat &baseProfile, const cv::Mat &hrscProfile,
                        const cv::Mat &spatialTransform, BrightnessCorrector &brightness)
{
  const int baseHeight = baseProfile.rows;
  const int hrscHeight = hrscProfile.rows;
  
  //const int spatialOffset = static_cast<int>(spatialTransform.at<float>(1, 2));
  
  cv::Mat gains, offsets;
  gains.create(  hrscHeight, 1, CV_32FC1);
  offsets.create(hrscHeight, 1, CV_32FC1);
  
  float rowGain = (float)baseHeight / (float)hrscHeight;
  
  // We compute one gain/offset pair per HRSC row.
  for (int r=0; r<hrscHeight; ++r)
  {
    // Compute the row in the base profile, it should be close to r.
    int baseRow = static_cast<int>((float)r*rowGain);
    //int baseRow = r + spatialOffset;
    //printf("HRSC row = %d --> Base row = %d\n", r, baseRow);
    
    
    // For now we don't use the offset.
    float newGain = 0.0;
    if (hrscProfile.at<float>(r,0) > 0) // Avoid divide by zero
      newGain = baseProfile.at<float>(baseRow,0) / hrscProfile.at<float>(r,0);
    gains.at<  float>(r,0) = newGain;
    offsets.at<float>(r,0) = 0.0;
  }
  brightness.set(gains, offsets);
}



//============================================================================


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc < 10)
  {
    printf("usage: WriteColorPairs <Base Image Path> <HRSC Red> <HRSC Green> <HRSC Blue> <HRSC NIR> <HRSC Nadir> <HRSC Mask> <Transform File Path> <Output Path>\n");
    return -1;
  }
  
  printf("Loading input images...\n");
  cv::Mat basemapImage, spatialTransform, hrscMask;
  std::string outputPath;
  std::vector<cv::Mat> hrscChannels(NUM_HRSC_CHANNELS);
  if (!loadInputImages(argc, argv, basemapImage, hrscChannels, hrscMask, spatialTransform, outputPath))
    return -1;

  printf("Computing brightness profiles...\n");
  
  // Squish the images down to vertical profiles
  cv::Mat basemapProfile, hrscProfile;
  rgbVertProfile(hrscChannels, spatialTransform, hrscMask, basemapImage,  basemapProfile);
  hrscVertProfile(hrscChannels, hrscMask, hrscProfile);
  
  /*
  printf("\n\nInput basemap profile\n");
  for (size_t r=0; r<30; ++r)
  {
    std::cout << 
        basemapProfile.at<float>(r+1400+static_cast<int>(spatialTransform.at<float>(1, 2)),0) 
        << std::endl;
  }

  printf("\n\nInput HRSC profile\n");
  for (size_t r=0; r<30; ++r)
  {
    std::cout << hrscProfile.at<float>(r+1400,0) << std::endl;
  }
  */
  
  printf("Smoothing profiles...\n");
  const int    PROFILE_SMOOTHING_SIZE = 31;
  const double SMOOTHING_SIGMA        = 21;  // TODO: Play with these values to improve quality
  cv::Mat basemapSmoothProfile, hrscSmoothProfile;
  
  // Apply a gaussian blur to the profiles
  cv::GaussianBlur(basemapProfile, basemapSmoothProfile, cv::Size(PROFILE_SMOOTHING_SIZE, PROFILE_SMOOTHING_SIZE), SMOOTHING_SIGMA, 0);
  cv::GaussianBlur(hrscProfile,    hrscSmoothProfile,    cv::Size(PROFILE_SMOOTHING_SIZE, PROFILE_SMOOTHING_SIZE), SMOOTHING_SIGMA, 0);

  // Compute a gain/offset to equalize the profiles
  printf("Computing gains...\n");
  BrightnessCorrector corrector;
  computeGainOffsets(basemapSmoothProfile, hrscSmoothProfile, spatialTransform, corrector);


  // Write the gain/offset numbers to a file.
  printf("Writing grayscale corrections...\n");
  corrector.writeProfileCorrection(outputPath);


  return 0;
}


