

#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <HrscCommon.h>

//=============================================================


bool loadInputData(int argc, char** argv,
                   std::vector<cv::Mat> &inputImages,
                   std::string &outputPath)
{
  if (argc < 3)
  {
    printf("Not enough input arguments passed in!\n");
    return false;
  }
  outputPath = argv[1];
  const int numInputImages = argc - 2;
  inputImages.resize(numInputImages);
  
  const int LOAD_GRAY = 0;
  const int LOAD_RGB  = 1;
  
  // Load all of the images
  for (int i=0; i<numInputImages; ++i)
  {
    std::string thisPath = argv[i+2];
    if (!readOpenCvImage(thisPath, inputImages[i], LOAD_GRAY))
      return false;
  }

  return true;
}


/// Make an output mask from a set of input images
bool generateSimpleImageMask(const std::vector<cv::Mat> &inputImages, cv::Mat &outputImage)
{
  const size_t numBands = inputImages.size();
    
  // Initialize the output image
  const size_t numRows = inputImages[0].rows;
  const size_t numCols = inputImages[0].cols;
  outputImage = cv::Mat(numRows, numCols, CV_8UC1);
 
  // Iterate over the pixels of the HRSC image
  for (int r=0; r<numRows; r+=1)
  {
    for (int c=0; c<numCols; c+=1)
    {
        
      // Check to see if this pixel should be masked out   
      // - Apply the brightness correction while we are at it.
      int count = 0;
      for (size_t i=0; i<numBands; ++i)
      {
        if (inputImages[i].at<unsigned char>(r,c) > 0)
          ++count;
      }
      // For HRSC we need all pixel values to be present!
      // - TODO: Some sort of BB system so internal pixels are not masked!
      if (count < numBands) 
        outputImage.at<unsigned char>(r, c) = 0;
      else
        outputImage.at<unsigned char>(r, c) = 255;

    } // End col loop
  } // End row loop

  return true;
}

//============================================================================


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc < 3)
  {
    printf("usage: makeSimpleImageMask <output path> <input image path>...\n");
    return -1;
  }
  
  printf("Loading input data...\n");
  
  // Load the input images  
  std::vector<cv::Mat> hrscChannels;
  std::string outputPath;
  if (!loadInputData(argc, argv, hrscChannels, outputPath))
    return -1;

  printf("Making output mask...\n");
  
  // Build the output mask
  cv::Mat outputMask;
  generateSimpleImageMask(hrscChannels, outputMask);

    
  // Write the output image
  cv::imwrite(outputPath, outputMask);

  return 0;
}


