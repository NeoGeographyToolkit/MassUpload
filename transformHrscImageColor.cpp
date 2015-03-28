

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <HrscCommon.h>


//=============================================================


bool loadInputImages(int argc, char** argv, cv::Mat &basemapImage, std::vector<cv::Mat> &hrscChannels, 
                     cv::Mat &transform, std::string &outputPath)
{
  std::vector<std::string> hrscPaths(NUM_HRSC_CHANNELS);
  std::string baseImagePath = argv[1];
  hrscPaths[0] = argv[2]; // R
  hrscPaths[1] = argv[3]; // G
  hrscPaths[2] = argv[4]; // B
  hrscPaths[3] = argv[5]; // NIR
  hrscPaths[4] = argv[6]; // NADIR
  std::string colorTransformPath = argv[7];
  outputPath  = argv[8];
  

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
  if (!readTransform(colorTransformPath, transform))
    return false;

  return true;
}

/// Apply a color transform matrix to the HRSC bands.
/// - This is a 5x3 matrix.
bool transformHrscColor(const cv::Mat &basemapImage, const cv::Mat &colorTransform, const std::vector<cv::Mat> hrscChannels,
                        cv::Mat &outputImage)
{
  // Initialize the output image
  const size_t numRows = hrscChannels[0].rows;
  const size_t numCols = hrscChannels[0].cols;
  outputImage = cv::Mat(numRows, numCols, CV_8UC3);

  // Iterate over the pixels of the HRSC image
  cv::Vec3b outputPixel;
  std::vector<unsigned char> hrscPixel(NUM_HRSC_CHANNELS);
  for (int r=0; r<numRows; r+=1)
  {
    for (int c=0; c<numCols; c+=1)
    {  
      // Check to see if this pixel should be masked out   
      bool maskOut = false;
      for (int i=0; i<NUM_HRSC_CHANNELS; ++i)
      {
        hrscPixel[i] = hrscChannels[i].at<unsigned char>(r,c);
        if (hrscPixel[i] == 0)
        {
          maskOut = true;
          break;
        }
      }
      // If any of input HRSC pixels are black, the output pixel is black.
      // TODO: A smarter mask method!
      if (maskOut)
      {
        outputImage.at<cv::Vec3b>(r, c) = cv::Vec3b(0,0,0);
        continue;
      }
    
      // Compute the output pixel values      
      for (int j=0; j<NUM_BASE_CHANNELS; ++j)
      {
        outputPixel[j] = 0;
        float temp = 0.0;
        for (int i=0; i<NUM_HRSC_CHANNELS; ++i)
           temp += static_cast<float>(hrscPixel[i])*colorTransform.at<float>(i,j);
         if (temp < 0.0)
           temp = 0;
         if (temp > 255.0)
           temp = 255.0;
        outputPixel[j] = static_cast<unsigned char>(temp);
      }
      // Store in output image
      outputImage.at<cv::Vec3b>(r, c) = outputPixel;

    } // End col loop
  } // End row loop

  return true;
}

//============================================================================


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc < 9)
  {
    printf("usage: transformHrscImageColor <Base Image Path> <HRSC Red> <HRSC Green> <HRSC Blue> <HRSC NIR> <HRSC Nadir> <Color Transform File Path> <Output Path>\n");
    return -1;
  }
  
  // Load the input images  
  cv::Mat basemapImage, colorTransform;
  std::string outputPath;
  std::vector<cv::Mat> hrscChannels(NUM_HRSC_CHANNELS);
  if (!loadInputImages(argc, argv, basemapImage, hrscChannels, colorTransform, outputPath))
    return -1;

  // The color transform is from HRSC to the basemap.

  // Generate the transformed color image
  cv::Mat newColorImage;
  transformHrscColor(basemapImage, colorTransform, hrscChannels, newColorImage);


  // TODO: Apply unsharp masking to the color image

  
  // Write the output image
  cv::imwrite(outputPath, newColorImage);

  return 0;
}


