

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <HrscCommon.h>



/**
  Program to generate a mosaic by adding new images to a base image.
  
  Currently the images are just pasted on but we can do something
  fancier in the future.
 
 */

//=============================================================

/// Load all the input files
bool loadInputImages(int argc, char** argv, cv::Mat &basemapImage,
                     std::vector<cv::Mat> &hrscImages,
                     std::vector<cv::Mat> &hrscMasks,
                     std::vector<cv::Mat> &spatialTransforms, std::string &outputPath)
{
  // Parse the input arguments
  const size_t numHrscImages = (argc - 3)/3;
  std::vector<std::string> hrscPaths(numHrscImages),
                           hrscMaskPaths(numHrscImages),
                           spatialTransformPaths(numHrscImages);
  std::string baseImagePath = argv[1];
  outputPath = argv[2];
  
  for (size_t i=0; i<numHrscImages; ++i)
  {
    hrscPaths[i]             = argv[3 + 2*i];
    hrscMaskPaths[i]         = argv[4 + 2*i];
    spatialTransformPaths[i] = argv[5 + 2*i];
  }

  const int LOAD_GRAY = 0;
  const int LOAD_RGB  = 1;
  
  // Load the base map
  if (!readOpenCvImage(baseImagePath, basemapImage, LOAD_RGB))
      return false;
  
  // Load all of the HRSC images and their spatial transforms
  cv::Mat tempTransform;
  hrscImages.resize(numHrscImages);
  hrscMasks.resize(numHrscImages);
  spatialTransforms.resize(numHrscImages);
  for (size_t i=0; i<numHrscImages; ++i)
  {
    if (!readOpenCvImage(hrscPaths[i], hrscImages[i], LOAD_RGB))
      return false;
    if (!readOpenCvImage(hrscMaskPaths[i], hrscMasks[i], LOAD_GRAY))
    {
      printf("Mask read error!\n");
      return false;
    }
    
    if (!readTransform(spatialTransformPaths[i], tempTransform))
    {
      printf("Failed to load HRSC spatial transform: %s\n", spatialTransformPaths[i].c_str());
      return false;
    }
    // Each transform is read in HRSC_to_basemap but we want basemap_to_HRSC so invert.
    double check = cv::invert(tempTransform, spatialTransforms[i]);
  }
  printf("Loaded %d images.\n", numHrscImages);

  return true;
}

//TODO: We need a mask image to handle the borders properly!

/// Just do a simple paste of one image on to another.
bool pasteImage(cv::Mat &outputImage,
                const cv::Mat imageToAdd, const cv::Mat imageMask, const cv::Mat &spatialTransform)
{
  const int tileSize = outputImage.rows; // Currently the code requires square tiles
    
  // Iterate over the pixels of the output image
  bool gotValue;
  cv::Vec3b pastePixel;
  for (int r=0; r<outputImage.rows; r++)
  {
    for (int c=0; c<outputImage.cols; c++)
    {        
      // Compute the equivalent location in the basemap image
      float interpX, interpY;
      affineTransform(spatialTransform, c, r, interpX, interpY);
      
      
      // Extract all of the basemap values at that location.
      // - Call the mirror version of the function so we retain all edges.
      pastePixel = interpPixelMirrorRgb(imageToAdd, imageMask, interpX, interpY, gotValue);

      // Skip masked pixels and out of bounds pixels
      if (!gotValue)
        continue;

      // If the interpolated pixel is good just overwrite the current value in the output image
      outputImage.at<cv::Vec3b>(r, c) = pastePixel;

    } // End col loop
  } // End row loop

  return true;
}

//============================================================================


int main(int argc, char** argv)
{
  // Check input arguments
  if ((argc - 3) % 3 != 0 )
  {
    printf("usage: hrscMosaic <Base Image Path> <Output Path> [<Hrsc Rgb Path> <HrscMaskPath> <Spatial transform path>]... \n");
    return -1;
  }

  printf("Loading input data...\n");
  
  // Load the input images  
  cv::Mat basemapImage;
  std::string outputPath;
  std::vector<cv::Mat> hrscImages, hrscMasks, spatialTransforms;
  if (!loadInputImages(argc, argv, basemapImage, hrscImages, hrscMasks, spatialTransforms, outputPath))
  {
    printf("Error reading input arguments!\n");
    return -1;
  }

  // The spatial transform is from the base map to HRSC


  // Initialize the output image to be identical to the input basemap image
  cv::Mat outputImage(basemapImage);

  printf("Pasting on HRSC images...\n");

  // For now, just dump all of the HRSC images in one at a time.  
  const size_t numHrscImages = hrscImages.size();
  for (size_t i=0; i<numHrscImages; ++i)
  {
    pasteImage(outputImage, hrscImages[i], hrscMasks[i], spatialTransforms[i]);
  }

  printf("Writing output file...\n");
  
  // Write the output image
  cv::imwrite(outputPath, outputImage);

  return 0;
}


