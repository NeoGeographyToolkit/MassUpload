

#include <stdio.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>


#include "opencv2/stitching/detail/autocalib.hpp"
#include "opencv2/stitching/detail/blenders.hpp"
#include "opencv2/stitching/detail/timelapsers.hpp"
#include "opencv2/stitching/detail/camera.hpp"
#include "opencv2/stitching/detail/exposure_compensate.hpp"
#include "opencv2/stitching/detail/matchers.hpp"
#include "opencv2/stitching/detail/motion_estimators.hpp"
#include "opencv2/stitching/detail/seam_finders.hpp"
#include "opencv2/stitching/detail/util.hpp"
#include "opencv2/stitching/detail/warpers.hpp"
#include "opencv2/stitching/warpers.hpp"


#include <HrscCommon.h>



/**
  Program to generate a mosaic by adding new images to a base image.
  
  Currently the images are just pasted on but we can do something
  fancier in the future.
 
 */

const int BLEND_DIST_GLOBAL = 64;

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
/*
// Without supporting classes this function is a mess
void getPasteBoundingBox(const cv::Mat &outputImage, const cv::Mat imageToAdd, const cv::Mat &spatialTransform,
                         int &minX, int &minY, int &maxX, int &maxY)
{
  // Init the bounds
  maxX = 0;
  maxY = 0;
  minX = outputImage.cols-1;
  minY = outputImage.rows-1;
  
  // Compute the four corners
  const int NUM_CORNERS = 4;
  float interpX[NUM_CORNERS], interpY[NUM_CORNERS];
  affineTransform(spatialTransform, 0,               0,               interpX[0], interpY[0]);
  affineTransform(spatialTransform, 0,               imageToAdd.rows, interpX[1], interpY[1]);
  affineTransform(spatialTransform, imageToAdd.cols, 0,               interpX[2], interpY[2]);
  affineTransform(spatialTransform, imageToAdd.cols, imageToAdd.rows, interpX[3], interpY[3]);
  
  // Adjust the bounds to match the corners
  for (int i=0; i<NUM_CORNERS; ++i)
  {
    if (floor(interpX[i]) < minX) minX = interpX[i];
    if (ceil( interpX[i]) > maxX) maxX = interpX[i];
    if (floor(interpY[i]) < minY) minY = interpY[i];
    if (ceil( interpY[i]) > maxY) maxY = interpY[i];
  }
}

// TODO: Need a bounding box class!

/// Just do a simple paste of one image on to another.
bool pasteImage(cv::Mat &outputImage,
                const cv::Mat &imageToAdd, const cv::Mat &imageMask, const cv::Mat &spatialTransform)
{
  const int tileSize = outputImage.rows; // Currently the code requires square tiles
    
  // Estimate the bounds of the new image so we do not have to iterate over the entire output image
  int minCol, minRow, maxCol, maxRow;
  cv::Mat newToOutput;
  cv::invert(spatialTransform, newToOutput);
  getPasteBoundingBox(outputImage, imageToAdd, newToOutput, minCol, minRow, maxCol, maxRow);
  
  // Restrict the paste ROI to valid bounds --> Should use a class function!
  if (minCol < 0)                minCol = 0;
  if (maxCol > outputImage.cols) maxCol = outputImage.cols;
  if (minRow < 0)                minRow = 0;
  if (maxRow > outputImage.rows) maxRow = outputImage.rows;
  //printf("minCol = %d, minRow = %d, maxCol = %d, maxRow = %d\n", minCol, minRow, maxCol, maxRow);
  
  // Iterate over the pixels of the output image
  bool gotValue;
  cv::Vec3b pastePixel;
  float interpX, interpY;
  for (int r=minRow; r<maxRow; r++)
  {
    for (int c=minCol; c<maxCol; c++)
    {        
      // Compute the equivalent location in the added image
      affineTransform(spatialTransform, c, r, interpX, interpY);
      
      // TODO: Don't use mirrored pixel over real pixel!
      
      // Extract all of the basemap values at that location.
      // - Call the mirror version of the function so we retain all edges.
      pastePixel = interpPixelMirrorRgb(imageToAdd, imageMask, interpX, interpY, gotValue);      
      
      // Skip masked pixels and out of bounds pixels
      if (!gotValue)
      {
        continue;
      }
      
      // If the interpolated pixel is good just overwrite the current value in the output image
      outputImage.at<cv::Vec3b>(r, c) = pastePixel;

    } // End col loop
  } // End row loop

  return true;
}


/// A weighted simple paste of one image on to another.
bool pasteWeightedImage(cv::Mat &outputImage,
                        const cv::Mat &imageToAdd, const cv::Mat &imageMask,
                        const cv::Mat &imageWeight, const cv::Mat &spatialTransform)
{
  const int tileSize = outputImage.rows; // Currently the code requires square tiles
    
  // Estimate the bounds of the new image so we do not have to iterate over the entire output image
  int minCol, minRow, maxCol, maxRow;
  cv::Mat newToOutput;
  cv::invert(spatialTransform, newToOutput);
  getPasteBoundingBox(outputImage, imageToAdd, newToOutput, minCol, minRow, maxCol, maxRow);
  
  // Restrict the paste ROI to valid bounds --> Should use a class function!
  if (minCol < 0)                minCol = 0;
  if (maxCol > outputImage.cols) maxCol = outputImage.cols;
  if (minRow < 0)                minRow = 0;
  if (maxRow > outputImage.rows) maxRow = outputImage.rows;
  //printf("minCol = %d, minRow = %d, maxCol = %d, maxRow = %d\n", minCol, minRow, maxCol, maxRow);
  
  // Iterate over the pixels of the output image
  bool gotValue;
  cv::Vec3b pastePixel, currentPixel;
  float interpX, interpY;
  const int NUM_RGB_CHANNELS = 3;
  for (int r=minRow; r<maxRow; r++)
  {
    for (int c=minCol; c<maxCol; c++)
    {        
      // Compute the equivalent location in the added image
      affineTransform(spatialTransform, c, r, interpX, interpY);
      
      // TODO: Don't use mirrored pixel over real pixel!
      
      // Extract all of the basemap values at that location.
      // - Call the mirror version of the function so we retain all edges.
      pastePixel = interpPixelMirrorRgb(imageToAdd, imageMask, interpX, interpY, gotValue);      
      
      // Skip masked pixels and out of bounds pixels
      if (!gotValue)
      {
        continue;
      }
      
      // If the interpolated pixel is good add the weighted value to the output image
      float weight = imageWeight.at<float>(r,c);
      currentPixel = outputImage.at<cv::Vec3b>(r, c);
      for (int chan=0; chan<NUM_RGB_CHANNELS; chan++)
      {
        float newVal = static_cast<float>(pastePixel[chan])*weight;
        currentPixel[chan] += static_cast<unsigned char>(newVal);
      }
      outputImage.at<cv::Vec3b>(r, c) = currentPixel;

    } // End col loop
  } // End row loop

  return true;
}


*/

/// A weighted simple paste of one image on to another.
bool pasteWeightedImage(cv::Mat &outputImage,
                         const cv::Mat &imageToAdd, const cv::Mat &imageMask,
                         const cv::Mat &imageWeight, const cv::Mat &spatialTransform)
{
  const int tileHeight = imageToAdd.rows;
  const int tileWidth  = imageToAdd.cols;
  
  // Use a simple integer translation to match the input masks!
  const int colOffset = static_cast<int>(spatialTransform.at<float>(0, 2));
  const int rowOffset = static_cast<int>(spatialTransform.at<float>(1, 2));
    
  // Compute the bounds of the new image so we do not have to iterate over the entire output image
  int minCol, minRow, maxCol, maxRow;
  minCol = -colOffset;
  minRow = -rowOffset;
  maxCol = minCol + tileWidth;
  maxRow = minRow + tileHeight;
  
  // Restrict the paste ROI to valid bounds --> Should use a class function!
  if (minCol < 0)                minCol = 0;
  if (maxCol > outputImage.cols) maxCol = outputImage.cols;
  if (minRow < 0)                minRow = 0;
  if (maxRow > outputImage.rows) maxRow = outputImage.rows;
  //printf("minCol = %d, minRow = %d, maxCol = %d, maxRow = %d\n", minCol, minRow, maxCol, maxRow);
  
  
  // TODO: Use OpenCV calls to replace this block of code!
  // TODO: Why does this produce better results than OpenCV?
  
  // Iterate over the pixels of the output image
  cv::Vec3b pastePixel, currentPixel;
  unsigned char maskPixel;
  const int NUM_RGB_CHANNELS = 3;
  for (int r=minRow; r<maxRow; r++)
  {
    for (int c=minCol; c<maxCol; c++)
    {
      const int smallR = r+rowOffset;
      const int smallC = c+colOffset;
      maskPixel = imageMask.at<unsigned char>(smallR, smallC);      
      
      if (maskPixel == 0) // Skip masked pixels and out of bounds pixels
        continue;

      pastePixel = imageToAdd.at<cv::Vec3b>(smallR, smallC);
      
      // If the interpolated pixel is good add the weighted value to the output image
      float weight = imageWeight.at<float>(smallR, smallC);
      currentPixel = outputImage.at<cv::Vec3b>(r, c);
      for (int chan=0; chan<NUM_RGB_CHANNELS; chan++)
      {
        float newVal = static_cast<float>(pastePixel[chan])*weight;
        currentPixel[chan] += static_cast<unsigned char>(newVal);
      }
      outputImage.at<cv::Vec3b>(r, c) = currentPixel;

    } // End col loop
  } // End row loop

  return true;
}









/// Sets up the base and paste masks so we get the desired images in the right places
void setImageMasks(const cv::Mat &baseMask,    const std::vector<cv::Mat> &pasteMasks,
                                               const std::vector<cv::Mat> &spatialTransforms,
                         cv::Mat &baseMaskOut)
{
  //  We want the base mask to be invalid underneath the paste mask except for the edges.
  //  Tile edges do not count as edges for this purpose.
  const int EDGE_SIZE = BLEND_DIST_GLOBAL;
  
  const size_t numMasks = pasteMasks.size();
  const cv::Rect baseRoi(0, 0, baseMask.cols, baseMask.rows);
  
  baseMask.copyTo(baseMaskOut); 
  
  cv::Mat kernel(EDGE_SIZE, EDGE_SIZE, CV_8UC1, 255);
  cv::Mat oneKernel(3, 3, CV_8UC1, 255);
  cv::Mat shrunkPasteMask, invertPasteMask, tempMat;
  for (size_t i=0; i<numMasks; ++i)
  {    
    // Generate a shrunk version to mask the base map for blending
    cv::erode(pasteMasks[i], shrunkPasteMask, kernel);
    
    // Subtract the the shrunk version from the base mask
    cv::absdiff(shrunkPasteMask, 255, invertPasteMask);
    
    
    // TODO: Improve image ROI handling so this does not break on certain images!
    
    // ROI of the small image in the base image
    cv::Rect pasteRoiInBase(static_cast<int>(-spatialTransforms[i].at<float>(0, 2)),
                            static_cast<int>(-spatialTransforms[i].at<float>(1, 2)),
                            invertPasteMask.cols, invertPasteMask.rows);
    
    // The portion of small ROI that is contained in the base image, base image coordinates.
    cv::Rect pasteRoiInBaseSafe = pasteRoiInBase & baseRoi;
    
    // The full ROI internal to the paste image
    cv::Rect pasteRoiInPaste(0, 0, invertPasteMask.cols, invertPasteMask.rows);
    // The full base ROI in the paste image
    cv::Rect baseRoiInPaste(static_cast<int>(spatialTransforms[i].at<float>(0, 2)),
                            static_cast<int>(spatialTransforms[i].at<float>(1, 2)),
                            baseMask.cols, baseMask.rows);
    // The portion of small ROI that is contained in the base image, paste image coordinates.
    cv::Rect pasteRoiInPasteSafe = pasteRoiInPaste & baseRoiInPaste;
    
    
    //std::cout << "Paste ROI: " << pasteRoiSafe << std::endl;
    baseMaskOut.copyTo(tempMat);
    cv::Mat outSection(baseMaskOut, pasteRoiInBaseSafe);
    cv::min(invertPasteMask(pasteRoiInPasteSafe), tempMat(pasteRoiInBaseSafe), outSection);
    
    // Expand base mask by a single pixel to help clean up strange
    //  pixel artifacts seen at the top and bottom of a pasted tile.
    cv::Mat tempBase = baseMaskOut;
    cv::dilate(tempBase, baseMaskOut, oneKernel);
    
    /*
    // DEBUG
    std::string path = "shrunkPasteMask" + itoa(i) + ".tif";
    cv::imwrite(path, shrunkPasteMask);
    path = "invertPasteMask" + itoa(i) + ".tif";
    cv::imwrite(path, invertPasteMask);
    path = "baseMaskOut_" + itoa(i) + ".tif";
    cv::imwrite(path, baseMaskOut);*/
  }
  
}

/// Paste new images using a graph cut to blend the seams.
bool pasteImagesGraphCut(const             cv::Mat  &baseImage,
                         const std::vector<cv::Mat> &pasteImages,
                         const std::vector<cv::Mat> &pasteMasks,
                         const std::vector<cv::Mat> &spatialTransforms,
                               std::vector<cv::Mat> &imageWeights,
                                           cv::Mat  &outputImage)
{
    const float FEATHER_BLEND_DIST = BLEND_DIST_GLOBAL;
    const int   NUM_BLEND_BANDS   = 1;
    
    // OpenCV provides the multi band blender and a feather blender
    // - For some reason the multi-band blender super smooths the entire image!
    const bool USE_MULTI_BLENDER = false;
    
    size_t numImages = pasteImages.size() + 1;
    
    // Set up the initial image masks
    // - Need to avoid feathering at the inter-tile boundaries.
    cv::Mat baseMaskTrue(baseImage.rows, baseImage.cols, CV_8UC1, 255);
    cv::Mat baseMaskShrunk;
    setImageMasks(baseMaskTrue, pasteMasks, spatialTransforms, baseMaskShrunk);
    
    printf("Converting data...\n");
    
    // Need to convert from Mat to UMat
    std::vector<cv::UMat >  umatImages(numImages);
    std::vector<cv::UMat >  seamMasks (numImages);
    std::vector<cv::Point> corners   (numImages);
    std::vector<cv::Size > sizes     (numImages);
    for (int i = 0; i < numImages-1; ++i)
    {
      pasteImages[i].convertTo(umatImages[i], CV_32F); // Seam Finder needs this format
      pasteMasks [i].copyTo(seamMasks [i]);
      
      // The input transform is just a translation in affine format
      // - TODO: Make this cleaner, get out inversion
      corners[i] = cv::Point(static_cast<int>(-spatialTransforms[i].at<float>(0, 2)),
                             static_cast<int>(-spatialTransforms[i].at<float>(1, 2)));
    }
    printf("Converting data 2...\n");
    // Add the base map to the list of input images with a mask
    const int baseIndex = numImages-1;
    baseImage.convertTo(umatImages[baseIndex], CV_32F);
    baseMaskShrunk.copyTo(seamMasks[baseIndex]);
    corners[baseIndex] = cv::Point(0,0);

    
    // Record the size of each input image
    for (int i = 0; i < numImages; ++i){
      sizes[i] = umatImages[i].size();
      std::cout << "sizeI = " << sizes[i] << std::endl;
      std::cout << "sizeM = " << seamMasks[i].size() << std::endl;
    }
/*
    int debugX = 1026;
    int debugY = 1972;
    std::cout << debugY+spatialTransforms[0].at<float>(1, 2) << ", "
              << debugX+spatialTransforms[0].at<float>(0, 2) << std::endl;
    std::cout << baseImage.at<cv::Vec3b>(debugY, debugX) << std::endl;
    std::cout << pasteImages[0].at<cv::Vec3b>(debugY+spatialTransforms[0].at<float>(1, 2),
                                              debugX+spatialTransforms[0].at<float>(0, 2)) << std::endl;
    */
    
    printf("Dumping input masks...\n");
      
    // Dump all the input masks to disk for debugging
    for (size_t i=0; i<numImages; ++i)
    {
      std::cout << "Corner: " << corners[i] << std::endl;
        
      std::string path = "pre_seam_mask" + itoa(i) + ".tif";
      cv::imwrite(path, seamMasks[i]);
      
      path = "image_in_" + itoa(i) + ".tif";
      cv::Mat temp;
      umatImages[i].convertTo(temp, CV_8U);
      cv::imwrite(path, temp);
    }
    
    // Initialize seam finder
    cv::Ptr<cv::detail::SeamFinder> seamFinder;
    
    // TODO: Pick one of these!
    //seamFinder = cv::makePtr<cv::detail::VoronoiSeamFinder>();
    //seamFinder = cv::makePtr<cv::detail::GraphCutSeamFinder>(cv::detail::GraphCutSeamFinderBase::COST_COLOR);
    float terminal_cost      = 100.f;//10000.f; // TODO Adjust these!
    float bad_region_penalty = 10.f;//1000.f;
    seamFinder = cv::makePtr<cv::detail::GraphCutSeamFinder>(cv::detail::GraphCutSeamFinderBase::COST_COLOR_GRAD,
                                                             terminal_cost, bad_region_penalty);
    if (!seamFinder)
    {
      std::cout << "Can't create the seam finder!'\n";
      return false;
    }

    // Call seam finder
    printf("Running seam finder...\n");
    //seamFinder->find(umatImages, corners, seamMasks); //TODO: This is wiping out all overlap!
    
    std::cout << "seamMasks type: " <<seamMasks[0].type() << std::endl;
    
    printf("Dumping output masks...\n");
    // Dump all the output masks to disk for debugging
    for (size_t i=0; i<numImages; ++i)
    {
      std::string path = "seam_mask" + itoa(i) + ".tif";
      cv::imwrite(path, seamMasks[i]);
    }
    
    // Initialize the blender
    printf("Initializing blender...\n");
    cv::Ptr<cv::detail::Blender> blender;
    bool TRY_GPU = false;
    
    if (USE_MULTI_BLENDER)
    {
      blender = cv::detail::Blender::createDefault(cv::detail::Blender::MULTI_BAND, TRY_GPU); // Feather or pyramid blend
      cv::detail::MultiBandBlender* mb = dynamic_cast<cv::detail::MultiBandBlender*>(blender.get());
      mb->setNumBands(NUM_BLEND_BANDS);
      //blender = cv::Ptr<cv::detail::Blender>(new cv::detail::MultiBandBlender(TRY_GPU, NUM_BLEND_BANDS));
      
    }
    else // Use the feather blender
    {
      blender = cv::detail::Blender::createDefault(cv::detail::Blender::FEATHER, TRY_GPU); // Feather or pyramid blend
      cv::detail::FeatherBlender* fb = dynamic_cast<cv::detail::FeatherBlender*>(blender.get());
      fb->setSharpness(1.0/FEATHER_BLEND_DIST);
    }
    
    
    blender->prepare(corners, sizes);
    
    // Feed all the tiles into the blender
    for (size_t i=0; i<numImages; ++i)
    {
      cv::Mat tempImg; // TODO: Have a seperate array of Mat images?
      umatImages[i].convertTo(tempImg, CV_16S); // Image must be CV_16SC3 and mask must be CV_8U
      //baseImage.convertTo(tempImg, CV_16S);
      blender->feed(tempImg, seamMasks[i], corners[i]);
    }
    
    if (!USE_MULTI_BLENDER)
    {
        std::vector<cv::UMat> imageWeightsUmat;
        cv::detail::FeatherBlender* fb = dynamic_cast<cv::detail::FeatherBlender*>(blender.get());
        cv::Rect r = fb->createWeightMaps(seamMasks, corners, imageWeightsUmat);
        //std::cout << "wmask type: " << imageWeightsUmat[0].type() << std::endl;
        
        imageWeights.resize(imageWeightsUmat.size());
        for (size_t i=0; i<imageWeightsUmat.size(); ++i)
        {
          imageWeightsUmat[i].convertTo(imageWeights[i], CV_32F);
            
          cv::Mat temp;
          imageWeightsUmat[i].convertTo(temp, CV_8U, 255.0);
          std::string path = "weight_mask_" + itoa(i) + ".tif";
          cv::imwrite(path, temp);
        }
    }
    /*
    // Blend the images
    printf("Running blender...\n");
    cv::Mat resultImage, resultMask;
    blender->blend(resultImage, resultMask);
    printf("Blended!\n");
    
    resultImage.convertTo(outputImage, CV_8U); // The result comes out as CV_16S
    //cv::imwrite("blendMask.tif", resultMask);
    */
    
    
    outputImage = cv::Mat::zeros(baseImage.rows, baseImage.cols, CV_8UC3);
    
    cv::Mat basemapMask(baseImage.rows, baseImage.cols, CV_8UC1, 255);
    cv::Mat basemapTransform = cv::Mat::eye(3, 3, CV_32F);
    pasteWeightedImage(outputImage, baseImage, baseMaskShrunk, imageWeights[numImages-1], basemapTransform);
  
    printf("Pasted the base.\n");
  
    // For now, just dump all of the HRSC images in one at a time.  
    for (size_t i=0; i<numImages-1; ++i)
    {
      pasteWeightedImage(outputImage, pasteImages[i], pasteMasks[i], imageWeights[i], spatialTransforms[i]);
    }    
    
    
    /*
    // Note: Images are stored in BGR format
    std::cout << "Output: " << outputImage.at<cv::Vec3b>(debugY, debugX) << std::endl;
    std::cout << "Output type: " <<outputImage.type() << std::endl;
    */
    cv::imwrite("blended.tif", outputImage);

   
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
  const size_t numHrscImages = hrscImages.size();

  // The spatial transform is from the base map to HRSC


  // Initialize the output image to be identical to the input basemap image
  

  printf("Pasting on HRSC images...\n");

  // Just use this function to compute the image weights
  // - The base image weights are in the last element
  std::vector<cv::Mat> imageWeights;
  cv::Mat outputImage, basemapFloat;
  pasteImagesGraphCut(basemapImage, hrscImages, hrscMasks, spatialTransforms, imageWeights, outputImage);
  
  /*
  cv::Mat outputImage =cv::Mat::zeros(basemapImage.rows, basemapImage.cols, CV_8UC3);
  
  printf("Ready to paste weighted...\n");
  
  cv::Mat basemapMask(basemapImage.rows, basemapImage.cols, CV_8UC1, 255);
  cv::Mat basemapTransform = cv::Mat::eye(3, 3, CV_32F);
  pasteWeightedImage(outputImage, basemapImage, basemapMask, imageWeights[numHrscImages], basemapTransform);
  
  printf("Pasted the base.\n");
  
  // For now, just dump all of the HRSC images in one at a time.
  cv::Mat outputImage(basemapImage);
  for (size_t i=0; i<numHrscImages; ++i)
  {
    //pasteWeightedImage(outputImage, hrscImages[i], hrscMasks[i], imageWeights[i], spatialTransforms[i]);
    pasteImage(outputImage, hrscImages[i], hrscMasks[i], spatialTransforms[i]);
  }
  */
  printf("Writing output file...\n");
  
  // Write the output image
  cv::imwrite(outputPath, outputImage);

  return 0;
}


