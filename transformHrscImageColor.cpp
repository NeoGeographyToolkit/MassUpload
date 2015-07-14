

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <HrscCommon.h>


//=============================================================


bool loadInputData(int argc, char** argv,
                   std::vector<cv::Mat> &hrscChannels, cv::Mat &hrscMask,
                   BrightnessCorrector &corrector, std::string &outputPath,
                   cv::Mat &mainTransform, double &mainWeight,
                   std::vector<cv::Mat>   &otherTransforms,
                   std::vector<double>    &otherWeights,
                   std::vector<cv::Vec2i> &otherOffsets)
{
  if (argc < 11)
  {
    printf("Not enough input arguments passed in!\n");
    return false;
  }
    
  std::vector<std::string> hrscPaths(NUM_HRSC_CHANNELS);
  hrscPaths[0] = argv[1]; // R
  hrscPaths[1] = argv[2]; // G
  hrscPaths[2] = argv[3]; // B
  hrscPaths[3] = argv[4]; // NIR
  hrscPaths[4] = argv[5]; // NADIR
  std::string hrscMaskPath   = argv[6];
  std::string brightnessPath = argv[7];
  outputPath = argv[8];
  
  std::string mainTransformPath = argv[9];
  mainWeight = atof(argv[10]);
  
  const int LOAD_GRAY = 0;
  const int LOAD_RGB  = 1;
  
  
  // Load all of the HRSC images
  for (size_t i=0; i<NUM_HRSC_CHANNELS; ++i)
  {
    if (!readOpenCvImage(hrscPaths[i], hrscChannels[i], LOAD_GRAY))
      return false;
  }
  // Read in the mask
  if (!readOpenCvImage(hrscMaskPath, hrscMask, LOAD_GRAY))
      return false;

  // Load brightness correction data
  if (!corrector.readProfileCorrection(brightnessPath))
    return false;
  
  // Load the main spatial transform
  //std::cout << "Reading transform from: " << mainTransformPath << std::endl;
  if (!readTransform(mainTransformPath, mainTransform))
    return false;

  // Compute the number of neighboring tiles
  const int numOtherTiles = (argc - 11) / 4;
  //printf("Loading information for %d neighboring tiles.\n", numOtherTiles);
    
  // Load all the information for neighboring tiles
  otherTransforms.resize(numOtherTiles);
  otherWeights.resize(numOtherTiles);
  otherOffsets.resize(numOtherTiles);
  for (int i=0; i<numOtherTiles; ++i)
  {
    int baseIndex = 11+4*i;
    //std::cout << "Reading transform from: " << argv[baseIndex] << std::endl;
    if (!readTransform(argv[baseIndex], otherTransforms[i]))
      return false;
    otherWeights[i]    = atof(argv[baseIndex+1]);
    otherOffsets[i][0] = atoi(argv[baseIndex+2]);
    otherOffsets[i][1] = atoi(argv[baseIndex+3]);
  }

  return true;
}



/// Given the squared distance, compute a smoothly decreasing weight.
double distWeightingFunction(const double distance)
{
  const double MAX_WEIGHT = 1.0;
  const double MIN_WEIGHT = 0.0;
  
  // Using the max tile size makes this reach just to the center of diagonal tiles
  // TODO: This needs to correspond with the HRSC tile size!
  const double MAX_TILE_SIZE = 5792; // Tile size 4096 //  1440;
  const double MAX_DIST = (MAX_TILE_SIZE);//*sqrt(2.0);
    
  // Simple linear drop
  double weight = (MAX_DIST - distance) / MAX_DIST;
  if (weight < MIN_WEIGHT)
    weight = MIN_WEIGHT;
  return weight;
}

/// Generate weighted blend of input pixels
/// - Each adjacent tile (4-connectivity) has a fixed weight.
/// - Absent tiles have zero weight.
/// - There is also a positional weight based on distance.
bool determineBlendedPixel(int height, int width, int row, int col,
                           const cv::Vec3b &mainPixel,  double mainWeightIn, // The main pixel
                           const std::vector<cv::Vec3b> &pixels,     // All the other pixels
                           const std::vector<cv::Vec2i> &offsets,    // Offset in tiles for the other pixels
                           const std::vector<double   > &weightsIn,  //
                           cv::Vec3b &outputPixel)
{
    
  int TEST_COL = 1010;
  int LOW_ROW  = 10;
  int HIGH_ROW = 1014;
  bool DEBUG_PIXEL = false;//((col == TEST_COL) && ((row < LOW_ROW) || (row >= HIGH_ROW)));
      
    
  // Each pixel is influenced by the main tile and each adjacent tile
  const size_t numOtherTiles = pixels.size();
  const size_t numInfluences = pixels.size() + 1;
  std::vector<double> influences(numInfluences);
  
  // The influence is a combination of tile size (weight) and distance to the tile.
  // - The squared distance is used to assess the penalty.  This is important because
  //   the influence of a tile needs to be about zero by the edge of adjacent tiles!
  double tileCenterY = static_cast<double>(height) / 2.0;
  double tileCenterX = static_cast<double>(width)  / 2.0;
  double mainDistY   = (double)row-tileCenterY;
  double mainDistX   = (double)col-tileCenterX;
  double mainDist    = sqrt(mainDistX*mainDistX + mainDistY*mainDistY);
  if (mainDist < 0.1) // Avoid divide by zero at center pixel
  {
    outputPixel = mainPixel;
    return true;
  }
  influences[0] = mainWeightIn * distWeightingFunction(mainDist);
  double influenceSum = influences[0];
  
  if (DEBUG_PIXEL) //DEBUG 
    printf("At row %d, distances: ", row); 
  
  for (size_t i=0; i<numOtherTiles; ++i)
  {
    // Tile distance is computed from the center of the other tile
    double thisTileCenterY = tileCenterY + offsets[i][1]*height;
    double thisTileCenterX = tileCenterX + offsets[i][0]*width;
    double tileDistY   = (double)row-thisTileCenterY;
    double tileDistX   = (double)col-thisTileCenterX;
    double tileDist    = sqrt(tileDistX*tileDistX + tileDistY*tileDistY); // Don't need to check divide by zero here!
    if (DEBUG_PIXEL)
      printf("%lf, ", tileDist);
    influences[i+1] = weightsIn[i] * distWeightingFunction(tileDist);
    influenceSum += influences[i+1];
  }
  
  if (DEBUG_PIXEL) //DEBUG 
    printf("\nAt row %d, influences: ", row); 
  
  // Now normalize all the influences so they total up to 1.0
  for (size_t i=0; i<numInfluences; ++i)
  {
    influences[i] /= influenceSum;
    
    if (DEBUG_PIXEL) //DEBUG 
      printf("%lf, ", influences[i]);
    
  }
  if (DEBUG_PIXEL) //DEBUG 
    printf("\n");
  
  // Compute the final value
  const size_t NUM_RGB_CHANNELS = 3;
  for (size_t c=0; c<NUM_RGB_CHANNELS; ++c)
  {
    // Accumulate the output pixel over all influences
    outputPixel[c] = influences[0]*mainPixel[c];
    for (size_t i=0; i<numOtherTiles; ++i)
    {
      outputPixel[c] += pixels[i][c] * influences[i+1]; // Incorporate input weights
    }    
  }
  
  return true;
}


/// Applies our color transform to a pixel
// TODO: Could just use the OpenCV matrix multiply code here
cv::Vec3b transformPixel(const std::vector<unsigned char> &hrscPixel, const cv::Mat &colorTransform)
{
  // Compute the output pixel values
  cv::Vec3b outputPixel;
  for (int j=0; j<NUM_BASE_CHANNELS; ++j)
  {
    outputPixel[j] = 0;
    float temp = 0.0;
    for (int i=0; i<NUM_HRSC_CHANNELS; ++i)
    {
      temp += static_cast<float>(hrscPixel[i])*colorTransform.at<float>(i,j);
    }
    if (temp < 0.0) // Clamp output to legal range.
      temp = 0;
    if (temp > 255.0)
      temp = 255.0;
    outputPixel[j] = static_cast<unsigned char>(temp);
  }
  return outputPixel;
}


/// Applies our color transform to a pixel
/// - This updated function applies the transform computed in "solveHrscColor.py".
///   See that file for a better description of the transform.
cv::Vec3b transformPixelYCC(const std::vector<unsigned char> &hrscPixel, const cv::Mat &colorTransform)
{
  // First apply the HRSC --> RGB transform
  cv::Vec3b outputPixel;
  for (int j=0; j<NUM_BASE_CHANNELS; ++j)
  {
    outputPixel[j] = 0;
    float temp = 0.0;
    for (int i=0; i<NUM_HRSC_CHANNELS; ++i)
    {
      temp += static_cast<float>(hrscPixel[i])*colorTransform.at<float>(i,j);
    }
    if (temp < 0.0) // Clamp output to legal range.
      temp = 0;
    if (temp > 255.0)
      temp = 255.0;
    outputPixel[j] = static_cast<unsigned char>(temp);
  }
  
  cv::Vec3b ycbcrPixel = rgb2ycbcr(outputPixel);
  
  // Replace the image intensity with the scaled NADIR channel.
  const int Y     = 0;
  const int NADIR = 4;
  float temp2 = static_cast<float>(hrscPixel[NADIR]) * colorTransform.at<float>(NUM_HRSC_CHANNELS, 0);
  if (temp2 < 0.0) // Clamp output to legal range.
    temp2 = 0.0;
  if (temp2 > 255.0)
    temp2 = 255.0;  
  ycbcrPixel[Y] = static_cast<unsigned char>(temp2);

  
  return ycbcr2rgb(ycbcrPixel);
  //return (outputPixel);
}





/// Apply a color transform matrix to the HRSC bands.
/// - This is a 5x3 matrix.
bool transformHrscColor(const std::vector<cv::Mat>   &hrscChannels, const cv::Mat &hrscMask,
                        const BrightnessCorrector &corrector, cv::Mat &outputImage,
                        const cv::Mat                &colorTransform,        double mainWeight,
                        const std::vector<cv::Mat>   &otherColorTransforms,  const std::vector<double> &otherWeights,
                        const std::vector<cv::Vec2i> &otherOffsets)
{
    
  const size_t numOtherTiles = otherOffsets.size();
    
  // Initialize the output image
  const size_t numRows = hrscMask.rows; // The mask is the smallest of the input images
  const size_t numCols = hrscMask.cols;
  outputImage = cv::Mat(numRows, numCols, CV_8UC3);

  //printf("numRows = %d, numCols = %d\n", numRows, numCols);
  
  //cv::Mat weightImage(numRows, numCols, CV_8UC1);
  
  // Iterate over the pixels of the HRSC image
  cv::Vec3b mainPixel, outputPixel;
  std::vector<unsigned char> hrscPixel(NUM_HRSC_CHANNELS);
  for (int r=0; r<numRows; r+=1)
  {
    for (int c=0; c<numCols; c+=1)
    {
      // Handle masked pixels
      if (hrscMask.at<unsigned char>(r, c) == 0)
      {
        outputImage.at<cv::Vec3b>(r, c) = cv::Vec3b(0,0,0);
        continue;
      }
/*
      std::cout << "Input HRSC channels:  ";
      for (int i=0; i<NUM_HRSC_CHANNELS; ++i) {
        std::cout << int(hrscChannels[i].at<unsigned char>(r,c)) << ", ";
      }
      std::cout << "\n";
        */
      // Build the HRSC pixel from the seperate channels with brightness correction
      //std::cout << "Brightness corrected HRSC channels:  ";
      for (int i=0; i<NUM_HRSC_CHANNELS; ++i) {
        hrscPixel[i] = corrector.correctPixel(hrscChannels[i].at<unsigned char>(r,c), r); // Correct brightness
        //std::cout << int(hrscPixel[i]) << ", ";
      }
      //std::cout << "\n";
    
      // Compute the main pixel transform
      mainPixel = transformPixelYCC(hrscPixel, colorTransform);
      //std::cout << "HRSC pixel: " << hrscPixel << std::endl;
      //std::cout << "Main pixel: " << mainPixel << std::endl;
      
      // Compute pixel transforms from all the adjacent tiles
      std::vector<cv::Vec3b> otherPixels(numOtherTiles);
      for (size_t i=0; i<numOtherTiles; ++i)
      {
        otherPixels[i] = transformPixelYCC(hrscPixel, otherColorTransforms[i]);
        //std::cout << "Other pixel " << i << " = " << otherPixels[i] << std::endl;
      }

      // Compute a blended composition of the different color transforms.
      const size_t tileSize = numRows;
      determineBlendedPixel(numRows, numCols, r, c,
                            mainPixel,   mainWeight,
                            otherPixels, otherOffsets, otherWeights,
                            outputPixel);
      
      //std::cout << "Output pixel: " << outputPixel << std::endl;
      
      //double calcWeight = computeMainWeight(numRows, numCols, r, c);
      //double a = calcWeight*255.0;
      //unsigned char b = static_cast<unsigned char>(a);
      //weightImage.at<unsigned char>(r, c) = b;
      //printf("%lf, %lf, %d ", calcWeight, a, b);
      //printf("%lf\n", calcWeight);
      
      // Store the result in output image
      outputImage.at<cv::Vec3b>(r, c) = outputPixel;

    } // End col loop
  } // End row loop

  //printf("Writing weightImage.tif\n");
  //cv::imwrite("weightImage.jpg", weightImage);
  
  return true;
}

//============================================================================


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc < 10)
  {
    printf("usage: transformHrscImageColor <HRSC Red> <HRSC Green> <HRSC Blue> <HRSC NIR> <HRSC Nadir> <HRSC Mask> <Brightness File Path> <Output Path> <Main Color Transform File Path> <Main Weight> [<Color Transform Path> <Weight> <xOffset> <yOffset>]...\n");
    return -1;
  }
  
  printf("Loading input data...\n");
  
  // Load the input images  
  std::vector<cv::Mat> hrscChannels(NUM_HRSC_CHANNELS);
  cv::Mat colorTransform, hrscMask;
  double  mainWeight;
  std::vector<cv::Mat>   otherColorTransforms;
  std::vector<double>    otherWeights;
  std::vector<cv::Vec2i> otherOffsets;
  BrightnessCorrector corrector;
  std::string outputPath;
  if (!loadInputData(argc, argv, hrscChannels, hrscMask, corrector, outputPath,
                     colorTransform, mainWeight,
                     otherColorTransforms, otherWeights, otherOffsets))
    return -1;

  // The color transform is from HRSC to the basemap.

  printf("Transforming image...\n");
  
  // Generate the transformed color image
  cv::Mat newColorImage;
  transformHrscColor(hrscChannels, hrscMask,
                     corrector, newColorImage,
                     colorTransform, mainWeight,
                     otherColorTransforms, otherWeights, otherOffsets);

  // TODO: Apply unsharp masking to the color image?

  
  // Write the output image
  cv::imwrite(outputPath, newColorImage);

  return 0;
}


