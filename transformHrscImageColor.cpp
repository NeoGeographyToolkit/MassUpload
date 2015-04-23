

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <HrscCommon.h>


//=============================================================


bool loadInputData(int argc, char** argv,
                   std::vector<cv::Mat> &hrscChannels,
                   BrightnessCorrector &corrector, std::string &outputPath,
                   cv::Mat &mainTransform, double &mainWeight,
                   std::vector<cv::Mat>   &otherTransforms,
                   std::vector<double>    &otherWeights,
                   std::vector<cv::Vec2i> &otherOffsets)
{
  if (argc < 10)
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
  std::string brightnessPath = argv[6];
  outputPath = argv[7];
  
  std::string mainTransformPath = argv[8];
  mainWeight = atof(argv[9]);
  
  const int LOAD_GRAY = 0;
  const int LOAD_RGB  = 1;
  
  // Load all of the HRSC images
  for (size_t i=0; i<NUM_HRSC_CHANNELS; ++i)
  {
    //std::cout << "Reading image from: " << hrscPaths[i] << std::endl;
    hrscChannels[i] = cv::imread(hrscPaths[i], LOAD_GRAY);
    if (!hrscChannels[i].data)
    {
      printf("Failed to load HRSC image: %s\n", hrscPaths[i].c_str());
      return false;
    }
  }

  // Load brightness correction data
  if (!corrector.readProfileCorrection(brightnessPath))
    return false;
  
  // Load the main spatial transform
  //std::cout << "Reading transform from: " << mainTransformPath << std::endl;
  if (!readTransform(mainTransformPath, mainTransform))
    return false;

  // Compute the number of neighboring tiles
  const int numOtherTiles = (argc - 10) / 4;
  //printf("Loading information for %d neighboring tiles.\n", numOtherTiles);
    
  // Load all the information for neighboring tiles
  otherTransforms.resize(numOtherTiles);
  otherWeights.resize(numOtherTiles);
  otherOffsets.resize(numOtherTiles);
  for (int i=0; i<numOtherTiles; ++i)
  {
    int baseIndex = 10+4*i;
    //std::cout << "Reading transform from: " << argv[baseIndex] << std::endl;
    if (!readTransform(argv[baseIndex], otherTransforms[i]))
      return false;
    otherWeights[i]    = atof(argv[baseIndex+1]);
    otherOffsets[i][0] = atoi(argv[baseIndex+2]);
    otherOffsets[i][1] = atoi(argv[baseIndex+3]);
  }

  return true;
}



/*
/// Generate weighted blend of input pixels
/// - Each adjacent tile (4-connectivity) has a fixed weight.
/// - Absent tiles have zero weight.
/// - There is also a positional weight based on distance.
bool determineBlendedPixel(int tileSize, int row, int col,
                           const cv::Vec3b &mainPixel,   double mainWeightIn,
                           const cv::Vec3b &leftPixel,   double leftWeightIn,
                           const cv::Vec3b &topPixel,    double topWeightIn,
                           const cv::Vec3b &rightPixel,  double rightWeightIn,
                           const cv::Vec3b &bottomPixel, double bottomWeightIn,
                           cv::Vec3b &outputPixel)
{
  // Compute a distance from the center
  double tileCenter = (double)tileSize / 2.0;
  double mainDistY  = abs((double)row-tileCenter);
  double mainDistX  = abs((double)col-tileCenter);
  double mainDist   = min(tileDistY, tileDistX);
  
  // At the edges of a tile, that tile's pixel has 50% weight.
  double otherTileWeight = (mainDist/tileCenter)/2.0;
  double mainWeight      = 1.0 - otherTileWeight
  
  // Compute the relative contribution of the other tiles
  // - There are never more than two other contributors
  double leftWeight   = tileCenter - (double)col; // Get other weights in 0-1 range
  double rightWeight  = -leftWeight;
  double topWeight    = tileCenter - (double)row;
  double bottomWeight = -topWeight;
  if (leftWeight   < 0) leftWeight   = 0; // Clamp to minimum weight of zero.
  if (rightWeight  < 0) rightWeight  = 0;
  if (topWeight    < 0) topWeight    = 0;
  if (bottomWeight < 0) bottomWeight = 0;
  
  // Normalize the weights
  double totalOtherWeight = leftWeight + rightWeight + topWeight + bottomWeight;
  if (totalOtherWeight > 0) // Handle case where the pixel is right on the tile center
  {
    leftWeight   /= totalOtherWeight;
    rightWeight  /= totalOtherWeight;
    topWeight    /= totalOtherWeight;
    bottomWeight /= totalOtherWeight;
  }
  
  // Compute the final value
  const size_t NUM_RGB_CHANNELS = 3;
  for (size_t i=0; i<NUM_RGB_CHANNELS; ++i)
  {
    outputPixel[i] = mainWeight      * mainPixel[i] +
                     otherTileWeight * (leftPixel  [i]*leftWeight   +
                                        rightPixel [i]*rightWeight  +
                                        topPixel   [i]*topWeight    +
                                        bottomPixel[i]*bottomWeight  );
  }
  
  return true;
}
*/


/// Generate weighted blend of input pixels
/// - Each adjacent tile (4-connectivity) has a fixed weight.
/// - Absent tiles have zero weight.
/// - There is also a positional weight based on distance.
bool determineBlendedPixel(int height, int width, int row, int col,
                           const cv::Vec3b &mainPixel,  double mainWeightIn, // The main pixel
                           const std::vector<cv::Vec3b> &pixels,     // All the other pixels
                           const std::vector<cv::Vec2i> &offsets,
                           const std::vector<double   > &weightsIn,
                           cv::Vec3b &outputPixel)
{
  // Compute a distance from the center
  double tileCenterY = static_cast<double>(width) / 2.0;
  double tileCenterX = static_cast<double>(height) / 2.0;
  double mainDistY   = abs((double)row-tileCenterY);
  double mainDistX   = abs((double)col-tileCenterX);
  double otherTileWeight = ((mainDistX/tileCenterX)+(mainDistY/tileCenterY))/4.0; // This could be improved
  double mainWeight = 1.0 - otherTileWeight;
  //std::cout << "Main distance weight = " << mainWeight << std::endl;
  // At the edges of a tile, that tile's pixel has 50% weight.
  
  //printf("row = %d, col = %d\n", row, col);
  //printf("tileCenterX = %lf, tileCenterY = %lf\n", tileCenterX, tileCenterY);
  
  
  
  // Compute the relative contribution of the other tiles
  // - There are never more than two other contributors
  const size_t numOtherTiles = pixels.size();
  double totalOtherPositionWeight = 0.0;
  std::vector<double> otherPositionWeights(numOtherTiles);
  for (size_t i=0; i<numOtherTiles; ++i) // Loop through the other tiles
  {
    // Compute a 0 to 1 fraction of the position relative to this tile
    double xDist    = ((double)col - tileCenterX) * offsets[i][0];
    double yDist    = ((double)row - tileCenterY) * offsets[i][1];
    
    //printf("Position %d, xDist = %lf, yDist = %lf\n", i, xDist, yDist);
    
    // TODO: Probably need to tweak this to handle diagonals better
    // Negative distances = zero weight!
    // - Tile weight does not carry past the center of the main tile.
    double thisWeight;
    if (xDist > yDist) // Use the closest distance to compute weight
      thisWeight = xDist / tileCenterX;
    else
      thisWeight = yDist / tileCenterY;
    //printf("This weight = %lf\n", thisWeight);
    if (thisWeight <= 0) // If on the opposite side of the tile center, no weight
      thisWeight = 0;
    thisWeight *= weightsIn[i]; // Incorporate the input weight for this tile
    otherPositionWeights[i]   = thisWeight;
    totalOtherPositionWeight += thisWeight;
  }
  // If none of the other tiles have weight, just return the main pixel.
  if (totalOtherPositionWeight == 0)
  {
    outputPixel = mainPixel;
    return true;
  }
  
  // Normalize the position weighting of the other tiles
  for (size_t i=0; i<numOtherTiles; ++i) // Loop through the other tiles
  {
    otherPositionWeights[i] /= totalOtherPositionWeight;
    //std::cout << "Normalized position weight " << i << " = " << otherPositionWeights[i] << std::endl;
  }
  
  
  // Compute the final value
  const size_t NUM_RGB_CHANNELS = 3;
  for (size_t c=0; c<NUM_RGB_CHANNELS; ++c)
  {
    double otherWeight = 0; // Accumulate value of the other tiles
    for (size_t i=0; i<numOtherTiles; ++i)
    {
      //std::cout << "Adding other weight: " << otherPositionWeights[i] << std::endl;
      otherWeight += pixels[i][c] * otherPositionWeights[i]; // Incorporate input weights
    }
    // Compute the final pixel value for this channel
    //std::cout << "Weighted other positions: " <<  otherWeight*otherTileWeight << std::endl;
    //std::cout << "Weighted main: " << mainWeight*mainPixel[c] << std::endl;
    outputPixel[c] = mainWeight*mainPixel[c] + otherWeight*otherTileWeight;
    
  }
  
  //std::cout << "outputPixel: " << outputPixel << std::endl;
  //int die = 5/0;
  
  return true;
}


cv::Vec3b transformPixel(const std::vector<unsigned char> &hrscPixel, const cv::Mat &colorTransform)
{
  //TODO: Speed this up using a built in OpenCV matrix transform call
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


/// Apply a color transform matrix to the HRSC bands.
/// - This is a 5x3 matrix.
bool transformHrscColor(const std::vector<cv::Mat>   &hrscChannels,
                        const BrightnessCorrector &corrector, cv::Mat &outputImage,
                        const cv::Mat                &colorTransform,        double mainWeight,
                        const std::vector<cv::Mat>   &otherColorTransforms,  const std::vector<double> &otherWeights,
                        const std::vector<cv::Vec2i> &otherOffsets)
{
    
  const size_t numOtherTiles = otherOffsets.size();
    
  // Initialize the output image
  const size_t numRows = hrscChannels[0].rows;
  const size_t numCols = hrscChannels[0].cols;
  outputImage = cv::Mat(numRows, numCols, CV_8UC3);

  //printf("numRows = %d, numCols = %d\n", numRows, numCols);
  
  // Iterate over the pixels of the HRSC image
  cv::Vec3b mainPixel, outputPixel;
  std::vector<unsigned char> hrscPixel(NUM_HRSC_CHANNELS);
  for (int r=0; r<numRows; r+=1)
  {
    for (int c=0; c<numCols; c+=1)
    {  
      // Check to see if this pixel should be masked out   
      // - Apply the brightness correction while we are at it.
      bool maskOut = false;
      for (int i=0; i<NUM_HRSC_CHANNELS; ++i)
      {
        hrscPixel[i] = hrscChannels[i].at<unsigned char>(r,c);
        if (hrscPixel[i] == 0)
        {
          maskOut = true;
          break;
        }
        hrscPixel[i] = corrector.correctPixel(hrscChannels[i].at<unsigned char>(r,c), r); // Correct brightness
      }
      // If any of input HRSC pixels are black, the output pixel is black.
      // TODO: A smarter mask method!
      if (maskOut)
      {
        outputImage.at<cv::Vec3b>(r, c) = cv::Vec3b(0,0,0);
        continue;
      }
    
      // Compute the main pixel transform
      mainPixel = transformPixel(hrscPixel, colorTransform);
      //std::cout << "Main pixel: " << mainPixel << std::endl;
      
      // Compute pixel transforms from all the adjacent tiles
      std::vector<cv::Vec3b> otherPixels(numOtherTiles);
      for (size_t i=0; i<numOtherTiles; ++i)
      {
        otherPixels[i] = transformPixel(hrscPixel, otherColorTransforms[i]);
        //std::cout << "Other pixel " << i << " = " << otherPixels[i] << std::endl;
      }

      // Compute a blended composition of the different color transforms.
      const size_t tileSize = numRows;
      determineBlendedPixel(numRows, numCols, r, c,
                            mainPixel,  mainWeight,
                            otherPixels, otherOffsets,otherWeights,
                            outputPixel);
      
      
      // Store the result in output image
      outputImage.at<cv::Vec3b>(r, c) = outputPixel;

    } // End col loop
  } // End row loop

  return true;
}

//============================================================================


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc < 10)
  {
    printf("usage: transformHrscImageColor <HRSC Red> <HRSC Green> <HRSC Blue> <HRSC NIR> <HRSC Nadir> <Brightness File Path> <Output Path> <Main Color Transform File Path> <Main Weight> [<Color Transform Path> <Weight> <xOffset> <yOffset>]...\n");
    return -1;
  }
  
  printf("Loading input data...\n");
  
  // Load the input images  
  std::vector<cv::Mat> hrscChannels(NUM_HRSC_CHANNELS);
  cv::Mat colorTransform;
  double  mainWeight;
  std::vector<cv::Mat>   otherColorTransforms;
  std::vector<double>    otherWeights;
  std::vector<cv::Vec2i> otherOffsets;
  BrightnessCorrector corrector;
  std::string outputPath;
  if (!loadInputData(argc, argv, hrscChannels, corrector, outputPath,
                     colorTransform, mainWeight,
                     otherColorTransforms, otherWeights, otherOffsets))
    return -1;

  // The color transform is from HRSC to the basemap.

  printf("Transforming image...\n");
  
  // Generate the transformed color image
  cv::Mat newColorImage;
  transformHrscColor(hrscChannels,
                     corrector, newColorImage,
                     colorTransform, mainWeight,
                     otherColorTransforms, otherWeights, otherOffsets);

  // TODO: Apply unsharp masking to the color image?

  
  // Write the output image
  cv::imwrite(outputPath, newColorImage);

  return 0;
}


