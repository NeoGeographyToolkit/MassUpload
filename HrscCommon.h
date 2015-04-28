

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <vw/Math/Geometry.h>
#include <vw/Math/RANSAC.h>
#include <vw/InterestPoint/InterestData.h>



const size_t NUM_HRSC_CHANNELS = 5;
const size_t NUM_BASE_CHANNELS = 3;

void affineTransform(const cv::Mat &transform, float xIn, float yIn, float &xOut, float &yOut)
{
  xOut = xIn*transform.at<float>(0,0) + yIn*transform.at<float>(0,1) + transform.at<float>(0,2);
  yOut = xIn*transform.at<float>(1,0) + yIn*transform.at<float>(1,1) + transform.at<float>(1,2);
}

/// Single channel image interpolation
template <typename T>
T interpPixel(const cv::Mat& img, const cv::Mat& mask, float xF, float yF, bool &gotValue)
{
  const int BORDER_SIZE = 1; // Stay away from border artifacts

  gotValue = false;
  int x = (int)xF;
  int y = (int)yF;

  // Get the bordering pixel coordinates, replacing out of bounds with zero.
  int minX = BORDER_SIZE; // Max legal pixel boundaries with the specified border.
  int minY = BORDER_SIZE;
  int maxX = img.cols-BORDER_SIZE;
  int maxY = img.rows-BORDER_SIZE;
  int x0 = x;   // The coordinates of the four bordering pixels
  int x1 = x+1;
  int y0 = y;
  int y1 = y+1;
  if ((x0 < minX) || (x0 >= maxX)) return 0; // Quit if we exceed any of the borders.
  if ((x1 < minX) || (x1 >= maxX)) return 0;
  if ((y0 < minY) || (y0 >= maxY)) return 0;
  if ((y1 < minY) || (y1 >= maxY)) return 0;
  
  // - Don't interpolate if any mask inputs are zero, this might indicate 
  //    that we are at a projection border.
  unsigned char i00 = mask.at<unsigned char>(y0, x0);
  unsigned char i01 = mask.at<unsigned char>(y0, x1);
  unsigned char i10 = mask.at<unsigned char>(y1, x0);
  unsigned char i11 = mask.at<unsigned char>(y1, x1);
  if ((i00 == 0) || (i01 == 0) || (i10 == 0) || (i11 == 0))
    return 0;


  float a = xF - (float)x;
  float c = yF - (float)y;
  
  float v00 = static_cast<float>(img.at<T>(y0, x0));
  float v01 = static_cast<float>(img.at<T>(y0, x1));
  float v10 = static_cast<float>(img.at<T>(y1, x0));
  float v11 = static_cast<float>(img.at<T>(y1, x1));

  T val = static_cast<T>( v00*(1-a)*(1-c)  + v10*a*(1-c) + v01*(1-a)*c + v11*a*c );

  gotValue = true;
  return val;
}

/// As interpPixel but specialized for RGB
cv::Vec3b interpPixelRgb(const cv::Mat& img, float xF, float yF, bool &gotValue)
{
  const size_t NUM_RGB_CHANNELS = 3;
  const int    BORDER_SIZE      = 1; // Stay away from border artifacts

  gotValue = false;
  int x = (int)xF;
  int y = (int)yF;

  // Get the bordering pixel coordinates, replacing out of bounds with zero.
  int minX = BORDER_SIZE;
  int minY = BORDER_SIZE;
  int maxX = img.cols-BORDER_SIZE;
  int maxY = img.rows-BORDER_SIZE;
  int x0 = x;
  int x1 = x+1;
  int y0 = y;
  int y1 = y+1;
  if ((x0 < minX) || (x0 >= maxX)) return 0;
  if ((x1 < minX) || (x1 >= maxX)) return 0;
  if ((y0 < minY) || (y0 >= maxY)) return 0;
  if ((y1 < minY) || (y1 >= maxY)) return 0;

  // Now interpolate each pixel channel

  float a = xF - (float)x;
  float c = yF - (float)y;
  
  cv::Vec3b outputPixel;
  for (size_t i=0; i<NUM_RGB_CHANNELS; ++i)
  {
    float v00 = static_cast<float>(img.at<cv::Vec3b>(y0, x0)[i]);
    float v01 = static_cast<float>(img.at<cv::Vec3b>(y0, x1)[i]);
    float v10 = static_cast<float>(img.at<cv::Vec3b>(y1, x0)[i]);
    float v11 = static_cast<float>(img.at<cv::Vec3b>(y1, x1)[i]);

    outputPixel[i] = static_cast<unsigned char>( v00*(1-a)*(1-c)  + v10*a*(1-c) + v01*(1-a)*c + v11*a*c );

  }

  gotValue = true;
  return outputPixel;
}

/// As interpPixelRgb but with pixels near the edges handled by mirroring
cv::Vec3b interpPixelMirrorRgb(const cv::Mat& img,  const cv::Mat& mask,
                               float xF, float yF, bool &gotValue)
{
  const size_t NUM_RGB_CHANNELS = 3;

  // Get the bounding pixel coordinates
  gotValue = false;
  int x = (int)xF;
  int y = (int)yF;
  int x0 = x;
  int x1 = x+1;
  int y0 = y;
  int y1 = y+1;
  /*
  // Mirror a border of one by adjusting the bounding coordinates
  if (x0 == -1)       x0 = 0;
  if (y0 == -1)       y0 = 0;
  if (x1 == img.cols) x1 = img.cols-1;
  if (y1 == img.rows) y1 = img.rows-1;
  */
  // Pixels past the border are still rejected
  if ((x0 < 0) || (x0 >= img.cols)) return 0;
  if ((x1 < 0) || (x1 >= img.cols)) return 0;
  if ((y0 < 0) || (y0 >= img.rows)) return 0;
  if ((y1 < 0) || (y1 >= img.rows)) return 0;

  // Check the mask
  // - Don't interpolate if any mask inputs are zero, this might indicate 
  //    that we are at a projection border.
  unsigned char i00 = mask.at<unsigned char>(y0, x0);
  unsigned char i01 = mask.at<unsigned char>(y0, x1);
  unsigned char i10 = mask.at<unsigned char>(y1, x0);
  unsigned char i11 = mask.at<unsigned char>(y1, x1);
  if ((i00 == 0) || (i01 == 0) || (i10 == 0) || (i11 == 0))
    return 0;
  
  // Now interpolate each pixel channel

  float a = xF - (float)x;
  float c = yF - (float)y;
  
  cv::Vec3b outputPixel;
  for (size_t i=0; i<NUM_RGB_CHANNELS; ++i)
  {
    float v00 = static_cast<float>(img.at<cv::Vec3b>(y0, x0)[i]);
    float v01 = static_cast<float>(img.at<cv::Vec3b>(y0, x1)[i]);
    float v10 = static_cast<float>(img.at<cv::Vec3b>(y1, x0)[i]);
    float v11 = static_cast<float>(img.at<cv::Vec3b>(y1, x1)[i]);  
    outputPixel[i] = static_cast<unsigned char>( v00*(1.0f-a)*(1.0f-c)  + v10*a*(1.0f-c) + v01*(1.0f-a)*c + v11*a*c );
  }
  
  gotValue = true;
  return outputPixel;
}


/// Computes the ROI of one image in another given the transform with bounds checking.
cv::Rect_<int> getboundsInOtherImage(const cv::Mat &imageA, const cv::Mat &imageB, const cv::Mat &transB_to_A)
{
  // Transform the four corners of imageB
  float x[4], y[4];
  affineTransform(transB_to_A, 0,             0,             x[0], y[0]);
  affineTransform(transB_to_A, imageB.cols-1, 0,             x[1], y[1]);
  affineTransform(transB_to_A, imageB.cols-1, imageB.rows-1, x[2], y[2]);
  affineTransform(transB_to_A, 0,             imageB.rows-1, x[3], y[3]);
  
  // Get the bounding box of the transformed points
  float xMin = x[0];
  float xMax = x[0];
  float yMin = y[0];
  float yMax = y[0];
  for (size_t i=0; i<4; ++i)
  {
    if (x[i] < xMin) xMin = x[i];
    if (x[i] > xMax) xMax = x[i];
    if (y[i] < yMin) yMin = y[i];
    if (y[i] > yMax) yMax = y[i];
  }
  
  if (xMin < 0) xMin = 0;
  if (yMin < 0) yMin = 0;
  if (xMax > imageA.cols-1) xMax = imageA.cols-1;
  if (yMax > imageA.rows-1) yMax = imageA.rows-1;

  // Return the results expanded to the nearest integer
  cv::Rect_<int> boundsInA(static_cast<int>(floor(xMin)), 
                           static_cast<int>(floor(yMin)),
                           static_cast<int>(ceil(xMax-xMin)), 
                           static_cast<int>(ceil(yMax-yMin)));
  return boundsInA;
}

/// Write a small matrix to a text file
bool writeTransform(const std::string &outputPath, const cv::Mat &transform)
{
  std::ofstream file(outputPath.c_str());
  file << transform.rows << ", " << transform.cols << std::endl;
  for (size_t r=0; r<transform.rows; ++r)
  {
    for (size_t c=0; c<transform.cols-1; ++c)
    {
      file << transform.at<float>(r,c) << ", ";
    }
    file << transform.at<float>(r,transform.cols-1) << std::endl;
  }
  file.close();
  
  return (!file.fail());
}

// Read a small matrix from a text file
bool readTransform(const std::string &inputPath, cv::Mat &transform)
{
  std::ifstream file(inputPath.c_str());
  char   comma;
  size_t numRows, numCols;
  file >> numRows >> comma >> numCols;
  transform.create(numRows, numCols, CV_32FC1);
  for (size_t r=0; r<transform.rows; ++r)
  {
    for (size_t c=0; c<transform.cols-1; ++c)
    {
      file >> transform.at<float>(r,c) >> comma;
    }
    file >> transform.at<float>(r,transform.cols-1);
  }
  file.close();
  
  return (!file.fail());
}

/// Try to load the image and then make sure we got valid data.
/// - The type must by 0 (gray) or 1 (RGB)
bool readOpenCvImage(const std::string &imagePath, cv::Mat &image, const int imageType)
{
  //printf("Reading image file: %s\n", imagePath.c_str());
  image = cv::imread(imagePath, imageType);
  if (!image.data)
  {
    printf("Failed to load image %s!\n", imagePath.c_str());
    return false;
  }
  return true;
}


/// Helper class for working with brightness information
class BrightnessCorrector
{
public:

  /// Data loading options
  BrightnessCorrector() {}
  BrightnessCorrector(cv::Mat &gain, cv::Mat &offset) : _gain(gain), _offset(offset) { }
  void set(cv::Mat &gain, cv::Mat &offset) { _gain = gain; _offset = offset; }

  /// Write a gain/offset pair to a CSV file
  bool writeProfileCorrection(const std::string &outputPath) const
  {
    std::ofstream file(outputPath.c_str());
    file << _gain.rows << std::endl;
    for (size_t r=0; r<_gain.rows; ++r)
    {
      file << _gain.at<float>(r,0) << ", " << _offset.at<float>(r,0) << std::endl;
    }
    file.close();
    return (!file.fail());
  }

  /// Write a gain/offset pair to a CSV file
  bool readProfileCorrection(const std::string &inputPath)
  {
    std::ifstream file(inputPath.c_str());
    char   comma;
    size_t numRows;
    file >> numRows;
    _gain.create(numRows, 1, CV_32FC1);
    _offset.create(numRows, 1, CV_32FC1);
    for (size_t r=0; r<numRows; ++r)
    {
      file >> _gain.at<float>(r,0) >> comma >> _offset.at<float>(r,0);
    }
    file.close();   
    return (!file.fail());
  }

  /// Get the corrected value of a single pixel
  unsigned char correctPixel(unsigned char inputPixel, int row) const
  {
    float result = static_cast<float>(inputPixel) * _gain.at<float>(row,0);
    if (result <   0.0) result = 0.0;  // Clamp the output value
    if (result > 255.0) result = 255.0;
    return static_cast<unsigned char>(result);
  }
  
private:

  cv::Mat _gain;
  cv::Mat _offset;
};





/// Replace the Value channel of the input HSV image
bool replaceValue(const cv::Mat &baseImageRgb, const cv::Mat &spatialTransform, const cv::Mat &nadir, cv::Mat &outputImage)
{
  printf("Converting image...\n");
  
  // Convert the input image to HSV
  cv::Mat hsvImage;
  cv::cvtColor(baseImageRgb, hsvImage, cv::COLOR_BGR2HSV);
 
  printf("Replacing value channel...\n");
 
  // TODO: There must be a better way to do this using OpenCV!
  // Replace the value channel
  //cv::Mat outputMask;
  bool gotValue;
  for (int r=0; r<baseImageRgb.rows; ++r)
  {
    for (int c=0; c<baseImageRgb.cols; ++c)
    {     
      float matchX = c*spatialTransform.at<float>(0,0) + r*spatialTransform.at<float>(0,1) + spatialTransform.at<float>(0,2);
      float matchY = c*spatialTransform.at<float>(1,0) + r*spatialTransform.at<float>(1,1) + spatialTransform.at<float>(1,2);
      
      unsigned char newVal = interpPixel<unsigned char>(nadir, nadir, matchX, matchY, gotValue);
      //hsvImage.at<unsigned char>(r,c, 2) = newVal;
      if (gotValue)
        hsvImage.at<cv::Vec3b>(r,c)[2] = newVal;
    }
  }
  cv::cvtColor(hsvImage, outputImage, cv::COLOR_HSV2BGR);
  
  cv::imwrite("value_replaced_image.jpeg", outputImage);
  
  //cv::namedWindow("Display Image", cv::WINDOW_AUTOSIZE );
  //cv::imshow("Display Image", outputImage);
  //cv::waitKey(0);
  
  printf("Finished replacing value\n");
  return true;

}






