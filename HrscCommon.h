

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <vw/Math/Geometry.h>
#include <vw/Math/RANSAC.h>
#include <vw/InterestPoint/InterestData.h>



const size_t NUM_HRSC_CHANNELS = 5;
const size_t NUM_BASE_CHANNELS = 3;


template <typename T>
T interpPixel(const cv::Mat& img, const cv::Mat& mask, float xF, float yF, bool &gotValue)
{
    const int BORDER_SIZE = 1; // Stay away from border artifacts

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

    T val = static_cast<short>( v00*(1-a)*(1-c)  + v10*a*(1-c) + v01*(1-a)*c + v11*a*c );

    gotValue = true;
    return val;
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
  transform = cv::Mat(numRows, numCols, CV_32FC1);
  for (size_t r=0; r<transform.rows; ++r)
  {
    for (size_t c=0; c<transform.cols-1; ++c)
    {
      file >> transform.at<float>(r,c) >> comma;
    }
    file >> transform.at<float>(r,transform.cols-1)
  }
  file.close();
  
  return (!file.fail());
}








