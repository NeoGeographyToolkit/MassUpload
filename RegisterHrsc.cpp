

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <vw/Math/Geometry.h>
#include <vw/Math/RANSAC.h>
#include <vw/InterestPoint/InterestData.h>

/**
  Program to compute the transform between an HRSC image and a base image using OpenCV
*/

/*
cv::Point templateMatchPoint(const cv::Mat &findImage, const cv::Mat &templateImage, 
                             const int      centerX,   const int      centerY)
{

//SQDIFF
//SQDIFF NORMED
//TM CCORR
//TM CCORR NORMED
//TM COEFF
//TM COEFF NORMED

  const int TEMPLATE_SIZE = 7;
  const int radius = (TEMPLATE_SIZE - 1) / 2;

  // Create the template image
  cv::Rect roi(centerX-radius, centerY-radius, TEMPLATE_SIZE, TEMPLATE_SIZE);
  cv::Mat templ(templateImage, roi);
  
  // Create the result matrix
  cv::Mat result;
  int result_cols = findImage.cols - templ.cols + 1;
  int result_rows = findImage.rows - templ.rows + 1;
  result.create( result_cols, result_rows, CV_32FC1 );  
  
  
  matchTemplate(findImage, templ, result, CV_TM_SQDIFF_NORMED );
  normalize( result, result, 0, 1, NORM_MINMAX, -1, Mat() );


  // Localizing the best match with minMaxLoc
  double    minVal, maxVal; 
  cv::Point minLoc, maxLoc, matchLoc;

  minMaxLoc( result, &minVal, &maxVal, &minLoc, &maxLoc, Mat() );

  // For SQDIFF and SQDIFF_NORMED, the best matches are lower values. 
  //   For all the other methods, the higher the better
  if( match_method  == CV_TM_SQDIFF || match_method == CV_TM_SQDIFF_NORMED )
    matchLoc = minLoc;
  else
    matchLoc = maxLoc;

  return matchLoc;
}
*/
/*
bool matchByTemplate(const cv::Mat &imageA, const cv::Mat &imageB, 
                     const std::vector<cv::KeyPoint> &keypointsA, 
                           std::vector<cv::KeyPoint> &keypointsB)
{

  // TODO: Handle edge cases?
  keypointsB.resize(keypointsA.size());
  for (size_t i=0; i<keypointsA.size(); ++i)
  {
    
    // Search for a match for keypoint A
    cv::Point foundPtB = templateMatchPoint(imageB, imageA, 
                                            keypointsA[i].pt.x, keypointsA[i].pt.y);


}
*/

/*
struct L2NormErrorMetric {
  double operator() (vw::math::TranslationFittingFunctor::result_type  const& H,
                     vw::Vector3 const& p1,
                     vw::Vector3 const& p2) const {
    std::cout<< "h = " << H << std::endl;
    return vw::math::norm_2( p2 - H * p1 );
  }
};
*/

/// Compute an affine transform for the given feature points using
///  the Vision Workbench RANSAC implementation
bool vwRansacAffine(const std::vector<cv::Point2f> &keypointsA, 
                    const std::vector<cv::Point2f> &keypointsB,
                    cv::Mat &transformMatrix)
{
  // Convert points to VW format
  const size_t numPoints = keypointsA.size();
  std::vector<vw::Vector3> vwPtsA(numPoints), vwPtsB(numPoints);
  for (size_t i=0; i<numPoints; ++i)
  {
    vwPtsA[i][0] = keypointsA[i].x;
    vwPtsA[i][1] = keypointsA[i].y;
    vwPtsA[i][2] = 1;
    vwPtsB[i][0] = keypointsB[i].x;
    vwPtsB[i][1] = keypointsB[i].y;
    vwPtsB[i][2] = 1;
    //printf("Pair: (%lf, %lf) ==> (%lf, %lf) \n", vwPtsA[i][0], vwPtsA[i][1], vwPtsB[i][0], vwPtsB[i][1]);
  }

  // Call VW RANSAC function
  typedef vw::math::TranslationFittingFunctor FittingFunctorType;
  typedef vw::math::InterestPointErrorMetric  ErrorFunctorType;
  int    num_iterations         = 100;
  double inlier_threshold       = 30; // Max point distance in pixels
  int    min_num_output_inliers = 10; // Min pixels to count as a match.
  bool   reduce_min_num_output_inliers_if_no_fit = true;
  vw::math::RandomSampleConsensus<FittingFunctorType, ErrorFunctorType> 
      ransac_instance(FittingFunctorType(), ErrorFunctorType(), 
                      num_iterations, inlier_threshold,
                      min_num_output_inliers, reduce_min_num_output_inliers_if_no_fit);
  FittingFunctorType::result_type vwTransform = ransac_instance(vwPtsA, vwPtsB);
  
  // Convert output to OpenCV format
  transformMatrix = cv::Mat(3, 3, CV_32FC1);
  for (int r=0; r<3; ++r)
    for (int c=0; c<3; ++c)
      transformMatrix.at<float>(r, c) = vwTransform[r][c];
  
  return true;
}


bool computeImageTransform(const cv::Mat &refImageIn, const cv::Mat &matchImageIn,
                                 cv::Mat &transform)
{

  
  // Preprocess the images to improve feature detection
  cv::Mat refImage, matchImage, temp;
  //GaussianBlur( src, src, Size(3,3), 0, 0, BORDER_DEFAULT );
  //cv::GaussianBlur( matchImageIn, matchImageIn, cv::Size(3,3), 0, 0, cv::BORDER_DEFAULT );
  
  
  const int kernel_size = 5;
  const int scale = 1;
  const int delta = 0;
  cv::Laplacian( refImageIn, temp, CV_16S, kernel_size, scale, delta, cv::BORDER_DEFAULT );
  cv::convertScaleAbs( temp, refImage, 0.3);
  cv::Laplacian( matchImageIn, temp, CV_16S, kernel_size, scale, delta, cv::BORDER_DEFAULT );
  cv::convertScaleAbs( temp, matchImage, 0.3 );

  
  std::vector<cv::KeyPoint> keypointsA, keypointsB;
  cv::Mat descriptorsA, descriptorsB;  

  cv::Ptr<cv::FeatureDetector    > detector  = cv::BRISK::create();
  cv::Ptr<cv::DescriptorExtractor> extractor = cv::BRISK::create();

  detector->detect(  refImage, keypointsA);
  extractor->compute(refImage, keypointsA, descriptorsA);

  detector->detect(  matchImage, keypointsB);
  extractor->compute(matchImage, keypointsB, descriptorsB);

  // Rule out obviously bad matches based on the known starting alignment accuracy
  cv::Mat mask(keypointsA.size(), keypointsB.size(), CV_8UC1);
  const float MAX_MATCH_PIXEL_DISTANCE = 3000; // Tune this for the base map resolution!
  cv::Point2f centerA(refImage.cols/2,   refImage.rows/2  );
  cv::Point2f centerB(matchImage.cols/2, matchImage.rows/2);
  for (size_t i=0; i<keypointsA.size(); ++i)
  {
    for (size_t j=0; j<keypointsB.size(); ++j)
    {
      cv::Point2f cA = keypointsA[i].pt - centerA;
      cv::Point2f cB = keypointsB[i].pt - centerB;
      float distance = cv::norm(cA - cB);
      if (distance < MAX_MATCH_PIXEL_DISTANCE)
        mask.at<uchar>(i,j) = 1;
      else // Too far, disallow a match
        mask.at<uchar>(i,j) = 0;
    } 
  }

  // Find the closest match for each feature
  //cv::FlannBasedMatcher matcher;
  cv::Ptr<cv::DescriptorMatcher> matcher = cv::DescriptorMatcher::create("BruteForce-Hamming");
  std::vector<cv::DMatch> matches;
  matcher->match(descriptorsA, descriptorsB, matches, mask);


  //-- Quick calculation of max and min distances between keypoints
  double max_dist = 0; double min_dist = 6000;
  for (int i=0; i<descriptorsA.rows; i++)
  { 
    if (matches[i].queryIdx < 0)
      continue;
    double dist = matches[i].distance;
    if (dist < min_dist) 
      min_dist = dist;
    if (dist > max_dist) 
      max_dist = dist;
  }

  printf("-- Max dist : %f \n", max_dist );
  printf("-- Min dist : %f \n", min_dist );

  //-- Pick out "good" matches
  float goodDist = 120;
  //if (argc > 3)
  //  goodDist = atof(argv[3]);
  std::vector< cv::DMatch > good_matches;
  for (int i=0; i<descriptorsA.rows; i++)
  { 
    // First verify that the match is valid
    if ( (matches[i].queryIdx < 0) ||
         (matches[i].trainIdx < 0) ||  
         (matches[i].queryIdx >= keypointsA.size()) || 
         (matches[i].trainIdx >= keypointsB.size()) )
      continue;
    // Now check the distance
    if (matches[i].distance < goodDist)
    {
      good_matches.push_back( matches[i]);
    }
  }
  //good_matches = matches;

  // Compute a transform between the images using RANSAC
  std::vector<cv::Point2f> refPts;
  std::vector<cv::Point2f> matchPts;
  for( int i = 0; i < good_matches.size(); i++ )
  {
    // Get the keypoints from the good matches
    refPts.push_back  (keypointsA[good_matches[i].queryIdx].pt);
    matchPts.push_back(keypointsB[good_matches[i].trainIdx].pt);
    //printf("Pair: HRSC(%lf, %lf) ==> NOEL(%lf, %lf),  DIFF = (%lf, %lf)\n", 
    //        matchPts[i].x, matchPts[i].y, refPts[i].x, refPts[i].y,
    //        refPts[i].x-matchPts[i].x, refPts[i].y-matchPts[i].y);
  }
  // Compute a transform using the Vision Workbench RANSAC tool
  cv::Mat H;
  vwRansacAffine(matchPts, refPts, H);
  std::cout << "H = \n" << H << std::endl;
   
  cv::perspectiveTransform(matchPts, refPts, H);
  for( int i = 0; i < good_matches.size(); i++ )
  {
    //printf("Pair: HRSC(%lf, %lf) ==> NOEL(%lf, %lf),  DIFF = (%lf, %lf)\n", 
    //        matchPts[i].x, matchPts[i].y, refPts[i].x, refPts[i].y,
    //        refPts[i].x-matchPts[i].x, refPts[i].y-matchPts[i].y);  
  }

  cv::Mat matches_image;
  cv::drawMatches( refImage, keypointsA, matchImage, keypointsB,
                       good_matches, matches_image, cv::Scalar::all(-1), cv::Scalar::all(-1),
                       std::vector<char>(),cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
                       

  // Draw the transform H
  std::vector<cv::Point2f> ptsIn, ptsOut;
  for (int r=25; r<matchImage.rows; r+=400)
  {
    for (int c=25; c<matchImage.cols; c+=160)
    {
      cv::Point2f ptIn(c,r); // Point in the HRSC image
      ptsIn.push_back(ptIn);
    }
  }
  cv::perspectiveTransform(ptsIn, ptsOut, H);
  for (size_t i=0; i<matchPts.size(); i+=4)
  {
    cv::Point2f matchPt(matchPts[i] + cv::Point2f(refImage.cols, 0));
    cv::Point2f refPt(refPts[i]); // Point in the reference image
    //printf("OUT Pair: HRSC(%lf, %lf) ==> NOEL(%lf, %lf) \n", ptsIn[i].x, ptsIn[i].y, ptsOut[i].x, ptsOut[i].y);
    //std::cout << "RefPt = " << refPt << std::endl;
    cv::line(matches_image, matchPt, refPt, cv::Scalar(0, 255, 0), 1);
  }
                       
  //cv::imshow( "BRISK All Matches", matches_image );
  //cv::waitKey(0);

  //cv::IplImage* outrecog = new cv::IplImage(matches_image);
  cv::imwrite( "match_image.jpeg", matches_image );


  return true;
}



short interpPixel(const cv::Mat& img, const cv::Mat& mask, float xF, float yF)
{
    const int BORDER_SIZE = 1; // Stay away from border artifacts

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
    
    float v00 = static_cast<float>(img.at<short>(y0, x0));
    float v01 = static_cast<float>(img.at<short>(y0, x1));
    float v10 = static_cast<float>(img.at<short>(y1, x0));
    float v11 = static_cast<float>(img.at<short>(y1, x1));

    short val = static_cast<short>( v00*(1-a)*(1-c)  + v10*a*(1-c) + v01*(1-a)*c + v11*a*c );
    
/*
    short val = cvRound((img.at<short>(y0, x0) * (1.f - a) 
                      +  img.at<short>(y0, x1) * a) * (1.f - c)
                      + (img.at<short>(y1, x0) * (1.f - a)
                      +  img.at<short>(y1, x1) * a) * c);
                      
*/
    //printf("=== %d == %f, %f, %f, %f ===\n", val, v00, v01, v10, v11);


    return val;
}

/// Sharpen the reference image using an aligned match image
bool sharpenReferenceImage(const cv::Mat &refImageIn, const cv::Mat &matchImageIn,
                           const cv::Mat &transform)
{

  // Compute the unsharp mask of the match image
  const float UNSHARP_SCALING = 4.0f; // Control strength of the enhancement
  const float MIN_ENHANCE     = 15;      // Only enhance pixels over this amount (keep out noise)
  cv::Mat blurredMatch, temp, sharpMask;
  cv::GaussianBlur( matchImageIn, blurredMatch, cv::Size(5,5), 0, 0, cv::BORDER_DEFAULT ); 
  cv::subtract(matchImageIn, blurredMatch, sharpMask, cv::noArray(), CV_16S);
  
  double min, max;
  cv::minMaxIdx(sharpMask, &min, &max);
  printf("min = %lf, max = %lf\n", min, max);
  cv::Mat displaySharp = (sharpMask + 128);
  //cv::convertScaleAbs(sharpMask, displaySharp, UNSHARP_SCALING);//255.0/(0.2*max));
  cv::imwrite("sharpMask.jpeg", displaySharp);
    
  // TODO: Show the sharp mask to make sure it looks ok!
  // TODO: Play around with the enhancement amount!
    
  // Now sharpen the input image
  cv::Mat outputImage, outputMask;
  outputImage = refImageIn;
  for (int r=0; r<refImageIn.rows; ++r)
  {
    for (int c=0; c<refImageIn.cols; ++c)
    {     
      float matchX = c*transform.at<float>(0,0) + r*transform.at<float>(0,1) + transform.at<float>(0,2);
      float matchY = c*transform.at<float>(1,0) + r*transform.at<float>(1,1) + transform.at<float>(1,2);
      //printf("%d, %d  =>  %lf, %lf\n", r, c, matchY, matchX);
      
      short temp = interpPixel(sharpMask, matchImageIn, matchX, matchY);
      float sharpVal = static_cast<float>(temp)*UNSHARP_SCALING;
      if (abs(sharpVal) < MIN_ENHANCE)
        continue; // Skip small values
      // Enforce limits
      short oldVal   = static_cast<short>(outputImage.at<unsigned char>(r,c));
      short newVal   = oldVal + (short)sharpVal;
      if (temp > max)
      {
        printf("%d, %d  =>  %lf, %lf\n", r, c, matchY, matchX);
        printf("%d, %d, %f\n", oldVal, temp, sharpVal);
        return false;
      }
      if (newVal > 255)
        newVal = 255;
      if (newVal < 0)
        newVal = 0;
      outputImage.at<unsigned char>(r,c) = (unsigned char)newVal;
    }
  }
  
  cv::imwrite("sharpened_image.jpeg", outputImage);
  
  //cv::namedWindow("Display Image", cv::WINDOW_AUTOSIZE );
  //cv::imshow("Display Image", outputImage);
  //cv::waitKey(0);
  
  printf("Finished sharpening\n");
  

}

//=============================================================

int main(int argc, char** argv )
{
  
  if (argc < 3)
  {
    printf("usage: DisplayImage.out <Ref Image Path> <Match Image Path>\n");
    return -1;
  }
  std::string refImagePath   = argv[1];
  std::string matchImagePath = argv[2];
  
  // Load the input image
  cv::Mat refImageIn = cv::imread(refImagePath, 0);
  if (!refImageIn.data)
  {
    printf("Failed to load reference image\n");
    return -1;
  }

  cv::Mat matchImageIn = cv::imread(matchImagePath, 0);
  if (!matchImageIn.data)
  {
    printf("Failed to load match image\n");
    return -1;
  }

  // Display the input image
  //namedWindow("Display Image", WINDOW_AUTOSIZE );
  //imshow("Display Image", image);
  //waitKey(0);
  
  // First compute the transform between the two images
  // - This could be moved to a seperate tool!
  cv::Mat transform(3, 3, CV_32FC1);
  //computeImageTransform(refImageIn, matchImageIn, transform);
  
  
  transform.at<float>(0, 0) = 1;  // Skip registration for now!
  transform.at<float>(0, 1) = 0;
  transform.at<float>(0, 2) = 49.44;//24.38;
  transform.at<float>(1, 0) = 0;
  transform.at<float>(1, 1) = 1;
  transform.at<float>(1, 2) = 1191.51;//32.46;
  transform.at<float>(2, 0) = 0;
  transform.at<float>(2, 1) = 0;
  transform.at<float>(2, 2) = 1;
  
  
  // For the next step the transform needs to be from reference to match!
  cv::Mat invTransform(transform);
  double check = cv::invert(transform, invTransform);
  std::cout << "H_inv = \n" << invTransform << std::endl;
  
  sharpenReferenceImage(refImageIn, matchImageIn, invTransform);


  return 0;
}








