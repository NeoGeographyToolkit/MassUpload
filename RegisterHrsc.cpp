

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <vw/Math/Geometry.h>
#include <vw/Math/RANSAC.h>
#include <vw/InterestPoint/InterestData.h>

#include <HrscCommon.h>

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
  float goodDist = 110;
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

//=============================================================




int main(int argc, char** argv )
{
  
  if (argc != 4)
  {
    printf("usage: RegisterHrsc <Base map path> <HRSC path> <Output path>\n");
    return -1;
  }
  std::string refImagePath   = argv[1];
  std::string matchImagePath = argv[2];
  std::string outputPath     = argv[3];
  
  const int LOAD_GRAY = 0;
  const int LOAD_RGB  = 1;
  
  // Load the input image
  cv::Mat refImageIn = cv::imread(refImagePath, LOAD_RGB);
  if (!refImageIn.data)
  {
    printf("Failed to load reference image\n");
    return -1;
  }

  cv::Mat matchImageIn = cv::imread(matchImagePath, LOAD_GRAY);
  if (!matchImageIn.data)
  {
    printf("Failed to load match image\n");
    return -1;
  }

  // First compute the transform between the two images
  // - This could be moved to a seperate tool!
  cv::Mat transform(3, 3, CV_32FC1);
  //computeImageTransform(refImageIn, matchImageIn, transform);
  
  
  transform.at<float>(0, 0) = 1;  // Skip registration for now!
  transform.at<float>(0, 1) = 0;
  transform.at<float>(0, 2) = 229-147;//49.44;//24.38;
  transform.at<float>(1, 0) = 0;
  transform.at<float>(1, 1) = 1;
  transform.at<float>(1, 2) = 3084-1892;//1191.51;//32.46;
  transform.at<float>(2, 0) = 0;
  transform.at<float>(2, 1) = 0;
  transform.at<float>(2, 2) = 1;
  
  
  // For the next step the transform needs to be from reference to match!
  cv::Mat invTransform(transform);
  double check = cv::invert(transform, invTransform);
  std::cout << "H_inv = \n" << invTransform << std::endl;

  writeTransform(outputPath, invTransform);
  
  return 0;
}








