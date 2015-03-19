

#include <stdio.h>
#include <opencv2/opencv.hpp>

/**
  Program to compute the transform between an HRSC image and a base image using OpenCV
*/

int main(int argc, char** argv )
{
  
  if (argc != 3)
  {
    printf("usage: DisplayImage.out <Ref Image Path> <Match Image Path>\n");
    return -1;
  }
  std::string refImagePath   = argv[1];
  std::string matchImagePath = argv[2];
  
  // Load the input image
  cv::Mat refImageIn = cv::imread(refImagePath, 1);
  if (!refImageIn.data)
  {
    printf("Failed to load reference image\n");
    return -1;
  }
  cv::Mat matchImageIn = cv::imread(matchImagePath, 1);
  if (!matchImageIn.data)
  {
    printf("Failed to load match image\n");
    return -1;
  }
  // Display the input image
  //namedWindow("Display Image", WINDOW_AUTOSIZE );
  //imshow("Display Image", image);

  //waitKey(0);
  
  
  
  // Preprocess the images to improve feature detection
  cv::Mat refImage, matchImage, temp;
  //GaussianBlur( src, src, Size(3,3), 0, 0, BORDER_DEFAULT );
  //cv::GaussianBlur( matchImageIn, matchImageIn, cv::Size(3,3), 0, 0, cv::BORDER_DEFAULT );
  
  
  const int kernel_size = 5;
  const int scale = 1;
  const int delta = 0;
  cv::Laplacian( refImageIn, temp, CV_16S, kernel_size, scale, delta, cv::BORDER_DEFAULT );
  cv::convertScaleAbs( temp, refImage );
  cv::Laplacian( matchImageIn, temp, CV_16S, kernel_size, scale, delta, cv::BORDER_DEFAULT );
  cv::convertScaleAbs( temp, matchImage );

  
  std::vector<cv::KeyPoint> keypointsA, keypointsB;
  cv::Mat descriptorsA, descriptorsB;  

  cv::Ptr<cv::FeatureDetector    > detector  = cv::ORB::create();
  cv::Ptr<cv::DescriptorExtractor> extractor = cv::BRISK::create();

  detector->detect(  refImage, keypointsA);
  extractor->compute(refImage, keypointsA, descriptorsA);

  detector->detect(  matchImage, keypointsB);
  extractor->compute(matchImage, keypointsB, descriptorsB);

  // Rule out obviously bad matches
  cv::Mat mask(keypointsA.size(), keypointsB.size(), CV_8UC1);
  const float MAX_MATCH_PIXEL_DISTANCE = 200; // Tune this for the base map resolution!
  for (size_t i=0; i<keypointsA.size(); ++i)
  {
    for (size_t j=0; j<keypointsB.size(); ++j)
    {
      float distance = cv::norm(keypointsA[i].pt - keypointsB[j].pt);
      if (distance < MAX_MATCH_PIXEL_DISTANCE)
        mask.at<uchar>(i,j) = 1;
      else // Too far, disallow a match
        mask.at<uchar>(i,j) = 0;
    } 
  }

  //cv::FlannBasedMatcher matcher;
  cv::Ptr<cv::DescriptorMatcher> matcher = cv::DescriptorMatcher::create("BruteForce-Hamming");
  std::vector<cv::DMatch> matches;
  matcher->match(descriptorsA, descriptorsB, matches, mask);



  //-- Quick calculation of max and min distances between keypoints
  double max_dist = 0; double min_dist = 100;
  for( int i = 0; i < descriptorsA.rows; i++ )
  { double dist = matches[i].distance;
    if( dist < min_dist ) min_dist = dist;
    if( dist > max_dist ) max_dist = dist;
  }

  printf("-- Max dist : %f \n", max_dist );
  printf("-- Min dist : %f \n", min_dist );

  //-- Pick out "good" matches
  std::vector< cv::DMatch > good_matches;
  for( int i = 0; i < descriptorsA.rows; i++ )
  { if( matches[i].distance < 200 )
     { good_matches.push_back( matches[i]); }
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
  }
  cv::Mat H = cv::findHomography( matchPts, refPts, cv::RANSAC );
  std::cout << "H = \n" << H << std::endl;
   

  cv::Mat matches_image;
  cv::drawMatches( refImage, keypointsA, matchImage, keypointsB,
                       good_matches, matches_image, cv::Scalar::all(-1), cv::Scalar::all(-1),
                       std::vector<char>(),cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
                       
                       /*
  // Draw the transform H
  std::vector<cv::Point2f> ptsIn, ptsOut;
  for (int r=25; r<matchImage.rows; r+=50)
  {
    for (int c=25; c<matchImage.cols; c+=50)
    {
      cv::Point2f ptIn(c,r);
      ptsIn.push_back(ptIn);
    }
  }
  cv::perspectiveTransform(ptsIn, ptsOut, H);
  for (size_t i=0; i<ptsIn.size(); ++i)
  {
    cv::Point2f matchPt(ptsIn[i] + cv::Point2f(refImage.cols, 0));
    cv::Point2f refPt(ptsOut[i]);
    std::cout << "RefPt = " << refPt << std::endl;
    cv::line(matches_image, matchPt, refPt, cv::Scalar(0, 255, 0), 2);
  }*/
                       
  cv::imshow( "BRISK All Matches", matches_image );
  cv::waitKey(0);

  //cv::IplImage* outrecog = new cv::IplImage(all_matches);
  //cv::imwrite( "BRISK_All_Matches.jpeg", all_matches );


  


  return 0;
}








