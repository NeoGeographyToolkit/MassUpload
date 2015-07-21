// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Filter.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>

#include <vector>

#include <HrscCommon.h>

#include <boost/program_options.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/foreach.hpp>
namespace po = boost::program_options;

using namespace vw;

/**
  This program takes in a set of images and returns a uint8 image
  which is 255 in pixels where all of the input images have data.

*/

//=================================================================================




// TODO: This could be handled with a generic class + a functor.
//      --> May it already exists somewhere?
/// Image view class which creates a binary mask image equalling
///  255 in locations where all of the input images are non-zero.
/// - The output mask is equal in size to the smallest input image.
template <class ImageT>
class ImageAndView : public ImageViewBase<ImageAndView<ImageT> >
{
private:
  std::vector<ImageT> m_image_vec;
  int m_num_rows;
  int m_num_cols;

public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::pixel_type result_type;
  
  typedef PixelGray<unsigned char> Uint8;

  // Constructor
  ImageAndView( std::vector<ImageT> & imageVec ) : m_image_vec(imageVec)
  {
    const size_t numImages = imageVec.size();
    VW_ASSERT( (numImages > 0), ArgumentErr() << "ImageAndView: One or more images required!." );

    // Determine the minimum image size
    m_num_rows = imageVec[0].rows();
    m_num_cols = imageVec[0].cols();
    for (size_t i=1; i<numImages; ++i)
    {
      if (imageVec[i].rows() < m_num_rows)  m_num_rows = imageVec[i].rows();
      if (imageVec[i].cols() < m_num_cols)  m_num_cols = imageVec[i].cols();
    }
  }

  inline int32 cols  () const { return m_num_cols; }
  inline int32 rows  () const { return m_num_rows; }
  inline int32 planes() const { return 1; }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const 
  { 
    return 0; // NOT IMPLEMENTED!
  }

  typedef ProceduralPixelAccessor<ImageAndView<ImageT> > pixel_accessor;
  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  typedef CropView<ImageView<result_type> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const 
  { 
    // Set up the output image tile
    ImageView<result_type> tile(bbox.width(), bbox.height());

    // Set up for pixel calculations
    size_t num_images = m_image_vec.size();

    // Rasterize all the input images at this particular tile
    std::vector<ImageView<Uint8> > input_tiles(num_images);
    for (size_t i=0; i<num_images; ++i)
    {
      // Make sure we don't go outside the bounds of this image
      BBox2i safe_bbox = bbox;
      safe_bbox.crop(bounding_box(m_image_vec[i]));
      input_tiles[i] = crop(m_image_vec[i], safe_bbox);
    }


    // Loop through each output pixel and compute each output value
    for (int c = 0; c < bbox.width(); c++)
    {
      for (int r = 0; r < bbox.height(); r++)
      {
        // Check all of the pixels at this location and perform the AND operation
        result_type thisPixel;
        result_type outputPixel;
        outputPixel[0] = 255; // UINT8 max
        for (size_t i=0; i<num_images; ++i) 
        {
          thisPixel = input_tiles[i](c,r);
          if (thisPixel[0] <= 0) {
            outputPixel[0] = 0;
            break; // Skip looking at the other pixels at this location
          }
        } // End image loop
        tile(c, r) = outputPixel;
      } // End row loop
    } // End column loop

  // Return the tile we created with fake borders to make it look the size of the entire output image
  return prerasterize_type(tile,
                           -bbox.min().x(), -bbox.min().y(),
                           cols(), rows() );

  } // End prerasterize function

 template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const 
 { 
   vw::rasterize( prerasterize(bbox), dest, bbox ); 
 }
  
}; // End class ImageAndView


struct Options {
  Options() : nodata(-1), feather_min(0), feather_max(255), filter("linear") {}
  // Input
  std::vector<std::string> input_files;

  // Settings
  double nodata;
  int feather_min, feather_max; // Currently these are always left at the defaults
  std::string filter;
  std::string output_filename;
};






// Handling input
void handle_arguments( int argc, char *argv[], Options& opt ) {

  po::options_description general_options("");
  general_options.add_options()
    ("nodata-value",      po::value(&opt.nodata), "Value that is nodata in the input image. Not used if input has alpha.")
    ("output-filename,o", po::value(&opt.output_filename), "Output file name.")
    ("help,h",            "Display this help message");

  po::options_description positional("");
  positional.add_options()
    ("input-files", po::value<std::vector<std::string> >(&opt.input_files));

  po::positional_options_description positional_desc;
  positional_desc.add("input-files", -1);

  po::options_description all_options;
  all_options.add(general_options).add(positional);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << general_options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <image-files>\n";
  boost::to_lower( opt.filter );

  if ( vm.count("help") )
    vw_throw( ArgumentErr() << usage.str() << general_options );
  if ( opt.input_files.empty() )
    vw_throw( ArgumentErr() << "Missing input files!\n" << usage.str() << general_options );

}

int main( int argc, char *argv[] ) {

  Options opt;
  handle_arguments( argc, argv, opt );

  // Use a default output path if none provided
  const std::string firstPath = opt.input_files[0];
  //vw_out() << "Loading: " << firstPath << "\n";
  size_t pt_idx = firstPath.rfind(".");
  std::string output;
  if (opt.output_filename.size() != 0)
    output = opt.output_filename;
  else {
    output = firstPath.substr(0,pt_idx)+"_union";
    output += firstPath.substr(pt_idx,firstPath.size()-pt_idx);
  }

  // Read the georef from the first file, they should all have the same value.
  cartography::GeoReference georef;
  cartography::read_georeference(georef, firstPath);

  // For now this is always the same
  typedef PixelGray<uint8> PixelT;
  
  const size_t numInputFiles = opt.input_files.size();
  std::vector< ImageViewRef<PixelT> > input_images(numInputFiles);

  // Loop through each input file
  for (size_t i=0; i<numInputFiles; ++i)
  {
    const std::string input = opt.input_files[i];

    // Determining the format of the input
    SrcImageResource *rsrc = DiskImageResource::open(input);
    ChannelTypeEnum channel_type = rsrc->channel_type();
    PixelFormatEnum pixel_format = rsrc->pixel_format();

    // Check for nodata value in the file
    if ( rsrc->has_nodata_read() ) {
      opt.nodata = rsrc->nodata_read();
      std::cout << "\t--> Extracted nodata value from file: " << opt.nodata << ".\n";
    }
    delete rsrc;

    DiskImageView<PixelT> this_disk_image(input);
    input_images[i] = this_disk_image;

  } // loop through input images


  vw_out() << "Writing: " << output << std::endl;
  block_write_gdal_image(output, 
                        ImageAndView< ImageViewRef<PixelT> >(input_images),
                        georef,
                        TerminalProgressCallback("bigMaskMaker","Writing:"));


}
