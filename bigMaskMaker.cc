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
// - Transfer functions that should be in a VW file somewhere

// Linear Transfer Function
class LinearTransFunc : public vw::UnaryReturnSameType {
public:
  LinearTransFunc() {}

  template <class ArgT>
  ArgT operator()( ArgT const& value ) const { return value; }
};


// Cosine Transfer Function (this tracks 180 degrees with level at high and low)
template <class PixelT>
class CosineTransFunc : public vw::UnaryReturnSameType {
  typedef typename CompoundChannelType<PixelT>::type channel_type;
  typedef ChannelRange<channel_type> range_type;
public:
  CosineTransFunc() {}

  template <class ArgT>
  inline typename boost::enable_if<typename boost::is_floating_point<ArgT>,ArgT>::type
  operator()( ArgT const& value ) const {
    return range_type::max()*((1.0-cos(value/float(range_type::max())*M_PI))/2.0);
  }

  template <class ArgT>
  inline typename boost::disable_if<typename boost::is_floating_point<ArgT>,ArgT>::type
  operator()( ArgT const& value ) const {
    ArgT result = ArgT(range_type::max()*((1.0-cos(float(value)/float(range_type::max())*M_PI))/2.0));
    if ( result == 0 && value != 0 )
      result = 1;
    return result;
  }
};

// 90 degree Cosine transfer function ( high slope at beginning and low slope at end )
template <class PixelT>
class Cosine90TransFunc : public vw::UnaryReturnSameType {
  typedef typename CompoundChannelType<PixelT>::type channel_type;
  typedef ChannelRange<channel_type> range_type;
public:
  Cosine90TransFunc() {}

  template <class ArgT>
  inline typename boost::enable_if<typename boost::is_floating_point<ArgT>,ArgT>::type
  operator()( ArgT const& value ) const {
    return range_type::max()*(-cos(value/float(range_type::max())*(M_PI/2.0) + M_PI/2.0));
  }

  template <class ArgT>
  inline typename boost::disable_if<typename boost::is_floating_point<ArgT>,ArgT>::type
  operator()( ArgT const& value ) const {
    ArgT result = ArgT(range_type::max()*(-cos(float(value)/float(range_type::max())*(M_PI/2.0) + M_PI/2.0)));
    if ( result == 0 && value != 0 )
      result = 1;
    return result;
  }
};

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

    // Loop through each output pixel and compute each output value
    for (int c = 0; c < bbox.width(); c++)
    {
      for (int r = 0; r < bbox.height(); r++)
      {
        // Check all of the pixels at this location and perform the AND operation
        result_type thisPixel;
        result_type outputPixel;
        outputPixel[0] = 255;
        for (size_t i=0; i<num_images; ++i) 
        {
          thisPixel = m_image_vec[i](c,r);
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


//======================================================================================================

/*

// Function for highlighting spots of data???
template<class PixelT>
class NotNoDataFunctor {
  typedef typename CompoundChannelType<PixelT>::type channel_type;
  channel_type m_nodata;
  typedef ChannelRange<channel_type> range_type;
public:
  NotNoDataFunctor( channel_type nodata ) : m_nodata(nodata) {}

  template <class Args> struct result {
    typedef channel_type type;
  };

  inline channel_type operator()( channel_type const& val ) const {
    if (val == m_nodata)
      return range_type::min();
    else // data
      return range_type::max();
  }
};
//???
template <class ImageT, class NoDataT>
UnaryPerPixelView<ImageT,UnaryCompoundFunctor<NotNoDataFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type>  >
inline notnodata( ImageViewBase<ImageT> const& image, NoDataT nodata ) {
  typedef UnaryCompoundFunctor<NotNoDataFunctor<typename ImageT::pixel_type>, typename ImageT::pixel_type> func_type;
  func_type func( nodata );
  return UnaryPerPixelView<ImageT,func_type>( image.impl(), func );
}

*/

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




// Operation code for data that uses nodata
template <class PixelT>
void grassfire_nodata( Options& opt,
                       cartography::GeoReference const& georef,
                       ImageView<PixelT> const& input_image,
                       std::string output_path ) 
{
/* DEBUG - No grassfire!
  vw_out() << "Writing: " << output_path << std::endl;
  cartography::write_georeferenced_image(output_path, 
                                         input_image, 
                                         georef,
                                         TerminalProgressCallback("bigMaskMaker","Writing:"));
  return;
*/

  typedef typename CompoundChannelType<PixelT>::type inter_type;
  typedef          ChannelRange<inter_type>          range_type;

  //ImageView<int32> distance = grassfire(notnodata(input_image, inter_type(opt.nodata)));
  //ImageView<int32> distance = grassfire(input_image);

  // Check to see if the user has specified a feather length.  If not,
  // then we send the feather_max to the max pixel value (which
  // results in a full grassfire blend all the way to the center of the image.)
  //if (opt.feather_max < 1)
  //  opt.feather_max = max_pixel_value(distance);
  vw_out() << "\t--> Distance range: [ " << opt.feather_min << " " << opt.feather_max << " ]\n";

/*
  ImageViewRef<inter_type> norm_dist;
  norm_dist = pixel_cast<inter_type>(range_type::max() / (opt.feather_max - opt.feather_min) *
                                     clamp(pixel_cast<float>(distance) - opt.feather_min,
                                           0.0, opt.feather_max - opt.feather_min));
*/

/*
  ImageViewRef<typename PixelWithAlpha<PixelT>::type> result;
  if        ( opt.filter == "linear"   ) { result = per_pixel_filter(norm_dist, LinearTransFunc              ());
  } else if ( opt.filter == "cosine"   ) { result = per_pixel_filter(norm_dist, CosineTransFunc  <inter_type>());
  } else if ( opt.filter == "cosine90" ) { result = per_pixel_filter(norm_dist, Cosine90TransFunc<inter_type>());
  } else {
    vw_throw( ArgumentErr() << "Unknown transfer function " << opt.filter );
  }
*/
  vw_out() << "Writing: " << output_path << std::endl;
  cartography::write_georeferenced_image(output_path, 
                                         //per_pixel_filter(norm_dist, LinearTransFunc()), 
                                         pixel_cast<PixelT>(clamp(grassfire(input_image), 0, 255)),
                                         georef,
                                         TerminalProgressCallback("bigMaskMaker","Writing:"));



/*
  GdalOptions gdal_options;

  // Set up the output image writer
  boost::scoped_ptr<DiskImageResourceGDAL> output_handle( build_gdal_rsrc( output_file,
                                                                           difference, gdal_options ) );
  //output_handle->set_nodata_write( opt.nodata_value );
  write_georeference( *output_handle, output_georef );

  // Do all the computations and write to file
  block_write_image( *output_handle, 
                      result,
                     TerminalProgressCallback("=)", "\t--> Generating Mask ") );

*/

}




// Handling input
void handle_arguments( int argc, char *argv[], Options& opt ) {
  size_t cache_size;

  po::options_description general_options("");
  general_options.add_options()
    ("nodata-value",      po::value(&opt.nodata), "Value that is nodata in the input image. Not used if input has alpha.")
    ("output-filename,o", po::value(&opt.output_filename), "Output file name.")
    ("cache",             po::value(&cache_size)->default_value(1024), "Source data cache size, in megabytes.")
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

  // Set the system cache size
  vw_settings().set_system_cache_size( cache_size*1024*1024 );

}

/*
// The GDAL output image options format
typedef vw::DiskImageResourceGDAL::Options GdalOptions;

/// Helper function to set up a GDAL image for writing
template <class ImageT>
vw::DiskImageResourceGDAL*
build_gdal_rsrc( const std::string &filename,
                 vw::ImageViewBase<ImageT> const& image,
                 GdalOptions const& gdal_options ) {
  return new vw::DiskImageResourceGDAL(filename, image.impl().format(), opt.raster_tile_size, gdal_options);
}
*/





int main( int argc, char *argv[] ) {

  Options opt;
  handle_arguments( argc, argv, opt );

  // Use a default output path if none provided
  const std::string firstPath = opt.input_files[0];
  vw_out() << "Loading: " << firstPath << "\n";
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

  // Pass the input images to the working function
  //ImageViewRef<PixelT> binary_image = ImageAndView(input_images);
  //grassfire_nodata<PixelGray<uint8> >(opt, georef, binary_image, output);

  // Convert the input images to a single binary mask and then pass to the grassfire function
  // - This call also writes the data to disk
  //ImageViewRef<PixelT> temp =  ImageAndView< ImageViewRef<PixelT> >(input_images);
  //grassfire_nodata<PixelT >(opt, georef, temp, output);
  grassfire_nodata<PixelT >(opt, georef, ImageAndView< ImageViewRef<PixelT> >(input_images), output);


/*
  switch (pixel_format) 
  {
    case VW_CHANNEL_UINT8:  grassfire_nodata<PixelGray<uint8  > >(opt, input_images, output); break;
    case VW_CHANNEL_INT16:  grassfire_nodata<PixelGray<int16  > >(opt, input, output); break;
    case VW_CHANNEL_UINT16: grassfire_nodata<PixelGray<uint16 > >(opt, input, output); break;
    default:                grassfire_nodata<PixelGray<float32> >(opt, input, output); break;
  }; // End switch
*/
/*

    GdalOptions gdal_options;

    // Set up the output image writer
    boost::scoped_ptr<DiskImageResourceGDAL> output_handle( build_gdal_rsrc( output_file,
                                                                             difference, gdal_options ) );
    //output_handle->set_nodata_write( opt.nodata_value );
    write_georeference( *output_handle, output_georef );

    // Do all the computations and write to file
    block_write_image( *output_handle, 
                        TODO_GIANT_IMAGE_EXPRESSION,
                       TerminalProgressCallback("=)", "\t--> Generating Mask ") );



*/










}
