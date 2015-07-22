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
#include <vector>

#include <HrscCommon.h>

#include <boost/program_options.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/foreach.hpp>
namespace po = boost::program_options;

using namespace vw;

/**
  This program performs a grassfire operation on an input binary image.

*/




//======================================================================================================

struct Options {
  Options() : nodata(-1), feather_min(0), feather_max(255), filter("linear") {}
  // Input
  std::string input_file;

  // Settings
  double nodata;
  int feather_min, feather_max; // Currently these are always left at the defaults
  std::string filter;
  std::string output_filename;
};





template <class ImageT>
class GrassfireView : public ImageViewBase<GrassfireView<ImageT> >
{
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::pixel_type result_type;
  

  typedef PixelGray<unsigned char> Uint8;
  typedef PixelGray<unsigned char> Uint16;

private:
  DiskImageView<Uint8> const& m_input_image;
  int m_num_rows;
  int m_num_cols;

public:

  // Constructor
  GrassfireView( DiskImageView<Uint8> const& input_image ) : m_input_image(input_image)
  {
    m_num_rows = input_image.rows();
    m_num_cols = input_image.cols();
  }

  inline int32 cols  () const { return m_num_cols; }
  inline int32 rows  () const { return m_num_rows; }
  inline int32 planes() const { return 1; }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const 
  { 
    return 0; // NOT IMPLEMENTED!
  }

  typedef ProceduralPixelAccessor<GrassfireView<ImageT> > pixel_accessor;
  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  typedef CropView<ImageView<result_type> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const 
  { 
    const int GRASSFIRE_DISTANCE = static_cast<int>(MASK_MAX);
  
    // Load up a section of the input image that is large enough to fully compute the
    // grassfile results for this tile.  That ensures only a single disk load per tile.
    BBox2 expanded_bbox = bbox;
    expanded_bbox.expand(GRASSFIRE_DISTANCE);
    
    // Figure out the section of the expanded ROI that corresponds to the input bbox
    // - This is normally constant but can change if the expanded ROI gets cropped by
    //   the boundaries of the input image.
    int outputX = GRASSFIRE_DISTANCE;
    int outputY = GRASSFIRE_DISTANCE;
    if (expanded_bbox.min().x() < 0) outputX += expanded_bbox.min().x();
    if (expanded_bbox.min().y() < 0) outputY += expanded_bbox.min().y();
    BBox2 output_section(outputX, outputY, bbox.width(), bbox.height());

    // Crop the expanded bbox to the image bounds and fetch that image data from disk    
    expanded_bbox.crop(bounding_box(m_input_image));
    ImageView<MASK_DATA_TYPE> expanded_tile = pixel_cast<MASK_DATA_TYPE>(crop(m_input_image, expanded_bbox));
    /*
    vw_out() << "Input tile    : " << bbox << std::endl;
    vw_out() << "Expanded tile : " << expanded_bbox << std::endl;
    vw_out() << "output section: " << output_section << std::endl;
    */
    const typename ImageView<Uint16>::pixel_type zero = typename ImageView<Uint16>::pixel_type();
    if (max_pixel_value(crop(expanded_tile, output_section)) == zero)
    {
        //vw_out() << "Skipping zero tile!\n";
        // If all pixels in the input bbox are zero, just return a blank image tile.
        return prerasterize_type(crop(expanded_tile, output_section),
                                 -bbox.min().x(), -bbox.min().y(),
                                 cols(), rows() );  
    }
    //vw_out() << "Normal tile!\n";
   
    // Grassfire computation
    ImageView<MASK_DATA_TYPE> grassfire_output(pixel_cast<MASK_DATA_TYPE>(clamp(grassfire(expanded_tile), 0, GRASSFIRE_DISTANCE-1)));
    
    // Return the tile we created with fake borders to make it look the size of the entire output image
    // - Only return the image data from the input bbox, not the entire expanded bbox
    return prerasterize_type(crop(grassfire_output, output_section),
                             -bbox.min().x(), -bbox.min().y(),
                             cols(), rows() );

  } // End prerasterize function

 template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const 
 { 
   vw::rasterize( prerasterize(bbox), dest, bbox ); 
 }
  
}; // End class ImageAndView



//--------------------------------------------------------




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
    ("input-file", po::value<std::string>(&opt.input_file));

  po::positional_options_description positional_desc;
  positional_desc.add("input-file", -1);

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
  //if ( opt.input_file.empty() )
  //  vw_throw( ArgumentErr() << "Missing input files!\n" << usage.str() << general_options );

  // Set the system cache size
  vw_settings().set_system_cache_size( cache_size*1024*1024 );

}



int main( int argc, char *argv[] ) {

  Options opt;
  handle_arguments( argc, argv, opt );

  // Use a default output path if none provided
  vw_out() << "Loading: " << opt.input_file << "\n";
  size_t pt_idx = opt.input_file.rfind(".");
  std::string output;
  if (opt.output_filename.size() != 0)
    output = opt.output_filename;
  else {
    output = opt.input_file.substr(0,pt_idx)+"_union";
    output += opt.input_file.substr(pt_idx,opt.input_file.size()-pt_idx);
  }

  // Read the georef from the first file, they should all have the same value.
  cartography::GeoReference georef;
  cartography::read_georeference(georef, opt.input_file);

  // The input binary mask is always 8 bit
  typedef PixelGray<uint8> PixelT;

  // Determining the format of the input
  SrcImageResource *rsrc = DiskImageResource::open(opt.input_file);
  ChannelTypeEnum channel_type = rsrc->channel_type();
  PixelFormatEnum pixel_format = rsrc->pixel_format();

  // Check for nodata value in the file
  if ( rsrc->has_nodata_read() ) {
    opt.nodata = rsrc->nodata_read();
    std::cout << "\t--> Extracted nodata value from file: " << opt.nodata << ".\n";
  }
  delete rsrc;

  DiskImageView<PixelT> disk_image(opt.input_file);


  vw_out() << "Writing: " << output << std::endl;
  block_write_gdal_image(output, 
                        GrassfireView<ImageView<PixelGray<MASK_DATA_TYPE> > >(disk_image),
                        georef,
                        TerminalProgressCallback("bigMaskGrassfire","Writing:"));
                              
}
