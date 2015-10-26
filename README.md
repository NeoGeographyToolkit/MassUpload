==========================================
                MassUpload
==========================================

------------------------------------------------------------------------------------------
Code contents
-----------------------------

Tools for mass upload of data to mapping servers.

There are two main components to this code.  The first set of files was for batch 
uploading lots of satellite imagery to the defunct Google Maps Engine tool.

- prepThemisMosaic.py = Small script to help generate a THEMIS mosaic.
- initDatabase.txt = Describes the database used to store uploaded file information.
- fix_jp2.cpp = A tool for fixing HiRISE JP2 files with a metadata error.
- common.py = General python utilities

- mapsEngineUpload.py = Lower level tool for uploading images to Maps Engine plus some other functions.
- unifiedDataLoader.py = Higher level tool for uploading large batches of files to Maps Engine.

- ctxDataLoader.py - Tool for finding, fetching, and preprocessing CTX images.
- hiriseDataLoader.py - Tool for finding, fetching, and preprocessing HiRISE images.
- hrscDataLoader.py - Tool for finding, fetching, and preprocessing HRSC images.

- remoteInstallScript.sh = Instructions for setting up a new upload environment.
- sampleBash.sh = Sample bashrc script for a new upload environment.

These tools updated a local database file with information about available
satellite data, fetched that data, converted problem formats, and uploaded it to Maps Engine.
The top level tool for doing this is unifiedDataLoader.py, though with Maps Engine going away
the tool as a whole is unlikely to be useful in the future.



The second set of files is for making a high resolution mosaic of HRSC satellite 
imagery spanning all of Mars.

- CMakeLists.txt = Build script for C++ tools.
- FindVisionWorkbench.cmake = Find the required VW installation.
- HrscCommon.h = Common C++ functions.
- MosaicUtilities.py = Supporting Python classes.
- RegisterHrsc.cpp = Improve the estimated registration of an HRSC image on the base map.
- badHrscSets.csv = List of HRSC image ID's which either contain artifacts or are not handled properly.
- bigMaskGrassfire.cc = Produce a blending mask for entire HRSC images.
- bigMaskMaker.cc = Produce a binary mask for entire HRSC images.
- computeBrightnessCorrection.cc = Compute a simple brightness correction for HRSC images.
- hrscFileCacher.py = Fetch HRSC images as they are requested and hold on to the most recent ones.
- hrscImageManager.py = Coordinates the processing for an HRSC image to get it ready to paste on the map.
- hrscMosaic.cpp = Paste processed HRSC image tiles on top of the output map tiles.
- makeSimpleImageMask.cpp = Simple binary bask for small images.
- marsColorMosaicCreator.py = The main program for generating the HRSC map.
- mosaicTileManager.py = Manage the output map tiles.
- sendToGoogleBucket.py = Standalone tool for sending up data to a Google Bucket using gsutil.
- solveHrscColor.py = Given color pairs, generate information for a good looking color transform.
- stackImagePyramid.py = Generate a simple kml tree to display a low res version of the map.
- transformHrscImageColor.cpp = Using a computed transform, generate a basemap-colored HRSC map.
- writeHrscColorPairs.cpp = Write pixel pairs from the HRSC image and the basemap.

Major external dependencies: Vision Workbench, OpenCV, GDAL, ImageMagick.

------------------------------------------------------------------------------------------
Implementation details
-----------------------------

The map generation process is heavily reliant on the use of an existing full-mars color image
as a basemap to enhance.  HRSC data does not contain the normal RGB spectrum, so the color of the
final map is a blended composite of HRSC spectral data and the RGB data from the existing map.
The output map will be much higher resolution than the input map as specified in the program.
For example, the map can be generated at 15 meters per pixel, resulting in a 3 Terrabyte output map.
The output map is composed of many large image tiles, each of which corresponds to a specific region
in the input full map.

A brief overview of the process used to generate the output map:
1 - A database of HRSC images, their footprints, and where to find them is needed.
    - This can be generated using the Maps Engine upload portion of the project.
2 - Generate a high resolution, tiled, and blurred version of the input map.
    - This is done using imagemagick but that program has a tendency to crash!
3 - Fetch one HRSC image at a time to process.  The red, green, blue, nadir, and near-IR
    channels are used.
4 - At low resolution (matching the full input map):
    - Register the image to the input map
    - Compute brightness corrections
5 - At high resolution (the output resolution):
    - Make blending masks
    - Split up the image into tiles for processing (no relation to the output tiles)
    - Compute a color transform for each tile
    - Generate color transformed tiles
6 - After high resolution processing is finished, all of the processed HRSC tiles
    are blended on to the output map tiles.
7 - Steps 3 to 6 are repeated for every HRSC image (there are over 3000 of them).
8 - After the map is created a low resolution KML tree can be created for debugging purposes
    and uploaded to a Google cloud bucket if desired.
    
The git repository has a "poles" branch which contains a modified version of the code
used to generate seperate stereographic projection maps for the north and south poles
of Mars.  The process used for these maps is very similar but many parts of the code
had to be converted from performing calculations in degrees to performing calculations
in projected coordinates.  Ideally the two branches would be merged to allow the same code
to generate all three of the maps.

Known issues with the code that if corrected would significantly improve the output map:
- The full resolution image blend map is the biggest processing bottleneck and is not fully
  parallelized.
- The image registration algorithm fails about 20% of the time, leading to gaps in coverage.
- There is a distinct vignetting effect that makes many of the HRSC images stand out from
  the rest of the output map.
- The north pole mosaic has trouble with changing ice locations.
- Some input images produce strange color artifacts.  Those images may just have to be isolated
  and excluded from processing.
- It would be nice to have a streamlined method of merging the main and polar mosaics.

------------------------------------------------------------------------------------------
Map projection details
-----------------------------

All of the maps generated by the tool are based on an existing full-coverage RGB map of mars.
The projection space used for this map is:

PROJ.4 string is:
'+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=3396200 +b=3376200 +units=m +no_defs '
Origin = (-10669477.099999999627471,5334738.599999999627471)
Pixel Size = (1852.339774305555466,-1852.339791666666542)

The repository is set up to generate three output maps: one equirectangular map for the
+/- 60 degree latitude region and one stereographic map for each pole of Mars.
The main map corresponds to the low resolution input map but is at 128x resolution.
It is broken up into tiles, each corresponding to a 45x45 pixel region of the input map,
meaning that there are 128 rows and 256 columns in the full map (fewer in the +/-60 degree region).
The output tiles can be assigned the correct georegistration data based on this information.
The output resolution is flexible if changed the blend settings used in some of the C++ files may
need to be updated.

In order to generate the two polar maps, polar versions of the input map had to be generated.
This was done by cropping the top and bottom sections of the map and using gdalwarp similar to
the following command:
gdalwarp map_upper_crop.tif map_north_pole.tif -t_srs "+proj=stere +lon_0=0 +lat_0=90 +a=3396200 +b=3376200 units=m" -wo sample_grid=yes, -wo sample_steps=1000, -wo source_extra=2 -wo CUTLINE_ALL_TOUCHED=true -r cubic

The projection info for the polar maps are:

PROJ.4 string is:
'+proj=stere +lat_0=90 +lat_ts=90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=3396200 +b=3376200 +units=m +no_defs '
Origin = (-1959288.086665658513084,1959261.357437750091776)
Pixel Size = (1852.339774305555466,-1852.339791666666542)

PROJ.4 string is:
'+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=3396200 +b=3376200 +units=m +no_defs '
Origin = (-1959439.740240827202797,1959418.630283087957650)
Pixel Size = (1852.339774305555466,-1852.339791666666542)

The polar maps were cropped to be exactly 2115x2115 in size so they could be evenly divided into 45x45
pixel tiles in the same manner as the main RGB map.  Since the polar images are larger than required,
losing some pixels at the edge is not an issue.  The output tiles are also generated at 128x
resolution.  The output tiles can be assigned the correct georegistration data based on this information.

In order to generate a single map, the polar map tiles will need to be converted into the same projection
as the main +/-60 degree map.












