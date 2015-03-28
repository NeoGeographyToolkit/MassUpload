
import os
import sys
import IrgGeoFunctions

"""
TODO:



- Generate an image 
- Algorithm to smooth out transitions.
- Generate a good result image.
- Generate a list of features Earth Engine would need to replicate all steps.

- Better mask handling

Existing tools:
- RegisterHrsc.cpp
    - Input  = Basemap, HRSC
    - Output = spatialTransform
- writeHrscColorPairs.cpp
    - Input  = Basemap, HRSC, spatialTransform
    - Output = File containing pixel color pairs
- transformHrscImageColor.cpp
    - Input  = HRSC, colorTransform
    - Output = Color transformed HRSC image
    - TODO   = Add cleanup/pansharp


"""

fullBasemapPath    = '/home/smcmich1/data/hrscMapTest/noel_basemap.tif'
cropBasemapPath = '/home/smcmich1/data/hrscMapTest/basemap_crop.tif'
hrscInputPaths = ['/home/smcmich1/data/hrscMapTest/h0022_0000_re3.tif',
                  '/home/smcmich1/data/hrscMapTest/h0022_0000_gr3.tif',
                  '/home/smcmich1/data/hrscMapTest/h0022_0000_bl3.tif',
                  '/home/smcmich1/data/hrscMapTest/h0022_0000_ir3.tif',
                  '/home/smcmich1/data/hrscMapTest/h0022_0000_nd3.tif']
HRSC_RED   = 0
HRSC_GREEN = 1
HRSC_BLUE  = 2
HRSC_NIR   = 3
HRSC_NADIR = 4


GDAL_DIR = '/home/smcmich1/programs/gdal-1.11.0-install/bin/'

print 'Starting basemap enhancement script...'

# Get the HRSC bounding box and expand it
HRSC_BB_EXPAND_DEGREES = 1.5
(minLon, maxLon, minLat, maxLat) = IrgGeoFunctions.getGeoTiffBoundingBox(hrscInputPaths[4])
minLon -= HRSC_BB_EXPAND_DEGREES
maxLon += HRSC_BB_EXPAND_DEGREES
minLat -= HRSC_BB_EXPAND_DEGREES
maxLat += HRSC_BB_EXPAND_DEGREES

print str((minLon, maxLon, minLat, maxLat))


# Convert the bounding box from degrees to the projected coordinate system (meters)
DEGREES_TO_PROJECTION_METERS = 59274.9
NOEL_MAP_METERS_PER_PIXEL = 1852.4
minX = minLon*DEGREES_TO_PROJECTION_METERS
maxX = maxLon*DEGREES_TO_PROJECTION_METERS
minY = minLat*DEGREES_TO_PROJECTION_METERS
maxY = maxLat*DEGREES_TO_PROJECTION_METERS

# Crop out the correct section of the base map
projCoordString = '%f %f %f %f' % (minX, maxLat, maxX, minY)
if False:#not os.path.exists(redCroppedPath):
    #cmd = (GDAL_DIR+'gdal_translate ' + redBasemapPath +' '+ redCroppedPath
    #                         +' -projwin '+ projCoordString)
    cmd = (GDAL_DIR+'gdal_translate ' + fullBasemapPath +' '+ redCroppedPath
                             +' -projwin '+ projCoordString)
    print cmd
    os.system(cmd)

# Increase the resolution of the cropped image
# TODO: Can this be done in one step?
RESOLUTION_INCREASE = 200 # In percent
if False:#not os.path.exists(redCropHighResPath):
    cmd = (GDAL_DIR+'gdal_translate ' + redCroppedPath +' '+ redCropHighResPath
                             +' -outsize '+str(RESOLUTION_INCREASE)+'% '+str(RESOLUTION_INCREASE)+'% ')
    print cmd
    os.system(cmd)



def warpHrscFile(sourcePath, metersPerPixel):
    warpedPath = sourcePath[:-4] + '_resample.tif'
    if not os.path.exists(warpedPath):
        cmd = (GDAL_DIR+'gdalwarp ' + sourcePath +' '+ warpedPath + ' -r cubicspline '
                 ' -t_srs "+proj=eqc +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m" -tr '
                 + str(metersPerPixel)+' '+str(metersPerPixel)+' -overwrite')
        print cmd
        os.system(cmd)
    return warpedPath

# Transform the HRSC image to the same projection/resolution as the upsampled base map crop
highResMetersPerPixel = NOEL_MAP_METERS_PER_PIXEL / (RESOLUTION_INCREASE/100.0)
hrscWarpedPaths       = [warpHrscFile(path, highResMetersPerPixel) for path in hrscInputPaths]


hrscPathString = ''
for path in hrscWarpedPaths:
    hrscPathString += path + ' '

# TODO: Don't hardcode this
spatialTransformPath = '/home/smcmich1/data/hrscMapTest/h0022_0000_spatial_transform.csv'

# Generate the color pairs
colorPairPath = '/home/smcmich1/data/hrscMapTest/h0022_0000_color_pairs.csv'
if not os.path.exists(colorPairPath):
    cmd = './writeHrscColorPairs ' + cropBasemapPath +' '+ hrscPathString +' '+ spatialTransformPath +' '+ colorPairPath
    print cmd
    os.system(cmd)


# Compute the color transform
colorTransformPath = '/home/smcmich1/data/hrscMapTest/h0022_0000_color_transform.csv'
if not os.path.exists(colorTransformPath):
    cmd = 'python /home/smcmich1/repo/MassUpload/solveHrscColor.py ' + colorTransformPath +' '+ colorPairPath 
    print cmd
    os.system(cmd)


# Transform the HRSC image color
hrscNewColorPath = '/home/smcmich1/data/hrscMapTest/h0022_0000_new_color.tif'
if not os.path.exists(hrscNewColorPath):
    cmd = './transformHrscImageColor ' + cropBasemapPath +' '+ hrscPathString +' '+ colorTransformPath +' '+ hrscNewColorPath
    print cmd
    os.system(cmd)


mosaicPath = '/home/smcmich1/data/hrscMapTest/outputMosaic.tif'
if True:#not os.path.exists(mosaicPath):
    cmd = './hrscMosaic ' + cropBasemapPath +' '+ mosaicPath +' '+ hrscNewColorPath +' '+ spatialTransformPath
    print cmd
    os.system(cmd)


raise Exception('DEBUG')

print 'Basemap enhancement script completed!'

"""
== pansharp prototype ==

Generate RED version of Noel's map (ONCE)
Equitorial circumference = 21338954.25548 / 2 = 10669477.1 <--- x span
Polar      circumference = 21338954.25548 / 4 =  5334738.6 <--- y span
gdal_translate RA_Albedo.tif noel_basemap.tif
               -a_srs "+proj=eqc +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m"
               -a_ullr -10669477.1 5334738.6 10669477.1 -5334738.6
gdal_translate noel_basemap.tif noel_basemap_red.tif -b 1
DONE

FETCH HRSC IMAGE

Pick out HRSC image of interest and download ND3 file
DONE

PREP INPUT IMAGES

--> TODO: Try these steps out a increased resolution from the base map!
            2x or 4x should be good.

Generate a reduced-resolution of the file to match Noel's map using GDAL_translate.
- Noel's map resolution is (circumference 21,344 km/ 11520 pixels) = 1852.4 meters per pixel.

Get a bounding box surrounding the HRSC image and extract that region from Noel's map.

Min latitude  = -50.1955
Max latitude  = -18.6077
Min longitude = 98.2023
Max longitude = 102.64

Convert from lat/lon to projected coordinates by multiplying degrees with these constants:
meters per degree lon = 59274.9
meters per degree lat = 59274.9

gdal_translate noel_basemap_red.tif red_crop.tif -projwin 97.2023 -17.6077 103.64 -51.1955
 5761656.61227 -1043694.65673  6143250.636 -3034608.14295


Transform the HRSC image to the same projection/resolution as the base map
~/programs/gdal-1.11.0-install/bin/gdalwarp h0022_0000_nd3.tif h0022_0000_nd3_resample.tif -t_srs "+proj=eqc +lat_ts=0 +lat_0=0 +a=3396200 +b=3376200 units=m"  -tr 1852  1852 -overwrite


ALIGN THE IMAGES

Compute a pixel to pixel match between those two images.
- HRSC to NOEL => 24.384001, 32.455997 (basemap res)


GENERATE SHARPENED BASE TILE

Generate an unsharp image from the HRSC pan image.

Iterate through the base tile and for each pixel use the transform to find the matching point
 in the unsharp image, then add the unsharp image value to the base tile.  Also create a mask
 showing which pixels were modified (many won't overlap with the HRSC image).






"""



"""
DB command to find a list of overlapping HRSC files:
select setname from Files where sensor=1 and minLon<102.64 and maxLon>98.2023 and minLat<-18.6077 and maxLat>-50.1955 and subtype="nd3";


Alignment target = /byss/mars/themis/dayir_100m/
- Noel's map goes down to +/-90 degrees lat, but matching will gets very hard outside 60 degrees!
- Alignment between Noel and Themis is close but not exact, definately possible with ASP tools.



A professionally made partial mosaic: http://maps.planet.fu-berlin.de/
    - This is a hard problem, our only advantages are Earth Engine and the color reference map!

Query sqlite3 database to get a list of all the HRSC data sets after a certain date
    3500 sets -> 18000 files!

From each of these, we will need the following:
    ir3
    nd3
    red
    green
    blue

Align the nd3 image to the THEMIS Day Time global map.
    Vision workbench?
    
Generate a grid of pixel tie points between the nd3 and THEMIS

Generate a grid of pixel tie points between the mars color reference and THEMIS
    Only do this once!

Make a list of pixel pairs: HRSC - COLOR REF

Compute 5x3 transform using pixel pairs

Transfrom HRSC bands to a new RGB representation using the transform

Apply any image enhancement algorithms

Transform the RGB image to align with the THEMIS base map

While we are at it, try applying PAN unsharp mask to the upsampled color base layer!
    - This might be the only way to easily make a good mosaic


OUTSTANDING ISSUES
- Dark side color handling?
    - Ensure fixed brightness across image?
        - What about dark spots and poles?
        - Maybe brightness profile changed to match the color reference map?
    - Transform changes gradually down image?
        - Think of color transform as a rotation that smoothly changes
        - Can this be done between images as well?
    HRSC_h0344_0000_nd3
    HRSC_h0543_0000_nd3
    HRSC_h2011_0001_nd3

- Streak handling and other random crap
    HRSC_h0334_0000_nd3
    HRSC_h1970_0000_nd3

- Alignment verification

- Seam blending --> Making sure the colors match!
    - Images are captured at all sorts of light levels and angles so even a grayscale
      mosaic is hard to get right (see http://maps.planet.fu-berlin.de/)


- Can we do everything in OpenCV?
    - The largest file is 13,500x150,000
    - Image alignment is done at a lower resolution so no problem.
    - Everything else can be done in dumb tiles.


List of overlapping images for h0022_0000_nd3:

Use these:
h0506_0000_nd3 <-- Feels like there should be overlap =)
h2411_0000_nd3 <-- Lots of overlap!
h6419_0000_nd3 <-- Higher resolution shot of middle of image


h0248_0000_nd3
h0300_0000_nd3
h0440_0000_nd3 <-- Big image, but no overlap =(
h0451_0000_nd3 <-- Messed up top!
h0462_0000_nd3
h0506_0000_nd3 <-- Feels like there should be overlap =)
h0637_0000_nd3
h0648_0000_nd3
h1592_0000_nd3 <-- Weird polar image
h1774_0000_nd3
h1786_0001_nd3
h1942_0000_nd3 <-- Small image in the center
h2345_0000_nd3
h2389_0000_nd3
h2400_0001_nd3 <-- Lots of overlap, but not all components =()
h2411_0000_nd3 <-- Lots of overlap!
h2466_0000_nd3 
h2510_0001_nd3
h2619_0000_nd3
h2652_0000_nd3
h2663_0001_nd3
h2726_0000_nd3
h2729_0000_nd3
h4261_0000_nd3 <-- Maybe aligned?
h4272_0000_nd3
h4294_0000_nd3
h4469_0001_nd3
h4642_0000_nd3
h4817_0000_nd3
h6411_0000_nd3
h6419_0000_nd3 <-- Higher resolution shot of middle of image
h6437_0000_nd3
h8429_0000_nd3
h8474_0000_nd3
h8562_0000_nd3
h8990_0000_nd3
ha526_0000_nd3
ha688_0000_nd3
ha776_0000_nd3
hc598_0017_nd3

"""




