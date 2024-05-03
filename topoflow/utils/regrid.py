
#  Copyright (c) 2019-2023, Scott D. Peckham
#
#  Aug. 2023.  K to C in create_rts_from_nc_files().
#  August, September, October, Nov. 2019
#  Jun. 2020.  Updated create_rts_from_nc_files() to suppport
#              non-precip variables.
#  July 2021.  Updates to CHIRPS routines.
#
#-------------------------------------------------------------------

#  test1()
#  test2()

#  download_data()       # commented out
#  get_raster_bounds()
#  bounds_disjoint()
#  get_raster_cellsize()

#  read_geotiff()
#  regrid_geotiff()
#  clip_geotiff()
#  save_grid_to_geotiff()      2019-10-31
#  regrid_rts_file()

#-----------------------
#  read_nc_grid()           # (unused)
#  gdal_open_nc_file()      # (unused)
#-----------------------

#  fix_gpm_raster_bounds()
#  fix_gpm_file_as_geotiff
#  gdal_regrid_to_dem_grid()
#  create_rts_from_nc_files()       ### Create a rainfall grid stack
#  create_rts_from_chirps_files()

#-----------------
#  Commented out
#-----------------
#  regrid_geotiff_to_dem()    # (special case of regrid_geotiff())
#  resave_grid_to_geotiff()

#-------------------------------------------------------------------
import numpy as np
from osgeo import gdal, osr  ## ogr
import glob, sys
import os, os.path
import netCDF4 as nc
from datetime import datetime, timedelta

from . import rti_files
from . import rtg_files
from . import ncgs_files
from . import time_utils
from . import import_grid
from . import met_utils

## from . import rts_files   # (avoid this extra dependency)

#-------------------------------------------------------------------
def test1():

    regrid_geotiff()

#   test1()
#-------------------------------------------------------------------
# def test2():
# 
# 
# #   test2()
#-------------------------------------------------------------------
# def download_data():
# 
#     from pydap.client import open_url
#     from pydap.cas.urs import setup_session
#     dataset_url = 'http://server.example.com/path/to/dataset'
#     session = setup_session(username, password, check_url=dataset_url)
#     dataset = open_url(dataset_url, session=session)
# 
# #    download_data()
#------------------------------------------------------------------- 
def get_raster_bounds( ds, VERBOSE=False):

    #-------------------------------------------------------------
    # Note:  The bounds depend on the map projection and are not
    # necessarily a Geographic bounding box of lons and lats.    
    #-------------------------------------------------------------
    # See:
    # https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html
    # and search on "geotransform".  An example of gdal.SetGeoTransform
    # gives: [xmin, pixel_size, 0, ymax, 0, -pixel_size].
    # Also says args are:
    # [ulx, xDist, rtnX, uly, yDist, rtnY]
    # This is consistent with information below.
    #-------------------------------------------------------------    
    # ulx = upper left x  = xmin
    # uly = upper left y  = ymax
    # lrx = lower right x = xmax
    # lry = lower right y = ymin
    #-----------------------------

    #----------------------------------------------------------  
    # Notice the strange order or parameters here is CORRECT.
    # It is not:  ulx, xres, xskew, uly, yres, yskew
    #----------------------------------------------------------
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    lrx = ulx + (ds.RasterXSize * xres)
    lry = uly + (ds.RasterYSize * yres)

    if (VERBOSE):
        print('ulx, uly   =', ulx, uly)
        print('lrx, lry   =', lrx, lry)
        print('xres, yres = ', xres, yres)
        print('xskew, yskew =', xskew, yskew)
        print('----------------------------------')

    #########################################################
    # Bounding box reported by gdal.info does not match
    # what the GES DISC website is saying.  The result is
    # that gdal.Warp gives all nodata values in output. 
    #########################################################
    return [ulx, lry, lrx, uly]  # [xmin, ymin, xmax, ymax]
 
    #########################################################
    # Bounding box reported by gdal.info does not match
    # what the GES DISC website is saying.  Reversing lats
    # and lons like this doesn't fix the problem.
    #########################################################    
    ## return [lry, ulx, uly, lrx]
    
#   get_raster_bounds()
#------------------------------------------------------------------- 
def bounds_disjoint( bounds1, bounds2, VERBOSE=False):
 
    #-----------------------------------------------------------
    # Note.  Assume both bounds are in same spatial reference
    #        system (SRS), e.g. Geographic lons and lats.
    #------------------------------------------------------------------
    # https://gamedev.stackexchange.com/questions/586/
    # what-is-the-fastest-way-to-work-out-2d-bounding-box-intersection
    #------------------------------------------------------------------    
    b1_xmin = bounds1[0]
    b1_xmax = bounds1[2]
    b2_xmin = bounds2[0]
    b2_xmax = bounds2[2]
#     x_overlap1 = (b1_xmin < b2_xmin) and (b2_xmin < b1_xmax)
#     x_overlap2 = (b2_xmin < b1_xmin) and (b1_xmin < b2_xmax)
#     x_overlap  = (x_overlap1 or x_overlap2)
    
    b1_ymin = bounds1[1]
    b1_ymax = bounds1[3]
    b2_ymin = bounds2[1]
    b2_ymax = bounds2[3]
#     y_overlap1 = (b1_ymin < b2_ymin) and (b2_ymin < b1_ymax) 
#     y_overlap2 = (b2_ymin < b1_ymin) and (b1_ymin < b2_ymax)
#     y_overlap  = (y_overlap1 or y_overlap2)
#     return not(x_overlap and y_overlap)

    disjoint = (b2_xmin > b1_xmax) or (b2_xmax < b1_xmin) or \
               (b2_ymax < b1_ymin) or (b2_ymin > b1_ymax)

    return disjoint
    
#   bounds_disjoint()
#-------------------------------------------------------------------
def get_raster_cellsize( gdal_unit ):

    geotransform = gdal_unit.GetGeoTransform()
    # ulx  = geotransform[0]
    xres = geotransform[1]
    # xrtn = geotransform[2]
    #-----------------------
    # uly  = geotransform[3]
    # yrtn = geotransform[4]  # (not yres !!)
    yres = geotransform[5]  # (not yrtn !!)
    
    return (xres, abs(yres))

#   get_raster_cellsize()
#-------------------------------------------------------------------
def read_geotiff(in_file=None, REPORT=True,
                 MAKE_RTG=False, rtg_file=None, GEO=True ):
 
    # Bounds = [ minlon, minlat, maxlon, maxlat ]

    if (in_file is None):
        in_dir  = '/Users/peckhams/Conflict/Data/GPW-v4/'  
        in_file = 'gpw_v4_population_count_rev11_2020_30_sec.tif'
        in_file = in_dir + in_file

    if (rtg_file is None):
        in_dir   = '/Users/peckhams/Conflict/Data/GPW-v4/'
        rtg_file = 'GPW-v4_global_pop_count.rtg'
        rtg_file = in_dir + rtg_file
          
    #-----------------------------------------    
    # Open the input GeoTIFF file & get info
    #-----------------------------------------
    print('Reading grid from GeoTIFF file...')
    in_unit     = gdal.Open( in_file, gdal.GA_ReadOnly )
    (dx, dy)    = get_raster_cellsize( in_unit )
    in_ncols    = in_unit.RasterXSize
    in_nrows    = in_unit.RasterYSize
    in_bounds   = get_raster_bounds( in_unit )   ######

    if (GEO):
        in_xres_deg = dx
        in_yres_deg = dy
        in_xres_sec = (in_xres_deg * 3600.0)
        in_yres_sec = (in_yres_deg * 3600.0)
        in_xres     = in_xres_sec
        in_yres     = in_yres_sec
    else:
        in_xres = dx
        in_yres = dy

    #------------------------------
    # Get min & max of input grid
    #------------------------------
    band = in_unit.GetRasterBand(1)
    in_nodata = band.GetNoDataValue()
    in_grid   = in_unit.ReadAsArray()
    in_dtype  = in_grid.dtype
    in_gmin   = in_grid.min()
    in_gmax   = in_grid.max()

    #----------------
    # Close in_file
    #----------------
    in_unit = None   # Close in_file

    #-----------------------------
    # Option to save as RTG file
    #-----------------------------
    if (MAKE_RTG):
        rtg_path = rtg_file
        #-------------------------------
        # Option to create an RTI file
        #-------------------------------
        if (GEO):
            pixel_geom = 0
        else:
            pixel_geom = 1
        rti = rti_files.make_info(grid_file=rtg_path,
                  ncols=in_ncols, nrows=in_nrows,
                  xres=in_xres, yres=in_yres,
                  #--------------------------------------
                  data_source='SEDAC',
                  data_type='FLOAT',  ## Must agree with dtype  #######
                  byte_order=rti_files.get_rti_byte_order(),
                  pixel_geom=pixel_geom,
                  zres=0.001, z_units='unknown',
                  y_south_edge=in_bounds[1],
                  y_north_edge=in_bounds[3],
                  x_west_edge=in_bounds[0],
                  x_east_edge=in_bounds[2],
                  box_units='DEGREES')
        rti_files.write_info( rtg_path, rti )    
        rtg = rtg_files.rtg_file() 
        OK  = rtg.open_new_file( rtg_path, rti )
        if not(OK):
            print('ERROR during open_new_file().')
            return       
        rtg.write_grid( in_grid, VERBOSE=True )
        rtg.close_file()
         
    #------------------
    # Optional report
    #------------------
    if (REPORT):
        print('GeoTIFF grid info:')
        print('   ' + in_file )
        print('   ncols  =', in_ncols )
        print('   nrows  =', in_nrows )
        if (GEO):
            print('   xres   =', in_xres, ' [arcsecs]' )
            print('   yres   =', in_yres, ' [arcsecs]' )
        else:
            print('   xres   =', in_xres, ' [meters]')
            print('   yres   =', in_yres, ' [meters]')
        print('   bounds =', in_bounds )
        print('   dtype  =', in_dtype )
        print('   nodata =', in_nodata )
        print('   gmin   =', in_gmin )
        print('   gmax   =', in_gmax ) 
        print('Finished.')
        print()

    #------------------------------
    # Return the grid as an array
    #------------------------------
    return in_grid

#   read_geotiff()
#-------------------------------------------------------------------
def regrid_geotiff(in_file=None, out_file=None, 
                   out_bounds=None,
                   out_xres_sec=None, out_yres_sec=None,
                   RESAMPLE_ALGO='bilinear', REPORT=True):
                   ##### in_nodata=None, out_nodata=None):

    #-----------------------------------   
    # Specify the resampling algorithm
    #-----------------------------------
    algo_dict = {
    'nearest'     : gdal.GRA_NearestNeighbour,
    'bilinear'    : gdal.GRA_Bilinear,
    'cubic'       : gdal.GRA_Cubic,
    'cubicspline' : gdal.GRA_CubicSpline,
    'lanczos'     : gdal.GRA_Lanczos,
    'average'     : gdal.GRA_Average,
    'min'         : gdal.GRA_Min,
    'max'         : gdal.GRA_Max,
    'mode'        : gdal.GRA_Mode,
    'med'         : gdal.GRA_Med }
    
    resample_algo = algo_dict[ RESAMPLE_ALGO ]
    
    #---------------------------------------------------------------
    # Note:  DEM_bounds = [dem_xmin, dem_ymin, dem_xmax, dem_ymax]
    #        Give xres, yres in *arcseconds* for Geographic.
    #        They will then be converted to decimal degrees.
    #        gdal.Warp() clips to a bounding box, and can also
    #        resample to a different resolution.
    #        gdal.Translate() is faster for simple clipping.
    #---------------------------------------------------------------
    if (in_file == None):
        #-----------------------------------------------------------
        # Use Pongo_30sec DEM as a test, which works well.
        # However,  the soil data has same resolution (xres, yres)
        # as the DEM, of 30 arcseconds.  In addition, grid cells
        # outside of South Sudan have NODATA values.
        #-----------------------------------------------------------
        test_dir   = '/Users/peckhams/TF_Tests/Regrid/'
        in_file    = test_dir + 'Baro_Gam_3sec_DEM.tif'
        out_file   = test_dir + 'TEST_11sec_DEM.tif'
        out_xres_sec = 11.0
        out_yres_sec = 11.0
        # Bounds = [ minlon, minlat, maxlon, maxlat ]
        ## out_bounds = [ 34.22125, 7.3704166666, 36.43791666666, 9.50375]
  
    #-----------------------------------------    
    # Open the input GeoTIFF file & get info
    #-----------------------------------------
    in_unit     = gdal.Open( in_file, gdal.GA_ReadOnly )
    (dx, dy)    = get_raster_cellsize( in_unit )
    in_xres_deg = dx
    in_yres_deg = dy
    in_xres_sec = (in_xres_deg * 3600.0)
    in_yres_sec = (in_yres_deg * 3600.0)
    in_ncols    = in_unit.RasterXSize
    in_nrows    = in_unit.RasterYSize
    in_bounds   = get_raster_bounds( in_unit )   ######

    #-----------------------------
    # Double-check the data type
    #--------------------------------------------
    # https://gdal.org/java/org/gdal/gdalconst/
    # gdalconstConstants.html
    # 'complex64'  = 2 32-bit floats
    # 'complex128' = 2 64-bit floats
    #--------------------------------------------
#     band = in_unit.GetRasterBand(1)
#     dtype_code = band.DataType
#     ## print('### dtype_code =', dtype_code)
#     code_to_dtype_map = {
#     0:'unknown', 1:'uint8', 2:'uint16', 3:'int16',
#     4:'uint32', 5:'int32',
#     6:'float32', 7:'float64',
#     8:'complexint32', 9:'complexint64',  # no corresponding dtype
#     10:'complex64',11:'complex128'} 
#     in_dtype2 = code_to_dtype_map[ dtype_code ]

    #-----------------------
    # Get the nodata value
    #-----------------------
    band = in_unit.GetRasterBand(1)
    in_nodata = band.GetNoDataValue()

    #------------------------------
    # Get min & max of input grid
    #------------------------------
    in_grid  = in_unit.ReadAsArray()
    in_dtype = in_grid.dtype
    in_gmin  = in_grid.min()
    in_gmax  = in_grid.max()
    
    #-------------------------------------------------------------
    # If a spatial resolution has not been specified for output,
    # then assume it is the same as the resolution of the input.
    #-------------------------------------------------------------
    if (out_xres_sec is not None):
        out_xres_deg = (out_xres_sec / 3600.0)  # [arcsec] -> [degrees]
    else:
        out_xres_deg = None    # (will default to in_yres_deg)
        ## out_xres_sec = in_xres_sec
        ## out_xres_deg = (out_xres_sec / 3600.0)  # [arcsec] -> [degrees]       
    #----------------------------------------------------------------------
    if (out_yres_sec is not None):
        out_yres_deg = (out_yres_sec / 3600.0)  # [arcsec] -> [degrees]
    else:
        out_yres_deg = None  # (will default to in_yres_deg)
        ## out_yres_sec = in_yres_sec
        ## out_yres_deg = (out_yres_sec / 3600.0)  # [arcsec] -> [degrees]     
                
    #------------------------------------------- 
    # Are out_bounds disjoint from in_bounds ?
    #-------------------------------------------
    # Specifying out_bounds is not required.
    #-------------------------------------------
    if (out_bounds is not None):   
        DISJOINT = bounds_disjoint( in_bounds, out_bounds, VERBOSE=False)
        if (DISJOINT):
            # print('ERROR:  Output bounds do not overlap input bounds.')
            print('ERROR:  Input & output bounding boxes are disjoint.')
            print( '       New grid would contain only nodata.')
            print('  in_file  =', in_file )
            print('  out_file =', out_file )
            in_unit  = None   # Close in_file
            return

    #-------------------------------------------------  
    # Resample & clip and write new grid to GeoTIFF
    #-------------------------------------------------
    # It seems SRTM and MERIT DEMs both use a nodata
    # value of -9999.0, but "in_nodata" is set to 
    # None in the source TIF DEM.  Next line fixes
    # an issue for NE corner of Awash River DEM, and
    # out_nodata becomes same as in_nodata.
    #-------------------------------------------------
    if (in_nodata is None):
        in_nodata = -9999.0
    ## print('##### in_nodata =', in_nodata)
    out_unit = gdal.Warp( out_file, in_unit,
        format = 'GTiff',  # (output format string)
        outputBounds=out_bounds,
        xRes=out_xres_deg, yRes=out_yres_deg,
        srcNodata = in_nodata, #### 2022-05-02
        # dstNodata = out_nodata,    ## equals in_nodata, by default
        resampleAlg = resample_algo )

    #-----------------------
    # Get info on new grid
    #-----------------------
    out_ncols  = out_unit.RasterXSize
    out_nrows  = out_unit.RasterYSize
    out_bounds = get_raster_bounds( out_unit )   ######

    #-------------------------------
    # Get min & max of output grid
    #-------------------------------
    out_band   = out_unit.GetRasterBand(1)
    out_nodata = out_band.GetNoDataValue()
    out_grid   = out_unit.ReadAsArray()
    out_dtype  = out_grid.dtype
    out_gmin   = out_grid.min()
    out_gmax   = out_grid.max()
        
    #----------------------------------        
    # Close both in_file and out_file
    #----------------------------------
    in_unit  = None   # Close in_file
    out_unit = None   # Close out_file

    #------------------
    # Optional report
    #------------------
    if (REPORT):
        print('Input grid file:')
        print('   ' + in_file )
        print('   ncols  =', in_ncols )
        print('   nrows  =', in_nrows )
        print('   xres   =', in_xres_sec, ' [arcsecs]' )
        print('   yres   =', in_yres_sec, ' [arcsecs]' )
        print('   bounds =', in_bounds )
        print('   dtype  =', in_dtype )
        # print('   dtype2 =', in_dtype2 )    # for testing
        print('   nodata =', in_nodata )
        print('   gmin   =', in_gmin )
        print('   gmax   =', in_gmax )
        print()
        print('Output grid file:')
        print('   ' + out_file )
        print('   ncols  =', out_ncols )
        print('   nrows  =', out_nrows )
        print('   xres   =', out_xres_sec, ' [arcsecs]' )
        print('   yres   =', out_yres_sec, ' [arcsecs]' ) 
        print('   bounds =', out_bounds )
        print('   dtype  =', out_dtype )
        print('   nodata =', out_nodata )
        print('   gmin   =', out_gmin )
        print('   gmax   =', out_gmax )      
        print('Finished regridding.')
        print()

    #--------------------------------------------------------  
    # This shows some of the other keywords to gdal.Warp.
    #--------------------------------------------------------      
    # WarpOptions(options=[], format=None, outputBounds=None,
    # outputBoundsSRS=None, xRes=None, yRes=None, targetAlignedPixels=False,
    # width=0, height=0, srcSRS=None, dstSRS=None, srcAlpha=False,
    # dstAlpha=False, warpOptions=None, errorThreshold=None,
    # warpMemoryLimit=None, creationOptions=None, outputType=GDT_Unknown,
    # workingType=GDT_Unknown, resampleAlg=None, srcNodata=None,
    # dstNodata=None, multithread=False, tps=False, rpc=False,
    # geoloc=False, polynomialOrder=None, transformerOptions=None,
    # cutlineDSName=None, cutlineLayer=None, cutlineWhere=None,
    # cutlineSQL=None, cutlineBlend=None, cropToCutline=False,
    # copyMetadata=True, metadataConflictValue=None,
    # setColorInterpretation=False, callback=None, callback_data=None)
    # 
    # Create a WarpOptions() object that can be passed to gdal.Warp()
    # Keyword arguments are : options --- can be be an array of strings,
    # a string or let empty and filled from other keywords.

#    regrid_geotiff()
#-------------------------------------------------------------------
def clip_geotiff(in_file=None, out_file=None, 
                 out_bounds=None, REPORT=True):

    #################################################
    # Note.  This is incomplete and not tested yet.
    #        regrid_geotiff() clips & resamples.
    #################################################
    if (in_file == None):
        #-----------------------------------------------------------
        # Use Pongo_30sec DEM as a test, which works well.
        # However,  the soil data has same resolution (xres, yres)
        # as the DEM, of 30 arcseconds.  In addition, grid cells
        # outside of South Sudan have NODATA values.
        #-----------------------------------------------------------
        test_dir   = '/Users/peckhams/TF_Tests/Regrid/'
        in_file    = test_dir + 'Baro_Gam_3sec_DEM.tif'
        out_file   = test_dir + 'TEST_3sec_DEM.tif'
        ## out_bounds = [24.079583333333,  6.565416666666, 27.379583333333, 10.132083333333 ]

    #-----------------------------------------    
    # Open the input GeoTIFF file & get info
    #-----------------------------------------
    in_unit     = gdal.Open( in_file, gdal.GA_ReadOnly )
    (dx, dy)    = get_raster_cellsize( in_unit )
    in_xres_deg = dx    # [degrees]
    in_yres_deg = dy    # [degrees]
    in_xres_sec = (in_xres_deg * 3600.0)
    in_yres_sec = (in_yres_deg * 3600.0)
    in_ncols    = in_unit.RasterXSize
    in_nrows    = in_unit.RasterYSize
    in_bounds   = get_raster_bounds( in_unit )   ######

    #------------------------------------------- 
    # Are out_bounds disjoint from in_bounds ?
    #-------------------------------------------
    # Specifying out_bounds is not required.
    #-------------------------------------------
    if (out_bounds is not None):   
        DISJOINT = bounds_disjoint( in_bounds, out_bounds, VERBOSE=False)
        if (DISJOINT):
            print('ERROR:  Output bounds do not overlap input bounds.')
            return

    #------------------------------------------------  
    # Resample & clip and write new grid to GeoTIFF
    #------------------------------------------------
#     out_unit = gdal.Warp( out_file, in_unit,
#         format = 'GTiff',  # (output format string)
#         outputBounds=out_bounds,
#         xRes=out_xres_deg, yRes=out_yres_deg,
#         resampleAlg = resample_algo )

    #----------------------------------------------------- 
    # Use gdal.Translate to clip to output bounding box.
    # See: https://gdal.org/programs/gdal_translate.html
    #-----------------------------------------------------
    # Convert "bounds array" to "projWin array".
    #    bounds  = [ulx, lry, lrx, uly] 
    #              [xmin, ymin, xmax, ymax]
    #    projWin = [ulx, uly, lrx, lry].
    #-----------------------------------------------------
    pWin = [out_bounds[0], out_bounds[3], out_bounds[2], out_bounds[1]]     
    out_unit = gdal.Translate(out_file, in_unit,
                    format = 'GTiff',  # (output format string)
                    projWin = pWin)
    
    #-----------------------
    # Get info on new grid
    #-----------------------
    out_ncols    = out_unit.RasterXSize
    out_nrows    = out_unit.RasterYSize
    (dx2, dy2)   = get_raster_cellsize( out_unit )
    out_xres_deg = dx2   # [degrees]
    out_yres_deg = dy2   # [degrees]
    out_bounds   = get_raster_bounds( out_unit )   ######
            
    #----------------------------------        
    # Close both in_file and out_file
    #----------------------------------
    in_unit  = None   # Close in_file
    out_unit = None   # Close out_file

    #------------------
    # Optional report
    #------------------
    if (REPORT):
        print('Input grid file:')
        print( in_file )
        print('   ncols  =', in_ncols )
        print('   nrows  =', in_nrows )
        print('   xres   =', in_xres_sec, ' [arcsecs]' )
        print('   yres   =', in_yres_sec, ' [arcsecs]' )
        print('   bounds =', in_bounds )
        print('Output grid file:')
        print( out_file )
        print('   ncols  =', out_ncols )
        print('   nrows  =', out_nrows )
        print('   xres   =', out_xres_sec, ' [arcsecs]' )
        print('   yres   =', out_yres_sec, ' [arcsecs]' ) 
        print('   bounds =', out_bounds )       
        print('Finished regridding.')
        print()      
           
#   clip_geotiff()
#-------------------------------------------------------------------  
def save_grid_to_geotiff( new_tif_file, grid,
                          ulx, uly, xres, yres,
                          dtype='float32', nodata=None):
                              
    #---------------------------------------------------------
    # Notes:
    #---------------------------------------------------------
#     new_nodata = -9999.0
#     grid[ grid <= nodata ] = new_nodata

    #--------------------------------------
    # These are the geotransform elements
    #--------------------------------------
    # ulx  = geotransform[0]
    # xres = geotransform[1]
    # xrtn = geotransform[2]  (rtn = rotation, if any)
    # uly  = geotransform[3]
    # yrtn = geotransform[4]  # (not yres !!)
    # yres = geotransform[5]  # (not yrtn !!) 
     
    if (dtype == 'float32'):
        data_type = gdal.GDT_Float32
    elif (dtype == 'int32'):
        data_type = gdal.GDT_Int32
    elif (dtype == 'int16'):
        data_type = gdal.GDT_Int16
    elif (dtype == 'float64'):
        data_type = gdal.GDT_Float64
    else:
        data_type = gdal.GDT_Float32
                        
    ncols  = grid.shape[1]
    nrows  = grid.shape[0]
    nbands = 1
    xrtn   = 0.0
    yrtn   = 0.0
    

    #-----------------------------------------
    # Use same EPSG code as source data that
    # was downloaded from GES DISC DAAC.
    #-----------------------------------------
    epsg_code = 4326  ## (WGS 84,  DOUBLE CHECK #####################)

    
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(new_tif_file, ncols, nrows, nbands, data_type)
    outRaster.SetGeoTransform((ulx, xres, xrtn, uly, yrtn, yres))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray( grid )
    outband.SetNoDataValue( nodata )
    outRasterSRS = osr.SpatialReference()
    #------------------------------------------
    # Need to specify Projection from scratch
    #------------------------------------------
    outRasterSRS.ImportFromEPSG( epsg_code )
    #------------------------------------------
    # Use Projection from another GDAL raster
    #------------------------------------------
    ###### outRasterSRS.ImportFromWkt(inRaster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

#   save_grid_to_geotiff()
#-------------------------------------------------------------------
def regrid_rts_file( rts_file_in, rts_file_out,
        xres_in_sec,  yres_in_sec, nx_in, ny_in,
        xres_out_sec, yres_out_sec, resample_algo='bilinear',
        BARO=True, LOL_KURU=False):

    #---------------------------------------
    # Note:  This has not been tested yet.
    #---------------------------------------
    
    #------------------------------------------
    # Use bounding boxes for Baro or Lol-Kuru
    #------------------------------------------
    if (BARO):
        #------------------------------------------
        # For the 30 arcsecond spatial resolution
        # (ncols, nrows) = (267, 257)
        #------------------------------------------
        minlat = 7.362083333332    # (south)
        maxlat = 9.503749999999    # (north)
        minlon = 34.221249999999   # (west)
        maxlon = 36.446249999999   # (east)
        
    if (LOL_KURU):
        #------------------------------------------
        # For the 30 arcsecond spatial resolution
        # (ncols, nrows) = (483, 364)
        #------------------------------------------
        minlat = 6.532916666667    # (south)
        maxlat = 9.566250000000    # (north)
        minlon = 23.995416666666   # (west)
        maxlon = 28.020416666666   # (east)

    #------------------------------------------------- 
    # Set lon & lat of upper-left corner (northwest)
    # ulx = upper left x = min longitude  (west)
    # uly = upper left y = max latitude   (north)
    #-------------------------------------------------
    ulx = minlon
    uly = maxlat
    bounds_in = [minlon, minlat, maxlon, maxlat]
    
    #--------------------------  
    # Open the input RTS file
    #-------------------------- 
    rts_unit_in = open( rts_file_in, 'rb')

    #---------------------------  
    # Open the output RTS file
    #--------------------------- 
    rts_unit_out = open( rts_file_out, 'wb')

    #-----------------------------------------
    # Prepare to read grids from RTS file in
    #-----------------------------------------
    nodata    = None
    dtype     = 'float32'    
    grid_in   = np.zeros( (nx_in, ny_in), dtype=dtype )
    n_values  = nx_in * ny_in
    grid_size = n_values * 4   # (float32)
    file_size = os.path.getsize( rts_file_in )
    n_grids   = (file_size / grid_size)
    tmp_file  = 'TEMP.tif'
    tmp_file2 = 'TEMP2.tif'
    count     = 0

    for k in range(n_grids):
        #----------------------------------------------
        # Read one grid (or "frame") from rts_file_in
        #----------------------------------------------
        grid = np.fromfile( rts_unit_in, count=n_values,
                            dtype=dtype )      
        grid = grid.reshape( ny_in, nx_in )

        #---------------------------------------------      
        # Save grid (an ndarray) to GeoTIFF tmp_file
        ####### xres_in_sec, or xres_deg ????
        ########################################
        #---------------------------------------------
        save_grid_to_geotiff( tmp_file, grid,
                              ulx, uly, xres_in_deg, yres_in_deg)

        #------------------------------------------
        # Resample tmp_file to xres_out, yres_out
        # and save to tmp_file2
        #------------------------------------------
#         out_xres_sec = (out_xres_deg * 3600.0)
#         out_yres_sec = (out_yres_deg * 3600.0)
#         regrid_geotiff(in_file=tmp_file, out_file=tmp_file2, 
#                    out_bounds=None,
#                    out_xres_sec=out_xres_sec, out_yres_sec=out_yres_sec,
#                    #### in_nodata, out_nodata, 
#                    RESAMPLE_ALGO='bilinear', REPORT=True):
                  
        #---------------------------------------------
        # Read ndarray from GeoTIFF into GDAL object
        #---------------------------------------------
        ds_in = gdal.Open( tmp_file, gdal.GA_ReadOnly )
        
        #----------------------------------------- 
        # Resample to grid to xres_out, yres_out
        #-----------------------------------------
        xres_out_deg = (xres_out_sec / 3600.0)
        yres_out_deg = (yres_out_sec / 3600.0)       
        grid2 = gdal_regrid_to_dem_grid( ds_in, tmp_file2,
                    bounds_out, xres_out_deg, yres_out_deg,
                    nodata=nodata, VERBOSE=False,
                    RESAMPLE_ALGO=resample_algo)
        #------------------------
        # Close the tmp tif_file
        #------------------------
        ds_in = None
###
        
        #---------------------------------------        
        # Write resampled grid to new RTS file
        #---------------------------------------
        grid2 = np.float32( grid2 )
        grid2.tofile( rts_unit_out )
        count += 1
        
    #------------------        
    # Close all files
    #------------------
    rts_unit_in.close()
    rts_unit_out.close()
 
    #------------------------------------  
    # Check the size of output RTS file
    #------------------------------------
    nx_out = grid2.shape[1]
    ny_out = grid2.shape[2]
    filesize = os.path.getsize( rts_file_out)
    expected = (nx_out * ny_out * 4 * count)
    if (filesize != expected):
        print('ERROR:')
        print('Output RTS filesize =', filesize)
        print('But expected size   =', expected)
        print()
          
    print('Finished regridding RTS file.')
    print('New RTS file is: ' + rts_file_out )
    print('Number of grids = ' + str(count) )
    if (n_grids != count):
        print('ERROR:')
        print('Input  RTS file has n_grids =', n_grids )
        print('Output RTS file has n_grids =', count )
    print()
                                     
#   regrid_rts_file()
#-------------------------------------------------------------------
#-------------------------------------------------------------------
def read_nc_grid( nc_file=None, var_name='precipitationCal',
                  REPORT=False):

    if (nc_file == None):
        nc_file = 'TEST.nc4'

    ds = gdal.Open("NETCDF:{0}:{1}".format(nc_file, layer_name))
    grid = ds.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    ds = None  # (close ds)
        
    if (REPORT):
        print( 'grid.min() =', grid.min() )
        print( 'grid.max() =', grid.max() )
        print ('grid.shape =', grid.shape )

    return grid
            
#   read_nc_grid()
#-------------------------------------------------------------------
# def read_nc_as_array( nc_file=None, var_name='precipitationCal',
#                       REPORT=False):
#                       
#     ds = gdal.Open( nc_file )
#     if (ds is None):
#         print( 'Open failed.')
#          sys.exit()
#     
#     if (ds.GetSubDatasets() >= 1):
#         subdataset = 'NETCDF:"' + nc_file + '":' + var_name
#         ds_sd = gdal.Open( subdataset )
#         NDV   = ds_sd.GetRasterBand(1).GetNoDataValue()
#         ncols = ds_sd.RasterXsize
#         nrows = ds_sd.RasterYsize
#         GeoT  = ds_sd.GetGeoTransform()
#         ds    = None
#         ds_sd = None  
#   
# #   read_nc_as_array()
#------------------------------------------------------------------- 
def gdal_open_nc_file( nc_file, var_name, VERBOSE=False):

    if (VERBOSE):  print('Opening netCDF file...')
    ### ds_in = gdal.Open("NETCDF:{0}:{1}".format(nc_file, var_name), gdal.GA_ReadOnly )
    ds_in  = gdal.Open("NETCDF:{0}:{1}".format(nc_file, var_name) )
    band   = ds_in.GetRasterBand(1)
    nodata = band.GetNoDataValue()

    if (VERBOSE):  print('Reading grid of values...')
    g1 = band.ReadAsArray()
    ## g1 = ds_in.ReadAsArray(0, 0, ds_in.RasterXSize, ds_in.RasterYSize)

    if (VERBOSE):
        print( 'grid: min =', g1.min(), ', max =', g1.max() )
        print( 'grid.shape =', g1.shape )
        print( 'grid.dtype =', g1.dtype )
        print( 'grid nodata =', nodata )
        w  = np.where(g1 == nodata)
        nw = w[0].size
        print( 'nodata count =', nw)
        print()

    return (ds_in, g1, nodata)

#   gdal_open_nc_file()
#------------------------------------------------------------------- 
def fix_gpm_raster_bounds( ds, VERBOSE=False):

    #------------------------------------------------------------
    # Note:  NetCDF files downloaded from the GES DISC website
    #        have corner coordinate lons and lats reversed.
    #        I checked with multiple files for which bounding
    #        box was known that when gdalinfo reports Corner
    #        Coordinates, it uses (lon, lat) vs. (lat, lon).
    #        Here, we use SetGeoTransform to fix the bounding
    #        box, so that gdal.Warp() and other gdal functions
    #        will work correctly.  (8/14/2019)
    #------------------------------------------------------------
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    lrx = ulx + (ds.RasterXSize * xres)
    lry = uly + (ds.RasterYSize * yres)
 
    ulx2  = lry
    uly2  = lrx
    lrx2  = uly
    lry2  = ulx
    #lrx2 = ulx2 + (ds.RasterXsize * xres)
    # Note:  (xres > 0, yres < 0)
  
    if (VERBOSE):
        #----------------------------------------------------- 
        # These print out correctly, but the reported corner
        # coordinates are now really messed up.
        # Need to close or flush to make new info "stick" ?
        #----------------------------------------------------- 
        print('in_bounds  =', ulx, lry, lrx, uly)      # (2,20,15,40)
        print('out_bounds =', ulx2, lry2, lrx2, uly2 ) # (20,2,40,15)
        print(' ')
     
    ds.SetGeoTransform( (ulx2, xskew, xres, uly2, yskew, yres) )
    
#   fix_gpm_raster_bounds()
#------------------------------------------------------------------- 
def fix_gpm_file_as_geotiff( nc_file, var_name, out_file, 
                             out_nodata=None, VERBOSE=False):

    if (VERBOSE):  print('Reading data from netCDF file...')
    ### raster = gdal.Open("NETCDF:{0}:{1}".format(nc_file, var_name), gdal.GA_ReadOnly )
    raster  = gdal.Open("NETCDF:{0}:{1}".format(nc_file, var_name) )
    band   = raster.GetRasterBand(1)
    ncols  = raster.RasterXSize
    nrows  = raster.RasterYSize
    proj   = raster.GetProjectionRef()
    bounds = get_raster_bounds( raster )   ######
    nodata = band.GetNoDataValue()
    array  = band.ReadAsArray() 
    ## array = raster.ReadAsArray(0, 0, ds_in.RasterXSize, ds_in.RasterYSize)
    #----------------------------------------------
    # Get geotransform for array in nc_file
    # Note:  These look strange, but are CORRECT.
    #----------------------------------------------
    geotransform = raster.GetGeoTransform()
    ulx  = geotransform[0]
    xres = geotransform[1]
    xrtn = geotransform[2]
    #-----------------------
    uly  = geotransform[3]
    yrtn = geotransform[4]  # (not yres !!)
    yres = geotransform[5]  # (not yrtn !!)
    raster = None    # Close the nc_file

    if (VERBOSE):
        print( 'grid: min   =', array.min(), ', max =', array.max() )
        print( 'grid.shape  =', array.shape )
        print( 'grid.dtype  =', array.dtype )
        print( 'grid nodata =', nodata )
        w  = np.where(array == nodata)
        nw = w[0].size
        print( 'nodata count =', nw)
        print( ' ' )

    #----------------------------------------------    
    # Rotate the array; column major to row major
    # a           = [[7,4,1],[8,5,2],[9,6,3]]
    # np.rot90(a) = [[1,2,3],[4,5,6],[7,8,9]]
    #----------------------------------------------
    if (VERBOSE):  print('Rotating grid (90 deg ccw)...')
    ### array2 = np.transpose( array )
    array2 = np.rot90( array )    ### counter clockwise
    ncols2 = nrows
    nrows2 = ncols

    #-----------------------------       
    # Change the nodata value ?
    # This may lead to problems.
    #-----------------------------
    if (out_nodata is not None):
        array2[ array2 <= nodata ] = out_nodata
    
    #-----------------------------------------    
    # Build new geotransform & projectionRef
    #-----------------------------------------
    if (VERBOSE):  print('Adjusting metadata...')
    lrx    = bounds[2]
    lry    = bounds[1]
    ulx2   = lry
    uly2   = lrx
    xres2  = -yres
    yres2  = -xres
    xrtn2  = yrtn
    yrtn2  = xrtn
    geotransform2 = (ulx2, xres2, xrtn2, uly2, yrtn2, yres2)
    proj2 = proj

    if (VERBOSE):
        print( 'geotransform  =', geotransform )
        print( 'geotransform2 =', geotransform2 )
    
    #------------------------------------
    # Write new array to a GeoTIFF file
    #------------------------------------
    if (VERBOSE):  print('Writing grid to GeoTIFF file...')
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(out_file, ncols2, nrows2, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform( geotransform2 )
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray( array2 )
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt( proj2 )
    outRaster.SetProjection( outRasterSRS.ExportToWkt() )
    outband.FlushCache()

    #---------------------    
    # Close the out_file
    #---------------------
    outRaster = None
    if (VERBOSE):
        print('Finished.')
        print()
    
#   fix_gpm_file_as_geotiff()
#-------------------------------------------------------------------
def gdal_regrid_to_dem_grid( ds_in, tmp_file, 
         DEM_bounds, DEM_xres_deg, DEM_yres_deg,
         nodata=None, RESAMPLE_ALGO='bilinear',
         VERBOSE=False, IN_MEMORY=False):

    #-----------------------------------   
    # Specify the resampling algorithm
    #-----------------------------------
    algo_dict = {
    'nearest'     : gdal.GRA_NearestNeighbour,
    'bilinear'    : gdal.GRA_Bilinear,
    'cubic'       : gdal.GRA_Cubic,
    'cubicspline' : gdal.GRA_CubicSpline,
    'lanczos'     : gdal.GRA_Lanczos,
    'average'     : gdal.GRA_Average,
    'min'         : gdal.GRA_Min,
    'max'         : gdal.GRA_Max,
    'mode'        : gdal.GRA_Mode,
    'med'         : gdal.GRA_Med }
    
    resample_algo = algo_dict[ RESAMPLE_ALGO ]

    #--------------------------------------------------
    # Use gdal.Warp to clip and resample to DEM grid
    # then save results to a GeoTIFF file (tmp_file).
    #--------------------------------------------------
    ds_tmp = gdal.Warp( tmp_file, ds_in,
        format = 'GTiff',  # (output format string)
        outputBounds=DEM_bounds,
        xRes=DEM_xres_deg, yRes=DEM_yres_deg,
        srcNodata=nodata,      ########
        ### dstNodata=nodata,  ########
        resampleAlg = resample_algo )

    grid = ds_tmp.ReadAsArray()
    
    ds_tmp = None   # Close tmp_file
    if (IN_MEMORY):
        gdal.Unlink( tmp_file )  ######## DOUBLE-CHECK THIS
                   
    return grid

#   gdal_regrid_to_dem_grid()       
#-------------------------------------------------------------------    
def create_rts_from_nc_files( rts_file='TEST.rts', NC4=False,
           var_name=None, GPM=True, resample_algo='bilinear',
           BARO_60=False, IN_MEMORY=False, VERBOSE=False,
           DEM_bounds=None, DEM_xres_sec=None, DEM_yres_sec=None,
           DEM_ncols=None, DEM_nrows=None, SILENT=False):

    VERBOSE = True
    #------------------------------------
    # Note: GPM and GLDAS are currently
    #       the only supported products
    #------------------------------------
    GLDAS = not(GPM)  #############
    
    #------------------------------------------------------
    # Note: See function above for resampling algorithms.
    #------------------------------------------------------
    # For info on GDAL constants, see:
    # https://gdal.org/python/osgeo.gdalconst-module.html
    #------------------------------------------------------
    if (BARO_60):
        DEM_bounds = [ 34.221249999999, 7.362083333332, 36.450416666666, 9.503749999999]
        ### DEM_bounds = [ 7.362083333332, 34.221249999999, 9.503749999999, 36.450416666666]
        DEM_xres_sec = 60.0   # (60 arcsecs = 60/3600 degrees)
        DEM_yres_sec = 60.0
        DEM_ncols  = 134
        DEM_nrows  = 129
        
    #--------------------------------------------------------    
    # See: Pongo_30sec DEM: MINT_2019/__Regions/South_Sudan
    #      Adjusted this later as the Lol-Kuru.
    #--------------------------------------------------------     
#     if (PONGO_30):
#         DEM_bounds = [24.079583333333,  6.565416666666, 27.379583333333, 10.132083333333]
#         ## DEM_bounds = [6.565416666666, 24.079583333333,  10.132083333333, 27.379583333333]
#         DEM_xres   = 1./120   # (30 arcsecs = 30/3600 degrees)
#         DEM_yres   = 1./120   # (30 arcsecs = 30/3600 degrees)
#         DEM_ncols  = 396
#         DEM_nrows  = 428

    if (DEM_bounds is None) or (DEM_xres_sec is None) or \
       (DEM_yres_sec is None) or (DEM_ncols is None) or (DEM_nrows is None):
        print('ERROR:  All DEM bounding box information is required.')
        return

    #---------------------------------------------
    # Convert DEM_xres_sec to DEM_xres_deg, etc.
    #---------------------------------------------
    DEM_xres_deg = (DEM_xres_sec / 3600.0)
    DEM_yres_deg = (DEM_yres_sec / 3600.0)
           
    #----------------------------------------
    # Set the name of the variable (precip)
    #----------------------------------------
    # IMERG Technical Documentation, Table 1,
    # p. 25 (Huffman et al. (2019)).
    if (GPM):
        #------------------------------------------
        # IMERG Technical Documentation, Table 1,
        # p. 25 (Huffman et al. (2019)).
        # Should be 'precipitationCal' and NOT
        #   HQprecipitation (microwave only).
        #------------------------------------------
        if (var_name is None):
            var_name = 'precipitationCal'    # [mmph]
            ## var_name = 'HQprecipitation'
    #------------------------------------------------  
    if (GLDAS):
        if (var_name is None):
            var_name = 'Rainf_f_tavg' # Total precip. rate [kg m-2 s-1] 
            # var_name = 'Rainf_tavg' # Rain  precip. rate [kg m-2 s-1] 

    #------------------------------------------------     
    if (var_name is None):
        print('ERROR: var_name is required.')

    #-----------------------------------------    
    # Use a temp file in memory or on disk ?
    #-----------------------------------------
    if (IN_MEMORY):
        tmp_file = '/vsimem/TEMP.tif'
    else:
        tmp_file = 'TEMP.tif'
  
    #-------------------------    
    # Open RTS file to write
    #-------------------------
    rts_unit = open( rts_file, 'wb' )

    #------------------------------------------------
    # Get list of all nc files in working directory
    #------------------------------------------------
    if (NC4):
        nc_file_list = glob.glob( '*.nc4' )
    else:
        nc_file_list = glob.glob( '*.nc' )

    #-------------------------------------       
    # IMPORTANT!  Sort the list of files
    # Bug fix: 08-04-2020.
    #-------------------------------------
    nc_file_list = sorted( nc_file_list)
    
    #-------------------------------------------------
    # It may not be wise to change the nodata values
    # because GDAL may not recognize them as nodata
    # during spatial interpolation, etc.
    # This could be the cause of "artifacts".
    #-------------------------------------------------
    # A nodata value of zero may be okay for rainfall
    # rates, but is not good in general.
    #-------------------------------------------------    
    #### rts_nodata = -9999.0    ######
    # rts_nodata = 0.0
    rts_nodata = None   # (do not change nodata values)
   
    count = 0
    vmax  = -1e8
    vmin  = 1e8
    bad_count = 0
    BAD_FILE  = False
    tif_file  = 'TEMP1.tif'
    
    print('Working...') 
    for nc_file in nc_file_list:
        if (GPM):
            #------------------------------------------
            # Option to fix problem with bounding box
            #------------------------------------------
            ### fix_gpm_raster_bounds( ds_in )

            #------------------------------------------
            # Fix GPM netCDF file, resave as GeoTIFF, 
            # then open the new GeoTIFF file
            #------------------------------------------
            if not(SILENT):  print('resaving as GeoTIFF...')
            fix_gpm_file_as_geotiff( nc_file, var_name, tif_file,
                                     out_nodata=rts_nodata )
            ds_in = gdal.Open( tif_file )
            grid1 = ds_in.ReadAsArray()
        elif (GLDAS):
            #-------------------------------
            # Open the original netCDF file
            #--------------------------------
            (ds_in, grid1, nodata) = gdal_open_nc_file( nc_file, var_name, VERBOSE=VERBOSE)

        #------------------------------------  
        # Note: This gives min & max before
        #       clipping to DEM bounds
        #------------------------------------
#         gmax  = grid1.max()
#         gmin  = grid1.min()
#         vmax  = max( vmax, gmax )
#         vmin  = min( vmin, gmin )
        band  = ds_in.GetRasterBand(1)
        nc_nodata = band.GetNoDataValue()  # same as nodata just above

        if (VERBOSE):
            print( '===============================================================')
            print( 'count =', (count + 1) )
            print( '===============================================================')
            print( 'grid1: min   =', grid1.min(), ', max =', grid1.max() )
            print( 'grid1.shape  =', grid1.shape )
            print( 'grid1.dtype  =', grid1.dtype )
            print( 'grid1 nodata =', nc_nodata, '(in original nc_file)' )
            ### w  = np.where(grid1 > nc_nodata)
            w  = np.where(grid1 != nc_nodata)
            nw = w[0].size
            print( 'grid1 # data =', nw)
            print( ' ' )
              
        #--------------------------------------        
        # Use gdal.Info() to print/check info
        #--------------------------------------
        ## print( gdal.Info( ds_in ) )
        ## print( '===============================================================')

        #-----------------------------------------------        
        # Check if the bounding boxes actually overlap
        #-----------------------------------------------
        if not(SILENT):  print('Comparing bounds...')
        ds_bounds = get_raster_bounds( ds_in, VERBOSE=False )
        if (bounds_disjoint( ds_bounds, DEM_bounds )):
            print( '###############################################')
            print( 'WARNING: Bounding boxes do not overlap.')
            print( '         New grid will contain only nodata.')
            print( '###############################################')
            print( 'count =', count )
            print( 'file  =', nc_file )
            print( 'ds_bounds  =', ds_bounds )
            print( 'DEM_bounds =', DEM_bounds )
            print( ' ')
            bad_count += 1
            BAD_FILE = True

        #-------------------------------------------
        # Clip and resample data to the DEM's grid
        # then save to a temporary GeoTIFF file.
        #-------------------------------------------
        if not(BAD_FILE):
            print('Regridding to DEM grid...')
            grid2 = gdal_regrid_to_dem_grid( ds_in, tmp_file,
                         DEM_bounds, DEM_xres_deg, DEM_yres_deg,
                         nodata=rts_nodata, IN_MEMORY=IN_MEMORY,
                         RESAMPLE_ALGO=resample_algo )

            # This was already done in gdal_regrid_to_dem_grid()
            ## ds_in = None   # Close the tmp_file
            ## if (IN_MEMORY):
            ##    gdal.Unlink( tmp_file )
            
            if (VERBOSE):
                print( 'grid2: min  =', grid2.min(), ', max =', grid2.max() )
                print( 'grid2.shape =', grid2.shape )
                print( 'grid2.dtype =', grid2.dtype )
                if (rts_nodata is not None):
                    ## w  = np.where(grid2 > rts_nodata)
                    w  = np.where(grid2 != rts_nodata)
                    nw = w[0].size
                    print( 'grid2 # data =', nw)
                    print()
        else:
            grid2 = np.zeros( (DEM_nrows, DEM_ncols), dtype='float32' )
            if (rts_nodata is not None):
                grid2 += rts_nodata

        #----------------------------------------------- 
        # Perform unit conversions for some GLDAS vars
        # Usually, nc_nodata = -9999.0.
        #-----------------------------------------------
        if (GLDAS):
            rain_vars = ['Rainf_f_tavg', 'Rainf_tavg']
            if (var_name in rain_vars):
                print('Converting units: [kg m-2 s-1] to [mmph]...')
                w = (grid2 != nc_nodata)  # (boolean array)
                grid2[w] *= 3600.0        # (preserve nodata)
            #--------------------------------------------------------------                
#             if (var_name == 'Snowf_tavg'):
#                 print('Converting units: [kg m-2 s-1] to [mmph]...')
#                 w = (grid2 != nc_nodata)  # (boolean array)
#                 grid2[w] *= 3600.0        # (preserve nodata)            
            #--------------------------------------------------------------
            temp_vars = ['Tair_f_inst', 'AvgSurfT_tavg', 'AvgSurfT_inst', 'SnowT_tavg']
            if (var_name in temp_vars):
                print('Converting units: Kelvin to Celsius...')
                w = (grid2 != nc_nodata)  # (boolean array)
                grid2[w] -= 273.15        # (preserve nodata)
            #--------------------------------------------------------------             
            if (var_name == 'Albedo_inst'):
                print('Converting units: % to [0,1]...')
                w = (grid2 != nc_nodata)  # (boolean array)
                grid2[w] /= 100.0         # (preserve nodata)
            #--------------------------------------------------------------             
            if (var_name == 'Psurf_f_inst'):
                print('Converting units: Pa to mbar...')
                w = (grid2 != nc_nodata)  # (boolean array)
                grid2[w] /= 100.0         # (preserve nodata)
                                
        #-------------------------  
        # Write grid to RTS file
        #-------------------------
        grid2 = np.float32( grid2 )
        vmax  = max( vmax, grid2.max() )
        vmin  = min( vmin, grid2.min() )
        grid2.tofile( rts_unit )
        count += 1
        if not(SILENT):  print('count =', count)

    #---------------------
    # Close the RTS file
    #---------------------
    rts_unit.close()

    #---------------------------------------
    # Create an RTI file for this RTS file
    #---------------------------------------
    rti = rti_files.make_info(grid_file=rts_file,
              ncols=DEM_ncols, nrows=DEM_nrows,
              xres=DEM_xres_sec, yres=DEM_yres_sec,
              #--------------------------------------
              data_source='TopoFlow 3.6 Regrid',
              data_type='FLOAT',
              byte_order=rti_files.get_rti_byte_order(),
              pixel_geom=0,
              zres=0.001, z_units='unknown',
              y_south_edge=DEM_bounds[1],
              y_north_edge=DEM_bounds[3],
              x_west_edge=DEM_bounds[0],
              x_east_edge=DEM_bounds[2],
              box_units='DEGREES')
    rti_files.write_info( rts_file, rti )

    #---------------------------
    # Determine variable units
    #---------------------------
    if (GPM):
        units = 'mmph'
    else:
        #--------------------------------------------------
        # Note: The GLDAS Noah-LSM Parameters (Table 3.3) 
        #       Other GLDAS products have other vars
        #         and some of those are included here.
        #       "_tavg" vars are backward 3-hour avg.
        #       "_acc" vars are backward 3-hour accum.
        #       "_inst" vars are instantaneous
        #       "_f" vars are forcing vars.
        #       Updated this map on 2023-08-28.
        #------------------------------------------------  
        umap = {
        'Swnet_tavg'         : 'W m-2',       # Net shortwave rad. flux
        'Lwnet_tavg'         : 'W m-2',       # Net longwave rad. flux
        'Qle_tavg'           : 'W m-2',       # Latent heat net flux
        'Qh_tavg'            : 'W m-2',       # Sensible heat net flux
        'Qg_tavg'            : 'W m-2',       # Ground heat flux
        'SWdown_f_tavg'      : 'W m-2',       # Downward SW rad. flux
        'LWdown_f_tavg'      : 'W m-2',       # Downward LW rad. flux
        'PotEvap_tavg'       : 'W m-2',       # Potential evap. rate
        'ECanop_tavg'        : 'W m-2',       # Canopy water evaporation
        'Tveg_tavg'          : 'W m-2',       # Transpiration
        'ESoil_tavg'         : 'W m-2',       # Direct evap. from bare soil
        #-----------------------------------------------------------------------
#       'Rainf_tavg'         : 'kg m-2 s-1',  # Rain precip. rate (see below)
#       'Rainf_f_tavg'       : 'kg m-2 s-1',  # Total precip. rate (see below)
        'Snowf_tavg'         : 'kg m-2 s-1',  # Snow precip. rate
        'Evap_tavg'          : 'kg m-2 s-1',  # Evapotranspiration
        'EvapSnow_tavg'      : 'kg m-2 s-1',  # Snow evaporation
        'Qs_tavg'            : 'kg m-2 s-1',  # Storm surface runoff
        'Qsb_tavg'           : 'kg m-2 s-1',  # Baseflow-gw runoff
        'Qsm_tavg'           : 'kg m-2 s-1',  # Snow melt
        #-----------------------------------------------------------------------
        'Qs_acc'             : 'kg m-2',      # Storm surface runoff
        'Qsb_acc'            : 'kg m-2',      # Baseflow-gw runoff
        'Qsm_acc'            : 'kg m-2',      # Snow melt
        'RootMoist_inst'     : 'kg m-2',      # Root zone soil moisture
        'CanopInt_inst'      : 'kg m-2',      # Plant canopy surface water
        #-----------------------------------------------------------------------
        'Rainf_tavg'         : 'mm h-1',      # Rain precip. rate (converted)
        'Rainf_f_tavg'       : 'mm h-1',      # Total precip. rate (converted)        
        #-----------------------------------------------------------------------
        'Tair_f_inst'        : 'C',           # Air temp. (K to C)
        'SnowT_tavg'         : 'C',           # Snow surf. temp. (K to C)
        'AvgSurfT_inst'      : 'C',           # Avg. surface skin temp. (K to C)
        'AvgSurfT_tavg'      : 'C',           # Avg. surface skin temp. (K to C)
        #-----------------------------------------------------------------------
        'SnowDepth_inst'     : 'm',           # Snow depth (inst.)
        'SnowDepth_tavg'     : 'm',           # Snow depth (time avg.)
        'Albedo_inst'        : 'none [0,1]',  # Albedo  (% to [0,1])
        'Psurf_f_inst'       : 'Pa',          # Surface pressure 
#       'Psurf_f_tavg'       : 'Pa',          # Surface pressure  
        'Qair_f_inst'        : 'kg kg-1',     # Specific humidity     
        #-----------------------------------------------------------------------
        'SWE_inst'              : 'kg m-2',   # Snow depth water equiv.
        'SWE_tavg'              : 'kg m-2',   # Snow depth water equiv.
        'Qsb_acc'               : 'kg m-2',   # Baseflow groundwater runoff
        'SoilMoi0_10cm_inst'    : 'kg m-2',   # Soil moisture (0-10 cm)
        'SoilMoi10_40cm_inst'   : 'kg m-2',   # Soil moisture (10-40 cm)
        'SoilMoi40_100cm_inst'  : 'kg m-2',   # Soil moisture (40-100 cm) 
        'SoilMoi100_200cm_inst' : 'kg m-2',   # Soil moisture (100-200 cm)  
#       'SoilMoist_S_tavg'      : 'kg m-2',      # Surface soil moisture
#       'SoilMoist_RZ_tavg'     : 'kg m-2',      # Root zone soil moisture
#       'SoilMoist_P_tavg'      : 'kg m-2',      # Profile soil moisture
        'CanopInt_inst'         : 'kg m-2',      # Plant canopy surface water
#       'CanopInt_tavg'         : 'kg m-2',      # Plant canopy surface water
        #-----------------------------------------------------------------------
        'ACond_tavg'         : 'm s-1',       # Aerodynamic conductance
        'Wind_f_inst'        : 'm s-1' }      # Wind speed

        try:
            units = umap[ var_name ]
        except:
            print("SORRY, Don't know units for this GLDAS var.")
            units = 'not available (see docs)'

    #----------------------     
    # Print final message
    #----------------------  
    print( ' ')
    print( 'Variable name  =', var_name )
    print( 'Variable units =', units)
    print( 'max(variable)  =', vmax)
    print( 'min(variable)  =', vmin, '(possible nodata)' )
    print( 'bad_count =', bad_count )
    print( 'n_grids   =', count )
    print( 'Finished saving data to rts file.')
    print( ' ')
    
#   create_rts_from_nc_files()
#-------------------------------------------------------------------    
def create_rts_from_chirps_files( rts_file='TEST.rts',
           time_interval_hours=6.0,   ################### 
           resample_algo='bilinear', BARO_60=False, 
           DEM_bounds=None, DEM_xres_sec=None, DEM_yres_sec=None,
           DEM_ncols=None, DEM_nrows=None,
           VERBOSE=False, SILENT=False, NON_NEGATIVE=True):

    #------------------------------------------------------
    # Note: See function above for resampling algorithms.
    #------------------------------------------------------
    # For info on GDAL constants, see:
    # https://gdal.org/python/osgeo.gdalconst-module.html
    #------------------------------------------------------
    if (BARO_60):
        DEM_bounds = [ 34.221249999999, 7.362083333332, 36.450416666666, 9.503749999999]
        ### DEM_bounds = [ 7.362083333332, 34.221249999999, 9.503749999999, 36.450416666666]
        DEM_xres_sec = 60.0   # (60 arcsecs = 60/3600 degrees)
        DEM_yres_sec = 60.0
        DEM_ncols  = 134
        DEM_nrows  = 129
        
    #--------------------------------------------------------    
    # See: Pongo_30sec DEM: MINT_2019/__Regions/South_Sudan
    #      Adjusted this later as the Lol-Kuru.
    #--------------------------------------------------------     
#     if (PONGO_30):
#         DEM_bounds = [24.079583333333,  6.565416666666, 27.379583333333, 10.132083333333]
#         ## DEM_bounds = [6.565416666666, 24.079583333333,  10.132083333333, 27.379583333333]
#         DEM_xres   = 1./120   # (30 arcsecs = 30/3600 degrees)
#         DEM_yres   = 1./120   # (30 arcsecs = 30/3600 degrees)
#         DEM_ncols  = 396
#         DEM_nrows  = 428

    if (DEM_bounds is None) or (DEM_xres_sec is None) or \
       (DEM_yres_sec is None) or (DEM_ncols is None) or (DEM_nrows is None):
        print('ERROR:  All DEM bounding box information is required.')
        return

    #---------------------------------------------
    # Convert DEM_xres_sec to DEM_xres_deg, etc.
    #---------------------------------------------
    DEM_xres_deg = (DEM_xres_sec / 3600.0)
    DEM_yres_deg = (DEM_yres_sec / 3600.0)
  
    #-------------------------    
    # Open RTS file to write
    #-------------------------
    rts_unit = open( rts_file, 'wb' )

    #-----------------------------------------------
    # Get list of all files in working directory
    #  e.g. rfe_gdas.bin.2014100100 (no extension)
    #-----------------------------------------------
    file_list = glob.glob( 'rfe*' )
    file_list = sorted(file_list)

    #-------------------------------------------------
    # It may not be wise to change the nodata values
    # because GDAL may not recognize them as nodata
    # during spatial interpolation, etc.
    # This could be the cause of "artifacts".
    #-------------------------------------------------
    # A nodata value of zero may be okay for rainfall
    # rates, but is not good in general.
    #-------------------------------------------------    
    #### rts_nodata = -9999.0    ######
    # rts_nodata = 0.0
    rts_nodata = None   # (do not change nodata values)
   
    count = 0
    vmax  = -1e8
    vmin  = 1e8
    bad_count = 0
    BAD_FILE  = False
    tmp_file  = 'TEMP.tif'
        
    #-----------------------------------    
    # CHIRPS Africa bounding box, etc.
    #-----------------------------------
    # chirps_minlat = -40.05
    # chirps_maxlat = 40.05
    # chirps_minlon = -20.05
    # chirps_maxlon = 55.05
    chirps_ulx      = -20.05
    chirps_uly      = 40.05
    chirps_xres_deg = 0.1
    chirps_yres_deg = -0.1     ###########   NEED NEGATIVE HERE
    chirps_ncols    = 751
    chirps_nrows    = 801
    chirps_dtype    = 'float32'
    chirps_nodata   = -1.0

    print('Working...') 
    for chirps_file in file_list:
        #-------------------------------------
        # Read CHIRPS binary file (rainfall)
        #-------------------------------------
        print()
        print('chirps_file =', chirps_file)
        chirps_unit = open(chirps_file, 'rb')
        grid = np.fromfile( chirps_unit, dtype=chirps_dtype )      
        grid = grid.reshape( chirps_nrows, chirps_ncols )
        chirps_unit.close()
        
        #------------------------------------        
        # Flip the y-axis to make row-major
        #------------------------------------
        # Confirmed correct: July 2021;
        # ocean is in lower-right corner
        #------------------------------------
        grid = np.flipud( grid )

        #---------------------------    
        # Byteswap from MSB to LSB
        #---------------------------
        grid = grid.byteswap()

        #------------------------------------------------------    
        # Convert units from mm/interval to mmph (2020-09-24)
        #------------------------------------------------------
        grid /= time_interval_hours

        #-------------------------------------        
        # Re-save grid as geotiff (for GDAL)
        #-------------------------------------
        if not(SILENT):  print('  resaving as GeoTIFF...')
        new_tif_file = 'tmp_' + chirps_file + '.tif'
        save_grid_to_geotiff( new_tif_file, grid,
                              chirps_ulx, chirps_uly,
                              chirps_xres_deg, chirps_yres_deg,
                              nodata=chirps_nodata)
                          
        #---------------------------------------------- 
        # Read grid & metadata from geotiff with GDAL
        #----------------------------------------------
        ds_in = gdal.Open( new_tif_file )
        grid1 = ds_in.ReadAsArray()
        #------------------------------------  
        # Note: This gives min & max before
        #       clipping to DEM bounds
        #------------------------------------
#         gmax  = grid1.max()
#         gmin  = grid1.min()
#         vmax  = max( vmax, gmax )
#         vmin  = min( vmin, gmin )
        band  = ds_in.GetRasterBand(1)
        ### band.SetNoDataValue( chirps_nodata )   ###### SET ABOVE
        ## chirps_nodata = band.GetNoDataValue()

        if (VERBOSE):
            print( '===============================================================')
            print( 'count =', (count + 1) )
            print( '===============================================================')
            print( 'grid1: min   =', grid1.min(), ', max =', grid1.max() )
            print( 'grid1.shape  =', grid1.shape )
            print( 'grid1.dtype  =', grid1.dtype )
            print( 'grid1 nodata =', chirps_nodata, '(in original file)' )
            ## w  = np.where(grid1 > chirps_nodata)
            w  = np.where(grid1 != chirps_nodata)
            nw = w[0].size
            print( 'grid1 # data =', nw)
            print( ' ' )
              
        #--------------------------------------        
        # Use gdal.Info() to print/check info
        #--------------------------------------
        ## print( gdal.Info( ds_in ) )
        ## print( '===============================================================')

        #-----------------------------------------------        
        # Check if the bounding boxes actually overlap
        #-----------------------------------------------
        if not(SILENT):  print('  Comparing bounds...')
        ds_bounds = get_raster_bounds( ds_in, VERBOSE=False )
        if (bounds_disjoint( ds_bounds, DEM_bounds )):
            print( '###############################################')
            print( 'WARNING: Bounding boxes do not overlap.')
            print( '         New grid will contain only nodata.')
            print( '###############################################')
            print( 'count =', count )
            print( 'file  =', new_tif_file )
            print( 'ds_bounds  =', ds_bounds )
            print( 'DEM_bounds =', DEM_bounds )
            print( ' ')
            bad_count += 1
            BAD_FILE = True
            sys.exit()

        #-------------------------------------------
        # Clip and resample data to the DEM's grid
        # then save to a temporary GeoTIFF file.
        #-------------------------------------------
        if not(BAD_FILE):
            print('  Regridding to DEM grid...')
            grid2 = gdal_regrid_to_dem_grid( ds_in, tmp_file,
                         DEM_bounds, DEM_xres_deg, DEM_yres_deg,
                         nodata=rts_nodata, IN_MEMORY=False,
                         RESAMPLE_ALGO=resample_algo )
            
            if (VERBOSE):
                print( 'grid2: min  =', grid2.min(), ', max =', grid2.max() )
                print( 'grid2.shape =', grid2.shape )
                print( 'grid2.dtype =', grid2.dtype )
                if (rts_nodata is not None):
                    ## w  = np.where(grid2 > rts_nodata)
                    w  = np.where(grid2 != rts_nodata)
                    nw = w[0].size
                    print( 'grid2 # data =', nw)
                    print()
        else:
            grid2 = np.zeros( (DEM_nrows, DEM_ncols), dtype='float32' )
            if (rts_nodata is not None):
                grid2 += rts_nodata
 
        #-------------------------------
        # Convert [kg m-2 s-1] to mmph
        #-------------------------------
        w = (grid2 != chirps_nodata)   # (boolean array)
        grid2[w] = grid2[w] * 3600.0   # (preserve nodata)

        #------------------------------------
        # Option to remove the new_tif_file
        #------------------------------------
        os.remove(new_tif_file)
        
        #-------------------------  
        # Write grid to RTS file
        #-------------------------
        grid2 = np.float32( grid2 )
        if (NON_NEGATIVE):
            grid2[ grid2 < 0 ] = 0.0
        vmax  = max( vmax, grid2.max() )
        vmin  = min( vmin, grid2.min() )
        grid2.tofile( rts_unit )
        count += 1
        if not(SILENT):  print('count =', count)

    #----------------------------------
    # Close the RTS file and tmp_file
    #----------------------------------
    rts_unit.close()
    if (os.path.exists( tmp_file )):
        os.remove( tmp_file )   # TEMP.tif

    #---------------------------------------
    # Create an RTI file for this RTS file
    #---------------------------------------
    rti = rti_files.make_info(grid_file=rts_file,
              ncols=DEM_ncols, nrows=DEM_nrows,
              xres=DEM_xres_sec, yres=DEM_yres_sec,
              #--------------------------------------
              data_source='TopoFlow 3.6 Regrid',
              data_type='FLOAT',
              byte_order=rti_files.get_rti_byte_order(),
              pixel_geom=0,
              zres=0.001, z_units='unknown',
              y_south_edge=DEM_bounds[1],
              y_north_edge=DEM_bounds[3],
              x_west_edge=DEM_bounds[0],
              x_east_edge=DEM_bounds[2],
              box_units='DEGREES')
    rti_files.write_info( rts_file, rti )

    #---------------------------
    # Determine variable units
    #---------------------------
    units = 'mmph'  # (converted from mass flux)

    #----------------------     
    # Print final message
    #----------------------  
    print( ' ')
    # print( 'Variable name  =', var_name )
    print( 'Variable units =', units)
    print( 'max(variable)  =', vmax)
    print( 'min(variable)  =', vmin, '(possible nodata)' )
    print( 'bad_count =', bad_count )
    print( 'n_grids   =', count )
    print( 'Finished saving CHIRPS data to rts file.')
    print( ' ')
    
#   create_rts_from_chirps_files()
#-------------------------------------------------------------------
# def regrid_geotiff_to_dem(in_file=None, out_file=None, 
#                           DEM_bounds=None, DEM_xres=None, DEM_yres=None ):
# 
#     #---------------------------------------------------------------
#     # Note:  DEM_bounds = [dem_xmin, dem_ymin, dem_xmax, dem_ymax]
#     #        Give xres, yres in decimal degrees for Geographic.
#     #        gdal.Warp() clips to a bounding box, and can also
#     #        resample to a different resolution.
#     #        gdal.Translate() is faster for simple clipping.
#     #---------------------------------------------------------------
#     if (in_file == None):
#         #-----------------------------------------------------------
#         # Use Pongo_30sec DEM as a test, which works well.
#         # However,  the soil data has same resolution (xres, yres)
#         # as the DEM, of 30 arcseconds.  In addition, grid cells
#         # outside of South Sudan have NODATA values.
#         #-----------------------------------------------------------
#         in_file    = 'SLTPPT_M_sl1_1km_South Sudan.tiff'
#         out_file   = 'Pongo_SLTPPT_sl1.tiff'
#         DEM_bounds = [24.079583333333,  6.565416666666, 27.379583333333, 10.132083333333 ]
#         DEM_xres   = 1./120   # (30 arcsecs = 30/3600 degrees)
#         DEM_yres   = 1./120   # (30 arcsecs = 30/3600 degrees)
#     
#     f1 = gdal.Open( in_file, gdal.GA_ReadOnly )
#     ## data_xres = f1.RasterXsize
#     ### data_yres = f1.RasterYsize
#     # print( f1.RasterCount )
#     # print( data_xres, data_yres )
#   
#     out_unit = gdal.Warp( out_file, f1,
#         format = 'GTiff',  # (output format string)
#         outputBounds=DEM_bounds, xRes=DEM_xres, yRes=DEM_yres,
#         resampleAlg = gdal.GRA_Bilinear )
#         ## resampleAlg = gdal.GRA_NearestNeighbour ) 
#         # (near, bilinear, cubic, cubicspline, lanczos, average, etc.)
#     out_unit = None   # Close out_file
# 
#     #--------------------------------------------------------  
#     # This shows some of the other keywords to gdal.Warp.
#     #--------------------------------------------------------      
#     # WarpOptions(options=[], format=None, outputBounds=None,
#     # outputBoundsSRS=None, xRes=None, yRes=None, targetAlignedPixels=False,
#     # width=0, height=0, srcSRS=None, dstSRS=None, srcAlpha=False,
#     # dstAlpha=False, warpOptions=None, errorThreshold=None,
#     # warpMemoryLimit=None, creationOptions=None, outputType=GDT_Unknown,
#     # workingType=GDT_Unknown, resampleAlg=None, srcNodata=None,
#     # dstNodata=None, multithread=False, tps=False, rpc=False,
#     # geoloc=False, polynomialOrder=None, transformerOptions=None,
#     # cutlineDSName=None, cutlineLayer=None, cutlineWhere=None,
#     # cutlineSQL=None, cutlineBlend=None, cropToCutline=False,
#     # copyMetadata=True, metadataConflictValue=None,
#     # setColorInterpretation=False, callback=None, callback_data=None)
#     # 
#     # Create a WarpOptions() object that can be passed to gdal.Warp()
#     # Keyword arguments are : options --- can be be an array of strings,
#     # a string or let empty and filled from other keywords.
# 
#      
# #    regrid_geotiff_to_dem()
#-------------------------------------------------------------------  
# def resave_grid_to_geotiff( ds_in, new_file, grid1, nodata=nodata ):
# 
#     if (nodata is not None):
#         new_nodata = -9999.0
#         grid1[ grid1 != nodata ] = new_nodata
#     
#     ##### raster = gdal.Open( nc_file )
#     raster = ds_in
#     ncols = raster.RasterXSize
#     nrows = raster.RasterYSize
# 
#     geotransform = raster.GetGeoTransform()
#     originX      = geotransform[0]
#     originY      = geotransform[3]
#     pixelWidth   = geotransform[1]
#     pixelHeight  = geotransform[5]
# 
#     driver = gdal.GetDriverByName('GTiff')
#     outRaster = driver.Create(new_file, ncols, nrows, 1, gdal.GDT_Float32)
#     outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
#     outband = outRaster.GetRasterBand(1)
#     outband.WriteArray( grid1 )
#     outRasterSRS = osr.SpatialReference()
#     outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
#     outRaster.SetProjection(outRasterSRS.ExportToWkt())
#     outband.FlushCache()
# 
# #   resave_grid_to_geotiff()
#-------------------------------------------------------------------  
# def rename_nc_files( txt_file_in=None, txt_file_out=None ):
# 
#     if (txt_file_in is None):
#         txt_file_in  = 'subset_GPM_3IMERGHH_06_20200604_191328.txt'
#     if (txt_file_out is None):
#         txt_file_out = 'subset_GPM_2014-08.txt'
# 
#     in_unit  = open(txt_file_in,  'r')
#     out_unit = open(txt_file_out, 'w')
#     
#     # Read a line, modify it, write new line to txt_file_out
#     while (True):
#         line  = in_unit.readline()
#         if (line == ''):
#             break
#         words = line.split('?')
#         ncols = len( words )
#         #-------------------------------
#         # Discard part after first "?"
#         #-------------------------------
#         new_line = words[0]
#         out_unit.write( new_line + '\n' )
#      
#     #------------------
#     # Close the files
#     #------------------
#     in_unit.close()
#     out_unit.close()
# 
# #   rename_nc_files()
#-------------------------------------------------------------------

def create_aorc_ncgs_forcings_files(nc_file_in=None, grid_info=None, time_info=None,
           resample_algo='bilinear', SILENT=False, VERBOSE=False, NC4=False):

    #-----------------------------------------    
    # Use a temp file in memory or on disk ?
    #-----------------------------------------
    IN_MEMORY = False
    if (IN_MEMORY):
        tmp_file = '/vsimem/TEMP.tif'
    else:
        tmp_file = 'TEMP.tif'
    
    #-------------------------------------------------
    # It may not be wise to change the nodata values
    # because GDAL may not recognize them as nodata
    # during spatial interpolation, etc.
    # This could be the cause of "artifacts".
    #-------------------------------------------------
    # A nodata value of zero may be okay for rainfall
    # rates, but is not good in general.
    #-------------------------------------------------    
    #### rts_nodata = -9999.0    ######
    # rts_nodata = 0.0
    rts_nodata = None   # (do not change nodata values)
    count = 0
    vmax  = -1e8
    vmin  = 1e8
    bad_count = 0
    BAD_FILE  = False
    tif_file  = 'TEMP1.tif'
    
    #-----------------------------------------    
    # Get time info from input nc file
    #-----------------------------------------
    ncgs_in = ncgs_files.ncgs_file()
    ncgs_in.open_file(file_name=nc_file_in)
    time_values = ncgs_in.ncgs_unit.variables['time'][:]
    time_units = ncgs_in.ncgs_unit.variables['time'].units
    ncgs_in.close_file()

    #-----------------------------------------    
    # Create time info for file
    #-----------------------------------------
    if time_units == "minutes since 2016-10-01 00:00:00":
        time_origin_datetime_obj = time_utils.datetime.datetime(2016,10,1,0,0,0)
        time_units = 'minutes'
    else:
        print('ERROR unrecognized units for time dimension in ',ncgs_in.ncgs_unit.filepath())
    start_datetime = time_utils.convert_times_from_minutes_to_datetime(times_min=[time_values[0]],origin_datetime_obj=time_origin_datetime_obj)
    end_datetime = time_utils.convert_times_from_minutes_to_datetime(times_min=[time_values[-1]],origin_datetime_obj=time_origin_datetime_obj)
    start_datetime = start_datetime[0]
    end_datetime = end_datetime[0]

    #-----------------------------------------    
    # Compare time info objects
    #-----------------------------------------
    start_datetime_obj_file = time_utils.get_datetime_obj_from_one_str(start_datetime)
    end_datetime_obj_file = time_utils.get_datetime_obj_from_one_str(end_datetime)
    start_datetime_obj_simulation = time_utils.get_datetime_obj_from_one_str(time_info.start_datetime)
    end_datetime_obj_simulation = time_utils.get_datetime_obj_from_one_str(time_info.end_datetime)
    if time_units != time_info.duration_units:
        print('ERROR specified time units do not match time units in ',nc_file_in)
    if start_datetime_obj_file > start_datetime_obj_simulation:
        print('ERROR specified start datetime is not within datetime range of ',nc_file_in)
    if end_datetime_obj_file < end_datetime_obj_simulation:
        print('ERROR specified end datetime is not within datetime range of ',nc_file_in)

    #-----------------------------------------    
    # Find first and second time step index values 
    #-----------------------------------------
    time_index_first = None
    time_index_second = None
    for i in range(len(time_values)):
        i_datetime_str = time_utils.convert_times_from_minutes_to_datetime(times_min=[time_values[i]],origin_datetime_obj=start_datetime_obj_file)
        i_datetime_obj = time_utils.get_datetime_obj_from_one_str(i_datetime_str[0])
        if i_datetime_obj == start_datetime_obj_simulation:
            time_index_first = i
        if i_datetime_obj == start_datetime_obj_simulation + timedelta(minutes=time_info.dt):
            time_index_second = i
        if time_index_first is not None and time_index_second is not None:
            break
    if time_index_first is None:
        print('ERROR could not find first simulation datetime in ',nc_file_in)
    if time_index_second is None:
        print('ERROR could not find second simulation datetime in ',ncgs_unit.filepath())

    #-----------------------------------------    
    # Determine needed time indices 
    #-----------------------------------------
    time_index = [i for i in range(time_index_first,len(time_values),time_index_second-time_index_first)]

    #-----------------------------------------    
    # Create array of datetime strings for netcdf datetime variable
    #-----------------------------------------
    var_datetimes = [time_utils.convert_times_from_minutes_to_datetime(times_min=[time_values[i]],origin_datetime_obj=start_datetime_obj_file)[0] for i in time_index]

    #---------------------------------------------
    # Define DEM spatial info for regridding
    #---------------------------------------------
    DEM_bounds = [grid_info.x_west_edge,grid_info.y_south_edge,grid_info.x_east_edge,grid_info.y_north_edge]
    DEM_xres_sec = grid_info.xres
    DEM_yres_sec = grid_info.yres
    DEM_xres_deg = (grid_info.xres / 3600.0)
    DEM_yres_deg = (grid_info.yres / 3600.0)
    DEM_ncols = grid_info.ncols
    DEM_nrows = grid_info.nrows

    #---------------------------------------------
    # Define dictionary to hold grids stacks and variable units
    #---------------------------------------------
    dt_grids = dict()
    dt_var_units = dict()

    #---------------------------------------------
    # Iterate over variables
    #---------------------------------------------
    var_names = ['Psurf_f_inst','Qair_f_inst','Rainf_tavg','Tair_f_inst','Wind_f_inst']
    for var_name in var_names:
        
        print('\nEvaluating variable: ',var_name)

        #-----------------------------------------    
        # Read var using gdal (to faciliate regridding)
        #-----------------------------------------
        if (VERBOSE):  print('Reading grid of values for variable: ',var_name)
        (ds_in, grid1, nc_nodata) = gdal_open_nc_file( nc_file_in, var_name, VERBOSE=VERBOSE)

        if (VERBOSE):
            print( '===============================================================')
            print( 'count =', (count + 1) )
            print( '===============================================================')
            print( 'grid1: min   =', grid1.min(), ', max =', grid1.max() )
            print( 'grid1.shape  =', grid1.shape )
            print( 'grid1.dtype  =', grid1.dtype )
            print( 'grid1 nodata =', nc_nodata, '(in original nc_file)' )
            ### w  = np.where(grid1 > nc_nodata)
            w  = np.where(grid1 != nc_nodata)
            nw = w[0].size
            print( 'grid1 # data =', nw)
            print( ' ' )
                
        #--------------------------------------        
        # Use gdal.Info() to print/check info
        #--------------------------------------
        ## print( gdal.Info( ds_in ) )
        ## print( '===============================================================')

        #-----------------------------------------------        
        # Check if the bounding boxes actually overlap
        #-----------------------------------------------
        if not(SILENT):  print('Comparing bounds...')
        ds_bounds = get_raster_bounds( ds_in, VERBOSE=False )
        if (bounds_disjoint( ds_bounds, DEM_bounds )):
            print( '###############################################')
            print( 'WARNING: Bounding boxes do not overlap.')
            print( '         New grid will contain only nodata.')
            print( '###############################################')
            print( 'count =', count )
            print( 'file  =', nc_file )
            print( 'ds_bounds  =', ds_bounds )
            print( 'DEM_bounds =', DEM_bounds )
            print( ' ')
            bad_count += 1
            BAD_FILE = True

        #-------------------------------------------
        # Clip and resample data to the DEM's grid
        # then save to a temporary GeoTIFF file.
        #-------------------------------------------
        if not(BAD_FILE):
            print('Regridding to DEM grid...')
            grid2 = gdal_regrid_to_dem_grid( ds_in, tmp_file,
                            DEM_bounds, DEM_xres_deg, DEM_yres_deg,
                            nodata=rts_nodata, IN_MEMORY=IN_MEMORY,
                            RESAMPLE_ALGO=resample_algo )
            # This was already done in gdal_regrid_to_dem_grid()
            ## ds_in = None   # Close the tmp_file
            ## if (IN_MEMORY):
            ##    gdal.Unlink( tmp_file )
            if (VERBOSE):
                print( 'grid2: min  =', grid2.min(), ', max =', grid2.max() )
                print( 'grid2.shape =', grid2.shape )
                print( 'grid2.dtype =', grid2.dtype )
                if (rts_nodata is not None):
                    w  = np.where(grid2 != rts_nodata)
                    nw = w[0].size
                    print( 'grid2 # data =', nw)
                    print()
        else:
            grid2 = np.zeros( (DEM_nrows, DEM_ncols), dtype='float32' )
            if (rts_nodata is not None):
                grid2 += rts_nodata

        #-----------------------------------------    
        # close the read-in dataset
        #-----------------------------------------
        ds_in = None  # (close ds)

        #-----------------------------------------    
        # forcings unit conversions
        #-----------------------------------------
        if (var_name == 'Rainf_tavg'):
            print('Converting units: [kg m-2 s-1] to [mmph]...')
            w = (grid2 != nc_nodata)  # (boolean array)
            grid2[w] *= 3600.0        # (preserve nodata)
            var_units = 'mmph'
        if (var_name == 'Tair_f_inst'):
            print('Converting units: Kelvin to Celsius...')
            w = (grid2 != nc_nodata)  # (boolean array)
            grid2[w] -= 273.15        # (preserve nodata)  
            var_units = 'C'               
        if (var_name == 'Psurf_f_inst'):
            print('Converting units: Pa to mbar...')
            w = (grid2 != nc_nodata)  # (boolean array)
            grid2[w] /= 100.0         # (preserve nodata)
            var_units = 'mbar'

        #-----------------------------------------    
        # name conversions
        #-----------------------------------------
        topoflow_var_name = None
        if (var_name == 'Rainf_tavg'):
            topoflow_var_name = 'P'
        elif (var_name == 'Tair_f_inst'):
            topoflow_var_name = 'T_air'
        elif (var_name == 'Psurf_f_inst'):
            topoflow_var_name = 'p0'
        elif (var_name == 'Lwnet_tavg'):
            topoflow_var_name = 'Qn_LW'
        elif (var_name == 'Swnet_tavg'):
            topoflow_var_name = 'Qn_SW'
        elif (var_name == 'Wind_f_inst'):
            topoflow_var_name = 'uz'
        if topoflow_var_name is None:
            topoflow_var_name = var_name

        #-----------------------------------------    
        # Slice the time dimension
        #-----------------------------------------
        grid2 = grid2[time_index,:,:]

        #-----------------------------------------    
        # Save grid and var units
        #-----------------------------------------
        dt_grids[topoflow_var_name] = grid2
        dt_var_units[topoflow_var_name] = var_units

    #-----------------------------------------    
    # Calculate grid for relative humidity
    #-----------------------------------------
    dt_grids['rh'] = met_utils.relative_humidity(dt_grids['Qair_f_inst'], dt_grids['T_air'], dt_grids['p0'])
    dt_var_units['rh'] = 'None'

    #-----------------------------------------    
    # Set minimum wind speed (topoflow may error when wind speed = 0)
    #-----------------------------------------
    dt_grids['uz'][dt_grids['uz'] == 0.] = numpy.nextafter(numpy.single(0), numpy.single(1))

    #-----------------------------------------    
    # iterate over variables and create output nc files
    #-----------------------------------------
    print('')
    for topoflow_var_name in ['P','T_air','p0','uz','rh']:
        print('Creating nc file for variable: ',topoflow_var_name)

        #-----------------------------------------    
        # initialize new file
        #-----------------------------------------
        ncgs_out = ncgs_files.ncgs_file()
        output_file_name = os.path.join(os.path.dirname(nc_file_in),'AORC_'+topoflow_var_name+'.nc')
        if os.path.isfile(output_file_name):
            os.remove(output_file_name)
        OK = ncgs_out.open_new_file(file_name=output_file_name,grid_info=grid_info,time_info=time_info,var_name=topoflow_var_name,
                                units_name=dt_var_units[topoflow_var_name],time_units='minutes',time_res=time_info.dt,
                                OVERWRITE_OK=True,MAKE_RTI=False)
        if not OK:
            print('ERROR unable to successfully create the new file ',output_file_name)

        #------------------------------------
        # Populate the variable
        #------------------------------------
        ncgs_out.ncgs_unit.variables[topoflow_var_name][:,:,:] = dt_grids[topoflow_var_name]

        #------------------------------------
        # Populate time variable
        #------------------------------------
        ncgs_out.ncgs_unit.variables['datetime'][:] = np.array(var_datetimes)
        ncgs_out.ncgs_unit.variables['time'][:] = [i+1 for i in range(len(time_index))]

        #------------------------------------
        # Set other attributes
        #------------------------------------
        ncgs_out.ncgs_unit.variables[topoflow_var_name].long_name    = topoflow_var_name
        ncgs_out.ncgs_unit.variables[topoflow_var_name].dx           = grid_info.xres
        ncgs_out.ncgs_unit.variables[topoflow_var_name].dy           = grid_info.yres 
        ncgs_out.ncgs_unit.variables[topoflow_var_name].y_south_edge = grid_info.y_south_edge
        ncgs_out.ncgs_unit.variables[topoflow_var_name].y_north_edge = grid_info.y_north_edge
        ncgs_out.ncgs_unit.variables[topoflow_var_name].x_west_edge  = grid_info.x_west_edge
        ncgs_out.ncgs_unit.variables[topoflow_var_name].x_east_edge  = grid_info.x_east_edge        

        #------------------------------------
        # Close the file
        #------------------------------------
        ncgs_out.close()