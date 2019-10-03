
#  Copyright (c) 2019, Scott D. Peckham
#  August 2019

#-------------------------------------------------------------------

#  regrid_geotiff_to_dem()
#  read_nc_grid()
#  gdal_open_nc_file()
#  get_raster_bounds()
#  bounds_disjoint()
#  gdal_regrid_to_dem_grid()
#  resave_grid_to_geotiff()
#  create_rts_from_nc_files()

#-------------------------------------------------------------------
#  Set up a "tf4" conda environment (TopoFlow 4.0)
#-------------------------------------------------------------------
#  % conda create --name tf4
#  % conda activate tf4
#  % conda install -c conda-forge gdal   (to read geotiff)
#  % conda install -c conda-forge scipy   (for gamma function)
#         (use conda-forge vs. anaconda; broken?)
#  % conda install -c conda-forge pydap
#  % conda install dask

#-------------------------------------------------------------------
import numpy as np
import gdal, osr  ## ogr
import glob

# import os.path
#-------------------------------------------------------------------
def regrid_geotiff_to_dem(in_file=None, out_file=None, 
                          DEM_bounds=None, DEM_xres=None, DEM_yres=None ):

    #---------------------------------------------------------------
    # Note:  DEM_bounds = [dem_xmin, dem_ymin, dem_xmax, dem_ymax]
    #        Give xres, yres in decimal degrees for Geographic.
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
        in_file    = 'SLTPPT_M_sl1_1km_South Sudan.tiff'
        out_file   = 'Pongo_SLTPPT_sl1.tiff'
        DEM_bounds = [24.079583333333,  6.565416666666, 27.379583333333, 10.132083333333 ]
        DEM_xres   = 1./120   # (30 arcsecs = 30/3600 degrees)
        DEM_yres   = 1./120   # (30 arcsecs = 30/3600 degrees)
    
    f1 = gdal.Open( in_file, gdal.GA_ReadOnly )
    ## data_xres = f1.RasterXsize
    ### data_yres = f1.RasterYsize
    # print( f1.RasterCount )
    # print( data_xres, data_yres )
  
    out_unit = gdal.Warp( out_file, f1,
        format = 'GTiff',  # (output format string)
        outputBounds=DEM_bounds, xRes=DEM_xres, yRes=DEM_yres,
        resampleAlg = gdal.GRA_Bilinear )
        ## resampleAlg = gdal.GRA_NearestNeighbour ) 
        # (near, bilinear, cubic, cubicspline, lanczos, average, etc.)
    out_unit = None   # Close out_file

    #-------------------------------------------------------- 
    # Example:  Use gdal.Translate to clip to bounding box.
    #--------------------------------------------------------   
    # ds = gdal.Open('original.tif')
    # ds = gdal.Translate('new.tif', ds, projWin = [-75.3, 5.5, -73.5, 3.7])
    # ds = None

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

     
#    regrid_geotiff_to_dem()
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
def read_nc_grid( nc_file=None, var_name='HQprecipitation',
                  REPORT=False):

    if (nc_file == None):
        nc_file = 'TEST.nc4'

    ds = gdal.Open("NETCDF:{0}:{1}".format(nc_file, layer_name))
    grid = ds.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    ds = None  # (close ds)
        
    if (REPORT):
        print(( 'grid.min() =', grid.min() ))
        print(( 'grid.max() =', grid.max() ))
        print(('grid.shape =', grid.shape ))

    return grid
    
    #--------------------
    # This doesn't work
    #--------------------
#     ds = gdal.Open( nc_file )
#     # print( ds.RasterCount )
#     # print( ds.RasterYSize, ds.RasterXsize )
#     data = ds.ReadAsArray()
#     # print( data.shape )
#     print( data.min() )
#     print( data.max() )
#     ds = None  # (close ds)
            
#   read_nc_grid()
#-------------------------------------------------------------------
# def read_nc_as_array( nc_file=None, var_name='HQprecipitation',
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

    ### ds_in = gdal.Open("NETCDF:{0}:{1}".format(nc_file, var_name), gdal.GA_ReadOnly )
    ds_in  = gdal.Open("NETCDF:{0}:{1}".format(nc_file, var_name) )
    band   = ds_in.GetRasterBand(1)
    nodata = band.GetNoDataValue()

    g1 = band.ReadAsArray()
    ## g1 = ds_in.ReadAsArray(0, 0, ds_in.RasterXSize, ds_in.RasterYSize)

    if (VERBOSE):
        print(( 'grid1: min =', g1.min(), 'max =', g1.max() ))
        print(( 'grid1.shape =', g1.shape ))
        print(( 'grid1.dtype =', g1.dtype ))
        print(( 'grid1 nodata =', nodata ))
        print( ' ' )

    return (ds_in, g1, nodata)

# gdal_open_nc_file()
#------------------------------------------------------------------- 
def get_raster_bounds( ds, VERBOSE=True):

    #-------------------------------------------------------------
    # Note:  The bounds depend on the map projection and are not
    # necessarily a Geographic bounding box of lons and lats.    
    #-------------------------------------------------------------
    # ulx = upper left x  = xmin
    # uly = upper left y  = ymax
    # lrx = lower right x = xmax
    # lry = lower right y = ymin
    #-----------------------------
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    lrx = ulx + (ds.RasterXSize * xres)
    lry = uly + (ds.RasterYSize * yres)

    if (VERBOSE):
        print(('ulx, uly   =', ulx, uly))
        print(('lrx, lry   =', lrx, lry))
        print(('xres, yres = ', xres, yres))
        print(('xskew, yskew =', xskew, yskew))
        print('----------------------------------')

    return [ulx, lry, lrx, uly]  # [xmin, ymin, xmax, ymax]

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
def gdal_regrid_to_dem_grid( ds_in, tmp_file, 
         nodata, DEM_bounds, DEM_xres, DEM_yres,
         RESAMPLE_ALGO='bilinear', VERBOSE=False):

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
        outputBounds=DEM_bounds, xRes=DEM_xres, yRes=DEM_yres,
        srcNodata=nodata,      ########
        ### dstNodata=nodata,  ########
        resampleAlg = resample_algo )

    grid = ds_tmp.ReadAsArray()
    
    ds_tmp = None   # Close tmp_file
   
    return grid

#   gdal_regrid_to_dem_grid()
#-------------------------------------------------------------------  
def resave_grid_to_geotiff( ds_in, new_file, grid1, nodata ):

    new_nodata = -9999.0
    grid1[ grid1 <= nodata ] = new_nodata
    
    ##### raster = gdal.Open( nc_file )
    raster = ds_in
    ncols = raster.RasterXSize
    nrows = raster.RasterYSize

    geotransform = raster.GetGeoTransform()
    originX      = geotransform[0]
    originY      = geotransform[3]
    pixelWidth   = geotransform[1]
    pixelHeight  = geotransform[5]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(new_file, ncols, nrows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray( grid1 )
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

#   resave_grid_to_geotiff()     
#-------------------------------------------------------------------    
def create_rts_from_nc_files( rts_file='TEST.rts',
                              IN_MEMORY=False, VERBOSE=True):

    #------------------------------------------------------
    # For info on GDAL constants, see:
    # https://gdal.org/python/osgeo.gdalconst-module.html
    #------------------------------------------------------  
    if (rts_file == 'TEST.rts'):
        #-----------------------------------------------------------
        # Use Pongo_30sec DEM as a test, which works well.
        # However,  the soil data has same resolution (xres, yres)
        # as the DEM, of 30 arcseconds.  In addition, grid cells
        # outside of South Sudan have NODATA values.
        #-----------------------------------------------------------
        DEM_bounds = [24.079583333333,  6.565416666666, 27.379583333333, 10.132083333333 ]
        DEM_xres   = 1./120   # (30 arcsecs = 30/3600 degrees)
        DEM_yres   = 1./120   # (30 arcsecs = 30/3600 degrees)

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
    ## nc_file_list = glob.glob( '*.nc4' )
    nc_file_list = glob.glob( '*.nc' )
    var_name = "HQprecipitation"    # HQ = high quality;  1/2 hourly
    count = 0
    bad_count = 0
  
    for nc_file in nc_file_list:
        #-------------------------------
        # Open the original netCDF file
        #--------------------------------
        (ds_in, grid1, nodata) = gdal_open_nc_file( nc_file, var_name, VERBOSE=True)
        print(( gdal.Info( ds_in ) ))

        #-----------------------------------------------        
        # Check if the bounding boxes actually overlap
        #-----------------------------------------------
        ds_bounds = get_raster_bounds( ds_in )
        if (bounds_disjoint( ds_bounds, DEM_bounds )):
            print( '###############################################')
            print( 'WARNING: Bounding boxes do not overlap.')
            print( '         New grid will contain only nodata.')
            print( '###############################################')
            print(( 'ds_bounds  =', ds_bounds ))
            print(( 'DEM_bounds =', DEM_bounds ))
            print( ' ')
            bad_count += 1

        #-------------------------------------------
        # Replace nodata value and save as GeoTIFF
        #-------------------------------------------
#         new_file = 'TEMP2.tif'
#         resave_grid_to_geotiff( ds_in, new_file, grid1, nodata )
#         ds_in = None  # Close the nc_file
#         ds_in = gdal.Open( new_file )   # Open the GeoTIFF file; new nodata

        #-------------------------------------------
        # Clip and resample data to the DEM's grid
        # then save to a temporary GeoTIFF file.
        #-------------------------------------------
        grid2 = gdal_regrid_to_dem_grid( ds_in, tmp_file,
                    nodata, DEM_bounds, DEM_xres, DEM_yres,
                    RESAMPLE_ALGO='bilinear' )
        print(( 'grid2: min =', grid2.min(), 'max =', grid2.max() ))
        print(( 'grid2.shape =', grid2.shape ))
        print(( 'grid2.dtype =', grid2.dtype ))
        print( ' ')
        ds_in = None   # Close the tmp_file
               
        #--------------------------------------------
        # Read resampled data from tmp GeoTIFF file
        #--------------------------------------------
        ds_tmp = gdal.Open( tmp_file )
        ## ds_tmp = gdal.Open( tmp_file, gdal.GA_ReadOnly  )
        ## print( gdal.Info( ds_tmp ) )
        grid3  = ds_tmp.ReadAsArray()
        print(( 'grid3: min, max =', grid3.min(), grid3.max() ))
        print(( 'grid3.shape =', grid3.shape))
        print(( 'grid3.dtype =', grid3.dtype))
        ds_tmp = None   # Close tmp file
   
        if (IN_MEMORY):
            gdal.Unlink( tmp_file )
       
        #-------------------------  
        # Write grid to RTS file
        #-------------------------
        grid3 = np.float32( grid3 )
        ## rts_unit.write( grid3 )
        grid3.tofile( rts_unit )
        count += 1
        if (VERBOSE):
            print(( 'count, min, max =', count, grid3.min(), grid3.max() ))   #######
            
        if (count == 300):  ##################################
            break

    #---------------------
    # Close the RTS file
    #---------------------
    rts_unit.close()

    print( ' ')
    print(( 'bad_count =', bad_count ))
    print(( 'n_grids   =', count ))
    print( 'Finished saving data to rts file.')
    print( ' ')
    
#   create_rts_from_nc_files()
#-------------------------------------------------------------------   

          
 