
#  Copyright (c) 2019, Scott D. Peckham
#  November 2019

#-------------------------------------------------------------------

#  test1()

#  read_from_geotiff()
#  read_from_netcdf()
#  read_from_rtg()

#  print_gdal_grid_info()

#  rtg_type_from_gdal_type()
#  make_rti_for_gdal_grid()

#  save_geotiff_as_rtg()

#-------------------------------------------------------------------

import gdal, osr  ## ogr
import os.path
# import glob

from . import rtg_files, rti_files

#-------------------------------------------------------------------
def test1():

    test_dir = '/Users/peckhams/TF_Tests/Import/'

    #-----------------------   
    # Read a GeoTIFF file.
    #----------------------- 
    tif_file = test_dir + 'Baro_Gam_1min_DEM.tif'
    rti_file = test_dir + 'Baro_Gam_1min.rti'
    grid = read_from_geotiff( tif_file, REPORT=True, rti_file=rti_file)
    print()
    
    #----------------------    
    # Read a netCDF file.
    #----------------------
    nc_file = test_dir + 'GPM_Rain_Grid_TEST.nc'
    grid = read_from_netcdf( nc_file, REPORT=True)    
    print()
    
    #------------------------------   
    # Read an RTG file (with RTI)
    #------------------------------
    rtg_file = test_dir + 'Baro_Gam_1min_DEM.rtg'
    grid = read_from_rtg( rtg_file, REPORT=True)
        
#   test1()
#-------------------------------------------------------------------
def read_from_geotiff( tif_file, REPORT=False, rti_file=None):

    #-----------------------------------------------------
    # Note:  tif_file may need to have ".tif" extension.
    #-----------------------------------------------------
    if (not(os.path.exists( tif_file ))):
        print('ERROR:  Wrong filename or file does not exist.')
        return -1

    ds = gdal.Open( tif_file, gdal.GA_ReadOnly )
    grid = ds.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
 
    if (REPORT):
        print_gdal_grid_info( tif_file, grid, ds )
        
    if (rti_file is not None):
        make_rti_for_gdal_grid(tif_file, grid, ds, rti_file=rti_file)

    #--------------------          
    # Close the tif_file
    #--------------------
    ds = None
    
    #------------------------------------
    # Return the DEM (or other 2D grid)
    #------------------------------------
    return grid

#   read_from_geotiff()
#-------------------------------------------------------------------
def read_from_netcdf( nc_file, layer_name='', REPORT=False,
                      rti_file=None):

    #--------------------------------------------------------
    # Note:  If layer_name is not specified it defaults to
    #        the null string.  GDAL will successfully read
    #        the file, at least if there is only one named
    #        variable or layer or band in the file.  It may
    #        simply use the first one it finds in the file.
    #--------------------------------------------------------
    if (not(os.path.exists( nc_file ))):
        print('ERROR:  Wrong filename or file does not exist.')
        return -1

    ds = gdal.Open("NETCDF:{0}:{1}".format(nc_file, layer_name))
    grid = ds.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    
    if (REPORT):
        print_gdal_grid_info( nc_file, grid, ds ) 

    if (rti_file is not None):
        make_rti_for_gdal_grid(nc_file, grid, ds, rti_file)

    #--------------------          
    # Close the nc_file
    #--------------------  
    ds = None

    #------------------------------------
    # Return the DEM (or other 2D grid)
    #------------------------------------
    return grid
    
#   read_from_netcdf()
#-------------------------------------------------------------------
def read_from_rtg( rtg_file, REPORT=False):

    if (not(os.path.exists( rtg_file ))):
        print('ERROR:  Wrong filename or file does not exist.')
        return -1
  
    info = rti_files.read_info( rtg_file, SILENT=False, REPORT=REPORT)
    
    grid = rtg_files.read_grid( rtg_file, info,
                     RTG_type=info.data_type,
                     ## RTG_type='FLOAT',
                     REPORT=False, SILENT=True )
  
    if (REPORT):
        print('Finished reading file:')
        print( rtg_file )
        print( 'grid.min() =', grid.min() )
        print( 'grid.max() =', grid.max() )
        print ('grid.shape =', grid.shape )
        print()
        
    #------------------------------------
    # Return the DEM (or other 2D grid)
    #------------------------------------                     
    return grid

#   read_from_rtg()
#-------------------------------------------------------------------
def print_gdal_grid_info( filename, grid, ds ):

    #----------------------
    # Get grid dimensions
    #----------------------
    shp = grid.shape
    ncols = shp[1]
    nrows = shp[0]
    #--------------------------
    # ncols = ds.RasterXSize
    # nrows = ds.RasterYSize
    
    #----------------------------------------------------------  
    # Notice the strange order or parameters here is CORRECT.
    # It is not:  ulx, xres, xskew, uly, yres, yskew
    #----------------------------------------------------------
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    lrx = ulx + (ds.RasterXSize * xres)
    lry = uly + (ds.RasterYSize * yres)

    print('Finished reading file:')
    print( filename )
    print('Grid info from GDAL:')
    print('ncols, nrows =', ncols, ',', nrows)
    print('xres, yres   =', xres, ',', abs(yres) )
    print('----------------------------------')
    print('ulx, uly     =', ulx, ',', uly)
    print('lrx, lry     =', lrx, ',', lry)
    print('xskew, yskew =', xskew, ',', yskew)
    print('----------------------------------')
    print( 'grid.min()  =', grid.min() )
    print( 'grid.max()  =', grid.max() )
    print()

#   print_gdal_grid_info()
#-------------------------------------------------------------------
def rtg_type_from_gdal_type( ds ):

#     gdal.GDT_Unknown = 0
#     gdal.GDT_Byte    = 1
#     gdal.GDT_UInt16  = 2	
#     gdal.GDT_Int16   = 3	
#     gdal.GDT_UInt32  = 4 	
#     gdal.GDT_Int32   = 5	
#     gdal.GDT_Float32 = 6	
#     gdal.GDT_Float64 = 7	
#     gdal.GDT_CInt16  = 8 	
#     gdal.GDT_CInt32  = 9 	
#     gdal.GDT_CFloat32 = 10	
#     gdal.GDT_CFloat64 = 11

    type_map = {
    0 : 'UNKNOWN',
    1 : 'BYTE',
    2 : 'UINT',
    3 : 'INTEGER',
    4 : 'ULONG',
    5 : 'LONG',
    6 : 'FLOAT',
    7 : 'DOUBLE' }

    gdal_data_type_code = ds.GetRasterBand(1).DataType
    rtg_type = type_map[ gdal_data_type_code ]
        
    return rtg_type
    
#   rtg_type_from_gdal_type()
#-------------------------------------------------------------------
def make_rti_for_gdal_grid( grid_file, grid, ds, rti_file,
                            z_units='METERS', box_units='DEGREES'):
    
    #----------------------
    # Get grid dimensions
    #----------------------
    shp = grid.shape
    ncols = shp[1]
    nrows = shp[0]
    #--------------------------
    # ncols = ds.RasterXSize
    # nrows = ds.RasterYSize

    #--------------------    
    # Get the data type
    #--------------------
    data_type = rtg_type_from_gdal_type( ds )

    #-----------------------------  
    # Get the machine byte order
    #-----------------------------
    byte_order = rti_files.get_rti_byte_order()

    data_source = 'TopoFlow & GDAL'
    zres        = 0.01     # [meters]    ###### ASSUMED            
    gmin = grid.min()
    gmax = grid.max()
        
    #----------------------------------------------------------  
    # Notice the strange order or parameters here is CORRECT.
    # It is not:  ulx, xres, xskew, uly, yres, yskew
    #----------------------------------------------------------
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    lrx = ulx + (ds.RasterXSize * xres)
    lry = uly + (ds.RasterYSize * yres)

    #-------------------------------------------------------
    # ASSUME for now that grid uses Geographic coordinates
    # and not a fixed-length projection like UTM
    #-------------------------------------------------------
    GEOGRAPHIC = True
    if (GEOGRAPHIC):
        pixel_geom   = 0   # (use 1 for fixed-length, like UTM)
        x_west_edge  = ulx
        x_east_edge  = lrx
        y_south_edge = lry
        y_north_edge = uly
        #-----------------------------------------------------------
        # Get UTM zone for the center, but may span multiple zones
        #-----------------------------------------------------------
        mid_lon = (x_west_edge + x_east_edge) / 2.0
        UTM_zone = rti_files.get_utm_zone( mid_lon )   # (string)
    
        xres_sec = xres * 3600       # [deg] -> [arcsec]
        yres_sec = abs(yres) * 3600  # [deg] -> [arcsec]
    else:
        #-------------------------------------------
        # Will need to do something different here
        #-------------------------------------------
        dum = 0
        return
    
    grid_info = rti_files.make_info(grid_file=grid_file, ncols=ncols, nrows=nrows,
                  xres=xres_sec, yres=yres_sec, data_source=data_source,
                  byte_order=byte_order, pixel_geom=pixel_geom,
                  data_type=data_type, zres=zres, z_units=z_units,
                  y_south_edge=y_south_edge, y_north_edge=y_north_edge,
                  x_west_edge=x_west_edge, x_east_edge=x_east_edge,
                  box_units=box_units, gmin=gmin, gmax=gmax,
                  UTM_zone=UTM_zone )
   
    #-----------------------------------------
    # Save georeferencing in RTI file format
    #-----------------------------------------
    if (rti_file is None):
        dot_pos  = grid_file.rfind('.')
        prefix   = grid_file[:dot_pos]
        rti_file = prefix + '.rti'
    rti_files.write_info( rti_file, grid_info )
    
    print('Finished making RTI file for:')
    print('  ' + grid_file )
    print()
 
#   make_rti_for_gdal_grid()
#-------------------------------------------------------------------
def save_geotiff_as_rtg(tif_file, rtg_file, rti_file, REPORT=False):

    grid = read_from_geotiff( tif_file, REPORT=REPORT, rti_file=rti_file )
    grid_info = rti_files.read_info( rti_file )
    out_type  = 'FLOAT'   # (regardless of grid_info.data_type)
    rtg_files.write_grid( grid, rtg_file, grid_info, RTG_type=out_type, SILENT=True)

    #------------------------------------------------------------    
    # Set SILENT=True above, then print a more specific message
    #------------------------------------------------------------
    print('Finished saving GeoTIFF file as RTG to:')
    print('  ' + rtg_file )
    print()
    
#   save_geotiff_as_rtg()
#-------------------------------------------------------------------


