
# Copyright (c) 2019, Scott D. Peckham
# August - November 2019
# See:
# https://csdms.colorado.edu/wiki/Model_help:TopoFlow-Soil_Properties_Page
#-------------------------------------------------------------------

#  test1()
#  test2()
#  test3()

#  get_nodata_values()
#  read_soil_grid_files()        # (as RTG)
#  check_soil_grid_values()
#  convert_soil_grid_values()    # ISRIC to Wosten PTF vars
#  read_isric_soil_grid_files()  # (as TIFF)

#  get_rtg_file_prefix()
#  regrid_isric_soil_grid_files()
#  transform_isric_soil_grid_files()

#  wosten_theta_s()
#  wosten_K_s()
#  wosten_alpha()
#  wosten_n()
#  wosten_l()

#  get_wosten_vars()
#  capillary_length_G()
#  get_tBC_from_vG_vars()
#  save_soil_hydraulic_vars()

#-------------------------------------------------------------------
import numpy as np
import gdal, osr  ## ogr
from scipy.special import gamma
from . import regrid as rg
from . import import_grid as ig
from . import rti_files
from . import rtg_files

## from . import soil_trans_BC as soil

import glob, os

# import os.path

#-------------------------------------------------------------------
def test1():

    directory = '/Users/peckhams/ISRIC_Files/'
    transform_isric_soil_grid_files( directory=directory,
        IN_MEMORY=False, VERBOSE=False, LOL_1MIN_TEST=True)
              
#   test1()
#-------------------------------------------------------------------
def test2(BARO_1MIN=True, LOL_1MIN=False):

    if (BARO_1MIN):
        site_prefix = 'Baro_Gam_1min'
        out_dir = '/Users/peckhams/ISRIC_Files/Ethiopia/'
        in_dir  = out_dir + '_isric_soil_data/'
        # Bounds = [ minlon, minlat, maxlon, maxlat ]   
        out_bounds = [ 34.22125, 7.3704166666, 36.43791666666, 9.50375]
        out_xres_sec = 60    # [arcsecs]
        out_yres_sec = 60    # [arcsecs]
        ## out_ncols  = 133  # (will be computed)
        ## out_nrows  = 128  # (will be computed)
               
    if (LOL_1MIN):
        site_prefix = 'Lol-Kuru_1min'
        out_dir = '/Users/peckhams/ISRIC_Files/South_Sudan/'
        in_dir  = out_dir + '_isric_soil_data/'
        # Bounds = [ minlon, minlat, maxlon, maxlat ] 
        out_bounds = [ 23.9954166666, 6.5329166666, 28.0120833333, 9.56625]
        out_xres_sec = 60   # [arcsecs]
        out_yres_sec = 60   # [arcsecs]
        ## out_ncols  = 241  # (will be computed)
        ## out_nrows  = 182  # (will be computed)
        
    save_soil_hydraulic_vars( site_prefix=site_prefix,
         in_dir=in_dir, out_dir=out_dir,
         out_bounds=out_bounds, VERBOSE=True,
         out_xres_sec=out_xres_sec, out_yres_sec=out_yres_sec,
         RESAMPLE_ALGO='bilinear')
           
#   test2()      
#-------------------------------------------------------------------
def test3():

    topsoil = 1
    
    print('Case of: C=90, S=5, OM=5, D=1.5')
    alpha = wosten_alpha( 0.9, 0.05, 0.05, 1.5, topsoil )
    # alpha = wosten_alpha( 90.0, 5.0, 5.0, 1.5, topsoil )
    print('alpha = ' + str(alpha) )
    print()
    #-----------------------------------------------------
    print('Case of: C=5, S=90, OM=5, D=1.5')
    alpha = wosten_alpha( 0.05, 0.9, 0.05, 1.5, topsoil )
    # alpha = wosten_alpha( 5.0, 90.0, 5.0, 1.5, topsoil )
    print('alpha = ' + str(alpha) )
    print()
    #-----------------------------------------------------
        
#   test3()
#-------------------------------------------------------------------
def get_nodata_values():

    #-------------------------------------------------
    # Replace zero values with another nodata value.
    # Zero values cause NaN or Inf in Wosten vars.
    # For %, 101.0 is an impossible value.
    # For D, use 10 [g / cm^3].
    #-------------------------------------------------
    # Density of water          = 1.00  [g / cm^3]
    # Density of silt           = 1.33  [g / cm^3]
    # Density of loose sand     = 1.442 [g / cm^3]
    # Density of dry sand       = 1.602 [g / cm^3]
    # Density of compacted clay = 1.746 [g / cm^3]
    # Density of most rocks     = 2.65  [g / cm^3]
    # Density of pure iron      = 7.87  [g / cm^3]
    # Density of pure actinium  = 10.07 [g / cm^3]  
    # Density of pure lead      = 11.34 [g / cm^3]
    # Density of pure gold      = 19.20 [g / cm^3]
    #-------------------------------------------------
    nd_values = {
    'C'  : 101.0,  # [%]
    'S'  : 101.0,  # [%]
    'X'  : 101.0,  # [%]
    'OM' : 101.0,  # [%]
    'D'  : 10.0,   # [g / cm^3]
    #------------------------------
    'theta' : -9999.0 }  # unitless, in [0,1]

#     nd_values = {
#     'C'  : -9999.0,  # [%]
#     'S'  : -9999.0,  # [%]
#     'X'  : -9999.0,  # [%]
#     'OM' : -9999.0,  # [%]
#     'D'  : -9999.0 }  # [g / cm^3]
#     #------------------------------
#     'theta' : -9999.0 }  # unitless, in [0,1]
         
#     nd_values = {
#     'C'  : 0.001,  # [%]
#     'S'  : 0.001,  # [%]
#     'X'  : 0.001,  # [%]
#     'OM' : 0.001,  # [%]
#     'D'  : 0.001 }  # [g / cm^3]
#     #------------------------------
#     'theta' : -9999.0 }  # unitless, in [0,1]
        
    return nd_values

#   get_nodata_values()
#-------------------------------------------------------------------
def read_soil_grid_files( site_prefix=None, directory=None,
         res_str='1km', layer=1, REPORT=True):

    layer_str = '_sl' + str(layer) + '_'
    if (directory is None):
        directory = '/Users/peckhams/ISRIC_Files/'
    os.chdir( directory )

    #-------------------------------------------------------------
    # Read soil property data from a set RTG files.  These
    # grids contain data from ISRIC - SoilGrids, but have been
    # clipped and resampled to a DEM grid to be used by TopoFlow.
    # These basic soil properties are needed to compute several
    # soil hydraulic properties using the pedotransfer functions
    # (PTFs) in Wosten (1998).
    #-------------------------------------------------------------
    # Read percent clay
    # Read percent silt
    # Read percent organic matter
    # Read bulk density
    #-------------------------------------------------------------
    # Values are available for 7 different depths:
    #   0.0, 0.05, 0.15, 0.30, 0.60, 1.00, 2.00
    #   sl1, sl2,  sl3,  sl4,  sl5,  sl6,  sl7
    #----------------------------------------------------- 
    suffix  = layer_str + res_str + '.rtg'
    C_file  = site_prefix + '_CLYPPT_M' + suffix  # (C)
    S_file  = site_prefix + '_SLTPPT_M' + suffix  # (S)
    OC_file = site_prefix + '_ORCDRC_M' + suffix  # (OC)
    D_file  = site_prefix + '_BLDFIE_M' + suffix  # (D)   

    #------------------------------------------------
    # Read header_file info for all of the RTG files
    #------------------------------------------------
    header_file = (site_prefix + '.rti')
    grid_info = rti_files.read_info( header_file, REPORT=False)
    
    #-------------------------------
    # Read the soil property grids
    #---------------------------------------------------
    # Read ISRIC clay mass fraction [kg/kg], as %.
    # This is the same as "C" in Wosten PTF equations.
    #--------------------------------------------------- 
    C = rtg_files.read_grid( C_file, grid_info, RTG_type='FLOAT' )
    
    #---------------------------------------------------
    # Read ISRIC silt mass fraction [kg/kg], as %.
    # This is the same as "S" in Wosten PTF equations.
    #---------------------------------------------------    
    S = rtg_files.read_grid( S_file, grid_info, RTG_type='FLOAT' )
    
    #----------------------------------------------------
    # Read ISRIC organic carbon content, units [g / kg]
    #----------------------------------------------------
    OC = rtg_files.read_grid( OC_file, grid_info, RTG_type='FLOAT' )    

    #-------------------------------------------------
    # Read ISRIC soil bulk density, units [kg / m^3]
    #-------------------------------------------------
    D = rtg_files.read_grid( D_file, grid_info, RTG_type='FLOAT' )

    #--------------------------------------
    # Check if grid values are reasonable
    #--------------------------------------
    check_soil_grid_values(C, S, OC, D, REPORT=REPORT)

    #------------------------------------------------
    # Convent units to those needed by Wosten PTFs.
    #------------------------------------------------
    (C,S,OM,D) = convert_soil_grid_values(C, S, OC, D)
        
    #-------------------------------------------------
    # Replace -9999.0 with another nodata value to
    # avoid NaN or Inf in Wosten PTF equations.
    # For %, 101.0 is an impossible value.
    # For D, use 10 [g / cm^3] (impossible for soil).
    #---------------------------------------------------
    # We have used -9999.0 (ISRIC) up to this point,
    # but that's not good for the wosten PTF equations.
    #---------------------------------------------------
    nodata = get_nodata_values()
    C[ C <= 0 ]   = nodata['C']   # [%]
    S[ S <= 0 ]   = nodata['S']   # [%]
    OM[ OM <= 0 ] = nodata['OM']  # [%]
    D[ D <= 0]    = nodata['D']   # [g / cm^3]
      
    #------------------------------------
    # Return Wosten input vars as tuple
    #------------------------------------   
    return (C, S, OM, D)
    
#   read_soil_grid_files()
#-------------------------------------------------------------------
def check_soil_grid_values(C, S, OC, D, REPORT=True):

    #-----------------------------
    # Get all the mins and maxes
    #-----------------------------
    Cmin  = C.min()
    Cmax  = C.max()
    Smin  = S.min()
    Smax  = S.max()
    OCmin = OC.min()
    OCmax = OC.max()
    Dmin  = D.min()
    Dmax  = D.max()

    #--------------------------------- 
    # Option to print mins and maxes
    #---------------------------------
    if (REPORT):
        print('min(C),  max(C)  = ' + str(Cmin)  + ', ' + str(Cmax) )
        print('min(S),  max(S)  = ' + str(Smin)  + ', ' + str(Smax) )
        print('min(OC), max(OC) = ' + str(OCmin) + ', ' + str(OCmax) )
        print('min(D),  max(D)  = ' + str(Dmin)  + ', ' + str(Dmax) )
        print()

    #--------------------------------------
    # Check if grid values are reasonable
    #--------------------------------------
    if ((Cmin < 0.0) or (Cmax > 100.0)):
        print('WARNING in read_soil_grid_files:')
        print('   Some values in C grid are out of range.')
        print('   min(C) = ' + str(Cmin) )
        print('   max(C) = ' + str(Cmax) )
        print()
    #--------------------------------------------------------
    if ((Smin < 0.0) or (Smax > 100.0)):
        print('WARNING in read_soil_grid_files:')
        print('   Some values in S grid are out of range.')
        print('   min(S) = ' + str(Smin) )
        print('   max(S) = ' + str(Smax) )
        print()
    #--------------------------------------------------------
    # ISRIC OC units are [g/kg]
    # ISRIC D units are [kg / m^3]
    #--------------------------------------------------------
    if ((OCmin < 0.0) or (OCmax > 1000.0)):
        print('WARNING in read_soil_grid_files:')
        print('   Some values in OC grid are out of range.')
        if (OCmin < 0.0):
            print('   min(OC) = ' + str(OCmin) )
            print('   Possible nodata value.')
        if (OCmax > 1000.0):
            print('   max(OC) = ' + str(OCmax) )
            ## print('   Possible nodata value.')
        print()
    #------------------------------------------------
    # Density of water          = 1000  [kg / m^3]
    # Density of silt           = 1330  [kg / m^3]
    # Density of loose sand     = 1442  [kg / m^3]
    # Density of dry sand       = 1602  [kg / m^3]
    # Density of compacted clay = 1746  [kg / m^3]  
    # Density of most rocks     = 2650  [kg / m^3]
    #------------------------------------------------         
    if ((Dmin < 0.0) or (Dmax > 2650)):
        print('WARNING in read_soil_grid_files:')
        print('   Some values in D grid are out of range.')
        if (Dmin < 0.0):
            print('   min(D) = ' + str(Dmin) )
            print('   Possible nodata value.')
        if (Dmax > 2650.0):
            print('   max(D) = ' + str(Dmax) )
            ## print('   Possible nodata value.')
        print()
     
#   check_soil_grid_values()
#-------------------------------------------------------------------
def convert_soil_grid_values(C, S, OC, D, VERBOSE=False):

    #----------------------------------------------
    # Convert ISIRC OC [g / kg] to Wosten OM [%].
    #--------------------------------------------------
    # OC = ISRIC organic carbon content in [g / kg]
    # However, Wosten PTFs need OM, where:
    # OM = organic mass fraction in [kg / kg], as %.
    #---------------------------------------------------
    # First convert [g/kg] to [kg/kg] (mass fraction).
    # Then multiply by 100 to convert [kg/kg] to %.
    # Note:  mass fraction [kg/kg] is in [0,1].
    #--------------------------------------------------
    #     g       1 kg        OM   kg
    # OM ---  *  ------  =   ----  --
    #     kg     1000 g      1000  kg
    #--------------------------------------------------
    # (OM/1000) [kg/kg] * 100 = (OM/10) %.
    #--------------------------------------------------
    w1 = (OC == -9999.0)
    OC = (OC / 10.0)
    OC[w1] = -9999.0   # (preserve nodata value)
    
    #---------------------------------------------------
    # Convert organic carbon to organic matter content
    #----------------------------------------------------------
    # NOTE:  ISRIC SoilGrids "ORCDRC" is "soil organic carbon
    #        content (fine earth fraction)", in [g / kg].
    #        But Wosten PTFs need soil organic *matter*.
    #        A conversion factor of 1.724 is given in 
    #        Wosten et al. (1998, p. 28) and is attributed to:
    #        Nelson and Sommers (1982).  Also found online.
    #        About 58% of organic matter is carbon, and
    #        100 / 58 = 1.724.
    #---------------------------------------------------------- 
    OM = OC * 1.724   # (Wosten et al (1998), p. 28)
    OM[w1] = -9999.0
    
    #-------------------------------------------------
    # Convert ISRIC [kg / m^3] to Wosten [g / cm^3]
    #-------------------------------------------------
    # Wosten D = bulk density, [g / cm^3]
    #-------------------------------------------------
    #    kg     1000 g      (1 m)^3       D    g
    # D ---- * -------- * ------------ = ---- ---
    #    m3      1 kg      (100 cm)^3    1000  cm3
    #-------------------------------------------------
    w1 = (D == -9999.0)
    D = (D / 1000.0)
    D[w1] = -9999.0    # (preserve nodata value)
    
    return (C, S, OM, D)

#   convert_soil_grid_values()
#-------------------------------------------------------------------
def read_isric_soil_grid_files( layer=1, region_str=None,
                                res_str='1km' ):

    if (region_str is None):
        region_str = '_South Sudan.tiff'
    layer_str  = '_sl' + str(layer) + '_'
    
    #-------------------------------------------------------------
    # Read soil property data from ISRIC - SoilGrids files,
    # as needed to compute Wosten (1998) pedotransfer functions.
    # SoilGrids files are in GeoTIFF format.
    #-------------------------------------------------------------
    # Another option is to use the rasterio package:
    #
    # import rasterio
    # with rasterio.open('sample.tif') as r:
    #     ar = r.read()    
    #-------------------------------------------------------------
    # Read percent clay
    # Read percent silt
    # Read percent organic matter
    # Read bulk density
    #-------------------------------------------------------------
    # Here we are using 1km grid cell data for the entire
    # country of South Sudan, downloaded from:
    #  https://soilgrids.org/#!/?layer=ORCDRC_M_sl2_250m&vector=1
    # Go to soilgrids.org.
    # Click on the "Download data" bitmap, on the right.
    # Scroll down to the section with Coverage ID (layer 1km)
    #   and Country droplist.  Choose "South Sudan".
    # Choose layers and click "Download" button.
    # Values are available for 7 different depths:
    #   0.0, 0.05, 0.15, 0.30, 0.60, 1.00, 2.00
    #   sl1, sl2,  sl3,  sl4,  sl5,  sl6,  sl7
    #-------------------------------------------------------------
    suffix  = layer_str + res_str + region_str
    C_file  = site_prefix + '_CLYPPT_M' + suffix  # (C, clay)
    S_file  = site_prefix + '_SLTPPT_M' + suffix  # (S, silt)
    OM_file = site_prefix + '_ORCDRC_M' + suffix  # (OM)
    D_file  = site_prefix + '_BLDFIE_M' + suffix  # (D)
    #-------------------------------------------------------  
    # This is not needed by Wosten pedotransfer functions.
    #-------------------------------------------------------   
#     X_file  = site_prefix + '_SNDPPT_M' + suffix  # (X, sand)
        
    #-------------------------------
    # Read the soil property grids
    #------------------------------------------------
    # Wosten C = clay mass fraction, [kg/kg], as %.
    # Same units in ISRIC and Wosten.
    #------------------------------------------------
    f1 = gdal.Open( C_file )
    # print( f1.RasterCount )
    # print( f1.RasterYSize, f1.RasterXsize )
    C  = f1.ReadAsArray()   # (byte type)
    f1 = None  # (close f1)
    #------------------------------------------------
    # Wosten S = silt mass fraction, [kg/kg], as %.
    # Same units in ISRIC and Wosten.
    #------------------------------------------------   
    f2 = gdal.Open( S_file )
    S  = f2.ReadAsArray()   # (byte type)
    f2 = None  # (close f2)
    #----------------------------------------------------------
    # Wosten OM = organic matter mass fraction, [kg/kg], as %
    # ISRIC units = [g / kg]  
    # Convert ISIRC [g / kg] to Wosten [%].
    # First convert [g/kg] to [kg/kg] (mass fraction),
    # Then multiply by 100 to convert [kg/kg] to %.
    # Note:  mass fraction [kg/kg] is in [0,1].
    #--------------------------------------------------
    # OM [g/kg] * (1 kg / 1000 g) = (OM/1000) [kg/kg]
    # (OM/1000) [kg/kg] * 100 = (OM/10) %.
    #--------------------------------------------------   
    f3 = gdal.Open( OM_file )
    OM  = f3.ReadAsArray()
    OM = (OM / 10.0)
    f3 = None  # (close f3)
    #--------------------------------------
    # Wosten D = bulk density, [g / cm^3]
    #----------------------------------------------------
    # Convert ISRIC [kg / m^3] to Wosten [g / cm^3]
    # D [kg / m^3] * (1000 g / 1 kg) * (1 m / 100 cm)^3
    # D [kg / m^3] * (1 / 1000) = D [g / cm^3]
    #---------------------------------------------------- 
    f4 = gdal.Open( D_file )
    D = f4.ReadAsArray()    # (bulk density, kg/m3)
    D = (D / 1000.0)
    f4 = None  # (close f4)

    #------------------------------------------------
    # Wosten X = sand mass fraction, [kg/kg], as %.
    # Same units in ISRIC and Wosten.
    #------------------------------------------------   
#     f5 = gdal.Open( X_file )
#     X  = f5.ReadAsArray()   # (byte type)
#     f5 = None  # (close f5)
    
    #------------------------------------
    # Return Wosten input vars as tuple
    #------------------------------------
    return (C, S, OM, D)
    
#   read_isric_soil_grid_files()
#-------------------------------------------------------------------  
def get_rtg_file_prefix( tiff_file, site_prefix ):

    p = tiff_file.split('_')
    n = len(p)
    #-------------------------------------------------------    
    out_suffix = ('_' + p[0])       # var name string
    if (n > 1):
        out_suffix += ('_' + p[1])  # unknown code string
    if (n > 2):
        out_suffix += ('_' + p[2])  # layer number string
    if (n > 3):
        out_suffix += ('_' + p[3])  # resolution string
#     if (n > 4):
#         out_suffix += ('_' + p[4])  # location string
    #-------------------------------------------------------
    rtg_prefix= site_prefix + out_suffix
    return rtg_prefix
    
#   get_rtg_file_prefix()
#-------------------------------------------------------------------    
def regrid_isric_soil_grid_files( site_prefix,
           in_dir=None, out_dir=None,
           out_bounds=None, out_xres_sec=None, out_yres_sec=None,
           RESAMPLE_ALGO='bilinear', REPORT=True):

    if (in_dir is None):
        return
    os.chdir( in_dir )
    if (out_dir[-1] != '/'):
        out_dir = out_dir + '/'
         
    #--------------------------------------------------
    # Get list of all TIFF files in working directory
    #--------------------------------------------------
    # Note:  glob.glob does NOT return a sorted list.
    # https://stackoverflow.com/questions/6773584/
    #    how-is-pythons-glob-glob-ordered
    #--------------------------------------------------
    ### tif_file_list = glob.glob( '*.tif' )
    tif_file_list = glob.glob( '*.tiff' )
    tif_file_list = sorted( tif_file_list )  #######
    n_grids = 0
    #----------------------------------------
    ## bad_box_count = 0
    ## out_nodata = -9999.0       #############
             
    for tif_file in tif_file_list:                                

        #----------------------
        # Construct filenames
        #----------------------
        rtg_prefix = get_rtg_file_prefix( tif_file, site_prefix )
        rtg_file   = out_dir + rtg_prefix + '.rtg' 
        rti_file   = out_dir + rtg_prefix + '.rti'
        #--------------------------------------------
        # This file gets overwritten multiple times
        # BAD IDEA BECAUSE data_types DIFFER.
        #--------------------------------------------
#         rti_file     = out_dir + site_prefix + '.rti'   ###### TEMPORARY FIX
        out_tif_file = out_dir + rtg_prefix + '.tif'
        
        rg.regrid_geotiff(in_file=tif_file,
               out_file=out_tif_file, 
               out_bounds=out_bounds,
               out_xres_sec=out_xres_sec,
               out_yres_sec=out_yres_sec,
               RESAMPLE_ALGO='bilinear', REPORT=True)

        #--------------------------------------------------        
        # Read regridded GeoTIFF & save as RTG (with RTI)
        #--------------------------------------------------
        ig.save_geotiff_as_rtg(out_tif_file, rtg_file, rti_file, REPORT=False)
        n_grids += 1

    print()
    print('Finished regridding ISRIC soil grid files.')
    print('   Soil grid files also saved to RTG format.')
    print('   Number of grids = ' + str(n_grids) )
    # print( '   Number outside of model domain = ' + str(bad_box_count) )
    print()       
 
#   regrid_isric_soil_grid_files()                 
#-------------------------------------------------------------------    
def transform_isric_soil_grid_files( site_prefix,
              in_dir=None, out_dir=None,
              out_bounds=None,
              out_xres_sec=None, out_yres_sec=None,
              RESAMPLE_ALGO='bilinear', 
              IN_MEMORY=False, REPORT=False ):

    #------------------------------------------------------
    # For info on GDAL constants, see:
    # https://gdal.org/python/osgeo.gdalconst-module.html
    #------------------------------------------------------ 
    if (in_dir is None):
        return
    os.chdir( in_dir )
    if (out_dir[-1] != '/'):
        out_dir = out_dir + '/'  
        
    #-----------------------------------------------------------  
    # Typically, out_bounds = DEM_bounds,
    # out_xres_sec = DEM_xres_sec, out_yres_sec = DEM_yres_sec
    #-----------------------------------------------------------
    out_xres_deg = out_xres_sec / 3600.0
    out_yres_deg = out_yres_sec / 3600.0
        
    #-----------------------------------------    
    # Use a temp file in memory or on disk ?
    #-----------------------------------------
    if (IN_MEMORY):
        tmp_file = '/vsimem/TEMP.tif'
    else:
        tmp_file = 'TEMP.tif'

    #--------------------------------------------------
    # Get list of all TIFF files in working directory
    #--------------------------------------------------
    tiff_file_list = glob.glob( '*.tiff' )
    tiff_file_list = sorted( tiff_file_list )   ######
    n_grids = 0
    bad_box_count = 0
    out_nodata = -9999.0       #############
             
    for tiff_file in tiff_file_list:
        #------------------------------------------------     
        # Print filename before open in case of failure
        #------------------------------------------------
        if (REPORT):
            print( '===============================================================')
            print( 'ISRIC File =', tiff_file)
            print( '===============================================================')                               
        ds_in = gdal.Open( tiff_file )
        grid1 = ds_in.ReadAsArray()
        band  = ds_in.GetRasterBand(1)
        tiff_nodata = band.GetNoDataValue()

        if (REPORT):
            print( 'grid1.min()  =', grid1.min() )
            print( 'grid1.max()  =', grid1.max() )
            print( 'grid1.shape  =', grid1.shape )
            print( 'grid1.dtype  =', grid1.dtype )
            print( 'grid1 nodata =', tiff_nodata )
            #-----------------------------------------
            w  = np.where(grid1 == tiff_nodata)
            nw = w[0].size
            print( 'nodata count =', nw)
            print()
              
        #--------------------------------------        
        # Use gdal.Info() to print/check info
        #--------------------------------------
        ## print( gdal.Info( ds_in ) )
        ## print( '===============================================================')

        #-----------------------------------------------        
        # Check if the bounding boxes actually overlap
        #-----------------------------------------------
        BAD_BOX = False
        in_bounds = rg.get_raster_bounds( ds_in, VERBOSE=False )
        if (rg.bounds_disjoint( in_bounds, out_bounds )):
            print( '###############################################')
            print( 'WARNING: Bounding boxes do not overlap.')
            print( '         New grid will contain only nodata.')
            print( '###############################################')
            print( 'count =', n_grids )
            print( 'file  =', nc_file )
            print( 'in_bounds  =', in_bounds )
            print( 'out_bounds =', out_bounds )
            print( ' ')
            bad_box_count += 1
            BAD_BOX = True

        #-------------------------------------------
        # Clip and resample data to the DEM's grid
        # then save to a temporary GeoTIFF file.
        #-------------------------------------------
        if not(BAD_BOX):
            grid2 = rg.gdal_regrid_to_dem_grid( ds_in, tmp_file,
                        out_nodata, 
                        out_bounds, out_xres_deg, out_yres_deg,
                        ### DEM_bounds, DEM_xres, DEM_yres,
                        RESAMPLE_ALGO=RESAMPLE_ALGO)
            if (REPORT):
                print( 'grid2.min()  =', grid2.min() )
                print( 'grid2.max()  =', grid2.max() )
                print( 'grid2.shape  =', grid2.shape )
                print( 'grid2.dtype  =', grid2.dtype )
                print( 'grid2 nodata =', out_nodata )
                #-----------------------------------------
                w  = np.where(grid2 == out_nodata)
                nw = w[0].size
                print( 'nodata count =', nw)
                print()

            ds_in = None   # Close the tmp_file
            if (IN_MEMORY):
                gdal.Unlink( tmp_file )
#         else:
#             grid2 = np.zeros( (out_nrows, out_ncols), dtype='float32' )
#             grid2 += out_nodata
 
        #--------------------------------  
        # Write grid to new output file
        #------------------------------------
        # Example ISRIC filename:
        # CLYPPT_M_sl2_1km_South_Sudan.tiff
        #------------------------------------
        if not(BAD_BOX):       
            rtg_prefix = get_rtg_file_prefix( tiff_file, site_prefix )
            rtg_file   = out_dir + rtg_prefix + '.rtg'
            rtg_unit   = open( rtg_file, 'wb' )
            grid2 = np.float32( grid2 )
            grid2.tofile( rtg_unit )
            rtg_unit.close()
        n_grids += 1

    os.remove( tmp_file )   # Delete the temp tif file.
    
    print()
    print( 'Finished transforming ISRIC soil grid files.')
    print( '   Number of grids = ' + str(n_grids) )
    print( '   Number outside of model domain = ' + str(bad_box_count) )
    print()
    
#   transform_isric_soil_grid_files() 
#-------------------------------------------------------------------
def wosten_theta_s( C, S, OM, D, topsoil, FORCE_RANGE=True):

    #--------------------------------
    # From Wosten (1998). R^2 = 76%
    #-------------------------------------------------------
    # Note that theta_s should be in [0, 1], unitless
    #-------------------------------------------------------
    # Equations in Wosten (1998), Wosten et al. (2001) and
    # Matula et al. (2007) all agree.
    #-------------------------------------------------------
    p1 = 0.7919 + 0.001691*C - 0.29619*D
    p2 = -0.000001491*S**2 + 0.0000821*OM**2  
    p3 = 0.02427 * (1/C) + 0.01113 * (1/S) + 0.01472 * np.log(S)
    p4 = (-0.0000733 * OM * C) - (0.000619 * D * C)
    p5 = (-0.001183 * D * OM) - (0.0001664 * topsoil * S)

    theta_s = (p1 + p2 + p3 + p4 + p5)

    #--------------------------  
    # Propagate nodata values
    #--------------------------
    nd_values    = get_nodata_values()
    C_nodata     = nd_values['C']
    S_nodata     = nd_values['S']
    OM_nodata    = nd_values['OM']
    D_nodata     = nd_values['D']
    theta_nodata = nd_values['theta']
    w1 = np.logical_or( OM == OM_nodata, D == D_nodata)
    w2 = np.logical_or( C == C_nodata, S == S_nodata)
    theta_s[w1] = theta_nodata
    theta_s[w2] = theta_nodata
    
    #--------------------------------------
    # Check if grid values are reasonable
    #--------------------------------------
    w1 = np.logical_or( theta_s < 0.0, theta_s > 1.0 )
    n1 = w1.sum()
    if (n1 > 0):
        qmin = theta_s.min()
        qmax = theta_s.max()
        print('WARNING in wosten_theta_s:')
        print('   Some values are not in [0, 1].')
        if (qmin < 0.0):
            print('   min(theta_s) = ' + str(qmin) )
            print('   Possible nodata value.' )
        if (qmax > 1.0):
            print('   max(theta_s) = ' + str(qmax) )

        #-------------------------------
        # Option to replace bad values
        #-------------------------------
        if (FORCE_RANGE):
            print('   Forcing bad values into range.')
            wneg = (theta_s < 0.0)
            wpos = (theta_s > 0.0)
            theta_s[ wneg ] = theta_s[ wpos ].min()
            w1 = (theta_s > 1.0)
            w2 = (theta_s < 1.0)
            theta_s[ w1 ] = theta_s[ w2 ].max()
        print()

    return theta_s
    
#   wosten_theta_s()
#-------------------------------------------------------------------
def wosten_K_s( C, S, OM, D, topsoil ):

    #-------------------------------------   
    # From Wosten (1998). R^2 = 19%.
    # K_s^* = ln(K_s), K_s > 0.
    # Units in Wosten (1998) are cm/day.
    # Note that K_s should be > 0.
    #-------------------------------------
    
    ######################################################
    # NOTE!!  In term p5, the coefficient is given by:
    #         0.02986 in Wosten (1998) and:
    #         0.2986 in Wosten et al. (2001).
    ######################################################
    p1 = 7.755 + 0.0352*S + 0.93*topsoil
    p2 = -0.967*D**2 - 0.000484*C**2 - 0.000322*S**2
    p3 = (0.001 / S) - (0.0748 / OM) - 0.643*np.log(S)
    p4 = (-0.01398 * D * C) - (0.1673 * D * OM)
    #----------------
    # Wosten (1998)
    #----------------
    p5 = (0.02986 * topsoil * C) - 0.03305 * topsoil * S
    #-----------------------    
    # Wosten et al. (2001)
    #-----------------------
    ### p5 = (0.2986 * topsoil * C) - 0.03305 * topsoil * S
    
    Ks =  np.exp( p1 + p2 + p3 + p4 + p5 )  # [cm /day]

    #--------------------------------------
    # Check if grid values are reasonable
    #----------------------------------------------------------------
    # A very high value of K_s is 10 [cm / s] or 0.1 [m / s] which
    # equals 864,000 [cm / day] or 8.64 [km / day].  K_s is usually
    # much smaller and ranges over many orders of magnitude.
    #----------------------------------------------------------------
    Ks_max = 864000.0
    w1 = np.logical_or( Ks < 0.0, Ks > Ks_max )
    n1 = w1.sum()
    if (n1 > 0):
        Ksmin = Ks.min()
        Ksmax = Ks.max()
        print('ERROR in wosten_K_s:')
        print('   Some values are out of range.')
        if (Ksmin < 0.0):
            print('   min(K_s) = ' + str(Ksmin) )
        if (Ksmax > Ks_max):
            print('   max(K_s) = ' + str(Ksmax) )
        print()
                    
    #--------------------------------------------
    # Convert units from cm/day to m/sec for TF
    #--------------------------------------------
    #  cm       1 m         1 day         1 hr
    #  ---  *  -------  *  -------  *  ----------
    #  day     100 cm       24 hr       3600 sec
    #----------------------------------------------
    Ks = Ks / (100 * 24.0 * 3600.0)   # [m/sec]
      
    return Ks
    
#   wosten_K_s()
#-------------------------------------------------------------------
def wosten_alpha( C, S, OM, D, topsoil, FORCE_RANGE=True ):

    #---------------------------------   
    # From Wosten (1998). R^2 = 20%.
    # a^* = ln(a), a > 0.
    # Units are [1 / cm] ??    ###############  CHECK
    #---------------------------------
    p1 = -14.96 + 0.03135*C + 0.0351*S + 0.646*OM + 15.29*D
    p2 = -0.192*topsoil - 4.671*D**2 - 0.000781*C**2 - 0.00687*OM**2
    #----------------------------------------------------------
    # Wosten (1998), Nemes et al. (2001), Matula et al. (2007)
    #----------------------------------------------------------
    p3 = (0.0449 / OM) + 0.0663*np.log(S) + 0.1482*np.log(OM)
    p4 = (-0.04546 * D * S) - (0.4852 * D * OM)
    #-----------------------    
    # Wosten et al. (2001)
    #-----------------------
    ## p3 = (0.449 / OM) + 0.0663*np.log(S) + 0.1482*np.log(OM)
    ## p4 = (-0.4546 * D * S) - (0.4852 * D * OM)
    #-------------------------------
    p5 = 0.00673 * topsoil * C

    #---------------------------------------------------------   
    # Wosten formula returns |alpha|, so multiply by -1.
    # The pressure head is negative in the unsaturated zone,
    # zero at the water table and positive below that.
    # alpha < 0 is also needed so that psi_B = 1/alpha < 0
    # and G > 0.
    #---------------------------------------------------------  
    alpha = -1.0 * np.exp( p1 + p2 + p3 + p4 + p5)
    
    #--------------------------------------
    # Check if grid values are reasonable
    #------------------------------------------------------------
    # According to van Genuchten et al. (1991, p. 50, Table 3),
    # average alpha values range from:
    # -0.027 (clay) to -0.138 (sand) [1/cm].
    # If psi_B = 1/alpha, then psi_B in (-37.04, -7.25) [cm].
    #------------------------------------------------------------
    # According to van Genuchten et al. (1991, p. 51, Table 4),
    # average alpha values range from:
    # -0.005 (silty clay) to -0.145 (sand) [1/cm].
    # If psi_B = 1/alpha, then psi_B in (-200.0, -6.90) [cm].
    #------------------------------------------------------------
    # According to Smith (2001) book, psi_B values range from
    # about -90 to -9 cm.
    # If alpha = 1 / psi_B, alpha in [-1/9, -1/90] =
    # [-0.1111, -0.01111] or roughly [-0.1, -0.01].
    #-----------------------------------------------------
    ## alpha_min = -1.0 / 9.0
    ## alpha_max = -1.0 / 90.0
    alpha_min = -0.15
    alpha_max = -0.004
    w1 = np.logical_or( alpha < alpha_min, alpha > alpha_max )
    n1 = w1.sum()
    if (n1 > 0):
        amin  = alpha.min()
        amax1 = alpha.max()
        amax2 = alpha[ alpha != 0 ].max()
        print('WARNING in wosten_alpha:')
        print('   Some values are out of typical range.')
        print('   Typical range: -0.15 to -0.004 [1/cm]')
        if (amin < alpha_min):
            print('   min(alpha) = ' + str(amin) )
        if (amax2 > alpha_max):
            print('   max(alpha) = ' + str(amax2) + ' (excl. zero)')
        if (amax1 == 0):
            print('   max(alpha) = ' + str(amax1) + ' (incl. zero)')

        #-------------------------------
        # Option to replace bad values
        #-------------------------------
        if (FORCE_RANGE):
            print('   Forcing bad values into range.')
            wneg = (alpha < 0)
            wbig = (alpha >= 0)
            alpha[ wbig ] = alpha[ wneg ].max()
            #------------------------------------
            wbad = (alpha < -0.5)   #######################
            w2   = np.invert( wbad )
            alpha[ wbad ] = alpha[ w2 ].max()
        print()

    #--------------
    # For testing
    #--------------
#     amin  = alpha.min()
#     amax  = alpha.max()
#     amin2 = np.min( alpha[ alpha != -np.inf ] )
#     ## ananmin = np.nanmin( alpha )
#     print('In wosten_alpha():')
#     print('   min(alpha)    = ' + str(amin) )
#     print('   min2(alpha)   = ' + str(amin2) ) 
#     ## print('   nanmin(alpha) = ' + str(ananmin) )
#     print('   max(alpha)    = ' + str(amax) )
#     print()
        
    #-----------------------------------------------       
    # Convert units from [1/cm] to [1/m].
    #    1     100 cm                 1
    # X ---- x -------  = (X * 100)  ---
    #    cm       m                   m
    #-------------------------------------------------
    # Note pressure heads, psi, considered by Wosten
    # are in the range:  0 to -16,000 cm.
    # Permanent wilting point = -15,000 cm.
    # Hygroscopic condition   = -31,000 cm.
    # Field capacity is at -340 cm
    #   (at which water can be held against gravity)
    #-------------------------------------------------  
    alpha *= 100.0  # [1/m]
    
    return alpha
    
#   wosten_alpha()
#-------------------------------------------------------------------
def wosten_n( C, S, OM, D, topsoil, FORCE_RANGE=True):

    #------------------------------------------   
    # From Wosten (1998). R^2 = 54%.
    # n^* = ln(n-1), n > 1, unitless.
    # Wosten (1998) assumes that m = 1 - 1/n.
    #-------------------------------------------------------------
    # Equations in Wosten (1998) and Wosten et al. (2001) agree.
    #-------------------------------------------------------------
    # In Matula et al. (2007), p2 has: (0.002885 * OM)
    #-------------------------------------------------------------
    p1 = -25.23 - 0.02195*C + 0.0074*S - 0.1940*OM + 45.5*D
    p2 = -7.24*D**2 + 0.0003658*C**2 + 0.002885*OM**2 - (12.81 / D)
    p3 = (-0.1524 / S) - (0.01958 / OM) - 0.2876 * np.log(S)
    p4 = (-0.0709 * np.log(OM)) - (44.6 * np.log(D))
    p5 = (-0.02264 * D * C) + (0.0896 * D * OM) + (0.00718 * topsoil * C)
    
    n = (np.exp(p1 + p2 + p3 + p4 + p5) + 1)

    #--------------------------------------
    # Check if grid values are reasonable
    #--------------------------------------
    ### w1 = np.logical_or( n < 1.0, n > 3.0 )
    w1 = np.logical_or( n <= 1.0, n > 3.0 )     ### (Use <= vs. <)
    n1 = w1.sum()
    if (n1 > 0):
        nmin = n.min()
        nmax = n.max()
        print('ERROR in wosten_n:')
        print('   Some values are out of range.')
        print('   Typical range: 1.0 to 3.0 [unitless].')
        if (nmin < 1.0):
            print('   min(n) = ' + str(nmin) )
        if (nmax > 3.0):
            print('   max(n) = ' + str(nmax) )
        #-------------------------------
        # Option to replace bad values
        #-------------------------------
        if (FORCE_RANGE):
            print('   Forcing bad values into range.')
            w1 = (n <= 1.0)
            n[ w1 ] = 1.001
            ## w2 = np.invert( w1 )
            ## n[ w1 ] = n[ w2 ].min()
            #------------------------------------
            w3 = (n > 3.0)   ########### MAYBE 2.0 ??
            n[ w3 ] = 1.3    ###########
            ## w4 = np.invert( w3 )
            ## n[ w3 ] = n[ w4 ].max()
        print()
                       
    return n

#   wosten_n()
#-------------------------------------------------------------------
def wosten_L( C, S, OM, D, topsoil ):

    #-------------------------------------------------------    
    # From Wosten (1998). R^2 = 12%. 
    # L^* = ln[(L+10)/(10-L)], -10 < L < +10, unitless.
    # Mualem (1976) says L should be about 0.5 on average.
    # Do modern van Genuchten formulas assume L = 1/2 ??
    # See:  Wosten et al. (2001) for more about L.
    #-------------------------------------------------------------
    # Equations in Wosten (1998) and Wosten et al. (2001) agree.
    #-------------------------------------------------------------
    p1 = 0.0202 + 0.0006193*C**2  - 0.001136*OM**2
    p2 = -0.2316 * np.log(OM) - (0.03544 * D * C)
    p3 = (0.00283 * D * S) + (0.0488 * D * OM)
    
    s1 = (p1 + p2 + p3)
    L = 10 * (np.exp(s1) - 1)/(np.exp(s1) + 1)

    #--------------------------------------
    # Check if grid values are reasonable
    #-----------------------------------------------------
    # Wosten constrains L to be in [-10, 10]
    #-----------------------------------------------------
    w1 = np.logical_or( L < -10, L > 10 )
    n1 = w1.sum()
    if (n1 > 0):
        Lmin = L.min()
        Lmax = L.max()
        print('ERROR in wosten_L:')
        print('   Some values are out of range.')
        print('   Typical range: -10 to 10 [unitless].')
        if (L < -10.0):
            print('   min(L) = ' + str(Lmin) )
        if (L > 10.0):
            print('   max(L) = ' + str(Lmax) )
        print()
            
    return L
        
#   wosten_L()
#-------------------------------------------------------------------
def get_wosten_vars(C, S, OM, D, topsoil ):

    #----------------------------------------------------------
    # Use the Wosten (1998) pedotransfer functions to compute
    # theta_s, K_s, and van Genuchten parameters, then save
    # them to files.
    #----------------------------------------------------------
    theta_s = wosten_theta_s( C, S, OM, D, topsoil )
    K_s     = wosten_K_s( C, S, OM, D, topsoil )
    alpha   = wosten_alpha( C, S, OM, D, topsoil )
    n       = wosten_n( C, S, OM, D, topsoil )
    L       = wosten_L( C, S, OM, D, topsoil )
    
    return (theta_s, K_s, alpha, n, L)
    
#   get_wosten_vars()
#-------------------------------------------------------------------
def capillary_length_G(c, eta, n, inv_alpha, TBC=True):

    #------------------------------------------
    # Compute the Green-Ampt parameter, G > 0
    #------------------------------------------
    # G = capillary length scale > 0
    # psi_B = bubbling pressure head < 0
    # Both have same units of:  ?????????
    #----------------------------------------------------------
    # This formula is for transitional Brooks-Corey, given in
    # Peckham et al. (2017) GPF paper.
    #----------------------------------------------------------
    # According to Table 8.1 (p. 136) in R.E. Smith (2002) book,
    # G ranges from about 82 mm for sand to 2230 mm for clay.
    #-------------------------------------------------------------
    if (TBC):
        psi_B = inv_alpha
        G = -psi_B * gamma(1 + 1/c) * gamma((eta-1)/c) / gamma(eta/c)

    #------------------------------------------
    # Compute the Green-Ampt parameter, G > 0
    #------------------------------------------
    # G = capillary length scale > 0
    #----------------------------------------------------------
    # (2019-10-22) Derived this formula to compute G directly
    # from van Genuchten m and alpha parameters.  Use the
    # definition of G from Smith, change variables to:
    # Z = [(alpha*psi)^n + 1], and then integrate symbolically
    # with Mathematica.
    #---------------------------------------------------------
    # Limit of G as m -> 0 = 0.
    # Limit of G as m -> 1 = -1/alpha.
    #-----------------------------------
    if not(TBC):
        #------------------------------------------------
        # Limit of G[m] as m -> 0 equals 0.
        # Note: gamma(0) = inf, gamma(-1) = complex inf
        # So treat the case of m == 0 separately.
        #------------------------------------------------
        m  = 1 - (1.0/n)
        G  = np.zeros( m.shape, dtype=m.dtype)
        w1 = np.logical_and(m != 0, m != 1)
        m1 = m[w1]
        #----------------------------------------
#         mmin = m.min()
#         mmax = m.max()
#         print('mmin = ' + str(mmin))
#         print('mmax = ' + str(mmax))
        #----------------------------------------
        coeff = -inv_alpha[w1] * (1-m1)   # (> 0)
        term1 = 4 / (2 - 3*m1)
        term2 = gamma(1-m1)/gamma(m1/2)  + gamma(1+m1)/gamma(5*m1/2)
        term3 = gamma((3*m1/2) - 1)
        G[w1] = coeff * (term1 + (term2 * term3))
        #--------------------------------------------
        ## w2 = (m == 0)
        ## G[w2] = 0.0   # (already the case)
        #--------------------------------------------
        w3 = (m == 1)
        G[w3] = -inv_alpha[ w3 ]
    
    #--------------------------------------
    # Check if grid values are reasonable
    #--------------------------------------------------
    # See Peckham_et_al (2017) GPF paper on TopoFlow.
    #--------------------------------------------------
    # Theoretical results:
    #  G[c, eta] is in [0, -2 * psi_B].
    #  Definitions imply:  eta > 2, c > 3.
    #  Ignoring the coefficient in G[c,eta] above:
    #     Limit[ g[c,eta], c-> Infinity] = (eta/(eta-1)).
    #     For eta =2, this equals 2.
    #-------------------------------------------------------------
    # According to Table 8.1 (p. 136) in R.E. Smith (2002) book,
    # G ranges from about 82 mm for sand to 2230 mm for clay.
    #-------------------------------------------------------------
    G_min = 0.08  # (sand, meters)
    G_max = 2.3   # (clay, meters)
    w1 = np.logical_or( G < G_min, G > G_max )
    ## w1 = np.logical_or( G < 0, G > -2*psi_B )  ### (psi_B is grid)
    n1 = w1.sum()
    if (n1 > 0):
        Gmin = G.min()
        Gmax = G.max()
        # Gmin = np.nanmin(G)
        # Gmax = np.nanmax(G)
        ## Gmin = np.nanmin(G[ np.isfinite(G) ])
        ## Gmax = np.nanmax(G[ np.isfinite(G) ])

        print('WARNING in get_tBC_from_vG_vars:')
        print('   Some values in G grid are out of range.')
        print('   Typical range: 0.08 to 2.3 [m].' )
        if (Gmin < G_min):
            print('   min(G) = ' + str(Gmin) + ' [m]' )
        if (Gmax > G_max):
            print('   max(G) = ' + str(Gmax) + ' [m]' )
        print()
     
    return G
   
#   capillary_length_G()
#-------------------------------------------------------------------
def get_tBC_from_vG_vars( alpha, n, L ):

    #----------------------------------------------------------
    # Note:  This function is based on results from the book:
    # Smith, R.E. (2002) Infiltration Theory for Hydrologic
    # Applications, AGU Water Resources Monograph 15, 212 pp.
    # Smith first discusses the well-known Brooks-Corey and
    # van Genuchten models, and then proposes what he calls
    # the "transitional Brooks-Corey" model that is used by
    # TopoFlow.
    #----------------------------------------------------------
    # Convert van Genuchten parameters to those of the
    # transitional Brooks-Corey model.  Although K
    # is computed quite differently for the transitional
    # Brooks-Corey and van Genuchten methods, the equations
    # for ψ are the same if we use the formulas below.
    # For more information, see:
    # https://csdms.colorado.edu/wiki/
    #   Model_help:TopoFlow-Soil_Properties_Page   
    #---------------------------------------------------------
    # NOTE: L is often fixed at 1/2, but Wosten lets it vary.
    # See p. 51 in Wosten (1998).
    # Also see p. 1 in Schaap and van Genuchten (2005).
    #---------------------------------------------------------
    # These equations appear to be general:
    #---------------------------------------------------------
    # (1)  m      = (1 - 1/n)     (from Mualem (1976))
    # (1') n      = 1 / (1 - m)
    # (1") n-1    = m / (1 - m)
    # (2)  eta    = 2 + (3 * lambda)
    #---------------------------------------------------------
    # These equations come from forcing the transitional
    # Brooks-Corey and van Genuchten equations for pressure
    # head (psi) to match exactly (eta is not involved):
    # tBC params = psi_B, c, lambda  (and eta)
    # vG  params = alpha, m, n, L
    # See R.E. Smith (2002, ch. 2, p. 22)
    #---------------------------------------------------------
    # NOTE:  We only need alpha and n (not L) to set the
    #        transitional Brooks-Corey parameters.  We
    #        cannot make the functional forms match for K.
    #---------------------------------------------------------
    # (3) psi_B    = 1 / alpha_g   (both < 0)
    # (4) c        = n             (both > 1)
    # (5) lambda   = m * c = m * n    (Note: c/lambda = 1/m)
    #
    # Using (4) and (5) and (1):
    # (6) lambda = m * c = (1 - 1/n)*c = (1 - 1/n)*n = (n-1)
    #
    # Using (6) and (2)
    # eta    = 2 + (3*lambda) = 2 + 3*(n-1) = 3*n - 1
    # eta    = 3/(1-m) - 1 = [3 - (3-m)]/(1-m) = (2+m)/(1-m)
    #---------------------------------------------------------
    # (n > 1) => 0 < m < 1  
    # (n > 1) => (lambda > 0)
    # (n > 1) => (eta > 2)
    #---------------------------------------------------------
    w1 = (alpha != 0)
    w2 = np.invert( w1 )
    inv_alpha = np.zeros( alpha.shape, alpha.dtype )
    inv_alpha[w1] = (1.0 / alpha[w1])    # (Do this only once here.)
    inv_alpha[w2] = -9999.0
    #--------------------------
    psi_B  = inv_alpha
    c      = n
    lam    = (n - 1)
    eta    = 3*n - 1
    
    TBC = False
    # TBC = True
    G = capillary_length_G(c, eta, n, inv_alpha, TBC=TBC)
           
    return ( psi_B, c, lam, eta, G )

#   get_tBC_from_vG_vars()
#-------------------------------------------------------------------
def save_soil_hydraulic_vars( site_prefix=None,
         in_dir=None, out_dir=None,
         out_bounds=None, REPORT=True,
         out_xres_sec=None, out_yres_sec=None,
         RESAMPLE_ALGO='bilinear'):

    os.chdir( in_dir )
    if (out_dir[-1] != '/'):
        out_dir = out_dir + '/'
 
    #----------------------------------------------------  
    # Apply transformations to all soil layers (*.tiff)
    #----------------------------------------------------
#     regrid_isric_soil_grid_files( site_prefix, 
#            in_dir=in_dir, out_dir=out_dir,
#            out_bounds=out_bounds,
#            out_xres_sec=out_xres_sec, out_yres_sec=out_yres_sec,
#            RESAMPLE_ALGO='bilinear', REPORT=REPORT)
    #----------------------------------------------------
    # Almost ready;  still need to create an RTI file
    #----------------------------------------------------               
    transform_isric_soil_grid_files( site_prefix,
          in_dir=in_dir, out_dir=out_dir,
           out_bounds=out_bounds,
           out_xres_sec=out_xres_sec, out_yres_sec=out_yres_sec,
           RESAMPLE_ALGO='bilinear', REPORT=REPORT)

    #-------------------------------------------------------     
    # Compute soil hydraulic variable grids for all layers
    #-------------------------------------------------------              
    for k in range(7):
        layer = k + 1
        
        if (REPORT):
            print('===============================================')
            print('Computing variables for soil layer = ' + str(layer) )
            print('===============================================')
              
        (C, S, OM, D) = read_soil_grid_files( layer=layer,
                             res_str='1km',   ################
                             directory=out_dir, REPORT=REPORT,
                             site_prefix=site_prefix )

        if (layer == 1):
            topsoil = 1
        else:
            topsoil = 0
 
        #----------------------------------------------------------   
        # Use the Wosten PTF equations to compute key
        # soil hydraulic variables, including van Genuchten vars.
        #----------------------------------------------------------
        (theta_s, K_s, alpha, n, L) = get_wosten_vars( C, S, OM, D, topsoil )
        
        #----------------------------------------------------------   
        # Compute "transitional Brooks-Corey" soil hydraulic vars
        # from the van Genuchten vars.
        #----------------------------------------------------------   
        (psi_B, c, lam, eta, G ) = get_tBC_from_vG_vars( alpha, n, L )

        #------------------------------------------
        # Get consistent estimates of:
        # theta_i = qi = initial soil moisture, &
        # theta_r = qr = residual soil moisture
        # K_i = initial hydraulic conductivity
        #------------------------------------------
        # NOT FINISHED YET.
     
        prefix = site_prefix + '_sl' + str(layer)
        ## prefix = case_prefix+ '_sl' + str(layer) 
        prefix = out_dir + prefix   ##################

        #---------------------------------   
        # Basic soil hydraulic variables
        #---------------------------------
        Ks_file   = prefix + '_2D-Ks.bin'
        qs_file   = prefix + '_2D-qs.bin'
        pB_file   = prefix + '_2D-pB.bin'
        #-------------------------------------------
        # The transitional Brooks-Corey parameters
        #-------------------------------------------
        c_file    = prefix + '_2D-c.bin'    
        lam_file  = prefix + '_2D-lam.bin'
        G_file    = prefix + '_2D-G.bin'
        ## eta_file  = prefix + '_2D-eta.bin'
    #     c_file    = prefix + '_2D-tBC-c.bin'    
    #     lam_file  = prefix + '_2D-tBC-lam.bin'
    #     G_file    = prefix + '_2D-tBC-G.bin'
    #     ## eta_file  = prefix + '_2D-tBC-eta.bin'
        #-------------------------------
        # The van Genuchten parameters
        #-------------------------------
        a_file    = prefix + '_2D-vG-alpha.bin'
        n_file    = prefix + '_2D-vG-n.bin'
        L_file    = prefix + '_2D-vG-L.bin'

        #------------------------------------------
        # Consistent estimates of:
        # theta_i = qi = initial soil moisture, &
        # theta_r = qr = residual soil moisture
        # K_i = initial hydraulic conductivity
        #------------------------------------------
        # qi_file = prefix + '_2D-qi.bin'
        # qr_file = prefix + '_2D-qr.bin'
        # Ki_file = prefix + '_2D-Ki.bin'
    
        #-------------------------------    
        # Write all variables to files
        #-------------------------------
        Ks_unit = open(Ks_file, 'wb')
        K_s = np.float32( K_s )  
        K_s.tofile( Ks_unit )
        Ks_unit.close()
        #----------------------------------
        qs_unit = open(qs_file, 'wb')
        theta_s = np.float32( theta_s )
        theta_s.tofile( qs_unit )
        qs_unit.close()
        #----------------------------------
        pB_unit = open(pB_file, 'wb')
        psi_B   = np.float32( psi_B )
        psi_B.tofile( pB_unit)
        pB_unit.close()
        #----------------------------------
        c_unit = open(c_file, 'wb')
        c = np.float32( c )
        c.tofile( c_unit )
        c_unit.close()
        #----------------------------------
        lam_unit = open(lam_file, 'wb')
        lam = np.float32( lam )
        lam.tofile( lam_unit)
        lam_unit.close()
        #----------------------------------
        G_unit = open(G_file, 'wb')
        G = np.float32( G )
        G.tofile( G_unit )
        G_unit.close()
        #----------------------------------
    #     eta_unit = open(eta_file, 'wb')
    #     eta = np.float32( eta )
    #     eta.tofile( eta_unit )
    #     eta_unit.close()
        #----------------------------------
        a_unit = open(a_file, 'wb')
        alpha = np.float32( alpha )
        alpha.tofile( a_unit )
        a_unit.close()
        #----------------------------------
        n_unit = open(n_file, 'wb')
        n = np.float32( n )
        n.tofile( n_unit)
        n_unit.close()
        #----------------------------------
        L_unit = open(L_file, 'wb')
        L = np.float32( L )
        L.tofile( L_unit )
        L_unit.close()
 
    print('Finished computing & saving soil hydraulic vars.')
    print()

#     if (REPORT):
#         print('Finished computing & saving soil hydraulic vars.')
#         print()
           
#   save_soil_hydraulic_vars()   
#-------------------------------------------------------------------
# def save_soil_hydraulic_vars0( site_prefix=None,
#                               in_dir=None, out_dir=None,
#                               layer=1, BARO_1MIN_TEST=False,
#                               LOL_1MIN_TEST=False):
# 
#     if (BARO_1MIN_TEST):
#         out_dir = '/Users/peckhams/ISRIC_Files/Ethiopia/'
#         in_dir  = out_dir + '_isric_soil_data/'    
#         site_prefix = 'Baro_Gam_1min'
#     if (LOL_1MIN_TEST):
#         out_dir = '/Users/peckhams/ISRIC_Files/South_Sudan/'
#         in_dir  = out_dir + '_isric_soil_data/' 
#         site_prefix = 'Lol-Kuru_1min'
#     if (directory is None):
#         in_dir = '/Users/peckhams/ISRIC_Files/Test/'
#         out_dir = '/Users/peckhams/ISRIC_Files/Test/'
#     os.chdir( directory )
#    
#     print('===============================================')
#     print('Computing variables for soil layer = ' + str(layer) )
#     print('===============================================')
#         
#     if (layer == 1): 
#         #----------------------------------------------------  
#         # Apply transformations to all soil layers (*.tiff)
#         #----------------------------------------------------        
#         transform_isric_soil_grid_files( directory=directory,
#                   site_prefix=site_prefix,
#                   BARO_1MIN_TEST=BARO_1MIN_TEST,
#                   LOL_1MIN_TEST=LOL_1MIN_TEST)
#               
#     (C, S, OM, D) = read_soil_grid_files( layer=layer, directory=directory,
#                                           site_prefix=site_prefix )
# 
#     if (layer == 1):
#         topsoil = 1
#     else:
#         topsoil = 0
#     
#     (theta_s, K_s, alpha, n, L) = get_wosten_vars( C, S, OM, D, topsoil )   
#     (psi_B, c, lam, eta, G ) = get_tBC_from_vG_vars( alpha, n, L )
# 
#     #------------------------------------------
#     # Get consistent estimates of:
#     # theta_i = qi = initial soil moisture, &
#     # theta_r = qr = residual soil moisture
#     # K_i = initial hydraulic conductivity
#     #------------------------------------------
#     # NOT FINISHED YET.
#      
#     prefix = site_prefix + '_sl' + str(layer)
#     ## prefix = case_prefix+ '_sl' + str(layer) 
#      
#     #---------------------------------   
#     # Basic soil hydraulic variables
#     #---------------------------------
#     Ks_file   = prefix + '_2D-Ks.bin'
#     qs_file   = prefix + '_2D-qs.bin'
#     pB_file   = prefix + '_2D-pB.bin'
#     #-------------------------------------------
#     # The transitional Brooks-Corey parameters
#     #-------------------------------------------
#     c_file    = prefix + '_2D-c.bin'    
#     lam_file  = prefix + '_2D-lam.bin'
#     G_file    = prefix + '_2D-G.bin'
#     ## eta_file  = prefix + '_2D-eta.bin'
# #     c_file    = prefix + '_2D-tBC-c.bin'    
# #     lam_file  = prefix + '_2D-tBC-lam.bin'
# #     G_file    = prefix + '_2D-tBC-G.bin'
# #     ## eta_file  = prefix + '_2D-tBC-eta.bin'
#     #-------------------------------
#     # The van Genuchten parameters
#     #-------------------------------
#     a_file    = prefix + '_2D-vG-alpha.bin'
#     n_file    = prefix + '_2D-vG-n.bin'
#     L_file    = prefix + '_2D-vG-L.bin'
# 
#     #------------------------------------------
#     # Consistent estimates of:
#     # theta_i = qi = initial soil moisture, &
#     # theta_r = qr = residual soil moisture
#     # K_i = initial hydraulic conductivity
#     #------------------------------------------
#     # qi_file = prefix + '_2D-qi.bin'
#     # qr_file = prefix + '_2D-qr.bin'
#     # Ki_file = prefix + '_2D-Ki.bin'
#     
#     #-------------------------------    
#     # Write all variables to files
#     #-------------------------------
#     Ks_unit = open(Ks_file, 'wb')
#     K_s = np.float32( K_s )  
#     K_s.tofile( Ks_unit )
#     Ks_unit.close()
#     #----------------------------------
#     qs_unit = open(qs_file, 'wb')
#     theta_s = np.float32( theta_s )
#     theta_s.tofile( qs_unit )
#     qs_unit.close()
#     #----------------------------------
#     pB_unit = open(pB_file, 'wb')
#     psi_B   = np.float32( psi_B )
#     psi_B.tofile( pB_unit)
#     pB_unit.close()
#     #----------------------------------
#     c_unit = open(c_file, 'wb')
#     c = np.float32( c )
#     c.tofile( c_unit )
#     c_unit.close()
#     #----------------------------------
#     lam_unit = open(lam_file, 'wb')
#     lam = np.float32( lam )
#     lam.tofile( lam_unit)
#     lam_unit.close()
#     #----------------------------------
#     G_unit = open(G_file, 'wb')
#     G = np.float32( G )
#     G.tofile( G_unit )
#     G_unit.close()
#     #----------------------------------
# #     eta_unit = open(eta_file, 'wb')
# #     eta = np.float32( eta )
# #     eta.tofile( eta_unit )
# #     eta_unit.close()
#     #----------------------------------
#     a_unit = open(a_file, 'wb')
#     alpha = np.float32( alpha )
#     alpha.tofile( a_unit )
#     a_unit.close()
#     #----------------------------------
#     n_unit = open(n_file, 'wb')
#     n = np.float32( n )
#     n.tofile( n_unit)
#     n_unit.close()
#     #----------------------------------
#     L_unit = open(L_file, 'wb')
#     L = np.float32( L )
#     L.tofile( L_unit )
#     L_unit.close()
#         
# #   save_soil_hydraulic_vars0()   
#-------------------------------------------------------------------

