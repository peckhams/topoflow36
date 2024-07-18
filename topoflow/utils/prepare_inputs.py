#!/usr/bin/env python
# coding: utf-8
 
# This module uses a collection of TopoFlow utilities to create
# required input files for the TopoFlow hydrologic model.
# Before you can run this code, you will need to install the
# TopoFlow 3.6 Python package. This notebook uses numerous
# utilities from the TopoFlow 3.6 package, in topoflow/utils.
# Detailed instructions and background information for how
# to install TopoFlow in a conda environment are given in
# a Jupyter notebook:
#  Appendix 1: Installing TopoFlow in a conda Environment.
# 
# This module addresses the important issue of portability.
# By simply specifying a geographic bounding box and a few other
# bits of information for any river basin (e.g. in Ethiopia) in
# the first few methods, you can create all of the input files that
# are required by TopoFlow -- including the space-time rainfall
# grids, if you have downloaded the source data files.  By changing
# the regional source files for the DEM, soils data and met data,
# you can use this module to create the input files required by
# TopoFlow to model any river basin on Earth.  This is because it
# uses data sets that are available globally, such as the
# [<b>MERIT DEM</b>](http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/)
# and the
# [<b>ISRIC SoilGrids</b>](https://www.isric.org/explore/soilgrids)
# database [and services](https://soilgrids.org/).
# 
# This module has been designed to create a set of TopoFlow input
# files for any river basin on Earth.  It assumes that you have the
# necessary source data files for a larger region that contains the
# basin of interest.  For the DEM, this could be a global DEM, or a
# DEM for a particular country, for example.  Source data files will
# be clipped to the river basin's geographic bounding box and also
# resampled to the desired spatial resolution.  The
# [<b>MERIT DEM</b>](http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/)
# tiles may be mosaicked and/or clipped to create the source DEM.
#
# It is also assumed that you have
# [<b>ISRIC SoilGrids</b>](https://soilgrids.org/)
# for a region that contains the basin of interest.  You must first
# download these files to your computer and then edit the method
# get_source_dirs() with the path.  The original SoilGrid data
# is from 2017 and the latest is from 2020.  The 2017 data is
#  available for the globe, aggregated to 1km spatial resolution,
# and is usually sufficient. That data is available here:
# [<b>ISRIC SoilGrids 2017 1km</b>]
# (https://files.isric.org/soilgrids/former/2017-03-10/aggregated/1km/).
# For TopoFlow, you need the files that begin with these 7 prefixes,
# for all soil layers:
# 
# BLDFIE = Bulk density [kg / m3]
# BDRICM = Absolute depth to bedrock [cm]
# BDTICM = (Total?) depth to bedrock [cm]
# CLYPPT = Mass fraction of clay [%]
# ORCDRC = Soil organic carbon content (fine earth fraction)  [g / kg]
# SLTPPT = Mass fraction of silt [%]
# SNDPPT = Mass fraction of sand [%]
# 
# The 2020 data is available here:
#     [<b>ISRIC SoilGrids 2020</b>](https://files.isric.org/soilgrids/latest/data/).
# You will also need to obtain sample CFG files for TopoFlow, as
# well as sample basin outlet files. (GIVE GitHub LINK HERE)  ###################
#
#----------------------------------------------------------------------------------
# Example use:
#
# Edit some info in this file like: paths & set_basin_info(), and then:
#
# >>> from topoflow.utils import prepare_inputs as prep
# >>> obj = prep.get_inputs()
# >>> obj.prepare_all_inputs()
#
#----------------------------------------------------------------------------------
#
#  class get_inputs()    ## Better name ???
#      __init__()
#      prepare_all_inputs()       # This one calls most of the others
#      get_topoflow_dirs()        #### NEW
#      get_source_dirs()
#      create_output_dirs()
#      set_output_filenames()     #### NEW
#      set_time_info()            #### NEW
#      set_basin_info()           ####  DEM_ncols, DEM_nrows, etc.
#      set_met_info()
#      -------------------------
#      copy_cfg_file_templates()
#      get_cfg_files()
#      copy_outlets_file()
#      get_outlets_file()         ##### Create good default
#      -------------------------
#      get_ngen_dem()
#      clip_dem_from_source()
#      read_dem_as_geotiff()
#      read_dem_as_rtg()
#      read_dem_as_netcdf()
#      write_dem_as_rtg()
#      get_profile_smoothed_dem()
#      replace_negative_vals_in_dem()
#      fill_pits_in_dem()
#      get_d8_flow_grid()
#      get_d8_area_grid()
#      get_d8_slope_grid()
#      get_d8_aspect_grid()
#      get_channel_width_grid()
#      get_manning_n_grid()
#      get_channel_sinuosity_grid()
#      get_bankfull_depth_grid()
#      get_init_channel_depth_grid()
#      ------------------------------
#      get_soil_hydraulic_grids()
#      ------------------------------
#      write_simple_rain_file()
#      get_rainfall_grid_stack()
#      subset_ngen_csv_file()
#
#----------------------------------------------------------------------------------

import numpy as np
import os, os.path
import glob, shutil

from topoflow.utils import cfg_templates as ct
from topoflow.utils import outlets  # (for write_outlet_file)
from topoflow.utils import regrid
from topoflow.utils import import_grid
from topoflow.utils import fill_pits
from topoflow.utils import rtg_files
from topoflow.utils import rti_files
from topoflow.utils import parameterize
from topoflow.utils import init_depth
from topoflow.utils import pedotransfer
from topoflow.utils import new_slopes

from topoflow.utils.ngen import hydrofab_utils as hfu

from topoflow.components import d8_global
from topoflow.components import smooth_DEM

#----------------------------------------------------------------------------------
class get_inputs():

    def __init__(self, SILENT=False):

        #----------------------------------------------------------   
        # Set various default options.  These can be changed
        # immediately after an instance of this class is created.
        #----------------------------------------------------------
        # If the NGEN_CSV option is set, self.test_dir should
        # be set to:  ngen_dir + 'data/topoflow/input_files/'
        #-------------------------------------------------------    
        self.home_dir    = os.path.expanduser("~") + os.sep
        self.test_dir    = self.home_dir + 'basins_TEST' + os.sep   ######
        self.output_dir  = self.home_dir + 'output' + os.sep
        self.ngen_dir    = self.home_dir + 'Dropbox/GitHub/ngen'
        self.src_ngen_met_dir = self.ngen_dir + 'data/topoflow/forcing/huc01/'
        ## self.src_ngen_met_dir = self.home_dir + 'TF_Data/NGEN_AORC/HUC01/csv/'        
        #-------------------------------------------
        self.site_prefix = 'Unknown'  ###########
        self.case_prefix = 'Test1'
        #----------------------------
        self.SILENT      = SILENT
        self.TEST_RAIN   = False
        #-------------------------------------------------------
        # Current options for preparing input files are:
        #   NGEN_CSV and EAST_AFRICA.
        # For EAST_AFRICA option, can also choose a rainfall
        #   product by setting CHIRPS, CHIRPS2, GLDAS, or GPM.
        #-------------------------------------------------------
        self.NGEN_CSV    = True
        self.EAST_AFRICA = False
        #-------------------------
        self.CHIRPS      = False
        self.CHIRPS2     = False
        self.GLDAS       = False
        self.GPM         = False

        #--------------------------------------------------------
        # For now, keep these the same for all the river basins
        # They are reasonable defaults.
        #--------------------------------------------------------
        self.channel_width_power = 0.5
        self.min_sinuosity   = 1.0    # (BY DEFINITION.  DO NOT CHANGE.)
        self.max_sinuosity   = 1.3    # [m/m]
        self.min_manning_n   = 0.03   # [m / s^(1/3)]
        self.max_manning_n   = 0.2    # [m / s^(1/3)]
        self.bankfull_depth_power = 0.4
        self.bank_angle = 30.0   # [degrees]
        self.DRAINS_TO_OCEAN = False   # (this is the default)

        #-----------------------------------------------------------   
        # These need to be set to reasonable values for each basin
        #-----------------------------------------------------------
        self.max_river_width = 10.0    # [meters]  from Google Maps or Google Earth
        self.A_out_km2       = 50.0    # total contributing area at basin outlet [km2]
        self.Qbase_out       = 1.0     # estimated baseflow discharge at basin outlet [m^3 / s]
        self.max_bankfull_depth = 2.0  # [meters]    # from literature or data

        #------------------------------------
        # Set the channel process time step
        #---------------------------------------------------------
        # The largest time step, dt, for which the model will be
        # numerically stable will depend on the grid cell size.
        # Larger dt values are possible with larger grid cells.
        # Making dt smaller than necessary increases run time.
        #---------------------------------------------------------
        self.chan_dt_str = '5.0'   # [secs]
        # self.chan_dt_str = '10.0'  # [secs]
        # self.chan_dt_str = '20.0'  # [secs]        
        # self.chan_dt_str = '30.0'  # [secs]
        # self.chan_dt_str = '100.0'  # [secs]
        # self.chan_dt_str = '200.0'  # [secs]  # (try larger values later)
        # chan_dt_str = '300.0'  # [secs]
        # if (out_xres_sec == 120.0):
        #     chan_dt_str = '900.0' # [secs]
        # elif (out_xres_sec == 60.0):
        #    chan_dt_str = '300.0' # [secs]
        # else:
        #    chan_dt_str = '300.0' # [secs]
        
    #   __init__()
    #-------------------------------------------------------------------------------
    def prepare_all_inputs(self, site_prefix='Tana_120sec',
                           case_prefix='Test1',
                           CLIP_DEM=True, NGEN_DEM=False,
                           NO_SOIL=True, NO_MET=True ):

        if not(self.NGEN_CSV) and not(self.EAST_AFRICA):
           print('ERROR in prepare_all_inputs:')
           print('To use this function, you must currently set')
           print('either the NGEN_CSV or EAST_AFRICA option to True.')
           print('This is needed by the set_time_info() function.')
           print()
           return
        self.site_prefix = site_prefix
        self.case_prefix = case_prefix
        ## self.input_directory = self.test_dir + self.site_prefix + '/'
        #-----------------------------
        self.NO_MET      = NO_MET
        self.NO_SOIL     = NO_SOIL
        #-----------------------------
        # self.get_topoflow_dirs()       
        self.get_source_dirs()         # calls get_topoflow_dirs()
        #-----------------------------
        self.set_basin_info()
        self.set_time_info()
        # Add info like DEM_ncols, DEM_nrows later.  ##############
        self.set_met_info()      ### Need here to set "met_dt_str" for met_cfg
        #------------------------------
        self.create_output_dirs()      # after set_time_info() & set_met_info()
        self.set_output_filenames()    # after set_basin_info()
        self.copy_cfg_file_templates()
        self.get_cfg_files()
        #------------------------------
        if (NGEN_DEM):
            self.get_ngen_dem( cat_id_str=self.site_prefix )  # 2024-07-05
        if (CLIP_DEM):
            # Otherwise, get DEM from S3 bucket/VRT
            # via multi_bmi.get_new_dem(), via hft.get_dem()      
            self.clip_dem_from_source()
        self.read_dem_as_geotiff()
        # self.read_dem_as_rtg()
        # self.read_dem_as_netcdf()
        self.write_dem_as_rtg()
        self.replace_negative_vals_in_dem()
        self.fill_pits_in_dem()
        #------------------------------ 
        self.get_d8_flow_grid()
        self.get_d8_area_grid()
        self.get_d8_slope_grid()
        self.get_d8_aspect_grid()
        self.get_outlets_file()       # After get_d8_flow_grid()
        #------------------------------
        self.get_channel_width_grid()
        self.get_manning_n_grid()
        self.get_channel_sinuosity_grid()
        self.get_bankfull_depth_grid()
        self.get_init_channel_depth_grid()
        #-------------------------------------
        self.write_simple_rain_file()
        if (self.NGEN_CSV):
            self.subset_ngen_csv_file()        #################
        if not(NO_SOIL):
            self.get_soil_hydraulic_grids()
        if not(NO_MET):
            self.get_rainfall_grid_stack()
        #--------------------------------------
        print('Finished preparing TopoFlow input files.')
        print()

    #   prepare_all_inputs()
    #-------------------------------------------------------------------------------
    def get_topoflow_dirs(self):
    
        #------------------------------------------------------
        # Get path to the current file (prepare_inputs.py),
        # as a way to determine the TopoFlow directories.
        # At top need: "#! /usr/bin/env python" ?
        # See: https://docs.python.org/2/library/os.path.html
        #------------------------------------------------------
        utils_dir = os.path.dirname( __file__ )
        tf_dir    = os.path.join( utils_dir, '..' )
        #--------------------------------------------------------
        self.utils_dir = os.path.abspath( utils_dir ) + os.sep
        self.tf_dir    = os.path.abspath( tf_dir )    + os.sep
        #------------------------------------------------------------------------
        # Include null string in join() call to add os.sep at the end.
        #------------------------------------------------------------------------
        self.components_dir   = os.path.join( self.tf_dir, 'components', '' )
        self.framework_dir    = os.path.join( self.tf_dir, 'framework', '' )
        self.cfg_template_dir = os.path.join( self.tf_dir, 'cfg_templates', 'Test1', '' )
        #----------------------------------------------------------------
        if not(self.SILENT):
           print('TopoFlow 3.6 package directory is:')
           print(self.tf_dir)
           print()
        #----------------------------------------------------------------
        if not(os.path.exists(self.cfg_template_dir)):
           print('ERROR:  TopoFlow CFG file template dir not found.')
           print('        Please upgrade to TopoFlow v3.65 or higher.')
           print()
            
    #   get_topoflow_dirs()
    #-------------------------------------------------------------------------------
    def get_source_dirs(self):
    
        """Get directories for regional DEM, soil data, cfg files and outlet file."""
        #----------------------------------------            
        # Change these directories as necessary      ############
        #----------------------------------------       
        self.data_dir = self.home_dir + 'TF_Data/'
        dd = self.data_dir

        #------------------------------------
        # Source dir for cfg template files
        #------------------------------------
        self.get_topoflow_dirs()  # will set self.cfg_template_dir

        #------------------------------
        # Source dir for outlet files
        #------------------------------
        self.src_outlet_dir = dd + '_Test1/outlet_files/'

        #--------------------------
        # Source dir for DEM data
        #--------------------------
        self.src_dem_dir  = self.home_dir + 'DEMs/Horn_of_Africa/'
        self.src_dem_file = 'Horn_of_Africa_MERIT_DEM.tif'      ########
               
        #---------------------------
        # Source dir for soil data
        #---------------------------
        self.src_soil_dir = dd + 'ISRIC_SoilGrids/__Global_2017_1km/'
                
        #---------------------------
        # Source dirs for met data
        #---------------------------
        self.src_chirps_dir   = dd + 'CHIRPS_Rain_Africa_6hr_2015-10_to_2018-10_onedir/'
        self.src_chirps2_dir  = dd + 'CHIRPS_Rain_Africa_6hr_2005-01_to_2015-01_onedir/'
        self.src_gldas_dir    = dd + 'GLDAS_Rain_Baro_2015-10_to_2018-10/'
        self.src_gpm_dir      = dd + 'GPM_Rain_Ethiopia_30min_2015-10_to_2018-10/'
        ### self.src_ngen_met_dir = dd + 'NGEN_AORC/HUC01/csv/'

        if not(self.SILENT):
            print( 'Home directory          =', self.home_dir )
            print( 'Source DEM directory    =', self.src_dem_dir )
            print( 'Source Soil directory   =', self.src_soil_dir )
            print( 'CFG template directory  =', self.cfg_template_dir )
            print( 'Source outlet directory =', self.src_outlet_dir )
            print()

    #   get_source_dirs()
    #-------------------------------------------------------------------------------
    def create_output_dirs(self):
 
        #----------------------------------------------   
        # Get the names of new directories.  home_dir
        # & test_dir are defined in __init__().
        #----------------------------------------------
        self.basin_dir = self.test_dir + self.site_prefix + '/'
        self.topo_dir  = self.basin_dir + '__topo/'
        self.soil_dir  = self.basin_dir + '__soil/'
        self.met_dir   = self.basin_dir + '__met/'
        self.misc_dir  = self.basin_dir + '__misc/'

        if (self.NO_MET):
            self.cfg_dir = self.basin_dir + self.case_prefix + '_cfg/'
        else:
            #------------------------------------------------
            # Include "date_range" & "rain_type" in cfg_dir
            #------------------------------------------------
            case_info_str = self.date_range_str + '_' + self.rain_type_str
            cp = self.case_prefix + '_'
            self.cfg_dir = self.basin_dir + cp + case_info_str + '_cfg/'

        #---------------------------------------        
        # Create the directories, if necessary
        #---------------------------------------       
        if not(os.path.exists( self.test_dir )):  os.mkdir( self.test_dir )
        if not(os.path.exists( self.basin_dir )): os.mkdir( self.basin_dir )
        if not(os.path.exists( self.topo_dir )):  os.mkdir( self.topo_dir )
        if not(os.path.exists( self.soil_dir )):  os.mkdir( self.soil_dir )
        if not(os.path.exists( self.met_dir )):   os.mkdir( self.met_dir )
        if not(os.path.exists( self.misc_dir )):  os.mkdir( self.misc_dir )
        if not(os.path.exists( self.cfg_dir )):   os.mkdir( self.cfg_dir )
            
        if not(self.SILENT):
            print( 'Home directory  =', self.home_dir )
            print( 'Basin directory =', self.basin_dir )
            print( 'Topo directory  =', self.topo_dir )
            print( 'Soil directory  =', self.soil_dir )
            print( 'Met directory   =', self.met_dir )
            print( 'Misc directory  =', self.misc_dir )
            print( 'CFG directory   =', self.cfg_dir )
            print()

    #   create_output_dirs()
    #-------------------------------------------------------------------------------
    def set_output_filenames(self):
    
        """Set all the output filenames"""
        sp = self.site_prefix
        cp = self.case_prefix
        #------------------------------------------
        # Set the "Topo" related output filenames
        #------------------------------------------
        self.DEM_tif_file   = self.topo_dir + sp + '_rawDEM.tif'
        self.DEM_file       = self.topo_dir + sp + '_rawDEM.rtg'
        self.new_DEM_file   = self.topo_dir + sp + '_DEM.rtg'
        self.rti_file       = self.topo_dir + sp + '.rti'
        self.d8_code_file   = self.topo_dir + sp + '_flow.rtg'
        self.d8_area_file   = self.topo_dir + sp + '_d8-area.rtg'
        self.d8_slope_file  = self.topo_dir + sp + '_slope.rtg'
        self.d8_aspect_file = self.topo_dir + sp + '_aspect.rtg'
        self.width_file     = self.topo_dir + sp + '_chan-w.rtg'
        self.manning_file   = self.topo_dir + sp + '_chan-n.rtg'
        self.sinu_file      = self.topo_dir + sp + '_sinu.rtg'
        self.dbank_file     = self.topo_dir + sp + '_d-bank.rtg'
        self.d0_file        = self.topo_dir + sp + '_d0.rtg'
                        
        #------------------------------------------
        # Set the "Soil" related output filenames
        #------------------------------------------ 
        # self.soil_file = self.soil_dir + sp + extension
         
        #-----------------------------------------
        # Set the "Met" related output filenames
        #-----------------------------------------
        self.test_rain_file = self.met_dir + cp + '_rain_rates.txt'
        # Note:  This is now set in: get_rainfall_grid_stack()
        # self.rain_rts_path  = self.met_dir + self.rain_rts_file
         
        if not(self.SILENT):
            print('DEM_tif_file   =', self.DEM_tif_file )
            print('DEM_file       =', self.DEM_file )
            print('new_DEM_file   =', self.new_DEM_file )
            print('rti_file       =', self.rti_file )
            print('d8_code_file   =', self.d8_code_file )
            print('d8_area_file   =', self.d8_area_file )
            print('d8_slope_file  =', self.d8_slope_file )
            print('d8_aspect_file =', self.d8_aspect_file )
            print('width_file     =', self.width_file )
            print('manning_file   =', self.manning_file )
            print('sinu_file      =', self.sinu_file )
            print('dbank_file     =', self.dbank_file )
            print('d0_file        =', self.d0_file )
            print()
                           
    #   set_output_filenames()
    #-------------------------------------------------------------------------------
    def set_time_info(self):
 
        """Set info for time period of interest (& met data)"""   
        #----------------------------------------------------------
        # Note that the number of minutes in a regular year
        #    = 60 * 24 * 365 = 525600,
        # and the number of minutes in a leap year
        #    = 60 * 24 * 366 = 527040.
        # The year 2016 is a leap year, so the number of minutes
        # in 2015-10-01 to 2018-10-01 = 2*(525600) + 527040
        #    = 1578240. (Note: 525600 * 3 = 1576800).
        #----------------------------------------------------------
        # Note: time steps are set elsewhere.
        #       met_dt_str is set in set_met_info().
        #----------------------------------------------------------
        if (self.EAST_AFRICA):
            if (self.CHIRPS or self.GLDAS or self.GPM):
                self.start_date_str  = '2015-10-01'
                self.end_date_str    = '2018-10-01'
                self.date_range_str  = '2015-10_to_2018-10'
                self.GMT_offset_str  = '3'    # Ethiopia and Kenya
                self.start_month_str = 'October'
                self.start_day_str   = '1'
                self.T_stop_str      = '1578240'  # [minutes]  See note above.
            elif (self.CHIRPS2):
                self.start_date_str  = '2005-01-01'
                self.end_date_str    = '2015-01-01'
                self.date_range_str  = '2005-01_to_2015-01'
                self.GMT_offset_str  = '3'    # Ethiopia and Kenya
                self.start_month_str = 'January'
                self.start_day_str   = '1'
                self.T_stop_str      = '5258880'  # [min] with 2 leap years: 2008, 2012
            else:
                print('## WARNING: Rainfall product not set (East Africa).')

        if (self.NGEN_CSV):
            # Wettest month for cat-2913 and cat-209 (headwater HUC01 catchments)
            self.start_date_str  = '2012-06-01'
            self.end_date_str    = '2012-07-01'
            self.date_range_str  = '2012-06_to_2012-07'
            self.GMT_offset_str  = '-5'    # Eastern Time  ############
            self.start_month_str = 'June'
            self.start_day_str   = '1'
            self.T_stop_str      = '43200'  # [min] (minutes in June 2012)

#         if (self.NGEN_CSV):
#             # Wettest month for cat-84 (headwater HUC01 catchment)
#             self.start_date_str  = '2011-05-01'
#             self.end_date_str    = '2011-06-01'
#             self.date_range_str  = '2011-05_to_2011-06'
#             self.GMT_offset_str  = '-5'    # Eastern Time  ############
#             self.start_month_str = 'May'
#             self.start_day_str   = '1'
#             self.T_stop_str      = '44640'  # [min] (minutes in May 2011)
                             
#         if (self.NGEN_CSV):
#             self.start_date_str  = '2007-01-01'
#             self.end_date_str    = '2007-02-01'
#             self.date_range_str  = '2007-01_to_2007-02'
#             self.GMT_offset_str  = '-5'    # Eastern Time  ############
#             self.start_month_str = 'January'
#             self.start_day_str   = '1'
#             self.T_stop_str      = '44640'  # [min] (minutes in January 2007)

    #   set_time_info() 
    #-------------------------------------------------------------------------------
    def set_basin_info(self):
#     def set_basin_info(self, site_prefix=None, case_prefix='Test1'):
#                        out_xres_sec=None, out_yres_sec=None,
#                        out_xres_m=None, out_yres_m=None,
#                        out_bounds=None, max_river_width=100.0,
#                        A_out_km2=None, Qbase_out=None,
#                        max_bankfull_depth=None, DRAINS_TO_OCEAN=False):
                
        """Set info for a specific basin of interest"""
        #### Could rename to "set_site_info()"  ???

        #----------------------------------------------------------------------- 
        # TopoFlow uses a "site_prefix" for all of the filenames in a
        # data set that pertain to the geographic location (the "site").
        # These files describe static properties of the location, such as
        # topography and soil.  The default site prefix in this module
        # is "Baro-Gam_60sec"".
        # 
        # Topoflow uses a "case_prefix" for all of the filenames
        # in a data set that describe a particular model scenario
        # (the "case" under consideration).  These files describe things
        # that can change from one model run to the next for the same site. 
        # The default case prefix in this module is "Test1"". 
        # Note: Component CFG filenames always start with the case prefix.
        # 
        # By simply changing the information in this function, this
        # notebook can generate TopoFlow input files for any river basin
        # in Ethiopia.  However, please heed the Important Warning below.
        # 
        # out_bounds = Geographic bounding box for the chosen river basin
        # out_xres_sec = the spatial grid cell xsize to use, in arcseconds
        # out_yres_sec = the spatial grid cell ysize to use, in arcseconds
        # 
        # Important Warning !!
        # Some of the input grid files generated by this notebook are
        # computed using empirical, power-law estimates, based on a grid
        # of total contributing area (TCA).
        # 
        # For these, it is necessary -- at a minimum  -- to know the values
        # of a few key variables at the basin outlet.  These must be
        # determined from the literature or other data sets. If these
        # are not set to reasonable values, the resulting predictions
        # will be meaningless.
        #-----------------------------------------------------------------------
        ## if (site_prefix is None):
        if (self.site_prefix is None):
            print('ERROR: site_prefix is required.')
            print()
            return
        # self.site_prefix = site_prefix    #######
        # self.case_prefix = case_prefix    #######

        #-----------------------        
        # site_prefix examples
        #-----------------------
        # Basins in Ethiopia
        # site_prefix = 'Baro-Gam_60sec'
        # site_prefix = 'Awash_120sec'
        # site_prefix = 'Omo_120sec'
        # site_prefix = 'Shebelle-Imi_60sec'
        # site_prefix = 'Ganale-Welmel_60sec'
        # site_prefix = 'Guder_30sec'
        # site_prefix = 'Muger_30sec'
        # Basins in Kenya
        # site_prefix = 'Tana_120sec'
        # site_prefix = 'Galana_120sec'
                
        if (self.site_prefix == 'Tana_120sec'):
            # Tana River, Kenya
            self.out_xres_sec = 120.0    # [arcseconds]
            self.out_yres_sec = 120.0    # [arcseconds]
            # Set the geographic bounding box and the grid cell size that
            # will be used for the TopoFlow model run, where
            #    Bounds = [ minlon, minlat, maxlon, maxlat ]
            # The bounding box MUST contain the entire watershed polygon.
            self.out_bounds = [ 36.4, -3.15, 40.75, 0.69]

            self.max_river_width = 120.0  # [meters]  from Google Maps or Google Earth; near Kipini
            self.A_out_km2 = 94284.4  # total contributing area at basin outlet [km2]
            self.Qbase_out = 40.0     # estimated baseflow discharge at basin outlet [m^3 / s]
            self.max_bankfull_depth = 8.0   # [meters]    # rough estimate
            self.DRAINS_TO_OCEAN = True
        #-----------------------------------------------------------------------
        if (self.site_prefix == 'Galana_120sec'):
            # Galana River, Kenya
            self.out_xres_sec = 120.0    # [arcseconds]
            self.out_yres_sec = 120.0    # [arcseconds]
            # Set the geographic bounding box and the grid cell size that
            # will be used for the TopoFlow model run, where
            #    Bounds = [ minlon, minlat, maxlon, maxlat ]
            # The bounding box MUST contain the entire watershed polygon.
            self.out_bounds = [ 36.3, -3.56, 40.23, -0.94]

            self.max_river_width = 140.0  # [meters]  from Google Maps or Google Earth
            self.A_out_km2 = 41822.0   # total contributing area at basin outlet [km2]
            self.Qbase_out = 40.0     # estimated baseflow discharge at basin outlet [m^3 / s]
            self.max_bankfull_depth = 8.0   # [meters]    # from literature or data
            self.DRAINS_TO_OCEAN = True
        #-----------------------------------------------------------------------
        if (self.site_prefix == 'Baro-Gam_60sec'):
            # Baro River at Gambella, Ethiopia
            self.out_xres_sec = 60.0    # [arcseconds]
            self.out_yres_sec = 60.0    # [arcseconds]
            # Set the geographic bounding box and the grid cell size that
            # will be used for the TopoFlow model run, where
            #    Bounds = [ minlon, minlat, maxlon, maxlat ]
            # The bounding box MUST contain the entire watershed polygon.
            # Next one gives ncols = 133, nrows = 128.
            ## out_bounds = [ 34.22125, 7.3704166666, 36.43791666666, 9.50375]
            # Next one gives ncols = 134, nrows = 129
            self.out_bounds = [ 34.221249999999, 7.362083333332, 36.450416666666, 9.503749999999]

            self.max_river_width = 140.0  # [meters]  from Google Maps or Google Earth
            self.A_out_km2 = 23567.7  # total contributing area at basin outlet [km2]
            self.Qbase_out = 40.0     # estimated baseflow discharge at basin outlet [m^3 / s]
            self.max_bankfull_depth = 8.0   # [meters]    # from literature or data
        #-----------------------------------------------------------------------
        if (self.site_prefix == 'Awash_120sec'):
            #-----------------------------------------------------------------
            # Note: The Awash river basin terminates in Lake Abbe (or Abhe)
            #       at an elevation near 240 m.  The "fill_pits" routine
            #       can handle this situation if we set the value of a grid
            #       cell in Lake Abhe to the "closed-basin code": -32000.
            #       Otherwise, the river is routed incorrectly to ocean.
            #       See the "Fill Pits" section for more info.
            #       According to Google Earth measuring tool, the Awash
            #       River reaches a max width of about 35 meters and then
            #       fans out in a distributary system of smaller channels
            #       that flow into Lake Abbe.
            #-----------------------------------------------------------------
            # Awash River, Ethiopia
            self.out_xres_sec = 120.0    # [arcseconds]
            self.out_yres_sec = 120.0    # [arcseconds]
            # Bounds = [ minlon, minlat, maxlon, maxlat ]
            # The bounding box MUST contain the entire watershed polygon.
            self.out_bounds = [37.92, 7.75, 43.36, 12.44]
            self.max_river_width = 35.0  # [meters]  from Google Maps or Google Earth
            self.A_out_km2 = 92638.1   # total contributing area at basin outlet [km2]
            self.Qbase_out = 5.0      # estimated baseflow discharge at basin outlet [m^3 / s]
            self.max_bankfull_depth = 8.0   # [meters]    # rough estimate; no data
            self.DRAINS_TO_OCEAN = False   # (Drains to Lake Abhe or Abbe; but routed wrong in DEM)
        #-----------------------------------------------------------------------
        if (self.site_prefix == 'Omo_120sec'):
            #---------------------------------------------------------------
            # Note: The Omo River terminates in Lake Turkana.
            # Lake Chew Bahir (Stephanie Wildlife Sanctuary), lies at the
            # bottom of a closed basin, with a lowest elevation around
            # 497.181 (col, row = 49, 140).  The lake level fluctuates.
            #---------------------------------------------------------------
            # Max river width is around 150+ meters before the first
            # distributary occurs just north of Lake Turkana.
            # Floodplain width may be greater than 5 near this point.
            #---------------------------------------------------------------
            # Omo-Gabe River, Ethiopia    
            self.out_xres_sec = 120.0    # [arcseconds]
            self.out_yres_sec = 120.0    # [arcseconds]
            # Bounds = [ minlon, minlat, maxlon, maxlat ]
            # The bounding box MUST contain the entire watershed polygon.
            self.out_bounds = [35.19, 4.00, 38.48, 9.41]
            self.max_river_width = 150.0  # [meters]  from Google Maps or Google Earth
            self.A_out_km2 = 72736.0   # total contributing area at basin outlet [km2]
            self.Qbase_out = 5.0      # estimated baseflow discharge at basin outlet [m^3 / s]
            self.max_bankfull_depth = 9.0   # [meters]    # rough estimate; no data
            self.DRAINS_TO_OCEAN = False   # (Drains to Lake Abhe or Abbe; but routed wrong in DEM)
        #-----------------------------------------------------------------------
        if (self.site_prefix == 'Shebelle-Imi_60sec'):
            # Shebelle River at Imi, Ethiopia
            self.out_xres_sec = 60.0    # [arcseconds]
            self.out_yres_sec = 60.0    # [arcseconds]
            # Bounds = [ minlon, minlat, maxlon, maxlat ]
            # The bounding box MUST contain the entire watershed polygon.
            # out_bounds = [ 38.159583333333, 6.324583333333, 43.559583333333, 9.899583333333]  # (for 30 arcseconds)
            self.out_bounds = [38.159583333333, 6.319583333333, 43.559583333333, 9.899583333333] 
            self.max_river_width = 130.0  # [meters]  from Google Maps or Google Earth
            self.A_out_km2 = 90662.1   # total contributing area at basin outlet [km2]
            self.Qbase_out = 50.0      # estimated baseflow discharge at basin outlet [m^3 / s]
            self.max_bankfull_depth = 7.0   # [meters]    # from literature or data
        #-----------------------------------------------------------------------
        # Actually, this tributary of the Ganale River is called the Welmel Shet River.
        # The outlet is at the border between Oromia and Somali Regions.
        if (self.site_prefix == 'Ganale-Welmel_60sec'):
            # Ganale-Welmel River at border, Ethiopia
            self.out_xres_sec = 60.0    # [arcseconds]
            self.out_yres_sec = 60.0    # [arcseconds]
            # Bounds = [ minlon, minlat, maxlon, maxlat ]
            # The bounding box MUST contain the entire watershed polygon.
            self.out_bounds = [39.174583333333, 5.527916666666, 41.124583333333, 7.098749999999]  # (3 arcseconds) 
            self.max_river_width = 40.0  # [meters]  from Google Maps or Google Earth
            self.A_out_km2 = 15241.7   # total contributing area at basin outlet [km2]
            self.Qbase_out = 3.0      # estimated baseflow discharge at basin outlet [m^3 / s]
            self.max_bankfull_depth = 2.0   # [meters]    # from literature or data    
         #-----------------------------------------------------------------------   
        if (self.site_prefix == 'Guder_30sec'):
            # Guder River at Blue Nile confluence, Ethiopia
            self.out_xres_sec = 30.0    # [arcseconds]
            self.out_yres_sec = 30.0    # [arcseconds]
            # Bounds = [ minlon, minlat, maxlon, maxlat ]
            # The bounding box MUST contain the entire watershed polygon.
            self.out_bounds = [37.149583333333, 8.596250000000, 38.266250000000, 9.904583333333]
            self.max_river_width = 20.0  # [meters]  from Google Maps or Google Earth
            self.A_out_km2 = 6487.8   # total contributing area at basin outlet [km2]
            self.Qbase_out = 2.0      # estimated baseflow discharge at basin outlet [m^3 / s]
            self.max_bankfull_depth = 2.0   # [meters]    # from literature or data
         #----------------------------------------------------------------------- 
        if (self.site_prefix == 'Muger_30sec'):
            # Muger River at Blue Nile confluence, Ethiopia
            self.out_xres_sec = 30.0    # [arcseconds]
            self.out_yres_sec = 30.0    # [arcseconds]
            # Bounds = [ minlon, minlat, maxlon, maxlat ]
            # The bounding box MUST contain the entire watershed polygon.
            self.out_bounds = [37.807916666667, 8.929583333333, 39.032916666667, 10.112916666667]
            self.max_river_width = 45.0  # [meters]  from Google Maps or Google Earth
            self.A_out_km2 = 6924.12    # total contributing area at basin outlet [km2]
            self.Qbase_out = 3.0      # estimated baseflow discharge at basin outlet [m^3 / s]
            self.max_bankfull_depth = 2.0   # [meters]    # from literature or data

    #   set_basin_info()
    #-------------------------------------------------------------------------------
    def set_met_info(self):
                ## CHIRPS=True, CHIRPS2=False, 
                ## GLDAS=False, GPM=False):

        """Set rainfall info to use for rainfall grid stack"""
        rts_extension = '_' + self.date_range_str + '.rts'
        sp = self.site_prefix

        if (self.NGEN_CSV):
            print('Will use NextGen CSV rainfall.')
            self.src_rain_dir  = None
            self.rain_rts_file = None
            ## self.rain_csv_file = self.site_prefix + '.csv'  ## in met_dir
            self.rain_csv_file = self.site_prefix + '_' + self.date_range_str + '.csv'
            self.met_dt_str    = '3600.0'  # [secs]  60 * 60 (secs in 1 hr)
            self.rain_type_str = 'NGEN_CSV'
  
        if (self.TEST_RAIN):
            print('Will use test rainfall.')
            self.src_rain_dir  = None
            self.rain_rts_file = None
            self.rain_csv_file = self.case_prefix + '_rain_rates.txt'  ## in met_dir
            ## self.rain_ts_file  = self.case_prefix + '_rain_rates.txt'  ## in met_dir
            ################################################################
            self.met_dt_str    = '3600.0'  # [secs]  60 * 60 (secs in 1 hr)
            self.rain_type_str = 'TEST_RAIN'
                                           
        if (self.CHIRPS):
            # src_rain_dir must contain only "rfe" files for specified time range
            print('Will use CHIRPS rainfall product.')
            self.src_rain_dir  = self.src_chirps_dir
            self.rain_rts_file = 'CHIRPS_Rain_' + sp + rts_extension
            self.met_dt_str    = '21600.0'  # [secs]  60 * 60 * 6 (secs in 6 hrs)
            self.rain_type_str = 'CHIRPS'

        if (self.CHIRPS2):
            # Use this for a 10-year baseline for stats
            print('Will use CHIRPS rainfall product.')           
            self.src_rain_dir  = self.src_chirps2_dir
            self.rain_rts_file = 'CHIRPS_Rain_' + sp + rts_extension
            self.met_dt_str    = '21600.0'  # [secs]  60 * 60 * 6 (secs in 6 hrs)
            self.rain_type_str = 'CHIRPS'
              
        if (self.GLDAS):
            # src_rain_dir must contain only "nc4" files for specified time range
            print('Will use GLDAS rainfall product.')
            self.src_rain_dir  = self.src_gldas_dir
            self.rain_rts_file = 'GLDAS_Rain_' + sp + rts_extension
            self.met_dt_str    = '10800.0'  # [secs]  60 * 60 * 3 (secs in 3 hrs)
            self.rain_type_str = 'GLDAS'
    
        if (self.GPM):
            # src_rain_dir must contain only "nc4" files for specified time range
            print('Will use GPM rainfall product.')
            self.src_rain_dir  = self.src_gpm_dir
            self.rain_rts_file = 'GPM_Rain_' + site_prefix + rts_extension
            self.met_dt_str    = '1800.0'  # [secs]  60 * 60 * 0.5 (secs in 1/2 hr)
            self.rain_type_str = 'GPM'

          # Note:  This is now set in: get_rainfall_grid_stack()
          # self.rain_rts_path = self.met_dir + self.rain_rts_file

    #   set_met_info()
    #-------------------------------------------------------------------------------
    def copy_cfg_file_templates(self):
    
        """Copy a set of default CFG files into new basin directory"""
        #---------------------------------------------------- 
        # Copy a set of default "CFG template" files from
        # "cfg_template_dir" to "cfg_dir".  Settings that
        # can currently be changed in the next block of this
        # notebook are wrapped in "double curly brackets",
        # e.g. "{{site_prefix}}".
        #----------------------------------------------------
        ctd = self.cfg_template_dir
        cfg_template_list = sorted(glob.glob( ctd + '*.cfg' ))
        
        for cfg_template_file in cfg_template_list:
            parts = cfg_template_file.split( os.sep )
            new_cfg_file = self.cfg_dir + parts[-1]
            shutil.copyfile( cfg_template_file, new_cfg_file)

        # Copy the default "provider_file"
        provider_file = self.case_prefix + '_providers.txt'
        shutil.copyfile( ctd + provider_file, self.cfg_dir + provider_file )
        
    #   copy_cfg_file_templates()
    #-------------------------------------------------------------------------------
    def get_cfg_files(self, REPORT=True):
    
        """Create CFG files from CFG templates"""
        in_dir   = self.basin_dir
        out_dir  = self.output_dir
        cp       = self.case_prefix
        prefix   = self.cfg_template_dir + cp    ###### RENAME #######
        cfg_list = []
        
#         print('### In get_cfg_files():')
#         print('### case_prefix =', cp)
#         print('### cfg_dir =', self.cfg_dir)

        #--------------------------------
        # Create new path_info CFG file
        #--------------------------------
        cfg_template = prefix + '_path_info.cfg'
        cfg_file     = self.cfg_dir + cp + '_path_info.cfg'
        vars = ['in_directory', 'out_directory', 'site_prefix', 'case_prefix' ]
        vals = [ in_dir, out_dir, self.site_prefix, self.case_prefix ]
        ct.make_new_cfg_file( vars, vals, cfg_template, cfg_file )
        cfg_list.append( cfg_file )

        #--------------------------------
        # Create new time_info CFG file
        #--------------------------------
        cfg_template = prefix + '_time_info.cfg'
        cfg_file     = self.cfg_dir + cp + '_time_info.cfg'
        vars = ['start_date', 'start_time', 'end_date', 'end_time' ]
        vals = [ self.start_date_str, '  00:00:00', self.end_date_str, '  00:00:00']
        ct.make_new_cfg_file( vars, vals, cfg_template, cfg_file )
        cfg_list.append( cfg_file )
        
        #--------------------------------------
        # Create new TopoFlow driver CFG file
        #--------------------------------------
        compute_stat_str       = '0'
        create_indicator_str   = '0'
        create_media_files_str = '1'
        if (self.CHIRPS2):
            compute_stat_str = '1'  #######################
        cfg_template = prefix + '_topoflow.cfg'
        cfg_file     = self.cfg_dir + cp + '_topoflow.cfg'
        vars = ['dt', 'T_stop_model', 'COMPUTE_STAT_GRIDS', 'CREATE_INDICATORS',
                'CREATE_MEDIA_FILES' ]
        vals = [self.chan_dt_str, self.T_stop_str, compute_stat_str,
                create_indicator_str, create_media_files_str ]
        ct.make_new_cfg_file( vars, vals, cfg_template, cfg_file )
        cfg_list.append( cfg_file )
        
        #---------------------------------------
        # Create new meteorology CFG file
        # GMT Offset is 3 for Ethiopia & Kenya
        #---------------------------------------
        cfg_template = prefix + '_meteorology.cfg'
        cfg_file     = self.cfg_dir + cp + '_meteorology.cfg'
        vars = ['dt', 'P_type', 'P', 'GMT_offset', 'start_month', 'start_day',
                'start_hour', 'save_grid_dt', 'save_pixels_dt', 'NGEN_CSV' ]
        if (self.NGEN_CSV):
            vals = [self.met_dt_str, 'Time_Series', self.rain_csv_file,  ######
                    self.GMT_offset_str, self.start_month_str,
                    self.start_day_str, '0.0', 
                    self.met_dt_str, self.met_dt_str, 'Yes']
        elif (self.TEST_RAIN):
            vals = [self.met_dt_str, 'Time_Series', self.rain_csv_file,  ######
                    self.GMT_offset_str, self.start_month_str,
                    self.start_day_str, '0.0', 
                    self.met_dt_str, self.met_dt_str, 'No']        
        else:        
            vals = [self.met_dt_str, 'Grid_Sequence', self.rain_rts_file,
                    self.GMT_offset_str, self.start_month_str,
                    self.start_day_str, '0.0', 
                    self.met_dt_str, self.met_dt_str, 'No']
        ct.make_new_cfg_file( vars, vals, cfg_template, cfg_file )
        cfg_list.append( cfg_file )
        
        #----------------------------------------------
        # Create new channels_kinematic_wave CFG file
        #----------------------------------------------
        cfg_template = prefix + '_channels_kinematic_wave.cfg'
        cfg_file     = self.cfg_dir + cp + '_channels_kinematic_wave.cfg'
        vars = ['dt', 'FLOOD_OPTION', 'ATTENUATE', 'save_grid_dt', 'save_pixels_dt' ]
        vals = [self.chan_dt_str, '0', '0', self.met_dt_str, self.met_dt_str] 
        ct.make_new_cfg_file( vars, vals, cfg_template, cfg_file )
        cfg_list.append( cfg_file )
        
        #----------------------------------------------
        # Create new channels_diffusive_wave CFG file
        #----------------------------------------------
        cfg_template = prefix + '_channels_diffusive_wave.cfg'
        cfg_file     = self.cfg_dir + cp + '_channels_diffusive_wave.cfg'
        vars = ['dt', 'FLOOD_OPTION', 'ATTENUATE', 'save_grid_dt', 'save_pixels_dt' ]
        vals = [self.chan_dt_str, '0', '0', self.met_dt_str, self.met_dt_str]  
        ct.make_new_cfg_file( vars, vals, cfg_template, cfg_file )
        cfg_list.append( cfg_file )
        
        if (REPORT):
            print('Created CFG files:')
            for name in cfg_list:
                print('  ' + name) 
            print()
               
    #   get_cfg_files()
    #-------------------------------------------------------------------------------
    def copy_outlets_file(self):
    
        """Copy an "outlets file" of grid cells to monitor"""
        #--------------------------------------------------------------
        # In TopoFlow component CFG files, flags can be set to tell
        # TopoFlow to write values of chosen gridded variables to a
        # file, to create a "grid stack"", indexed by time.
        # Grids are saved at a time interval set by "save_grid_dt"".
        # 
        # Other flags in a CFG file can be set to tell TopoFlow to
        # write values of chosen gridded variables to a file, but
        # only at a specified set of grid cells.  These "monitored
        # grid cells" or "virtual gauges" are set in an "outlets
        # file" named ""[case_prefix]_outlets.txt"".  It is a simple,
        # multi-column text file, like this:
        #
        # ------------------------------------------------------------
        #  Monitored Grid Cell (Outlet) Information
        # -------------------------------------------------------------------------------------
        #     Column       Row     Area [km^2]      Relief [m]      Lon [deg]     Lat [deg]
        # -------------------------------------------------------------------------------------
        #         4         79        24647.3         2433.62       34.296250     8.1787500
        #        14         76        24011.0         2425.90       34.462917     8.2287500
        #        29         79        22972.0         2393.16       34.712917     8.1787500
        #        48         59          507.3         1434.04       35.029583     8.5120833
        #
        # It is not necessary for the Area and Relief columns to
        # contain valid values (they are for reference, but are unused).
        # However, the column, row, longitude and latitude of each grid
        # cell to be monitored must be specified.  They must match
        # columns and rows in the DEM that is being used for the model
        # run.  Creating an outlets file therefore requires a "human in
        # the loop"", using interactive GIS (Geographic Information
        # System) software and making intelligent choices.  It cannot
        # be automated.  Here we just copy an existing outlets file for
        # one the selected basin.
        #-----------------------------------------------------------------
        outlet_file = self.case_prefix + '_outlets.txt'
        src_dir     = self.src_outlet_dir + self.site_prefix + os.sep
        src_file    = src_dir + outlet_file
        new_file    = self.cfg_dir + outlet_file

        if (os.path.exists( src_file )):
            shutil.copyfile( src_file, new_file )
        else:
            print('WARNING: No outlets file found for: ' + self.site_prefix)
            print('         You will need to create one manually.')
            print()
            
    #   copy_outlets_file()
    #-------------------------------------------------------------------------------
    def get_outlets_file(self, STEP_BACKS=0):
    
        """Create a default "outlets file" of grid cells to monitor"""
        #--------------------------------------------------------------
        # In TopoFlow component CFG files, flags can be set to tell
        # TopoFlow to write values of chosen gridded variables to a
        # file, to create a "grid stack"", indexed by time.
        # Grids are saved at a time interval set by "save_grid_dt"".
        # 
        # Other flags in a CFG file can be set to tell TopoFlow to
        # write values of chosen gridded variables to a file, but
        # only at a specified set of grid cells.  These "monitored
        # grid cells" or "virtual gauges" are set in an "outlets
        # file" named ""[case_prefix]_outlets.txt"".  It is a simple,
        # multi-column text file, like this:
        #
        # ------------------------------------------------------------
        #  Monitored Grid Cell (Outlet) Information
        # -------------------------------------------------------------------------------------
        #     Column       Row     Area [km^2]      Relief [m]      Lon [deg]     Lat [deg]
        # -------------------------------------------------------------------------------------
        #         4         79        24647.3         2433.62       34.296250     8.1787500
        #        14         76        24011.0         2425.90       34.462917     8.2287500
        #        29         79        22972.0         2393.16       34.712917     8.1787500
        #        48         59          507.3         1434.04       35.029583     8.5120833
        #
        # It is not necessary for the Area and Relief columns to
        # contain valid values (they are for reference, but are unused).
        # However, the column, row, longitude and latitude of each grid
        # cell to be monitored must be specified.  They must match
        # columns and rows in the DEM that is being used for the model
        # run.  Here, we attempt to create an outlets file automatically
        # by using the grid cell with the largest total contributing area
        # or a grid cell a bit upstream from that one.  While it is
        # better to have a "human in the loop", this should work in most
        # cases.
        #-----------------------------------------------------------------
        outlet_file = self.case_prefix + '_outlets.txt'
        new_file    = self.cfg_dir + outlet_file
        
        if not(hasattr(self, 'd8')):
            self.get_d8_flow_grid()
        A_grid  = self.d8.A
        d8_grid = self.d8.d8_grid
            
        #---------------------------------------------------        
        # Or open and read the D8 area and flow code files
        #---------------------------------------------------
        # A_grid  = rtg_files.read_grid( self.d8_area_file, self.grid_info)
        # d8_grid = rtg_files.read_grid( self.d8_code_file, self.grid_info) 

        #--------------------------------------------------
        # Find the grid cell that has the largest D8 area
        #--------------------------------------------------
        Amax = A_grid.max()
        w = np.where( A_grid == Amax )
        row = w[0][0]
        col = w[1][0]
        old_row_list = list()
        old_col_list = list()
        old_row_list.append(row)
        old_col_list.append(col)
        good_indices = np.ones(8, dtype='bool')
        
        #---------------------------------------------------------        
        # Find the D8 child (or just D8 neighbor) of this cell
        # that has the largest TCA  (OPTIONAL, but maybe better)
        #---------------------------------------------------------
        # Successful tests on Tana River, for SB in {0,1,2}.
        #--------------------------------------------------------
        ## STEP_BACKS = 2
        if (STEP_BACKS > 0):
            for k in range(STEP_BACKS):
                # Order:  NW, N, NE, W, E, SW, S, SE
                rows = [row-1, row-1, row-1, row, row, row+1, row+1, row+1]
                cols = [col-1, col, col+1, col-1, col+1, col-1, col, col+1]
                rows = np.array( rows )
                cols = np.array( cols )
                for j in range(8):
                    if (rows[j] in old_row_list) and (cols[j] in old_col_list):
                        good_indices[j] = False
                rows = rows[ good_indices ]
                cols = cols[ good_indices ]
                #-------------------------------
                areas = A_grid[ [rows, cols] ]
                w     = np.argmax( areas )
                Amax  = areas[ w[0] ]
                col   = cols[ w[0] ]
                row   = rows[ w[0] ]
                old_row_list.append(row)
                old_col_list.append(col)
    
        #-----------------------------
        # Get the outlet information
        #-----------------------------
        area   = Amax
        relief = -9999

        #---------------------------------------        
        # Compute lon and lat from col and row
        #--------------------------------------------
        # Assume that pixel geometry is fixed-angle
        # Lat and lon is for center of grid cell
        #--------------------------------------------
        if (self.grid_info.pixel_geom == 0):
            minlon   = self.grid_info.x_west_edge
            maxlat   = self.grid_info.y_north_edge
            ## minlat   = self.grid_info.y_south_edge
            xres_deg = self.grid_info.xres / 3600.0
            yres_deg = self.grid_info.yres / 3600.0
            lon = minlon + (xres_deg * col) + (xres_deg/2)
            lat = maxlat - (yres_deg * row) - (yres_deg/2)
            ## lat = minlat + (yres_deg * row) + (yres_deg/2)
        else:
            print('ERROR: Cannot yet compute outlet lat and lon')
            print('       for DEMs with fixed-length grid cells.')
            print()
            lon = -999.0
            lat = -999.0
        #------------------------------------
        # lon = 35.029583  # (for testing)
        # lat = 8.5120833 
        
        outlets.write_outlet_file(new_file, col, row, area, relief, lon, lat)
          
    #   get_outlets_file()
    #-------------------------------------------------------------------
    def get_ngen_dem( self, DEM_out_file=None,
                     cat_id_str='cat-67',
                     GEO=True, ALBERS=False ):

        if (DEM_out_file is None):
            DEM_out_file = self.topo_dir + cat_id_str + '_rawDEM.tif'

        c = hfu.catchment()
        # c.print_info( cat_id_str )
        dem = c.get_dem(cat_id_str=cat_id_str,   #########
                RESAMPLE_ALGO='bilinear',
                CLIP_GEO=GEO, CLIP_ALBERS=ALBERS,
                zfactor=100, REPORT=True, SAVE_DEM=True,
                DEM_out_file=DEM_out_file)
                
        self.out_bounds   = c.DEM_bounds
        self.out_xres_sec = c.DEM_xres  ####################
        self.out_yres_sec = c.DEM_yres  ####################     
        # self.out_bounds = c.get_bounding_box( cat_id_str )

    #   get_ngen_dem()      
    #-------------------------------------------------------------------------------
    def clip_dem_from_source(self):
    
        """Clip a source DEM to a bounding box and resample"""
        #------------------------------------------------------------------
        # Here, we use the TopoFlow regrid utility to clip a DEM for
        # a larger region to the geographic bounding box of a particular
        # river basin in that region.
        #
        # This utility uses the gdal.warp() function in the GDAL Python
        # package.  At the same time, we resample (via spatial bilinear
        # interpolation) the resulting DEM to a different (e.g. coarser)
        # spatial resolution.  The source DEM may have any grid cell
        # size (e.g. 3 arcseconds, roughly 93 meters), and the new DEM
        # can have a different grid cell size (e.g. 60 arcseconds,
        # roughly 1800 meters).  Both the source DEM and new DEM are
        # stored in GeoTIFF format.  Resampling typically causes the
        # bounding box to change slightly.
        #------------------------------------------------------------------
        ### in_nodata = ????
        ### out_nodata = -9999.0

        in_file  = self.src_dem_dir + self.src_dem_file
        out_file = self.DEM_tif_file

        regrid.regrid_geotiff(in_file=in_file, out_file=out_file, 
                           out_bounds=self.out_bounds,
                           out_xres_sec=self.out_xres_sec,
                           out_yres_sec=self.out_yres_sec,
                           ### in_nodata=None, out_nodata=None, 
                           RESAMPLE_ALGO='bilinear', REPORT=True)
                           
    #   clip_dem_from_source()
    #-------------------------------------------------------------------------------
    def read_dem_as_geotiff(self, tif_file=None):

        """Read new, clipped DEM from a GeoTIFF file"""
        #---------------------------------------------------------------
        # Here we import a DEM in GeoTIFF format.
        # See: read_dem_as_rtg() for RiverTools Grid (RTG) format.
        # See: read_dem_as_netcdf() for NetCDF (.nc) format.
        # 
        # Most of the TopoFlow utilities use grids saved in the
        # RiverTools Grid (RTG) format, which is a generic, binary,
        # row-major format.  Georeferencing information for the grid
        # is stored in a small, separate text file in RiverTools Info
        # (RTI) format.  When the rti_file argument is specified,
        # georeferencing information is also saved in the RTI file
        # format for later use.
        #---------------------------------------------------------------       
        if (tif_file is None):
            tif_file = self.topo_dir + self.site_prefix + '_rawDEM.tif'

#         print('### In read_dem_as_geotiff:')
#         print('### rti_file =', self.rti_file )

        self.DEM = import_grid.read_from_geotiff( tif_file, REPORT=True,
                                                  rti_file=self.rti_file)
        self.grid_info = rti_files.read_info( self.rti_file )
        
    #   read_dem_as_geotiff()
    #-------------------------------------------------------------------------------
    def read_dem_as_rtg(self, rtg_file=None):

        """Read DEM from an RTG format file (with RTI file)"""
        if (rtg_file is None):
            rtg_file = self.topo_dir + self.site_prefix + '_rawDEM.rtg'
        self.DEM = import_grid.read_from_rtg( rtg_file, REPORT=True)
        self.grid_info = rti_files.read_info( rtg_file )

    #   read_dem_as_rtg()
    #-------------------------------------------------------------------------------
    def read_dem_as_netcdf(self, nc_file=None):

        """Read DEM from a netCDF format file"""
        if (nc_file is None):
            nc_file = self.topo_dir + self.site_prefix + '_rawDEM.nc'
        self.DEM = import_grid.read_from_netcdf( nc_file, REPORT=True)
        ## rti_file = self.site_prefix + '.rti'
        ## self.grid_info = rti_files.read_info( rtg_file )
        
    #   read_dem_as_netcdf()
    #-------------------------------------------------------------------------------
    def write_dem_as_rtg(self, rtg_file=None):
    
        """Write DEM to a RiverTools Grid (RTG) format file"""
        if (rtg_file is None):
            rtg_file = self.topo_dir + self.site_prefix + '_rawDEM.rtg' 
        out_type  = 'FLOAT'   # (regardless of grid_info.data_type)
        ## out_type = grid_info.data_type
        rtg_files.write_grid( self.DEM, rtg_file, self.grid_info,
                              RTG_type=out_type, SILENT=True)

        if not(self.SILENT):
            print('Finished saving GeoTIFF file as RTG to:')
            print('  ' + rtg_file )
            print()
   
    #   write_dem_as_rtg()
    #-------------------------------------------------------------------------------
#     def get_profile_smoothed_dem(self):
#     
#         """Create a DEM with smoother slopes"""
#         #------------------------------------------------------------
#         # The new_slopes.py utility has a different method,
#         # and is currently used instead.
#         #-------------------------------------------------------------
#         # For so-called "mature" landscapes this "profile smoothing"
#         # algorithm works well, and results in a DEM with smoothly
#         # decreasing, nonzero channel slopes everywhere.  However,
#         # the landscape of the Baro River basin is not a good
#         # candidate because it is not mature.
#         #-------------------------------------------------------------
#         os.chdir( self.topo_dir )
#         c = smooth_DEM.DEM_smoother()
#         c.DEBUG = True
#         cfg_file = self.cfg_dir + self.case_prefix + '_dem_smoother.cfg'
#         c.initialize( cfg_file=cfg_file, mode='driver')
#         c.update()
# 
#     #   get_profile_smoothed_dem()
    #-------------------------------------------------------------------------------
    def replace_negative_vals_in_dem(self):

        """Replace negative values in DEM with zeros"""
        #--------------------------------------------------------
        # Use this to prevent problems that can occur for DEMs
        # where the main river basin drains to the ocean.  But
        # some DEMs have valid negative elevations.
        #-------------------------------------------------------- 
        if (self.DRAINS_TO_OCEAN):
            self.DEM[ self.DEM < 0 ] = 0.0

    #   replace_negative_values_in_dem()
    #-------------------------------------------------------------------------------
    def fill_pits_in_dem(self):

        """Fill depressions in the DEM and save""" 
        # This may be necessary to make the DEM hydrologically sound.

        # shp   = self.DEM.shape
        # nrows = shp[0]
        # ncols = shp[1]

        data_type = self.grid_info.data_type   # e.g. "INTEGER" or "FLOAT"
        ncols = self.grid_info.ncols
        nrows = self.grid_info.nrows

        if (self.site_prefix == 'Awash_120sec'):
            #-----------------------------------------------------------------
            # Note: The Awash River basin terminates in Lake Abhe (or Abbe)
            #       at an elevation near 240 m.  The "fill_pits" routine
            #       can handle this situation if we set the value of a grid
            #       cell in Lake Abhe to the "closed-basin code": -32000.
            #       Otherwise, the river is routed incorrectly to ocean.
            #-----------------------------------------------------------------
            self.DEM[36, 115] = -32000
            ## self.DEM[37, 115] = -32000
            ## self.DEM[37, 114] = -32000   # elev: 239.942

        if (self.site_prefix == 'Omo_120sec'):
            #-----------------------------------------------------------------
            # Note: The Omo River basin terminates in Lake Turkana
            #       at an elevation near 360.4 m.  Lake Chew Bahir lies to
            #       the northeast and is the bottom of a closed basin.
            #       The "fill_pits" routine can handle this situation if we
            #       set the value of a grid cell in Lake Chew Bahir to the
            #      "closed-basin code": -32000.
            #       Otherwise, the river is routed incorrectly over a ridge.
            #-----------------------------------------------------------------
            self.DEM[140, 49] = -32000
    
        fill_pits.fill_pits( self.DEM, data_type, ncols, nrows, 
                             SILENT=False)
    
        rtg_files.write_grid( self.DEM, self.new_DEM_file,
                              self.grid_info, SILENT=False)

    #   fill_pits_in_dem()
    #-------------------------------------------------------------------------------
    def get_d8_flow_grid(self):
    
        """Compute & save the D8 flow direction grid"""
        # TopoFlow includes a component called d8_global that can
        # compute a grid of D8 flow direction codes (Jenson 1984
        # convention), as well as several additional, related grids
        # such as a grid of total contributing area (TCA).  TopoFlow
        # components are configured through the use of configuration
        # files, which are text files with the extension ".cfg".
        # Therefore, we now need to make use of the CFG file for the 
        # D8-Global component, Test1_d8_global.cfg.

        d8 = d8_global.d8_component()
        d8.DEBUG = False
        
        ######################################
        # CHECK IF fill_pits IS CALLED AGAIN
        ### d8.FILL_PITS_IN_Z0 = 0
        ######################################
        cfg_file = self.cfg_dir + self.case_prefix + '_d8_global.cfg' # (need full path)
        time = 0.0
        d8.initialize( cfg_file=cfg_file, SILENT=False, REPORT=True )

        # If elevation <= nodata, set D8 flow code to 0 (undefined)
        if (self.DRAINS_TO_OCEAN):
            d8.DEM_nodata = 0.0
        else:
            d8.DEM_nodata = -9999.0

        d8.update( time, REPORT=True )
        self.d8 = d8   ##################
        
        grid_info = rti_files.read_info( self.rti_file, REPORT=True )  ##########

        rtg_files.write_grid(d8.d8_grid, self.d8_code_file,
                             grid_info, RTG_type='BYTE')

    #   get_d8_flow_grid()
    #-------------------------------------------------------------------------------
    def get_d8_area_grid(self):

        """Compute and save the D8 total contributing area (TCA) grid"""

        rtg_files.write_grid( self.d8.A, self.d8_area_file, self.grid_info,
                              RTG_type='FLOAT', SILENT=False)

    #   get_d8_area_grid()
    #-------------------------------------------------------------------------------
    def get_d8_slope_grid(self):
        """Compute and save the D8 slope grid"""

        #---------------------------------------------
        # Method 1:  Standard D8, cell-to-cell slope
        #---------------------------------------------
        # self.d8.update_slope_grid()
        # rtg_files.write_grid( self.d8.S, d8_slope_file, self.grid_info,
        #                        RTG_type='FLOAT', SILENT=False)

        #-----------------------------------------------------------------
        # Method 2:  Better method for handling cells with slope of zero
        #-----------------------------------------------------------------
        # (Recommended)
        new_slopes.get_new_slope_grid(site_prefix=self.site_prefix,
                                      case_prefix=self.case_prefix,
                                      cfg_dir=self.cfg_dir,
                                      slope_file=self.d8_slope_file)

    #   get_d8_slope_grid()
    #-------------------------------------------------------------------------------
    def get_d8_aspect_grid(self):

        """Compute and save the D8 aspect grid"""        
        self.d8.update_aspect_grid()
        rtg_files.write_grid(self.d8.aspect, self.d8_aspect_file, self.grid_info,
                             RTG_type='FLOAT', SILENT=False)

    #   get_d8_aspect_grid()
    #-------------------------------------------------------------------------------
    def get_channel_width_grid(self):

        """Compute the estimated channel width grid"""
        #-----------------------------------------------------------------
        # First, use Google Maps or Google Earth to estimate the
        # width of the river at the outlet to your river basin,
        # in meters.  Here, we'll assume that width equals 140 meters.
        # 
        # The idea is to estimate the channel widths throughout the
        # basin (as a grid with the same dimensions as the DEM),
        # using an empirical power law of the form:  $w = c \, A^p$
        # where A is the total contributing area (TCA) that we
        # computed as a grid above and saved into "d8_area_file".
        # A typical value of p is 0.5.  The value that w should have
        # where A is maximum (e.g. the river outlet) is specified as g1. 
        #-----------------------------------------------------------------
        # Note:  max_river_width and channel_width_power should be
        #        specified in the set_basin_info() method.
        #-----------------------------------------------------------------       
        parameterize.get_grid_from_TCA(site_prefix=self.site_prefix,
                     topo_dir=self.topo_dir,
                     area_file=self.d8_area_file, out_file=self.width_file,
                     g1=self.max_river_width, p=self.channel_width_power)


    #   get_channel_width_grid()
    #-------------------------------------------------------------------------------
    def get_manning_n_grid(self):
    
        """Compute the estimated "Manning's n" grid"""
        #-----------------------------------------------------------------
        # In order to compute grids of river flow velocity and
        # discharge (volume flow rate), a very well-known, empirical
        # formula known as Manning's formula (see Wikipedia) is the
        # method used by default within TopoFlow.  This formula
        # includes a parameter called "Manning's n"", that
        # characterizes the roughness of the channel bed and
        # resulting frictional loss of momentum.  Typical values
        # in larger river channels range between 0.03 and 0.05.
        # Manning's formula can also be used for non-channelized,
        # overland flow, but then a much larger value of 0.2 to 0.3
        # should be used.
        # 
        # The following code uses a power-law estimate of the form:
        # $n = c \, A^p$, where A is the total contributing area
        # (TCA) grid, to create a grid of Manning's n values.  The
        # value that n should have where A is maximum (e.g. the
        # river outlet) is set as <b>g1</b>.  Similarly, the value
        # that n should have where A is minimum (e.g. on a ridge)
        # is set as <b>g2</b>.  The coefficient, c, and power, p,
        # are then set to match these constraints.
        #-----------------------------------------------------------------
        # Note:  min_manning_n and max_manning_n should be
        #        specified in the set_basin_info() method.
        #-----------------------------------------------------------------        
        parameterize.get_grid_from_TCA(site_prefix=self.site_prefix,
                     topo_dir=self.topo_dir,
                     area_file=self.d8_area_file,
                     out_file=self.manning_file,
                     g1=self.min_manning_n, g2=self.max_manning_n )

    #   get_manning_n_grid()
    #-------------------------------------------------------------------------------
    def get_channel_sinuosity_grid(self):
    
        """Compute and save the estimated channel sinuosity grid"""
        #------------------------------------------------------------------------- 
        # There are different definitions of channel sinuosity.  Here we
        # are referring to the "absolute sinuosity", defined as the ratio
        # of the "along-channel flow distance" between the two endpoints
        # of a channel and the "straight-line distance" between those endpoints.
        # 
        # By this definition, sinuosity is dimensionless \[km/km\], with a
        # minimum possible value of 1.0.  It tends to increase slowly from
        # 1 where TCA is small to a larger value where TCA is big, but
        # typically does not exceed 1.3.
        # 
        # The following code uses a power-law estimate of the form:
        # $s = c \, A^p$, where A is the total contributing area (TCA)
        # grid, to create a grid of sinuosity values.  The value that
        # s should have where A is maximum (e.g. the river outlet) is
        # set as g1.  Similarly, the value that n should have
        # where A is minimum (e.g. near a ridge) is set as g2.
        # The coefficient, c, and power, p, are then set to match these
        # constraints.
        #-------------------------------------------------------------------------
        # Note:  max_sinuosity should be
        #        specified in the set_basin_info() method.
        #-------------------------------------------------------------------------
        parameterize.get_grid_from_TCA(site_prefix=self.site_prefix,
                     topo_dir=self.topo_dir,
                     area_file=self.d8_area_file, out_file=self.sinu_file,
                     g1=self.max_sinuosity, g2=self.min_sinuosity )

    #   get_channel_sinuosity_grid()
    #-------------------------------------------------------------------------------
    def get_bankfull_depth_grid(self):
    
        """Compute and save the estimated bankfull depth grid"""
        #-----------------------------------------------------------------
        # The "bankfull depth" is the maximum in-channel water depth of
        # a river at a given location.  (It varies throughout a river
        # basin.)  When the depth of water in a river exceeds this depth,
        # "overbank flow" occurs and water enters the flood plain adjacent
        # to the channel. "Overbank flow depth", "inundation depth" or
        # simply "flooding depth" are terms that refer to the depth of
        # water on land outside of the river channel.  It is important
        # to know the bankfull depth in order to more accurately predict
        # the flooding depth. 
        # 
        # While remote sensing images can be used to estimate a river's
        # bankfull width, the river bed typically cannot be "seen"
        # through the water.  Moreover, bankfull depth is typically only
        # measured at a few locations (e.g. at gauging stations) within
        # a river basin, so accurate values of bankfull depth are
        # difficult to obtain.
        # 
        # The following code uses a power-law estimate of the form:
        # $d_b = c \, A^p$, where A is the total contributing area (TCA)
        # grid, to create a grid of bankfull depth values.  The value
        # that $d_b$ should have where A is maximum (e.g. the river
        # outlet) is set as g1.  A typical, empirical value for p is 0.4.
        # The coefficient, c, is then set to match these constraints.
        #-----------------------------------------------------------------
        # Note:  max_bankfull_depth and bankfull_depth_power should be
        #        specified in the set_basin_info() method.
        #-----------------------------------------------------------------        
        parameterize.get_grid_from_TCA(site_prefix=self.site_prefix,
                topo_dir=self.topo_dir,
                area_file=self.d8_area_file,  out_file=self.dbank_file,
                g1=self.max_bankfull_depth, p=self.bankfull_depth_power )

    #   get_bankfull_depth_grid()
    #-------------------------------------------------------------------------------
    def get_init_channel_depth_grid(self):
    
        """Compute the estimated initial channel water depth grid"""
        #----------------------------------------------------------------- 
        # Here we attempt to estimate the initial depth of water for
        # every channel in the river network.  This is supposed to be
        # the "normal depth" of the river that is maintained by baseflow
        # from groundwater (i.e. due to the groundwater table intersecting
        # the channel bed) and is not attributed to a recent rainfall
        # event. This is the starting or initial condition for a model run.
        # 
        # This routine uses a grid-based Newton-Raphson iterative scheme
        # to solve a transcendental equation (see Wikipedia) for the
        # initial depth of water in a channel network that results from
        # groundwater baseflow.  The variables involved are:
        # 
        # w = bed bottom width, trapezoid [m]<br>
        # A = upstream area [$km^2$]<br>
        # S = downstream slope [m/m]<br>
        # n = Manning roughness parameter  [$s/m^{1/3}$]<br>
        # $\theta$ = bank angle [degrees]<br>
        # d = water depth in channel [m]<br>
        # $A_c$ = wetted cross-section area [$m^2$]<br>
        # P  = wetted cross-section perimeter [m]<br>
        # $R_h = (A_c / P)$ = hydraulic radius [m]<br>
        # B = spatially-uniform baseflow volume flux [$m s^{-1}$]<br>
        # 
        # The equations used here are: <br>
        # $Q = v \, A_c = B \,A$    [$m^3 s^{-1}$] (steady-state) <br>
        # $v = (1/n) \, {R_h}^{2/3} \, S^{1/2} \,\,\,$  [SI units] <br>
        # $R_h = A_c / P$ <br>
        # $A_c = d \, [w + (d \, \tan(\theta))]$ <br>
        # $P = w + [2 \, d \, / \cos(\theta)]$ <br>
        # 
        # Note that B can be estimated from a baseflow discharge
        # measured at the basin outlet.
        # 
        # If we are given w, n, theta, A, S and B, then we get an
        # equation for d that cannot be solved in closed form.
        # However, we can write the equation $v \, A_c = B \, A$
        # in the form needed to solve for d (in every grid cell)
        # by Newton's method, i.e.:
        # $F(d) = [v(d) \, A_c(d)] - (B \, A) = 0$.
        #-----------------------------------------------------------------
        # Note:  A_out_km2 (TCA at basin outlet) and Qbase_out
        #        (estimated baseflow discharge at basin outlet,
        #        m3/s, should be specified in set_basin_info() method.
        #-----------------------------------------------------------------        
        B_mps = init_depth.get_baseflow_volume_flux( self.A_out_km2,
                               self.Qbase_out, REPORT=True)
        init_depth.compute_initial_depth( site_prefix=self.site_prefix,
                   topo_dir=self.topo_dir,
                   SILENT=False, baseflow_rate=B_mps,
                   bank_angle=self.bank_angle,    #################
                   # angle_file=self.angle_file,
                   area_file=self.d8_area_file,
                   slope_file=self.d8_slope_file,
                   width_file=self.width_file,
                   manning_file=self.manning_file,
                   sinu_file=self.sinu_file,
                   d0_file=self.d0_file)

    #   get_init_channel_depth_grid()
    #-------------------------------------------------------------------------------
    def get_soil_hydraulic_grids(self):

        """Compute and save soil hydraulic property grids (via pedotransfer)"""
        #----------------------------------------------------------------- 
        # Here, we first read a set of ISRIC soil property grids for
        # each of 7 soil layers (see above) in GeoTIFF format.  We then
        # clip them to a chosen geographic bounding box and resample
        #  (or regrid) them to a new spatial resolution (grid cell size).
        # (Typically, we regrid to the TopoFlow model grid.)  Finally,
        # we use a
        # [<b>pedotransfer functions</b>]
        # (https://en.wikipedia.org/wiki/Pedotransfer_function)
        # (Wosten et al., 1998, 2001 ) to compute a corresponding set
        # of soil hydraulic property grids that are used to compute
        # infiltration. (This uses pedotransfer.py in topoflow/utils.
        # These grids are referenced in the CFG files for the TopoFlow
        # infiltration components.  The REPORT flag can be set to True
        # to see more detailed information about each of the soil
        # property grids.  Any warnings or errors are printed regardless.
        # In some cases, out-of-range values are generated and by default
        # are forced into the valid range (from above or below).
        # Spurious values can cause the infiltration model to become
        # [<b>numerically unstable</b>]
        # (https://en.wikipedia.org/wiki/Numerical_stability).
        # 
        # Files in the <b>src_soil_dir</b> folder span the entire
        # country of Ethiopia and were downloaded from the
        # [<b>ISRIC SoilGrids website</b>](https://soilgrids.org/).
        # ISRIC SoilGrids provides <b>global</b> coverage at two different
        # spatial resolutions:  250m and 1km.  (The products actually use
        # Geographic coordinates and resolutions are 7.5 and 30 arcseconds.
        # The north-south dimension of grid cells therefore varies with
        # latitude, being roughly 250m or 1km at the equator, but smaller
        # at other latitudes.)  They contain soil property variables for
        # each of 7 soil layers, as follows.
        # 
        # Variables:
        # BLDFIE = Bulk density [kg / m3]
        # BDRICM = Absolute depth to bedrock [cm]
        # CLYPPT = Mass fraction of clay [%]
        # ORCDRC = Soil organic carbon content (fine earth fraction)  [g / kg]
        # SLTPPT = Mass fraction of silt [%]
        # SNDPPT = Mass fraction of sand [%]
        # 
        # Layer 1 = sl1 = 0.00 to 0.05 m   (0  to  5 cm)
        # Layer 2 = sl2 = 0.05 to 0.15 m   (5  to 15 cm)
        # Layer 3 = sl3 = 0.15 to 0.30 m   (15 to 30 cm)
        # Layer 4 = sl4 = 0.30 to 0.60 m   (30 to 60 cm)
        # Layer 5 = sl5 = 0.60 to 1.00 m   (60 to 100 cm)
        # Layer 6 = sl6 = 1.00 to 2.00 m   (100 to 200 cm)
        # Layer 7 = sl7 = 2.00 to ??   m   (200 to ??? cm)
        # 
        # ISRIC Data Citation
        # Hengl T, Mendes de Jesus J, Heuvelink GBM, Ruiperez Gonzalez M,
        # Kilibarda M, Blagotić A, et al. (2017) SoilGrids250m:
        # Global gridded soil information based on machine learning.
        # PLoS ONE 12(2): e0169748. doi:10.1371/journal.pone.0169748
        #-----------------------------------------------------------------
        
        #--------------------------------------------
        # Copy RTI file from topo to soil directory
        #--------------------------------------------
        topo_rti_file = self.rti_file
        soil_rti_file = self.soil_dir + self.site_prefix + '.rti'
        shutil.copyfile( topo_rti_file, soil_rti_file )
        ### REPORT = not(self.SILENT)
        if not(self.SILENT):
            print('Clipping/regridding soil hydraulic property grids...')

        pedotransfer.save_soil_hydraulic_vars( site_prefix=self.site_prefix,
                 in_dir=self.src_soil_dir, out_dir=self.soil_dir,
                 out_bounds=self.out_bounds, REPORT=False,
                 out_xres_sec=self.out_xres_sec,
                 out_yres_sec=self.out_yres_sec,
                 RESAMPLE_ALGO='bilinear')

    #   get_soil_hydraulic_grids()
    #-------------------------------------------------------------------------------
    def write_simple_rain_file(self):
    
        """Write a simple, uniform "rainrates file" to "met" directory"""
 
        # In this step we create a text file with a simple time series of
        # spatially uniform rainrates for testing purposes.
        # The get_rainfall_grid_stack() method can be used to create
        # space-time rainfall grid stacks (e.g. CHIRPS).

        self.test_rain_file = self.met_dir + self.case_prefix + '_rain_rates.txt'
        test_rain_unit = open( self.test_rain_file, 'w')
        for k in range(5):
            test_rain_unit.write('3.0\n')
        for k in range(15):
            test_rain_unit.write('10.0\n')
        for k in range(5):
            test_rain_unit.write('3.0\n')
        test_rain_unit.close()
        
    #   write_simple_rain_file()
    #-------------------------------------------------------------------------------
    def get_rainfall_grid_stack(self):
    
        """Create a space-time rainfall rate grid stack"""
        #-------------------------------------------------------------------------------
        # Note: Downloading all of the required files can take hours to days.
        # However, assuming all required files have already been downloaded
        # to your computer, this can also be performed using TopoFlow's
        # regrid utility.  The specific functions for this are called:
        # create_rts_from_nc_files() and create_rts_from_chirps_files().
        # 
        # For the [<b>MINT Project</b>](http://mint-project.info/), we
        # experimented with three different global, remote-sensing-based
        # products for rainfall, called:  GPM, GLDAS and CHIRPS.
        # 
        # [Version 06 of the GPM GPM IMERG Final Precipitation L3]
        # (https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGHH_06/summary?keywords=GPM)
        # has a spatial resolution of <b>0.1 x 0.1 degrees</b> (360 x 360 arcseconds)
        # and a temporal resolution of <b>30 minutes</b> (1800 seconds). 
        # This data is available for dates between: <b>2000-06-01</b> and <b>2020-01-31</b>.
        # This is the highest spatial and temporal resolution available for GPM data.
        # 
        # [Version 2.1 of the GLDAS Noah Land Surface Model L4]
        # (https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_3H_2.1/summary?keywords=GLDAS)
        # has a spatial resolution of <b>0.25 x 0.25 degrees</b> (900 x 900 arcseconds)
        # and a temporal resolution of <b>3 hours</b> (180 minutes).
        # This data is available for dates between: <b>2000-01-01</b> and <b>2020-02-29</b>.
        # This is the highest spatial and temporal resolution available for GLDAS data.
        # 
        # [Version 2.0 of the UCSB CHIRPS Rainfall Product]
        # (https://data.chc.ucsb.edu/products/CHIRPS-2.0/).
        # CHIRPS data (Climate Hazards Group InfraRed Precipitation with Station) is
        # available for various spatial and temporal resolutions as well as regions.
        # Here we use the Africa 6-hourly data set, which
        # has a spatial resolution of <b>0.1 x 0.1 degrees</b> (360 x 360 arcseconds)
        # and a temporal resolution of <b>6 hours</b> (21600 seconds). 
        # This data is available for dates between: <b>1981-01</b> and <b>2021-05</b>.
        # (Note: Dates beyond 2015-12 are in the extra_step folder.)
        # See the document README-CHIRPS.txt for more info.
        # The FTP site is:
        # (ftp://chg-ftpout.geog.ucsb.edu/pub/org/chg/products/CHIRPS-2.0/).
        #------------------------------------------------------------------------------- 
        # In addition to space-time rainfall rates, some of the TopoFlow
        # components can make use of other meteorological variables such as:
        # air temperature, soil temperature, relative humidity, surface wind
        # speed, shortware radiation and longwave radiation.  The procedure
        # for preparing these other space-time datasets for use in TopoFlow
        # is very similar to what is done for precipitation.
        #-------------------------------------------------------------------------------  
        rti_file  = self.topo_dir + self.site_prefix + '.rti'
        grid_info = rti_files.read_info( rti_file )  # needed for DEM ncols & nrows
        self.rain_rts_path = self.met_dir + self.rain_rts_file  ######
        
        if (self.CHIRPS or self.CHIRPS2):
            resample_algo = 'bilinear'
            ## resample_algo = 'nearest'
            time_interval_hours = 6.0    # For the 6-hourly rainfall product
            non_negative = True          # zero out any negative values (e.g. nodata)

            os.chdir( self.src_rain_dir )
            regrid.create_rts_from_chirps_files( rts_file=self.rain_rts_path,
                   time_interval_hours=time_interval_hours,
                   resample_algo=resample_algo, NON_NEGATIVE=non_negative,
                   VERBOSE=False, DEM_bounds=self.out_bounds,
                   DEM_xres_sec=self.out_xres_sec,
                   DEM_yres_sec=self.out_yres_sec,
                   DEM_ncols=grid_info.ncols, DEM_nrows=grid_info.nrows )

        if (self.GLDAS or self.GPM):
            NC4 = True
            resample_algo = 'bilinear'
            ## resample_algo = 'nearest'
            if (self.GLDAS):
                time_interval_hours = 3.0   # not used yet
            else:
                time_interval_hours = 0.5   # not used yet
    
            os.chdir( self.src_rain_dir )
            regrid.create_rts_from_nc_files( rts_file=self.rain_rts_path, NC4=NC4,
                   GPM=self.GPM, VERBOSE=False, DEM_bounds=self.out_bounds,
                   DEM_xres_sec=self.out_xres_sec,
                   DEM_yres_sec=self.out_yres_sec,
                   DEM_ncols=grid_info.ncols, DEM_nrows=grid_info.nrows,
                   resample_algo=resample_algo)

    #   get_rainfall_grid_stack()
    #-------------------------------------------------------------------------------
    def subset_ngen_csv_file(self):
    
        """Subset a NextGen CSV file using a specified date range"""
        sp = self.site_prefix   # (e.g. 'cat-84')
        in_file  = self.src_ngen_met_dir + sp + '.csv'
        out_file = self.met_dir + self.rain_csv_file
        ## out_file = self.met_dir + sp + '_' + self.date_range_str + '.csv'
        #--------------------------------------
        csv_unit     = open( in_file, 'r')
        new_csv_unit = open( out_file, 'w')

        #----------------------------------------        
        # Read header line and copy to new file
        #----------------------------------------
        line = csv_unit.readline()
        new_csv_unit.write( line )
 
        while (True):
            line = csv_unit.readline()
            if (line.startswith( self.start_date_str)):
                new_csv_unit.write( line  )
                break

        while (True):
             line = csv_unit.readline()
             if (line.startswith( self.end_date_str)):
                # Don't write out this line
                break
             new_csv_unit.write( line )
             
        #-------------------                
        # Close both files
        #-------------------
        csv_unit.close()
        new_csv_unit.close()
              
    #   subset_ngen_csv_file()   
    #-------------------------------------------------------------------------------
    
    
    