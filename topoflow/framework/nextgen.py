#! /usr/bin/env python
"""
This module initializes settings for running TopoFlow in original NextGen.
"""
#      
#-----------------------------------------------------------------------      
# Copyright (c) 2022-2025, Scott D. Peckham
#
# 2025-10.  Moved code here from topoflow/framework/multi_bmi.py
#           to reduce impact of NextGen changes on multi_bmi.py.
# 2024-09.  In multi_bmi.py, renamed prepare_tf_inputs() to
#           prepare_inputs_for_ngen_basin() and added NGEN_CSV
#           keyword to prepare_inputs_for_ngen_basin().
# 2022-10.  Wrote multi_bmi.py starting from emeli.py.
#
#-----------------------------------------------------------------------
#
# initialize_ngiab()
# initialize_nextgen()
# prepare_tf_inputs_for_ngen()
#
#-----------------------------------------------------------------------

import os, os.path, sys
from topoflow.utils import prepare_inputs as prep

#-----------------------------------------------------------------------
def initialize_ngiab(self):

    #------------------------------------------------------------
    # Note:  For a NGIAB Docker image, the directory names are
    #        standardized and don't need to be set in CFG file.
    #        These include ngen_dir, spatial_dir, forcing_dir,
    #        and rc_dir.
    #------------------------------------------------------------
    self.cfg_prefix    = filename.replace('_ngiab.cfg', '')
    self.cfg_directory = self.cfg_file.replace(filename, '')
    self._att_map['cfg_extension'] = '_ngiab.cfg'
    ## self.case_prefix  = self.cfg_prefix
    
    self.read_config_file()

    #--------------------------------------------------    
    # Set the "ngen" directory for NGIAB Docker image
    #--------------------------------------------------
    self.ngen_dir = '/ngen/ngen/'

    #------------------------------------------------------
    # Set the "spatial" dir with hydrofabric GPKG file(s)
    #------------------------------------------------------
    self.spatial_dir = (self.ngen_dir + '/data/config/')
    #------------------------------------------------------    
    self.catchment_gpkg_file = self.catchment_file
    self.nexus_gpkg_file     = self.nexus_file
    ###################################################
    
    #------------------------------------------------
    # Set the "rc" dir with realization config file
    #------------------------------------------------
    self.rc_dir = (self.ngen_dir + '/data/config/')
      
    #---------------------------------------------
    # Set the "forcing" dir with forcing file(s)    
    #---------------------------------------------
    self.forcing_dir = (self.ngen_dir + '/data/forcing/')
    #-------------------------------------------------------
    self.NGEN_CSV = False  # (NetCDF is always used?)

    #--------------------------------------------------
    # Expand things like ".." and "~" in dir names ??
    #--------------------------------------------------
    self.ngen_dir    = os.path.expanduser( self.ngen_dir )
    self.spatial_dir = os.path.expanduser( self.spatial_dir )
    self.forcing_dir = os.path.expanduser( self.forcing_dir )
    self.rc_dir      = os.path.expanduser( self.rc_dir )       

#   initialize_ngiab()
#-----------------------------------------------------------------------
def initialize_nextgen(self):

    #------------------------------------------------------------
    # Note:  For the original NextGen, directory names are not
    #        standardized and must be set in the CFG file.
    #        These include ngen_dir, spatial_dir, forcing_dir,
    #        and rc_dir.
    #------------------------------------------------------------
    self.cfg_prefix    = filename.replace('_ngiab.cfg', '')
    self.cfg_directory = self.cfg_file.replace(filename, '')
    self._att_map['cfg_extension'] = '_ngiab.cfg'
    ## self.case_prefix  = self.cfg_prefix
    
    self.read_config_file()

    #--------------------------------------------------    
    # Set the "ngen" directory for NGIAB Docker image
    #--------------------------------------------------
    self.ngen_dir = '/ngen/ngen/'

    #------------------------------------------------------
    # Set the "spatial" dir with hydrofabric GPKG file(s)
    #------------------------------------------------------
    if (self.spatial_dir[:5] == 'ngen/'):
        self.spatial_dir = self.spatial_dir[5:]
    if (self.spatial_dir[-1] != '/'):
        self.spatial_dir += '/'           
    self.spatial_dir = (self.ngen_dir + self.spatial_dir)
    #------------------------------------------------------    
    self.catchment_gpkg_file = self.catchment_file
    self.nexus_gpkg_file     = self.nexus_file
    ###################################################
    
    #------------------------------------------------
    # Set the "rc" dir with realization config file
    #------------------------------------------------
    if (self.rc_dir[:5] == 'ngen/'):
        self.rc_dir = self.rc_dir[5:]
    if (self.rc_dir[-1] != '/'):
        self.rc_dir += '/'  
    self.rc_dir = (self.ngen_dir + self.rc_dir)
      
    #---------------------------------------------
    # Set the "forcing" dir with forcing file(s)    
    #---------------------------------------------
    if (self.forcing_dir[:5] == 'ngen/'):
        self.forcing_dir = self.forcing_dir[5:]
    if (self.forcing_dir[-1] != '/'):
        self.forcing_dir += '/'  
    self.forcing_dir = (self.ngen_dir + self.forcing_dir)
    #-------------------------------------------------------
    self.NGEN_CSV = (self.forcing_type == 'local_aorc_csv')  # OR ngen_aorc_csv

    #--------------------------------------------------
    # Expand things like ".." and "~" in dir names ??
    #--------------------------------------------------
    self.ngen_dir    = os.path.expanduser( self.ngen_dir )
    self.spatial_dir = os.path.expanduser( self.spatial_dir )
    self.forcing_dir = os.path.expanduser( self.forcing_dir )
    self.rc_dir      = os.path.expanduser( self.rc_dir )       

#   initialize_nextgen()
#-----------------------------------------------------------------------
def prepare_tf_inputs_for_ngen(self):

    #--------------------------------------------
    # Note: DEM must already exist as GeoTIFF.
    # This is now called from initialize().
    #--------------------------------------------------------
    # site_prefix & case_prefix are read from path_info.cfg
    # file and set into self by initialize_config_vars(),
    # defined in BMI_base.py.
    #--------------------------------------------------------        
    p = prep.get_inputs()
    p.test_dir = self.ngen_dir + 'data/topoflow/input_files/'    ##########
    p.in_directory = self.in_directory
    #-----------------------------------------------
    # Copy values from this object into object "p"
    # These are now read from "_nextgen.cfg" file
    # via read_config_file() in initialize().
    #-----------------------------------------------
    p.ngen_dir    = self.ngen_dir
    p.NGEN_CSV    = self.NGEN_CSV  # Use NGEN CSV files?
    p.forcing_dir = self.forcing_dir
    p.spatial_dir = self.spatial_dir
    p.catchment_file = self.catchment_file  #######  GPKG
    p.nexus_file     = self.nexus_file      #######  GPKG
    #----------------------------------------------------                
    p.prepare_all_inputs( site_prefix=self.site_prefix,
         case_prefix=self.case_prefix,
         NGEN_DEM=True, CLIP_DEM=False,
         NO_SOIL=True, NO_MET=True)

    #---------------------------------------------------------
    # Next line is for when initialize calls read_grid_info.
    # Otherwise, self.topo_directory = self.in_directory.
    #---------------------------------------------------------         
    # self.topo_directory = p.topo_dir  # Not needed now.

#   prepare_tf_inputs_for_ngen()
#-------------------------------------------------------------------

