
## Copyright (c) 2009-2012, Scott D. Peckham
##
## October 2012   (CSDMS Standard Names and BMI)
## January 2009  (converted from IDL)
## December 2009 (starting from evap_priestley_taylor.py)
## May 2010 (changes to unit_test() and read_cfg_file()

#################################################
#  NB!  See note at "THIS MAY BE COSTLY"
#################################################

#-----------------------------------------------------------------------
#  NOTES:  This file defines a "Read from file" ET component
#          and related functions.  It inherits from the ET
#          "base class" in "evap_base.py".
#-----------------------------------------------------------------------
#
#  class evap_component
#
#      get_component_name()
#      get_attribute()           # (10/26/11)
#      get_input_var_names()     # (10/23/12)
#      get_output_var_names()    # (10/23/12)
#      get_var_name()            # (10/23/12)
#      get_var_units()           # (10/23/12)
#---------------------------
#      check_input_types()
#      update_ET_rate()
#      update_water_balance()   ###### Overrides one in evap_base.py
#      open_input_files()
#      read_input_files()
#      close_input_files()
#
#  Functions:
#      RTG_to_RTS()
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.components import evap_base

from topoflow.utils import model_input
from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
class evap_component( evap_base.evap_component ):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Evaporation_Read_File',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'EvapReadFile',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Evap_Read_File.cfg.in',
        'cfg_extension':      '_evap_read_file.cfg',
        'cmt_var_prefix':     '/EvapReadFile/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Evap_Read_File.xml',
        'dialog_title':       'Evaporation: Read File Parameters',
        'time_units':         'seconds' }

    _input_var_names = [
        'channel_water_x-section__mean_depth' ]    # channels comp

    _output_var_names = [
        'land_surface_water__evaporation_volume_flux',
        'land_surface_water__domain_time_integral_of_evaporation_volume_flux',
        'model_grid_cell__area',
        'model__time_step' ]

    _var_name_map = {
        'channel_water_x-section__mean_depth' : 'd',
        'land_surface_water__evaporation_volume_flux': 'ET',
        'land_surface_water__domain_time_integral_of_evaporation_volume_flux': 'vol_ET',
        'model_grid_cell__area': 'da',
        'model__time_step': 'dt' }

    _var_units_map = {
        'channel_water_x-section__mean_depth' : 'm',
        'land_surface_water__evaporation_volume_flux': 'm s-1',
        'land_surface_water__domain_time_integral_of_evaporation_volume_flux': 'm3',
        'model_grid_cell__area': 'm2',
        'model__time_step': 's' }

    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Evaporation_Read_File'

    #   get_component_name()     
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print('###################################################')
            print(' ERROR: Could not find attribute: ' + att_name)
            print('###################################################')
            print(' ')

    #   get_attribute()
    #-------------------------------------------------------------------
    def get_input_var_names(self):

        #--------------------------------------------------------
        # Note: These are currently variables needed from other
        #       components vs. those read from files or GUI.
        #--------------------------------------------------------   
        return self._input_var_names
    
    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):
 
        return self._output_var_names
    
    #   get_output_var_names()
    #-------------------------------------------------------------------
    def get_var_name(self, long_var_name):
            
        return self._var_name_map[ long_var_name ]

    #   get_var_name()
    #-------------------------------------------------------------------
    def get_var_units(self, long_var_name):

        return self._var_units_map[ long_var_name ]
   
    #   get_var_units()
    #-------------------------------------------------------------------
##    def get_var_type(self, long_var_name):
##
##        #---------------------------------------
##        # So far, all vars have type "double",
##        # but use the one in BMI_base instead.
##        #---------------------------------------
##        return 'float64'
##    
##    #   get_var_type()
    #-------------------------------------------------------------------
##    def check_input_types(self):
##
##        #----------------------------------------------------
##        # Notes: Usually this will be overridden by a given
##        #        method of computing ET.
##        #----------------------------------------------------
##        are_scalars = np.array([
##                         self.is_scalar('T_soil_x'),
##                         self.is_scalar('soil_x'),
##                         self.is_scalar('K_soil'),
##                         self.is_scalar('alpha'),
##                         #-------------------------------
##                         self.mp.is_scalar('Qn_SW'),
##                         self.mp.is_scalar('Qn_LW'),
##                         self.mp.is_scalar('T_air'),
##                         self.mp.is_scalar('T_surf'),
##                         #-------------------------------
##                         self.cp.is_scalar('d'),
##                         #-------------------------------
##                         self.gp.is_scalar('h_table') ])
##
##        self.ALL_SCALARS = np.all(are_scalars)
##        
##    #   check_input_types()
    #-------------------------------------------------------------------
    def update_ET_rate(self):

        #----------------------------------------------------------
        # Note:  This function should return without doing
        #        anything.  For this component, the ET rate
        #        is not computed but is read directly as the
        #        only "input" variable.  Note that "initialize()"
        #        calls "read_input_files()" so that ET has been
        #        defined before first call to evap_base.update().
        #        Recall that evap_base.update() will call
        #        read_input_files() to get next ET grid.
        #----------------------------------------------------------
        #        We may want to make sure ET values read from
        #        input file are positive before using them.
        #        If so, uncomment 2 lines below.
        #----------------------------------------------------------
        return 
##        if (self.ET is not None):
##            self.ET = np.maximum(self.ET, np.float64(0))

        ##########################################
        #  THIS MAY BE COSTLY.  BETTER WAY OR
        #  ALLOW ET TO BE A SCALAR ??
        ##########################################
##        if (size(self.ET) == 1):
##            self.ET += np.zeros((self.ny, self.nx), dtype='Float64')
            
    #   update_ET_rate()
    #-------------------------------------------------------------------  
    def update_water_balance(self):

        #----------------------------------------------------------
        # Notes: Return without modifying the surface water depth
        #        or soil moisture in the top soil layer.
        #        Depending on how the model is run, this may
        #        lead to an inconsistent result where the mass
        #        of water in the system is not conserved.
        #----------------------------------------------------------
        print('-----------------------------------------------------')
        print('WARNING: ET_read_file component assumes that the')
        print('         ET rates read from files are actual vs.')
        print('         potential rates.  It does not consume')
        print('         surface or subsurface water.')
        print('-----------------------------------------------------')
        return
    
    #   update_water_balance()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        self.ET_file = self.in_directory + self.ET_file

        self.ET_unit = model_input.open_file(self.ET_type, \
                                             self.ET_file)
        print('self.ET_type =', self.ET_type)
        print('self.ET_file =', self.ET_file)
        
##        self.duration_unit  = model_input.open_file(self.duration_type,  \
##                                                    self.duration_file)
##        print 'self.duration_type (ET) =', self.duration_type
##        print 'self.duration_file (ET) =', self.duration_file

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        #--------------
        # For testing
        #--------------
        ## print 'self.ET_type =', self.ET_type
        ## print 'self.ET_file =', self.ET_file
        
        rti = self.rti
        
        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        ET = model_input.read_next(self.ET_unit, self.ET_type, rti,
                                   factor=self.mmph_to_mps)
        if (ET is not None):
            self.ET = ET
            print('min(ET) =', ET.min() * self.mps_to_mmph, ' [mmph]')
            print('min(ET) =', ET.max() * self.mps_to_mmph, ' [mmph]')
##            print 'min(ET) =', ET.min(), ' [mps]'
##            print 'max(ET) =', ET.max(), ' [mps]'
            
    #   read_input_files()
    #-------------------------------------------------------------------
    def close_input_files(self):
 
        try:
            self.ET_unit.close()
        except:
            print('Could not close input file "ET_unit".')
            
        ## if (self.ET_file != ''): self.ET_unit.close()
        
    #   close_input_files()
    #-------------------------------------------------------------------
#-----------------------------------------------------------------------
def RTG_to_RTS( RTG_prefix, new_RTS_file,
                nx=None, ny=None, dtype='float32' ):

    #--------------------------------------------------------------
    # Notes:  This function "bundles" a set of binary grid files
    #         (generic RTG format) with the same file_name prefix
    #         and embedded time index into a single RTS file.
    #--------------------------------------------------------------
    import glob
    import os.path
    import numpy
    import rti_files

    #--------------------
    # Open new RTS file
    #--------------------
    RTS_EXISTS = os.path.exists( new_RTS_file )
    if (RTS_EXISTS):
        print('SORRY, An RTS file with the name')
        print(new_RTS_file)
        print('already exists.')
        return

    RTG_list = glob.glob( RTG_prefix + '*' )
    RTG_list = np.sort( RTG_list )
    print('Number of RTG files to merge =', len(RTG_list))
    
    #--------------------------------
    # Try to get info from RTI file
    #--------------------------------
    if (nx == None) and (ny == None):
        RTG_file1 = RTG_list[0]
        RTI_file  = rti_files.try_to_find_rti_file( RTG_file1 )
        if (RTI_file != 'none'):
            info = rti_files.read_info( RTI_file )
            nx   = info.ncols
            ny   = info.nrows
        else:
            print('ERROR: Could not find RTI file and nx and ny not provided.')
            return
##    else:
##        #--------------------------------------------
##        # For case when used as a script under Unix
##        # (used in the file rtg2rts.py)
##        #--------------------------------------------
##        nx = eval(nx)
##        ny = eval(ny)
        
    RTS_unit = open( new_RTS_file, 'wb' )

    #------------------------------
    # Write RTG grids to RTS file
    #------------------------------    
    for RTG_file in RTG_list:
        print('Reading values from:', RTG_file)
        RTG_unit = open( RTG_file, 'rb')
        grid     = np.fromfile( RTG_unit, count=nx*ny, dtype=dtype )
        RTG_unit.close()
        grid.tofile( RTS_unit )

    #---------------------
    # Close new RTS file
    #---------------------
    RTS_unit.close()
    print('Finished writing RTG files to RTS format.')
    print(' ')
    
#   RTG_to_RTS()
#-----------------------------------------------------------------------
    


