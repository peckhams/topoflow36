#
#  Copyright (c) 2001-2014, Scott D. Peckham
#
#  Sep 2014.  Cleanup and testing.
#             Own versions of input file routines, at end.
#  Aug 2014.  Customized initialize_computed_vars(). 
#             Updates to standard names and BMI.
#  Nov 2013.  Converted TopoFlow to Python package.
#  Jan 2013.  Revised handling of input/output names.
#  Oct 2012.  CSDMS Standard Names and BMI.
#  Aug 2009.  Updates.
#  Jul 2009.  Updates.
#  May 2010.  Changes to unit_test() and read_cfg_file().
#  Jan 2009.  Converted from IDL.
#
#-----------------------------------------------------------------------
#  NOTES:  This file defines a Degree-Day snowmelt component
#          and related functions.  It inherits from the snowmelt
#          "base class" in "snow_base.py".
#-----------------------------------------------------------------------
#
#  class snow_degree_day
#
#      get_component_name()
#      get_attribute()          # (10/26/11)
#      get_input_var_names()    # (5/14/12)
#      get_output_var_names()   # (5/14/12)
#      get_var_name()           # (5/16/12, Bolton)
#      get_var_units()          # (5/16/12, Bolton)
#      ----------------------
#      check_input_types()
#      update_meltrate()
#      ----------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#
#-----------------------------------------------------------------------

import numpy as np

from topoflow.utils import model_input
from topoflow.components import snow_base

#-----------------------------------------------------------------------
class snow_component( snow_base.snow_component ):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Snowmelt_Degree_Day',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'SnowDegreeDay',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Snow_Degree_Day.cfg.in',
        'cfg_extension':      '_snow_degree_day.cfg',
        'cmt_var_prefix':     '/SnowDegreeDay/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Snow_Degree_Day.xml',
        'dialog_title':       'Snowmelt: Degree-Day Parameters',    # (Snowmelt ?)
        'time_units':         'seconds' }

    #----------------------------------------------------------
    # The only input variable that the Degree-Day method
    # needs from the Met component is "T_air".  The others
    # are used for consistent use of scalars vs. grids;
    # see "check_input_types()" below.
    #----------------------------------------------------------
    # Some input vars, like c0 and T0 are read from the CFG
    # file and not from other components.  They are therefore
    # included with the output_vars.
    #------------------------------------------------------------
    # Note: snowfall_leq-volume_flux *must* be nonliquid precip.
    #------------------------------------------------------------   
    _input_var_names = [
        'atmosphere_bottom_air__temperature',
        'atmosphere_water__snowfall_leq-volume_flux',
        'land_surface__temperature',        
        'water-liquid__mass-per-volume_density' ]
        #------------------------------------------------------------        
        # None of these should be needed for the Degree-Day method.
        #------------------------------------------------------------
        # 'atmosphere_bottom_air__mass-per-volume_density',
        # 'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity',
        # 'atmosphere_bottom_air__pressure',  # (not needed)
        # 'atmosphere_bottom_air_flow__log_law_roughness_length',
        # 'atmosphere_bottom_air_flow__speed_reference_height',
        # 'atmosphere_bottom_air_flow__reference-height_speed',        
        # 'atmosphere_bottom_air_water-vapor__relative_saturation',
                        
    #------------------------------------------------------------
    # Note: The main output of the Degree-Day method is SM.
    #       Functions in snow_base.py compute additional output
    #       vars, such as:
    #       update_SM_integral(), update_depth(), update_swe().
    #       They need things like: rho_H2O and rho_snow.
    #------------------------------------------------------------
    # Note: Cp_snow is a constant set in the "set_constants()"
    #       function in snow_base.py.
    #------------------------------------------------------------
    #       vol_SM was "basin_cumulative_snow_meltwater_volume"
    #------------------------------------------------------------
    _output_var_names = [
        'model__time_step',                                   # dt 
        'snowpack__degree-day_coefficient',                   # c0   (read from CFG)
        'snowpack__degree-day_threshold_temperature',         # T0   (read from CFG)
        'snowpack__domain_time_integral_of_melt_volume_flux', # vol_SM
        'snowpack__depth',                                    # h_snow
        'snowpack__initial_depth',                            # h0_snow
        'snowpack__initial_liquid-equivalent_depth',          # h0_swe
        'snowpack__liquid-equivalent_depth',                  # h_swe
        'snowpack__melt_volume_flux',                         # SM (MR is used for ice)
        'snowpack__z_mean_of_mass-per-volume_density' ]       # rho_snow
        #------------------------------------------------------------        
        # None of these should be needed for the Degree-Day method.
        #------------------------------------------------------------
        # 'snowpack__z_mean_of_mass-specific_isobaric_heat_capacity',       # Cp_snow 
                    
    _var_name_map = {
        'atmosphere_bottom_air__temperature':                 'T_air',
        'atmosphere_water__snowfall_leq-volume_flux':         'P_snow',
        'land_surface__temperature':                          'T_surf',        
        'water-liquid__mass-per-volume_density':              'rho_H2O',
        #------------------------------------------------------------------
        'model__time_step':                                   'dt',
        'snowpack__domain_time_integral_of_melt_volume_flux': 'vol_SM',
        'snowpack__degree-day_coefficient':                   'c0',
        'snowpack__degree-day_threshold_temperature':         'T0',
        'snowpack__depth':                                    'h_snow',
        'snowpack__initial_depth':                            'h0_snow',
        'snowpack__initial_liquid-equivalent_depth':          'h0_swe',
        'snowpack__liquid-equivalent_depth':                  'h_swe',
        'snowpack__melt_volume_flux':                         'SM',
        'snowpack__z_mean_of_mass-per-volume_density':        'rho_snow' }
        #------------------------------------------------------------        
        # None of these should be needed for the Degree-Day method.
        #------------------------------------------------------------
        # 'atmosphere_bottom_air__mass-per-volume_density':              'rho_air',
        # 'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'Cp_air',
        # 'atmosphere_bottom_air__pressure':                             'p0',
        # 'atmosphere_bottom_air_flow__log_law_roughness_length':        'z0_air',
        # 'atmosphere_bottom_air_flow__speed_reference_height':          'z',
        # 'atmosphere_bottom_air_flow__reference-height_speed':          'uz',
        # 'atmosphere_bottom_air_water-vapor__relative_saturation':      'RH',                
        # 'snowpack__z_mean_of_mass-specific_isobaric_heat_capacity':    'Cp_snow',
                    
    #-------------------------------------------------------------------
    # Note: We need to be careful with whether units are C or K,
    #       for all "thermal" quantities (e.g. thermal_capacity).
    #       Use [deg C] for actual temperatures and K everywhere else.
    #       This will usually be K-1, and a degree C or K is same size.
    #       Use "deg_C" because "C" = Coulombs in SI units.
    #-------------------------------------------------------------------       
    _var_units_map = {
        'atmosphere_bottom_air__temperature':                 'deg_C',  # (see Notes above)
        'atmosphere_water__snowfall_leq-volume_flux':         'm s-1',
        'land_surface__temperature':                          'deg_C',        
        'water-liquid__mass-per-volume_density':              'kg m-3',
        #----------------------------------------------------------------
        'model__time_step':                                   's',
        'snowpack__domain_time_integral_of_melt_volume_flux': 'm3',
        'snowpack__degree-day_coefficient':                   'mm day-1 K-1 ',
        'snowpack__degree-day_threshold_temperature':         'deg_C',
        'snowpack__depth':                                    'm',
        'snowpack__initial_depth':                            'm',
        'snowpack__initial_liquid-equivalent_depth':          'm',
        'snowpack__liquid-equivalent_depth':                  'm',
        'snowpack__melt_volume_flux':                         'm s-1',
        'snowpack__z_mean_of_mass-per-volume_density':        'kg m-3' }
        #------------------------------------------------------------        
        # None of these should be needed for the Degree-Day method.
        #------------------------------------------------------------
        # 'atmosphere_bottom_air__mass-per-volume_density':              'kg m-3',
        # 'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'J kg-1 K-1 ', # (see Notes above)
        # 'atmosphere_bottom_air__pressure':                             'mbar',
        # 'atmosphere_bottom_air_flow__log_law_roughness_length':        'm',
        # 'atmosphere_bottom_air_flow__speed_reference_height':          'm',
        # 'atmosphere_bottom_air_flow__reference-height_speed':          'm s-1',
        # 'atmosphere_bottom_air_water-vapor__relative_saturation':      '1',
        # 'snowpack__z_mean_of_mass-specific_isobaric_heat_capacity':    'J kg-1 K-1 ',
                                
    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )
 
    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Snow_Degree_Day'

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
    def check_input_types(self):

        #--------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are
        #        currently always scalars.
        #--------------------------------------------------        
        are_scalars = np.array([
                          self.is_scalar('P_snow'),
                          self.is_scalar('rho_H2O'),
                          self.is_scalar('T_air'),
                          ## self.is_scalar('rho_air'),  # (not needed)
                          ## self.is_scalar('Cp_air'),   # (not needed)
                          #-------------------------------
                          self.is_scalar('rho_snow'),
                          ## self.is_scalar('Cp_snow'),  # (not needed)
                          self.is_scalar('h0_snow'),
                          self.is_scalar('h0_swe'),
                          self.is_scalar('c0'),
                          self.is_scalar('T0') ])

        self.ALL_SCALARS = np.all(are_scalars)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def update_meltrate(self):

        #------------------------------------------------------------
        # Notes: This computes a "potential" meltrate, which can't
        #        be realized unless there is enough snow.
        #        See snow_base.enforce_max_meltrate().
        #------------------------------------------------------------
        # Notes: Arguments are assumed to be references to a scalar
        #        or grid so that:
        #        c0     = degree-day melt factor [mm/day/deg_C]
        #        T0     = threshold temperature [deg_C]
        #        M      = water equivalent of snowmelt [m/s]
        #        M_max  = max possible meltrate if all snow melts
        #        T_air  = air temperature [deg_C]

        #        Model must start when snow is isothermal.
        #        Cooling of the snowpack is not considered.

        #        This is a simple snowmelt model that can be used
        #        when there is insufficient data to use the energy
        #        balance method.  If c0 and T0 are established for a
        #        range of conditions, Kane et al. (1997) showed that
        #        this method gives comparable results to the energy
        #        balance method.

        #        86400000d = 1000 [mm/m] * 60 [sec/min] *
        #                    60 [min/hr] * 24 [hrs/day]

        #        rho_snow is not needed in formula here, but is
        #        needed to convert snowmelt to snow water
        #        equivalent (swe) in update_swe() in "snow_base.py".
        #-------------------------------------------------------------
        T_air = self.T_air  # (2/3/13, new framework)
        
        #------------------------------
        # Compute degree-day meltrate
        #------------------------------
        M = (self.c0 / np.float64(8.64E7)) * (T_air - self.T0)   #[m/s]

        # This is really an "enforce_min_meltrate()"
        self.SM = np.maximum(M, np.float64(0))
   
        #-------------------------------------------------------
        # Note: enforce_max_meltrate() method is always called
        #       by the base class to make sure that meltrate
        #       does not exceed the max possible.
        #-------------------------------------------------------
        #  self.enforce_max_meltrate()
        
    #   update_meltrate()
    #------------------------------------------------------------------- 
    def open_input_files(self):

        self.c0_file       = self.in_directory + self.c0_file
        self.T0_file       = self.in_directory + self.T0_file
        self.rho_snow_file = self.in_directory + self.rho_snow_file
        self.h0_snow_file  = self.in_directory + self.h0_snow_file
        self.h0_swe_file   = self.in_directory + self.h0_swe_file

        self.c0_unit       = model_input.open_file(self.c0_type,       self.c0_file)
        self.T0_unit       = model_input.open_file(self.T0_type,       self.T0_file)
        self.rho_snow_unit = model_input.open_file(self.rho_snow_type, self.rho_snow_file)
        self.h0_snow_unit  = model_input.open_file(self.h0_snow_type,  self.h0_snow_file)
        self.h0_swe_unit   = model_input.open_file(self.h0_swe_type,   self.h0_swe_file)

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        rti = self.rti

        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        c0 = model_input.read_next(self.c0_unit, self.c0_type, rti)
        if (c0 is not None): self.c0 = c0

        T0 = model_input.read_next(self.T0_unit, self.T0_type, rti)
        if (T0 is not None): self.T0 = T0

        rho_snow = model_input.read_next(self.rho_snow_unit, self.rho_snow_type, rti)
        if (rho_snow is not None): self.rho_snow = rho_snow

        h0_snow = model_input.read_next(self.h0_snow_unit, self.h0_snow_type, rti)
        if (h0_snow is not None): self.h0_snow = h0_snow
        
        h0_swe = model_input.read_next(self.h0_swe_unit, self.h0_swe_type, rti)
        if (h0_swe is not None): self.h0_swe = h0_swe
        
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):

        if (self.c0_type       != 'Scalar'): self.c0_unit.close()        
        if (self.T0_type       != 'Scalar'): self.T0_unit.close()
        if (self.rho_snow_type != 'Scalar'): self.rho_snow_unit.close()
        if (self.h0_snow_type  != 'Scalar'): self.h0_snow_unit.close()
        if (self.h0_swe_type   != 'Scalar'): self.h0_swe_unit.close()
        
##        if (self.c0_file       != ''): self.c0_unit.close()        
##        if (self.T0_file       != ''): self.T0_unit.close()
##        if (self.rho_snow_file != ''): self.rho_snow_unit.close()
##        if (self.h0_snow_file  != ''): self.h0_snow_unit.close()
##        if (self.h0_swe_file   != ''): self.h0_swe_unit.close()

    #   close_input_files()    
    #-------------------------------------------------------------------
