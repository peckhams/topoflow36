"""
This class defines a hydrologic evaporation component that implements
the Priestley-Taylor method.  It uses the same notation as the paper
Zhang et al. (2000).

This class inherits from the infiltration "base class" in "evap_base.py".

See: Zhang et al. (2000) Development and application of a spatially-
-distributed Arctic hydrological and thermal process model (ARHYTHM),
Hydrological Processes, 14, 1017-1044.
"""
#
#  NB!  If PRECIP_ONLY in met component, then Qn_SW and Qn_LW
#       both equal 0 and evaporation rate is zero.
#-----------------------------------------------------------------------
#
#  Copyright (c) 2001-2020, Scott D. Peckham
#
#  May 2020.  Testing with Jupyter notebook.
#             Get Q_net from met component directly?
#             In-place assignments in update_ET_rate().
#  Sep 2014.  New standard names and BMI updates and testing.
#  Aug 2014.  Updates to standard names and BMI.
#             Wrote latent_heat_of_evaporation(); not used yet.
#             Moved update_water_balance() to satzone_base.py.
#  Nov 2013.  Converted TopoFlow to Python package.
#  Feb 2013.  Adapted to use EMELI framework.
#  Jan 2013.  Revised handling of input/output names.
#  Oct 2012.  CSDMS Standard Names and BMI.
#  May 2010.  Changes to initialize() and read_cfg_file().
#  Aug 2009.  Updates.
#  Jul 2009.  Updates.
#  May 2009.  Updates.
#  Apr 2009.  Updates.
#  Jan 2009.  Converted from IDL to Python with I2PY.
#
#-----------------------------------------------------------------------
#  NOTES:  This file defines a Priestley-Taylor ET component
#          and related functions.  It inherits from the ET
#          "base class" in "evap_base.py".
#-----------------------------------------------------------------------
#
#  class evap_component
#
#      get_component_name()
#      get_attribute()             # (10/26/11)
#      get_input_var_names()       # (5/15/12)
#      get_output_var_names()      # (5/15/12)
#      get_var_name()              # (5/16/12), Bolton
#      get_var_units()             # (5/16/12), Bolton
#      ------------------------
#      check_input_types()
#      update_ET_rate()
#      ------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#
#  Functions:
#      Priestley_Taylor_ET_Rate   # Not used now.
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.components import evap_base

from topoflow.utils import model_input

#-----------------------------------------------------------------------
class evap_component( evap_base.evap_component ):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Evaporation_Priestley_Taylor',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #------------------------------------------------------
        'comp_name':          'EvapPriestleyTaylor',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Evap_Priestley_Taylor.cfg.in',
        'cfg_extension':      '_evap_priestley_taylor.cfg',
        'cmt_var_prefix':     '/EvapPriestleyTaylor/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Evap_Priestley_Taylor.xml',
        'dialog_title':       'Evaporation: Priestley-Taylor Parameters',
        'time_units':         'seconds' }

    #----------------------------------------------------------------
    # Note that the "meteorology" component uses the following to
    # compute Q_sum and Qe, but they aren't needed directly here:
    #     uz, z, z0_air, rho_air, Cp_air, Qn_SW, Qn_LW 
    #----------------------------------------------------------------
    _input_var_names = [
        'atmosphere_bottom_air__temperature',                # (meteorology)
        'land_surface_net-longwave-radiation__energy_flux',  # (meteorology)
        'land_surface_net-shortwave-radiation__energy_flux', # (meteorology)
        'land_surface__temperature' ]                        # (meteorology)
        #----------------------------------------------
        # These are no longer needed here. (9/25/14)
        #----------------------------------------------        
#         'channel_water_x-section__mean_depth',           # (@channels)       
#         'soil_top-layer__porosity',                      # (@satzone)
#         'soil_top-layer__saturated_thickness',           # (@satzone)
#         'soil_water_sat-zone_top_surface__elevation' ]   # (@satzone)

         #-------------------------------------------------
         # These are currently obtained from the GUI/file
         # and are not obtained from other components.
         #-------------------------------------------------
##        'land_surface__elevation',               # (GUI, DEM)
##        'land_surface_water__priestley-taylor_alpha_coefficient'  # (GUI, alpha)
##        'soil__reference_depth_temperature',     # (GUI, T_soil_x)
##        'soil_surface__temperature',             # (GUI, T_surf)
##        'soil__temperature_reference_depth',     # (GUI, soil_x)
##        'soil__thermal_conductivity' :           # (GUI, K_soil)
        
        #----------------------------------------------------
        # These could be added in the future; not used yet.
        #----------------------------------------------------
##        'soil_model_top_layer__saturated_water_content', # (satzone comp)
##        'land_surface_water_potential_evaporation_rate':'PET' }

    _output_var_names = [
        'land_surface_soil__conduction_heat_flux',     # (Qc)
        'land_surface_water__domain_time_integral_of_evaporation_volume_flux',  # (vol_ET)
        'land_surface_water__evaporation_volume_flux', # (ET)
        'model__time_step' ] # (dt)

        #-----------------------------------------------------
        # These are read from GUI/file, but can be returned.
        #-----------------------------------------------------       
        #'land_surface__elevation',
        #'land_surface_water__priestley-taylor_alpha_coefficient',
        #'soil__reference_depth_temperature',
        ## 'soil_surface__temperature',
        #'soil__temperature_reference_depth',
        #'soil__thermal_conductivity' ]
        
    #----------------------------------------------------------------
    # Should we use "ponded_water__depth" or "surface_water__depth"
    # instead of "channel_water__depth" in this case ?
    #----------------------------------------------------------------
    # Should we use "soil_surface__temperature" or
    # "land_surface__temperature" here ?  (Both, for now.)
    #----------------------------------------------------------------   
    _var_name_map = {   
        'atmosphere_bottom_air__temperature' :               'T_air',
        'land_surface__temperature':                         'T_surf',
        'land_surface_net-longwave-radiation__energy_flux':  'Qn_LW',
        'land_surface_net-shortwave-radiation__energy_flux': 'Qn_SW',
        #---------------------------------------------------------------
        'land_surface_soil__conduction_heat_flux' :           'Qc',   # (computed)
        'land_surface_water__domain_time_integral_of_evaporation_volume_flux': 'vol_ET',
        'land_surface_water__evaporation_volume_flux' :       'ET',
        'model__time_step':                                   'dt',
        #-----------------------------------------------------
        # These are read from GUI/file, but can be returned.
        #-----------------------------------------------------       
        'land_surface__elevation' :                       'DEM',
        'land_surface_water__priestley-taylor_alpha_coefficient': 'alpha',
        'soil__reference_depth_temperature' :             'T_soil_x',
        # 'soil_surface__temperature' :                   'T_surf',    # (from met)
        'soil__temperature_reference_depth':              'soil_x',
        'soil__thermal_conductivity' :                    'K_soil' }   # (thermal !)
        #----------------------------------------------
        # These are no longer needed here. (9/25/14)
        #---------------------------------------------- 
#         'channel_water_x-section__mean_depth' :           'depth', 
#         'soil_top-layer__porosity':                       'p0',
#         'soil_top-layer__saturated_thickness' :           'y0',
#         'soil_water_sat-zone_top_surface__elevation' :    'h_table' }    
        
    #------------------------------------------------
    # What is the correct unit string for "deg_C" ?
    #------------------------------------------------
    _var_units_map = {
        'atmosphere_bottom_air__temperature' :               'deg_C',
        'land_surface__temperature':                         'deg_C',      
        'land_surface_net-longwave-radiation__energy_flux':  'W m-2',
        'land_surface_net-shortwave-radiation__energy_flux': 'W m-2',
        #--------------------------------------------------------------
        'land_surface_soil__conduction_heat_flux' :      'W m-2',
        'land_surface_water__evaporation_volume_flux' :  'm s-1',
        'land_surface_water__domain_time_integral_of_evaporation_volume_flux': 'm3',
        'model__time_step' :                     's',
        #-----------------------------------------------------
        # These are read from GUI/file, but can be returned.
        #-----------------------------------------------------
        'land_surface__elevation' :                       'm',
        'land_surface_water__priestley-taylor_alpha_coefficient': '1',
        'soil__reference_depth_temperature' :             'deg_C',
        # 'soil_surface__temperature' :                   'deg_C',
        'soil__temperature_reference_depth':              'm',
        'soil__thermal_conductivity' :                    'W m-1 K-1]' } 
        #----------------------------------------------
        # These are no longer needed here. (9/25/14)
        #---------------------------------------------- 
#         'channel_water_x-section__mean_depth' :           'm',
#         'soil_top-layer__porosity':                       '1',
#         'soil_top-layer__saturated_thickness' :           'm',
#         'soil_water_sat-zone_top_surface__elevation' :    'm' }    
        
    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Evaporation_Priestley-Taylor'

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

        #--------------------------------------------------------
        # As of 7/9/10, Qn_SW and Qn_LW are computed internally
        # from other vars, including slope and aspect grids.
        # So they'll always be grids and so will self.ET
        # unless PRECIP_ONLY = True.
        #--------------------------------------------------------
        are_scalars = np.array([
                         self.is_scalar('T_soil_x'),
                         self.is_scalar('soil_x'),
                         self.is_scalar('K_soil'),
                         self.is_scalar('alpha'),
                         #-------------------------------
                         self.is_scalar('ET'),       # @evap
                         self.is_scalar('Qn_SW'),    # @met
                         self.is_scalar('Qn_LW'),    # @met
                         self.is_scalar('T_air'),    # @met
                         self.is_scalar('T_surf') ]) # @met
                         #---------------------------------
#                          self.is_scalar('depth'),     # d@chan
#                          self.is_scalar('h_table') ]) # satzone
                         #-------------------------------
##                         Qn_SW_IS_SCALAR,
##                         Qn_LW_IS_SCALAR,
##                         self.is_scalar('T_air'),
##                         self.is_scalar('T_surf') ])


        self.ALL_SCALARS = np.all(are_scalars)

        ## self.ALL_SCALARS = False
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def update_ET_rate(self):

        #--------------------------------------------------------------
        # Notes: Qet   = energy used for ET of water from surface
        #        Qn_SW = net shortwave irradiation flux (solar)
        #        Qn_LW = net longwave irradiation flux (air, surface)
        #        Qh    = sensible heat flux from turbulent convection
        #                between snow surface and air
        #        Qc    = energy transferred from surface to subsurface

        #        All of the Q's have units of [W/m^2].

        #        T_air    = air temperature [deg_C]
        #        T_surf   = soil temp at the surface [deg_C]
        #        T_soil_x = soil temp at depth of x meters [deg_C]

        #        Ks   = thermal conductivity of soil [W m-1 K-1]
        #        Ks = 0.45   ;[W m-1 K-1] (thawed soil; moisture content
        #                     near field capacity)
        #        Ks = 1.0    ;[W m-1 K-1] (frozen soil)

        #        alpha = evaporation parameter
        #        alpha = 0.95   ;(average found by Rouse)
        #        alpha = 1.26   ;(Jackson et al. (1996), at saturation)

        #        Modification of alpha:  alpha = (a1 * R) + a2
        #        R  = 1.0d   ;(equals 1 for saturation; R in [0,1])
        #        a1 = 1.0d   ;(accounts for moisture content of soil)
        #        a2 = 0.2d   ;(accounts for vegetation effect)
        #--------------------------------------------------------------
        Qn_SW  = self.Qn_SW   # (2/3/13, new framework)
        Qn_LW  = self.Qn_LW   # (2/3/13, new framework)
        T_air  = self.T_air   # (2/3/13, new framework)
        T_surf = self.T_surf  # (2/3/13, new framework)
        Q_net  = Qn_SW + Qn_LW
        
        #---------------------------------------------
        # Compute the conductive energy between the
        # surface and subsurface using Fourier's law
        #--------------------------------------------------------
        # This is now computed by update_Qc() in evap_base.py,
        # for both Priestley-Taylor and Energy Balance methods.
        # update() method calls update_Qc() & update_ET_rate().
        # Zhang et al. (2000) used -delta_T and (Q_net - Qc).
        # See detailed notes there.        
        #--------------------------------------------------------
        Qc = self.Qc
        # delta_T = (T_surf - self.T_soil_x)  
        # Qc  = self.K_soil * delta_T / self.soil_x
        # self.Qc[:] = Qc  # (in-place)
        #---------------------------------------------
        # In Qet formula, the constant 0.011 has
        # units of 1/[deg_C] to cancel T_air units.
        #---------------------------------------------
        Qet = self.alpha * (np.float64(0.406) + (np.float64(0.011) * T_air)) * (Q_net + Qc)
    
        #-----------------------------------
        # Convert ET energy to a loss rate
        #------------------------------------------
        # Lf = latent heat of fusion [J/kg]
        # Lv = latent heat of vaporization [J/kg]
        # ET = (Qet / (rho_w * Lv))
        #------------------------------------------
        # rho_w = 1000d       ;[kg/m^3]
        # Lv    = -2500000d   ;[J/kg]
        # So (rho_w * Lv) = -2.5e+9  [J/m^3]
        #-------------------------------------
        ET = (Qet / np.float64(2.5E+9))  #[m/s]  (A loss, but returned as positive.)

        #------------------------------- 
        # Save new ET values, in-place
        #-------------------------------
        np.maximum(ET, np.float64(0), self.ET)     # in-place     
        ## self.ET = np.maximum(ET, np.float64(0))

        TEST = False
        if (TEST):
            print('min(Qn_SW), max(Qn_SW) =', Qn_SW.min(), Qn_SW.max() )
            print('min(Qn_LW), max(Qn_LW) =', Qn_LW.min(), Qn_LW.max() )
            print('min(Qc),    max(Qc)    =', Qc.min(), Qc.max() )
            print('min(ET),    max(ET)    =', ET.min(), ET.max() )
            print('vol_ET =', self.vol_ET )
            print()
                       
        #----------------------------------------------------------
        # ET is 2D array in evap_base.initialize_computed_vars().
        # Should we ever allow ET to be a scalar ?
        # THIS WOULD BE COSTLY
        #----------------------------------------------------------
        # if (np.size(self.ET) == 1):
        #     self.ET += np.zeros((self.ny, self.nx), dtype='float64')
            
    #   update_ET_rate()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        #----------------------------------------------------
        # Note: Priestley-Taylor method needs alpha but the
        #       energy balance method doesn't. (2/5/13)
        #----------------------------------------------------
        self.alpha_file    = self.soil_directory + self.alpha_file
        self.K_soil_file   = self.soil_directory + self.K_soil_file
        self.soil_x_file   = self.soil_directory + self.soil_x_file
        self.T_soil_x_file = self.soil_directory + self.T_soil_x_file

        self.alpha_unit    = model_input.open_file(self.alpha_type,    self.alpha_file)
        self.K_soil_unit   = model_input.open_file(self.K_soil_type,   self.K_soil_file)
        self.soil_x_unit   = model_input.open_file(self.soil_x_type,   self.soil_x_file)
        self.T_soil_x_unit = model_input.open_file(self.T_soil_x_type, self.T_soil_x_file)
        
    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        rti = self.rti
        
        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        alpha = model_input.read_next(self.alpha_unit, self.alpha_type, rti)
        if (alpha is not None): self.alpha = alpha

        K_soil = model_input.read_next(self.K_soil_unit, self.K_soil_type, rti)
        if (K_soil is not None): self.K_soil = K_soil

        soil_x = model_input.read_next(self.soil_x_unit, self.soil_x_type, rti)
        if (soil_x is not None): self.soil_x = soil_x

        T_soil_x = model_input.read_next(self.T_soil_x_unit, self.T_soil_x_type, rti)
        if (T_soil_x is not None): self.T_soil_x = T_soil_x
        
    #   read_input_files()        
    #-------------------------------------------------------------------  
    def close_input_files(self):

        if (self.alpha_type    != 'Scalar'): self.alpha_unit.close()        
        if (self.K_soil_type   != 'Scalar'): self.K_soil_unit.close()
        if (self.soil_x_type   != 'Scalar'): self.soil_x_unit.close()
        if (self.T_soil_x_type != 'Scalar'): self.T_soil_x_unit.close()
        
    #   close_input_files()
    #-------------------------------------------------------------------


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------   
# def Priestley_Taylor_ET_Rate(alpha, Ks, T_soil_x, soil_x, \
#                              Qn_SW, Qn_LW, T_air, T_surf):
# 
#     #--------------------------------------------------------------
#     # Notes: Qet   = energy used for ET of water from surface
#     #        Qn_SW = net shortwave radiation flux (solar)
#     #        Qn_LW = net longwave radiation flux (air, surface)
#     #        Qh    = sensible heat flux from turbulent convection
#     #                between snow surface and air
#     #        Qc    = energy transferred from surface to subsurface
# 
#     #        All of the Q's have units of [W/m^2].
# 
#     #        T_air    = air temperature [deg_C]
#     #        T_surf   = soil temp at the surface [deg_C]
#     #        T_soil_x = soil temp at depth of x meters [deg_C]
# 
#     #        Ks   = thermal conductivity of soil [W m-1 K-1]
#     #        Ks = 0.45   ;[W m-1 K-1] (thawed soil; moisture content
#     #                     near field capacity)
#     #        Ks = 1.0    ;[W m-1 K-1] (frozen soil)
# 
#     #        alpha = evaporation parameter
#     #        alpha = 0.95   ;(average found by Rouse)
#     #        alpha = 1.26   ;(Jackson et al. (1996), at saturation)
# 
#     #        Modification of alpha:  alpha = (a1 * R) + a2
#     #        R  = 1.0d   ;(equals 1 for saturation; R in [0,1])
#     #        a1 = 1.0d   ;(accounts for moisture content of soil)
#     #        a2 = 0.2d   ;(accounts for vegetation effect)
#     #--------------------------------------------------------------
#     
#     #---------------------------------------------
#     # Compute the conductive energy between the
#     # surface and subsurface using Fourier's law
#     #---------------------------------------------
#     # soil_x is converted from [cm] to [m] when
#     # it is read from the GUI and then stored
#     #---------------------------------------------
#     # In Qet formula, the constant 0.011 has
#     # units of 1/[deg_C] to cancel T_air units.
#     #---------------------------------------------  
#     Qc   = Ks * (T_soil_x - T_surf) / (soil_x)
#     Qnet = Qn_SW + Qn_LW
#     Qet  = alpha * (np.float64(0.406) + (np.float32(0.011) * T_air)) * (Qnet - Qc)
#     
#     #-----------------------------------
#     # Convert ET energy to a loss rate
#     #------------------------------------------
#     # Lf = latent heat of fusion [J/kg]
#     # Lv = latent heat of vaporization [J/kg]
#     # ET = (Qet / (rho_w * Lv))
#     #------------------------------------------
#     # rho_w = 1000d       ;[kg/m^3]
#     # Lv    = -2500000d   ;[J/kg]
#     # So (rho_w * Lv) = -2.5e+9  [J/m^3]
#     #-------------------------------------
#     ET = (Qet / np.float32(2.5E+9))  #[m/s]  (A loss, but returned as positive.)
#     
#     return np.maximum(ET, np.float64(0))
# 
# #   Priestley_Taylor_ET_Rate
#-----------------------------------------------------------------------

