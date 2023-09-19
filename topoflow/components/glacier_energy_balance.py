"""
This file defines an "energy balance" glacier melt and snowmelt component and related
functions.  It inherits from the glacier "base class" in
"glacier_base.py".
"""
#-----------------------------------------------------------------------
#
#  Copyright (c) 2001-2023, Scott D. Peckham
#
#  Aug 2023.  Removed trailing space in 'J kg-1 K-1 '.
#             Added missing "vol_swe" in initialize_computed_vars.
#
#-----------------------------------------------------------------------
#
#  class snow_energy_balance
#
#      get_component_name()
#      get_attribute()          
#      get_input_var_names()   
#      get_output_var_names()   
#      get_var_name()          
#      get_var_units()         
#      ----------------------
#      check_input_types()
#      initialize_input_file_vars()   
#      initialize_snow_cold_content() 
#      initialize_ice_cold_content()   
#      initialize_computed_vars()  # from glacier_base.py
#      update_snow_meltrate()
#      update_ice_meltrate()
#      ----------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#
#-----------------------------------------------------------------------

import numpy as np

from topoflow.utils import model_input
from topoflow.components import glacier_base

#-----------------------------------------------------------------------
class glacier_component( glacier_base.glacier_component ):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Glacier_Energy_Balance',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham, Lauren A. Bolotin',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'GlacierEnergyBalance',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Glacier_Energy_Balance.cfg.in',
        'cfg_extension':      '_glacier_energy_balance.cfg',
        'cmt_var_prefix':     '/GlacierEnergyBalance/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Snow_Energy_Balance.xml', # LB: ??
        'dialog_title':       'Snow and Glacier Melt: Energy Balance Parameters', 
        'time_units':         'seconds' }

    #----------------------------------------------------------
    # The Energy-Balance method needs several vars from the
    # Met component, such as T_air, Q_sum.  The others
    # are used for consistent use of scalars vs. grids;
    # see "check_input_types()" below. (e.g. p0 and RH)
    #----------------------------------------------------------
    # Some input vars, like ****** are read from the CFG
    # file and not from other components.  They are therefore
    # included with the output_vars.
    #------------------------------------------------------------
    # Note: snowfall_leq-volume_flux *must* be nonliquid precip.
    #------------------------------------------------------------  
    _input_var_names = [
        'atmosphere_bottom_air__mass-per-volume_density',
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity',
        'atmosphere_bottom_air__temperature',
        'atmosphere_water__snowfall_leq-volume_flux',
        'land_surface__temperature',    # (used to initialize Eccs, Ecci)
        'land_surface_net-total-energy__energy_flux',
        'water-liquid__mass-per-volume_density' ]

        #------------------------------------------------------------
        # These are used by Meteorology component to compute Q_sum.
        # but are not needed directly by this component. 
        #------------------------------------------------------------
        # 'atmosphere_bottom_air__pressure',
        # 'atmosphere_bottom_air_flow__log_law_roughness_length',
        # 'atmosphere_bottom_air_flow__speed_reference_height',
        # 'atmosphere_bottom_air_flow__reference-height_speed',
        # 'atmosphere_bottom_air_water-vapor__relative_saturation',
        # 'land_surface_net-longwave-radiation__energy_flux',
        # 'land_surface_net-shortwave-radiation__energy_flux',
                        
                                  
    #-------------------------------------------------------------
    # Note: The main output of the Energy-Balance method is SM.
    #       Functions in glacier_base.py compute additional output
    #       vars, such as:
    #       update_SM_integral(), update_snow/ice_depth(), 
    #       update_swe/iwe().
    #       They need things like: rho_H2O, rho_snow, and rho_ice.
    #-------------------------------------------------------------
    # Note: Cp_snow and Cp_ice are constants set in the 
    # "set_constants()" function in glacier_base.py.
    #------------------------------------------------------------

    _output_var_names = [
        'model__time_step',                            # dt   
        'snowpack__domain_time_integral_of_melt_volume_flux',   # vol_SM
        'snowpack__initial_domain_integral_of_liquid-equivalent_depth', # vol_swe_start'
        'snowpack__domain_integral_of_liquid-equivalent_depth',         # vol_swe
        'snowpack__energy-per-area_cold_content',      # Eccs
        'snowpack__depth',                             # h_snow
        'snowpack__initial_depth',                     # h0_snow
        'snowpack__initial_liquid-equivalent_depth',   # h0_swe
        'snowpack__liquid-equivalent_depth',           # h_swe
        'snowpack__melt_volume_flux',                  # SM
        'snowpack__z_mean_of_mass-per-volume_density', # rho_snow 
        'snowpack__z_mean_of_mass-specific_isobaric_heat_capacity', # Cp_snow   
        'glacier_ice__domain_time_integral_of_melt_volume_flux', # vol_IM
        'glacier__initial_domain_integral_of_liquid-equivalent_depth', # vol_iwe_start
        'glacier__domain_integral_of_liquid-equivalent_depth',         # vol_iwe
        'glacier__energy-per-area_cold_content',      # Ecci
        'glacier_ice__thickness', # h_ice
        'glacier_ice__initial_thickness', #h0_ice
        'glacier__initial_liquid_equivalent_depth', # h0_iwe
        'glacier_ice__melt_volume_flux', # IM
        'glacier_ice__mass-per-volume_density', # rho_ice
        'glacier_ice__mass-specific_isobaric_heat_capacity' ] # Cp_ice 
    
    _var_name_map = {
        'atmosphere_bottom_air__mass-per-volume_density': 'rho_air',
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'Cp_air',
        'atmosphere_bottom_air__temperature': 'T_air',
        'atmosphere_water__snowfall_leq-volume_flux': 'P_snow',
        'land_surface_net-total-energy__energy_flux': 'Q_sum',
        'land_surface__temperature': 'T_surf',        
        'water-liquid__mass-per-volume_density': 'rho_H2O',
        #------------------------------------------------------------
        # These are used by Meteorology component to compute Q_sum.
        # but are not needed directly by this component. 
        #------------------------------------------------------------
        #'atmosphere_bottom_air__pressure': 'p0',
        #'atmosphere_bottom_air_flow__log_law_roughness_length': 'z0_air',
        #'atmosphere_bottom_air_flow__speed_reference_height': 'z',
        #'atmosphere_bottom_air_flow__reference-height_speed': 'uz',
        #'atmosphere_bottom_air_water-vapor__relative_saturation': 'RH',        
        #'land_surface_net-longwave-radiation__energy_flux': 'Qn_LW',
        #'land_surface_net-shortwave-radiation__energy_flux': 'Qn_SW',        
        #----------------------------------------------------------
        'model__time_step': 'dt',     
        'snowpack__domain_time_integral_of_melt_volume_flux':   'vol_SM',
        'snowpack__initial_domain_integral_of_liquid-equivalent_depth': 'vol_swe_start',
        'snowpack__domain_integral_of_liquid-equivalent_depth':         'vol_swe',           
        'snowpack__depth': 'h_snow',
        'snowpack__energy-per-area_cold_content': 'Eccs',
        'snowpack__initial_depth': 'h0_snow',
        'snowpack__initial_liquid-equivalent_depth': 'h0_swe',
        'snowpack__liquid-equivalent_depth': 'h_swe',
        'snowpack__melt_volume_flux': 'SM',
        'snowpack__z_mean_of_mass-per-volume_density': 'rho_snow',
        'snowpack__z_mean_of_mass-specific_isobaric_heat_capacity': 'Cp_snow',
        'glacier_ice__domain_time_integral_of_melt_volume_flux': 'vol_IM', 
        'glacier__initial_domain_integral_of_liquid-equivalent_depth': 'vol_iwe_start',
        'glacier__domain_integral_of_liquid-equivalent_depth': 'vol_iwe',
        'glacier__energy-per-area_cold_content': 'Ecci',      
        'glacier_ice__thickness': 'h_ice',
        'glacier_ice__initial_thickness': 'h0_ice',
        'glacier__initial_liquid_equivalent_depth': 'h0_iwe',
        'glacier_ice__melt_volume_flux': 'IM',
        'glacier_ice__mass-per-volume_density': 'rho_ice',
        'glacier_ice__mass-specific_isobaric_heat_capacity': 'Cp_ice' }
    
    #-----------------------------------------------------------------
    # Note: We need to be careful with whether units are C or K,
    #       for all "thermal" quantities (e.g. thermal_capacity).
    #-----------------------------------------------------------------       
    _var_units_map = {
        'atmosphere_bottom_air__mass-per-volume_density': 'kg m-3',
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'J kg-1 K-1', # (see Notes above)
        'atmosphere_bottom_air__temperature': 'deg_C',  # (see Notes above)
        'atmosphere_water__snowfall_leq-volume_flux': 'm s-1',
        'land_surface_net-total-energy__energy_flux': 'W m-2',
        'land_surface__temperature': 'deg_C',        
        'water-liquid__mass-per-volume_density': 'kg m-3',
        #------------------------------------------------------------
        # These are used by Meteorology component to compute Q_sum.
        # but are not needed directly by this component. 
        #------------------------------------------------------------
        # 'atmosphere_bottom_air__pressure': 'mbar',
        # 'atmosphere_bottom_air_flow__log_law_roughness_length': 'm',
        # 'atmosphere_bottom_air_flow__speed_reference_height': 'm',
        # 'atmosphere_bottom_air_flow__reference-height_speed': 'm s-1',
        # 'atmosphere_bottom_air_water-vapor__relative_saturation': '1',                
        # 'land_surface_net-longwave-radiation__energy_flux': 'W m-2',
        # 'land_surface_net-shortwave-radiation__energy_flux': 'W m-2',        
        #--------------------------------------------------------------
        'model__time_step': 's',
        'snowpack__domain_time_integral_of_melt_volume_flux': 'm3',
        'snowpack__initial_domain_integral_of_liquid-equivalent_depth': 'm3',  
        'snowpack__domain_integral_of_liquid-equivalent_depth': 'm3',          
        'snowpack__depth': 'm',
        'snowpack__energy-per-area_cold_content': 'J m-2',
        'snowpack__initial_depth': 'm',
        'snowpack__initial_liquid-equivalent_depth': 'm',
        'snowpack__liquid-equivalent_depth': 'm',
        'snowpack__melt_volume_flux': 'm s-1',
        'snowpack__z_mean_of_mass-per-volume_density': 'kg m-3',
        'snowpack__z_mean_of_mass-specific_isobaric_heat_capacity': 'J kg-1 K-1',
        'glacier_ice__domain_time_integral_of_melt_volume_flux': 'm3', 
        'glacier__initial_domain_integral_of_liquid-equivalent_depth': 'm3',
        'glacier__domain_integral_of_liquid-equivalent_depth': 'm3',
        'glacier__energy-per-area_cold_content': 'J m-2',
        'glacier_ice__thickness': 'm',
        'glacier_ice__initial_thickness': 'm',
        'glacier__initial_liquid_equivalent_depth': 'm',
        'glacier_ice__melt_volume_flux': 'm s-1',
        'glacier_ice__mass-per-volume_density': 'kg m-3',
        'glacier_ice__mass-specific_isobaric_heat_capacity': 'J kg-1 K-1' 
          }

    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Glacier_Energy_Balance'

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
        #       components vs. those read from files.
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
    def check_input_types(self):

        #--------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_air, Cp_air,
        # rho_ice, and Cp_ice are currently always scalars.
        #--------------------------------------------------        
        are_scalars = np.array([
                          self.is_scalar('P_snow'),
                          self.is_scalar('rho_H2O'),
                          self.is_scalar('rho_air'),
                          self.is_scalar('Cp_air'),
                          self.is_scalar('T_air'),   ##### CHECK THIS ONE.
                          self.is_scalar('T_surf'),
                          self.is_scalar('Q_sum'),
                          self.is_scalar('rho_snow'),
                          self.is_scalar('Cp_snow'),
                          self.is_scalar('h0_snow'),
                          self.is_scalar('h0_swe'),
                          self.is_scalar('rho_ice'),
                          self.is_scalar('Cp_ice'),
                          self.is_scalar('h0_ice'),
                          self.is_scalar('h0_iwe'),
                          self.is_scalar('h_active_layer')  ])


        self.ALL_SCALARS = np.all(are_scalars)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def initialize_input_file_vars(self):     
    
        #---------------------------------------------------
        # Initialize vars to be read from files 
        #---------------------------------------------------
        # Need this in order to use "bmi.update_var()".
        #----------------------------------------------------------
        # NOTE: read_config_file() sets these to '0.0' if they
        #       are not type "Scalar", so self has the attribute.
        #----------------------------------------------------------
        dtype = 'float64'
        if (self.Cp_snow_type.lower() != 'scalar'):
            self.Cp_snow = self.initialize_var(self.Cp_snow_type, dtype=dtype)
        if (self.rho_snow_type.lower() != 'scalar'):
            self.rho_snow = self.initialize_var(self.rho_snow_type, dtype=dtype)
        if (self.T0_type.lower() != 'scalar'):
            self.T0 = self.initialize_var(self.T0_type, dtype=dtype)
        if (self.h0_snow_type.lower() != 'scalar'):
            self.h0_snow = self.initialize_var(self.h0_snow_type, dtype=dtype)    
        if (self.h0_swe_type.lower() != 'scalar'):
            self.h0_swe = self.initialize_var(self.h0_swe_type, dtype=dtype)    
        if (self.Cp_ice_type.lower() != 'scalar'):
            self.Cp_ice = self.initialize_var(self.Cp_ice_type, dtype=dtype)
        if (self.rho_ice_type.lower() != 'scalar'):
            self.rho_ice = self.initialize_var(self.rho_ice_type, dtype=dtype)
        if (self.h0_ice_type.lower() != 'scalar'):
            self.h0_ice = self.initialize_var(self.h0_ice_type, dtype=dtype)
        if (self.h0_iwe_type.lower() != 'scalar'):
            self.h0_iwe = self.initialize_var(self.h0_iwe_type, dtype=dtype) 
        if (self.h_active_layer_type.lower() != 'scalar'):
            self.h_active_layer = self.initialize_var(self.h_active_layer_type, dtype=dtype)
    
    #   initialize_input_file_vars()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        #----------------------------------------------
        # NOTE:  This function overrides the version
        #        that is inherited from glacier_base.py.
        #----------------------------------------------
        # If T_air or precip are grids, then make
        # sure that h_snow, h_swe, h_ice, and h_iwe 
        # are grids
        #----------------------------------------------
        T_IS_GRID = self.is_grid('T_air')
        P_IS_GRID = self.is_grid('P_snow')
        H0_SNOW_IS_SCALAR = self.is_scalar('h0_snow')
        H0_SWE_IS_SCALAR  = self.is_scalar('h0_swe') 
        H0_ICE_IS_SCALAR  = self.is_scalar('h0_ice')
        H0_IWE_IS_SCALAR  = self.is_scalar('h0_iwe') 

        #------------------------------------------------------
        # If h0_snow, h0_swe, h0_ice, or h0_iwe are scalars,
        # the use of copy() here requires they were converted 
        # to numpy scalars. Using copy() may not be necessary 
        # for scalars.
        #------------------------------------------------------
        h_snow = self.h0_snow.copy()    # [meters]
        h_swe  = self.h0_swe.copy()     # [meters]
        h_ice  = self.h0_ice.copy()     # [meters]
        h_iwe  = self.h0_iwe.copy()     # [meters]

        #--------------------------------------------------------------       
        # For the Energy Balance method, SM, IM, h_snow, h_swe, 
        # and h_ice are always grids because Q_sum is always a grid.
        #---------------------------------------------------------------------
        # Convert h_snow, h_swe, h_ice and h_iwe to grids if not already grids
        #---------------------------------------------------------------------
        if (H0_SNOW_IS_SCALAR):
            self.h_snow = h_snow + np.zeros([self.ny, self.nx], dtype='float64')
        else:
            self.h_snow = h_snow  # (is already a grid)
        if (H0_SWE_IS_SCALAR):
            self.h_swe = h_swe + np.zeros([self.ny, self.nx], dtype='float64')
        else:
            self.h_swe = h_swe    # (is already a grid) 
        if (H0_ICE_IS_SCALAR):
            self.h_ice = h_ice + np.zeros([self.ny, self.nx], dtype='float64')
        else:
            self.h_ice = h_ice # (is already a grid)         
        if (H0_IWE_IS_SCALAR):
            self.h_iwe = h_iwe + np.zeros([self.ny, self.nx], dtype='float64')
        else:
            self.h_iwe = h_iwe    # (is already a grid)

        self.SM      = np.zeros([self.ny, self.nx], dtype='float64')
        self.IM      = np.zeros([self.ny, self.nx], dtype='float64')
        self.vol_SM  = self.initialize_scalar( 0, dtype='float64') # (m3)
        self.vol_IM  = self.initialize_scalar( 0, dtype='float64') # (m3)
        #-------------------------------------------
        # 2023-08-28.  Added next line to fix bug.
        #-------------------------------------------
        self.vol_swe = self.initialize_scalar( 0, dtype='float64') # (m3)         
        self.vol_swe_start = self.initialize_scalar( 0, dtype='float64') # (m3)  
        self.vol_iwe = self.initialize_scalar( 0, dtype='float64') # (m3)         
        self.vol_iwe_start = self.initialize_scalar( 0, dtype='float64') # (m3)  
        #-----------------------------------------------------
        # Compute density ratio for water to snow/ice
        # rho_H2O is for liquid water close to 0 degrees C.
        # Water is denser than snow/ice, so density_ratio > 1.
        #-----------------------------------------------------
        self.ws_density_ratio = (self.rho_H2O / self.rho_snow)
        self.wi_density_ratio = (self.rho_H2O / self.rho_ice)

        #-------------------------------------------------------------
        # Initialize the cold content of snowpack/ice column
        #-------------------------------------------------------------
        # This is the only difference from initialize_computed_vars()
        # method in glacier_base.py.
        #-------------------------------------------------------------
        self.initialize_snow_cold_content()
        self.initialize_ice_cold_content()
        
    #   initialize_computed_vars()
    #---------------------------------------------------------------------
    def initialize_snow_cold_content( self ):

        #----------------------------------------------------------------
        # NOTES: This function is used to initialize the cold content
        #        of a snowpack.
        #        The cold content has units of [J m-2] (_NOT_ [W m-2]).
        #        It is an energy (per unit area) threshold (or deficit)
        #        that must be overcome before melting of snow can occur.
        #        Cold content changes over time as the snowpack warms or
        #        cools, but must always be non-negative.
        #
        #        K_snow is between 0.063 and 0.71  [W m-1 K-1]
        #        All of the Q's have units of W m-2 = J s-1 m-2).
        #
        #        T0 is read from the config file.  This is a different
        #        T0 from the one used by the Degree-Day method.
        #---------------------------------------------------------------

        #--------------------------------------------
        # Compute initial cold content of snowpack
        #--------------------------------------------
        T_snow    = self.T_surf
        del_T     = (self.T0 - T_snow)
        self.Eccs  = (self.rho_snow * self.Cp_snow) * self.h0_snow * del_T
        # print('Snow Cold Content:')
        # print(self.Eccs)
        # print("T_Surf:")
        # print(self.T_surf)
        # print('Cp_snow')
        # print(self.Cp_snow)
        # print('rho_snow')
        # print(self.rho_snow)
        # print('h0_snow')
        # print(self.h0_snow)

        #------------------------------------        
        # Cold content must be nonnegative.
        #----------------------------------------------
        # Eccs > 0 if (T_snow < T0).  i.e. T_snow < 0.
        #----------------------------------------------
        self.Eccs = np.maximum( self.Eccs, np.float64(0))
        ### np.maximum( self.Ecc, np.float64(0), self.Ecc)  # (in place)
        # print('Updated Snow Cold Content:')
        # print(self.Eccs)
        # print('Max initial Eccs:')
        # print(np.max(self.Eccs))

    #   initialize_snow_cold_content()
    #-------------------------------------------------------------------
    def initialize_ice_cold_content( self ):

        #----------------------------------------------------------------
        # NOTES: This function is used to initialize the cold content
        #        of glacier ice.
        #        The cold content has units of [J m-2] (_NOT_ [W m-2]).
        #        It is an energy (per unit area) threshold (or deficit)
        #        that must be overcome before melting of ice can occur.
        #        Cold content changes over time as the ice warms or
        #        cools, but must always be non-negative.
        #
        #        K_snow is between 0.063 and 0.71  [W m-1 K-1]
        #        All of the Q's have units of W m-2 = J s-1 m-2).
        #
        #        T0 is read from the config file.  This is a different
        #        T0 from the one used by the Degree-Day method.
        #---------------------------------------------------------------

        #--------------------------------------------
        # Compute initial cold content of ice
        #--------------------------------------------
        T_ice    = self.T_surf
        del_T     = (self.T0 - T_ice)
        self.Ecci  = (self.rho_ice * self.Cp_ice) * self.h_active_layer * del_T 
        # print('Ice Cold Content:')
        # print(self.Ecci)
        # print('Cp_ice')
        # print(self.Cp_ice)
        # print('rho_ice')
        # print(self.rho_ice)
        # print('h0_ice')
        # print(self.h0_ice)

        #------------------------------------        
        # Cold content must be nonnegative.
        #----------------------------------------------
        # Ecci > 0 if (T_ice < T0).  i.e. T_ice < 0.
        #----------------------------------------------
        self.Ecci = np.maximum( self.Ecci, np.float64(0))
        ### np.maximum( self.Ecci, np.float64(0), self.Ecci)  # (in place)
        # print('Updated Ice Cold Content:')
        # print(self.Ecci)
        # print('Max initial Ecci:')
        # print(np.max(self.Ecci))

    #   initialize_ice_cold_content()
    #-------------------------------------------------------------------
    def update_snow_meltrate(self):

        #------------------------------------------------------------
        # Notes: This computes a "potential" meltrate, which can't
        #        be realized unless there is enough snow.
        #        See snow_base.enforce_max_meltrate().
        #------------------------------------------------------------        
        # Notes: See notes in "met_base.py" for the method called
        #        "update_net_energy_flux()".

        #        This version uses "Q_sum" which is computed as a
        #        state variable for a meteorology component
        #        (e.g. met_base.py).

        #        Arguments are assumed to be references to a scalar
        #            or grid so that:

        #        M  = water equivalent of snowmelt [m/s]
        #        M_max  = max possible meltrate if all snow melts
        #        T_air  = air temperature [deg_C]

        #        Model must start when snow is isothermal. (CHECK)
        #        Cooling of the snowpack is not considered.

        #        86400000d = 1000 [mm/m] * 60 [sec/min] *
        #                    60 [min/sec] * 24 [hrs/day]

        #        rho_snow is not needed in formula here, but is
        #        needed to convert snowmelt to water equivalent?
        #-------------------------------------------------------------
     
        #----------------------------------
        # Compute energy-balance meltrate   
        #------------------------------------------------------
        # Eccs is initialized by initialize_snow_cold_content().
        #------------------------------------------------------
        # The following pseudocode only works for scalars but
        # is otherwise equivalent to that given below and
        # clarifies the logic:
        #------------------------------------------------------
        #  if (Q_sum gt 0) then begin
        #      if ((Q_sum * dt) gt Eccs) then begin
        #          ;-------------------------------------------
        #          ; Snow is melting.  Use some of Q_sum to
        #          ; overcome Eccs, and remainder to melt snow
        #          ;-------------------------------------------
        #          Qm  = Q_sum - (Eccs/dt)
        #          Eccs = 0
        #          M   = (Qm / (rho_w * Lf))
        #      endif else begin
        #          ;------------------------------
        #          ; Snow is warming; reduce Eccs
        #          ;------------------------------
        #          Eccs = (Eccs - (Q_sum * dt))
        #          M   = 0d
        #      endelse
        #  endif else begin
        #      ;--------------------------------
        #      ; Snow is cooling; increase Eccs
        #      ;--------------------------------
        #      Eccs = Eccs - (Q_sum * dt)
        #      M   = 0d
        #  endelse
        #---------------------------------------------------------
        # Q_sum = Qn_SW + Qn_LW + Qh + Qe + Qa + Qc    # [W m-2]
        #---------------------------------------------------------
        
        #-----------------------------------------------        
        # New approach; easier to understand
        #----------------------------------------------- 
        # E_in  = energy input over one time step
        # E_rem = energy remaining in excess of Eccs
        #-----------------------------------------------
        E_in  = (self.Q_sum * self.dt)
        E_rem = np.maximum( E_in - self.Eccs, np.float64(0) )
        Qm    = (E_rem / self.dt)  # [W m-2]
        print(E_in)
        
        #-------------------------------------
        # Convert melt energy to a melt rate
        #------------------------------------------
        # Lf = latent heat of fusion [J/kg]
        # Lv = latent heat of vaporization [J/kg]
        # M  = (Qm / (rho_w * Lf))
        #------------------------------------------
        # rho_w = 1000d       ;[kg/m^3]
        # Lf    = 334000d     ;[J/kg = W*s/kg]
        #------------------------------------------
        M       = (Qm / (self.rho_H2O * self.Lf))   #[m/s]
        self.SM = np.maximum(M, np.float64(0))

        #--------------------------------------------------
        # Update the cold content of the snowpack [J m-2]
        # If this is positive, there was no melt so far.
        #--------------------------------------------------
        self.Eccs = np.maximum((self.Eccs - E_in), np.float64(0))
        # print('Max Eccs:')
        # print(np.max(self.Eccs))

        #-----------------------------------------------------------
        # Note: enforce_max_snow_meltrate() method is always called
        #       by the base class to make sure that meltrate
        #       does not exceed the max possible.
        #-----------------------------------------------------------
        #  self.enforce_max_meltrate()

    #   update_snow_meltrate()
    #---------------------------------------------------------------------
    def update_ice_meltrate(self):

        #------------------------------------------------------------
        # Notes: This computes a "potential" meltrate, which can't
        #        be realized unless there is enough ice.
        #        See glacier_base.enforce_max_ice_meltrate().
        #------------------------------------------------------------        
        # Notes: See notes in "met_base.py" for the method called
        #        "update_net_energy_flux()".

        #        This version uses "Q_sum" which is computed as a
        #        state variable for a meteorology component
        #        (e.g. met_base.py).

        #        Arguments are assumed to be references to a scalar
        #            or grid so that:

        #        M  = water equivalent of ice melt [m/s]
        #        M_max  = max possible meltrate if all ice melts
        #        T_air  = air temperature [deg_C]

        #        Model must start when ice is isothermal. (CHECK) # LB: full ice column?
        #        Cooling of ice is not considered.

        #        86400000d = 1000 [mm/m] * 60 [sec/min] *
        #                    60 [min/sec] * 24 [hrs/day]

        #        rho_ice is not needed in formula here, but is
        #        needed to convert ice melt to water equivalent?
        #-------------------------------------------------------------
     
        #----------------------------------
        # Compute energy-balance meltrate   
        #------------------------------------------------------
        # Ecci is initialized by initialize_ice_cold_content().
        #------------------------------------------------------
        # The following pseudocode only works for scalars but
        # is otherwise equivalent to that given below and
        # clarifies the logic:
        #------------------------------------------------------
        #  if (Q_sum gt 0) then begin
        #      if ((Q_sum * dt) gt Ecci) then begin
        #          ;-------------------------------------------
        #          ; Ice is melting.  Use some of Q_sum to
        #          ; overcome Ecci, and remainder to melt ice
        #          ;-------------------------------------------
        #          Qm  = Q_sum - (Ecci/dt)
        #          Ecci = 0
        #          M   = (Qm / (rho_w * Lf))
        #      endif else begin
        #          ;------------------------------
        #          ; Ice is warming; reduce Ecci
        #          ;------------------------------
        #          Ecci = (Ecci - (Q_sum * dt))
        #          M   = 0d
        #      endelse
        #  endif else begin
        #      ;--------------------------------
        #      ; Ice is cooling; increase Ecci
        #      ;--------------------------------
        #      Ecci = Ecci - (Q_sum * dt)
        #      M   = 0d
        #  endelse
        #---------------------------------------------------------
        # Q_sum = Qn_SW + Qn_LW + Qh + Qe + Qa + Qc    # [W m-2]
        #---------------------------------------------------------
        
        #-----------------------------------------------        
        # New approach; easier to understand 
        #----------------------------------------------- 
        # E_in  = energy input over one time step
        # E_rem = energy remaining in excess of Ecci
        #-----------------------------------------------
        E_in  = (self.Q_sum * self.dt)
        E_rem = np.maximum( E_in - self.Ecci, np.float64(0) )
        Qm    = (E_rem / self.dt)  # [W m-2]
        
        #-------------------------------------
        # Convert melt energy to a melt rate
        #------------------------------------------
        # Lf = latent heat of fusion [J/kg]
        # Lv = latent heat of vaporization [J/kg]
        # M  = (Qm / (rho_w * Lf))
        #------------------------------------------
        # rho_w = 1000d       ;[kg/m^3]
        # Lf    = 334000d     ;[J/kg = W*s/kg]
        #------------------------------------------
        M       = (Qm / (self.rho_H2O * self.Lf))   #[m/s]
        self.IM = np.maximum(M, np.float64(0))

        #--------------------------------------------------
        # Update the cold content of the ice [J m-2]
        # If this is positive, there was no melt so far.
        #--------------------------------------------------
   
        self.Ecci = np.maximum((self.Ecci - E_in), np.float64(0))
        # print('Max Ecci:')
        # print(np.max(self.Ecci))

        #----------------------------------------------------------
        # Note: enforce_max_ice_meltrate() method is always called
        #       by the base class to make sure that meltrate
        #       does not exceed the max possible.
        #----------------------------------------------------------
        #  self.enforce_max_ice_meltrate()
            
    #   update_ice_meltrate()
    #---------------------------------------------------------------------
    def open_input_files(self):

        self.Cp_snow_file  = self.in_directory + self.Cp_snow_file
        self.rho_snow_file = self.in_directory + self.rho_snow_file
        self.T0_file       = self.in_directory + self.T0_file
        self.h0_snow_file  = self.in_directory + self.h0_snow_file
        self.h0_swe_file   = self.in_directory + self.h0_swe_file
        self.Cp_ice_file   = self.in_directory + self.Cp_ice_file
        self.rho_ice_file  = self.in_directory + self.rho_ice_file
        self.h0_ice_file   = self.in_directory + self.h0_ice_file
        self.h0_iwe_file   = self.in_directory + self.h0_iwe_file
        self.h_active_layer_file = self.in_directory + self.h_active_layer_file

        self.Cp_snow_unit  = model_input.open_file(self.Cp_snow_type,  self.Cp_snow_file)
        self.rho_snow_unit = model_input.open_file(self.rho_snow_type, self.rho_snow_file)
        self.T0_unit       = model_input.open_file(self.T0_type,       self.T0_file)
        self.h0_snow_unit  = model_input.open_file(self.h0_snow_type,  self.h0_snow_file)
        self.h0_swe_unit   = model_input.open_file(self.h0_swe_type,   self.h0_swe_file)
        self.Cp_ice_unit   = model_input.open_file(self.Cp_ice_type,   self.Cp_ice_file)
        self.rho_ice_unit  = model_input.open_file(self.rho_ice_type,  self.rho_ice_file)
        self.h0_ice_unit   = model_input.open_file(self.h0_ice_type,   self.h0_ice_file)
        self.h0_iwe_unit   = model_input.open_file(self.h0_iwe_type,   self.h0_iwe_file)
        self.h_active_layer_unit = model_input.open_file(self.h_active_layer_type, self.h_active_layer_file)

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        rti = self.rti

        #--------------------------------------------------------
        # All grids are assumed to have a data type of float32.
        #--------------------------------------------------------
        Cp_snow = model_input.read_next(self.Cp_snow_unit, self.Cp_snow_type, rti)
        if (Cp_snow is not None):
            self.update_var( 'Cp_snow', Cp_snow )
        
        rho_snow = model_input.read_next(self.rho_snow_unit, self.rho_snow_type, rti)
        if (rho_snow is not None):
            self.update_var( 'rho_snow', rho_snow )

        T0 = model_input.read_next(self.T0_unit, self.T0_type, rti)
        if (T0 is not None):
            self.update_var( 'T0', T0 )
        
        h0_snow = model_input.read_next(self.h0_snow_unit, self.h0_snow_type, rti)
        if (h0_snow is not None):
            self.update_var( 'h0_snow', h0_snow )
        
        h0_swe = model_input.read_next(self.h0_swe_unit, self.h0_swe_type, rti)
        if (h0_swe is not None):
            self.update_var( 'h0_swe', h0_swe )

        Cp_ice = model_input.read_next(self.Cp_ice_unit, self.Cp_ice_type, rti)
        if (Cp_ice is not None):
            self.update_var( 'Cp_ice', Cp_ice)

        rho_ice = model_input.read_next(self.rho_ice_unit, self.rho_ice_type, rti)
        if (rho_ice is not None):
            self.update_var( 'rho_ice', rho_ice)
        
        h0_ice = model_input.read_next(self.h0_ice, self.h0_ice_type, rti)
        if (h0_ice is not None):
            self.update_var( 'h0_ice', h0_ice)

        h0_iwe = model_input.read_next(self.h0_iwe_unit, self.h0_iwe_type, rti)
        if (h0_iwe is not None):
            self.update_var( 'h0_iwe', h0_iwe )

        h_active_layer = model_input.read_next(self.h_active_layer_unit, self.h_active_layer_type, rti)
        if (h_active_layer is not None):
            self.update_var( 'h_active_layer', h_active_layer )
        
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):
     
        if (self.Cp_snow_type  != 'Scalar'): self.Cp_snow_unit.close()
        if (self.rho_snow_type != 'Scalar'): self.rho_snow_unit.close()
        if (self.T0_type       != 'Scalar'): self.T0_unit.close()
        if (self.h0_snow_type  != 'Scalar'): self.h0_snow_unit.close()
        if (self.h0_swe_type   != 'Scalar'): self.h0_swe_unit.close()
        if (self.Cp_ice_type   != 'Scalar'): self.Cp_ice_unit.close()
        if (self.rho_ice_type  != 'Scalar'): self.rho_ice_unit.close()
        if (self.h0_ice_type   != 'Scalar'): self.h0_ice_unit.close()
        if (self.h0_iwe_type   != 'Scalar'): self.h0_iwe_unit.close()
        if (self.h_active_layer_type != 'Scalar'): self.h_active_layer_unit.close()

    #   close_input_files()    
    #-------------------------------------------------------------------    
#-------------------------------------------------------------------------  