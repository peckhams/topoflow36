"""
This file defines an "energy balance" glacier melt and snowmelt component and related
functions.  It inherits from the glacier "base class" in
"glacier_base.py".
"""
#-----------------------------------------------------------------------
#
#  Copyright (c) 2001-2023, Scott D. Peckham
#
#  Oct 2023. Updated the original snow cold content routine
#            so that new cold content is added with new snow in 
#            update_snowfall_cold_content(), and still accounts
#            for land surface energy fluxes in update_snowpack_cold_content(), 
#            and there is no cold content if there is no snow.
#  Sep 2023. Checked sign in initialize_cold_content().
#            Separate function: update_cold_content().
#            Moved initialize_cold_content() back to base class.
#            Renamed T0 to avoid conflict w/ degree-day method.
#            Removed separate version of initialize_computed_vars().
#  Aug 2023. Renamed snow functions to have 'snow' in the name
#            Duplicated existing snow functions for 'ice'
#
#-----------------------------------------------------------------------
#
#  class glacier_energy_balance
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
#      initialize_computed_vars()  # from glacier_base.py
#      update_snow_meltrate()
#      update_snowfall_cold_content()
#      update_snowpack_cold_content()
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
        'atmosphere_bottom_air_water-vapor__relative_saturation',
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
    # Note: The main output of the Energy-Balance method is SM/IM.
    #       Functions in glacier_base.py compute additional output
    #       vars, such as:
    #       update_SM/IM_integral(), update_snow/ice_depth(), 
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
        'glacier_ice__mass-specific_isobaric_heat_capacity', # Cp_ice 
        'cryosphere__melt_volume_flux', # M_total
        'cryosphere__domain_time_integral_of_melt_volume_flux' ] #vol_M_total 
    
    _var_name_map = {
        'atmosphere_bottom_air__mass-per-volume_density': 'rho_air',
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'Cp_air',
        'atmosphere_bottom_air__temperature': 'T_air',
        'atmosphere_bottom_air_water-vapor__relative_saturation': 'RH',
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
        'glacier_ice__mass-specific_isobaric_heat_capacity': 'Cp_ice',
        'cryosphere__melt_volume_flux': 'M_total',
        'cryosphere__domain_time_integral_of_melt_volume_flux' : 'vol_M_total'}
    
    #-----------------------------------------------------------------
    # Note: We need to be careful with whether units are C or K,
    #       for all "thermal" quantities (e.g. thermal_capacity).
    #-----------------------------------------------------------------       
    _var_units_map = {
        'atmosphere_bottom_air__mass-per-volume_density': 'kg m-3',
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'J kg-1 K-1', # (see Notes above)
        'atmosphere_bottom_air__temperature': 'deg_C',  # (see Notes above)
        'atmosphere_bottom_air_water-vapor__relative_saturation': '1',
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
        'glacier_ice__mass-specific_isobaric_heat_capacity': 'J kg-1 K-1',
        'cryosphere__melt_volume_flux': 'm s-1' 
          }

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
                          self.is_scalar('RH'),
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

        #------------------------------------------------        
        # enforce_max_meltrate() is always called after
        # update_meltrate() and also enforces min=0
        #------------------------------------------------
        ### np.maximum(M, np.float64(0), M)

        #------------------------------------------
        # Here, the "fill" method works whether M
        # is a 0D array or a scalar.
        #------------------------------------------
        # Could use update_var() in BMI_base also
        #------------------------------------------
        if (np.size(self.SM) == 1):
            M = np.float64(M)  # avoid type change
            self.SM.fill( M )
        else:
            self.SM[:] = M

        #-----------------------------------------------------------
        # Note: enforce_max_snow_meltrate() method is always called
        #       by update() in the base class to make sure that
        #       meltrate does not exceed the max possible.
        #-----------------------------------------------------------
        #  self.enforce_max_snow_meltrate()

    #   update_snow_meltrate()
    #---------------------------------------------------------------------
    def update_snowfall_cold_content(self):
        new_h_snow = (self.P_snow * self.dt) * self.ws_density_ratio
        
        #----------------------------------------------------
        # Copy previous timestep's CC and adjust from here
        #----------------------------------------------------
        Eccs = self.Eccs 

        #----------------------------------------------------
        # Prepare to adjust CC for land surface energy fluxes
        #----------------------------------------------------
        E_in = (self.Q_sum * self.dt)  # [J m-2]

        #--------------------------------------------------------------------
        # For newly fallen snow, add cold content using
        # the same equation used for initializing cold content,
        # but use wet bulb temperature as T_snow for new snow:
        #--------------------------------------------------------------------
        # Wet bulb temp. equation from Stull 2011. Adapted from R code:
        # https://github.com/SnowHydrology/humidity/blob/master/R/humidity.R
        #--------------------------------------------------------------------
        T_wb = self.T_air * np.arctan(0.151977 * ( (self.RH + 8.313659) ** 0.5)) + \
            np.arctan(self.T_air + self.RH) - \
            np.arctan(self.RH - 1.676331) + \
            ((0.00391838 * (self.RH ** 1.5)) * np.arctan(0.023101 * self.RH)) - \
            4.86035
        
        del_T    = (self.T0_cc - T_wb)

        #----------------------------------------------------
        # Only where NEW snow has fallen (P_snow > 0), ADD 
        # cold content for the new snow AND account for land 
        # surface energy fluxes
        #----------------------------------------------------
        Eccs = np.where((self.P_snow > 0), 
                        (np.maximum((Eccs + ((self.rho_snow * self.Cp_snow) * new_h_snow * del_T) - E_in), np.float64(0))),
                        Eccs) # make sure signs check out
        # Eccs = np.maximum(Eccs, np.float64(0)) # make sure signs check out

        if (np.size(self.Eccs) == 1):
            Eccs = np.float64(Eccs)  # avoid type change
            self.Eccs.fill( Eccs )
        else:
            self.Eccs[:] = Eccs

        # update_snowfall_cold_content()
    #---------------------------------------------------------------------
    def update_snowpack_cold_content(self):

        #-----------------------------------------------------
        # Copy the CC that was only adjusted in places WITH 
        # new snowfall before adjusting in places WITHOUT new 
        # snowfall 
        #-----------------------------------------------------   
        Eccs = self.Eccs 

        #------------------------------------------------------
        # In places where no new snow fell (P_snow !> 0), ONLY
        # adjust CC for land surface energy fluxes
        #------------------------------------------------------
        E_in = (self.Q_sum * self.dt)  # [J m-2]
        # Eccs  = np.maximum((self.Eccs - E_in), np.float64(0))
        Eccs = np.where((self.P_snow <= 0), 
                        np.maximum((Eccs - E_in), np.float64(0)), 
                        Eccs) 

        #--------------------------------------------------
        # Where there is no snow, there is no cold content.
        # Set Eccs to 0:
        #--------------------------------------------------
        Eccs = np.where((self.h_snow == 0), np.float64(0), Eccs)       
        
        if (np.size(self.Eccs) == 1):
            Eccs = np.float64(Eccs)  # avoid type change
            self.Eccs.fill( Eccs )
        else:
            self.Eccs[:] = Eccs

    #   update_snowpack_cold_content()
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

        #------------------------------------------------        
        # enforce_max_ice_meltrate() is always called after
        # update_ice_meltrate() and also enforces min=0
        #------------------------------------------------
        ### np.maximum(M, np.float64(0), M)

        #------------------------------------------
        # Here, the "fill" method works whether M
        # is a 0D array or a scalar.
        #------------------------------------------
        # Could use update_var() in BMI_base also
        #------------------------------------------
        IM = np.maximum(M, np.float64(0))
        self.IM = np.where((self.h_swe == 0) & (self.previous_swe == 0), IM, np.float64(0))

        #--------------------------------------------------
        # Update the cold content of the ice [J m-2]
        # If this is positive, there was no melt so far.
        #--------------------------------------------------
        Ecci = np.maximum((self.Ecci - E_in), np.float64(0))

        #--------------------------------------------------
        # Where there is no ice, there is no cold content.
        # Set Ecci to 0:
        #--------------------------------------------------
        Ecci = np.where((self.h_ice == 0), np.float64(0), Ecci)   

        if (np.size(self.Ecci) == 1):
            Ecci = np.float64(Ecci)  # avoid type change
            self.Ecci.fill( Ecci )
        else:
            self.Ecci[:] = Ecci

        #-----------------------------------------------------------
        # Note: enforce_max_ice_meltrate() method is always called
        #       by update() in the base class to make sure that
        #       meltrate does not exceed the max possible.
        #-----------------------------------------------------------
        #  self.enforce_max_ice_meltrate()
            
    #   update_ice_meltrate()
    #---------------------------------------------------------------------
    # def update_combined_meltrate(self):
    #     #-------------------------------------------------------------
    #     # We want to feed combined snow and ice melt to GIUH for 
    #     # runoff, so combine the IM and SM variables to create Mtotal.
    #     #-------------------------------------------------------------

    #     M_total = self.IM + self.SM

    #     self.M_total = M_total

    # #   update_combined_meltrate()
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

        T0 = model_input.read_next(self.T0_unit, self.T0_type, rti)
        if (T0 is not None):
            self.update_var( 'T0', T0 )
        
        h0_snow = model_input.read_next(self.h0_snow_unit, self.h0_snow_type, rti)
        if (h0_snow is not None):
            self.update_var( 'h0_snow', h0_snow )

        rho_snow = model_input.read_next(self.rho_snow_unit, self.rho_snow_type, rti)
        if (rho_snow is not None):
            self.update_var( 'rho_snow', rho_snow )
        
        h0_swe = model_input.read_next(self.h0_swe_unit, self.h0_swe_type, rti)
        if (h0_swe is not None):
            self.update_var( 'h0_swe', h0_swe )

        Cp_ice = model_input.read_next(self.Cp_ice_unit, self.Cp_ice_type, rti)
        if (Cp_ice is not None):
            self.update_var( 'Cp_ice', Cp_ice)
        
        h0_ice = model_input.read_next(self.h0_ice_unit, self.h0_ice_type, rti)
        if (h0_ice is not None):
            self.update_var( 'h0_ice', h0_ice)

        rho_ice = model_input.read_next(self.rho_ice_unit, self.rho_ice_type, rti)
        if (rho_ice is not None):
            self.update_var( 'rho_ice', rho_ice)
            
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