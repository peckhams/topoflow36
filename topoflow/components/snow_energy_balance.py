"""
This file defines an "energy balance" snowmelt component and related
functions.  It inherits from the snowmelt "base class" in
"snow_base.py".
"""
#-----------------------------------------------------------------------
#
#  Copyright (c) 2001-2023, Scott D. Peckham
#
#  Sep 2023.  Checked sign in initialize_cold_content().
#             Separate function: update_cold_content().
#             Moved initialize_cold_content() back to base class.
#             Renamed T0 to avoid conflict w/ degree-day method.
#             Removed separate version of initialize_computed_vars().
#  Aug 2023.  Removed trailing space in 'J kg-1 K-1 '.
#             Added "vol_swe_start".
#             Added missing "vol_swe" in initialize_computed_vars.
#             Bug fix: "domain_time_integral_of_liquid-equivalent_depth"
#               "domain_integral_of_liquid-equivalent_depth"
#  Sep 2014.  Cleanup and testing.
#             Own versions of input file routines, at end.
#  Aug 2014.  Customized initialize_computed_vars(), which calls
#             initialize_cold_content(). 
#             Updates to standard names and BMI.
#  Nov 2013.  Converted TopoFlow to Python package.
#  Jan 2013.  Revised handling of input/output names.
#  Oct 2012.  CSDMS Standard Names and BMI.
#  May 2010.  Changes to unit_test() and read_cfg_file().
#  Aug 2009.  Updates.
#  Jul 2009.  Updates.
#  Jan 2009.  Converted from IDL.
#
#-----------------------------------------------------------------------
#
#  class snow_energy_balance
#
#      get_component_name()
#      get_attribute()          # (10/26/11)
#      get_input_var_names()    # (5/14/12)
#      get_output_var_names()   # (5/14/12)
#      get_var_name()           # (5/16/12, Bolton)
#      get_var_units()          # (5/16/12, Bolton)
#      ----------------------
#      check_input_types()
#      initialize_input_file_vars()   # (7/3/20)
#      update_meltrate()
#      update_cold_content()
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
        'model_name':         'TopoFlow_Snowmelt_Energy_Balance',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'SnowEnergyBalance',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Snow_Energy_Balance.cfg.in',
        'cfg_extension':      '_snow_energy_balance.cfg',
        'cmt_var_prefix':     '/SnowEnergyBalance/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Snow_Energy_Balance.xml',
        'dialog_title':       'Snowmelt: Energy Balance Parameters',  # (Snowmelt ?)
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
        'land_surface__temperature',    # (used to initialize Ecc)
        'land_surface_net-total-energy__energy_flux',
        'water-liquid__mass-per-volume_density' ]

        #------------------------------------------------------------
        # These are used by Meteorology component to compute Q_sum.
        # but are not needed directly by this component. (9/14/14)
        #------------------------------------------------------------
        # 'atmosphere_bottom_air__pressure',
        # 'atmosphere_bottom_air_flow__log_law_roughness_length',
        # 'atmosphere_bottom_air_flow__speed_reference_height',
        # 'atmosphere_bottom_air_flow__reference-height_speed',
        # 'atmosphere_bottom_air_water-vapor__relative_saturation',
        # 'land_surface_net-longwave-radiation__energy_flux',
        # 'land_surface_net-shortwave-radiation__energy_flux',
                        
                                  
    #------------------------------------------------------------
    # Note: The main output of the Energy-Balance method is SM.
    #       Functions in snow_base.py compute additional output
    #       vars, such as:
    #       update_SM_integral(), update_depth(), update_swe().
    #       They need things like: rho_H2O and rho_snow.
    #------------------------------------------------------------
    # Note: Cp_snow is a constant set in the "set_constants()"
    #       function in snow_base.py.
    #------------------------------------------------------------
    _output_var_names = [
        'model__time_step',                            # dt   
        'snowpack__domain_time_integral_of_melt_volume_flux',   # vol_SM
        'snowpack__initial_domain_integral_of_liquid-equivalent_depth', # vol_swe_start'
        'snowpack__domain_integral_of_liquid-equivalent_depth',         # vol_swe
        'snowpack__energy-per-area_cold_content',      # Ecc
        'snowpack__depth',                             # h_snow
        'snowpack__initial_depth',                     # h0_snow
        'snowpack__initial_liquid-equivalent_depth',   # h0_swe
        'snowpack__liquid-equivalent_depth',           # h_swe
        'snowpack__melt_volume_flux',                  # SM
        'snowpack__z_mean_of_mass-per-volume_density', # rho_snow 
        'snowpack__z_mean_of_mass-specific_isobaric_heat_capacity' ]  # Cp_snow   
    
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
        # but are not needed directly by this component. (9/14/14)
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
        'snowpack__energy-per-area_cold_content': 'Ecc',
        'snowpack__initial_depth': 'h0_snow',
        'snowpack__initial_liquid-equivalent_depth': 'h0_swe',
        'snowpack__liquid-equivalent_depth': 'h_swe',
        'snowpack__melt_volume_flux': 'SM',
        'snowpack__z_mean_of_mass-per-volume_density': 'rho_snow',
        'snowpack__z_mean_of_mass-specific_isobaric_heat_capacity': 'Cp_snow' }
    
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
        # but are not needed directly by this component. (9/14/14)
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
        'snowpack__z_mean_of_mass-specific_isobaric_heat_capacity': 'J kg-1 K-1' }

    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Snow_Energy_Balance'

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
                          self.is_scalar('rho_air'),
                          self.is_scalar('Cp_air'),
                          self.is_scalar('T_air'),   ##### CHECK THIS ONE.
                          self.is_scalar('T_surf'),
                          self.is_scalar('Q_sum'),
                          #--------------------------------                          
#                           self.is_scalar('RH'),
#                           self.is_scalar('p0'),
#                           self.is_scalar('uz'),
#                           self.is_scalar('z'),
#                           self.is_scalar('z0_air'),
#                           self.is_scalar('Qn_SW'),
#                           self.is_scalar('Qn_LW'),
                          #--------------------------------
                          self.is_scalar('rho_snow'),
                          self.is_scalar('Cp_snow'),
                          self.is_scalar('h0_snow'),
                          self.is_scalar('h0_swe') ])


        self.ALL_SCALARS = np.all(are_scalars)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def initialize_input_file_vars(self):     
    
        #---------------------------------------------------
        # Initialize vars to be read from files (11/16/16)
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
    
    #   initialize_input_file_vars()
    #-------------------------------------------------------------------
    def update_meltrate(self):

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
        # Ecc is initialized by initialize_cold_content().
        #------------------------------------------------------
        # The following pseudocode only works for scalars but
        # is otherwise equivalent to that given below and
        # clarifies the logic:
        #------------------------------------------------------
        #  if (Q_sum gt 0) then begin
        #      if ((Q_sum * dt) gt Ecc) then begin
        #          ;-------------------------------------------
        #          ; Snow is melting.  Use some of Q_sum to
        #          ; overcome Ecc, and remainder to melt snow
        #          ;-------------------------------------------
        #          Qm  = Q_sum - (Ecc/dt)
        #          Ecc = 0
        #          M   = (Qm / (rho_w * Lf))
        #      endif else begin
        #          ;------------------------------
        #          ; Snow is warming; reduce Ecc
        #          ;------------------------------
        #          Ecc = (Ecc - (Q_sum * dt))
        #          M   = 0d
        #      endelse
        #  endif else begin
        #      ;--------------------------------
        #      ; Snow is cooling; increase Ecc
        #      ;--------------------------------
        #      Ecc = Ecc - (Q_sum * dt)
        #      M   = 0d
        #  endelse
        #---------------------------------------------------------
        # Q_sum = Qn_SW + Qn_LW + Qh + Qe + Qa + Qc    # [W m-2]
        #---------------------------------------------------------
        
        #-----------------------------------------------        
        # New approach; easier to understand (9/14/14)
        #----------------------------------------------- 
        # E_in  = energy input over one time step
        # E_rem = energy remaining in excess of Ecc
        #-----------------------------------------------
        # 1 Watt = 1 Joule / second
        # W m-2 = J s-1 m-2
        #-----------------------------------------------        
        E_in  = (self.Q_sum * self.dt)  # [J m-2]
        E_rem = np.maximum( E_in - self.Ecc, np.float64(0) )
        Qm    = (E_rem / self.dt)  # [W m-2]
       
        ##################################
        # Used before 9/14/14/.
        ##################################        
        # Q_sum = self.Q_sum  # (2/3/13, new framework)
        # Qcc   = (self.Ecc / self.dt)                   # [W m-2]
        # Qm    = np.maximum((Q_sum - Qcc), float64(0))  # [W m-2]
        
        #-------------------------------------
        # Convert melt energy to a melt rate
        #------------------------------------------
        # Lf = latent heat of fusion [J/kg]       (solid to liquid)
        # Lv = latent heat of vaporization [J/kg] (liquid to gas)
        # M  = (Qm / (rho_w * Lf))
        # [m s-1] = [W m-2] / ([kg m-3]*[J kg-1])
        # [m s-1] = [J s-1 m-2] * m3 J-1
        #------------------------------------------
        # rho_w = 1000d       ;[kg/m^3]
        # Lf    = 334000d     ;[J/kg = W*s/kg]
        # So (rho_w * Lf) = 3.34e+8  [J/m^3]
        #------------------------------------------
        M = (Qm / np.float64(3.34E+8))   #[m/s]

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

        #--------------------------------------------------------
        # Note: enforce_max_meltrate() method is always called
        #       by update() in the base class to make sure that
        #       meltrate does not exceed the max possible.
        #--------------------------------------------------------
        # self.enforce_max_meltrate()
            
    #   update_meltrate()
    #-------------------------------------------------------------------
    def update_cold_content(self):

        #--------------------------------------------------
        # Update the cold content of the snowpack [J m-2]
        # If this is positive, there was no melt so far.
        #--------------------------------------------------
        # (9/13/14) Bug fix: Ecc wasn't stored into self.
        #-----------------------------------------------------
        # Recall that Ecc was initialized as:
        #    Ecc  = (rho_snow * Cp_snow) * h0_snow * del_T
        # Why doesn't Ecc still depend on h0_snow and del_T?
        #-----------------------------------------------------        
        ## self.Ecc = np.maximum((self.Ecc - E_in), np.float64(0))
        E_in = (self.Q_sum * self.dt)  # [J m-2]
        Ecc  = np.maximum((self.Ecc - E_in), np.float64(0))
        if (np.size(self.Ecc) == 1):
            Ecc = np.float64(Ecc)  # avoid type change
            self.Ecc.fill( Ecc )
        else:
            self.Ecc[:] = Ecc

    #   update_cold_content()
    #---------------------------------------------------------------------
    def open_input_files(self):

        self.Cp_snow_file  = self.in_directory + self.Cp_snow_file
        self.rho_snow_file = self.in_directory + self.rho_snow_file
        self.T0_file       = self.in_directory + self.T0_file
        self.h0_snow_file  = self.in_directory + self.h0_snow_file
        self.h0_swe_file   = self.in_directory + self.h0_swe_file

        self.Cp_snow_unit  = model_input.open_file(self.Cp_snow_type,  self.Cp_snow_file)
        self.rho_snow_unit = model_input.open_file(self.rho_snow_type, self.rho_snow_file)
        self.T0_unit       = model_input.open_file(self.T0_type,       self.T0_file)
        self.h0_snow_unit  = model_input.open_file(self.h0_snow_type,  self.h0_snow_file)
        self.h0_swe_unit   = model_input.open_file(self.h0_swe_type,   self.h0_swe_file)

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
                    
        rho_snow = model_input.read_next(self.rho_snow_unit, self.rho_snow_type, rti)
        if (rho_snow is not None):
            self.update_var( 'rho_snow', rho_snow )
  
        h0_snow = model_input.read_next(self.h0_snow_unit, self.h0_snow_type, rti)
        if (h0_snow is not None):
            self.update_var( 'h0_snow', h0_snow )
        
        h0_swe = model_input.read_next(self.h0_swe_unit, self.h0_swe_type, rti)
        if (h0_swe is not None):
            self.update_var( 'h0_swe', h0_swe )
        
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):
     
        if (self.Cp_snow_type  != 'Scalar'): self.Cp_snow_unit.close()
        if (self.rho_snow_type != 'Scalar'): self.rho_snow_unit.close()
        if (self.T0_type       != 'Scalar'): self.T0_unit.close()
        if (self.h0_snow_type  != 'Scalar'): self.h0_snow_unit.close()
        if (self.h0_swe_type   != 'Scalar'): self.h0_swe_unit.close()

    #   close_input_files()    
    #-------------------------------------------------------------------  

