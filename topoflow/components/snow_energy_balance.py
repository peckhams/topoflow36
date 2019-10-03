#
#  Copyright (c) 2001-2014, Scott D. Peckham
#
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
#  NOTES:  This file defines an Energy-Balance snowmelt component
#          and related functions.  It inherits from the snowmelt
#          "base class" in "snow_base.py".
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
#      initialize_cold_content()   # (9/13/14, from snow_base.py)
#      initialize_computed_vars()  # (9/13/14, from snow_base.py)
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
    #       vol_SM was "basin_cumulative_snow_meltwater_volume"
    #------------------------------------------------------------
    _output_var_names = [
        'model__time_step',                            # dt   
        'snowpack__domain_time_integral_of_melt_volume_flux',   # vol_SM
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
        'snowpack__domain_time_integral_of_melt_volume_flux': 'vol_SM',
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
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'J kg-1 K-1 ', # (see Notes above)
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
    def initialize_computed_vars(self):

        #------------------------------------------
        # If T_air or precip are grids, then make
        # sure that h_snow and h_swe are grids
        #------------------------------------------
        T_IS_GRID = self.is_grid('T_air')
        P_IS_GRID = self.is_grid('P_snow')
        H0_SNOW_IS_SCALAR = self.is_scalar('h0_snow')
        H0_SWE_IS_SCALAR  = self.is_scalar('h0_swe') 

        #------------------------------------------------------
        # If h0_snow or h0_swe are scalars, the use of copy()
        # here requires they were converted to numpy scalars.
        # Using copy() may not be necessary for scalars.
        #------------------------------------------------------
        h_snow = self.h0_snow.copy()    # [meters]
        h_swe  = self.h0_swe.copy()     # [meters]

        #------------------------------------------------------       
        # For the Energy Balance method, SM, h_snow and h_swe
        # are always grids because Q_sum is always a grid.
        #--------------------------------------------------------------
        # Convert both h_snow and h_swe to grids if not already grids
        #--------------------------------------------------------------
        if (H0_SNOW_IS_SCALAR):
            self.h_snow = h_snow + np.zeros([self.ny, self.nx], dtype='float64')
        else:
            self.h_snow = h_snow  # (is already a grid)
        if (H0_SWE_IS_SCALAR):
            self.h_swe = h_swe + np.zeros([self.ny, self.nx], dtype='float64')
        else:
            self.h_swe = h_swe    # (is already a grid)          

        self.SM     = np.zeros([self.ny, self.nx], dtype='float64')
        self.vol_SM = self.initialize_scalar( 0, dtype='float64') # (m3)

        #----------------------------------------------------
        # Compute density ratio for water to snow.
        # rho_H2O is for liquid water close to 0 degrees C.
        # Water is denser than snow, so density_ratio > 1.
        #----------------------------------------------------
        self.density_ratio = (self.rho_H2O / self.rho_snow)
                                                       
        #----------------------------------------------------
        # Initialize the cold content of snowpack (2/21/07)
        #-------------------------------------------------------------
        # This is the only difference from initialize_computed_vars()
        # method in snow_base.py.
        #-------------------------------------------------------------
        self.initialize_cold_content()
        
    #   initialize_computed_vars()
    #---------------------------------------------------------------------
    def initialize_cold_content( self ):

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
        self.Ecc  = (self.rho_snow * self.Cp_snow) * self.h0_snow * del_T

        #------------------------------------        
        # Cold content must be nonnegative.
        #----------------------------------------------
        # Ecc > 0 if (T_snow < T0).  i.e. T_snow < 0.
        #----------------------------------------------
        self.Ecc = np.maximum( self.Ecc, np.float64(0))
        ### np.maximum( self.Ecc, np.float64(0), self.Ecc)  # (in place)
        
    #   initialize_cold_content()
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
        E_in  = (self.Q_sum * self.dt)
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
        # Lf = latent heat of fusion [J/kg]
        # Lv = latent heat of vaporization [J/kg]
        # M  = (Qm / (rho_w * Lf))
        #------------------------------------------
        # rho_w = 1000d       ;[kg/m^3]
        # Lf    = 334000d     ;[J/kg = W*s/kg]
        # So (rho_w * Lf) = 3.34e+8  [J/m^3]
        #------------------------------------------
        M       = (Qm / np.float64(3.34E+8))   #[m/s]
        self.SM = np.maximum(M, np.float64(0))

        #--------------------------------------------------
        # Update the cold content of the snowpack [J m-2]
        # If this is positive, there was no melt so far.
        #--------------------------------------------------
        # (9/13/14) Bug fix: Ecc wasn't stored into self.
        #--------------------------------------------------        
        self.Ecc = np.maximum((self.Ecc - E_in), np.float64(0))

        #-------------------------------------------------------
        # Note: enforce_max_meltrate() method is always called
        #       by the base class to make sure that meltrate
        #       does not exceed the max possible.
        #-------------------------------------------------------
        #  self.enforce_max_meltrate()
            
    #   update_meltrate()
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
        # All grids are assumed to have a data type of Float32.
        #--------------------------------------------------------
        Cp_snow = model_input.read_next(self.Cp_snow_unit, self.Cp_snow_type, rti)
        if (Cp_snow is not None): self.Cp_snow = Cp_snow
        
        rho_snow = model_input.read_next(self.rho_snow_unit, self.rho_snow_type, rti)
        if (rho_snow is not None): self.rho_snow = rho_snow

        T0 = model_input.read_next(self.T0_unit, self.T0_type, rti)
        if (T0 is not None): self.T0 = T0
        
        h0_snow = model_input.read_next(self.h0_snow_unit, self.h0_snow_type, rti)
        if (h0_snow is not None): self.h0_snow = h0_snow
        
        h0_swe = model_input.read_next(self.h0_swe_unit, self.h0_swe_type, rti)
        if (h0_swe is not None): self.h0_swe = h0_swe
        
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):
     
        if (self.Cp_snow_type  != 'Scalar'): self.Cp_snow_unit.close()
        if (self.rho_snow_type != 'Scalar'): self.rho_snow_unit.close()
        if (self.T0_type       != 'Scalar'): self.T0_unit.close()
        if (self.h0_snow_type  != 'Scalar'): self.h0_snow_unit.close()
        if (self.h0_swe_type   != 'Scalar'): self.h0_swe_unit.close()
            
##        if (self.T0_file       != ''): self.T0_unit.close()
##        if (self.rho_snow_file != ''): self.rho_snow_unit.close()
##        if (self.h0_snow_file  != ''): self.h0_snow_unit.close()
##        if (self.h0_swe_file   != ''): self.h0_swe_unit.close()

    #   close_input_files()    
    #-------------------------------------------------------------------    
#-------------------------------------------------------------------------  
###-------------------------------------------------------------------------
##def Energy_Balance_Meltrate(Qn_SW, Qn_LW, T_air, T_surf, RH, p0, \
##                            uz, z, z0_air, rho_air, Cp_air, Ecc, \
##                            h_snow, rho_snow, Cp_snow, dt, \
##                            e_air, e_surf):  #(returned)
##
##    #-----------------------------------------------------------------
##    # Notes: 3/13/07.  This function used to have vapor pressure
##    #        arguments e_air and e_surf.  However, e_air is always
##    #        computed as a function of T_air and RH and e_surf is
##    #        computed as a function of T_surf (saturated vap. press.)
##    #        So it makes sense to remove these two arguments and add
##    #        RH (relative humidity).  This change only affects the
##    #        Latent_Heat_Flux function call, which calls a new
##    #        function called Vapor_Pressure.
##    #-----------------------------------------------------------------
##    #        Qm    = energy used to melt snowpack (if > 0)
##    #        Qn_SW = net shortwave radiation flux (solar)
##    #        Qn_LW = net longwave radiation flux (air, surface)
##    #        Qh    = sensible heat flux from turbulent convection
##    #                between snow surface and air
##    #        Qe    = latent heat flux from evaporation, sublimation,
##    #                and condensation
##    #        Qa    = energy advected by moving water (i.e. rainfall)
##    #                (ARHYTHM assumes this to be negligible; Qa=0.)
##    #        Qc    = energy flux via conduction from snow to soil
##    #                (ARHYTHM assumes this to be negligible; Qc=0.)
##    #        Ecc   = cold content of snowpack = amount of energy
##    #                needed before snow can begin to melt [J/m^2]
##
##    #        All Q's here have units of [W/m^2].
##    #        Are they all treated as positive quantities ?
##
##    #        rho_air  = density of air [kg/m^3]
##    #        rho_snow = density of snow [kg/m^3]
##    #        Cp_air   = specific heat of air [Jkg-1 K-1]
##    #        Cp_snow  = heat capacity of snow [Jkg-1 K-1]
##    #                 = ???????? = specific heat of snow
##    #        Kh       = eddy diffusivity for heat [m^2/s]
##    #        Ke       = eddy diffusivity for water vapor [m^2/s]
##    #        Lv       = latent heat of vaporization [J/kg]
##    #        Lf       = latent heat of fusion [J/kg]
##    #        ------------------------------------------------------
##    #        Dn       = bulk exchange coeff for the conditions of
##    #                   neutral atmospheric stability [m/s]
##    #        Dh       = bulk exchange coeff for heat
##    #        De       = bulk exchange coeff for vapor
##    #        ------------------------------------------------------
##    #        T_air    = air temperature [deg_C]
##    #        T_surf   = surface temperature [deg_C]
##    #        T_snow   = average snow temperature [deg_C]
##    #        RH       = relative humidity [none] (in [0,1])
##    #        e_air    = air vapor pressure at height z [mbar]
##    #        e_surf   = surface vapor pressure [mbar]
##    #        ------------------------------------------------------
##    #        h_snow   = snow depth [m]
##    #        z        = height where wind speed is uz [m]
##    #        uz       = wind speed at height z [m/s]
##    #        P0       = atmospheric pressure [mbar]
##    #        T0       = snow temperature when isothermal [deg_C]
##    #                   (This is usually 0.)
##    #        z0_air   = surface roughness length scale [m]
##    #                   (includes vegetation not covered by snow)
##    #                   (Values from page 1033: 0.0013, 0.02 [m])
##    #        kappa    = von Karman's constant [unitless] = 0.41
##    #        dt       = snowmelt timestep [seconds]
##    #----------------------------------------------------------------
##   
##    # FORWARD_FUNCTION Richardson_Number
##    
##    #---------------------------------
##    #Some required physical constants
##    #are defined in the functions:
##    #e.g. Lv, Lf
##    #---------------------------------
##    
##    #------------------------------
##    #Compute the Richardson number
##    #------------------------------
##    Ri = Richardson_Number(z, uz, T_air, T_surf)
##    
##    #-------------------------------------------------
##    #Compute bulk exchange coeffs (neutral stability)
##    #-------------------------------------------------
##    Dn = Bulk_Exchange_Coeff(uz, z, h_snow, z0_air, T_air, T_surf)
##    Dh = Dn
##    De = Dn
##    
##    #---------------------------
##    #Compute sensible heat flux
##    #---------------------------
##    Qh = Sensible_Heat_Flux(rho_air, Cp_air, Dh, T_air, T_surf)
##    #Formula:  Qh = rho_air * Cp_air * Dh * (T_air - T_surf)
##    #print,'Dh = ', Dh
##    #print,'Qh = ', Qh
##    
##    #-------------------------
##    #Compute latent heat flux
##    #-------------------------
##    Qe = Latent_Heat_Flux(rho_air, De, T_air, T_surf, RH, p0,
##                          e_air, e_surf)  #(these 2 returned)
##    #Formula:  Qe = rho_air * Lv * De * (0.662/p0) * (e_air - e_surf)
##    
##    #print,'Qe = ', Qe
##    
##    #-----------------------------
##    #Compute conduction heat flux
##    #-----------------------------
##    Qc = Conduction_Heat_Flux()
##    #Formula:  Qc = 0d
##    
##    #-----------------------------
##    #Compute advective heat flux
##    #-----------------------------
##    Qa = Advection_Heat_Flux()
##    #Formula:  Qa = 0d
##    
##    #---------------------------------
##    #Qn_SW, Qn_SW & Ecc are pointers,
##    #others are local variables
##    #----------------------------------------------------
##    #Ecc is initialized with the Initial_Cold_Content
##    #function by Initialize_Snow_Vars function (2/21/07)
##    #----------------------------------------------------
##    #The following pseudocode only works for scalars but
##    #is otherwise equivalent to that given below and
##    #clarifies the logic:
##    #----------------------------------------------------
##    #  if (Q_sum gt 0) then begin
##    #      if ((Q_sum * dt) gt Ecc) then begin
##    #          ;-----------------------------------------
##    #          ;Snow is melting.  Use some of Q_sum to
##    #          ;overcome Ecc, and remainder to melt snow
##    #          ;-----------------------------------------
##    #          Qm  = Q_sum - (Ecc/dt)
##    #          Ecc = 0
##    #          M   = (Qm / (rho_w * Lf))
##    #      endif else begin
##    #          ;----------------------------
##    #          ;Snow is warming; reduce Ecc
##    #          ;----------------------------
##    #          Ecc = (Ecc - (Q_sum * dt))
##    #          M   = 0d
##    #      endelse
##    #  endif else begin
##    #      ;------------------------------
##    #      ;Snow is cooling; increase Ecc
##    #      ;------------------------------
##    #      Ecc = Ecc - (Q_sum * dt)
##    #      M   = 0d
##    #  endelse
##    #-------------------------------------------------------
##    Q_sum = Qn_SW + Qn_LW + Qh + Qe + Qa + Qc    #[W/m^2]
##    Qcc  = (Ecc / dt)                            #[W/m^2]
##    Qm   = maximum((Q_sum - Qcc), float64(0))                      #[W/m^2]
##    Ecc  = maximum((Ecc - (Q_sum * dt)), float64(0))              #[J/m^2]
##    #print,'Qm = ', Qm
##    #print,' '
##    
##    #-----------------------------------
##    #Convert melt energy to a melt rate
##    #----------------------------------------
##    #Lf = latent heat of fusion [J/kg]
##    #Lv = latent heat of vaporization [J/kg]
##    #M  = (Qm/ (rho_w * Lf))
##    #----------------------------------------
##    #rho_w = 1000d       ;[kg/m^3]
##    #Lf    = 334000d     ;[J/kg = W*s/kg]
##    #So (rho_w * Lf) = 3.34e+8  [J/m^3]
##    #-------------------------------------
##    M = (Qm / float32(3.34E+8))   #[m/s]
##    
##    return maximum(M, float32(0.0))
##
###   Energy_Balance_Meltrate
###-------------------------------------------------------------------------
