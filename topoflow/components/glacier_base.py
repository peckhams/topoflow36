"""
This file defines a "base class" for combined glacier melt and snowmelt components as well
as functions used by most or all glacier melt and snowmelt methods.  That is, all
glacier melt/snowmelt components inherit methods from this class.  The methods
of this class should be over-ridden as necessary (especially the
"update_<snow/ice>meltrate() method) for different methods of modeling
glacier melt/snowmelt.  This class, in turn, inherits from the "BMI base class"
in BMI_base.py.

See glacier_degree_day.py and glacier_energy_balance.py.
"""
#-----------------------------------------------------------------------
#
#  Copyright (c) 2001-2023, Scott D. Peckham
#
# Oct 2023.  Updated the original snow cold content routine
#            so that new cold content is added with new snow in 
#            update_snowfall_cold_content(), and still accounts
#            for land surface energy fluxes in update_snowpack_cold_content(), 
#            and there is no cold content if there is no snow.
#  Sep 2023. Added update_density_ratio().
#            Fixed bug in enforce_max_meltrate().
#            Moved initialize_cold_content() from snow_energy_balance.py
#            back to this base class.  Note that T0 plays 2 roles in
#            the degree-day component, in meltrate and cold_content.
#            Degree-day comp now has option to set T0_cc in its CFG.
#  Aug 202., Renamed snow functions to have 'snow' in the name
#            Duplicated existing snow functions for 'ice'
#
#-----------------------------------------------------------------------
#  NOTES:  update_snow_vars() in precip.py sets values here. # LB: is this outdated? Is met_base.py meant here? 
#-----------------------------------------------------------------------
#
#  class glacier_component    (inherits from BMI_base.py)
#
#      set_constants()
#      -----------------------
#      initialize()
#      update()
#      finalize()
#      set_computed_input_vars()
#      --------------------------
#      set_missing_cfg_options()   
#      check_input_types()
#      initialize_computed_vars()
#      initialize_snow_cold_content()   
#      initialize_ice_cold_content() 
#      ----------------------------
#      update_snow_meltrate()
#      update_ice_meltrate()
#      enforce_max_snow_meltrate()
#      enfore_max_ice_meltrate()
#      update_SM_integral()
#      update_IM_integral()
#      update_snowfall_cold_content()
#      update_snowpack_cold_content()
#      extract_previous_swe()
#      update_swe()
#      update_density_ratio()
#      update_swe_integral()   
#      update_iwe()
#      update_iwe_integral() 
#      extract_previous_snow_depth()
#      update_snow_depth()
#      update_ice_depth()
#      update_total_snowpack_water_volume()
#      update_total_ice_water_volume()
#      -----------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      -----------------------
#      update_outfile_names()
#      disable_all_output()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.utils import BMI_base
# from topoflow.utils import model_input  # (not used here)
from topoflow.utils import model_output

#-----------------------------------------------------------------------
class glacier_component( BMI_base.BMI_component ):

    def set_constants(self):

        #-----------------------------------
        # Constants not changeable by user
        #---------------------------------------------------------
        # Cp_snow = mass-specific isobaric heat capacity of snow
        # Cp_ice = mass-specific isobaric heat capacity of ice
        #---------------------------------------------------------
        # Lf = latent heat of fusion for water [J kg -1]
        #---------------------------------------------------------        
        self.Cp_snow  = np.float64( 2090.0 )  # [J kg-1 K-1]
        self.Cp_ice       = np.float64(2060.0)     # [J/(kg * K)]
        self.Lf       = np.float64( 334000 )  # [J kg-1] 
    
        #--------------------------------------
        # Not a constant; read from CFG file.
        #--------------------------------------
        ## self.rho_snow = np.float64(300) # [kg m-3]
        ## self.rho_H2O  = np.float64(1000)  # (See initialize() method.)
        ## self.rho_ice      = np.float64(917) # [kg m-3]
                
    #   set_constants()       
    #-------------------------------------------------------------------
    def latent_heat_of_sublimation(self):

        #----------------------------------------------------------    
        # Notes:  See:  http://en.wikipedia.org/wiki/Latent_heat
        #         Valid for T in [-40, 0] deg C.
        #----------------------------------------------------------
        # sublimation/deposition.  What about fusion/melting?
        # deposition, desublimation and resublimation are synonyms
        #----------------------------------------------------------
        a = np.float64( 2834.1 )
        b = np.float64( -0.29 )
        c = np.float64( -0.004 )
        T = self.T_air
        
        self.Lf = a + (b * T) + (c * (T**2))  # [J g-1]
        self.Lf *= np.float64( 1000 ) # [J kg-1]
        
    #   latent_heat_of_sublimation()  
    # #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        self.SILENT = SILENT

        #---------------------------------------------------------
        # Notes:  Need to make sure than h_swe matches h_snow 
        #         User may have entered incompatible values.
        #         The same goes for h_iwe and h_ice.
        #---------------------------------------------------------
        if not(self.SILENT):
            print(' ')
            print('Glacier component: Initializing...')
            
        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode
        self.cfg_file   = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------       
        self.set_constants()
        self.initialize_config_vars() 
        self.set_missing_cfg_options()
        self.initialize_basin_vars()  # (5/14/10)
        #-----------------------------------------
        # This must come before "Disabled" test.
        #-----------------------------------------
        self.initialize_time_vars()
   
        if (self.comp_status.lower() == 'disabled'):
            if not(self.SILENT):
                print('Glacier component: Disabled in CFG file.')
            self.disable_all_output()
            self.h_snow  = self.initialize_scalar(0, dtype='float64')
            self.h_ice = self.initialize_scalar(0, dtype='float64')
            self.h_swe   = self.initialize_scalar(0, dtype='float64')
            self.h_iwe   = self.initialize_scalar(0, dtype='float64')
            self.SM      = self.initialize_scalar(0, dtype='float64')
            self.IM      = self.initialize_scalar(0, dtype='float64')
            self.M_total = self.initialize_scalar(0, dtype='float64')
            self.vol_SM  = self.initialize_scalar(0, dtype='float64') # [m3]
            self.vol_IM = self.initialize_scalar(0, dtype='float64') 
            self.vol_M_total = self.initialize_scalar(0, dtype='float64')
            self.vol_swe = self.initialize_scalar(0, dtype='float64') # [m3]
            self.vol_swe_start = self.initialize_scalar(0, dtype='float64')
            self.vol_iwe = self.initialize_scalar(0, dtype='float64')
            self.vol_iwe_start = self.initialize_scalar(0, dtype='float64')
            self.DONE    = True
            self.status  = 'initialized'
            return
 
        #----------------------------------------
        # Initialize vars to be read from files
        #----------------------------------------
        self.initialize_input_file_vars()
 
        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()

        #---------------------------
        # Initialize computed vars
        #---------------------------
        self.check_input_types()  # (maybe not used yet)
        self.initialize_computed_vars()  # (h_snow, h_swe, etc.)

        self.open_output_files()
        self.status = 'initialized'        
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, dt=-1.0):

        #----------------------------------------------------------
        # Note: The read_input_files() method is first called by
        #       the initialize() method.  Then, the update()
        #       method is called one or more times, and it calls
        #       other update_*() methods to compute additional
        #       variables using input data that was last read.
        #       Based on this pattern, read_input_files() should
        #       be called at end of update() method as done here.
        #       If the input files don't contain any additional
        #       data, the last data read persists by default.
        #----------------------------------------------------------

        #--------------------------------
        # Has component been disabled ?
        #-------------------------------------------------
        # Note: self.SM already set to 0 by initialize()
        #-------------------------------------------------
        if (self.comp_status.lower() == 'disabled'):
            # Note: self.status should be 'initialized'.
            return

        #-----------------------------------------
        # Read next snow vars from input files ?
        #-----------------------------------------------------       
        # Note: read_input_files() is called by initialize()
        # and those values must be used for the "update"
        # calls before reading new ones.
        #-----------------------------------------------------
        self.status = 'updating'
        if (self.time_index > 0):
            self.read_input_files()
                   
        #-------------------------
        # Update computed values 
        #-------------------------
        self.extract_previous_swe()
        self.extract_previous_snow_depth() # used in update_snow_cold_content()

        self.update_snow_meltrate()       # (meltrate = SM)
        self.enforce_max_snow_meltrate()  # (before SM integral!)
        self.update_SM_integral()

        #----------------------------------------------------------
        # Call update_swe/iwe() and before update_snow/ice_depth()
        #----------------------------------------------------------   
        self.update_swe()
        self.update_snowfall_cold_content()
        # self.update_swe_integral() 

        self.update_ice_meltrate()
        self.enforce_max_ice_meltrate()
        self.update_IM_integral()

        self.update_combined_meltrate()

        self.update_iwe() # relies on previous timestep's swe value
        # self.update_iwe_integral() # relies on previous timestep's iwe value   

        self.update_density_ratio()
        # self.update_ice_density_ratio()

        self.update_snow_depth()
        self.update_ice_depth()
        #--------------------------------------------------------------
        # Call update_snowpack_cold_content() after update_snow_depth()
        # so lack of snow or added cold content from new snow can
        # be accounted for
        #--------------------------------------------------------------   
        self.update_snowpack_cold_content() 
        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # Components use own self.time_sec by default.
        #-----------------------------------------------
        self.write_output_files()
        ## self.write_output_files( time_seconds )

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        self.update_time( dt )
        self.status = 'updated'  # (OpenMI)
        
    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        #--------------------------------
        # Has component been disabled ?
        #--------------------------------
        if (self.comp_status.lower() == 'disabled'):
            # Note: self.status should be 'initialized'.
            return

        self.status = 'finalizing' 
        self.close_input_files()   ##  TopoFlow input "data streams"
        self.close_output_files()
        if not(self.SILENT):
            self.print_final_report(comp_name='Glacier component')
        self.status = 'finalized'
               
    #   finalize()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        #---------------------------------------------------------
        np.maximum(self.save_grid_dt,   self.dt, out=self.save_grid_dt)
        np.maximum(self.save_pixels_dt, self.dt, out=self.save_pixels_dt)
        
    #   set_computed_input_vars()        
    #-------------------------------------------------------------------
    def set_missing_cfg_options(self):

        #-------------------------------------------------------  
        # Note: This is called in initialize() AFTER calling
        #       initialize_config_vars().  It is used to set
        #       newer toggles, etc. that may not have been
        #       set in the CFG file.
        #-------------------------------------------------------    
        pass
        
    #   set_missing_cfg_options()            
    #-------------------------------------------------------------------    
    def check_input_types(self):

        #----------------------------------------------------
        # Notes: Usually this will be overridden by a given
        #        method of computing snow meltrate.
        #----------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_ice, Cp_ice, rho_air 
        # and Cp_air are currently always scalars.
        #----------------------------------------------------        
        are_scalars = np.array([
                          self.is_scalar('P_snow'),
                          self.is_scalar('rho_H2O'),
                          self.is_scalar('rho_air'),
                          self.is_scalar('Cp_air'),
                          self.is_scalar('RH'),
                          #----------------------------------
                          self.is_scalar('rho_snow'),
                          self.is_scalar('Cp_snow'),
                          self.is_scalar('h0_snow'),
                          self.is_scalar('h0_swe'),
                          #----------------------------------
                          self.is_scalar('Cp_ice'),
                          self.is_scalar('rho_ice'),
                          self.is_scalar('h0_ice'),
                          self.is_scalar('h0_iwe'),
                          self.is_scalar('h_active_layer') ])


        self.ALL_SCALARS = np.all(are_scalars)
  
    #   check_input_types()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        #------------------------------------------
        # If T_air or precip are grids, then make
        # sure that h_snow, h_swe, h_ice, and h_iwe 
        # are grids
        #------------------------------------------
        T_IS_GRID = self.is_grid('T_air')
        P_IS_GRID = self.is_grid('P_snow')
        RH_IS_GRID = self.is_grid('RH')
        H0_SNOW_IS_SCALAR = self.is_scalar('h0_snow')
        H0_SWE_IS_SCALAR  = self.is_scalar('h0_swe') 
        H0_ICE_IS_SCALAR = self.is_scalar('h0_ice')
        H0_IWE_IS_SCALAR = self.is_scalar('h0_iwe')
        comp_name = self.get_component_name()
        EN_BAL_COMP = (comp_name == 'TopoFlow_Glacier_Energy_Balance')

        #------------------------------------------------------
        # If h0_snow, h0_swe, h0_ice, or h0_iwe are scalars, 
        # the use of copy() here requires they were converted 
        # to numpy scalars. Using copy() may not be necessary 
        # for scalars.
        #------------------------------------------------------
        h_snow = self.h0_snow.copy()    # [meters]
        h_swe  = self.h0_swe.copy()     # [meters]
        h_ice = self.h0_ice.copy()
        h_iwe = self.h0_iwe.copy()
        
        if (T_IS_GRID or P_IS_GRID or EN_BAL_COMP):
            self.SM = np.zeros([self.ny, self.nx], dtype='float64')
            self.IM = np.zeros([self.ny, self.nx], dtype='float64')
            self.M_total = np.zeros([self.ny, self.nx], dtype='float64')
            self.T_surf = np.full((self.ny, self.nx), self.T_surf)
            #-----------------------------------------
            # Convert h_snow, h_swe, h_ice, and h_iwe 
            # to grids if not already grids.
            # For the Energy Balance method, SM, h_snow and h_swe
            # are always grids because Q_sum is always a grid.
            #-----------------------------------------
            if (H0_SNOW_IS_SCALAR):
                self.h_snow = h_snow + np.zeros([self.ny, self.nx], dtype='float64')
            else:
                self.h_snow = h_snow  # (is already a grid)
            #------------------------------------------------
            if (H0_SWE_IS_SCALAR):
                self.h_swe = h_swe + np.zeros([self.ny, self.nx], dtype='float64')
            else:
                self.h_swe = h_swe    # (is already a grid)         
            #------------------------------------------------
            if (H0_ICE_IS_SCALAR):
                self.h_ice = h_ice + np.zeros([self.ny, self.nx], dtype='float64')
            else: 
                self.h_ice = h_ice 
            if (H0_IWE_IS_SCALAR):
                self.h_iwe = h_iwe + np.zeros([self.ny, self.nx], dtype='float64')
            else: 
                self.h_iwe = h_iwe

        else:
            #---------------------------------------------------
            # Both are scalars and that's for snow_degree_day.py
            #---------------------------------------------------
            self.SM     = self.initialize_scalar( 0, dtype='float64')
            self.IM     = self.initialize_scalar( 0, dtype='float64')
            self.M_total= self.initialize_scalar( 0, dtype='float64')
            self.h_snow = self.initialize_scalar( h_snow, dtype='float64')
            self.h_swe  = self.initialize_scalar( h_swe,  dtype='float64')
            self.h_ice  = self.initialize_scalar( h_ice,  dtype='float64')
            self.h_iwe  = self.initialize_scalar( h_iwe,  dtype='float64')

        # vol_swe is to track volume of water in the snowpack.
        # vol_iwe is to track volume of water in the ice column.
        self.vol_SM  = self.initialize_scalar( 0, dtype='float64') # (m3)
        self.vol_swe = self.initialize_scalar( 0, dtype='float64') # (m3)

        self.update_total_snowpack_water_volume()
        self.vol_swe_start = self.vol_swe.copy()
        
        self.vol_IM  = self.initialize_scalar( 0, dtype='float64') # (m3)
        self.vol_iwe = self.initialize_scalar( 0, dtype='float64') # (m3)

        self.update_total_ice_water_volume()
        self.vol_iwe_start = self.vol_iwe.copy()

        self.vol_M_total = self.initialize_scalar( 0, dtype='float64')
        #----------------------------------------------------
        # Compute density ratio for water to snow.
        # rho_H2O is for liquid water close to 0 degrees C.
        # Water is denser than snow, so density_ratio > 1.
        #----------------------------------------------------
        self.ws_density_ratio = (self.rho_H2O / self.rho_snow)

        #----------------------------------------------------
        # Compute density ratio for water to ice.
        # rho_H2O is for liquid water close to 0 degrees C.
        # Water is denser than ice, so density_ratio > 1.
        #----------------------------------------------------
        self.wi_density_ratio = (self.rho_H2O / self.rho_ice)

        #------------------------------------------
        # Initialize the cold content of snowpack
        #------------------------------------------
        self.initialize_snow_cold_content()
        self.initialize_ice_cold_content()
                
    #   initialize_computed_vars()
    #-------------------------------------------------------------------
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
        #        T0 is read from the config file.
        #        The degree-day component has another T0 for meltrate,
        #        but now T0_cc can also be specified in its CFG file.
        #        See it's "set_missing_cfg_options()" method.
        #        For the energy-balance comp, it's T0 is T0_cc.
        #        This is the last var to be set in initialize().
        #---------------------------------------------------------------
        if not(hasattr(self, 'T0_cc')):
           self.T0_cc = self.T0   # synonyms

        #--------------------------------------------
        # Compute initial cold content of snowpack
        # See equation (10) in Zhang et al. (2000).
        #--------------------------------------------
        T_snow   = self.T_surf
        del_T    = (self.T0_cc - T_snow)
        self.Eccs = (self.rho_snow * self.Cp_snow) * self.h0_snow * del_T

        #------------------------------------        
        # Cold content must be nonnegative.
        #----------------------------------------------
        # Ecc > 0 if (T_snow < T0).  i.e. T_snow < 0.
        #----------------------------------------------
        np.maximum( self.Eccs, np.float64(0), out=self.Eccs)  # (in place)
        
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
        #        T0 is read from the config file.
        #        The degree-day component has another T0 for meltrate,
        #        but now T0_cc can also be specified in its CFG file.
        #        See it's "set_missing_cfg_options()" method.
        #        For the energy-balance comp, it's T0 is T0_cc.
        #        This is the last var to be set in initialize().
        #---------------------------------------------------------------
        if not(hasattr(self, 'T0_cc')):
           self.T0_cc = self.T0   # synonyms

        #--------------------------------------------
        # Compute initial cold content of snowpack
        # See equation (10) in Zhang et al. (2000).
        #--------------------------------------------
        T_ice   = self.T_surf
        del_T    = (self.T0_cc - T_ice)
        self.Ecci = (self.rho_ice * self.Cp_ice) * self.h_active_layer * del_T

        #------------------------------------        
        # Cold content must be nonnegative.
        #----------------------------------------------
        # Ecc > 0 if (T_snow < T0).  i.e. T_snow < 0.
        #----------------------------------------------
        np.maximum( self.Ecci, np.float64(0), out=self.Ecci)  # (in place)
        
    #   initialize_ice_cold_content()
    #-------------------------------------------------------------------
    def update_snow_meltrate(self):
        #---------------------------------------------------------
        # Notes: This is for a "potential" meltrate, which can't
        #        be realized unless there is enough snow.
        #        See glacier_base.enforce_max_snow_meltrate().   
        #---------------------------------------------------------
        # Note: We don't need to update any variables if
        #       the snowmelt method is None.  But we need
        #       to make sure that self.SM = 0.0.
        #       This "method" will be over-ridden by a
        #       particular snowmelt method.
        #--------------------------------------------------
        print('ERROR: update_snow_meltrate() method for Glacier component')
        print('       has not been implemented.')
    #-------------------------------------------------------------------
    def update_ice_meltrate(self):
        #---------------------------------------------------------
        # Notes: This is for a "potential" meltrate, which can't
        #        be realized unless there is enough ice.
        #        See glacier_base.enforce_max_ice_meltrate().   
        #---------------------------------------------------------
        print('ERROR: update_ice_meltrate() method for Glacier component')
        print('       has not been implemented.')       
    #   update_snow_meltrate()
    #-------------------------------------------------------------------
    def update_combined_meltrate(self):
        #---------------------------------------------------------
        # We want to feed combined snow and ice melt to GIUH for 
        # runoff, so combine the IM and SM variables to create Mtotal.
        #---------------------------------------------------------
        M_total = self.IM + self.SM

        self.M_total = M_total 
        
    #   update_combined_meltrate()
    #-------------------------------------------------------------------
    def enforce_max_snow_meltrate(self):
    
        #-------------------------------------------------------
        # The max possible meltrate would be if all snow (given
        # by snow depth, h_snow, were to melt in the one time
        # step, dt.  Meltrate should never exceed this value.
        # Recall that: (h_snow / h_swe) = (rho_H2O / rho_snow)
        #                               = density_ratio > 1
        # So h_swe = h_snow / density_ratio.
        # Previous version had a bug; see below. 
        # Now also using "out" keyword for "in-place".  
        #------------------------------------------------------- 
        SM_max = self.h_swe / self.dt
        self.SM = np.minimum(self.SM, SM_max, out=self.SM)  # [m s-1]

        #------------------------------------------------------
        # Make sure meltrate is positive, while we're at it ?
        # Is already done by "Energy-Balance" component.
        #------------------------------------------------------
        np.maximum(self.SM, np.float64(0), out=self.SM)
   
    #   enforce_max_snow_meltrate()
    #-------------------------------------------------------------------
    def enforce_max_ice_meltrate(self):
        
        #-------------------------------------------------------
        # The max possible meltrate would be if all ice (given
        # by ice depth, h_ice, were to melt in the one time
        # step, dt.  Meltrate should never exceed this value.
        #------------------------------------------------------- 
        IM_max = self.h_iwe / self.dt
        self.IM = np.minimum(self.IM, IM_max, out=self.IM)  # [m s-1]

        #------------------------------------------------------
        # Make sure meltrate is positive, while we're at it ?
        # Is already done by "Energy-Balance" component.
        #------------------------------------------------------
        np.maximum(self.IM, np.float64(0), out=self.IM)

    #   enforce_max_ice_meltrate()
    #-------------------------------------------------------------------
    def update_SM_integral(self):

        #------------------------------------------------
        # Update mass total for SM, sum over all pixels
        #------------------------------------------------   
        volume = np.float64(self.SM * self.da * self.dt)  # [m^3]
        if (np.size(volume) == 1):
            self.vol_SM += (volume * self.rti.n_pixels)
        else:
            self.vol_SM += np.sum(volume)  #### np.sum vs. sum ???
    #   update_SM_integral()
    #-------------------------------------------------------------------
    def update_IM_integral(self):

        #------------------------------------------------
        # Update mass total for IM, sum over all pixels
        #------------------------------------------------  
        volume = np.float64(self.IM * self.da * self.dt)
        if (np.size(volume) == 1):
            self.vol_IM += (volume * self.rti.n_pixels)
        else:
            self.vol_IM += np.sum(volume)
            
    #   update_IM_integral()
    #-------------------------------------------------------------------
    def update_snowfall_cold_content(self):

        #-----------------------------------------------------------    
        # 2023-09-25. This is overridden in glacier_energy_balance.py
        # and could do the same in glacier_degree_day.py using some
        # other (temperature-based) method.
        # This is only needed for snow because snow can accumulate
        # and thus add cold content, whereas ice only melts (in 
        # this model)
        #-----------------------------------------------------------
        pass
    #-------------------------------------------------------------------
    def update_snowpack_cold_content(self):

        #-----------------------------------------------------------    
        # 2023-09-25. This is overridden in glacier_energy_balance.py
        # and could do the same in glacier_degree_day.py using some
        # other (temperature-based) method.
        # This is only needed for snow because snow can accumulate
        # and thus add cold content, whereas ice only melts (in 
        # this model)
        #-----------------------------------------------------------
        pass

    #   update_cold_content()
    #-------------------------------------------------------------------   
    def extract_previous_swe(self):
        #------------------------------------------------
        # Extract swe from previous timestep for use in 
        # toggling between ice/snow routines
        #------------------------------------------------
        self.previous_swe = self.h_swe.copy()

    #   extract_previous_swe()
    #-------------------------------------------------------------------
    def update_swe(self):

        #--------------------------------------------------------
        # Note: The Meteorology component uses air temperature
        # to compute P_rain (precip that falls as liquid) and
        # P_snow (precip that falls as snow or ice) separately.
        # P_snow = (self.P * (self.T_air <= 0)) 
        #----------------------------------------------------------
        # Note: This method must be written to work regardless
        # of whether P_rain and T are scalars or grids. (3/14/07)
        #------------------------------------------------------------
        # If P or T_air is a grid, then h_swe and h_snow are grids.
        # This is set up in initialize_computed_vars().
        #------------------------------------------------------------
      
        #------------------------------------------------
        # Increase snow water equivalent due to snowfall
        #------------------------------------------------
        # Meteorology and Channel components may have
        # different time steps, but then self.P_snow
        # will be a time-interpolated value.
        #------------------------------------------------
        dh1_swe  = (self.P_snow * self.dt)
        self.h_swe += dh1_swe
                
        #------------------------------------------------
        # Decrease snow water equivalent due to melting
        # Note that SM depends partly on h_snow.
        #------------------------------------------------
        dh2_swe    = self.SM * self.dt
        self.h_swe -= dh2_swe
        np.maximum(self.h_swe, np.float64(0), self.h_swe)  # (in place)
        
    #   update_swe()
    #-------------------------------------------------------------------
    def update_iwe(self):
                
        #------------------------------------------------
        # Decrease ice water equivalent due to melting 
        #------------------------------------------------
        dh2_iwe    = self.IM * self.dt
        self.h_iwe -= dh2_iwe
        np.maximum(self.h_iwe, np.float64(0), self.h_iwe)  # (in place)

    #   update_iwe()
    #-------------------------------------------------------------------
    def update_density_ratio(self):

        #-----------------------------------------------    
        # Return if density_ratio is constant in time.
        #-----------------------------------------------
        if (self.rho_snow_type.lower() in ['scalar', 'grid']):
            return
        density_ratio = self.rho_H2O / self.rho_snow

        #-------------------------------------             
        # Save updated density ratio in self
        #-------------------------------------
        if (np.ndim( self.density_ratio ) == 0):
            density_ratio = np.float64( density_ratio )  ### (from 0D array to scalar)
            self.density_ratio.fill( density_ratio )     ### (mutable scalar)
        else:
            self.density_ratio[:] = density_ratio

    #   update_density_ratio()
    #-------------------------------------------------------------------
    def update_swe_integral(self):

        #------------------------------------------------
        # Update mass total for water in the snowpack,
        # sum over all pixels.
        #------------------------------------------------   
        volume = np.float64(self.h_swe * self.da)  # [m^3]
        if (np.size(volume) == 1):
            self.vol_swe += (volume * self.rti.n_pixels)
        else:
            self.vol_swe += np.sum(volume)
            
    #   update_swe_integral()   
    #-------------------------------------------------------------------
    def update_iwe_integral(self):

        #------------------------------------------------
        # Update mass total for water in the ice column,
        # sum over all pixels.
        #------------------------------------------------   
        volume = np.float64(self.h_iwe * self.da)  # [m^3]
        if (np.size(volume) == 1):
            self.vol_iwe += (volume * self.rti.n_pixels)
        else:
            self.vol_iwe += np.sum(volume)
            
    #   update_iwe_integral()   
    #-------------------------------------------------------------------
    def extract_previous_snow_depth(self):
        #------------------------------------------------
        # Extract swe from previous timestep for use in 
        # toggling between ice/snow routines
        #------------------------------------------------
        self.previous_h_snow = self.h_snow.copy()

    #   extract_previous_snow_depth()
    #-------------------------------------------------------------------
    def update_snow_depth(self):

        #--------------------------------------------------------
        # Note: The Meteorology component uses air temperature
        # to compute P_rain (precip that falls as liquid) and
        # P_snow (precip that falls as snow or ice) separately.
        # P_snow = (self.P * (self.T_air <= 0)) 
        #----------------------------------------------------------
        # Note: This method must be written to work regardless
        # of whether P_rain and T are scalars or grids.
        #------------------------------------------------------------
        # If P or T_air is a grid, then h_swe and h_snow are grids.
        # This is set up in initialize_computed_vars().
        #------------------------------------------------------------
        # Note that for a region of area, A:
        #     rho_snow = (mass_snow / (h_snow * A))
        #     rho_H2O  = (mass_H20  / (h_swe * A))
        # Since mass_snow = mass_H20 (for SWE):
        #     rho_snow * h_snow = rho_H2O * h_swe
        #     (h_snow / h_swe)  = (rho_H2O / rho_snow)
        #      h_snow = h_swe * density_ratio
        # Since liquid water is denser than snow:
        #      density_ratio > 1 and
        #      h_snow > h_swe
        # self.density_ratio = (self.rho_H2O / self.rho_snow)
        # rho_H2O is for liquid water close to 0 degrees C.
        #------------------------------------------------------------

        #-------------------------------------------------        
        # Change snow depth due to melting or falling snow
        #-------------------------------------------------
        # This assumes that update_swe() is called
        # before update_snow_depth().
        #-------------------------------------------
        h_snow = self.h_swe * self.ws_density_ratio
        
        #----------------------------------             
        # Save updated snow depth in self
        #----------------------------------
        if (np.ndim( self.h_snow ) == 0):
            h_snow = np.float64( h_snow )  ### (from 0D array to scalar)
            self.h_snow.fill( h_snow )     ### (mutable scalar)
        else:
            self.h_snow[:] = h_snow

    def update_ice_depth(self):

        #---------------------------------  
        # Change ice depth due to melting 
        #---------------------------------
        # This assumes that update_iwe() is called
        # before update_ice_depth().
        #-------------------------------------------
        h_ice = self.h_iwe * self.wi_density_ratio

        #---------------------------------             
        # Save updated ice depth in self
        #---------------------------------
        if (np.ndim( self.h_ice ) == 0):
            h_ice = np.float64( h_ice )
            self.h_ice.fill( h_ice )
        else:
            self.h_ice[:] = h_ice
        
    #   update_ice_depth() 
    #------------------------------------------------------------------- 
    def update_total_snowpack_water_volume(self):

        #--------------------------------------------------------   
        # Note:  Compute the total volume of water stored
        #        in the snowpack for all grid cells in the DEM.
        #        Use this in the final mass balance reporting.
        #        (2023-08-31)
        #--------------------------------------------------------        
        # Note:  This is called from initialize() & finalize().
        #--------------------------------------------------------
  
        #----------------------------------------------------
        # Update total volume of liquid water stored in the
        # current snowpack, sum over all grid cells but no
        # integral over time.  (2023-08-31)
        #----------------------------------------------------   
        volume = np.float64(self.h_swe * self.da)  # [m^3]
        if (np.size(volume) == 1):
            vol_swe = (volume * self.rti.n_pixels)
        else:
            ## volume[ self.edge_IDs ] = 0.0  # (not needed)
            vol_swe = np.sum(volume)

        self.vol_swe.fill( vol_swe )
    
    #   update_total_snowpack_water_volume() 
    #-------------------------------------------------------------------   
    def update_total_ice_water_volume(self):

        #--------------------------------------------------------   
        # Note:  Compute the total volume of water stored
        #        in the ice for all grid cells in the DEM.
        #        Use this in the final mass balance reporting.
        #        (2023-08-31)
        #--------------------------------------------------------        
        # Note:  This is called from initialize() & finalize().
        #--------------------------------------------------------
  
        #----------------------------------------------------
        # Update total volume of liquid water stored in the
        # current ice, sum over all grid cells but no
        # integral over time.  (2023-08-31)
        #----------------------------------------------------   
        volume = np.float64(self.h_iwe * self.da)  # [m^3]
        if (np.size(volume) == 1):
            vol_iwe = (volume * self.rti.n_pixels)
        else:
            ## volume[ self.edge_IDs ] = 0.0  # (not needed)
            vol_iwe = np.sum(volume)

        self.vol_iwe.fill( vol_iwe )
    
    #   update_total_snowpack_water_volume() 
    #-------------------------------------------------------------------   
    def open_input_files(self):

        #------------------------------------------------------
        # Each component that inherits from glacier_base.py must
        # implement its own versions of these.
        #------------------------------------------------------
        print('ERROR: open_input_files() for Glacier component')
        print('       has not been implemented.')

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        print('ERROR: read_input_files() for Glacier component')
        print('       has not been implemented.')
        
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):

        print('ERROR: close_input_files() for Glacier component')
        print('       has not been implemented.')

    #   close_input_files()
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.smr_gs_file = (self.out_directory + self.smr_gs_file)
        self.hs_gs_file = (self.out_directory + self.hs_gs_file)
        self.sw_gs_file = (self.out_directory + self.sw_gs_file)
        self.cc_gs_file = (self.out_directory + self.cc_gs_file)
        self.imr_gs_file = (self.out_directory + self.imr_gs_file)
        self.hi_gs_file = (self.out_directory + self.hi_gs_file)
        self.cci_gs_file = (self.out_directory + self.cci_gs_file)
        self.iw_gs_file = (self.out_directory + self.iw_gs_file)
        self.tmr_gs_file = (self.out_directory + self.tmr_gs_file)
        #---------------------------------------------------------
        self.smr_ts_file = (self.out_directory + self.smr_ts_file)
        self.hs_ts_file = (self.out_directory + self.hs_ts_file)
        self.sw_ts_file = (self.out_directory + self.sw_ts_file)
        self.cc_ts_file = (self.out_directory + self.cc_ts_file)
        self.imr_ts_file = (self.out_directory + self.imr_ts_file)
        self.hi_ts_file = (self.out_directory + self.hi_ts_file)
        self.cci_ts_file =(self.out_directory + self.cci_ts_file)
        self.iw_ts_file =(self.out_directory + self.iw_ts_file)
        self.tmr_ts_file =(self.out_directory + self.tmr_ts_file)

    #   update_outfile_names()
    #-------------------------------------------------------------------  
    def disable_all_output(self):
    
        self.SAVE_SMR_GRIDS  = False
        self.SAVE_HS_GRIDS  = False
        self.SAVE_SW_GRIDS  = False
        self.SAVE_CC_GRIDS  = False
        self.SAVE_IMR_GRIDS = False
        self.SAVE_HI_GRIDS = False
        self.SAVE_CCI_GRIDS = False
        self.SAVE_TMR_GRIDS = False
        #-------------------------------
        self.SAVE_SMR_PIXELS = False
        self.SAVE_HS_PIXELS  = False
        self.SAVE_SW_PIXELS  = False
        self.SAVE_CC_PIXELS  = False
        self.SAVE_IMR_PIXELS = False
        self.SAVE_HI_PIXELS = False
        self.SAVE_CCI_PIXELS = False
        self.SAVE_TMR_PIXELS = False
        
    #   disable_all_output()  
    #-------------------------------------------------------------------  
    def open_output_files(self):

        model_output.check_netcdf( SILENT=self.SILENT )
        self.update_outfile_names()
        
        #----------------------------------
        # Open files to write grid stacks
        #----------------------------------
        if (self.SAVE_SMR_GRIDS):
            model_output.open_new_gs_file( self, self.smr_gs_file, self.rti,
                                           ## var_name='MR',
                                           var_name='smr',
                                           long_name='snow_meltrate',
                                           units_name='m/s')
            
        if (self.SAVE_HS_GRIDS):
            model_output.open_new_gs_file( self, self.hs_gs_file, self.rti,
                                           ## var_name='h_snow',
                                           var_name='hs',
                                           long_name='snow_depth',
                                           units_name='m')
            
        if (self.SAVE_SW_GRIDS):
            model_output.open_new_gs_file( self, self.sw_gs_file, self.rti,
                                           ## var_name='SWE',
                                           var_name='sw',
                                           long_name='snow_water_equivalent',
                                           units_name='m')
            
        if (self.SAVE_CC_GRIDS):
            model_output.open_new_gs_file( self, self.cc_gs_file, self.rti,
                                           ## var_name='SCC',
                                           var_name='cc',
                                           long_name='snow_cold_content',
                                           units_name='J/m^2')
            
        if (self.SAVE_IMR_GRIDS):
            model_output.open_new_gs_file( self, self.imr_gs_file, self.rti,
                                          var_name = 'imr',
                                          long_name = 'ice_meltrate',
                                          units_name = 'm/s')
        if (self.SAVE_HI_GRIDS):
            model_output.open_new_gs_file( self, self.hi_gs_file, self.rti,
                                          var_name='hi',
                                          long_name='ice_depth',
                                          units_name='m')
            
        if (self.SAVE_IW_GRIDS):
            model_output.open_new_gs_file( self, self.iw_gs_file, self.rti,
                                           var_name='iw',
                                           long_name='ice_water_equivalent',
                                           units_name='m')
        if (self.SAVE_CCI_GRIDS):
            model_output.open_new_gs_file( self, self.cci_gs_file, self.rti,
                                           var_name='cci',
                                           long_name='ice_cold_content',
                                           units_name='J/m^2')
            
        if (self.SAVE_TMR_GRIDS):
            model_output.open_new_gs_file( self, self.tmr_gs_file, self.rti,
                                           ## var_name='MR',
                                           var_name='tmr',
                                           long_name='combined_snow_ice_meltrate',
                                           units_name='m/s')
        #---------------------------------------
        # Open text files to write time series
        #---------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_SMR_PIXELS):
            model_output.open_new_ts_file( self, self.smr_ts_file, IDs,
                                           ## var_name='MR',
                                           var_name='smr',
                                           long_name='snow_meltrate',
                                           units_name='m/s')

        if (self.SAVE_HS_PIXELS):
            model_output.open_new_ts_file( self, self.hs_ts_file, IDs,
                                           ## var_name='h_snow',
                                           var_name='hs',
                                           long_name='snow_depth',
                                           units_name='m')

        if (self.SAVE_SW_PIXELS):
            model_output.open_new_ts_file( self, self.sw_ts_file, IDs,
                                           ## var_name='SWE',
                                           var_name='sw',
                                           long_name='snow_water_equivalent',
                                           units_name='m')
            
        if (self.SAVE_CC_PIXELS):
            model_output.open_new_ts_file( self, self.cc_ts_file, IDs,
                                           ## var_name='SCC',
                                           var_name='cc',
                                           long_name='snow_cold_content',
                                           units_name='J/m^2')
        
        if (self.SAVE_IMR_PIXELS):
            model_output.open_new_ts_file( self, self.imr_ts_file, IDs,
                                          var_name='imr',
                                          long_name='ice_meltrate',
                                          units_name='m/s')
        
        if (self.SAVE_HI_PIXELS):
            model_output.open_new_ts_file( self, self.hi_ts_file, IDs,
                                          var_name='hi',
                                          long_name='ice_depth',
                                          units_name='m')

        if (self.SAVE_IW_PIXELS):
            model_output.open_new_ts_file( self, self.iw_ts_file, IDs,
                                           var_name='iw',
                                           long_name='ice_water_equivalent',
                                           units_name='m')
        if (self.SAVE_CCI_PIXELS):
            model_output.open_new_ts_file( self, self.cci_ts_file, IDs,
                                           var_name='cci',
                                           long_name='ice_cold_content',
                                           units_name='J/m^2')
            
        if (self.SAVE_TMR_PIXELS):
            model_output.open_new_ts_file( self, self.tmr_ts_file, IDs,
                                           ## var_name='MR',
                                           var_name='tmr',
                                           long_name='combined_snow_ice_meltrate',
                                           units_name='m/s')
            
    #   open_output_files()
    #-------------------------------------------------------------------
    def write_output_files(self, time_seconds=None):

        #-----------------------------------------
        # Allows time to be passed from a caller
        #-----------------------------------------
        if (time_seconds is None):
            time_seconds = self.time_sec
        model_time = int(time_seconds)
        
        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
        if (model_time % int(self.save_grid_dt) == 0):
            self.save_grids()
        if (model_time % int(self.save_pixels_dt) == 0):
            self.save_pixel_values()

        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
##        if ((self.time_index % self.grid_save_step) == 0):
##             self.save_grids()
##        if ((self.time_index % self.pixel_save_step) == 0):
##             self.save_pixel_values()
        
    #   write_output_files()
    #-------------------------------------------------------------------
    def close_output_files(self):
    
        if (self.SAVE_SMR_GRIDS):model_output.close_gs_file( self, 'smr')   
        if (self.SAVE_HS_GRIDS): model_output.close_gs_file( self, 'hs')   
        if (self.SAVE_SW_GRIDS): model_output.close_gs_file( self, 'sw')   
        if (self.SAVE_CC_GRIDS): model_output.close_gs_file( self, 'cc')
        if (self.SAVE_IMR_GRIDS):model_output.close_gs_file( self, 'imr')
        if (self.SAVE_HI_GRIDS): model_output.close_gs_file( self, 'hi')
        if (self.SAVE_CCI_GRIDS): model_output.close_gs_file( self, 'cci')
        if (self.SAVE_TMR_GRIDS):model_output.close_gs_file( self, 'tmr')   
        #-----------------------------------------------------------------        
        if (self.SAVE_SMR_PIXELS):model_output.close_ts_file( self, 'smr')  
        if (self.SAVE_HS_PIXELS): model_output.close_ts_file( self, 'hs')   
        if (self.SAVE_SW_PIXELS): model_output.close_ts_file( self, 'sw')   
        if (self.SAVE_CC_PIXELS): model_output.close_ts_file( self, 'cc')
        if (self.SAVE_IMR_PIXELS):model_output.close_ts_file( self, 'imr')
        if (self.SAVE_HI_PIXELS): model_output.close_ts_file( self, 'hi')
        if (self.SAVE_CCI_PIXELS): model_output.close_ts_file( self, 'cci')
        if (self.SAVE_TMR_PIXELS):model_output.close_ts_file( self, 'tmr')  
        
    #-------------------------------------------------------------------  
    def save_grids(self):
     
        if (self.SAVE_SMR_GRIDS):
            model_output.add_grid( self, self.SM, 'smr', self.time_min )
            
        if (self.SAVE_HS_GRIDS):
            model_output.add_grid( self, self.h_snow, 'hs', self.time_min )
            
        if (self.SAVE_SW_GRIDS):
            model_output.add_grid( self, self.h_swe, 'sw', self.time_min )

        if (self.SAVE_CC_GRIDS):
            model_output.add_grid( self, self.Eccs, 'cc', self.time_min )

        if (self.SAVE_IMR_GRIDS):
            model_output.add_grid(self, self.IM, 'imr', self.time_min )

        if (self.SAVE_HI_GRIDS):
            model_output.add_grid(self, self.h_ice, 'hi', self.time_min )

        if (self.SAVE_IW_GRIDS):
            model_output.add_grid( self, self.h_iwe, 'iw', self.time_min )

        if (self.SAVE_CCI_GRIDS):
            model_output.add_grid( self, self.Ecci, 'cci', self.time_min )

        if (self.SAVE_TMR_GRIDS):
            model_output.add_grid( self, self.M_total, 'tmr', self.time_min )

    #   save_grids()     
    #-------------------------------------------------------------------  
    def save_pixel_values(self):

        IDs  = self.outlet_IDs
        time = self.time_min   ###
        
        if (self.SAVE_SMR_PIXELS):
            model_output.add_values_at_IDs( self, time, self.SM, 'smr', IDs )
            
        if (self.SAVE_HS_PIXELS):
            model_output.add_values_at_IDs( self, time, self.h_snow, 'hs', IDs )
            
        if (self.SAVE_SW_PIXELS):
            model_output.add_values_at_IDs( self, time, self.h_swe, 'sw', IDs )
            
        if (self.SAVE_CC_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Eccs, 'cc', IDs )

        if (self.SAVE_IMR_PIXELS):
            model_output.add_values_at_IDs( self, time, self.IM, 'imr', IDs )
        
        if (self.SAVE_HI_PIXELS):
            model_output.add_values_at_IDs( self, time, self.h_ice, 'hi', IDs )

        if (self.SAVE_IW_PIXELS):
            model_output.add_values_at_IDs( self, time, self.h_iwe, 'iw', IDs )

        if (self.SAVE_CCI_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Ecci, 'cci', IDs )

        if (self.SAVE_TMR_PIXELS):
            model_output.add_values_at_IDs( self, time, self.M_total, 'tmr', IDs )
            
    #   save_pixel_values()
    #-------------------------------------------------------------------