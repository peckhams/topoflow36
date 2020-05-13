"""
This file defines a "base class" for infiltration components as well
as functions used by most or all infiltration methods.  That is, all
infiltration components inherit methods from this class.  The methods
of this class should be over-ridden as necessary (especially the
"update_infil_rate() method) for different methods of modeling
infiltration.  This class, in turn, inherits from the "BMI base class"
in BMI_base.py.

See infil_beven.py, infil_green_ampt.py, infil_smith_parlange.py and
infil_richards_1D.py.
"""
#-----------------------------------------------------------------------
#
#  Copyright (c) 2001-2020, Scott D. Peckham
#
#  May 2020.  Added update_total_storage(), vol_soil.
#  Apr 2020.  Added set_new_defaults(), disable_all_output().
#  Jan 2013.  Removed CCA "get_port_data" calls, etc.
#  May 2010.  Changes to initialize() and read_cfg_file()
#  Aug 2009   Updates
#  May 2009.  Updates
#  Jan 2009.  Converted from IDL
#
#-----------------------------------------------------------------------
#  Notes:  This file defines a "base class" for infiltration
#          components as well as functions used by most or
#          all infiltration methods.  The methods of this class
#          (especially "update_infil_rate") should be over-ridden as
#          necessary for different methods of modeling infiltration.
#          See infil_green_ampt.py, infil_smith_parlange.py and
#          infil_richards_1D.py.
#-----------------------------------------------------------------------
#
#  class infil_component
#
#      set_constants()
#      set_new_defaults()    # (2020-04-29)
#      initialize()
#      update()
#      update_nondrivers()          ####### (OBSOLETE ??)
#      finalize()
#      initialize_layer_vars()      # (5/11/10)
#      build_layered_var()
#      build_layer_z_vector()       # Moved here, 11/13/16.
#      set_computed_input_vars()
#      ---------------------------
#      check_input_types()
#      initialize_computed_vars()    #####
#      ---------------------------
#      update_surface_influx()       #####
#      update_infil_rate()
#      adjust_infil_rate()           #####
#      update_v0_integral()
#      update_Rg()
#      update_Rg_integral()
#      update_I()
#      update_q0()
#      update_total_storage()
#      check_infiltration()
#      -------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      ------------------------
#      update_outfile_names()
#      disable_all_output()    # (2020-04-29)
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#      save_profiles()
#      save_cubes()
#      ------------------------
#      save_profiles_old()
#      save_cubes_old()

#-----------------------------------------------------------------------

import numpy as np
import os, sys

from topoflow.utils import BMI_base
from topoflow.utils import model_input
from topoflow.utils import model_output

#-----------------------------------------------------
# Moved these imports to: "embed_child_components()"
# because first one tries to "import infil_base" and
# this leads to a problem.
#-----------------------------------------------------
## from topoflow.components import channels_kinematic_wave
## from topoflow.components import snow_degree_day
## from topoflow.components import evap_priestley_taylor
## from topoflow.components import satzone_darcy_layers
## from topoflow.components import met_base

#-----------------------------------------------------------------------
class infil_component( BMI_base.BMI_component):

    #---------------------------------------------------------
    # Notes: Default settings are average for 'Loam', as
    #        returned by the Get_Soil_Params routine.
    #        Vars in 2nd line are set fairly arbitrary.

    #        eta = (2 + (3*lambda)) and needs to be
    #        updated whenever lambda is updated, as in
    #        Update_Infil_Vars.  We want to avoid
    #        recomputing it everywhere it is used in
    #        order to save time.

    #        The vars computed by Richards' method are
    #        set to scalar 0d below, but they are set to
    #        arrays by Route_Flow according to user choices.
    #---------------------------------------------------------
    # NB!    soil types are only used to pass defaults to
    #        the droplists in the GUI.  "soil_type" field
    #        is only used by the non-Richards routines.

    # NB!    dz is a pointer, but dz1, dz2 & dz3 are scalars.
    #        For multiple soil layers, we build a 1D array
    #        for dz from the others.
    #---------------------------------------------------------

    #-------------------------------------------------------------------
    def set_constants(self):       

        #------------------------
        # Define some constants
        #--------------------------------------------
        # See Figure 6-13, p. 237 in Dingman (2002)
        #----------------------------------------------------
        # Psi_field_capacity is often used for psi_init.
        # See initialize_theta_i() in infil_richards_1D.py.
        #----------------------------------------------------
        self.psi_oven_dry    = np.float64(-1e8)      # [m], oven dry
        self.psi_air_dry     = np.float64(-1e4)      # [m], air dry
        self.psi_min         = np.float64(-1e4)      # [m], air dry
        self.psi_hygro       = np.float64(-310)      # [m], hygroscopic
        self.psi_wilt        = np.float64(-150)      # [m], wilting pt.
        self.psi_field       = np.float64(-3.4)      # [m], field cap.
        #---------------------------------------------------------------
        self.psi_oven_dry_cm = np.float64(-1e10)     # [cm], oven dry
        self.psi_air_dry_cm  = np.float64(-1e6)      # [cm], air dry
        self.psi_min_cm      = np.float64(-1e6)      # [cm], air dry
        self.psi_hygro_cm    = np.float64(-31000)    # [cm], hygroscopic
        self.psi_wilt_cm     = np.float64(-15000)    # [cm], wilting pt.
        self.psi_field_cm    = np.float64(-340)      # [cm], field cap.
        
        #-------------------------------------
        # Why isn't this used anywhere yet ?
        #-------------------------------------
        self.g = np.float64(9.81)   # (gravitation const.)
        
    #   set_constants()
    #-------------------------------------------------------------------
    def set_new_defaults(self):
 
        #------------------------------------------------------
        # Notes:  Set defaults for flags or vars that may not
        #         be found in older CFG files. (2020-04-29)
        #------------------------------------------------------
        # This is much better than putting them in the
        # save_computed_input_vars() function.
        #------------------------------------------------------
        self.SINGLE_PROFILE  = True  # (must be set)
        self.CHECK_STABILITY = False

        #----------------------------------------------
        # These may not be defined in older CFG files
        # and are needed by open_output_files(), etc.
        #----------------------------------------------
        self.SAVE_Q0_GRIDS  = False
        self.SAVE_ZW_GRIDS  = False
        #-----------------------------
        self.SAVE_Q0_PIXELS = False
        self.SAVE_ZW_PIXELS = False
        #-----------------------------
        self.SAVE_Q_CUBES   = False
        self.SAVE_P_CUBES   = False
        self.SAVE_K_CUBES   = False
        self.SAVE_V_CUBES   = False

        #----------------------------------------
        # These may no longer be needed;  check.
        #----------------------------------------
        self.q_cs_file = ''
        self.p_cs_file = ''
        self.K_cs_file = ''
        self.v_cs_file = ''
        self.save_cube_dt = np.float64( 60 )
        
    #   set_new_defaults()         
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        ## self.DEBUG = True
                
        #---------------------------------------------------------
        # Notes:  Need to make sure than h_swe matches h_snow ?
        #         User may have entered incompatible valueself.
        #---------------------------------------------------------
        # (3/14/07) If the Energy Balance method is used for ET,
        # then we must initialize and track snow depth even if
        # there is no snowmelt method because the snow depth
        # affects the ET rate.  Otherwise, return to caller.
        #---------------------------------------------------------
        if not(SILENT):
            print(' ')
            print('Infiltration component: Initializing...')
            
        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file

        #-------------------------------------------------
        # Richards' method is special, so check for it
        # But also set now in set_computed_input_vars().
        #-------------------------------------------------
        cfg_extension = self.get_attribute( 'cfg_extension' )
        self.RICHARDS = ('richards' in cfg_extension.lower())
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-------------------------------------------------------
        # NOTE!  In BMI_base.py: function calls are nested as:
        #        initialize_config_vars()
        #            read_config_file()
        #                # read n_layers...
        #                initialize_layer_vars()
        #            set_computed_input_vars()
        #---------------------------------------------------------- 
        # HOWEVER, read_input_files() is called after
        # initialize_config_vars().
        # Should we move set_computed_input_vars() to AFTER
        # read_input_files() & before initialize_computed_vars()?
        #----------------------------------------------------------              
        self.set_constants()
        self.set_new_defaults()
        # Note: initialize_config_vars() calls initialize_layer_vars()
        self.initialize_config_vars() 
        ### self.set_missing_cfg_options() #### MAYBE LATER
        ### self.read_grid_info()  # NOW IN initialize_config_vars()
        self.initialize_basin_vars()  # (5/14/10)
        self.initialize_time_vars()

        #----------------------------------
        # Has component been turned off ?
        #----------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print('Infiltration component: Disabled in CFG file.')
            self.disable_all_output()
            #-------------------------------------------------
            # IN = infiltration rate at land surface
            # v0 = Darcy velocity at bottom of top cell
            # Rg = vertical flow rate just above water table
            #-------------------------------------------------
            self.IN     = self.initialize_scalar( 0, dtype='float64' )
            self.v0     = self.initialize_scalar( 0, dtype='float64' )
            self.Rg     = self.initialize_scalar( 0, dtype='float64' )
            self.vol_IN = self.initialize_scalar( 0, dtype='float64' )
            self.vol_v0 = self.initialize_scalar( 0, dtype='float64' )
            self.vol_Rg = self.initialize_scalar( 0, dtype='float64' )
            self.vol_soil = self.initialize_scalar( 0, dtype='float64' )
            self.DONE   = True
            self.status = 'initialized'  # (OpenMI 2.0 convention)
            return
   
        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        self.open_input_files()
        print('INFIL calling read_input_files()...')
        self.read_input_files()

        #----------------------------------------------
        # Must come before initialize_computed_vars()
        # because it uses ALL_SCALARS.
        #----------------------------------------------
        self.check_input_types()

        #--------------------------------------------
        # Set any input variables that are computed
        #--------------------------------------------------
        # NOTE:  Must be called AFTER read_input_files().
        #--------------------------------------------------
        print('INFIL calling set_computed_input_vars()...')
        self.set_computed_input_vars()
        
        print('INFIL calling initialize_computed_vars()...')
        self.initialize_computed_vars()
        
        # Get and save initial total storage in soil
        self.update_total_storage(REPORT=False)
        self.storage_init = self.storage.copy()

        self.open_output_files()
        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, dt=-1.0):

        #-------------------------------------------------
        # Note: self.v0 already set to 0 by initialize()
        #-------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        self.status = 'updating'  # (OpenMI 2.0 convention)

        #----------------------------------------
        # Read next met vars from input files ?
        #-----------------------------------------------------       
        # Note: read_input_files() is called by initialize()
        # and those values must be used for the "update"
        # calls before reading new ones.
        #-----------------------------------------------------
        if (self.time_index > 0):
            self.read_input_files()
                                      
        #-------------------------
        # Update computed values 
        #-------------------------
        self.update_surface_influx()    # (P_total = P + SM)
        self.update_infil_rate()
        self.adjust_infil_rate()
        self.update_v0_integral()
        self.update_Rg()
        self.update_Rg_integral()
        self.update_I()   # (total infiltrated depth)
        self.update_q0()  # (soil moisture at surface)

        #---------------------------------
        # Check for NaNs in infiltration
        #---------------------------------    
        self.check_infiltration()

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
        self.status = 'updated'  # (OpenMI 2.0 convention)
        
    #   update()  
    #-------------------------------------------------------------------
    def finalize(self):

        if (self.comp_status.lower() == 'enabled'):
            self.update_total_storage() # (2020-05-07)
  
        self.status = 'finalizing'  # (OpenMI 2.0 convention)
        self.close_input_files()   ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'  # (OpenMI 2.0 convention)

        self.print_final_report(comp_name='Infiltration component')
    
    #   finalize()
    #-------------------------------------------------------------------
#     def initialize_layer_vars(self):
# 
#         pass
#         
#     #   initialize_layer_vars()
    #-------------------------------------------------------------------
    def build_layer_z_vector(self):

        #----------------------------------------------------------
        # Most infil components don't have nz_val & dz_val set in
        # their CFG file, so add next few lines to be safe.
        #----------------------------------------------------------
        if not(self.RICHARDS):
            #-------------------------------------
            # Need these for "build_layered_var"
            #-------------------------------------
            self.nz_val = np.zeros(self.n_layers, dtype='int32') + 1
            self.nz = np.sum( self.nz_val )
            return
#         if (self.nz_val is None):
#             self.nz_val = 1
#         if (self.dz_val is None):
#             self.dz_val = 1

        #---------------------------------------------------------
        # Compute "total nz".  Need this for "build_layered_var"
        #---------------------------------------------------------
        self.nz = np.sum( self.nz_val )

        #-------------------------------        
        # First, compute the dz vector
        #-------------------------------
        dz = np.repeat(self.dz_val[0], self.nz_val[0])  # (1D ndarray)
        for j in range(1, self.n_layers):
            layer_dz = self.dz_val[j]
            layer_nz = self.nz_val[j]
            dz_j = np.repeat(layer_dz, layer_nz)  # (1D ndarray)
            dz = np.concatenate( (dz, dz_j) )

        #------------------------------------------
        # If all dz's are equal, make it a scalar
        #------------------------------------------
        dz_min = self.dz_val.min()
        dz_max = self.dz_val.max()
        if (dz_min == dz_max):
            self.dz = self.dz_val[0]
        else:
            self.dz = dz
            
        #---------------------------------------            
        # Compute z vector as a cumulative sum
        #---------------------------------------
        self.z = np.cumsum(dz)
     
    #   build_layer_z_vector()
    #------------------------------------------------------------------- 
    def build_layered_var(self, var_list_for_layers):

        #--------------------------------------------------------------
        # Notes:  This routine examines the parameters for each soil
        #         layer as set in the CFG file.  If all layers have
        #         a scalar value for a given parameter, then a 1D
        #         array (z-profile) is constructed for that variable
        #         and used for all grid cells.
        #         However, if any layer has a 2D value, then a 3D
        #         array is constructed for the variable.
        #         Due to Python's dynamic data typing, functions that
        #         use the variable will work correctly whether it is
        #         a 1D or 3D array.

        #         Note that self.nz was previously set to the sum:
        #            long(total(self.nz_val))
        #-----------------------------------------------------------------
        # Note:   Code moved here from infil_richards_1D (11/13/16).
        #-----------------------------------------------------------------

        #-----------------------------
        # Create an array of indices
        #-----------------------------
        #i[0] = 0
        #i[1] = self.nz_val[0]
        #i[2] = self.nz_val[0] + self.nz_val[1]
        #etc.
        #----------------------------------------------
        i = np.concatenate(([np.int32(0)], np.int32(np.cumsum(self.nz_val))) )
        
        #----------------------------------------
        # Do all layers have a scalar parameter
        # value for this particular variable ??
        #----------------------------------------
        nmax = np.int16(1)
        for j in range(self.n_layers):
            nj = var_list_for_layers[j].size
            nmax = np.maximum( nmax, nj )
        ALL_SCALARS = (nmax == 1)
        
        #--------------------------------------------------
        # Build either a 1D or 3D array for this variable
        #--------------------------------------------------
        if (ALL_SCALARS):
            #------------------------------------------
            # All layers have a scalar value for this
            # variable, so build a 1D array.
            #------------------------------------------
            var = np.zeros([self.nz], dtype='float64')
            for j in range(self.n_layers):
                var[i[j]: i[j + 1]] = var_list_for_layers[j]
        else:
            #------------------------------------------------
            # One or more layers have a grid value for this
            # variable, so build a 3D array or "data cube".    
            #--------------------------------------------------------
            # Note that all nz "levels" in a *given* layer can be
            # initialized to a grid, but not yet to different grids
            #--------------------------------------------------------
            var = np.zeros([self.nz, self.ny, self.nx], dtype='float64')
            for j in range(self.n_layers):
                for k in range(i[j], i[j + 1]):
                    var[k,:,:] = var_list_for_layers[j]

        return var

    #   build_layered_var
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        # self.nz = 1  # (needed by self.save_profiles() ??)
    
        #------------------------------------------    
        # Note: Moved some things from here to:
        #       set_new_defaults() on 2020-04-29.
        #------------------------------------------
        
        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt    = np.maximum(self.save_grid_dt,    self.dt)
        self.save_pixels_dt  = np.maximum(self.save_pixels_dt,  self.dt)
        self.save_profile_dt = np.maximum(self.save_profile_dt, self.dt)
        self.save_cube_dt    = np.maximum(self.save_cube_dt,    self.dt)
        
    #   set_computed_input_vars()
    #-------------------------------------------------------------------
    def check_input_types(self):

        #------------------------------------------------------
        # Notes: Usually this will be overridden by a given
        #        method of computing IN.  But this one should
        #        work for Green-Ampt and Smith-Parlange.
        #------------------------------------------------------
        are_scalars = np.array([
                         self.is_scalar('P_rain'),
                         self.is_scalar('SM'),
                         self.is_scalar('h_table'),
                         #----------------------------
                         self.is_scalar('Ks_list[0]'),
                         self.is_scalar('Ki_list[0]'),
                         self.is_scalar('qs_list[0]'),
                         self.is_scalar('qi_list[0]'),
                         self.is_scalar('G_list[0]')  ])
#                          self.is_scalar('Ks'),
#                          self.is_scalar('Ki'),
#                          self.is_scalar('qs'),
#                          self.is_scalar('qi'),
#                          self.is_scalar('G')  ])

        self.ALL_SCALARS = np.all(are_scalars)
        
    #   check_input_types()
    #-------------------------------------------------------------------
#     def initialize_computed_vars(self):
#   
#         #-----------------------------------------------------------
#         # Note: This is implemented separately by each of the
#         #       infiltration components. 
#         #-----------------------------------------------------------
#         pass
# 
#     #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def update_surface_influx(self):

        if (self.DEBUG):
            print('Calling update_surface_influx()...')

        P_rain = self.P_rain    # (2/3/13, new framework)
        SM     = self.SM        # (2/3/13, new framework)
        ## ET  = self.ET        # (2/3/13, new framework)

        self.P_total = (P_rain + SM)
        ## self.P_total = (P_rain + SM) - ET
        
    #   update_surface_influx()
    #-------------------------------------------------------------------
    def update_infil_rate(self):

        if (self.DEBUG):
            print('Calling update_infil_rate()...')

        #---------------------------------------------------------
        # Note: Do nothing now unless this method is overridden
        #       by a particular method of computing infil rate.
        #---------------------------------------------------------
        print("Warning: update_infil_rate() method is inactive.")
        
        #------------------------------------------------------------
        # Note: P  = precipitation rate [m/s]
        #       SM = snowmelt rate [m/s]
        #       r  = (P + SM)  [m/s]
        #       ET = evapotranspiration rate [m/s]
        #       IN = infiltration rate with units of [m/s]
        #       Rg = rate at which water reaches water table [m/s]
        #            (Rg default is 0, but changed by Richards)
        #       h  = water table elevation [m]
        #       z  = land surface elevation [m]
        #       I  = total infiltrated depth [m]
        #       n  = time step (for computing t_start & t_end)
        #------------------------------------------------------------
##        r = self.P_total   # (P_total, not R=runoff)
##         
##        if (self.method == 0):    
##            self.v0 = np.float64(0)
##            ## self.r_last = ???
##            return
##        elif (self.method == 1):    
##            self.v0 = r
##            #--------------------------------------------------------
##            # These next two are not totally correct but are stable
##            # and give results similar to the correct method
##            #--------------------------------------------------------
##        elif (self.method == 2):
##            self.v0 = Green_Ampt_Infil_Rate_v1(self, r)
##        elif (self.method == 3):    
##            self.v0 = Smith_Parlange_Infil_Rate_v1(self, r)
##            #-------------------------------------------------
##            # These next two use the correct approach but so
##            # far have convergence and "jump" issues
##            #------------------------------------------------------------
##            #** 2 : IN = Green_Ampt_Infil_Rate_v2(self, r, r_last, n)
##            #** 3 : IN = Smith_Parlange_Infil_Rate_v2(self, r, r_last, n)
##            #------------------------------------------------------------
##        elif (self.method == 4):
##            P  = self.mp.P_rain
##            SM = self.sp.SM
##            ET = self.ep.ET
##            self.v0 = Richards_Infil_Rate(self, P, SM, ET, self.Rg)
##            #########################################################
##            ##  Richards_Infil_Rate should also return and save Rg
##            #########################################################
##            #** 5 : IN = Beven_Exp_K_Infil_Rate_v1(self, r)    ;(no GUI yet)
##        else:    
##            self.v0 = np.float64(0)
##        
##        #---------------------------
##        # Print min and max values
##        #---------------------------
##        #nI = np.size(Iself.N)
##        #if (nI == 1):
##        #    print 'IN =', self.v0
##        #else:
##        #    imin = self.v0.min()
##        #    imax = self.v0.max()
##        #    print '(imin, imax) =', imin, imax
##        
##        #--------------------------
##        # For debugging & testing
##        #--------------------------
##        #print 'min(IN), max(IN) = ', self.v0.min(), self.v0.max()
##        #print 'self.dt =', self.dt

    #   update_infil_rate()
    #-------------------------------------------------------------
    def adjust_infil_rate(self):

        if (self.DEBUG):
            print('Calling adjust_infil_rate()...')
    
        #-------------------------------------
        # Is P_total less than Ks anywhere ?
        # If so, set IN = P_total there.
        #-------------------------------------
        CHECK_LOW_RAIN = not(self.RICHARDS)
        if (CHECK_LOW_RAIN):
            self.check_low_rainrate()

        #####################################################
        # (10/8/10) The rest of this routine doesn't work
        # if IN is a scalar.  Need to look at this more.
        #####################################################
        ## if (self.SINGLE_PROFILE):
        if (self.IN.size == 1):
            return
        
        #-------------------------------------------
        # Get water table and land surface from gp
        #-------------------------------------------
        ## H_IS_GRID = self.gp.is_grid('h_table')
        ## Z_IS_GRID = self.gp.is_grid('elev')
        
        h = self.h_table  # (2/3/13, new framework)
        z = self.elev     # (2/3/13, new framework)
        
        ##### if (h or z is undefined): return
        
        #----------------------------------------
        # Can't infiltrate where water table
        # is above ground, i.e. where (h ge z)
        # If z & h are given, IN is a grid too.            
        #--------------------------------------------------
        # Note: h  = water table elevation [m]
        #       z  = land surface elevation [m]
        #
        #       Currently, h and z are always grids,
        #           so (nh > 1) and (nz > 1)
        #--------------------------------------------------
        w  = np.where( h >= z )  # (changed on 8/19/09)
        ### w  = np.where( h == z )
        nw = np.size(w[0])
        if (nw != 0):
            self.IN[w] = np.float64(0)  
            ## self.v0[w] = np.float64(0)
            
        ##########################################
        #  ADD SIMILAR THING FOR GC2D
        ##########################################
        ##########################################
      
    #   adjust_infil_rate()
    #-------------------------------------------------------------------
    def update_v0_integral(self):

        if (self.DEBUG):
            print('Calling update_v0_integral()...')
            
        #------------------------------------------------
        # Update mass total for v0, sum over all pixels
        #------------------------------------------------   
        volume = np.double(self.v0 * self.da * self.dt)  # [m^3]
        if (np.size( volume ) == 1):
            self.vol_v0 += (volume * self.rti.n_pixels)
        else:
            self.vol_v0 += np.sum(volume)
 
        #-----------------------------------------------------       
        # NOTE!  v0 is downward darcy velocity at bottom of
        #        topmost layer, with thickness dz.  So must
        #        add the amount of water that was added to
        #        that layer as well, in order to know how
        #        much water infiltrated.   (2020-05-09)
        #-----------------------------------------------------
        vol2 = np.double(self.IN * self.da * self.dt)  # [m^3]
        if (np.size( vol2 ) == 1):
            self.vol_IN += (vol2 * self.rti.n_pixels)
        else:
            self.vol_IN += np.sum(vol2)
        
        # Obsolete.  v0 is not same as IN.    
        #--------------------------------------------------         
        # Create a synonym (2nd reference to same values)
        #--------------------------------------------------  
        ## self.vol_IN = self.vol_v0
    
    #   update_v0_integral()
    #-------------------------------------------------------------------
    def update_Rg(self):

        if (self.DEBUG):
            print('Calling update_Rg()...')
            
        #------------------------------------------------
        # This works for Green-Ampt and Smith_Parlange,
        # but should be overridden for Richards 1D.
        #------------------------------------------------
        # Set groundwater recharge rate to IN ?
        # Save last value of r for next time.
        #----------------------------------------   
        self.Rg = self.IN.copy()
        P_rain  = self.P_rain   # (2/3/13, new framework)
        SM      = self.SM       # (2/3/13, new framework)
        #---------------------
        self.r_last = (P_rain + SM)

    #   update_Rg()
    #-------------------------------------------------------------------
    def update_Rg_integral(self):

        if (self.DEBUG):
            print('Calling update_Rg_integral()...')
            
        #------------------------------------------------
        # Update mass total for Rg, sum over all pixels
        #------------------------------------------------   
        volume = np.double(self.Rg * self.da * self.dt)  # [m^3]
        if (np.size( volume ) == 1):
            self.vol_Rg += (volume * self.rti.n_pixels)
        else:
            self.vol_Rg += np.sum(volume)
            
    #   update_Rg_integral()
    #-------------------------------------------------------------------
    def update_I(self):

        if (self.DEBUG):
            print('Calling update_I()...')
            
        #---------------------------------------
        # Update the total infiltrated depth
        #---------------------------------------
        # Do this for all methods ?  I is not
        # used for Richards 1D method.
        #
        # This becomes a grid if IN is a grid.
        # v0 = Darcy velocity at bottom of top cell
        # IN = Infiltration rate at top of top cell
        #--------------------------------------------
        self.I += (self.IN * self.dt)  # [meters]

        #--------------
        # For testing
        #--------------
        #if (np.size(self.I) == 1):
        #    print '    Tot. infil. depth =' + str(self.I)
        
        #------------------------------------------
        # Reset the total infiltrated depth after
        # each "event";  not sure how to do this.
        #------------------------------------------
        #if (????): self.I = np.minimum(self.I, 0.0)

    #   update_I()
    #-------------------------------------------------------------------
    def update_q0(self):

        if (self.DEBUG):
            print('Calling update_q0()...')
            
        #-----------------------------------------------
        # Note: self.q0 = np.float64(0) in __init__().
        # Most infil methods don't compute q0.
        # This method is over-ridden for Richards 1D
        #-----------------------------------------------
        pass

    #   update_q0()
    #-------------------------------------------------------------------
    def update_total_storage(self, REPORT=True):

        #------------------------------------------------------
        # Note:  Update total volume of water stored in soil.
        #        Call this from finalize() for mass balance.
        #        Vw = volume of water, Vc = volume of cell
        #        theta = (Vw/Vc) = q = soil water content
        #------------------------------------------------------
        # https://docs.scipy.org/doc/numpy/reference/
        #       routines.linalg.html
        #------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_total_storage()...')
            
        if (self.SINGLE_PROFILE):
            #-------------------------------
            # theta = q is a 1D array in z
            #-------------------------------
            Vc = self.da * self.dz  # 0D or 1D array
            Vw = self.q * Vc        # 1D array
            Vw_profile = Vw.sum()
            self.storage = (Vw_profile * self.rti.n_pixels)
        else:
           #--------------------------
           # theta = q is a 3D array
           #--------------------------
           # da can be 0D or 2D; dz can be 0D or 1D
           da_size = self.da.size   # (1 if 0D, >1 if 2D)
           dz_size = self.dz.size   # (1 if 0D, >1 if 1D)
                  
           if (da_size == 1):
               #-------------------------
               # Vc is a 0D or 1D array
               #-------------------------
               Vc = self.da * self.dz
           else:
               #-------------------
               # Vc is a 3D array
               #-------------------
               if (dz_size == 1):
                   # Convert dz to a 1D array
                   dz = self.dz * np.ones( self.nz )
               # Yes, this really works.
               Vc = self.da[None,:,:] * self.dz[:,None,None]
               #----------------------------------------------
               # This also works.  Need dz as first arg.
               # Vc = np.outer(self.dz, self.da)
               # Vc = Vc.reshape(self.nz, self.ny, self.nx)

           if (Vc.ndim == 1):
               Vw = self.q * Vc[:,None,None]
           else:
               Vw = self.q * Vc   # (3D * 1D, or 3D * 3D, elementwise)
    
           self.storage = Vw.sum()

        #-------------------------------------------------
        # Save water added to vol_soil, for final report
        #-------------------------------------------------
        if not(hasattr(self, 'storage_init')):
            self.storage_init = np.float64(0)
        vol_added = (self.storage - self.storage_init)
        self.vol_soil.fill( vol_added )
        
        if (REPORT):
            print()
            print('Total volume of water in soil at start =')
            print( self.storage_init, '[m3]')
            print('Total volume of water in soil at end =')
            print( self.storage, '[m3]')
            print('Difference in volume =')
            print( vol_added, '[m3]')
            print()

    #   update_total_storage()
    #-------------------------------------------------------------------
    def check_infiltration(self):

        if (self.DEBUG):
            print('Calling check_infiltration()...')
            
        #--------------------------------------
        # Check for NaNs in infiltration rate
        #--------------------------------------
        # NB!  Don't set DONE = False, it may
        # already have been set to True
        #--------------------------------------
        if (np.size( self.IN ) == 1):
            OK = np.isfinite( self.IN )
            nbad = 1
        else:
            # Should be faster than using WHERE
            wbad = np.logical_not( np.isfinite( self.IN ) )
            nbad = wbad.sum()
#             wbad = np.where( np.logical_not(np.isfinite( self.IN )) )
#             nbad = np.size( wbad[0] )
            OK = (nbad == 0)
        if (OK):
            return
        
        #------------------------------------------
        # Issue warning message and abort the run
        #------------------------------------------
        msg = np.array(['ERROR:  Aborting model run.', \
                        '        NaNs found in infiltration rates.', \
                        '        Number of NaN values = ' + str(nbad) ])

        print('##############################################')
        for line in msg:
            print(line)
        print('##############################################')
        print(' ')
        
        self.status = 'failed'

        #--------------------------------------------------------        
        # This won't stop the run unless infiltration component
        # is the driver component.  So use sys.exit().
        #--------------------------------------------------------
        self.DONE = True
        sys.exit()
        
    #   check_infiltration
    #-----------------------------------------------------------------------
    def check_low_rainrate(self):

        #------------------------------------------------------------
        # Notes:  If (P_total < Ks), then we need to set the infil
        #         rate to P_total.  P_total = (P + SM).
        #
        #         This needs to be called by Green-Ampt and Smith-
        #         Parlange methods for computing IN;  perhaps by
        #         any method based on total infiltrated depth, I.
        #         This isn't needed for Richards' equation method.
        #------------------------------------------------------------

        #--------------------------------------
        # Is P_total less than Ks anywhere ?
        # If so, set IN = P_total there.
        #--------------------------------------
        nPt = np.size( self.P_total )
        nK  = np.size( self.Ks[0] )
        if ((nPt == 1) and (nK == 1)):    
            #----------------------------------
            # P_total and Ks are both scalars
            #----------------------------------
            if (self.P_total < self.Ks[0]):
                self.IN = self.P_total    
                ## self.v0 = self.P_total
        else:    
            #---------------------------------
            # Either P_total or Ks is a grid
            # so IN will have become a grid
            #---------------------------------
            w  = np.where( self.P_total < self.Ks[0] )
            nw = np.size( w[0] )
            
            if (nw != 0):    
                if (nPt > 1):
                    self.IN[w] = self.P_total[w]   
                    ## self.v0[w] = self.P_total[w]
                else:
                    self.IN[w] = self.P_total   
                    ## self.v0[w] = self.P_total

    #   check_low_rainrate
    #-------------------------------------------------------------------  
    def open_input_files(self):

        #-------------------------------------------------------
        # This method works for Green-Ampt and Smith-Parlange
        # but must be overridden for Richards 1D.
        #-------------------------------------------------------
        # NB! Green-Ampt and Smith-Parlange currently only
        #     support ONE layer (n_layers == 1).
        #-------------------------------------------------------
        self.Ks_unit  = []  # (empty lists to hold file objects)
        self.Ki_unit  = []
        self.qs_unit  = []
        self.qi_unit  = []
        self.G_unit   = []
        self.gam_unit = []

        for k in range(self.n_layers):
            self.Ks_file[k] = self.soil_directory + self.Ks_file[k]
            self.Ki_file[k] = self.soil_directory + self.Ki_file[k]
            self.qs_file[k] = self.soil_directory + self.qs_file[k]
            self.qi_file[k] = self.soil_directory + self.qi_file[k]
            self.G_file[k]  = self.soil_directory + self.G_file[k]
            self.gam_file[k]= self.soil_directory + self.gam_file[k]

            self.Ks_unit.append(  model_input.open_file(self.Ks_type[k],  self.Ks_file[k]) )
            self.Ki_unit.append(  model_input.open_file(self.Ki_type[k],  self.Ki_file[k]) )
            self.qs_unit.append(  model_input.open_file(self.qs_type[k],  self.qs_file[k]) )
            self.qi_unit.append(  model_input.open_file(self.qi_type[k],  self.qi_file[k]) )
            self.G_unit.append(   model_input.open_file(self.G_type[k],   self.G_file[k])  )
            self.gam_unit.append( model_input.open_file(self.gam_type[k], self.gam_file[k]) )
        
    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        if (self.DEBUG):
            print('Calling read_input_files()...')
            
        rti = self.rti

        #-------------------------------------------------------
        # All grids are assumed to have data type of Float32.
        #-------------------------------------------------------
        # This method works for Green-Ampt and Smith-Parlange
        # but must be overridden for Richards 1D.
        #-------------------------------------------------------
        # NB! Green-Ampt and Smith-Parlange currently only
        #     support ONE layer (n_layers == 1).
        #------------------------------------------------------- 
        for k in range(self.n_layers):
            Ks = model_input.read_next(self.Ks_unit[k], self.Ks_type[k], rti)
            if (Ks is not None): self.Ks[k] = Ks

            Ki = model_input.read_next(self.Ki_unit[k], self.Ki_type[k], rti)
            if (Ki is not None): self.Ki[k] = Ki

            qs = model_input.read_next(self.qs_unit[k], self.qs_type[k], rti)
            if (qs is not None): self.qs[k] = qs

            qi = model_input.read_next(self.qi_unit[k], self.qi_type[k], rti)
            if (qi is not None): self.qi[k] = qi
            
            G  = model_input.read_next(self.G_unit[k], self.G_type[k], rti)
            if (G is not None): self.G[k] = G

            gam = model_input.read_next(self.gam_unit[k], self.gam_type[k], rti)
            if (gam is not None): self.gam[k] = gam
          
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):

        #-------------------------------------------------------
        # This method works for Green-Ampt and Smith-Parlange
        # but must be overridden for Richards 1D.
        #-------------------------------------------------------
        # NB! Green-Ampt and Smith-Parlange currently only
        #     support ONE layer (n_layers == 1).
        #-------------------------------------------------------        
        for k in range(self.n_layers):
            if (self.Ks_type[k]  != 'Scalar'): self.Ks_unit[k].close()        
            if (self.Ki_type[k]  != 'Scalar'): self.Ki_unit[k].close()
            if (self.qs_type[k]  != 'Scalar'): self.qs_unit[k].close()
            if (self.qi_type[k]  != 'Scalar'): self.qi_unit[k].close()
            if (self.G_type[k]   != 'Scalar'): self.G_unit[k].close()
            if (self.gam_type[k] != 'Scalar'): self.gam_unit[k].close()
            #------------------------------------------------------------
##            if (self.Ks_file[k]  != ''): self.Ks_unit[k].close()        
##            if (self.Ki_file[k]  != ''): self.Ki_unit[k].close()
##            if (self.qs_file[k]  != ''): self.qs_unit[k].close()
##            if (self.qi_file[k]  != ''): self.qi_unit[k].close()
##            if (self.G_file[k]   != ''): self.G_unit[k].close()
##            if (self.gam_file[k] != ''): self.gam_unit[k].close()
          
    #   close_input_files()
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.v0_gs_file = (self.out_directory + self.v0_gs_file)
        self.I_gs_file  = (self.out_directory + self.I_gs_file)
        self.q0_gs_file = (self.out_directory + self.q0_gs_file)
        self.Zw_gs_file = (self.out_directory + self.Zw_gs_file)
        #-------------------------------------------------------------
        self.v0_ts_file = (self.out_directory + self.v0_ts_file)
        self.I_ts_file  = (self.out_directory + self.I_ts_file)
        self.q0_ts_file = (self.out_directory + self.q0_ts_file)
        self.Zw_ts_file = (self.out_directory + self.Zw_ts_file)
        #-----------------------------------------------------------------
        self.q_ps_file  = (self.out_directory + self.q_ps_file)
        self.p_ps_file  = (self.out_directory + self.p_ps_file)
        self.K_ps_file  = (self.out_directory + self.K_ps_file)
        self.v_ps_file  = (self.out_directory + self.v_ps_file)
        #-------------------------------------------------------------
        self.q_cs_file  = (self.out_directory + self.q_cs_file)
        self.p_cs_file  = (self.out_directory + self.p_cs_file)
        self.K_cs_file  = (self.out_directory + self.K_cs_file)
        self.v_cs_file  = (self.out_directory + self.v_cs_file)

  
##        self.v0_gs_file = (self.case_prefix + '_2D-v0.rts')
##        self.q0_gs_file = (self.case_prefix + '_2D-q0.rts')
##        self.I_gs_file  = (self.case_prefix + '_2D-I.rts')
##        self.Zw_gs_file = (self.case_prefix + '_2D-Zw.rts')
##        #---------------------------------------------------------
##        self.v0_ts_file = (self.case_prefix + '_0D-v0.txt')
##        self.q0_ts_file = (self.case_prefix + '_0D-q0.txt')
##        self.I_ts_file  = (self.case_prefix + '_0D-I.txt')
##        self.Zw_ts_file = (self.case_prefix + '_0D-Zw.txt')
##        #---------------------------------------------------------
##        self.q_cs_file = (self.case_prefix + '_3D-q.rt3')
##        self.p_cs_file = (self.case_prefix + '_3D-p.rt3')
##        self.K_cs_file = (self.case_prefix + '_3D-K.rt3')
##        self.v_cs_file = (self.case_prefix + '_3D-v.rt3')
##        #---------------------------------------------------------
##        self.q_ps_file = (self.case_prefix + '_1D-q.txt')
##        self.p_ps_file = (self.case_prefix + '_1D_p.txt')
##        self.K_ps_file = (self.case_prefix + '_1D_K.txt')
##        self.v_ps_file = (self.case_prefix + '_1D_v.txt')

    #   update_outfile_names()
    #-------------------------------------------------------------------
    def disable_all_output(self):
    
        self.SAVE_V0_GRIDS = False
        self.SAVE_Q0_GRIDS = False
        self.SAVE_I_GRIDS  = False
        self.SAVE_ZW_GRIDS = False
        #-----------------------------
        self.SAVE_V0_PIXELS = False
        self.SAVE_Q0_PIXELS = False
        self.SAVE_I_PIXELS  = False
        self.SAVE_ZW_PIXELS = False
        #-----------------------------
        self.SAVE_Q_PROFILES = False
        self.SAVE_P_PROFILES = False
        self.SAVE_K_PROFILES = False
        self.SAVE_V_PROFILES = False
        #-----------------------------
        self.SAVE_Q_CUBES = False
        self.SAVE_P_CUBES = False
        self.SAVE_K_CUBES = False
        self.SAVE_V_CUBES = False
                               
    #   disable_all_output()  
    #-------------------------------------------------------------------  
    def open_output_files(self):

        #-------------------------------------------------
        # Notes:  v0 = infiltration rate at surface
        #         q0 = soil moisture at surface
        #          I = total infiltrated depth
        #         Zw = wetting front
        #          q = soil moisture
        #          p = pressure head
        #          K = hydraulic conductivity
        #          v = vertical flow rate (see v0)
        #-------------------------------------------------
        model_output.check_netcdf()
        self.update_outfile_names()
        
        #--------------------------------------
        # Open new files to write grid stacks
        #----------------------------------------
        # NOTE:  var_name = 'IN' vs. 'v0' here.
        # Change all to 'IN' later.
        #----------------------------------------
        if (self.SAVE_V0_GRIDS):
            model_output.open_new_gs_file( self, self.v0_gs_file, self.rti,
                                           var_name='IN',
                                           long_name='infiltration_rate_at_surface',
                                           units_name='m/s')

        if (self.SAVE_I_GRIDS):
            model_output.open_new_gs_file( self, self.I_gs_file, self.rti,
                                           var_name='I',
                                           long_name='total_infiltrated_depth',
                                           units_name='m')

        if (self.SAVE_Q0_GRIDS):
            model_output.open_new_gs_file( self, self.q0_gs_file, self.rti,
                                           var_name='q0',
                                           long_name='soil_moisture_at_surface',
                                           units_name='none')
            
        if (self.SAVE_ZW_GRIDS):
            model_output.open_new_gs_file( self, self.Zw_gs_file, self.rti,
                                           var_name='Zw',
                                           long_name='depth_to_wetting_front',
                                           units_name='m')                               
                                
        #--------------------------------------
        # Open new files to write time series
        #--------------------------------------
        # NOTE:  var_name = 'IN' vs. 'v0' here.
        # Change all to 'IN' later.
        #----------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_V0_PIXELS):
            model_output.open_new_ts_file( self, self.v0_ts_file, IDs,
                                           var_name='IN',
                                           long_name='infiltration_rate_at_surface',
                                           units_name='m/s')

        if (self.SAVE_I_PIXELS):
            model_output.open_new_ts_file( self, self.I_ts_file, IDs,
                                           var_name='I',
                                           long_name='total_infiltrated_depth',
                                           units_name='m')
            

        if (self.SAVE_Q0_PIXELS):
            model_output.open_new_ts_file( self, self.q0_ts_file, IDs,
                                           var_name='q0',
                                           long_name='soil_moisture_at_surface',
                                           units_name='none')

        if (self.SAVE_ZW_PIXELS):
            model_output.open_new_ts_file( self, self.Zw_ts_file, IDs,
                                           var_name='Zw',
                                           long_name='depth_to_wetting_front',
                                           units_name='m')

        #-----------------------------------------------------
        # Remaining parts are only valid for Richards method
        #-----------------------------------------------------
        if not(self.RICHARDS):
            return
        
        #--------------------------------------------------
        # Open "profile files" to write vertical profiles
        #--------------------------------------------------
        if (self.SAVE_Q_PROFILES):
            model_output.open_new_ps_file( self, self.q_ps_file, IDs,
                                           z_values=self.z, z_units='m',
                                           var_name='q',
                                           long_name='soil_water_content',
                                           units_name='none')

        if (self.SAVE_P_PROFILES):    
            model_output.open_new_ps_file( self, self.p_ps_file, IDs,
                                           z_values=self.z, z_units='m',
                                           var_name='p',
                                           long_name='pressure_head',
                                           units_name='m')

        #################################################################
        # NOTE:  Should we convert these units from "m/s" to "mm/hr" ??
        #################################################################
        if (self.SAVE_K_PROFILES):    
            model_output.open_new_ps_file( self, self.K_ps_file, IDs,
                                           z_values=self.z, z_units='m',
                                           var_name='K',
                                           long_name='hydraulic_conductivity',
                                           units_name='m/s')

        if (self.SAVE_V_PROFILES):    
            model_output.open_new_ps_file( self, self.v_ps_file, IDs,
                                           z_values=self.z, z_units='m',
                                           var_name='v',
                                           long_name='vertical_flow_rate',
                                           units_name='m/s')

        #---------------------------------------------
        # Open "cube files" to write 3D grid "cubes"
        #---------------------------------------------
        if (self.SAVE_Q_CUBES):
            model_output.open_new_cs_file( self, self.q_cs_file, self.rti,
                                           var_name='q',
                                           long_name='soil_water_content',
                                           units_name='none')

        if (self.SAVE_P_CUBES):    
            model_output.open_new_cs_file( self, self.p_cs_file, self.rti,
                                           var_name='p',
                                           long_name='pressure_head',
                                           units_name='m')

        #################################################################
        # NOTE:  Should we convert these units from "m/s" to "mm/hr" ??
        #################################################################
        if (self.SAVE_K_CUBES):    
            model_output.open_new_cs_file( self, self.K_cs_file, self.rti,
                                           var_name='K',
                                           long_name='hydraulic_conductivity',
                                           units_name='m/s')

        if (self.SAVE_V_CUBES):    
            model_output.open_new_cs_file( self, self.v_cs_file, self.rti,
                                           var_name='v',
                                           long_name='vertical_flow_rate',
                                           units_name='m/s')

            
        #--------------------------------------------------
        # Open "profile files" to write vertical profiles
        #--------------------------------------------------
##        if (self.SAVE_Q_PROFILES):
##            self.q_profile_unit = open(self.q_ps_file, 'w')
##            write_ps_file_header(self.q_profile_unit, IDs, var_name='q')
##
##        if (self.SAVE_P_PROFILES):    
##            self.p_profile_unit = open(self.p_ps_file, 'w')
##            write_ps_file_header(self.p_profile_unit, IDs, var_name='p')
##
##        if (self.SAVE_K_PROFILES):    
##            self.K_profile_unit = open(self.K_ps_file, 'w')
##            write_ps_file_header(self.K_profile_unit, IDs, var_name='K')
##
##        if (self.SAVE_V_PROFILES):    
##            self.v_profile_unit = open(self.v_ps_file, 'w')
##            write_ps_file_header(self.v_profile_unit, IDs, var_name='v')

        #---------------------------------------
        # Open RT3 files to write grid "cubes"
        #---------------------------------------
##        if (self.SAVE_Q_STACKS):    
##            self.q_stack_unit = open(self.q_cs_file, 'wb')
##        if (self.SAVE_P_STACKS):    
##            self.p_stack_unit = open(self.p_cs_file, 'wb')
##        if (self.SAVE_K_STACKS):    
##            self.K_stack_unit = open(self.K_cs_file, 'wb')
##        if (self.SAVE_V_STACKS):    
##            self.v_stack_unit = open(self.v_cs_file, 'wb')
            
    #   open_output_files()
    #-------------------------------------------------------------------
    def write_output_files(self, time_seconds=None):

        #---------------------------------------------------------
        # Notes:  This function was written to use only model
        #         time (maybe from a caller) in seconds, and
        #         the save_grid_dt and save_pixels_dt parameters
        #         read by read_cfg_file().
        #
        #         read_cfg_file() makes sure that all of
        #         the "save_dts" are larger than or equal to the
        #         process dt.
        #---------------------------------------------------------
        if (self.DEBUG):
            print('Calling write_output_files()...')

        #-----------------------------------------
        # Allows time to be passed from a caller
        #-----------------------------------------
        if (time_seconds is None):
            time_seconds = self.time_sec
        model_time = int(time_seconds)
 
#         print('self.time = ' + str(self.time) )
#         print('self.time_sec = ' + str(self.time_sec) )
#         print('self.time_units = ' + self.time_units )        
#         print('model_time = ' + str(model_time) )
#         print('save_grid_dt = ' + str(self.save_grid_dt) )
#         print(' ')
        
        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
        if (model_time % int(self.save_grid_dt) == 0):
            ### print('Saving grids at time = ' + str(model_time) )
            self.save_grids()
            
        if (model_time % int(self.save_pixels_dt) == 0):
            self.save_pixel_values()

        if not(self.RICHARDS):
            return
        
        if (model_time % int(self.save_profile_dt) == 0):
            self.save_profiles()

        if (model_time % int(self.save_cube_dt) == 0):
            self.save_cubes()
      
    #   write_output_files()           
    #-------------------------------------------------------------------
    def close_output_files(self):
    
        if (self.SAVE_V0_GRIDS): model_output.close_gs_file( self, 'v0')   
        if (self.SAVE_I_GRIDS):  model_output.close_gs_file( self, 'I') 
        if (self.SAVE_Q0_GRIDS): model_output.close_gs_file( self, 'q0')   
        if (self.SAVE_ZW_GRIDS): model_output.close_gs_file( self, 'Zw') 

        if (self.SAVE_V0_PIXELS): model_output.close_ts_file( self, 'v0')      
        if (self.SAVE_I_PIXELS):  model_output.close_ts_file( self, 'I') 
        if (self.SAVE_Q0_PIXELS): model_output.close_ts_file( self, 'q0')   
        if (self.SAVE_ZW_PIXELS): model_output.close_ts_file( self, 'Zw') 

        if not(self.RICHARDS):
            return

        if (self.SAVE_Q_PROFILES): model_output.close_ps_file( self, 'q')    
        if (self.SAVE_P_PROFILES): model_output.close_ps_file( self, 'p')   
        if (self.SAVE_K_PROFILES): model_output.close_ps_file( self, 'K')   
        if (self.SAVE_V_PROFILES): model_output.close_ps_file( self, 'v')       

        if (self.SAVE_Q_CUBES): model_output.close_cs_file( self, 'q')    
        if (self.SAVE_P_CUBES): model_output.close_cs_file( self, 'p')    
        if (self.SAVE_K_CUBES): model_output.close_cs_file( self, 'K')     
        if (self.SAVE_V_CUBES): model_output.close_cs_file( self, 'v') 
        
    #   close_output_files()
    #-------------------------------------------------------------------  
    def save_grids(self):
        
        #-----------------------------------
        # Save grid stack to a netCDF file
        #---------------------------------------------
        # Note that add_grid() methods will convert
        # var from scalar to grid now, if necessary.
        #---------------------------------------------
        if (self.DEBUG):
            print('Calling save_grids()...')
            
        if (self.SAVE_V0_GRIDS):
            model_output.add_grid( self, self.IN, 'IN', self.time_min )
            ## model_output.add_grid( self, self.v0, 'v0', self.time_min )
            
        if (self.SAVE_I_GRIDS):
            model_output.add_grid( self, self.I, 'I', self.time_min )

        if (self.SAVE_Q0_GRIDS):
            model_output.add_grid( self, self.q0, 'q0', self.time_min )

        if (self.SAVE_ZW_GRIDS):
            model_output.add_grid( self, self.Zw, 'Zw', self.time_min )

    #   save_grids()                             
    #-------------------------------------------------------------------  
    def save_pixel_values(self):

        if (self.DEBUG):
            print('Calling save_pixel_values()...')
            
        IDs  = self.outlet_IDs
        time = self.time_min   ########
        
        #--------------------------------------------
        # Save a subsequence of IN var pixel values
        #--------------------------------------------  
        if (self.SAVE_V0_PIXELS):
            model_output.add_values_at_IDs( self, time, self.IN, 'IN', IDs )
            ## model_output.add_values_at_IDs( self, time, self.v0, 'v0', IDs )
        if (self.SAVE_I_PIXELS):
            model_output.add_values_at_IDs( self, time, self.I, 'I', IDs )
                                 
        #----------------------------------------
        # Richards' equation for infiltration ?
        #----------------------------------------
        if not(self.RICHARDS):
            return
        
        if (self.SAVE_Q0_PIXELS):
            model_output.add_values_at_IDs( self, time, self.q0, 'q0', IDs )
        if (self.SAVE_ZW_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Zw, 'Zw', IDs )
              
        #-------------------------------------
        # This should no longer be necessary
        #--------------------------------------                                
##            if (self.ALL_SCALARS):    
##                if (self.SAVE_Q0_PIXELS):
##                    write_ts_file_line(self.q0_ts_unit, time_min, self.q0, IDs)
##                if (self.SAVE_ZW_PIXELS):    
##                    write_ts_file_line(self.Zw_ts_unit, time_min, self.Zw, IDs)
##            else:    
##                if (self.SAVE_Q0_PIXELS):
##                    write_ts_file_line(self.q0_ts_unit, time_min, self.q0, IDs)
##                if (self.SAVE_ZW_PIXELS):    
##                    write_ts_file_line(self.Zw_ts_unit, time_min, self.Zw, IDs)

    #   save_pixel_values()
    #-------------------------------------------------------------------
    def save_profiles(self):

        #-----------------------------------------------------------
        # Notes:  A "profile" is a 1D array, in this case a set of
        #         values that vary with depth below the surface
        #         (z-axis) and that are obtained by "skewering" a
        #         "stack variable" (see above) at a prescribed set
        #         of pixel or grid cell IDs (outlet_IDs).
        #-----------------------------------------------------------
        if (self.DEBUG):
            print('Calling save_profiles()...')
            
        IDs  = self.outlet_IDs
        time = self.time_min
        
        #------------------------------------------
        # Save a subsequence of vertical profiles
        #------------------------------------------   
        if (self.SAVE_Q_PROFILES):
            model_output.add_profiles_at_IDs( self, self.q, 'q', IDs, time ) 
            
        if (self.SAVE_P_PROFILES):
            model_output.add_profiles_at_IDs( self, self.p, 'p', IDs, time ) 
            
        if (self.SAVE_K_PROFILES):
            model_output.add_profiles_at_IDs( self, self.K, 'K', IDs, time ) 
            
        if (self.SAVE_V_PROFILES):
            model_output.add_profiles_at_IDs( self, self.v, 'v', IDs, time ) 
            
    #   save_profiles()
    #-------------------------------------------------------------------
    def save_cubes(self):

        #---------------------------------------------------------
        # Notes:  A "cube" is a 3D array, in this case for a set
        #         of subsurface values such as K or v that vary
        #         in 3 space dimensions and with time.  This
        #         function saves a "snapshot" of such a 3D array
        #         to a "cube file" whenever it is called.
        #---------------------------------------------------------
        if (self.DEBUG):
            print('Calling save_cubes()...')
            
        time = self.time_min
        
        #------------------------------------------
        # Save a subsequence of vertical profiles
        #------------------------------------------   
        if (self.SAVE_Q_CUBES):
            model_output.add_cube( self, time, self.q, 'q' ) 
            
        if (self.SAVE_P_CUBES):
            model_output.add_cube( self, time, self.p, 'p' ) 
            
        if (self.SAVE_K_CUBES):
            model_output.add_cube( self, time, self.K, 'K' ) 
            
        if (self.SAVE_V_CUBES):
            model_output.add_cube( self, time, self.v, 'v' ) 
            
    #   save_cubes()
    #-------------------------------------------------------------------
##    def save_profiles_old(self):
##
##        #-----------------------------------------------------------
##        # Notes:  A "profile" is a 1D array, in this case a set of
##        #         values that vary with depth below the surface
##        #         (z-axis) and that are obtained by "skewering" a
##        #         "stack variable" (see above) at a prescribed set
##        #         of pixel or grid cell IDs (outlet_IDs).
##        #-----------------------------------------------------------
##        nz  = self.nz
##        IDs = self.outlet_IDs
##        
##        #----------------------------------
##        # Construct a "time stamp" string
##        #----------------------------------
##        tmstr = '***********  Time = '
##        tmstr = tmstr + ('%f8.1' % self.time_min)   ######
##        tmstr = tmstr + '  [minutes]'
##        
##        #---------------------------------------
##        # Save a subsequence of IN var profiles
##        #---------------------------------------   
##        if (self.SAVE_Q_PROFILES):    
##            Write_Profile(self.q_profile_unit, self.q, IDs, nz, tmstr)
##        if (self.SAVE_P_PROFILES):    
##            Write_Profile(self.p_profile_unit, self.p, IDs, nz, tmstr)
##        if (self.SAVE_K_PROFILES):    
##            Write_Profile(self.K_profile_unit, self.K, IDs, nz, tmstr)
##        if (self.SAVE_V_PROFILES):    
##            Write_Profile(self.v_profile_unit, self.v, IDs, nz, tmstr)
##
##    #   save_profiles_old()
    #-------------------------------------------------------------------  
##    def save_cubes_old(self):
##
##        #---------------------------------------------------------
##        # Notes:  A "stack" is a 3D array, in this case for a set
##        #         of subsurface values such as K or v that vary
##        #         in 3 space dimensions and with time.  This
##        #         function saves a "snapshot" of such a 3D array
##        #         to a "stack file" whenever it is called.
##        #
##        #         The Stack function will work whether argument
##        #         is a 1D profile or already a 3D array.
##        #         The Profile_Var function will work whether its
##        #         its argument is a 1D profile or a 3D array.
##        #         (It is called by Write_Profile.)
##        #---------------------------------------------------------
##        nx = self.nx
##        ny = self.ny
##        SWAP_ENDIAN = self.rti.SWAP_ENDIAN
##
##        #--------------------------------------
##        # Save a subsequence of IN var stacks
##        #--------------------------------------   
##        if (self.SAVE_Q_STACKS):    
##            if (SWAP_ENDIAN): self.q.byteswap(True)
##            Stack(self.q, nx, ny).tofile(self.q_stack_unit)
##        if (self.SAVE_P_STACKS):    
##            if (SWAP_ENDIAN): self.p.byteswap(True)
##            Stack(self.p, nx, ny).tofile(self.p_stack_unit)
##        if (self.SAVE_K_STACKS):    
##            if (SWAP_ENDIAN): self.K.byteswap(True)
##            Stack(self.K, nx, ny).tofile(self.K_stack_unit)
##        if (self.SAVE_V_STACKS):    
##            if (SWAP_ENDIAN): self.v.byteswap(True)
##            Stack(self.v, nx, ny).tofile(self.v_stack_unit)
##            
##    #   save_cubes_old()
##    #-------------------------------------------------------------------

 


