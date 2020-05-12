"""
This file defines a "base class" for evaporation components as well
as functions used by most or all evaporation methods.  That is, all
evaporation components inherit methods from this class.  The methods
of this class should be over-ridden as necessary (especially the
update_ET_rate() method) for different methods of modeling evaporation.
This class, in turn, inherits from the "BMI base class" in BMI_base.py.

See evap_priestley_taylor.py, evap_energy_balance.py, evap_read_file.py
"""
#-----------------------------------------------------------------------
#
#  Copyright (c) 2001-2020, Scott D. Peckham
#
#  May 2020.  Added disable_all_output()
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
#  class evap_component    (inherits from BMI_base)
#
#      (see non-base components for BMI functions)
#
#      ------------------------
#      set_constants()
#      latent_heat_of_evaporation()  # (not used yet)
#      initialize()
#      update()
#      finalize()
#      set_computed_input_vars()
#      -----------------------------
#      check_input_types()
#      check_if_types_match()
#      initialize_computed_vars()
#      -----------------------------
#      update_Qc()                   # (not used yet)
#      update_ET_rate()
#      update_ET_integral()
#      update_water_balance()        # (OBSOLETE, commented out)
#      ------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      ------------------------
#      update_outfile_names()
#      disable_all_output()      # 2020-05-09
#      open_output_files()
#      write_output_files()     #####
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.utils import BMI_base
from topoflow.utils import cfg_files as cfg
from topoflow.utils import model_input
from topoflow.utils import model_output

#-----------------------------------------------------------------------
class evap_component( BMI_base.BMI_component):

    #-------------------------------------------------------------------
    def set_constants(self):

        #---------------------------------
        # From Bob Bolton (Nov. 3, 2009)
        #--------------------------------------------
        # Lv = latent heat of vaporization [J kg-1]
        #--------------------------------------------
        self.mps_to_mmph = np.float64(3600000)
        self.mmph_to_mps = (np.float64(1) / np.float64(3600000))
        self.forever     = np.float64(999999999)  # [minutes]

        #--------------------------------------------
        # Lv = latent heat of vaporization [J kg-1]
        #--------------------------------------------
        # self.Lv = np.float64( 2260000 )    # (at T = 100 C)
        self.Lv = np.float64( 2500000 )              

        
    #   set_constants()
    #-------------------------------------------------------------------
    def latent_heat_of_evaporation(self):

        #----------------------------------------------------------    
        # Notes:  See:  http://en.wikipedia.org/wiki/Latent_heat
        #         Valid for T in [-25, 40] deg C.
        #----------------------------------------------------------
        # latent heat of condensation/evaporation.
        # What about latent heat of vaporization (boiling)?
        #----------------------------------------------------------
        a = np.float64( 2500.8 )
        b = np.float64( -2.36 )
        c = np.float64( 0.0016 )
        d = np.float64( -0.00006 )
        T = self.T_air
        
        self.Lv = a + (b * T) + (c * (T**2)) + (d * (T**3)) # [J g-1]
        self.Lv *= np.float64( 1000 ) # [J kg-1]
                
    #   latent_heat_of_evaporation()
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        if not(SILENT):
            print(' ')
            print('Evaporation component: Initializing...')
        
        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode
        self.cfg_file   = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()    # (12/3/09)
        self.initialize_config_vars()
        # self.read_grid_info()    # NOW IN initialize_config_vars()
        self.initialize_basin_vars()  # (5/14/10)
        #-----------------------------------------
        # This must come before "Disabled" test.
        #-----------------------------------------
        self.initialize_time_vars()
        
        #------------------------------------------------------
        # NB! "Sample steps" must be defined before we return
        #     Check all other process modules.
        #------------------------------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print('Evaporation component: Disabled in CFG file.')
            self.disable_all_output()
            self.ET     = self.initialize_scalar(0, dtype='float64')
            self.vol_ET = self.initialize_scalar(0, dtype='float64')
            self.DONE   = True
            self.status = 'initialized'
            return

        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()

        #-----------------------
        # Initialize variables
        #-----------------------
        self.initialize_computed_vars()  # (such as 'ET')
        self.check_input_types()   # (Uses "mp" vars)
        
        self.open_output_files()
        self.status = 'initialized'
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, dt=-1.0):
        
        #-------------------------------------------------
        # Note: self.ET already set to 0 by initialize()
        #-------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        self.status = 'updating'  # (OpenMI)

        #---------------------------------------
        # Read next ET vars from input files ?
        #-------------------------------------------
        # Note that read_input_files() is called
        # by initialize() and these values must be
        # used for "update" calls before reading
        # new ones.
        #-------------------------------------------
        if (self.time_index > 0):
            self.read_input_files() 

        #---------------------------------------------
        # Qc is used by both the Priestly-Taylor and
        # Energy Balance methods, and both have the
        # required parameters in their CFG files
        # Started using this on 2020-05-10.
        #---------------------------------------------
        self.update_Qc()
                            
        #-------------------------
        # Update computed values 
        #-------------------------
        self.update_ET_rate()
        self.update_ET_integral()

        #----------------------------------------------
        # Write user-specified data to output files ?
        #-----------------------------------------------
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

        self.status = 'finalizing'  # (OpenMI)
        if (self.comp_status == 'Enabled'):
            self.close_input_files()   ##  TopoFlow input "data streams"
            self.close_output_files()
        self.status = 'finalized'  # (OpenMI)

        self.print_final_report(comp_name='Evaporation component')
  
    #   finalize()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt   = np.maximum(self.save_grid_dt,   self.dt)
        self.save_pixels_dt = np.maximum(self.save_pixels_dt, self.dt)
        
    #   set_computed_input_vars()
    #-------------------------------------------------------------------
    def check_input_types(self):

        #----------------------------------------------------
        # Notes: Usually this will be overridden by a given
        #        method of computing ET.
        #----------------------------------------------------
        are_scalars = np.array([
                         # self.is_scalar('d'),
                         #---------------------------------
                         # self.is_scalar('h_table'),
                         #---------------------------------
                         self.is_scalar('T_air'),
                         self.is_scalar('T_surf') ])

        self.ALL_SCALARS = np.all(are_scalars)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        dtype = 'float64'
        self.ET = self.initialize_grid(0, dtype=dtype)
        self.Qc = self.initialize_grid(0, dtype=dtype)
        #---------------------------------------------------------
        self.vol_ET = self.initialize_scalar(0, dtype=dtype)
        
        #------------------------------------------
        # h_table = water table height
        # Assume h_table is always a grid.
        # h_table, dzw and ET must be compatible.
        #------------------------------------------
##        H_IS_GRID = self.is_grid('h_table')
##        if (H_IS_GRID):
##            self.ET = np.zeros([self.ny, self.nx], dtype='float64')            
##        else:
##            self.ET = np.float64(0)
##            print '********* WARNING: water table is not a grid'
##            print '                   but it should be.'

    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def update_Qc(self):
  
        #---------------------------------------------
        # Compute the conductive energy between the
        # surface and subsurface using Fourier's law
        #------------------------------------------------------------
        # Heat flow can be in either direction;
        # If (T_surf > T_soil), then Qc > 0, heat flow is downward.
        # If (T_surf < T_soil), then Qc < 0, heat flow is upward.
        # If Qc > 0, heat flows down into soil and adds to Q_net.
        # Q_net and Qc use the same sign convention.
        #------------------------------------------------------------
        # soil_x is converted from [cm] to [m] when
        # it is read from the GUI and then stored
        #---------------------------------------------
        delta_T = (self.T_surf - self.T_soil_x)
        self.Qc[:] = self.K_soil * delta_T / self.soil_x

        #-----------------------------------------------------------
        # In Zhang et al. (2000), Qc was defined with the opposite
        # sign, but then (Q_net - Qc) was used in the Qet equation.
        # So it seems Qc was defined with a different sign
        # convention that Qn_SW, Qn_LW, etc. (i.e. incoming > 0).
        #-----------------------------------------------------------
        ## delta_T = (self.T_soil_x - self.T_surf)  # Zhang sign.
        ## self.Qc[:] = self.K_soil * delta_T / self.soil_x
                       
    #   update_Qc()
    #-------------------------------------------------------------------
    def update_ET_rate(self):

        #------------------------------------------------------
        # Each component that inherits from evap_base.py must
        # implement its own version of this method.
        #------------------------------------------------------
        print('ERROR: update_ET_rate() method for Evaporation')
        print('       component has not been implemented.')
        
    #   update_ET_rate()
    #-------------------------------------------------------------------
    def update_ET_integral(self):

        #------------------------------------------------
        # Update mass total for ET, sum over all pixels
        #------------------------------------------------   
        volume = np.double(self.ET * self.da * self.dt)  # [m^3]
        if (np.size( volume ) == 1):
            self.vol_ET += (volume * self.rti.n_pixels)
        else:
            self.vol_ET += np.sum(volume)
            
    #   update_ET_integral()
    #-------------------------------------------------------------------
#     def update_water_balance(self):
# 
#         #-------------------------------------------------------
#         # Note: Computed ET values are generally taken to be
#         #       "potential" values which may not be achieved
#         #       if there is not enough water at or near the
#         #       surface.  This function first tries to consume
#         #       the required water from surface water (depth)
#         #       and then goes on to extract water from the
#         #       top soil layer (subsurface).
#         #-------------------------------------------------------
#         # Note: ET = ET rate with units of [m/s].
#         #        d = depth of surface water [m]
#         #        h = water table height above datum
#         #        y = thicknesses [m] of all soil layers
#         #            when using Darcy subsurface flow
#         #-------------------------------------------------------
# 
#         #-------------------------------------------
#         # If Richards' equation is being used for
#         # infiltration, then don't need to remove
#         # water from layers as done in remainder
#         # and y (wetted thicknesses) is not needed
#         #-------------------------------------------
#         # But still need to remove surface water
#         # first !!  This isn't done yet. ********
#         #-------------------------------------------
#         ## if (self.time < 5*self.dt):
#         ##     print 'RICHARDS_EQN =', self.RICHARDS_EQN
#         if (self.RICHARDS_EQN.lower() in ['yes', 'true']):
#             ## print '### Returning from update_water_balance()...'
#             return
#         
#         #-------------------------------------
#         # Depth of water to be removed by ET
#         #-------------------------------------------
#         # (8/25/09) Does it make sense to allow ET
#         # ET and dzw to be scalars ??
#         #-------------------------------------------
#         dzw = (self.dt * self.ET)
#         ## print 'size(dzw) =', np.size(dzw)
#         
#         #----------------
#         # For debugging
#         #----------------
#         #if (np.size(dzw) == 1) then begin
#         #    msg = [' ','ERROR: dzw is not an array. ', ' ']
#         #    result = GUI_Message(msg, /INFO)
#         #    STOP
#         #endif
# 
#         depth = self.depth    # (2/3/13, "d@channel")
#         UPDATE_DEPTH = False
#         
#         wL  = np.where( dzw <= depth )
#         nwL = np.size( wL[0] )
#         wG  = np.where( dzw > depth )
#         nwG = np.size( wG[0] )
# 
#         if (nwL != 0):    
#             #---------------------------------
#             # Reduce the surface water depth
#             #---------------------------------
#             depth[wL]    = (depth[wL] - dzw[wL])
#             UPDATE_DEPTH = True
#             dzw[wL]      = np.float64(0)
#         
#         if (nwG != 0):    
#             #-----------------------------
#             # Save a copy of initial dzw
#             #-----------------------------
#             dzw0 = dzw.copy()
#             
#             #-------------------------------------
#             # Consume all surface water first
#             # This doesn't account for channels.
#             #-------------------------------------
#             dzw[wG]      = dzw[wG] - depth[wG]
#             depth[wG]    = np.float64(0)
#             UPDATE_DEPTH = True
#             
#             #---------------------------------------
#             # Try to take remainder from top layer
#             # Compute water content of top layer
#             #---------------------------------------
#             # Used before 7/13/06
#             #----------------------
#             # p  = gv.soil_P[0]  ;(top layer porosity)
#             # y0 = y[*,*,0]
#             # content_1 = (y0[wG] * p)
#             #---------------------------------------------
#             # self.gp.qs is a 1D array of doubles that
#             # gives theta_sat for each soil layer.
#             # This is taken equal to porosity here.
#             #---------------------------------------------
#             # self.gp.y[0,:,:] is a grid of doubles that
#             # gives the "wetted thickness" of top layer
#             #---------------------------------------------
#             p0 = self.p0       # (2/3/13, new framework)
#             y0 = self.y0       # (2/3/13, new framework)
#             h  = self.h_table  # (2/3/13, new framework)
# 
#             SCALAR_POROSITY = (np.size(p0) == 1)  # (Always True now)
#             if (SCALAR_POROSITY):    
#                 content_1 = (y0[wG] * p0)
#             else:    
#                 content_1 = (y0[wG] * p0[wG])
#             
#             wwL  = np.where( dzw[wG] <= content_1 )
#             nwwL = np.size( wwL[0] )
#             wwG  = np.where( dzw[wG] > content_1 )
#             nwwG = np.size( wwG[0] )
# 
#             #####################################################
#             # See Notes at top regarding "nested WHERE calls".
#             #####################################################
# 
#             #---------------------------------------------
#             # Can get all remaining water from top layer
#             # Reduce the water table height
#             #---------------------------------------------
#             if (nwwL != 0):    
#                 if (SCALAR_POROSITY):
#                     dh = dzw.flat[wwL] / p0
#                     #### dh = dzw[wG][wwL] / p0
#                 else:
#                     dh = dzw.flat[wwL] / p0.flat[wwL]
#                     #### dh = dzw[wG][wwL] / p0[wG][wwL]
# 
#                 h.flat[wwL]   = h.flat[wwL]  - dh
#                 y0.flat[wwL]  = y0.flat[wwL] - dh
#                 dzw.flat[wwL] = np.float64(0)     # (not really needed ?)
#                 
# ##                h[wG][wwL]   = h[wG][wwL] - dh
# ##                y0[wG][wwL]  = y0[wG][wwL] - dh
# ##                dzw[wG][wwL] = np.float64(0)   # (not really needed ?)
#             
#             #-----------------------------------------------
#             # Can't get all remaining water from top layer
#             #-----------------------------------------------
#             # Get what is available, and then redefine ET
#             # for mass balance consistency
#             #-----------------------------------------------
#             if (nwwG != 0):
#                 dh = y0.flat[wwG]
#                 h.flat[wwG]   = h.flat[wwG] - dh
#                 y0.flat[wwG]  = np.float64(0)
#                 dzw.flat[wwG] = dzw.flat[wwG] - content_1[wwG]
#                 #################################################
#                 ##### Is there a problem in above line with
#                 ##### content_1[wwG] part ???
#                 #------------------------------------------------
#                 dzw_used    = dzw0.flat[wwG] - dzw.flat[wwG]
#                 self.ET.flat[wwG] = (dzw_used / self.dt)
#                 
# ##                dh = y0[wG][wwG]
# ##                h[wG][wwG]   = h[wG][wwG] - dh
# ##                y0[wG][wwG]  = np.float64(0)
# ##                dzw[wG][wwG] = dzw[wG][wwG] - content_1[wwG]
# ##                #--------------------------------------------
# ##                dzw_used    = dzw0[wG][wwG] - dzw[wG][wwG]
# ##                self.ET[wG][wwG] = (dzw_used / self.dt)
#       
#             #-------------------------
#             # Replace top layer in y
#             #-------------------------
#             print '    ET component changing "h_table" in GW component.'
# ##            print '       type(y0) =', type(y0)
# ##            print '       type(h)  =', type(h)
#             self.set_port_data('y[0,:,:]', y0, self.gp)
#             self.set_port_data('h_table', h, self.gp)
#             
#     #   update_water_balance()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        #------------------------------------------------------
        # Each component that inherits from evap_base.py must
        # implement its own versions of these.
        #------------------------------------------------------
        print('ERROR: open_input_files() for Evaporation component')
        print('       has not been implemented.')

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        print('ERROR: read_input_files() for Evaporation component')
        print('       has not been implemented.')
        
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):

        print('ERROR: close_input_files() for Evaporation component')
        print('       has not been implemented.')

    #   close_input_files()
    #-------------------------------------------------------------------
    def update_outfile_names(self):

        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.ET_gs_file = (self.out_directory + self.er_gs_file)
        #---------------------------------------------------------
        self.ET_ts_file = (self.out_directory + self.er_ts_file)

    #   update_outfile_names()
    #-------------------------------------------------------------------
    def disable_all_output(self):

        self.SAVE_ER_GRIDS  = False
        self.SAVE_ER_PIXELS = False
            
    #   disable_all_output()   
    #-------------------------------------------------------------------  
    def open_output_files(self):

        #---------------------------------------------------------
        # Note:  Qc (conduction heat flux), from Priestly-Taylor
        #        component, could also be saved.
        #---------------------------------------------------------
        model_output.check_netcdf()
        self.update_outfile_names()

        #--------------------------------------
        # Open new files to write grid stacks
        #--------------------------------------
        if (self.SAVE_ER_GRIDS):
            model_output.open_new_gs_file( self, self.ET_gs_file, self.rti,
                                           var_name='ET',
                                           long_name='evaporation_rate',
                                           units_name='mm/hr')
                                           ### units_name='m/s')
            
        #--------------------------------------
        # Open new files to write time series
        #--------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_ER_PIXELS):
            model_output.open_new_ts_file( self, self.ER_ts_file, IDs,
                                           var_name='ET',
                                           long_name='evaporation_rate',
                                           units_name='mm/hr')

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
##        if (SAVE_ER_GRIDS  == False) and  \
##           (SAVE_ER_PIXELS == False): return
           
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
    #---------------------------------------------------------------------
    def close_output_files(self):
    
        if (self.SAVE_ER_GRIDS):  model_output.close_gs_file( self, 'ET')
        #-----------------------------------------------------------------
        if (self.SAVE_ER_PIXELS): model_output.close_gs_file( self, 'ET')  

    #   close_output_files()   
    #---------------------------------------------------------------------  
    def save_grids(self):
        
        #-----------------------------------
        # Save grid stack to a netCDF file
        #---------------------------------------------
        # Note that add_grid() methods will convert
        # var from scalar to grid now, if necessary.
        #--------------------------------------------- 
        if (self.SAVE_ER_GRIDS):
            ET_mmph = self.ET * self.mps_to_mmph    # (Bolton 28 Aug)
            model_output.add_grid( self, ET_mmph, 'ET', self.time_min )

    #   save_grids()            
    #---------------------------------------------------------------------  
    def save_pixel_values(self):

        IDs  = self.outlet_IDs
        time = self.time_min   ########
         
        if (self.SAVE_ER_PIXELS):
            ET_mmph = self.ET * self.mps_to_mmph    # (Bolton 28 Aug)
            model_output.add_values_at_IDs( self, time, ET_mmph, 'ET', IDs )

    #   save_pixel_values()
    #---------------------------------------------------------------------

    
        
