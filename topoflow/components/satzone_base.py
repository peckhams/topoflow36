
# SEARCH FOR THIS!  It doesn't work for Python lists! (2/19/13)
# self.qs_layer_0 = self.qs[0]
##########################################


## Copyright (c) 2001-2016, Scott D. Peckham
## January 2009   (converted from IDL)
## May, July, August 2009
## May 2010 (changes to initialize() and read_cfg_file()

###########################################################################

## Clean up initialize() for case where method=0 ??

## NB!  "adjust_flow_depths()" has a "conceptual bug".  See notes there.

## NB!  There is a fair amount of IDL "where subscripting" used here in
##      update_water_table().  I2PY does not handle this correctly yet,
##      so it was fixed by hand using ".flat", etc.  CHECK MORE.

###########################################################################

#-----------------------------------------------------------------------
#  Notes:  This file defines a "base class" for groundwater
#          components as well as functions used by most or
#          all groundwater methods.  The methods of this class
#          should be over-ridden as necessary for different
#          methods of modeling groundwater flow.
#-----------------------------------------------------------------------

# (5/7/09)  "Nested WHERE calls" work differently
# in numpy than in IDL. For a 1D array called "a":
#
#     >>> a = np.arange(11)-5
#     >>> w = np.where(a < 0)
#     >>> w2 = np.where(a[w] > -3)
#     >>> print a[w[w2]]  # (this gives an error)
#     >>> print a[w2]     # (this works)
#     >>> print a[w][w2]  # (this works, too)

# For a 2D array called "a":
#
#     >>> a = np.arange(9) - 4
#     >>> a = a.reshape(3,3)
#     >>> w = np.where(a < 0)
#     >>> w2 = np.where(a[w] > -3)
#     >>> print a[w[w2]]    # TypeError: tuple indices must be integers
#     >>> print a[w2]       # IndexError: index (3) out of range (0<=index<=2) in dimension 0
#     >>> a[w][w2] = 99     # No error, but this doesn't work.
#     >>> print a[w][w2]    # (this works)
#     >>> print a.flat[w2]  # (this works, same result as last line)
#     >>> a.flat[w2] = 99   # (this works)
#     >>> a.flat[w2] = [-2,-1]  # (this works)
#     >>> np.put(a, w2, 99)  # (this works)
               
#-----------------------------------------------------------------------
#
#  class satzone_component
#
#      set_constants()
#      initialize()
#      update()
#      finalize()
#      ---------------------------
#      check_input_types()           ## (not written; needed?)
#      set_computed_input_vars()
#      --------------------------
#      initialize_layer_vars()        # (5/11/10)
#      initialize_computed_vars()
#      initialize_water_table()
#      initialize_wetted_thicknesses()
#      initialize_d8_vars()
#      adjust_flow_depths()                   ## (has a bug)
#      ----------------------------------
#      update_Sh()
#      update_Q_gw()
#      update_top_layer_for_ET()     # (9/25/14)
#      update_water_table()
#      update_seep_rate()
#      update_GW_integral()          # (GW = seeprate)
#      ------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      ------------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
#      get_total_darcy_layer_flow()   # (OBSOLETE NOW ??)
#
#  Functions:
#      Darcy_Flow()
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.utils import cfg_files as cfg
from topoflow.utils import model_input
from topoflow.utils import model_output
from topoflow.utils import tf_utils
from topoflow.utils import rti_files

## from topoflow.utils import BMI_base
from topoflow.components import infil_base   # (for build_layered_var())

#-------------------------------------------------------
# NOTE:  Do not import "d8_base" itself, it won't work
#-------------------------------------------------------
from topoflow.components import d8_global as d8_base    # (11/16/16)
## from topoflow.utils import tf_d8_base as d8_base

#-----------------------------------------------------------------------
### class satzone_component( BMI_base.BMI_component ):
# (11/16/16) Inherit from infil_base to get build_layered_var()
class satzone_component( infil_base.infil_component ):


    #-----------------------------------------------------------
    # Notes:  h0_table = init. elevation of water table [m]
    #         d_thaw   = depth to nonfrozen soil [m]
    #                    (everything below is thawed)
    #
    # NB!    (3/16/07) dt is now set to zero initially, so
    #        that Read_Run_Info_Panel can determine if user
    #        has set the value (e.g. by loading a saved value).
    #        If unset, then it is set by Read_Run_Info_Panel.
    #        If GUI is not used, it is set in the input file.
    #-----------------------------------------------------------
    # Notes: h_snow is needed by the Bulk_Exchange_Coeff
    #        function to adjust reference height, z.
    #        Do we need h0_snow here ??
    #-----------------------------------------------------

    #-------------------------------------------------------------------
    def set_constants(self):

        #------------------------
        # Define some constants
        #------------------------
        self.nodata   = np.float64(-9999)
        self.RICHARDS = False

    #   set_constants()
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        if not(SILENT):
            print(' ')
            print('Groundwater component: Initializing...')
            
        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-------------------------------------------------------------
        # NOTE!  initialize_config_vars() calls read_config_file(),
        #        which now calls initialize_layer_vars(). (11/15/16)
        #-------------------------------------------------------------
        self.set_constants()  
        ## self.initialize_layer_vars()  # (5/11/10)
        self.initialize_config_vars() 
        # self.read_grid_info()    # NOW IN initialize_config_vars()
        self.initialize_basin_vars()  # (5/14/10)
        #-----------------------------------------
        # This must come before "Disabled" test.
        #-----------------------------------------
        self.initialize_time_vars()
        
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print('Groundwater component: Disabled in CFG file.')
            
            #-------------------------------------------------------
            # Other processes, such as evap, may still need DEM.
            #-------------------------------------------------------
            # The default data type for model_input.read_next() is
            # Float32, but the DEM may have another type.
            #-------------------------------------------------------
            self.elev_file = self.topo_directory + self.elev_file
            self.elev_unit = model_input.open_file(self.elev_type,
                                                   self.elev_file)
            DEM_dtype = rti_files.get_numpy_data_type( self.rti.data_type )
            elev = model_input.read_next(self.elev_unit, self.elev_type,
                                         self.rti, dtype=DEM_dtype)
            model_input.close_file(self.elev_unit)
            if (elev is not None): self.elev = elev

            #---------------------------------------------------
            # Initialize GW, vol_GW, h_table, h_last, y, etc.
            #--------------------------------------------------
            self.initialize_computed_vars()  # (2/19/13)
            ##################################################
      
            self.DONE = True
            self.status = 'initialized'  # (OpenMI 2.0 convention)
            return
        
        #----------------------------
        # Initialize more variables
        #----------------------------
        self.initialize_d8_vars()
        
        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()
  
        #----------------------------------------------
        # Must come before initialize_computed_vars()
        # because it uses ALL_SCALARS.
        #-----------------------------------------------------
        # Used in other components, but not here. (11/16/16)
        #-----------------------------------------------------
        ## self.check_input_types()
      
        #------------------------------------------------
        # DEM was read by read_input_files().
        # Data type of DEM need not be 'FLOAT'.
        #------------------------------------------------
        # Note that infiltration.update() computes Rg,
        # the rate at which water reaches water table.
        #------------------------------------------------
        self.initialize_computed_vars()       # (vol_GW, GW, etc.)
     
        self.open_output_files()
        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, dt=-1.0):
         
        #-------------------------------------------------
        # Note: self.GW already set to 0 by initialize()
        #-------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        self.status = 'updating'  # (OpenMI)

        #---------------------------------------
        # Read next GW vars from input files ?
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
        self.update_Sh()
        self.update_Q_gw()
        # print 'CALLING update_top_layer_for_ET()...'
        self.update_top_layer_for_ET()
        # print 'CALLING update_water_table()...'
        self.update_water_table()
        # print 'CALLING update_seep_rate()...'
        self.update_seep_rate()
        # print 'CALLING update_GW_integral()...'
        self.update_GW_integral()

        #----------------------------------------------
        # Write user-specified data to output files ?
        #-----------------------------------------------
        # Components use own self.time_sec by default.
        #-----------------------------------------------
        # print 'CALLING write_output_files()...'
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

        self.status = 'finalizing'  # (OpenMI 2.0 convention)
        self.close_input_files()   ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'  # (OpenMI 2.0 convention)

        self.print_final_report(comp_name='Groundwater component')
    
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
    def initialize_layer_vars(self):

        #----------------------------------------------------------
        # Notes: initialize_config_vars() calls read_cfg_file().
        #        If read_cfg_file() finds a variable "n_layers",
        #        then it calls initialize_layer_vars() so that
        #        subsequent layer variables - as indicated by a
        #        subscript in the CFG file - can be read directly
        #        into a list or array.
        #----------------------------------------------------------
        ## n_layers = 6    # (before 11/16/16)
        n_layers = self.n_layers
        
        #---------------------------------------
        # Get arrays for soil layer parameters
        #---------------------------------------
        self.Ks_type = np.zeros(n_layers, dtype='<U200')
        self.qs_type = np.zeros(n_layers, dtype='<U200')
        self.th_type = np.zeros(n_layers, dtype='<U200')
        #------------------------------------------------------
        self.Ks_file = np.zeros(n_layers, dtype='<U200')
        self.qs_file = np.zeros(n_layers, dtype='<U200')
        self.th_file = np.zeros(n_layers, dtype='<U200')
        #---------------------------------------------------------
        # Note: self.Ks is a Python list.  Initially, each entry
        # is a numpy scalar (type 'np.float64').  However, we
        # can later change any list entry to a scalar or grid
        # (type 'np.ndarray'), according to its "Ks_type".
        #---------------------------------------------------------
        # While CFG file for Richards 1D uses "Ks_val[0]", the
        # CFG file for Green-Ampt, etc. just uses "Ks[0]".
        # Too late to change it to be consistent.
        # Later, build_layered_vars() will use these lists to
        # build ndarrays *with the same names*.
        #---------------------------------------------------------  
        self.Ks  = list( np.zeros(n_layers, dtype='float64') )
        self.qs  = list( np.zeros(n_layers, dtype='float64') )
        self.th  = list( np.zeros(n_layers, dtype='float64') )

        #----------------------------------------------------------------
        # Note: The variables qs, th and y are ndarrays.  If we define
        #       another variable as a slice or subset of these, such as
        #       qs_top = qs[0], or y_top = y[0,:,:], then they will
        #       also change whenever the main ndarray changes.
        #----------------------------------------------------------------
        #       This doesn't work for Python lists, however! (2/19/13)
        #----------------------------------------------------------------

        #--------------------------------------------------
        # These are not used anywhere;  for illustration?
        #--------------------------------------------------
#         self.qs_layer_0 = self.qs[0]
#         self.qs_layer_1 = self.qs[1]
#         self.qs_layer_2 = self.qs[2]
#         #-----------------------------
#         self.th_layer_0 = self.th[0]
#         self.th_layer_1 = self.th[1]
#         self.th_layer_2 = self.th[2]
        
    #   initialize_layer_vars()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):
        
        self.initialize_water_table()         # (h_table)
        self.initialize_wetted_thicknesses()  # (y)
        
        self.h_last = self.h_table.copy()
   
        #----------------------------------
        # Do this here?  See update_Q_gw.
        #----------------------------------
        self.vol_GW = self.initialize_scalar(0, dtype='float64')
        self.GW     = self.initialize_grid(0, dtype='float64')
        self.Q_gw   = self.initialize_grid(0, dtype='float64')
 
        #-----------------------------------------------------
        # Compute dz as 1D array from scalars in self.dz_val      
        #-----------------------------------------------------
        # Compute the z-vector, for plotting profiles
        #----------------------------------------------
        self.dz_val = self.th   ###### (11/16/16)
        self.build_layer_z_vector()

        #------------------------------------------------
        # Now build a 1D or 3D array for each input var
        #----------------------------------------------------------
        # (3/12/08) Same code should work if (self.n_layers eq 1)
        #----------------------------------------------------------
        # Convert from lists to arrays; same name. (11/15/16)
        #----------------------------------------------------------
        self.Ks  = self.build_layered_var( self.Ks )
        self.qs  = self.build_layered_var( self.qs )
        self.th  = self.build_layered_var (self.th )
       
        #---------------------------------
        # If water table > land surface,
        # increment the flow depth grid.
        #---------------------------------
        # self.adjust_flow_depths()         #### NOT READY YET ####

    #   initialize_computed_vars()   
    #-------------------------------------------------------------------
    def initialize_water_table(self):

        #----------------------------------------
        # Convert h_table from scalar to grid ?
        #----------------------------------------
        nh = np.size( self.h0_table )
        if (nh == 1):
            h0_scalar = self.h0_table.copy() 
            self.h_table = self.initialize_grid(h0_scalar, dtype='float64')   
        else:    
            self.h_table = self.h0_table.copy()

    #   initialize_water_table()
    #-------------------------------------------------------------------
    def initialize_wetted_thicknesses(self):

        #-----------------------------------------------------------
        # Notes: z    = elevation of land surface [m]
        #        h    = elevation of water table [m]
        #        (z-h) = depth to water table [m]
        #        diff = (partial sum of soil thicknesses -
        #               depth to water table)

        #        Let ti be the thickness of the ith soil layer,
        #        starting from the surface, and yi the wetted
        #        thickness of this layer, which must be between
        #        0 and ti.  If (z gt h), then we have:
        #        sum_{i=1}^n yi = [sum_{i=1}^n ti] ?? (z-h)

        #        e.g.
        #        y1 = (t1 ?? (z-h)) > 0 < t1
        #        y2 = ((t1+t2) ?? y1 - (z-h)) > 0 < t2
        #        y3 = ((t1+t2+t3) ?? (y1+y2) ?? (z-h)) > 0 < t3

        #        This function was tested by defining some vars
        #        as follows, and both methods seemed to work.

        #        IDL> gv = {n_layers:5, soil_thick:(findgen(5)+1)*0.1}
        #        IDL> z  = fltarr(3,3) + 2.0
        #        IDL> h  = fltarr(3.3) + 1.85
        #        IDL> y  = Wetted_Thicknesses(gv, z, h)
        #        IDL> print, y[*,*,0]
        #        IDL> print, y[*,*,1]  ;(etc.)

        #        If (h gt z), then testing shows that (yk eq tk)
        #        for all layers, as it should.
        #-----------------------------------------------------------
        nx = self.nx
        ny = self.ny
        self.y = np.zeros([self.n_layers, ny, nx], dtype='float64')

        #----------------------------------------------------------------
        # Note: The variables qs, th and y are ndarrays.  If we define
        #       another variable as a slice or subset of these, such as
        #       qs_top = qs[0], or y_top = y[0,:,:], then they will
        #       also change whenever the main ndarray changes. (2/19/13)
        #----------------------------------------------------------------    
        self.y_layer_0 = self.y[0,:,:]
        self.y_layer_1 = self.y[1,:,:]
        self.y_layer_2 = self.y[2,:,:]
        
        #--------------------------------
        # Original way: Prior to 3/1/04
        # Seems to agree with new way?
        #--------------------------------
        # Compute wetted thickness, y,
        # for each soil layer
        #-------------------------------
        tsum = (self.h_table - self.elev)
        for k in range(self.n_layers):
            #-------------------------------
            # Before 7/7/06. (tk = scalar)
            #-------------------------------
            # tk   = gv.soil_thick[k]
            # tsum = (tsum + tk)
            # y[*,*,k] = (tsum > 0.0) < tk)
            
            #--------------------------------
            # 7/7/06.  Allow scalar or grid
            #--------------------------------
            tsum  += self.th[k]
            grid_k = np.maximum(tsum, np.float64(0))
            self.y[k,:,:] = np.minimum(grid_k, self.th[k])
        
        #-------------------------------
        # Compute wetted thickness, y,
        # for each soil layer
        #-------------------------------
        # Slightly more costly method.
        #-------------------------------
        #ysum = zeros([ny, nx], dtype='Float32')
        #tsum = (self.h_table - self.elev)            ;(negative, if z > h)
        #for k in xrange(self.n_layers):
        #    tk = self.th[k]            ;(is now a scalar)
        #    tsum += tk
        #    self.y[k,:,:] = np.minimum(np.maximum(tsum - ysum,  0.0), tk)
        #    ysum += y[k,:,:]
    
    #   initialize_wetted_thicknesses()
    #-------------------------------------------------------------------
    def initialize_d8_vars(self):

        #---------------------------------------------
        # Compute and store a variety of (static) D8
        # flow grid variables.  Embed structure into
        # the "channel_base" component.
        #---------------------------------------------
        self.d8 = d8_base.d8_component()
 
        #--------------------------------------------------         
        # D8 component builds its cfg filename from these  
        #-------------------------------------------------------------
        # Note:  This D8 component is serving a satzone component
        #        that has already been instantiated and knows its
        #        directory and prefix information.  So we can build
        #        the correct D8 cfg_file name from that info.  It
        #        will then read path_info CFG file to get other info.
        #-------------------------------------------------------------
        cfg_file = (self.case_prefix + '_d8_global.cfg')
        cfg_file = (self.cfg_directory + cfg_file)
        self.d8.initialize( cfg_file=cfg_file, SILENT=self.SILENT, \
                            REPORT=self.REPORT )
        #---------------------------------------------------------     
#         self.d8.site_prefix  = self.site_prefix
#         self.d8.in_directory = self.in_directory
#         self.d8.initialize( cfg_file=None, SILENT=self.SILENT, \
#                             REPORT=self.REPORT )

        #---------------------------------------------------
        # The next 2 "update" calls are needed when we use
        # the new "d8_base.py", but are not needed when
        # using the older "tf_d8_base.py".      
        #---------------------------------------------------
        self.d8.update(self.time, SILENT=False, REPORT=True)

        #----------------------------------------------------------- 
        # Note: This is also needed, but is not done by default in
        #       d8.update() because it hurts performance of Erode.
        #----------------------------------------------------------- 
        self.d8.update_noflow_IDs()       

    #   initialize_d8_vars()
    #-------------------------------------------------------------------
    def adjust_flow_depths(self):

        #---------------------------------------------------------        
        # Note: This is meant to be called from the method
        #       initialize_computed_vars(), but isn't called yet
        #       because it is not implemented correctly.
        #       It is only needed if there are places where the
        #       water table is above the land surface.
        #---------------------------------------------------------
        
        #------------------------------------------------------------
        # If water table > land surface, increment flow depth grid.
        #------------------------------------------------------------
        diff = (self.h_table - self.elev)
        ### w = np.where(diff > 0)
        w    = np.where(np.logical_and(np.logical_and((diff > 0), \
                                           (self.elev > self.nodata)), \
                                           (np.isfinite(self.elev))))
        nw = np.size(w[0])
        if (nw != 0):
            d = self.d   ## (2/3/13, new framework)
            #######################################
            #  THIS IS NOT CORRECT FOR CHANNELS
            #######################################            
            d[w] = d[w] + diff[w]
            self.cp.set_grid_double('d', d)  ###########
            
        #--------------------------
        # For debugging & testing
        #--------------------------
        dmin = d.min()
        dmax = d.max()
        dstr = str(dmin) + ', ' + str(dmax) + ')'
        print('Initial depths due to water table: (dmin, dmax) = (' + dstr)

    #   adjust_flow_depths()
    #-------------------------------------------------------------------
    def update_Sh(self):

        #-----------------------------------
        # Compute water table slope from h
        #-----------------------------------
        # NB!  h is assumed to be a grid.
        #-----------------------------------
        # NB!  Q is zero where Sh is zero.
        #-------------------------------------------
        # NB!  Flow direction is still assumed to
        #      be given by the DEM's D8 flow grid.
        #-------------------------------------------
        # Sh = Free_Surface_Slope(np.float64(0), self.h_table, \
        #                         self.d8.ds, self.d8.pIDs)

        #-----------------------------------------------------------
        # Notes:  It is assumed that the flow directions don't
        #         change even though the free surface is changing.
        #-----------------------------------------------------------
        # Notes:  Check that ds includes channel sinuosity.
        #-----------------------------------------------------------
        delta_h = (self.h_table - self.h_table[self.d8.parent_IDs])
        self.Sh =  delta_h / self.d8.ds
        
    #   update_Sh()
    #-------------------------------------------------------------------
    def update_Q_gw(self):

        #-----------------------------------------
        # NB!  initialize_water_table makes sure
        #      that h_table is a grid.
        #-----------------------------------------
        # self.Q_gw = np.zeros([self.ny, self.nx], dtype='Float64')
        self.Q_gw = np.minimum(self.Q_gw, 0)  # (zero out Q_gw ??)
        
        for k in range(self.n_layers):
            #--------------------------------------
            # Add Q for this layer, via Darcy law
            #--------------------------------------
            self.Q_gw += (self.Ks[k] * self.Sh * self.d8.dw * self.y[k,:,:])
            
    #   update_Q_gw()
    #-------------------------------------------------------------------
    def update_seep_rate(self):

        #--------------------------------------------------------
        # Notes: h_table  = elevation of water table [m]
        #        h0_table = init. elev. of water table
        #        dt       = GW timestep [sec]
        #--------------------------------------------------------
        #        h  = elevation of water table [m]
        #        z  = elevation of bed [m]
        #        y  = wetted thicknesses [m] of all soil
        #             layers
        #        Rg = rate at which water from the surface
        #             arrives at the water table [m/s]
        #             Now set to infiltration rate.
        #             ########################################
        #        Sh = slope of the water table [m/m]
        #        dw = flow width grid [m]
        #        ds = flow length grid [m]
        #        da = pixel area [m^2]  (grid or scalar)

        #        GW = seepage rate (to/from surface) [m/s]
        #--------------------------------------------------------
        #        METHOD 1:  (Maybe the only method)

        #        Where (h lt z) there is no contribution,
        #        but where (h gt z), an amount (h - h_last)
        #        must be contributed over a time dt.
        #--------------------------------------------------------
        # NB!  This function is a bit different than the
        #      others in that GW is computed in the same way
        #      each time, but Q_gw may be computed by different
        #      methods.
        #--------------------------------------------------------

        #-----------------------------------------------
        # Compute the "seep rate" as the rate at which
        # the water table has risen over a time step.
        # Allow sign to be positive or negative?
        #-----------------------------------------------
        # Is this independent of method ?
        # If not, don't do it here.
        #-----------------------------------
        dh_dt   = (self.h_table - self.h_last) / self.dt       # [m/s]
        self.GW = (self.h_table > self.elev) * dh_dt
        ## np.maximum(self.GW, 0.0, self.GW)    ;(force to be positive)
        
        #--------------------------------
        # Redefine h_last for next time
        #--------------------------------
        self.h_last = self.h_table.copy()
          
    #   update_seep_rate()
    #-------------------------------------------------------------------
    def update_GW_integral(self):

        #------------------------------------------------
        # Update mass total for GW, sum over all pixels
        #------------------------------------------------   
        volume = np.double(self.GW * self.da * self.dt)  # [m^3]
        if (np.size( volume ) == 1):
            self.vol_GW += (volume * self.rti.n_pixels)
        else:
            self.vol_GW += np.sum(volume)

    #   update_GW_integral()
    #-------------------------------------------------------------------
    def update_top_layer_for_ET(self):

        #-------------------------------------------------------
        # Note: Computed ET values are generally taken to be
        #       "potential" values which may not be achieved
        #       if there is not enough water at or near the
        #       surface.  The Channels component includes ET
        #       in excess rainrate, R, and ET may therefore
        #       offset contributions from P, SM, MR and GW.
        #       If R < 0, the Channels component will remove
        #       water from the volume in the channel before
        #       computing the channel flow depth.
        #-------------------------------------------------------
        #       This function attempts to consume water from
        #       the top soil layer (subsurface).
        #-------------------------------------------------------
        # Note: ET = ET rate with units of [m/s].
        #        d = depth of surface water [m]
        #        h = water table height above datum
        #        y = thicknesses [m] of all soil layers
        #            when using Darcy subsurface flow
        #-------------------------------------------------------

        #############        
        return          # (not ready yet)
        #############
        
        #-----------------------------------------------------
        # Potential depth of water that can be removed by ET
        #-----------------------------------------------------
        # (8/25/09) Does it make sense to allow ET
        # and dzw to be scalars ??
        #-------------------------------------------
        dzw = (self.dt * self.ET)
        ## print 'size(dzw) =', np.size(dzw)
        
        #----------------
        # For debugging
        #----------------
        #if (np.size(dzw) == 1) then begin
        #    msg = [' ','ERROR: dzw is not an array. ', ' ']
        #    result = GUI_Message(msg, /INFO)
        #    STOP
        #endif

        depth = self.depth    # (2/3/13, "d@channel")
        UPDATE_DEPTH = False
        
        wL  = np.where( dzw <= depth )
        nwL = np.size( wL[0] )
        wG  = np.where( dzw > depth )
        nwG = np.size( wG[0] )

        if (nwL != 0):    
            #---------------------------------
            # Reduce the surface water depth
            #---------------------------------
            depth[wL]    = (depth[wL] - dzw[wL])
            UPDATE_DEPTH = True
            dzw[wL]      = np.float64(0)
        
        if (nwG != 0):    
            #-----------------------------
            # Save a copy of initial dzw
            #-----------------------------
            dzw0 = dzw.copy()
            
            #-------------------------------------
            # Consume all surface water first
            # This doesn't account for channels.
            #-------------------------------------
            dzw[wG]      = dzw[wG] - depth[wG]
            depth[wG]    = np.float64(0)
            UPDATE_DEPTH = True
            
            #---------------------------------------
            # Try to take remainder from top layer
            # Compute water content of top layer
            #---------------------------------------
            # Used before 7/13/06
            #----------------------
            # p  = gv.soil_P[0]  ;(top layer porosity)
            # y0 = y[*,*,0]
            # content_1 = (y0[wG] * p)
            #---------------------------------------------
            # self.gp.qs is a 1D array of doubles that
            # gives theta_sat for each soil layer.
            # This is taken equal to porosity here.
            #---------------------------------------------
            # self.gp.y[0,:,:] is a grid of doubles that
            # gives the "wetted thickness" of top layer
            #---------------------------------------------
            p0 = self.p0       # (2/3/13, new framework)
            y0 = self.y0       # (2/3/13, new framework)
            h  = self.h_table  # (2/3/13, new framework)

            SCALAR_POROSITY = (np.size(p0) == 1)  # (Always True now)
            if (SCALAR_POROSITY):    
                content_1 = (y0[wG] * p0)
            else:    
                content_1 = (y0[wG] * p0[wG])
            
            wwL  = np.where( dzw[wG] <= content_1 )
            nwwL = np.size( wwL[0] )
            wwG  = np.where( dzw[wG] > content_1 )
            nwwG = np.size( wwG[0] )

            #####################################################
            # See Notes at top regarding "nested WHERE calls".
            #####################################################

            #---------------------------------------------
            # Can get all remaining water from top layer
            # Reduce the water table height
            #---------------------------------------------
            if (nwwL != 0):    
                if (SCALAR_POROSITY):
                    dh = dzw.flat[wwL] / p0
                    #### dh = dzw[wG][wwL] / p0
                else:
                    dh = dzw.flat[wwL] / p0.flat[wwL]
                    #### dh = dzw[wG][wwL] / p0[wG][wwL]

                h.flat[wwL]   = h.flat[wwL]  - dh
                y0.flat[wwL]  = y0.flat[wwL] - dh
                dzw.flat[wwL] = np.float64(0)     # (not really needed ?)
                
##                h[wG][wwL]   = h[wG][wwL] - dh
##                y0[wG][wwL]  = y0[wG][wwL] - dh
##                dzw[wG][wwL] = np.float64(0)   # (not really needed ?)
            
            #-----------------------------------------------
            # Can't get all remaining water from top layer
            #-----------------------------------------------
            # Get what is available, and then redefine ET
            # for mass balance consistency
            #-----------------------------------------------
            if (nwwG != 0):
                dh = y0.flat[wwG]
                h.flat[wwG]   = h.flat[wwG] - dh
                y0.flat[wwG]  = np.float64(0)
                dzw.flat[wwG] = dzw.flat[wwG] - content_1[wwG]
                #################################################
                ##### Is there a problem in above line with
                ##### content_1[wwG] part ???
                #------------------------------------------------
                dzw_used    = dzw0.flat[wwG] - dzw.flat[wwG]
                
                ##############################################
                # self.ET.flat[wwG] = (dzw_used / self.dt)
                ##############################################
                                
##                dh = y0[wG][wwG]
##                h[wG][wwG]   = h[wG][wwG] - dh
##                y0[wG][wwG]  = np.float64(0)
##                dzw[wG][wwG] = dzw[wG][wwG] - content_1[wwG]
##                #--------------------------------------------
##                dzw_used    = dzw0[wG][wwG] - dzw[wG][wwG]
##                self.ET[wG][wwG] = (dzw_used / self.dt)
      
            #-----------------------------------------------
            # Replace top layer in y (saturated thickness)
            #-----------------------------------------------
            self.y[0,:,:]   = y0
            self.h_table[:] = h
##            print '       type(y0) =', type(y0)
##            print '       type(h)  =', type(h)  
    
    #   update_top_layer_for_ET()
    #-------------------------------------------------------------------
    def update_water_table(self):

        #----------------------------------------------------
        # Notes: h  = elevation of water table [m]
        #        h2 = temp version of h
        #        Q_gw = total subsurface flux [m^3/s]
        #        Rg = rate at which water from the surface
        #             arrives at the water table [m/s]
        #        da = pixel area [m^2]
        #        dt = GW timestep [sec]
        #        w1 = IDs of pixels that flow in direction 1
        #        p1 = IDs of parent pixels for "w1 pixels"

        # Note:  h and wetted-depths, y, are updated
        #-------------------------------------------------------------
        # Notes: There seems to be an implicit assumption here
        #        that Ks is nondecreasing towards the surface.
        #        Once an element is saturated there is no internal
        #        storage and the amount of water flowing in through
        #        its faces must equal the amount that is flowing out.
        #        So if the amount flowing in from the sides exceeds
        #        the amount flowing out, (which will occur when the
        #        flow is convergent) then the excess must flow
        #        through the upper or lower faces.  With saturated
        #        soil all the way down to an impermeable bedrock
        #        boundary, this means that it must flow through the
        #        upper face. But all that enters in a given time step
        #        can only flow through the upper face if Ks in the
        #        element above is high enough to accommodate it.
        #-------------------------------------------------------------

        #------------------------------------------
        # Compute wetted-depth, y, for each layer
        #------------------------------------------
        # Now passed as a variable in & out
        #------------------------------------------
        #dims  = np.size(z, /DIM)
        #ncols = dims[0]
        #nrows = dims[1]
        #y     = fltarr(ncols, nrows)
        #diff  = -(z - h)
        #for k=0, (n_layers - 1) do begin
        #    diff = (diff + t[k])
        #    y[*,*,k] = (diff > 0.0) < t[k]
        #endfor
        
        #--------------------------------------------
        # Compute dzw = total amount of water to be
        # added to or removed from the soil column
        # during the subsurface flow timestep
        #--------------------------------------------
        # Initialize dzw with outflow term.
        # Doesn't involve neighbor pixels.
        #------------------------------------
        Rg  = self.Rg   # (using new framework, 5/18/12)
        
        dzw = self.dt * (Rg - self.Q_gw / self.da)    ### USE gp.dt vs. main_dt !! ###
        
        #----------------
        # For debugging
        #----------------
##        if (np.size( dzw ) == 1):    
##            msg = array([' ', 'ERROR: dzw is not an array. ', ' '])
##            result = GUI_Message(msg, INFO=True)
##            sys.exit()

        #----------------
        # For debugging
        #----------------
        #print 'dt =', dt
        #print 'Rg =', Rg
        #------------------------------------
        #print 'Qg_min  =', self.Q_gw.min()
        #print 'Qg_max  =', self.Q_gw.max()
        #------------------------------------            
        dz_min = dzw.min()
        dz_max = dzw.max()  #***********************
        print('   dz_min = ' + str(dz_min))
        print('   dz_max = ' + str(dz_max))
        ## print ' '

        #-------------------------------------------
        # Local synonyms  (Any performance hit ??)
        #-------------------------------------------
        p1 = self.d8.p1   ;  w1 = self.d8.w1
        p2 = self.d8.p2   ;  w2 = self.d8.w2
        p3 = self.d8.p3   ;  w3 = self.d8.w3
        p4 = self.d8.p4   ;  w4 = self.d8.w4
        p5 = self.d8.p5   ;  w5 = self.d8.w5
        p6 = self.d8.p6   ;  w6 = self.d8.w6
        p7 = self.d8.p7   ;  w7 = self.d8.w7
        p8 = self.d8.p8   ;  w8 = self.d8.w8
        
        #-----------------------------------------
        # Add contributions from neighbor pixels
        #-----------------------------------------
        dt = self.dt          ############# CHECK dt #########
        if (np.size( self.da ) == 1):    
            factor = (dt / self.da)
            if (self.d8.p1_OK):    
                dzw[p1] += (self.Q_gw[w1] * factor)
            if (self.d8.p2_OK):    
                dzw[p2] += (self.Q_gw[w2] * factor)
            if (self.d8.p3_OK):    
                dzw[p3] += (self.Q_gw[w3] * factor)
            if (self.d8.p4_OK):    
                dzw[p4] += (self.Q_gw[w4] * factor)
            if (self.d8.p5_OK):    
                dzw[p5] += (self.Q_gw[w5] * factor)
            if (self.d8.p6_OK):    
                dzw[p6] += (self.Q_gw[w6] * factor)
            if (self.d8.p7_OK):    
                dzw[p7] += (self.Q_gw[w7] * factor)
            if (self.d8.p8_OK):    
                dzw[p8] += (self.Q_gw[w8] * factor)
        else:    
            if (self.d8.p1_OK):    
                dzw[p1] += (dt * self.Q_gw[w1] / self.da[p1])
            if (self.d8.p2_OK):    
                dzw[p2] += (dt * self.Q_gw[w2] / self.da[p2])
            if (self.d8.p3_OK):    
                dzw[p3] += (dt * self.Q_gw[w3] / self.da[p3])
            if (self.d8.p4_OK):    
                dzw[p4] += (dt * self.Q_gw[w4] / self.da[p4])
            if (self.d8.p5_OK):    
                dzw[p5] += (dt * self.Q_gw[w5] / self.da[p5])
            if (self.d8.p6_OK):    
                dzw[p6] += (dt * self.Q_gw[w6] / self.da[p6])
            if (self.d8.p7_OK):    
                dzw[p7] += (dt * self.Q_gw[w7] / self.da[p7])
            if (self.d8.p8_OK):    
                dzw[p8] += (dt * self.Q_gw[w8] / self.da[p8])
        
        #--------------------------------------------------
        # Find pixels where water table will rise or fall
        # Note: R = Rising, F = Falling
        #--------------------------------------------------
        wR = np.where( dzw > 0 )
        n_rising  = np.size( wR[0] )
        wF = np.where( dzw <= 0 )
        n_falling = np.size( wF[0] )
        print('   n_rising  = ' + str(n_rising))
        print('   n_falling = ' + str(n_falling))
        
        #-----------------------------------------
        # For debugging: save initial value of h     #*****************
        #-----------------------------------------
        start_h = self.h_table.copy()
        
        #-----------------------------------------
        # Process pixels where water table rises
        #-----------------------------------------
        if (n_rising > 0):    
            #---------------------------------
            # Must work from bottom layer up
            #----------------------------------------
            # Note:  xrange(start, stop, step), and
            #        last value is 0 if stop == -1.
            #----------------------------------------
            for k in range((self.n_layers - 1), -1, -1):
                
                #-------------------------------------
                # Compute unused capacity of layer k
                # Used before 7/13/06.
                #-------------------------------------
                #** yk = y[*,*,k]
                #** capacity_k = (t[k] - yk[wR]) * p[k]
                
                #-------------------------------------
                # Compute unused capacity of layer k
                #-------------------------------------
                yk = self.y[k,:,:]
                tk = self.th[k]    # (thickness of layer)
                pk = self.qs[k]    # ("porosity" of layer)
                SCALAR_THICKNESS = (np.size( tk ) == 1)
                SCALAR_POROSITY  = (np.size( pk ) == 1)

                # print 'k =', k
                # print 'IN SCALAR_THICKNESS test block...'
                ##############################################
                
                if (SCALAR_THICKNESS):    
                    if (SCALAR_POROSITY):    
                        capacity_k = (tk - yk[wR]) * pk
                    else:    
                        capacity_k = (tk - yk[wR]) * pk[wR]
                else:    
                    if (SCALAR_POROSITY):    
                        capacity_k = (tk[wR] - yk[wR]) * pk
                    else:    
                        capacity_k = (tk[wR] - yk[wR]) * pk[wR]
                
                w = np.where(capacity_k > 0)
                n_not_full = np.size( w[0] )

                ######################################################
                #  See Notes at top regarding "nested WHERE calls".
                ######################################################
            
                # print 'n_not_full =', n_not_full
                ######################################
                
                if (n_not_full != 0):
                    #-------------------------------------------
                    # Raise water table, update y, consume dzw
                    # Note: Both nwG and nwL may be nonzero.
                    #-------------------------------------------
                    wG  = np.where( dzw.flat[w] > capacity_k[w] )
                    nwG = np.size( wG[0] )
                    wL  = np.where( dzw.flat[w] <= capacity_k[w] )
                    nwL = np.size( wL[0] )
                    
##                    wG  = np.where( dzw[wR[w]] > capacity_k[w] )
##                    nwG = np.size( wG[0] )
##                    wL  = np.where( dzw[wR[w]] <= capacity_k[w] )
##                    nwL = np.size( wL[0] )
                    
##                    print 'nwG =', nwG  ############
##                    print 'nwL =', nwL  ############
                    
                    if (nwG != 0):

                        if (SCALAR_THICKNESS):    
                            self.h_table.flat[wG] += (tk - yk.flat[wG])
                            yk.flat[wG] = tk
                        else:    
                            self.h_table.flat[wG] += (tk.flat[wG] - yk.flat[wG])
                            yk.flat[wG] = tk.flat[wG]             
                        dzw.flat[wG] -= capacity_k[wG]
                        
##                        IDs = wR[w[wG]]         
##                        if (SCALAR_THICKNESS):    
##                            self.h_table[IDs] += (tk - yk[IDs])
##                            yk[IDs] = tk
##                        else:    
##                            self.h_table[IDs] += (tk[IDs] - yk[IDs])
##                            yk[IDs] = tk[IDs]             
##                        dzw[IDs] = dzw[IDs] - capacity_k[w[wG]]
                    
                    if (nwL != 0):
                        #--------------------------------
                        # Note: p[k]=0 => capacity_k=0,
                        # so okay to divide by p[k]
                        #--------------------------------
                        if (SCALAR_POROSITY):    
                            dh = dzw.flat[wL] / pk
                        else:    
                            dh = dzw.flat[wL] / pk.flat[wL]
                        self.h_table.flat[wL] += dh
                        yk.flat[wL] += dh                       
                        dzw.flat[wL] = np.float64(0)
                        
##                        IDs = wR[w[wL]]                        
##                        if (SCALAR_POROSITY):    
##                            dh = dzw[IDs] / pk
##                        else:    
##                            dh = dzw[IDs] / pk[IDs]
##                        self.h_table[IDs] += dh
##                        yk[IDs] += dh                       
##                        dzw[IDs] = np.float64(0)
                        
                    self.y[k,:,:] = yk   # (replace a layer in y)
                else:
                    pass
                    # dum = np.int16(0)
                    #kstr = str(k)
                    #print '   Layer ' + kstr + ' is full.'

                #--------------
                # For testing
                #--------------
                #hmin = self.h_table.min()
                #hmax = self.h_table.max()
                #print 'hmin = ' + str(hmin)
                #print 'hmax = ' + str(hmax)
            
            #------------------------------------------------
            # Where dzw is still gt 0, we must add it to h
            # since we have exhausted the capacity of the
            # soil layers.  This will bring h above the DEM
            # surface, z.  The increase in h will result in
            # surface runoff via a positive seepage rate.
            #------------------------------------------------
            # NB!  Not dzw.flat; wR is a tuple, see above.
            #------------------------------------------------            
            ## print 'IN dzw[wR] block...'  ###########
            w    = np.where( dzw[wR] > 0 )
            nIDs = np.size( w[0] )
            if (nIDs > 0):
                self.h_table.flat[w] += dzw.flat[w]
                dzw.flat[w] = np.float64(0)   # (need this)
                
##            w = np.where(dzw[wR] > np.float64(0))
##            nIDs = np.size(w[0])
##            if (nIDs > 0):    
##                IDs = wR[w]
##                self.h_table[IDs] += dzw[IDs]
##                dzw[IDs] = np.float64(0)   # (need this)
        
        #-----------------------------------------
        # Process pixels where water table falls
        #-----------------------------------------     
        if (n_falling > 0):    
            
            #--------------------------------
            # Must work from top layer down
            #--------------------------------
            for k in range(self.n_layers):
                #-----------------------------------
                # Compute water content of layer k
                #-----------------------------------
                yk = self.y[k,:,:]
                tk = self.th[k]    #  (thickness of layer)
                pk = self.qs[k]    # ("porosity" of layer)
                SCALAR_THICKNESS = (np.size( tk ) == 1)
                SCALAR_POROSITY  = (np.size( pk ) == 1)
                
                if (SCALAR_POROSITY):    
                    content_k = (yk[wF] * pk)
                else:    
                    content_k = (yk[wF] * pk[wF])        
                w = np.where( content_k > 0 )
                n_not_empty = np.size( w[0] )
                
                if (n_not_empty != 0):    
                    #-------------------------------------------
                    # Lower water table, update y, consume dzw
                    # Now dzw is negative or zero.
                    # Note: Both nwG and nwL may be nonzero.
                    #-------------------------------------------
                    wG  = np.where(np.absolute(dzw.flat[w])  > content_k[w])
                    nwG = np.size( wG[0] )
                    wL  = np.where(np.absolute(dzw.flat[w]) <= content_k[w])
                    nwL = np.size( wL[0] )
                    
##                    wG  = np.where(np.absolute(dzw[wF[w]])  > content_k[w])
##                    nwG = np.size( wG[0] )
##                    wL  = np.where(np.absolute(dzw[wF[w]]) <= content_k[w])
##                    nwL = np.size( wL[0] )
                    
                    #print 'nwG =', nwG 
                    #print 'nwL =', nwL
                    
                    if (nwG != 0):
                        self.h_table.flat[wG] -= yk.flat[wG]
                        yk.flat[wG]            = np.float64(0)
                        dzw.flat[wG]          += content_k[wG]  #(neg + pos)
                        
##                        IDs = wF[w[wG]]
##                        self.h_table[IDs] -= yk[IDs]
##                        yk[IDs]  = np.float64(0)
##                        dzw[IDs] += content_k[w[wG]]  #(neg + pos)
                    
                    if (nwL != 0):      #(3/8/04: BUG FIX)
                        #-------------------------------
                        # Note: p[k]=0 => content_k=0,
                        # so okay to divide by p[k]
                        #-------------------------------
                        if (SCALAR_POROSITY):    
                            dh = dzw.flat[wL] / pk
                        else:    
                            dh = dzw.flat[wL] / pk.flat[wL]
                        self.h_table.flat[wL] += dh
                        yk.flat[wL]           += dh
                        dzw.flat[wL]           = np.float64(0)
                        
##                        IDs = wF[w[wL]]
##                        if (SCALAR_POROSITY):    
##                            dh = dzw[IDs] / pk
##                        else:    
##                            dh = dzw[IDs] / pk[IDs]
##                        self.h_table[IDs] += dh
##                        yk[IDs]           += dh
##                        dzw[IDs]           = np.float64(0)
                        
                        #dz_min = dzw.min()
                        #dz_max = dzw.max()
                        #print 'dz_min = ',dz_min
                        #print 'dz_max = ',dz_max
                        #dh_min = self.h_table.min()
                        #dh_max = self.h_table.max()
                        #print 'pk[0]  =', pk[0]
                        #print 'dh_min =', dh_min
                        #print 'dh_max =', dh_max
                    
                    self.y[k,:,:] = yk    # (replace a layer in y)
                else:
                    pass
                    # dum = int16(0)
                    #print '   Layer ' + str(k) + ' is full.'
                
                #hmin = self.h_table.min()
                #hmax = self.h_table.max()
                #print '(hmin, hmax) =', hmin, hmax
            
            #------------------------------------------------
            # Where dzw is still lt 0, we must subtract it
            # from h;  all soil layers are now empty.  This
            # will bring h below the depth of the lowest
            # soil layer.  Should we assume that porosity
            # is the same as for the lowest layer or should
            # bottom of bottom layer be impermeable?
            #------------------------------------------------
            # This is where we should use depth to bedrock.
            #------------------------------------------------
            # NB!  Not dzw.flat; wF is a tuple, see above.
            #------------------------------------------------ 
            w    = np.where( dzw[wF] < 0 )
            nIDs = np.size( w[0] )
            
            if (nIDs > 0):             
                #--------------------------
                # pk = gp.qs[n_layers - 1]
                #--------------------------
                if (SCALAR_POROSITY):    
                    dh = dzw.flat[w] / pk
                else:    
                    dh = dzw.flat[w] / pk.flat[w] 
                self.h_table.flat[w] += dh
                dzw.flat[w] = np.float64(0)  # (need this ?)
                
##                IDs = wF[w]
##                if (SCALAR_POROSITY):    
##                    dh = dzw[IDs] / pk
##                else:    
##                    dh = dzw[IDs] / pk[IDs] 
##                self.h_table[IDs] += dh
##                dzw[IDs] = np.float64(0)  # (need this ?)
        
        #--------------------------------------------
        # For debugging: How much has h changed ?       ;****************
        # NB!  h may have NaNs if derived from DEM.
        # Exclude these from min & max via /NAN.
        #--------------------------------------------
        #dh2    = (h - start_h)
        #dh_min = dh2.min()
        #dh_max = dh2.max()  # (exclude NaNs as in IDL ??)
        #print,'(dh_min, dh_max) = ', dh_min, dh_max
        
        #---------------------------------
        # Option to print mins and maxes
        #--------------------------------- 
        #print ' '
        #print 'hmin =', self.h_table.min()
        #print 'hmax =', self.h_table.max()
        #------------------------------------
        #print 'ymin =', self.y.min()
        #print 'ymax =', self.y.max()
        
        #--------------
        # For debugging
        #--------------
        #if (np.size(h) == 1):
        #    msg = [' ', 'ERROR: h is not a 2D array.', ' ']
        #    result = GUI_Message(msg, /INFO)
        #    sys.exit()
     
    #   update_water_table()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        #------------------------------------------------------
        # This method uses prepend_directory() in BMI_base.py
        # which uses both eval and exec.
        #------------------------------------------------------
#         in_files = ['elev_file', 'h0_table_file', 'd_freeze_file',
#                     'd_thaw_file' ]
#         self.prepend_directory( in_files, INPUT=True )

        #------------------------------------------------------
        # This avoids eval/exec, but is brute-force
        # 2020-05-03. Changed in_directory to topo_directory.
        # See set_directories() in BMI_base.py.
        #------------------------------------------------------
        self.elev_file     = self.topo_directory + self.elev_file
        self.h0_table_file = self.soil_directory + self.h0_table_file
        ### self.d_bedrock_file = self.soil_directory + self.d_bedrock_file
        self.d_freeze_file = self.soil_directory + self.d_freeze_file
        self.d_thaw_file   = self.soil_directory + self.d_thaw_file

        #----------------------------------------------        
        # Open all input files and store file objects
        #----------------------------------------------
        self.elev_unit      = model_input.open_file(self.elev_type,      self.elev_file)
        self.h0_table_unit  = model_input.open_file(self.h0_table_type,  self.h0_table_file)
##        self.d_bedrock_unit = model_input.open_file(self.d_bedrock_type, self.d_bedrock_file)
        self.d_freeze_unit  = model_input.open_file(self.d_freeze_type,  self.d_freeze_file)
        self.d_thaw_unit    = model_input.open_file(self.d_thaw_type,    self.d_thaw_file)
        
        self.Ks_unit = []  # (empty lists to hold file objects)
        self.qs_unit = []
        self.th_unit = []
        
        for j in range(self.n_layers):
            self.Ks_file[j] = self.soil_directory + self.Ks_file[j]
            self.qs_file[j] = self.soil_directory + self.qs_file[j]
            self.th_file[j] = self.soil_directory + self.th_file[j]

            self.Ks_unit.append( model_input.open_file(self.Ks_type[j], self.Ks_file[j]) )
            self.qs_unit.append( model_input.open_file(self.qs_type[j], self.qs_file[j]) )
            self.th_unit.append( model_input.open_file(self.th_type[j], self.th_file[j]) )

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        rti = self.rti
      
        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        elev = model_input.read_next(self.elev_unit, self.elev_type, rti)
        if (elev is not None): self.elev = elev

        h0_table = model_input.read_next(self.h0_table_unit, self.h0_table_type, rti)
        if (h0_table is not None): self.h0_table = h0_table

        # This is not used yet. ####################
##        d_bedrock = model_input.read_next(self.d_bedrock_unit, self.d_bedrock_type, rti)
##        if (d_bedrock is not None): self.d_bedrock = d_bedrock

        #-----------------------------------------------------
        # These are computed variables, not input variables,
        # although we may read their initial value from file
        #-----------------------------------------------------
##        d_freeze = model_input.read_next(self.d_freeze_unit, self.d_freeze_type, rti)
##        if (d_freeze is not None): self.d_freeze = d_freeze
##
##        d_thaw = model_input.read_next(self.d_thaw_unit, self.d_thaw_type, rti)
##        if (d_thaw is not None): self.d_thaw = d_thaw
        
        #----------------------------------------------------
        # These are used by the new, more general GW method
        #----------------------------------------------------
        for j in range(self.n_layers):
            Ks = model_input.read_next(self.Ks_unit[j], self.Ks_type[j], rti)
            if (Ks is not None): self.Ks[j] = Ks

            qs = model_input.read_next(self.qs_unit[j], self.qs_type[j], rti)
            if (qs is not None): self.qs[j] = qs

            th = model_input.read_next(self.th_unit[j], self.th_type[j], rti)
            if (th is not None): self.th[j] = th            
        
    #   read_input_files()        
    #-------------------------------------------------------------------  
    def close_input_files(self):

        #############################################################
        # NOTE:  These first 3 files have single grids and should
        #        be closed right after they are read the 1st time.
        #############################################################
        #if (self.elev_file      != ''): self.elev_unit.close()        
        if (self.h0_table_file  != ''): self.h0_table_unit.close()
##        if (self.d_bedrock_file != ''): self.d_bedrock_unit.close()
        if (self.d_freeze_file  != ''): self.d_freeze_unit.close()
        if (self.d_thaw_file    != ''): self.d_thaw_unit.close()
        
        for j in range(self.n_layers):
            if (self.Ks_type[j] != 'Scalar'): self.Ks_unit[j].close()        
            if (self.qs_type[j] != 'Scalar'): self.qs_unit[j].close()
            if (self.th_type[j] != 'Scalar'): self.th_unit[j].close()
            #-----------------------------------------------------------
##            if (self.Ks_file[j] != ''): self.Ks_unit[j].close()        
##            if (self.qs_file[j] != ''): self.qs_unit[j].close()
##            if (self.th_file[j] != ''): self.th_unit[j].close()
        
    #   close_input_files()
    #-------------------------------------------------------------------
    def update_outfile_names(self):

        #----------------------
        # Filename extensions
        #----------------------
        # '_2D-htable.rts', '_2D-dfreeze.rts', '_2D-dthaw.rts'
        # '_0D-htable.txt', '_0D-dfreeze.txt', '_0D-dthaw.txt'
        
        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.ht_gs_file = (self.out_directory + self.ht_gs_file)
        self.df_gs_file = (self.out_directory + self.df_gs_file)
        self.dt_gs_file = (self.out_directory + self.dt_gs_file)
        #----------------------------------------------------------
        self.ht_ts_file = (self.out_directory + self.ht_ts_file)
        self.df_ts_file = (self.out_directory + self.df_ts_file)
        self.dt_ts_file = (self.out_directory + self.dt_ts_file)

    #   update_outfile_names()   
    #-------------------------------------------------------------------  
    def open_output_files(self):

        model_output.check_netcdf()
        self.update_outfile_names()
        
        #--------------------------------------
        # Open new files to write grid stacks
        #--------------------------------------
        if (self.SAVE_HT_GRIDS):
            model_output.open_new_gs_file( self, self.ht_gs_file, self.rti,
                                           var_name='ht',
                                           long_name='water_table_height',
                                           units_name='m')
            
        if (self.SAVE_DF_GRIDS):
            model_output.open_new_gs_file( self, self.df_gs_file, self.rti,
                                           var_name='df',
                                           long_name='depth_of_frozen_soil',
                                           units_name='m')
            
        if (self.SAVE_DT_GRIDS):
            model_output.open_new_gs_file( self, self.dt_gs_file, self.rti,
                                           var_name='dt',
                                           long_name='depth_of_thaw',
                                           units_name='m')
            
        #--------------------------------------
        # Open new files to write time series
        #--------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_HT_PIXELS):
            model_output.open_new_ts_file( self, self.ht_ts_file, IDs,
                                           var_name='ht',
                                           long_name='water_table_height',
                                           units_name='m')

        if (self.SAVE_DF_PIXELS):
            model_output.open_new_ts_file( self, self.df_ts_file, IDs,
                                           var_name='df',
                                           long_name='depth_of_frozen_soil',
                                           units_name='m')

        if (self.SAVE_DT_PIXELS):
            model_output.open_new_ts_file( self, self.dt_ts_file, IDs,
                                           var_name='dt',
                                           long_name='depth_of_thaw',
                                           units_name='m')

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
    
        if (self.SAVE_HT_GRIDS):  model_output.close_gs_file( self, 'ht')   
        if (self.SAVE_DF_GRIDS):  model_output.close_gs_file( self, 'df') 
        if (self.SAVE_DT_GRIDS):  model_output.close_gs_file( self, 'dt') 
        #-----------------------------------------------------------------
        if (self.SAVE_HT_PIXELS): model_output.close_ts_file( self, 'ht')  
        if (self.SAVE_DF_PIXELS): model_output.close_ts_file( self, 'df')
        if (self.SAVE_DT_PIXELS): model_output.close_ts_file( self, 'dt')

    #   close_output_files()        
    #-------------------------------------------------------------------  
    def save_grids(self):
        
        #------------------------------------------------------
        # Notes:  The RTG function will work whether argument
        #         is a scalar or already a 2D grid.
        #------------------------------------------------------
        if (self.SAVE_HT_GRIDS):
            model_output.add_grid( self, self.h_table, 'ht', self.time_min )

        if (self.SAVE_DF_GRIDS):
            model_output.add_grid( self, self.d_freeze, 'df', self.time_min )

        if (self.SAVE_DT_GRIDS):
            model_output.add_grid( self, self.d_thaw, 'dt', self.time_min )
            
    #   save_grids()            
    #-------------------------------------------------------------------  
    def save_pixel_values(self):

        IDs  = self.outlet_IDs
        time = self.time_min    ######
        
        if (self.SAVE_HT_PIXELS):
            model_output.add_values_at_IDs( self, time,
                                            self.h_table, 'ht', IDs )

        if (self.SAVE_DF_PIXELS):
            model_output.add_values_at_IDs( self, time,
                                            self.d_freeze, 'df', IDs )

        if (self.SAVE_DT_PIXELS):
            model_output.add_values_at_IDs( self, time,
                                            self.d_thaw, 'dt', IDs )
            
    #   save_pixel_values()
    #-------------------------------------------------------------------
    def get_total_darcy_layer_flow(self, cv):

        #---------------------------------------------------------
        # Notes: gv = gw_vars = structure
        #        z  = elevation of land surface [m]
        #        h  = elevation of water table [m]
        #        (z-h) = depth to water table [m]
        #        Sh = water table slope [unitless]
        #        Ks = soil_K = sat. hydraulic conductivity [m/s]
        #        dw = element width [m]
        #        ds = hor. Dist. between pixel and parent [m]
        #        y  = wetted flow depth in each layer [m]
        #             (could be a recycled local variable)
        #        Qg = total Darcy-flow discharge [m^3/s]
        #             (summed over all layers)
        #        diff = (partial sum of soil thicknesses -
        #                depth to water table)

        #        See Notes for Update_Water_Table routine.
        #---------------------------------------------------------
        
        #-------------------------------------------
        # Compute water table slope from h_table
        #-------------------------------------------
        # NB!  h_table is assumed to be a grid.
        #-------------------------------------------
        # NB!  Q is zero where Sh is zero.
        #-------------------------------------------
        # NB!  Flow direction is still assumed to
        #      be given by the DEM's D8 flow grid.
        #-------------------------------------------
        # Sh = Free_Surface_Slope(np.float64(0), self.h_table, self.d8.ds, self.d8.pIDs)

        #-----------------------------------------------------------
        # Notes:  It is assumed that the flow directions don't
        #         change even though the free surface is changing.
        #-----------------------------------------------------------
        delta_h = (self.h_table - self.h_table[self.d8.parent_IDs])
        Sh      = (delta_h / self.d8.ds)
    
        #------------------------------------------
        # Compute wetted-depth, y, for each layer
        # Now passed by caller.
        #------------------------------------------
        #** diff = -(z - h)
        
        #----------------------------------
        # NB!  h is assumed to be a grid.
        # It is converted, if necessary,
        # by Initialize_GW_Vars.
        #----------------------------------
        nx = self.nx
        ny = self.ny
        self.Q_gw = zeros([ny, nx], dtype='Float64')
        
        for k in range(self.n_layers):
            #--------------------------------------
            # Add Q for this layer, via Darcy law
            #--------------------------------------
            self.Q_gw += (self.Ks[k] * Sh * self.d8.dw * self.y[k,:,:])   
        
    #   get_total_darcy_layer_flow
    #-------------------------------------------------------------------


#-------------------------------------------------------------------
#   Outside of the class    
#-------------------------------------------------------------------
def Darcy_Flow(K, S, dw, y):

    #-----------------------------------------------------
    # Notes: K  = hydraulic conductivity [m/s]
    #        S  = water table slope (perhaps topo slope)
    #        dw = width of element [m]
    #        y  = wetted flow depth in a layer [m]
    #        Q  = discharge through a layer [m^3 / s]
    #-----------------------------------------------------
    Q = K * S * dw * y
    
    return Q

#   Darcy_Flow   
#-------------------------------------------------------------------        
