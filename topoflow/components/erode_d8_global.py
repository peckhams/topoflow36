
## Copyright (c) 2001-2013, Scott D. Peckham
##
## January 2013   (Revised handling of input/output names).
## January 2009  (converted from IDL)
## September, October, November 2009
## February, March 2010  (separate from local time-stepping)
## November 2010  (added update_dt_grid_method1().)
## October 2011 (changes to finalize())

#-------------------------------------------------------------
## NB!  LINK_FLATS is set to False in initialize_d8_vars().
##
#-----------------------------------------------------------------------
#
#  class erosion_component    # (inherits from erode_base.py)
#
#      get_component_name()
#      get_attribute()            ## (10/27/11)
#      get_input_var_names()
#      get_output_var_names()
#      get_var_name()
#      get_var_units()
#      get_var_type()             ## (not ready yet)
#-------------------------------
#      update()
#      finalize()
#-------------------------------
#      initialize_d8_vars()
#      update_d8_vars()
#      update_slope_grid()
#      update_Q_grid()
#      update_Qs_grid()
#      update_dz_dt_grid()         ## (2/19/10)
#-------------------------------
#      update_dt_grid()
#      update_dt_grid_method0()    ## (Uses Qs[k] - Qs[p])
#      update_dt_grid_method1()    ## (Uses Qs[k], 11/9/10)
#-------------------------------
#      update_DEM()                ## (2/22/10)
#-------------------------------
#      save_final_grids()          ## (10/17/11)
#
#-------------------------------
#      update_dt_grid_method2()    ## (in erode_base.py)
#      update_min_dz_up_grid()     ## (in erode_base.py)
#
#-----------------------------------------------------------------------

import numpy as np
import os.path
import time
#----------------------
# Could also do this.
#----------------------
# from numpy import where, logical_and, logical_or

from . import d8_global   ### (9/19/14.  Attempt to fix issue.)
from . import erode_base
from topoflow.utils import rtg_files

# from topoflow.components import d8_global
# from topoflow.components import erode_base
# from topoflow.utils      import rtg_files

#-----------------------------------------------------------------------
class erosion_component( erode_base.erosion_component ):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'Erode_D8_Global',
        'version':            '0.5',
        'author_name':        'Scott Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'adaptive',
        'step_method':        'explicit',
        #------------------------------------------------------        
        'comp_name':          'ErodeGlobal',
        'model_family':       'Erode',
        'cfg_template_file':  'Erode_Global.cfg.in',
        'cfg_extension':      '_erode_global.cfg',
        'cmt_var_prefix':     '/ErodeGlobal/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/erode/0.5/src/share/cmt/gui/Erode_Global.xml',
        'dialog_title':       'LEM: Erode Global-Timestep Parameters',
        'time_units':         'years' }

    #-------------------------------------------------------------
    # Later on, we may want to get initial DEM (z0), uplift rate
    # (U) and geomorphic rain rate (R) from another component.
    #-------------------------------------------------------------
    _input_var_names = []
#         'atmosphere_water__geomorphic_precipitation_leq-volume_flux', 
#         'bedrock_uplift_rate',
#         'land_surface__initial_elevation']

    #-------------------------------------------------------------
    # Later, maybe add: nx, ny, codes, time_sec, time_min.
    #-------------------------------------------------------------
    _output_var_names = [
        'atmosphere_water__geomorphic_precipitation_leq-volume_flux',  # R
        'bedrock__uplift_rate',                                        # U
        'channel_water_x-section__volume_flow_rate',                   # Q
        'channel_water_x-section__volume_flow_rate_law_area_exponent', # p
        ##'channel_water_x-section__volume_flow_rate_law_coefficient', # R
        'channel_water_total-sediment__volume_flow_rate',                    # Qs
        'channel_water_total-sediment__volume_flow_rate_law_area_exponent',  # m
        'channel_water_total-sediment__volume_flow_rate_law_coefficient',    # K
        'channel_water_total-sediment__volume_flow_rate_law_slope_exponent', # n
        'land_surface__elevation',               # z
        'land_surface__increment_of_elevation',  # dz
        'land_surface__initial_elevation',       # z0
        'land_surface__domain_max_of_increment_of_elevation',     # dz_max
        'land_surface__slope',                                    # S
        'land_surface__time_derivative_of_elevation',             # dz_dt
        'land_surface_contour-segment__total_contributing_area',  # A
        'model__time_step',                 # dt
        'model_domain_boundary__lowering_rate', # BLR
        'model_grid_cell__area',            # da
        'model_grid_cell__d8_flow_width',   # dw
        'model_grid_cell__d8_flow_length',  # ds
        'model_grid_cell__diameter',        # dd
        'model_grid_cell__x_length',        # dx
        'model_grid_cell__y_length' ]       # dy

    #-------------------------------------------------------
    # Qs = annual sed. discharge [m^3/ yr]
    # Q  = annual discharge [m^3 / yr]
    # A  = contributing area [meters^2]
    # S  = slope [unitless]
    # R  = geomorphically effective rainrate [meters / yr]
    #         (unless p ne 0, then R = coefficient)
    #-------------------------------------------------------    
    # K = coefficient [(m^3 / yr)^(1 - m)]
    # m = discharge exponent (usually in [1,2])
    # n = slope exponent (usually in [1,2])
    #-------------------------------------------------------
    _var_name_map = {
        'atmosphere_water__geomorphic_precipitation_leq-volume_flux':  'R',
        'bedrock__uplift_rate':                                        'U',
        'channel_water_x-section__volume_flow_rate':                   'Q',
        'channel_water_x-section__volume_flow_rate_law_area_exponent': 'p',
        ##'channel_water_x-section__volume_flow_rate_law_coefficient': 'R',
        'channel_water_total-sediment__volume_flow_rate':                    'Qs',
        'channel_water_total-sediment__volume_flow_rate_law_area_exponent':  'm',
        'channel_water_total-sediment__volume_flow_rate_law_coefficient':    'K',
        'channel_water_total-sediment__volume_flow_rate_law_slope_exponent': 'n',
        'land_surface__elevation':                               'z',
        'land_surface__increment_of_elevation':                  'dz',
        'land_surface__initial_elevation':                       'z0',
        'land_surface__domain_max_of_increment_of_elevation':    'dz_max',
        'land_surface__slope':                                   'S',
        'land_surface__time_derivative_of_elevation':            'dz_dt',
        'land_surface_contour-segment__total_contributing_area': 'A',
        'model__time_step':                 'dt',
        'model_domain_boundary__lowering_rate':                  'BLR',
        'model_grid_cell__area':            'da',
        'model_grid_cell__d8_flow_width':   'dw',
        'model_grid_cell__d8_flow_length':  'ds',
        'model_grid_cell__diameter':        'dd',
        'model_grid_cell__x_length':        'dx',
        'model_grid_cell__y_length':        'dy' }

    #---------------------------------------------------------
    # Note that the units of "K" depend on the exponent "m".
    # How do we handle this in a general way ?? (2/6/13)
    #---------------------------------------------------------
    # self.U has units of [mm/yr], but
    # self.U_mpyr has units of [m/yr].
    # self.BLR_mpyr has units of [m/yr].
    #---------------------------------------------------------
    # Maybe:  geomorphic -> domain_long-time_average_of
    #---------------------------------------------------------    
    _var_units_map = {
        'atmosphere_water__geomorphic_precipitation_leq-volume_flux':  'm yr-1',
        'bedrock__uplift_rate':                                        'mm yr-1',
        'channel_water_x-section__volume_flow_rate':                   'm3 yr-1',
        'channel_water_x-section__volume_flow_rate_law_area_exponent': '1',
        ##'channel_water_x-section__volume_flow_rate_law_coefficient': 'm(1-p) yr(p-1)',        
        'channel_water_total-sediment__volume_flow_rate':              'm3 yr-1',
        'channel_water_total-sediment__volume_flow_rate_law_area_exponent':  '1',
        'channel_water_total-sediment__volume_flow_rate_law_coefficient':    'm3(1-m) yr(m-1)',
        'channel_water_total-sediment__volume_flow_rate_law_slope_exponent': '1',
        'land_surface__domain_max_of_increment_of_elevation': 'm',
        'land_surface__elevation':                  'm',
        'land_surface__increment_of_elevation':     'm',
        'land_surface__initial_elevation':          'm',
        'land_surface__slope':                      '1',
        'land_surface__time_derivative_of_elevation': 'm yr-1',
        'land_surface_contour-segment__total_contributing_area':   'm2',
        'model__time_step':                 'yr', 
        'model_domain_boundary__lowering_rate': 'mm yr-1',
        'model_grid_cell__area':            'm2',
        'model_grid_cell__d8_flow_width':   'm',
        'model_grid_cell__d8_flow_length':  'm',
        'model_grid_cell__diameter':        'm',
        'model_grid_cell__x_length':        'm',
        'model_grid_cell__y_length':        'm' }
        
    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Erode_D8_Global'

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
    def update(self, dt=-1.0):

##        if (SILENT == None): SILENT=self.SILENT
##        if (REPORT == None): REPORT=self.REPORT

        SILENT = self.SILENT
        REPORT = self.REPORT
        
        #### if not(SILENT) and (self.time_index == 0):
        if (self.time_index == 0):
            print('Erosion component: Processing...')
        self.status = 'updating'  # (OpenMI 2.0 convention)

        #---------------------------------
        # Print dz_max to track progress
        #---------------------------------
        if (self.mode == 'driver'):
            self.print_time_and_value(self.dz_max, 'dz_max', '[m]',
                                      interval=5.0 ) ### PRINT_INDEX=True)
            
        #--------------------------------------------
        # Write DEM, etc. before anything happens
        # Added this here (like local) on 11/15/11.
        # Moved down to fix bug on 2/22/13.
        #--------------------------------------------
##        if (self.time_index == 0):
##            self.write_output_files(0)
            
        #--------------------------------------
        # Update values from other components
        #--------------------------------------
        self.update_R()
        self.update_R_integral()
        self.update_U()
        self.update_U_integral()
        
        #--------------------------------------------------------
        # (11/8/10)  As written, update_base_level() in
        # erode_base.py makes all base_IDs have same elevation.
        #--------------------------------------------------------
        # update_DEM_edge_values() is only for improved use of
        # color in plots with certain boundary conditions.
        #--------------------------------------------------------
        ## self.update_base_level()
        self.update_DEM_edge_values()
        if (self.FILL_PITS):
            self.fill_pits_in_DEM( SILENT=SILENT )

        #--------------------------------------------
        # Update the D8 flow grid and all vars that
        # depend on it, including D8 area grid.
        #--------------------------------------------
        self.update_d8_vars( SILENT=SILENT, REPORT=REPORT )
        self.update_slope_grid( SILENT=SILENT, REPORT=REPORT )
        self.update_Q_grid( SILENT=SILENT, REPORT=REPORT)
        self.update_Qs_grid( SILENT=SILENT, REPORT=REPORT )
        self.update_dz_dt_grid( SILENT=SILENT, REPORT=REPORT )

        #---------------------------------------------
        # One version of update_dt_grid() will call
        # a function called update_min_dz_up_grid().
        #---------------------------------------------
        self.update_dt_grid( SILENT=SILENT, REPORT=REPORT )

        #---------------------------------------------
        # For testing;  comment out later. (4/14/10)
        # Save stack of n_grids for viewing.
        #---------------------------------------------
        # self.update_n_grid(SILENT=SILENT, REPORT=REPORT)   ###########

##        #----------------------------------------------
##        # Write initial DEM, etc. before update_DEM() 
##        #----------------------------------------------
##        # Moved from top to here on 2/22/13 to fix 
##        # bug, but must happen before update_DEM(). 
##        #--------------------------------------------
##        if (self.time_index == 0):
##            self.write_output_files(0)

        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # (2/22/13) Moved up to just before update_DEM()
        # so that all values correspond to the same time.
        # This also causes writing of initial DEM, etc.
        #----------------------------------------------
        # self.write_output_files(0)  # (force writing every step)
        # self.write_output_files( time )
        self.write_output_files()
        
        self.update_DEM( SILENT=SILENT, REPORT=REPORT )
        self.update_DEM_min_and_max()

##        print 'min(A),     max(A)     =', self.d8.A.min(), self.d8.A.max()
##        print 'min(S),     max(S)     =', self.S.min(), self.S.max()
##        print 'min(Q),     max(Q)     =', self.Q.min(), self.Q.max()
##        print 'min(Qs),    max(Qs)    =', self.Qs.min(), self.Qs.max()
##        print 'min(dz/dt), max(dz/dt) =', self.dz_dt.min(), self.dz_dt.max()
##        print 'min(z),     max(z)     =', self.DEM.min(), self.DEM.max()
##        print 'min(dz),    max(dz)    =', self.dz.min(), self.dz.max()
##        print 'min(dt),    max(dt)    =', self.dt_grid.min(), self.dt_grid.max()
        
        if not(SILENT):
            print(' ')
            
        #----------------------------
        # Print some mins and maxes
        #----------------------------
        step = self.time_index
        if ((step % 1000) == 0):
            # self.print_mins_and_maxes(step, DEM=True)
            self.print_mins_and_maxes(step, DZ_DT=True)

    
        ########################################
        # CAN THE UPDATED DEM HAVE PITS ??
        # IF NOT, DON'T CALL FILL_PITS.
        ########################################
        
        #------------------------
        # Check computed values
        #------------------------
        OK = self.check_stability()
        if (OK):
            self.status = 'updated'  # (OpenMI 2.0 convention)
        else:
            self.status = 'failed'
            self.DONE   = True

        #-------------------------------------------
        # Read from files as needed to update vars 
        #-----------------------------------------------------
        # NB! This is currently not needed for the "erosion
        # process" because values don't change over time and
        # read_input_files() is called by initialize().
        #-----------------------------------------------------
        # if (self.time_index > 0):
        #     self.read_input_files()

##        #----------------------------------------------
##        # Write user-specified data to output files ?
##        #----------------------------------------------
##        # self.write_output_files(0)  # (force writing every step)
##        # self.write_output_files( time )
##        self.write_output_files()

        #------------------------
        # Update internal clock
        #------------------------
        self.update_time( dt )
        ## print 'time_index =', self.time_index
        
        #-------------------------------------------
        # Check for steady-state condition instead
        #-------------------------------------------
        self.check_finished()   ######################

        #-----------------------------------------
        # Print some info before final report?
        # Final report is printed by finalize().
        #-----------------------------------------
        ## if (self.DONE):
            ## self.write_output_files(0)
            ###########################################
        
    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        #--------------------------------------------------------------
        # Note: This overrides finalize() in erode_base.py. (10/4/11)
        #--------------------------------------------------------------
        self.status = 'finalizing'  # (OpenMI)   
        self.close_input_files()    # Input "data streams"
        self.close_output_files()

        if (self.mode == 'driver'):
            self.print_time_and_value(self.dz_max, 'dz_max', '[m]',
                                      interval=0.0 ) ### PRINT_INDEX=True)

        #----------------------------------------
        # Save all other final grids (10/17/11)
        #----------------------------------------
        self.save_final_grids()
        
        #----------------------------------------------
        # Added on 10/4/11 to help gauge performance.
        # See CSDMS_base.print_final_report.
        # This could eventually go there because it
        # would be useful to any model.  But not all
        # of them have self.time in "years".
        #----------------------------------------------
        finish         = time.time()
        run_time_secs  = (finish - self.start_time)
        sim_time_years = self.time
        years_per_sec  = (sim_time_years / run_time_secs)
        print(' ')
        print('Years per second =', years_per_sec)
        
        #----------------------------------------------
        # Print report info common with Global version
        #-----------------------------------------------
        print(' ')
        print('CFL_factor =', self.CFL_factor)
        print(' ')
        print('min(dt), max(dt) =', self.dt_grid.min(), self.dt_grid.max())
        print(' ')
        print('min(z), max(z)   =', self.DEM.min(), self.DEM.max())
        print(' ')
        print('dz_max_vec[0]    =', self.dz_max_vec[0])
        print('dz_max_vec[nt-1] =', self.dz_max_vec[-1])
        print('dz_max_vec.min() =', self.dz_max_vec.min())
        print('dz_max_vec.max() =', self.dz_max_vec.max())
        print(' ')
        self.print_final_report(comp_name='Erode-D8-Global 3.1 (2/7/13)',
                                mode='driver')  ## NEED THIS !!
        
        self.status = 'finalized'  # (OpenMI)
          
    #   finalize()
    #-------------------------------------------------------------------
    def initialize_d8_vars(self):
        
        #---------------------------------------------
        # Compute and store a variety of (static) D8
        # flow grid variables.  Embed structure into
        # the "erosion_base" component.
        #---------------------------------------------
        self.d8 = d8_global.d8_component()   ### NOTE: GLOBAL HERE ###

        self.d8.DEBUG  = False  ###############
        self.d8.SILENT = True
        self.d8.REPORT = False
      
        #----------------------------------------------------
        # Copy all of this from the "host" component, just
        # in case the D8 component doesn't have a CFG file.
        # If it does, these variables will get overwritten.
        # Something similar is done in GW_base.py.
        #----------------------------------------------------
        # (1/23/12) Note that d8_base.py now has a new
        # method called: set_default_config_vars()
        # that is used to intialize vars in cases
        # where there is no "*_d8_global.cfg" file.
        # It is called in d8_base.initialize().
        # However, it currently sets:
        #    self.in_directory  = '~/CMT_Output/'  AND
        #    self.out_directory = '~/CMT_Output/'
        # which may be worse defaults than these.
        #-------------------------------------------------------------
        # (2/11/2017) The initialize() method in d8_base.py now
        # uses case_prefix (vs. site_prefix) for its CFG file:
        # <site_prefix>_d8_global.cfg.  This is to prevent confusion
        # since this was the only CFG file that used site_prefix.
        #-------------------------------------------------------------  
        self.d8.site_prefix   = self.site_prefix
        self.d8.case_prefix   = self.case_prefix
        self.d8.in_directory  = self.in_directory
        self.d8.out_directory = self.out_directory
        #-------------------------------------------
        self.d8.FILL_PITS_IN_Z0 = 0                   # (1/23/12)
        self.d8.A_units         = 'm^2'               # (1/23/12) May be needed.

        #---------------------------------------------         
        # D8 component builds its cfg filename from
        # in_directory, site_prefix and extension. 
        #---------------------------------------------     
        self.d8.initialize( cfg_file=None,
                            SILENT=self.SILENT,
                            REPORT=self.REPORT )
##                            SILENT=not(self.DEBUG),
##                            REPORT=self.DEBUG )

        #---------------------------------------------------
        # This overrides settings from D8_Local CFG file.
        # We don't need to link flats in the Erode models
        # because pits are filled dynamically.  Note that
        # S=0 for flats, so outfluxes Q and Qs will also
        # be zero and deposition will therefore occur.
        #---------------------------------------------------      
        self.d8.LINK_FLATS = False
        self.d8.BREAK_TIES = True
        print('In erode_d8_global.initialize_d8_vars():')
        print('   LINK_FLATS =', self.d8.LINK_FLATS)
        print('   BREAK_TIES =', self.d8.BREAK_TIES)

    #   initialize_d8_vars()
    #-------------------------------------------------------------------
    # def update_d8_vars(self):
    #
    #     (Inherited from erode_base.py)
    #
    #     update_d8_vars()
    #-------------------------------------------------------------------
    def update_slope_grid(self, SILENT=True, REPORT=False):

        #----------------------------------------------------
        # Notes: Make sure that d8 component is initialized
        #        with the same directory, site_prefix, etc.
        #        as this component.  Otherwise, the "shape"
        #        of DEM and ds, A, etc. won't match.
        #----------------------------------------------------
        if not(SILENT):    
            print('Updating slope grid...')
   
        #--------------------------------------------------
        # Compute slope (rise/run) toward D8 parent pixel
        #--------------------------------------------------
        # pIDs gives indices in the "np.where" style
        #--------------------------------------------------        
        pIDs   = self.d8.parent_IDs
        self.S = ((self.DEM - self.DEM[pIDs]) / self.d8.ds)

##        w = np.where(self.DEM < self.DEM[pIDs])
##        n_bad = w[0].size
##        if (n_bad != 0):
##            print '   Number of uphill parents =', n_bad
            
        #---------------------------------------------------
        # Set slope to zero wherever the D8 flow direction
        # is undefined.  This makes Qs(Q,S) = 0.
        #---------------------------------------------------
        ## w  = np.where(self.d8.parent_ID_grid == 0)
        w  = np.where(self.d8.d8_grid == 0)
        nw = w[0].size
        if (nw != 0):    
            self.S[w] = 0

        #-------------------------
        # Check for bad S values
        #-------------------------
        w2  = np.where(self.S < 0)
        nw2 = w2[0].size
        if (nw2 != 0):    
            # self.S[w2] = 1e-8
            print('   Negative values found in slope grid.')
            print('   Found ' + str(nw2) + ' negative values.')
            
            #---------------------------
            # Return the cols and rows
            #---------------------------
            rows = w2[0]
            cols = w2[1]
            print('   cols =')
            print(cols)
            print('   rows =')
            print(rows)
            sys.exit()   #####################
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            ## S_str = str(self.S.min())  + ', '  + str(self.S.max())
            S_str = str(np.nanmin(self.S)) + ', ' + str(np.nanmax(self.S))
            print('    min(S), max(S) = ' + S_str + ' [m/m]')

    #   update_slope_grid()
    #-------------------------------------------------------------------
    def update_Q_grid(self, SILENT=True, REPORT=False):

        #--------------------------------------------------------------
        # NOTES:  Q = annual discharge [m^3 / yr] (*leaving* a cell)
        #         R = geomorphically effective rainrate [meters / yr]
        #             (unless p ne 0, then R = coefficient)
        #         A = contributing area [meters^2]
        #--------------------------------------------------------------
        # Note: Q is zero wherever A is zero, and A will be zero
        #       wherever the D8 flow code is zero.
        #--------------------------------------------------------------        
        if not(SILENT):    
            print('Updating discharge grid...')
        
        if (self.p != 1):    
            self.Q = self.R * (self.d8.A ** self.p)
        else:    
            #----------------------
            # Don't convert units
            #----------------------
            #** R2 = R / self.secs_per_year   ;[m/yr -> m/s]
            self.Q = self.R * self.d8.A

        #--------------------------------------------------------
        # Test whether there are any pixels where flow code is
        # zero but area grid isn't.  This would be bad because
        # flow code is zero in pits, and A > 0 would lead to a
        # nonzero discharge (water and sed.) *out* of the cell.
        #--------------------------------------------------------
        # Tested for 100 timesteps on 2/22/10.
        #--------------------------------------------------------
##        w = np.where(np.logical_and(self.d8.d8_grid == 0, self.d8.A != 0))
##        nw = w[0].size
##        if (nw != 0):
##            print 'WARNING: There are places where flow grid is'
##            print '         zero and area grid is nonzero.'
##            print ' '
            
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            ## Q_str = str(self.Q.min())  + ', '  + str(self.Q.max())
            Q_str = str(np.nanmin(self.Q)) + ', ' + str(np.nanmax(self.Q))
            print('    min(Q), max(Q) = ' + Q_str + ' [m^3/yr]')

    #   update_Q_grid()
    #-------------------------------------------------------------------
    def update_Qs_grid(self, SILENT=True, REPORT=False):
     
        #--------------------------------------------------------
        #NOTES:  Qs = annual sed. discharge [m^3/ yr]
        #        Q  = annual discharge [m^3 / yr]
        #        S  = slope [unitless]
        #        K  = coefficient [(m^3 / yr)^(1 - m)]
        #        m  = discharge exponent (usually in [1,2])
        #        n  = slope exponent (usually in [1,2])

        #        The standard formula Qs = K * Q^m * S^n is for
        #        the case where Q and Qs have units of m^3/sec.
        #        The extra factor below makes the formula valid
        #        for the case where units of both are m^3/yr.
        #        That is, Qs' = fac * K * (Q')^m * S^n.

        #        The Slope_Grid function checks for negative
        #        slopes and either adjusts them or aborts.
        #--------------------------------------------------------
        if not(SILENT):
            print('Updating sed. discharge grid...')
        
        fac1    = self.secs_per_year ** (np.float64(1) - self.m)
        fac2    = self.K * (self.Q ** self.m)
        fac3    = (self.S ** self.n)
        self.Qs = fac1 * fac2 * fac3     #[m^3 / year]
        
        #---------------------------------------
        # Impose no-flux boundary condition
        # on all four edges of DEM ?
        # NB!  If S=0 on edges, then Qs=0.
        #---------------------------------------
        # Is it better to put a "wall" around
        # the outside to impose this B.C. ?
        #---------------------------------------
        # nx = self.nx
        # ny = self.ny
        # self.Qs[0, :]    = 0.0
        # self.Qs[nx-1, :] = 0.0
        # self.Qs[:, 0]    = 0.0
        ### self.Qs[:, ny-1] = 0.0    ;(Let sediment get out)
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            ## Qs_str = str(self.Qs.min())  + ', '  + str(self.Qs.max())
            Qs_str = str(np.nanmin(self.Qs)) + ', ' + str(np.nanmax(self.Qs))
            print('    min(Qs), max(Qs) = ' + Qs_str + ' [m^3/yr]')

    #   update_Qs_grid()
    #-------------------------------------------------------------------
    def update_dz_dt_grid(self, SILENT=True, REPORT=False):
        
        #-------------------------------------------------------------
        # Notes: dV_dt is the net rate at which the volume in a
        #        grid cell is changing.  It is computed as the
        #        difference between all of the fluxes into the grid
        #        cell and flux out of the grid cell. (D8-based here)
        #        dz_dt = (dV_dt / da)  [m/yr].
        #
        #        (dz_dt > 0) => deposition at that cell
        #        (dz_dt < 0) => erosion at that cell
        #
        #        Qs    = K * Q^m * S^n = sed. discharge [m^3/yr]
        #        da    = pixel area grid  [m^2]
        #-------------------------------------------------------------        
        if not(SILENT):    
            print('Updating dz/dt grid...')

        #----------------------------------------------
        # Initialize dz_dt with the sediment flux out
        # of each pixel divided by the pixel area.
        # Pixel area, self.da, may be scalar or grid.
        # Should be 0 for pixels with flow code of 0.
        #----------------------------------------------
        # Note that in timestep, dt,  dz = dz_dt * dt
        #----------------------------------------------
        flux_out   = self.Qs / self.da
        self.dz_dt = -flux_out   # [m^3 / year]

        #-----------------------------------------
        # Add the contribution to to uplift, U.
        # self.U has units of [mm/yr], but
        # self.U_mpyr has units of [m/yr].
        #-----------------------------------------
        self.dz_dt += self.U_mpyr
        
        #-----------------------------------------
        # Add contributions from neighbor pixels
        #-----------------------------------------
        if (self.d8.n1 != 0):    
            self.dz_dt[ self.d8.p1 ] += flux_out[ self.d8.w1 ]
        if (self.d8.n2 != 0):    
            self.dz_dt[ self.d8.p2 ] += flux_out[ self.d8.w2 ]
        if (self.d8.n3 != 0):    
            self.dz_dt[ self.d8.p3 ] += flux_out[ self.d8.w3 ]
        if (self.d8.n4 != 0):    
            self.dz_dt[ self.d8.p4 ] += flux_out[ self.d8.w4 ]
        if (self.d8.n5 != 0):    
            self.dz_dt[ self.d8.p5 ] += flux_out[ self.d8.w5 ]
        if (self.d8.n6 != 0):    
            self.dz_dt[ self.d8.p6 ] += flux_out[ self.d8.w6 ]
        if (self.d8.n7 != 0):    
            self.dz_dt[ self.d8.p7 ] += flux_out[ self.d8.w7 ]
        if (self.d8.n8 != 0):    
            self.dz_dt[ self.d8.p8 ] += flux_out[ self.d8.w8 ]
     
        #----------------------------------------------
        # This is only necessary for "completeness"
        # since elevations at base_IDs are set by the
        # update_base_level() function.
        #----------------------------------------------        
        # Use Base-level Lowering Rate (BLR) to set
        # dz_dt on edges.
        # self.BLR_mpyr has units of [m/yr].
        #----------------------------------------------
        # base_IDs are defined in the function:
        # erode_base.initialize_boundary_conditions()
        #----------------------------------------------        
        self.dz_dt.flat[ self.base_IDs ] = -self.BLR_mpyr

        #------------------
        # Optional report
        #------------------
        ## REPORT = True  #########
        if (REPORT):
            dz_dt_min = self.dz_dt.min()
            dz_dt_max = self.dz_dt.max()
            print('(dz_dt_min, dz_dt_max) =', dz_dt_min, dz_dt_max)
            print(' ')
            
        #--------------------------------------
        # For testing.  Check where dz_dt = 0
        #--------------------------------------
##        w  = np.where(self.dz_dt == 0)
##        nw = w[0].size
##        if (nw != 0):
##            print 'dz_dt grid is zero at', nw, 'pixel(s).'
##            if (nw <= 20):
##                for k in xrange(nw):
##                    print '(row, col) =', w[0][k], w[1][k]
##            self.dz_dt[w] = self.BLR_mpyr  # (make positive)
        
    #   update_dz_dt_grid()
    #-------------------------------------------------------------------
    def update_dt_grid(self, SILENT=True, REPORT=False,
                       SAVE_RTG=False):

##        self.update_dt_grid_method0( SILENT=SILENT, REPORT=REPORT,
##                                     SAVE_RTG=SAVE_RTG )

        self.update_dt_grid_method1( SILENT=SILENT, REPORT=REPORT,
                                     SAVE_RTG=SAVE_RTG )
        
        #--------------------------------
        # Defined in erode_base.py
        # Uses min_dz_up, experimental.
        #--------------------------------
##        self.update_dt_grid_method2( SILENT=SILENT, REPORT=REPORT,
##                                     SAVE_RTG=SAVE_RTG )
        
    #   update_dt_grid()
    #-------------------------------------------------------------------
    def update_dt_grid_method0(self, SILENT=True, REPORT=False,
                               SAVE_RTG=False):

        #------------------------------------------------------------
        # Notes: Compute dt value that each pixel needs to
        #        satisfy an experimental stability condition.
        #
        #        Idea is that no D8 "kid pixel" can contribute
        #        enough sediment to its D8 parent to make its
        #        parent have a higher elevation than it does.
        #
        #        dt < fac * da * (del_z / del_Qs)
        #
        #        fac    = factory of safety
        #        DEM    = current elevation grid  [m]
        #        Qs     = K * Q^m * S^n = sed. discharge [m^3/yr]
        #        da     = pixel area grid  [m^2]
        #        del_z  = elevation drops to D8 parent pixels [m]
        #        del_Qs = sed. discharge "drops" to D8 parents
        #------------------------------------------------------------
        # Notes: fac is a factory of safety.  Decreased fac from
        #        0.8 to 0.4 on 8/10/04.  Tried factor of 0.6 on
        #        8/12/04, but for m = n = 1, one or more pits would
        #        form near the outlet so decreased back to 0.4.
        #------------------------------------------------------------     
        if not(SILENT):    
            print('Updating dt_grid...')

        #-----------------------
        # Initialize some vars
        #-----------------------
        dt_too_small = np.float64(1E-2)
        # dt_too_big = 5d  # (should there be a max?)

        #----------------------------------------------------------- 
        # (8/17/10) RT Time Profiles at random pixels shows that
        # state vars are still noisy (oscillating) for some pixels
        # when we use 0.4, but much less so than when using 1.0.
        # Reduced default to 0.3.
        #-----------------------------------------------------------         
        # fac = np.float64(0.4)   # (less than 1)
        fac = np.float64(0.3)   # (less than 1)
        
        #------------------------------------------------------------- 
        # (8/17/10) With default parameters and n_steps=5000,
        # using fac = 1.0 also works.  Despite very similar
        # drainage patterns the base level drops to -56.9028
        # vs. -24.5874.  Max elevations are 2.02733 and 2.9543,
        # respectively, and distribution of elevations is different.
        # Total simulated time is 541213 vs. 218060.
        # These comments are for an erode_d8_global simulation.
        #-------------------------------------------------------------
        # fac = np.float64(1)
        
        #---------------------------
        # Compute downstream drops
        #----------------------------------------------------
        # Get a nonzero value if a pixel's elevation is
        # greater than than that of its parent pixel.
        #----------------------------------------------------
        # #####  THIS SHOULD TRUE EXCEPT FOR "FLATS" #######
        #----------------------------------------------------
        # Computing pIDs as: "self.d8.parent_IDs" works for
        # erode_d8_global but doesn't for erode_d8_local.
        #----------------------------------------------------        
        pIDs  = divmod(self.d8.parent_ID_grid, self.nx)
        del_z = np.maximum((self.DEM - self.DEM[ pIDs ]), 0)
        
        #------------------------
        # Experiment 2: 8/10/04
        #----------------------------------------------------
        # Get a nonzero value if the amount leaving a pixel
        # is more than the amount leaving its parent pixel
        #----------------------------------------------------
        del_Qs = np.maximum((self.Qs - self.Qs[pIDs]), 0)
        
        #--------------------------------------------------
        # Initialize dt_grid to a very long dt value.
        # This becomes the value at non-deposition sites.
        #--------------------------------------------------
        ###########################################################
        ###########################################################
        # NB! Original version had default of 1 year vs. 1e9 !!!
        ###########################################################
        ###########################################################
        self.dt_grid = np.zeros([self.ny, self.nx], dtype='float64')  # [years]
        self.dt_grid += 1e9  # [years]
        ## self.dt_grid += 10000.0 # [years]  # (this is a mid-range value)

        ## self.dt_grid = np.ones([self.ny, self.nx], dtype='float64')  # [years]
        
        #----------------------------
        # Find the deposition sites
        #----------------------------
        wd = np.where(np.logical_and((del_Qs > 0), (del_z > 0)))
        nd = wd[0].size 
        if (nd != 0):
            #------------------------------
            # Compute stable dt, in years
            #------------------------------
            term = fac * (del_z[wd] / del_Qs[wd])
            if (self.da.size == 1):
                self.dt_grid[wd] = term * self.da 
            else:
                self.dt_grid[wd] = term * self.da[wd]
        else:
            if not(SILENT):
                print('-------------------------------------------')
                print(' WARNING: There are no deposition sites.')
                print('-------------------------------------------')
                print(' ')
            
        #--------------------------
        # Save the min and max dt
        #--------------------------
        self.dt_min = self.dt_grid.min()
        self.dt_max = self.dt_grid.max()
        ## if not(SILENT):
        if (self.DEBUG):
            print('#########################################')
            print(' dt_min =', np.around(self.dt_min, 2))
            print(' dt_max =', np.around(self.dt_max, 2))
            print(' in update_dt_grid() (del_z method)')
            print('#########################################')
        
        #-------------------------------
        # Don't let dt get too small ?
        #-----------------------------------------------------
        # dt_min must be less than 1e-4 for the case m=1.5,
        # n=1.0, K=1.0. Making K smaller allows dt_too_small
        # to be bigger. Now K is smaller.
        #-----------------------------------------------------
        # Fixed units issue in Qs_Grid so should be able to
        # reduce dt_too_small.
        #-----------------------------------------------------
        dt_grid_min = np.nanmin(self.dt_grid)
        if (dt_grid_min < dt_too_small):   
            print('******************************************')
            print(' Aborting: Stable dt is too small.')
            print(' Computed dt = ' + str(dt_grid_min))
            print('******************************************')
            print(' ')
            sys.exit()

        #------------------
        # Optional report
        #------------------
        if (REPORT):
            # dt_grid_min = np.nanmin(self.dt_grid)
            dt_grid_max = np.nanmax(self.dt_grid)
            dt_str     = str(dt_grid_min)  + ', ' + str(dt_grid_max)
            del_z_str  = str(del_z.min())  + ', ' + str(del_z.max())
            del_Qs_str = str(del_Qs.min()) + ', ' + str(del_Qs.max())
            print('    min(dt),     max(dt)     = ' + dt_str + ' [yrs]')
            print('    min(del_z),  max(del_z)  = ' + del_z_str)
            print('    min(del_Qs), max(del_Qs) = ' + del_Qs_str)
            print(' ')

        #-----------------------------
        # Save grid as an RTG file ?
        #-----------------------------
        if (SAVE_RTG):
            RTG_file = (self.case_prefix + '_dt.rtg')
            rtg_files.write_grid( self.dt_grid, RTG_file, self.rti )
            
    #   update_dt_grid_method0()
    #-------------------------------------------------------------------
    def update_dt_grid_method1(self, SILENT=True, REPORT=False,
                               SAVE_RTG=False):

        #---------------------------------------------------------------
        # Notes: This version is based on the following stability
        #        condition:
        #
        #        dt < dx * dy * [z(k) - z(p)] / [2 * Qs(k)]
        #
        #        This condition can be derived from the fact that
        #        the rate of sediment transport between cells k and
        #        p, Qs(k), must be less than the rate one gets by
        #        taking half of the sediment volume at cell k that
        #        is higher than cell p and moving it to cell p in
        #        one timestep.
        #
        #        This condition ensures a well-defined timestep
        #        everywhere except at pits and flats.  It may,
        #        however, allow the elevation at cell p to exceed
        #        the elevation at cell k during a (global) timestep
        #        if there is a sediment flux from other neighbors
        #        into cell k and sufficiently low flux out of cell p.
        #
        #        It could be that with global timesteps, z(k) will
        #        be raised by inflow from its neighbors instead of
        #        lowered by the amount used here.  But the result
        #        should still be valid.
        #---------------------------------------------------------------      
        if not(SILENT):    
            print('Updating dt_grid...')
        
        #------------------------------------
        # Compute downstream elevation drop
        #------------------------------------------------
        # del_z will always be positive as long as D8
        # flow code is not zero and LINK_FLATS = False.
        #------------------------------------------------
        # NB! del_z will be wrong wherever pID=0 since
        # it will be computed from DEM[0].  Places with
        # pID=0 are fixed further down.
        #------------------------------------------------
        pIDs  = self.d8.parent_ID_grid
        del_z = (self.DEM - self.DEM.flat[ pIDs ])

        #------------------------------
        # Compute stable dt, in years
        #---------------------------------------------
        # NOTE: CFL_factor is set in set_constants()
        #       function within erode_base.py.
        #--------------------------------------------------
        # See the Wikipedia article for Courant number;
        # it differs between 1D and 2D models.
        # CFL condition is necessary, but not sufficient.
        # Also, this is actually something different.
        #--------------------------------------------------
        self.dt_grid = np.zeros( [self.ny, self.nx], dtype='float64' )
 
        w1 = np.where( self.Qs > 0 )
        if (w1[0].size > 0):
            if (self.da.size == 1):
                da = self.da
            else:
                da = self.da[ w1 ]
            ## da = (self.rti.xres * self.rti.yres)
            dt_equal = (da / 2) * del_z[w1] / self.Qs[w1]
            self.dt_grid[ w1 ] = self.CFL_factor * dt_equal

        #--------------------------------------
        # This should never happen. (11/1/10)
        #--------------------------------------
        w2 = np.where( del_z < 0)
        if (w2[0].size > 0):
            self.dt_grid[ w2 ] = self.dt_limit
            
        #----------------------------------------------------
        # Deactivate all pixels that have a flow code of 0.
        # Their z values can then only be changed when they
        # are synchronized by one of their neighbors.
        #----------------------------------------------------
        w4 = np.where( pIDs == 0 )
        if (w4[0].size > 0):
            self.dt_grid[ w4 ] = self.dt_limit
            
        #--------------------------
        # Save the min and max dt
        #--------------------------
        self.dt_min = np.nanmin( self.dt_grid )
        self.dt_max = np.nanmax( self.dt_grid )
##        self.dt_min = self.dt_grid.min()
##        self.dt_max = self.dt_grid.max()
        ## if not(SILENT):
        if (self.DEBUG):
            print('#########################################')
            print(' dt_min =', np.around(self.dt_min, 2))
            print(' dt_max =', np.around(self.dt_max, 2))
            print(' in update_dt_grid() (del_z/Qs method)')
            print('#########################################')
        
        #-------------------------------
        # Don't let dt get too small ?
        #-----------------------------------------------------
        # dt_min must be less than 1e-4 for the case m=1.5,
        # n=1.0, K=1.0. Making K smaller allows dt_too_small
        # to be bigger. Now K is smaller.
        #-----------------------------------------------------
        # Fixed units issue in Qs_Grid so should be able to
        # reduce dt_too_small.
        #-----------------------------------------------------
        ## if (self.dt_min <= 0):
        if (self.dt_min < self.dt_too_small):
            print('### WARNING: dt < ' + str(self.dt_too_small))
            print('        dt = ' + str(self.dt_min))
            print(' ')
            if (self.dt_min <= 0):
                sys.exit()

        #------------------
        # Optional report
        #------------------
        if (REPORT):
            print('min(del_z), max(del_z) =', del_z.min(), del_z.max())
            print('min(Qs), max(Qs) =', self.Qs.min(), self.Qs.max())
            print('min(dt), max(dt) =', self.dt_grid.min(), self.dt_grid.max())
            wm = np.where( self.dt_grid != self.dt_grid.max() )
            print('max(dt_used)     =', self.dt_grid[wm].max())
            print(' ')

        #-----------------------------
        # Save grid as an RTG file ?
        #-----------------------------
        if (SAVE_RTG):
            RTG_file = (self.case_prefix + '_dt.rtg')
            rtg_files.write_grid( self.dt_grid, RTG_file, self.rti )
            
    #   update_dt_grid_method1()
    #-------------------------------------------------------------------
    def update_DEM(self, SILENT=True, REPORT=False):

        #------------------------------------------------------------
        # Notes: This version computes dt dynamically to ensure
        #        stability.
        #
        #        DEM   = current elevation grid  [m]
        #        Qs    = K * Q^m * S^n = sed. discharge [m^3/yr]
        #        da    = pixel area grid  [m^2]
        #        dt    = channel flow timestep  [s]  (returned)
        #        del_z = elevation drops to parent pixels [m]
        #        dz    = net elevation change due to sed. flux [m]
        #        U     = tectonic uplift rate [mm/year]
        #        w1    = IDs of pixels that...
        #        p1    = IDs of parent pixels that...
        #
        #        NB!  Don't want elevations to change at pixels
        #        where base level is fixed, as with flow to sea.
        #------------------------------------------------------------
        if not(SILENT):    
            print('Updating elevations...')

        #--------------------------------------------------------
        # The update_dt_grid() function was called previously
        # and saved dt_min (using min vs. nanmin).  The dt_grid
        # already includes any factor of safety.
        #--------------------------------------------------------
        self.dt = self.dt_min
        
        #-----------------------------------------------------------
        # Compute dz and the largest dz obtained during this
        # timestep. This is used to check stability & convergence.
        #-----------------------------------------------------------
        dz = self.dz_dt * self.dt
        self.dz_max.fill( np.nanmax(dz) )  ### (2/7/13)
        # dz_max_vec has size of n_steps.  (2/13/12)
        if (self.stop_code == 0):
            self.dz_max_vec[ self.time_index - 1 ] = self.dz_max
        
        #----------------
        # Add dz to DEM
        #----------------
        self.DEM += dz
        self.dz = dz      # (for option to save it)
        
##        self.DEM += np.float32(dz)
##        self.dz   = np.float32(dz)   # (for option to save it, etc.)

        #------------------
        # Optional report
        #------------------
        if (REPORT):
            dz_str = str(np.nanmin(dz))   + ', ' + str(self.dz_max)
            print('    min(dz),     max(dz)     = ' + dz_str + ' [m]')
            
    #   update_DEM()
    #-------------------------------------------------------------------
##    def update_DEM2(self, SILENT=True, REPORT=False):
##
##        #------------------------------------------------------------
##        # Notes: This version computes dt dynamically to ensure
##        #        stability.
##
##        #        DEM   = current elevation grid  [m]
##        #        Qs    = kf * Q^mf * S^nf = sed. discharge [m^3/yr]
##        #        da    = pixel area grid  [m^2]
##        #        dt    = channel flow timestep  [s]  (returned)
##        #        del_z = elevation drops to parent pixels [m]
##        #        dz    = net elevation change due to sed. flux [m]
##        #        U     = tectonic uplift rate [mm/year]
##        #        w1    = IDs of pixels that...
##        #        p1    = IDs of parent pixels that...
##
##        #        NB!  Don't want elevations to change at pixels
##        #        where base level is fixed, as with flow to sea.
##        #------------------------------------------------------------
##        if not(SILENT):    
##            print 'Updating elevations...'
##
##        #--------------
##        # For testing
##        #--------------
##        # self.dt = self.dt_grid.min()            ###################
##        # w  = where(self.dt_grid == self.dt)
##        #---------------------------------------        
##        self.dt = 1.01 * self.dt_grid.min()       ###################
##        # self.dt = 2.0 * self.dt_grid.min()
##        w  = where(self.dt_grid <= self.dt)
##        #---------------------------------------
##        nw = w[0].size
##        print 'number of pixels with min dt =', nw
##        dt_min_row = w[0][0]
##        dt_min_col = w[1][0]
##        print '   (row, col)         = ', dt_min_row, dt_min_col
##        print '   dz_dt(row, col)    =', self.dz_dt[dt_min_row, dt_min_col]
##        print '   code(row,col)      =', self.d8.d8_grid[dt_min_row, dt_min_col]
##        print '   min_dz_up(row,col) =', self.min_dz_up_grid[dt_min_row, dt_min_col]
##        print '   dz_dt(row,col)* dt =', self.dz_dt[dt_min_row, dt_min_col] * self.dt
##        print '   DEM(row,col)       =', self.DEM[dt_min_row, dt_min_col]
##        print '   S(row,col)         =', self.S[dt_min_row, dt_min_col]
##        print '   Qs(row,col)        =', self.Qs[dt_min_row, dt_min_col]
##        print 'self.dt      =', self.dt
##        print 'min(dt_grid) =', self.dt_grid.min()
##        print '------------------------------------------'
##        
##        #-----------------------------------
##        # Should there be a max timestep ?
##        #-----------------------------------
##        #*** dt_max = 5d
##        #*** self.dt = (self.dt < dt_max)
##        #*** print,'dt = ' + str(self.dt)
##        
##        #-------------------------------
##        # Don't let dt get too small ?
##        #-----------------------------------
##        # dt_min must be less than 1e-4
##        # for case mf=1.5, nf=1.0, kf=1.0.
##        # Making kf smaller allows dt_min
##        # to be bigger. Now kf is smaller.
##        #-----------------------------------
##        # Fixed units issue in Qs_Grid so
##        # should be able to reduce dt_min.
##        #-----------------------------------
##        dt_min = float32(1E-2)
##        if (self.dt < dt_min):    
##            print '******************************************'
##            print 'Aborting: Stable dt is too small.'
##            print 'Computed dt = ' + str(self.dt)
##            print '******************************************'
##            print ' '
##            sys.exit()
##        
##        #-----------------------------------------------------------
##        # Compute dz and the largest dz obtained during this
##        # timestep. This is used to check stability & convergence.
##        #-----------------------------------------------------------
##        dz = self.dz_dt * self.dt
##        self.dz_max = nanmax(dz)
##        self.dz_max_vec[ self.time_index - 1 ] = self.dz_max
##        
##        #------------------
##        # Optional report
##        #------------------
##        if (REPORT):
##            dz_str     = str(nanmin(dz))   + ', ' + str(self.dz_max)
##            del_z_str  = str(del_z.min())  + ', ' + str(del_z.max())
##            del_Qs_str = str(del_Qs.min()) + ', ' + str(del_Qs.max())
##            print '    min(dz),     max(dz)     = ' + dz_str + ' [m]'
##            print '    min(del_z),  max(del_z)  = ' + del_z_str
##            print '    min(del_Qs), max(del_Qs) = ' + del_Qs_str
##        
##        #----------------
##        # Add dz to DEM
##        #----------------
##        self.DEM += float32(dz)
##        self.dz   = float32(dz)   # (save for retrieval by caller ?)
## 
##    #   update_DEM()   
    #-------------------------------------------------------------------
    def save_final_grids(self):

        #---------------------------
        # Save "update count" grid
        #---------------------------
        # self.save_update_count_grid()  # (now called in finalize().)

        #-----------------------------------
        # Save last D8 flow direction grid
        #-----------------------------------
        d8_file = (self.out_directory + self.case_prefix + '_lastflow.rtg')
        rtg_files.write_grid( self.d8.d8_grid, d8_file, self.rti,
                              RTG_type='BYTE')

        #--------------------------------------
        # Save last D8 contributing area grid
        #--------------------------------------        
        area_file = (self.out_directory + self.case_prefix + '_lastA.rtg')
        rtg_files.write_grid( self.d8.A, area_file, self.rti)
        
        #---------------------
        # Save the final DEM
        #---------------------
        DEM_file = (self.out_directory + self.case_prefix + '_lastZ.rtg')
        rtg_files.write_grid( self.DEM, DEM_file, self.rti,
                              RTG_type='FLOAT')

        #----------------------------
        # Save the final slope grid
        #----------------------------       
        S_file = (self.out_directory + self.case_prefix + '_lastS.rtg')
        rtg_files.write_grid( self.S, S_file, self.rti,
                              RTG_type='FLOAT')

        #--------------------------------
        # Save the final discharge grid
        #--------------------------------       
        Q_file = (self.out_directory + self.case_prefix + '_lastQ.rtg')
        rtg_files.write_grid( self.Q, Q_file, self.rti,
                              RTG_type='FLOAT')

        #-------------------------------------
        # Save the final sed. discharge grid
        #-------------------------------------       
        Qs_file = (self.out_directory + self.case_prefix + '_lastQs.rtg')
        rtg_files.write_grid( self.Qs, Qs_file, self.rti,
                              RTG_type='FLOAT')

        #------------------------
        # Save the final D grid
        #------------------------       
##        D_file = (self.out_directory + self.case_prefix + '_lastD.rtg')
##        rtg_files.write_grid( self.D, D_file, self.rti,
##                              RTG_type='FLOAT')

        #----------------------------
        # Save the final dz_dt grid
        #----------------------------       
        dz_dt_file = (self.out_directory + self.case_prefix + '_lastdzdt.rtg')
        rtg_files.write_grid( self.dz_dt, dz_dt_file, self.rti,
                              RTG_type='FLOAT')

        #-------------------------
        # Save the final dt grid
        #-------------------------   
##        dt_grid_file = (self.out_directory + self.case_prefix + '_lastdt.rtg')
##        rtg_files.write_grid( self.dt_grid, dt_grid_file, self.rti,
##                              RTG_type='FLOAT')
        
    #   save_final_grids()
    #-------------------------------------------------------------------
    
