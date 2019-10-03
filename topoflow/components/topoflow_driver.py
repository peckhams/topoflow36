
###### 
###### channels_base.py called diversions.update() COMMENTED OUT !!!!!

## Copyright (c) 2001-2014, Scott D. Peckham
##
## Jan 2013. Revised handling of input/output names.
##
## Oct 2012. CSDMS Standard Names and BMI.
##
## May 2012. Q_outlet -> Q_outlet[0]   (5/19/12)
##
## May 2010. Changes to initialize() and read_cfg_file()
##
## Feb 2010. Channels comp. now calls dp.update() itself.)
##
## April, May, July, August 2009
##
## Jan 2009. Converted from IDL.
##

###############################################################
#  Note: We could have the channels component make most of
#        the process update() calls.  This may clean up the
#        appearance in the CMT.

#  NB! T Some of the  "print_*" functions are not ready yet.

#      In CCA mode, all update() method calls require the
#      "time_sec" argument.
###############################################################

#--------------------------------------------------------------------
#  Notes:  This is an object-oriented, Python version of the
#          TopoFlow spatial, hydrologic model, originally developed
#          in IDL.  Every process component now has a Basic Model
#          Interface (BMI) and uses CSDMS Standard Names.
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#  Notes:  Later go through code line-by-line to find places where
#          numpy calculations (e.g. where) could be done more
#          efficiently.
#--------------------------------------------------------------------
#
#  class topoflow_driver      (inherits from BMI_base.py)
#
#      get_component_name()
#      get_attribute()          # (10/26/11)
#      get_input_var_names()    # (5/16/12, Bolton)
#      get_output_var_names()   # (5/16/12, Bolton)
#      get_var_name()           # (5/16/12, Bolton)
#      get_var_units()          # (5/16/12, Bolton)
#      ---------------------
#      set_constants()
#      run_model()
#      initialize()
#      update()
#      finalize()
#      -----------------------------------
#      check_finished()
#      check_steady_state()   ###
#      check_interrupt()
#      -----------------------------------
#      initialize_stop_vars()
#      initialize_mass_totals()       # (OBSOLETE ??)
#      initialize_GUI()               # (commented out)
#      initialize_hydrograph_plot()   # (commented out)
#      -----------------------------------
#      update_mass_totals()           # (OBSOLETE ??)
#      update_hydrograph_plot()
#      -----------------------------------
#      print_final_report()
#      print_mins_and_maxes()
#      print_uniform_precip_data()
#      print_dimless_number_data()
#      print_mass_balance_report()

#-----------------------------------------------------------------------

import numpy as np
import os
import time

from topoflow.utils import BMI_base
from topoflow.utils import cfg_files as cfg

from topoflow.utils.tf_utils import TF_Print, TF_String, TF_Version
from topoflow.utils import tf_utils   # (must come after the "from" line ???)

## from tf_utils import TF_Print, TF_String
## import tf_utils  # (must come after the "from" line ???)

#-----------------------------------------------------------------------
class topoflow_driver( BMI_base.BMI_component ):

    #-------------------------------------------------------------------
    # Don't define an __init__() function.
    # We need to inherit the BMI_base.__init__()
    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Driver',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'TopoFlow',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'TopoFlow.cfg.in',
        'cfg_extension':      '_topoflow.cfg',
        'cmt_var_prefix':     '/TopoFlow/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/TopoFlow.xml',
        'dialog_title':       'TopoFlow: Driver Parameters',
        'time_units':         'seconds' }

    #------------------------------------------------------------------------------
    # (2/3/13) Added "channel_model__time_step" since it is needed by
    # both the TopoFlow driver (topoflow.py) and Diversions (diversion_base.py).
    # But source and sink files provide "dt" for Diversions, so check.
    #------------------------------------------------------------------------------
    _input_var_names = [
        'atmosphere_water__domain_time_integral_of_precipitation_leq-volume_flux', # vol_P@meteorology
        'atmosphere_water__domain_time_max_of_precipitation_leq-volume_flux',    # P_max@meteorology
        'basin_outlet_water_flow__half_of_fanning_friction_factor',              # f_outlet@channels
        'basin_outlet_water_x-section__mean_depth',                              # d_outlet@channels
        'basin_outlet_water_x-section__peak_time_of_depth',                      # Td_peak@channels
        'basin_outlet_water_x-section__peak_time_of_volume_flow_rate',           # T_peak@channels
        'basin_outlet_water_x-section__peak_time_of_volume_flux',                # Tu_peak@channels
        'basin_outlet_water_x-section__time_integral_of_volume_flow_rate',       # vol_Q@channels
        'basin_outlet_water_x-section__time_max_of_mean_depth',                  # d_peak@channels
        'basin_outlet_water_x-section__time_max_of_volume_flow_rate',            # Q_peak@channels
        'basin_outlet_water_x-section__time_max_of_volume_flux',                 # u_peak@channels
        'basin_outlet_water_x-section__volume_flow_rate',                        # Q_outlet@channels
        'basin_outlet_water_x-section__volume_flux',                             # u_outlet@channels
        'channel_bottom_water_flow__domain_max_of_log_law_roughness_length',     # z0val_max@channels
        'channel_bottom_water_flow__domain_min_of_log_law_roughness_length',     # z0val_min@channels
        ## 'channel_model__time_step',  ####################### (no longer needed?)
        'channel_water_flow__domain_max_of_manning_n_parameter',                 # nval_max@channels
        'channel_water_flow__domain_min_of_manning_n_parameter',                 # nval_min@channels
        #-----------------------------------------------------
        # These might only be available at the end of run ??
        # These are now over the entire domain (or DEM).
        #-----------------------------------------------------
        'channel_water_x-section__domain_max_of_mean_depth',                     # d_max
        'channel_water_x-section__domain_max_of_volume_flow_rate',               # Q_max
        'channel_water_x-section__domain_max_of_volume_flux',                    # u_max
        'channel_water_x-section__domain_min_of_mean_depth',                     # d_min
        'channel_water_x-section__domain_min_of_volume_flow_rate',               # Q_min
        'channel_water_x-section__domain_min_of_volume_flux',                    # u_min
#         'snowpack__domain_max_of_depth',                                       # hs_max
#         'snowpack__domain_min_of_depth',                                       # hs_min
        #-----------------------------------------------------------       
        'glacier_ice__domain_time_integral_of_melt_volume_flux',                 # vol_MR@ice
        'land_surface_water__baseflow_volume_flux',                              # GW@satzone
        'land_surface_water__domain_time_integral_of_baseflow_volume_flux',      # vol_GW@satzone
        'land_surface_water__domain_time_integral_of_evaporation_volume_flux',   # vol_ET@evap
        'land_surface_water__domain_time_integral_of_runoff_volume_flux',        # vol_R@channels
        'land_surface_water__runoff_volume_flux',                                # R@channels
        'snowpack__domain_time_integral_of_melt_volume_flux',                    # vol_SM
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux',  # vol_IN@infil
        'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux', # vol_Rg
        #-----------------------------------------------------------   
        'network_channel_water__volume',                                         # vol_chan@channels
        'land_surface_water__area_integral_of_depth' ]                           # vol_land@channels
        
        #----------------------------------------------------------------
        # The TopoFlow driver no longer needs to get the time_steps of
        # the other model components; this is now the framework's job.
        # That means it can't include them in its final report.
        #----------------------------------------------------------------
##        'channel:model__time_step',     # dt@channels
##        'diversion:model__time_step',
##        'evap:model_time_step',
##        'ice:model_time_step',
##        'infil:model__time_step',
##        'meteorology:model__time_step',
##        'satzone:model__time_step',
##        'snow:model__time_step' ]

        ###################################################################
        #### Bolton comments, 5/12/2012  ---
        ####          Not sure what to do with these missing /unknow vars
         ####          cp.get_status()
        ####          save_pixels_dt@channels    'model__save_pixels_flag' ? 
        ####          MANNING@channels           'model__manning_flag' ?
        ####          LAW_OF_WALL@channels       'model__wall_law_flag ?
        ####          RICHARDS@infiltration      'model__richards_flag' ?
        ###################################################################

    _output_var_names = [
        'model__time_step' ]   # dt

    _var_name_map = {
        'atmosphere_water__domain_time_integral_of_precipitation_leq-volume_flux': 'vol_P',
        'atmosphere_water__domain_time_max_of_precipitation_leq-volume_flux':      'P_max',
        'basin_outlet_water_flow__half_of_fanning_friction_factor':          'f_outlet', 
        'basin_outlet_water_x-section__mean_depth':                          'd_outlet',
        'basin_outlet_water_x-section__peak_time_of_depth':                  'Td_peak',
        'basin_outlet_water_x-section__peak_time_of_volume_flow_rate':       'T_peak',
        'basin_outlet_water_x-section__peak_time_of_volume_flux':            'Tu_peak',
        'basin_outlet_water_x-section__time_integral_of_volume_flow_rate':   'vol_Q',
        'basin_outlet_water_x-section__time_max_of_mean_depth':              'd_peak',
        'basin_outlet_water_x-section__time_max_of_volume_flow_rate':        'Q_peak',
        'basin_outlet_water_x-section__time_max_of_volume_flux':             'u_peak',
        'basin_outlet_water_x-section__volume_flow_rate':                    'Q_outlet',
        'basin_outlet_water_x-section__volume_flux':                         'u_outlet', 
        'channel_bottom_water_flow__domain_max_of_log_law_roughness_length': 'z0val_max',
        'channel_bottom_water_flow__domain_min_of_log_law_roughness_length': 'z0val_min',   
        ## 'channel_model__time_step':                                          'channel_dt', ## (2/3/13)
        'channel_water_flow__domain_max_of_manning_n_parameter':             'nval_max',
        'channel_water_flow__domain_min_of_manning_n_parameter':             'nval_min',        
        #-------------------------------------------------------
        # These 6 might only be available at the end of run ??
        #-------------------------------------------------------
        'channel_water_x-section__domain_max_of_mean_depth':                     'd_max',
        'channel_water_x-section__domain_max_of_volume_flow_rate':               'Q_max',
        'channel_water_x-section__domain_max_of_volume_flux':                    'u_max',
        'channel_water_x-section__domain_min_of_mean_depth':                     'd_min',
        'channel_water_x-section__domain_min_of_volume_flow_rate':               'Q_min',
        'channel_water_x-section__domain_min_of_volume_flux':                    'u_min',
#         'snowpack__domain_max_of_depth':                                       'hs_max',
#         'snowpack__domain_min_of_depth':                                       'hs_min',
        #------------------------------------------------------------                
        'glacier_ice__domain_time_integral_of_melt_volume_flux':                 'vol_MR',
        'land_surface_water__baseflow_volume_flux':                              'GW',
        'land_surface_water__domain_time_integral_of_baseflow_volume_flux':      'vol_GW',
        'land_surface_water__domain_time_integral_of_evaporation_volume_flux':   'vol_ET',
        'land_surface_water__domain_time_integral_of_runoff_volume_flux':        'vol_R',
        'land_surface_water__runoff_volume_flux':                                'R',        
        'snowpack__domain_time_integral_of_melt_volume_flux':                    'vol_SM',       
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux':  'vol_IN',
        'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux': 'vol_Rg',
        #---------------------
        'model__time_step':                                      'dt',
        #----------------------------------------------------------------------------   
        'network_channel_water__volume':                         'vol_chan',
        'land_surface_water__area_integral_of_depth':            'vol_land' }       

            
    _var_units_map = {
        'atmosphere_water__domain_time_integral_of_precipitation_leq-volume_flux': 'm3',
        'atmosphere_water__domain_time_max_of_precipitation_leq-volume_flux':      'm s-1',
        'basin_outlet_water_flow__half_of_fanning_friction_factor':                '1',
        'basin_outlet_water_x-section__mean_depth':                                'm',
        'basin_outlet_water_x-section__peak_time_of_depth':                        'min',
        'basin_outlet_water_x-section__peak_time_of_volume_flow_rate':             'min',
        'basin_outlet_water_x-section__peak_time_of_volume_flux':                  'min',
        'basin_outlet_water_x-section__time_integral_of_volume_flow_rate':         'm3',
        'basin_outlet_water_x-section__time_max_of_mean_depth':                    'm',
        'basin_outlet_water_x-section__time_max_of_volume_flow_rate':              'm3 s-1',
        'basin_outlet_water_x-section__time_max_of_volume_flux':                   'm s-1',
        'basin_outlet_water_x-section__volume_flow_rate':                          'm3 s-1',
        'basin_outlet_water_x-section__volume_flux':                               'm s-1',
        'channel_bottom_water_flow__domain_max_of_log_law_roughness_length':       'm',
        'channel_bottom_water_flow__domain_min_of_log_law_roughness_length':       'm',        
        ## 'channel_model__time_step':                                                's', ### (2/3/13)
        'channel_water_flow__domain_max_of_manning_n_parameter':                   'm-1/3 s',        
        'channel_water_flow__domain_min_of_manning_n_parameter':                   'm-1/3 s',
        #-----------------------------------------------------
        # These might only be available at the end of run ??
        #-----------------------------------------------------
        'channel_water_x-section__domain_max_of_mean_depth':                       'm',
        'channel_water_x-section__domain_max_of_volume_flow_rate':                 'm3 s-1',
        'channel_water_x-section__domain_max_of_volume_flux':                      'm s-1',
        'channel_water_x-section__domain_min_of_mean_depth':                       'm',
        'channel_water_x-section__domain_min_of_volume_flow_rate':                 'm3 s-1',
        'channel_water_x-section__domain_min_of_volume_flux':                      'm s-1',
#         'snowpack__domain_max_of_depth':                                         'm',
#         'snowpack__domain_min_of_depth':                                         'm',
        #------------------------------------------------------------                
        'glacier_ice__domain_time_integral_of_melt_volume_flux':                   'm3',
        'land_surface_water__baseflow_volume_flux':                                'm s-1',
        'land_surface_water__domain_time_integral_of_baseflow_volume_flux':        'm3',
        'land_surface_water__domain_time_integral_of_evaporation_volume_flux':     'm3',
        'land_surface_water__domain_time_integral_of_runoff_volume_flux':          'm3',
        'land_surface_water__runoff_volume_flux':                                  'm s-1',
        'network_channel_water__volume':                                           'm3',
        'snowpack__domain_time_integral_of_melt_volume_flux':                      'm3',        
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux':    'm3',
        'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux':   'm3',
        #----------------------------
        'model__time_step': 's',
        #----------------------------------------------------------------------------   
        'network_channel_water__volume':                         'm3',
        'land_surface_water__area_integral_of_depth':            'm3'   }  
    
    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Driver'  ##### TopoFlow_Run_Monitor,  Report_Maker?

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
##    def get_output_var_names(self):
##
##        names = ['hydro_model_time_step_size' ]   # dt   
##        
##        return names
####      return np.array( names )    # (string array vs. list)
##    
##    #   get_output_var_names()
    #-------------------------------------------------------------------
##    def get_var_name(self, long_var_name):
##
##        #-------------------------------------------------
##        # Define this map just once in "__init__()"  ??
##        #-------------------------------------------------
##        name_map = {
##            'basin_cumulative_discharged_water_volume':'vol_Q',
##            'basin_cumulative_evaporated_water_volume':'vol_ET',
##            'basin_cumulative_ice_meltwater_volume':'vol_MR',
##            'basin_cumulative_infiltrated_water_volume':'vol_IN',
##            'basin_cumulative_lwe_precipitated_water_volume':'vol_P',
##            'basin_cumulative_runoff_water_volume':'vol_R',
##            'basin_cumulative_saturated_zone_infiltrated_water_volume':'vol_Rg',
##            'basin_cumulative_snow_meltwater_volume':'vol_SM',
##            'basin_cumulative_subsurface_to_surface_seeped_water_volume':'vol_GW',
##            'basin_outlet_water_discharge':'Q_outlet',
##            'basin_outlet_water_mean_depth':'d_outlet',   ####           
##            'basin_outlet_water_mean_speed':'u_outlet',   ####
##            'channel_bed_max_manning_roughness_parameter':'nval_max',
##            'channel_bed_max_roughness_length':'z0val_max',
##            'channel_bed_min_manning_roughness_parameter':'nval_min',
##            'channel_bed_min_roughness_length':'z0val_min',
##            'channel_outgoing_peak_water_discharge':'Q_peak',
##            'channel_outgoing_peak_water_discharge_time':'T_peak',
##            'channel_time_step_size':'channel_dt',                     ####### Need "_" vs. ".".
##            'channel_water_peak_mean_depth':'d_peak',
##            'channel_water_peak_mean_depth_time':'Td_peak',
##            'channel_water_peak_mean_speed':'u_peak',
##            'channel_water_peak_mean_speed_time':'Tu_peak',
##            'diversion_time_step_size':'diversion_dt',
##            'evap_time_step_size':'evap_dt',      ####         
##        'hydro_model_time_step_size':'dt',
##            'ice_time_step_size':'ice_dt',        ####
##            'infil_time_step_size':'infil_dt',    ####         
##            'lwe_max_precipitation_rate':'P_max',
##            'meteorology_time_step_size':'meteorology_dt',   ####
##            'satzone_time_step_size':'satzone_dt',           ####
##            'snow_time_step_size':'snow_dt' }                ####
##            
##        return name_map[ long_var_name ]
##
##    #   get_var_name()
##    #-------------------------------------------------------------------
##    def get_var_units(self, long_var_name):
##
##        #-------------------------------------------------
##        # Define this map just once in "__init__()"  ??
##        #-------------------------------------------------
##        units_map = {
##            'vol_ET'         : 'm3',
##            'vol_MR'         : 'm3',
##            'vol_IN'         : 'm3',
##            'vol_P'          : 'm3',
##            'vol_R'          : 'm3',
##            'vol_Rg'         : 'm3',
##            'vol_SM'         : 'm3',
##            'vol_GW'         : 'm3',
##            'Q_outlet'       : 'm3 s-1',
##            'd_outlet'       : 'm',
##            'u_outlet'       : 'm s-1',
##            'nval_max'       : 's m-1/3',   ##### CHECK
##            'z0val_max'      : 'm',
##            'nval_max'       : 's m-1/3',
##            'nval_min'       : 's m-1/3',
##            'z0val_min'      : 'm',
##            'Q_peak'         : 'm3 s-1',
##            'T_peak'         : 'min',
##            'd_peak'         : 'm',
##            'Td_peak'        : 'min',
##            'u_peak'         : 'm s-1',
##            'Tu_peak'        : 'min',
##            # 'z0val'     : 'm',
##            'P_max'          : 'm s-1',
##            'dt'             : 's',
##            'channel_dt'     : 's',
##            'diversion_dt'   : 's',
##            'evap_dt'        : 's',
##            'ice_dt'         : 's',
##            'infil_dt'       : 's',
##            'satzone_dt'     : 's',
##            'snow_dt'        : 's',
##            'meteorology_dt' : 's' }
##            
##        var_name = self.get_var_name( long_var_name )
##        
##        return units_map[ var_name ]
##    
##    #   get_var_units()
    #-------------------------------------------------------------------
    def set_constants(self):

        #------------------------
        # Define some constants
        #------------------------
        self.mps_to_mmph = np.float64(3600000)   # [m/s] to [mm/hr]

        #------------------------------------------
        # Moved these from __init__() on 5/17/12.
        #------------------------------------------
        self.OK           = True
        self.comment      = None
        self.WRITE_LOG    = True
        self.VERBOSE      = False
        self.PLOT         = False
       
    #   set_constants()       
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        #------------------------------------------------------
        # Note:  If using as a CCA component, then we need to
        #        call "get_cca_ports()" before calling this
        #        "initialize()" method in the component's CCA
        #        Impl file.
        #------------------------------------------------------
        if not(SILENT):
            print('TopoFlow component: Initializing...')
        
        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode
        self.cfg_file   = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()
        self.initialize_config_vars() 
        self.read_grid_info()
        self.initialize_basin_vars()  # (5/14/10)

        #----------------------------------
        # Has component been turned off ?
        #----------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print('TopoFlow Main component: Disabled in CFG file.')
            self.DONE = True
            self.status = 'initialized.'  # (OpenMI 2.0 convention)
            return
        
        dc = (self.out_directory + self.case_prefix)
        self.comment_file = dc + '_comments.txt'
        self.log_file     = dc + '.log'

        #-----------------------
        # Initialize variables
        #-----------------------
        self.initialize_time_vars()  # (uses cp.dt from above)
        self.initialize_stop_vars()   #### (do this with CFG file ?)
                          
        #### self.nz = self.ip.get_scalar_long('nz')   #######

        #------------------------------------
        # Check if output options are valid
        #------------------------------------
##        outlet_IDs, n_outlets, outlet_ID, OK = \
##            Check_Output_Options(grid_vars, rv, cv, sv, ev, iv, gv, nx, ny)
      
        #---------------------
        # Open the logfile ?      *** (Append or overwrite ??)
        #---------------------
        if (self.WRITE_LOG):
            if (self.log_file == None):
                self.log_file = (self.case_prefix + '_LOG.txt')
            print('Opening log file:')
            print('    log_file = ' + self.log_file)
            self.log_unit = open(self.log_file, 'w')

        #----------------------
        # Print start message
        #----------------------
        TF_Print(' ')
        TF_Print('Starting TopoFlow model run...')
        hline = ''.ljust(60, '-')
        TF_Print(hline)
        self.status = 'initialized'
        
    #   initialize()        
    #-------------------------------------------------------------            
    ## def update(self, dt=-1.0, time_seconds=None):
    def update(self, dt=-1.0):

        #--------------------------------------------------
        # Note:  This method no longer calls the update()
        # methods of other components; this is now taken
        # care of by the framework, regardless of which
        # component is the "driver". (2/4/13)
        #--------------------------------------------------
        
        #-------------------------------
        # Check for interrupt by user ?
        #-------------------------------
        # OK = self.check_interrupt()
        # if not(OK):
        #    self.status = 'stopped'
        #    return
        
        # self.DEBUG = True  ########################
        
        self.status = 'updating'
        OK = True
        ## if (self.VERBOSE):
        if (self.mode == 'driver'):
            self.print_time_and_value(self.Q_outlet, 'Q_out', '[m^3/s]',
                                      interval=0.5)  # [seconds]
            
        ## self.update_hydrograph_plot()

        ## print '##### In topoflow.update(), time_sec =', self.time_sec
        
##        #------------------------------------------------
##        # Update channel process first, using the first
##        # values that were read for each process.
##        #------------------------------------------------
##        self.cp.update(self.time_sec)  # (uses others, so after them)
##        if (self.DEBUG): print 'TF UPDATED Channels Component...'
##
##        OK = (self.cp.get_status() != 'failed')
##        if not(OK):
##            if (self.DEBUG): print 'CHANNELS component failed.'
##            self.Q_outlet = np.float64(0)   # (why is this here again?)
##            self.DONE = True
##            return
##
##        if (self.DEBUG): print 'MADE IT PAST fail check...'
##        
##        #------------------------------------
##        # Step each process forward in time
##        #------------------------------------
##        UPDATE_MP = (self.time_index % self.mp_update_step) == 0
##        UPDATE_SP = (self.time_index % self.sp_update_step) == 0
##        UPDATE_EP = (self.time_index % self.ep_update_step) == 0
##        UPDATE_IP = (self.time_index % self.ip_update_step) == 0
##        UPDATE_GP = (self.time_index % self.gp_update_step) == 0
##        UPDATE_IIP= (self.time_index % self.iip_update_step)== 0
##        ## UPDATE_DP = (self.time_index % self.dp_update_step) == 0
##
##        if (self.DEBUG): print 'MADE IT PAST update checks...'
##        
##        # NB! Precip process now uses fixed dt. (August 09)
##        # NB! Precip has been absorbed into Met (Sept 09)
##
##        if (UPDATE_MP):
##            self.mp.update(self.time_sec)    # update met_vars
##            if (self.DEBUG): print 'TF UPDATED Met Component...'
##        if (UPDATE_SP):
##            self.sp.update(self.time_sec)    # update snow_vars
##            if (self.DEBUG): print 'TF UPDATED Snow Component...'
##        if (UPDATE_EP):
##            self.ep.update(self.time_sec)    # update ET
##            if (self.DEBUG): print 'TF UPDATED Evap Component...'
##        if (UPDATE_IP):
##            self.ip.update(self.time_sec)    # update infil
##            if (self.DEBUG): print 'TF UPDATED Infil Component...'
##        if (UPDATE_IIP):
##            self.iip.update(self.time_sec)   # update ice_vars
##            if (self.DEBUG): print 'TF UPDATED Ice Component...'
##            ## move update q0 into here from merged2.py
##            ## check_infiltration()
##        if (UPDATE_GP):
##            self.gp.update(self.time_sec)   # update GW ?
##            if (self.DEBUG): print 'TF UPDATED GW Component...'
##
##        if (self.DEBUG): print 'MADE IT PAST updates ...'
        
        #---------------------------------------------
        # Note that u and d from previous time step
        # must be used on RHS of the equations here.
        #---------------------------------------------
        # self.cp.update(self.time_sec)  # (uses others, so after them)
        # if (DEBUG): print 'UPDATED Channels Component...'
        
##        OK = (self.cp.get_status() != 'failed')
##        if not(OK):
##            self.Q_outlet = np.float64(0)   # (why is this here?)
##            self.DONE = True
##            return

        #---------------------------------------------
        # (2/1/10) dp.update() is now called by the
        # channels component in update_flow_volume()
        #---------------------------------------------
        # dp.update() uses a setter to modify cp.vol
        #---------------------------------------------
        ## if (UPDATE_DP):
##        self.dp.update(self.time_sec)
##        if (self.DEBUG): print 'TF UPDATED Diversions Component...'

        
        #-------------------------
        # Increment the timestep
        #-------------------------------------------------------
        # Note that the update_time() method in BMI_base.py
        # calls "check_finished()".  There is a simple version
        # of "check_finished()" in BMI_base.py and another
        # version in this file that supports more general
        # stopping conditions, including "steady state".
        #------------------------------------------------------   
        self.update_time( dt )
        self.status = 'updated'
        
    #   update()          
    #-------------------------------------------------------------            
    def finalize(self):

        self.status = 'finalizing'

        #--------------------------------------------
        # This is called by the update() method.
        # We'd like to force another output, but it
        # is now based on elapsed time.  10/29/11.
        #--------------------------------------------        
        # if (self.mode == 'driver'):
        #     self.print_time_and_value(self.Q_outlet, 'Q_out', '[m^3/s]')
 
        print('=======================')
        print('Simulation complete.')
        print('=======================')
        print(' ')
        
        #----------------
        # Print reports
        #----------------
        self.print_final_report()
##        self.print_mins_and_maxes( FINAL=True )   # (Moved into final report.)
##        self.print_uniform_precip_data()  # (not ready yet)
##        self.print_mass_balance_report()  # (not ready yet)

        #--------------------
        # Close the logfile
        #--------------------
        if (self.WRITE_LOG):
            ## print '###  Closing log file.'
            self.log_unit.close()

        #----------------------
        # Print final message
        #----------------------
        TF_Print('Finished.' + '  (' + self.case_prefix + ')')
        TF_Print(' ')
        self.status = 'finalized'
   
    #   finalize()
    #-------------------------------------------------------------            
    def check_finished(self):

        #------------------------------------------------------ 
        # Notes: This version of "check_finished()" overrides
        #        the simpler one that all components inherit
        #        from BMI_base.py.
        #------------------------------------------------------ 
        
        #---------------------------------
        # Methods to stop the simulation
        #---------------------------------
        if (self.stop_method == 'Q_peak_fraction'):
            #----------------------------------------------------
            # Run until the outlet hydrograph drops to less
            # than "Qp_fraction" of the peak value before that.
            #----------------------------------------------------  
            FALLING_LIMB = (self.Q_last > self.Q_outlet)
            
##            if (FALLING_LIMB):  print "ON FALLING LIMB..."
##            print 'self.Q_last  =', self.Q_last
##            print 'self.Q_outlet =', self.Q_outlet
##            print ' '
            
            #--------------------------------------------------------
            # With DYNAMIC_WAVE, it is possible for some reason for
            # Q_outlet to drop back to zero early in the simulation
            # so that model run ends right away.  (2/13/07)
            # Uncomment the debugging section below.
            #--------------------------------------------------------
            if (FALLING_LIMB):
                Q_stop   = (self.Q_peak * self.Qp_fraction)
                self.DONE = (self.Q_outlet <= Q_stop) and \
                            (self.Q_outlet > 0)

            if (self.DONE):
                stop_str = str(self.Qp_fraction) + '.\n'
                print('Stopping: Reached Q_peak fraction = ' + stop_str)
                
##            print 'FALLING_LIMB   =', FALLING_LIMB
##            print 'Q_last         =', self.Q_last
##            print 'Q_peak         =', self.Q_peak
##            print 'Qpeak_fraction =', self.Qp_fraction
##            print ' '
            
                #--------------
                # For testing
                #--------------
                #if (DONE):
                #    print 'Q_last         =', self.Q_last
                #    print 'Q_peak         =', self.Q_peak
                #    print 'Qpeak_fraction =', self.Qp_fraction
                #    print 'Q[self.outlet_ID]   =', Q_outlet
                #    print 'Q_stop         =', Q_stop
                #    print ' '

        elif (self.stop_method == 'Until_model_time'):
            #--------------------------------------------------
            # Run until specified "stopping time", in minutes
            #--------------------------------------------------
            self.DONE = (self.time_min >= self.T_stop)  # [minutes]
            if (self.DONE):
                stop_str = str(self.T_stop) + '.\n'
                print('Stopping: Reached stop time = ' + stop_str)
        elif (self.stop_method == 'Until_n_steps'):
            #--------------------------------------
            # Run for a specified number of steps
            #--------------------------------------
            self.DONE = (self.time_index >= self.n_steps)
            if (self.DONE):
                stop_str = str(self.n_steps) + '.\n'
                print('Stopping: Reached number of steps = ' + stop_str)
        else:
            raise RuntimeError('No match found for expression')

        #--------------------------------------------------------------
        # This model run is finished if the user-selected stopping
        # condition has been reached (above) OR if the model appears
        # to have reached a steady-state condition with discharge,
        # OR if the channel component has failed for some reason.
        # (2/4/13)
        #--------------------------------------------------------------
        FINISHED = self.DONE
        STEADY   = self.check_steady_state()
        # FAILED = (self.cp.get_status() == 'failed')  ###### FIX SOON (5/18/12)
        FAILED   = False
        if (FAILED):
            if (self.DEBUG): print('CHANNELS.update() failed.')
            self.status   = 'failed'
            ### self.Q_outlet = np.float64(0)   # (why is this here?)
        ### self.Q_last = self.Q_outlet
        self.Q_last = self.Q_outlet.copy()  ## (2/7/13)
        self.DONE   = (FINISHED or STEADY or FAILED)

        return self.DONE
      
    #   check_finished()
    #-------------------------------------------------------------------
    def check_steady_state(self):

        #-------------------------------------------------------
        # Notes:  See "initialize_stop_vars()" for definitions
        #         of steady_tol, nonzero_tol, n_same_max and
        #         Q_last.
        #-------------------------------------------------------
        STEADY = False
    
        #------------------------------------------
        # Count number of steps with same Q-value
        #------------------------------------------
        delta_Q = np.absolute(self.Q_outlet - self.Q_last)
        if ( delta_Q <= self.steady_tol ):
            ## print '(time_index, dQ) =', self.time_index, delta_Q
            self.n_same += 1
        else:
            self.n_same  = 0
        
        #------------------------------------
        # Check for steady-state, with Q > 0
        #------------------------------------
        if (self.stop_method == 'Q_peak_fraction') and \
           (self.Q_outlet > self.nonzero_tol) and \
           (self.n_same   > self.n_same_max):

            STEADY = True
            if not(self.DONE):   # (5/19/12. No message if already done.)
                msg = ['-----------------------------------------------------------', \
                       'WARNING: It appears that discharge, Q, has reached', \
                       '         a steady-state condition.', \
                       '         Discharge at outlet near: ' + str(self.Q_outlet), \
                       '             for ' + str(self.n_same) + ' timesteps.', \
                       '         Aborting model run.', \
                       '-----------------------------------------------------------', \
                       ' ']
                       ### 'Do you want to continue anyway ?', ' '])
                for line in msg:
                    print(line)

            
            ## answer = GUI_Message(msg, QUESTION=True)
            ## DONE = (answer.lower() == 'no')
            ## if not(DONE):    
            ##     n_same = int32(0)
            ##
            #TF_Print,'****************************************************'
            #TF_Print,'Aborting model run: '
            #TF_Print,'Discharge, Q, has reached steady-state.'
            #TF_Print,'****************************************************'
            #msg = [ $
            #'WARNING:  Route_Flow aborted.', ' ',$
            #'Discharge, Q, has reached steady-state. ', ' ']
            #GUI_Error_Message, msg
            #STEADY = True
        
        #-----------------------------------------------
        # (3/20/07) Commented out, since user can now
        # press any key to stop the model run.  Note
        # that Q-value can remain zero for a long time
        # if everything is infiltrating or snow depth
        # is building, etc.
        #-----------------------------------------------
        # Check for unchanging Q-value of zero ?
        #-----------------------------------------
        #if (STOP_METHOD eq 0b) AND (Q_peak eq 0.0) AND (n_same gt nn) then begin
        #    msg = [' ', $
        #    'ERROR: ', ' ', $
        #    'Discharge at outlet is zero for all times. ', ' ',$
        #    'Is there a runoff-producing process ? ', ' ']
        #    TF_Print,'*************************************************'
        #    TF_Print, msg[1]
        #    TF_Print, msg[3]
        #    TF_Print,'*************************************************'
        #    TF_Print,' '
        #    GUI_Error_Message, msg
        #    DONE = 1b
        #endif
          
        return STEADY
    
    #   check_steady_state
    #-------------------------------------------------------------
    def check_interrupt(self):

        #------------------------------------------------------------
        #( 3/21/07) This method checks for a keyboard interrupt
        # after a fixed amount of real time has elapsed.  This
        # works much better than using something like (n mod n_check)
        # since it avoids checking fast runs too often (which slows
        # down performance) or slow runs too seldom (which means the
        # user can't interrupt the run).  It only checks the elapsed
        # time once each time through the loop, however, so the time
        # required for one pass imposes a lower bound.
        #------------------------------------------------------------
        elapsed_time = (time.time() - self.last_check_time)
        if (elapsed_time > 2):    
            #print,'****************************************************'
            #print,'Checking interrupt: n = ' + TF_String(n)
            #print,'****************************************************'

            ########################################
            ##  REPLACE WITH THE CODE ITSELF
            ########################################
            ## Check_Interrupt(STOP_ID, OK)
            OK = True  # (over-ridden, for now)
            
            if not(OK):
                self.finalize()
                TF_Print(' ')
                TF_Print('----------------------------------')
                TF_Print(' Simulation terminated by user.')
                TF_Print('----------------------------------')
                return

            self.last_check_time = time.time()

    #   check_interrupt()
    #-------------------------------------------------------------
    def initialize_time_vars(self):

        #------------------
        # Start the clock
        #------------------
        self.start_time = time.time()

        #--------------------------------
        # Initialize the time variables
        #--------------------------------        
        self.time_units = 'seconds'
        self.time_index = np.int32(0)
        self.time       = np.float64(0)
        self.DONE       = False
        
        self.time_sec   = np.float64(0)
        self.time_min   = np.float64(0)

        self.last_check_time  = time.time()  # (for check_interrupt() )
        self.last_print_time  = time.time()  # (for print_time_and_value() )
        self.last_plot_time   = np.float64(0)   # (for update_hydrograph_plot() )

        #---------------------------------------
        # Set the model run timestep to that
        # of the "channel_flow" process [secs]
        #---------------------------------------
        ### self.dt = self.channel_dt  # (5/17/12.  New framework approach.)
        
    #   initialize_time_vars()
    #-------------------------------------------------------------
    def initialize_stop_vars(self):

        TF_Print('Setting stop method to: ' + self.stop_method)

        #---------------------------------------------------------
        # Define some constants for "check_steady_state() method
        #----------------------------------------------------------------
        # Best value of tolerance also depends on the time step.
        # For "plane" case, result changed with timestep = 2 or 4 secs.
        #----------------------------------------------------------------
        # "Optimal" value of steady_tol was found by trial and error.
        #----------------------------------------------------------------        
        self.steady_tol  = np.float64(1E-5)
        self.nonzero_tol = np.float64(1E-5)
        #self.nonzero_tol = 1e-6      #(worked better for "plane" case with step=2s)
        self.n_same      = np.int32(0)
        self.n_same_max  = np.int32(499) # (a number of time steps)

        #-------------------------------------------------------------- 
        # Note: Q_last will be compared to Q_outlet later, so we
        #       must use copy() method for them to ever be different.
        #       Q_outlet is a mutable scalar reference from channels.
        #--------------------------------------------------------------        
        self.Q_last = self.Q_outlet.copy()
        
        if (self.stop_method == 'Q_peak_fraction'):
            #----------------------------------------------------
            # Run until the outlet hydrograph drops to less
            # than "Qp_fraction" of the peak value before that.
            #----------------------------------------------------
            T_stop = 0
            Tstr   = '  [min]'
        elif (self.stop_method == 'Until_model_time'):
            #--------------------------------------------------
            # Run until specified "stopping time", in minutes
            #--------------------------------------------------
            T_stop = self.T_stop_model
            mstr   = ('%10.2f' % T_stop)
            Tstr   = ' of ' + mstr + '  [min]'
        elif (self.stop_method == 'Until_n_steps'):
            #--------------------------------------
            # Run for a specified number of steps
            #--------------------------------------
            n_steps = self.n_steps
            T_stop  = (n_steps * self.dt / np.float64(60))   #[sec] -> [min]
            mstr    = ('%10.2f' % T_stop)
            Tstr   = ' of ' + mstr + '  [min]'
        else:
            print('ERROR: Invalid stopping method.')
            return

        self.T_stop = T_stop
        self.Tstr   = Tstr
            
    #   initialize_stop_vars()       
    #-------------------------------------------------------------
##    def initialize_gui(self):
##        
##        #-----------------------------------------
##        # Set various widget IDs and info
##        # (should this be called by __init__ ??)
##        #-----------------------------------------
##        self.leader_ID  = np.int32(0)
##        self.base_ID    = int32(0)
##        self.start_ID   = int32(0)
##        self.stop_ID    = int32(0)
##        self.draw_ID    = int32(0)
##        self.win_num    = int32(0)   ########
##        self.temp_winID = int32(-1)
##        self.npanels    = npanels
##        self.panel_IDs  = zeros([npanels], dtype='Int32')
##        
##        self.stop_method_ID = int32(0)  #####
##
##        #----------------------------------------
##        # Option to plot hydrograph dynamically
##        #----------------------------------------
##        if (self.PLOT):
##            self.initialize_hydrograph_plot()
##
##    #   initialize_gui()
    #-------------------------------------------------------------
##    def initialize_hydrograph_plot(self):
##        
##        #------------------------------------------------------
##        # NB! Get the window number for the draw widget.
##        #     This assumes that a small draw window has been
##        #     added in the output log panel.  See the function
##        #     called GUI_Stopping in GUI_main.pro.
##        #------------------------------------------------------
##        # NB! Rainbow + white color table was loaded earlier
##        #     by GUI_Stopping, so black=0 and white=255.
##        #------------------------------------------------------
##        ## Initialize_Hydrograph_Window(DRAW_ID, win_num)
##        
##        self.nQ = int32(0)
##        self.nQ_init = int16(1000)
##        nQ_max = self.nQ_init
##        self.tvals = zeros([nQ_max], dtype='Float32')
##        self.Qvals = zeros([nQ_max], dtype='Float32')  
##
##    #   initialize_hydrograph_plot()
    #-------------------------------------------------------------
##    def initialize_mass_totals(self):
## 
##        #-------------------------------------------------------
##        # Prepare to track total mass of each process for the
##        # entire DEM.  This will actually be a volume, in m^3.
##        #-------------------------------------------------------
##        # Volume into and out of main basin is tracked by
##        # data members of the "basins" class
##        #-------------------------------------------------------        
##        self.vol_P  = np.float64(0)
##        self.vol_SM = np.float64(0)
##        self.vol_IN = np.float64(0)
##        self.vol_ET = np.float64(0)
##        self.vol_GW = np.float64(0)
##        self.vol_R  = np.float64(0)
##        self.vol_Rg = np.float64(0)
##
##    #   initialize_mass_totals()       
##    #-------------------------------------------------------------
##    def update_mass_totals(self):
##
##        N_P = N_S = N_I = N_E = N_G = N_R = N_Rg = int32(1)
##        
##        if (np.size(self.pp.P)  == 1): N_P  = self.n_pixels
##        if (np.size(self.sp.SM) == 1): N_S  = self.n_pixels
##        if (np.size(self.ip.IN) == 1): N_I  = self.n_pixels
##        if (np.size(self.ip.Rg) == 1): N_Rg = self.n_pixels
##        if (np.size(self.ep.ET) == 1): N_E  = self.n_pixels
##        if (np.size(self.gp.GW) == 1): N_G  = self.n_pixels
##        if (np.size(self.R)     == 1): N_R  = self.n_pixels
##
##        mfac = self.dt * self.da
##        
##        if (size(self.da) == 1):    
##            self.vol_P  += np.sum(np.double(self.pp.P   * mfac)) * N_P
##            self.vol_SM += np.sum(np.double(self.sp.SM  * mfac)) * N_S
##            self.vol_IN += np.sum(np.double(self.ip.IN  * mfac)) * N_I
##            self.vol_Rg += np.sum(np.double(self.ip.Rg  * mfac)) * N_Rg            
##            self.vol_ET += np.sum(np.double(self.ep.ET  * mfac)) * N_E
##            self.vol_GW += np.sum(np.double(self.gp.GW  * mfac)) * N_G
##            self.vol_R  += np.sum(np.double(self.R      * mfac)) * N_R
##
##        else:    
##            #----------------------------------------
##            # Note:  da is a grid so mfac is a grid
##            #----------------------------------------
##            self.vol_P  += sum(double(self.pp.P   * mfac))
##            self.vol_SM += sum(double(self.sp.SM  * mfac))
##            self.vol_IN += sum(double(self.ip.IN  * mfac))
##            self.vol_Rg += sum(double(self.ip.Rg  * mfac))
##            self.vol_ET += sum(double(self.ep.ET  * mfac))
##            self.vol_GW += sum(double(self.gp.GW  * mfac))
##            self.vol_R  += sum(double(self.R      * mfac))
##
##    #   update_mass_totals
    #-------------------------------------------------------------
    def update_hydrograph_plot(self):
        
        #-----------------------------------------
        # Plot hydrograph up until now (3/22/07)
        #-----------------------------------------
        # plus sign (psym=1), point (psym=3)
        #-----------------------------------------
        elapsed_plot_time = (self.time_min - self.last_plot_time)
        if (self.PLOT and (elapsed_plot_time > 1.0)):    
            #------------------------------------
            # Report an "instantaneous" Q value
            #------------------------------------
            ########################
            nQ      = self.nQ
            nQ_init = self.nQ_init
            tvals   = self.tvals
            Qvals   = self.Qvals
            ########################
            Q_main_out = Pixel_Var(self.cp.Q, self.outlet_ID)   ########
            self.tvals[nQ] = self.time_min
            self.Qvals[nQ] = Q_main_out
            matplotlib.pyplot.figure( self.win_num + 1 )
            matplotlib.pyplot.plot(tvals[0:nQ+1], Qvals[0:nQ+1], \
                                   color='k', marker='.')
            matplotlib.pyplot.axes(axisbg='w')
            matplotlib.pyplot.ylim(np.array(0, np.float64(1.03) * Qvals.max()))
            matplotlib.pyplot.show()  #**** -1
            time.sleep( np.float64(0.005) )  #[seconds]
            nQ = (nQ + np.int32(1))
            if (nQ == self.nQ_max):
                ## Use np.concatenate() here ??
                tvals = array([tvals, np.zeros([nQ_init], dtype='Float32')])
                Qvals = array([Qvals, np.zeros([nQ_init], dtype='Float32')])
                self.nQ_max = np.size(tvals)
            self.last_plot_time = self.time_min
            ########################
            self.nQ     = nQ
            self.tvals  = tvals
            self.Qvals  = Qvals
            ########################

    #   update_hydrograph_plot()
    #-------------------------------------------------------------
    def print_final_report(self, comp_name='TopoFlow',
                           mode='nondriver'):

        #------------------------------------------------------
        # NB! This overrides BMI_base.print_final_report(),
        # so it must have the same arguments. (10/27/11)
        #------------------------------------------------------
        
        #------------------------------------
        # Gather information for the report
        #------------------------------------
        NAME      = self.site_prefix
        COMMENT   = self.comment
        T_stop    = self.T_stop
        Q_final   = self.Q_outlet
        T_final   = self.time_min
        outlet_ID = self.outlet_ID
        n_steps   = self.time_index
        # outlet_col = (outlet_ID % self.nx)   ## (may need this)
        # outlet_row = (outlet_ID / self.nx)

        #---------------------------------------
        # New framework method with 0-d numpy
        # arrays for mutable scalars (2/7/13).
        #---------------------------------------
        Q_peak     = self.Q_peak
        T_peak     = self.T_peak
        u_peak     = self.u_peak
        Tu_peak    = self.Tu_peak
        d_peak     = self.d_peak
        Td_peak    = self.Td_peak
        #-----------------------------
        P_max      = self.P_max

        vol_P  = self.vol_P
        vol_Q  = self.vol_Q 
        vol_SM = self.vol_SM
        vol_MR = self.vol_MR
        vol_ET = self.vol_ET
        vol_IN = self.vol_IN
        vol_Rg = self.vol_Rg
        vol_GW = self.vol_GW
        vol_R  = self.vol_R
        vol_chan = self.vol_chan    # (2019-09-17)
        vol_land = self.vol_land    # (2019-09-17)
        
        #-----------------------------------------------------------------
        # (2/6/13) BMI_base.py has an "initialize_basin_vars()" method
        # that embeds an instance of the "basin" component in the caller
        # as "bp".  So this block of code should still work.
        #-----------------------------------------------------------------
        basin_area = self.basin_area
        ##################################################################
        # We need an "update_basin_vars()" method to be called from each
        # component's update() method in order to compute "volume_in".
        # Then we can set TRACK_VOLUME to True.
        ##################################################################
        ## volume_in  = self.bp.volume_in
        TRACK_VOLUME = False    ##### (since not ready yet) ######
        volume_in    = self.initialize_scalar( 0, dtype='float64')

        #----------------------------
        # Construct run time string
        #----------------------------
        run_time = (time.time() - self.start_time)
        if (run_time > 60):
            run_time = run_time / np.float64(60)
            rt_units = ' [min]'
        else:
            rt_units = ' [sec]'
        run_time_str = str(run_time) + rt_units
        
        if (TRACK_VOLUME):  
            if (volume_in != 0):    
                percent_out = np.float64(100) * vol_Q / volume_in
            else:    
                percent_out = np.float64(0)
  
        #-----------------------------------------------      
        # Prepare to save report as a list of strings
        #-----------------------------------------------
        report = list()  ############# (NEW: 11/15/16)

        #-------------------
        # Build the report
        #-------------------
        hline = ''.ljust(60, '-')
        report.append( hline )
        report.append( TF_Version() )
        report.append( time.asctime() )  #####
        report.append(' ')
        # report.append('Simulated Hydrograph for ' + NAME)
        report.append('Input directory:      ' + self.in_directory)
        report.append('Output directory:     ' + self.out_directory)
        report.append('Site prefix:          ' + self.site_prefix)
        report.append('Case prefix:          ' + self.case_prefix)
        if (COMMENT is not None):
            report.append(' ')
            report.append( COMMENT )
        
        report.append(' ')
        report.append('Simulated time:      ' + str(T_final) + ' [min]')
        report.append('Program run time:    ' + run_time_str)
        report.append(' ')
        report.append('Number of timesteps: ' + str(n_steps))
        report.append('Number of columns:   ' + str(self.nx) )
        report.append('Number of rows:      ' + str(self.ny) )
        report.append(' ')

        #------------------------------------------------------------
        # (2/6/13) With the new framework and use of CSDMS Standard
        # Names there is no simple way for the TopoFlow driver to
        # get the time steps from the other components.  That is,
        # it can request "model__time_step" from the framework, but
        # cannot specify which component(s) that it wants it from.
        # The framework has no trouble getting and printing this
        # info however, and this is now done instead.  But current
        # approach does not write time steps to the log_file.
        #------------------------------------------------------------
##        TF_Print('Channel timestep:    ' + str(cp_dt) + ' [s]')
##        TF_Print('Precip/met timestep: ' + str(mp_dt) + ' [s]')
##        TF_Print('Snowmelt timestep:   ' + str(sp_dt) + ' [s]')
##        TF_Print('ET timestep:         ' + str(ep_dt) + ' [s]')
##        TF_Print('Infiltr. timestep:   ' + str(ip_dt) + ' [s]')
##        TF_Print('Sat. zone timestep:  ' + str(gp_dt) + ' [s]')
##        TF_Print('Icemelt timestep:    ' + str(iip_dt) + ' [s]')        
        # TF_Print('Overland timestep:   ' + str(op_dt) + ' [s]')
        # TF_Print('Sampled timestep:    ' + str(sample_dt) + ' [s]')

        if (self.stop_method == 'Until_model_time'):    
            report.append('T_stop:            ' + str(T_stop) + ' [min]')
            report.append(' ')
        report.append('Main outlet ID:    ' + str(outlet_ID) + ' (row, col)')
        report.append('Basin_area:        ' + str(basin_area) + ' [km^2] ')
        #*** report.append('Basin_length:      ' + TF_String(basin_length) + ' [m]')
        report.append(' ')
            
        if (hasattr(self, 'nval_min')):
            if (self.nval_min != -1):
                report.append("Min Manning's n:   " + str(self.nval_min))
                report.append("Max Manning's n:   " + str(self.nval_max))

        if (hasattr(self, 'z0val_min')):
            if (self.z0val_min != -1):
                report.append("Min z0 value:      " + str(self.z0val_min) + ' [m]')
                report.append("Max z0 value:      " + str(self.z0val_max) + ' [m]')
            
        report.append(' ')
        report.append('Q_final:           ' + str(Q_final) + ' [m^3/s]')
        report.append('Q_peak:            ' + str(Q_peak)  + ' [m^3/s]')
        report.append('Q_peak_time:       ' + str(T_peak)  + ' [min]')
        report.append('u_peak:            ' + str(u_peak)  + ' [m/s]')
        report.append('u_peak_time:       ' + str(Tu_peak) + ' [min]')
        report.append('d_peak:            ' + str(d_peak)  + ' [m]')
        report.append('d_peak_time:       ' + str(Td_peak) + ' [min]')
        report.append(' ')
        
        ##############################################################################
        if (TRACK_VOLUME):    
            report.append('Total volume out:  ' + str(vol_Q) + ' [m^3]')
            if (basin_area != 0):    
                report.append('Rain volume in:    ' + str(volume_in)   + ' [m^3]')
                report.append('Percentage out:    ' + str(percent_out) + ' [%]')
                report.append(' ')
        ##############################################################################
                
        #--------------------------------
        # Print the maximum precip rate
        #--------------------------------
        MPR = (P_max * self.mps_to_mmph)   # ([m/s] -> [mm/hr])
        report.append('Max(precip rate):  ' + str(MPR) + ' [mm/hr]')
        report.append(' ')

        #--------------------------------------------
        # Print the area_time integrals over domain
        #--------------------------------------------
        report.append('Total accumulated volumes over entire DEM:')
        report.append('vol_P    (rainfall):      ' + str(vol_P)  + ' [m^3]   (snowfall excluded)')
        report.append('vol_Q    (discharge):     ' + str(vol_Q)  + ' [m^3]   (main basin outlet)')
        report.append('vol_SM   (snowmelt):      ' + str(vol_SM) + ' [m^3]')
        report.append('vol_MR   (icemelt):       ' + str(vol_MR) + ' [m^3]')
        report.append('vol_ET   (evaporation):   ' + str(vol_ET) + ' [m^3]')
        report.append('vol_IN   (infiltration):  ' + str(vol_IN) + ' [m^3]')
        report.append('vol_Rg   (recharge):      ' + str(vol_Rg) + ' [m^3]   (to water table)')
        report.append('vol_GW   (baseflow):      ' + str(vol_GW) + ' [m^3]')
        report.append('vol_R    (runoff):        ' + str(vol_R)  + ' [m^3]   R = (P+SM+MR) - (ET+IN)')
        report.append('vol_chan (channels):      ' + str(vol_chan) + ' [m^3]')
        report.append('vol_land (surface):       ' + str(vol_land) + ' [m^3]')
        report.append(' ')

        #----------------------------------        
        # Print the report to the console
        #----------------------------------
        for line in report:
            print(line)

        #----------------------------------
        # Print the report to a logfile ?
        #----------------------------------
        if (self.WRITE_LOG):
            for line in report:
                self.log_unit.write( line + '\n' )

        self.print_mins_and_maxes( FINAL=True )

        if (self.WRITE_LOG):
            #------------------------------------------------
            # This line is printed to console in finalize()
            #------------------------------------------------
            self.log_unit.write( 'Finished. (' + self.case_prefix + ')\n' )
            self.log_unit.write( '\n' )
      
    #   print_final_report()
    #-------------------------------------------------------------
    def print_mins_and_maxes(self, FINAL=False):

        #-------------------------------
        # New framework method, 2/6/13
        #-------------------------------
        Q_min = self.Q_min
        Q_max = self.Q_max
        u_min = self.u_min
        u_max = self.u_max
        d_min = self.d_min
        d_max = self.d_max
        
        f1 = '(F14.4)'   #(2/12/07)
        
        Qstr = TF_String( Q_min, FORMAT=f1 )
        Qstr = Qstr + ', ' + TF_String( Q_max, FORMAT=f1 )
        #----------------------------------------------------
        ustr = TF_String( u_min, FORMAT=f1 )
        ustr = ustr + ', ' + TF_String( u_max, FORMAT=f1 )
        #----------------------------------------------------
        dstr = TF_String( d_min, FORMAT=f1 )
        dstr = dstr + ', ' + TF_String( d_max, FORMAT=f1 )

        #-----------------------------------------------      
        # Prepare to save report as a list of strings
        #-----------------------------------------------
        report = list()  ############# (NEW: 11/15/16)
        
        if (FINAL):    
            report.append('Final grid mins and maxes:')
        else:    
            report.append('------------------------------------------')
        report.append('Min(Q), Max(Q):   ' + Qstr + ' [m^3/s]')
        report.append('Min(u), Max(u):   ' + ustr + ' [m/s]')
        report.append('Min(d), Max(d):   ' + dstr + ' [m]')
        report.append(' ')

        #----------------------------------        
        # Print the report to the console
        #----------------------------------
        for line in report:
            print(line)

        #----------------------------------
        # Print the report to a logfile ?
        #----------------------------------
        if (self.WRITE_LOG):
            for line in report:
                self.log_unit.write( line + '\n' )
        
    #   print_mins_and_maxes()
    #-------------------------------------------------------------
    def print_uniform_precip_data(self):

        ## precip_method = self.pp.method
        precip_method = 2
        
        #---------------------------------------------------
        # Precip method 1 is special and spatially uniform
        #---------------------------------------------------
        if (precip_method == 1):
            rates = (self.pp.method1_rates * self.mps_to_mmph)    #[m/s] -> [mm/hr]
            durs  = self.pp.method1_durations      
            nr    = np.size(rates)

            if (nr > 1): 
                rstr = str(rates[0])
                dstr = str(durs[0])
                for m in range(1, nr):
                    rstr += ('  ' + str(rates[m]))
                    dstr += ('  ' + str(durs[m]))
            else:
                rstr = TF_String(rates)
                dstr = TF_String(durs)
        elif (precip_method == 2) and \
             (self.pp.rate_type == 0) and \
             (self.pp.duration_type == 0):
                rstr = TF_String(self.pp.rate * self.mps_to_mmph)
                dstr = TF_String(self.pp.duration)
        else:
            #------------------------------------------------
            # Could have uniform precip with method 2 where
            # values are stored in a text file and read in
            # one by one.  Could read that file here.
            #------------------------------------------------
            return
            
        #-------------------------------------
        # This is too verbose in most cases?
        #------------------------------------- 
        #TF_Print,'Uniform precip. rate information: '
        #TF_Print,'Precip. rate:     ' + rstr + ' [mm/hr]'
        #TF_Print,'Duration:         ' + dstr + ' [min]'
        #TF_Print,' '
        
        #----------------------
        # Write to log file ?
        #----------------------
        if (self.WRITE_LOG):
            log_unit = self.log_unit
            log_unit.write("Uniform precip. rate information: \n")
            log_unit.write('Precip. rate:     ' + rstr + " [mm/hr]\n")
            log_unit.write('Duration:         ' + dstr + " [min]\n")
            log_unit.write("\n")
        
    #   print_uniform_precip_data()
    #-------------------------------------------------------------
    def print_dimless_number_data(self, basin_length):
 
##        rates     = self.pp.rates
##        durations = self.pp.durations     
##        rates     = self.mp.P
##        durations = self.mp.dt

        ###########################################
        # Updated for new framework (2/5/13),
        # but still need dt@met for durations !!
        ###########################################
        rates     = self.P
        durations = self.dt
        
        nd        = np.size(durations)
        T_peak    = self.T_peak
        log_unit  = self.log_unit
        
        #-------------------------------------
        # Compute some dimensionless numbers
        # that characterize the hydrograph
        #-------------------------------------
        TAU_P = (self.T_peak / durations)
        TAU_F = (self.T_final / durations)
        PSI   = (basin_length / (rates[0] * durations[0] * np.float64(60)))
        
        tpstr = str(TAU_P[0])
        tfstr = str(TAU_F[0])
        for m in range(1, nd):
            tpstr = (tpstr + '  ' + str(TAU_P[m]))
            tfstr = (tfstr + '  ' + str(TAU_F[m]))
        
        #----------------------------------------
        # Make predictions for Q_peak and T_peak
        #-----------------------------------------
        Q_peak_pred = (np.float64(0.2) * rates[0] * (basin_length) ** np.float64(2) / np.float64(3))
        T_peak_pred = (np.float64(3) * durations[0])
        
        TF_Print('Dimensionless number information:\n')
        TF_Print('T_peak /Duration: ' + tpstr + "\n")
        TF_Print('T_final/Duration: ' + tfstr + "\n") 
        vTF_Print('Psi=L/(R*TD):     ' + str(PSI) + ' [unitless]\n')
        TF_Print(' ')
        TF_Print('Q_peak predicted: ' + str(Q_peak_pred) + ' [m^3/s]\n')
        TF_Print('T_peak predicted: ' + str(T_peak_pred) + ' (min]\n')
        TF_Print(' ')
        
        #----------------------
        # Write to log file ?
        #----------------------
        if (self.WRITE_LOG):    
            self.log_unit.write('Dimensionless number information:\n')   
            self.log_unit.write('T_peak /Duration: ' + tpstr + "\n")
            self.log_unit.write('T_final/Duration: ' + tfstr + "\n")
            self.log_unit.write("\n")
            #printf,LU,'Psi=L/(R*TD):     ' + strPSI) + ' [unitless]'
            #printf,LU,' '
            #printf,LU,'Q_peak predicted: ' + str(Q_peak_pred) + ' [m^3/s]'
            #printf,LU,'T_peak predicted: ' + str(T_peak_pred) + ' [min]'
            #printf,LU,' '
                  
    #   print_dimless_number_data()
    #-------------------------------------------------------------
    def print_mass_balance_report(self):            

        #--------------------------------------
        # Updated for new framework. (2/5/13)
        #--------------------------------------
        vol_P  = self.vol_P
        vol_SM = self.vol_SM
        vol_IN = self.vol_IN
        vol_Rg = self.vol_Rg
        vol_ET = self.vol_ET
        vol_GW = self.vol_GW
        vol_R  = self.vol_R
        vol_chan = self.vol_chan    # (2019-09-17)
        
        TF_Print('Volume rain         = ' + str(vol_P)  + ' [m^3]')
        TF_Print('Volume snowmelt     = ' + str(vol_SM) + ' [m^3]')
        TF_Print('Volume infiltrated  = ' + str(vol_IN) + ' [m^3]')
        TF_Print('Volume evaporated   = ' + str(vol_ET) + ' [m^3]')
        TF_Print('Volume seepage      = ' + str(vol_GW) + ' [m^3]')
        TF_Print('Volume runoff       = ' + str(vol_R)  + ' [m^3]')
        TF_Print('Volume bottom loss  = ' + str(vol_Rg) + ' [m^3]')
        TF_Print('Volume all channels = ' + str(vol_chan) + ' [m^3]')
        TF_Print(' ')

        RICHARDS = self.ip.get_scalar_long('RICHARDS')
        print('In topoflow.print_mass_balance_report(),')
        print('   RICHARDS =', RICHARDS)
        print(' ')
        if (RICHARDS):    
            #----------------------------------------
            # Richards' equation for infiltration
            # Ground water losses not counted above
            #----------------------------------------
            da         = self.da
            n_pixels   = self.n_pixels   #####
            vol_stored = np.float64(0)
            SCALAR_dz  = (np.size( self.ip.dz ) == 1)
            SCALAR_da  = (np.size(da) == 1)     #(assume da is scalar or grid)
            #-----------------------------------
            dim_q  = np.ndim( self.ip.q )
            dim_qi = np.ndim( self.ip.qi )
            #----------------------------------
            for k in range(self.ip.nz):
                #--------------------------------------
                # In this loop, dz is always a scalar
                #--------------------------------------
                if (SCALAR_dz):    
                    dz = self.ip.dz
                else:    
                    dz = self.ip.dz[k]
                if (dim_q == 3):    
                    qq = self.ip.q[k,:,:]
                else:    
                    qq = self.ip.q[k]
                if (dim_qi == 3):    
                    qqi = self.ip.qi[k,:,:]
                else:    
                    qqi = self.ip.qi[k]
                #----------------------------------------
                # Get increase in q over initial values
                # (but q may have decreased due to ET)
                #----------------------------------------
                dq = (qq - qqi)   #(grid or scalar)
                da
                SCALAR_dq = (np.size(dq) == 1)
                if (SCALAR_da and SCALAR_dq):    
                    dm = dq * (dz * da * n_pixels)      #(da is a SCALAR)
                else:    
                    dm = np.sum(np.double(dq * (da * dz)))    #(da is a GRID)
                vol_stored += dm
            mass_error = np.float64(100) * (vol_IN - vol_stored) / vol_IN
            err_str = ('%6.2f' % mass_error) + ' % '
            TF_Print('Volume stored      = ' + TF_String(vol_stored) + ' [m^3]')
            TF_Print('Volume error       = ' + err_str)
        TF_Print(' ')

    #   print_mass_balance_report()
    #-------------------------------------------------------------

    


