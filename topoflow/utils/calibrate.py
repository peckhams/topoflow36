
#------------------------------------------------------------------------
#  Copyright (c) 2020, Scott D. Peckham
#
#  Jun 2020.  Modified for use with remote-sensing data.
#  May 2020.  Created for TopoFlow calibration notebook.
#
#---------------------------------------------------------------------
#
#  test()
#
#  class calibrator():
#      __init__()
#      read_cfg_file()
#      set_file_info()
#      print_file_info()
#      print_cfg_info()
#
#      read_time_series()
#      regularize_obs_time_series()
#
#      get_observed_values()
#      get_simulated_values()
#      modify_cfg_file()
#      backup_input_file()
#      modify_channel_grid_file()
#      compute_cost()
#      get_parameter_values()
#      calibrate()
#
#      plot_time_series()
#      update_sim_data_file()
#      delete_output_files()    ## commented out; dangerous
#
#--------------------------------------------
#  These functions are outside of the class
#--------------------------------------------
#  read_any_time_series()
#  compute_rain_rates_from_accumulations()
#  spearman_metric()
#
#---------------------------------------------------------------------

import glob, os, os.path, shutil, time

import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import spearmanr

from topoflow.framework import emeli  # (for examples dir)
from topoflow.utils import parameterize
from topoflow.utils import time_utils as tu
from topoflow.utils import visualize as tfvis
from topoflow import main

#------------------------------------------------------------------------
def test():

    c = calibrator()
    c.print_file_info()
    c.get_simulated_values( SILENT=False )
    c.get_observed_values( SILENT=False  )
    print('obs_times[0:20]  =', c.obs_times[0:8])
    print('obs_values[0:20] =', c.obs_values[0:8])
    print('sim_times[0:20]  =', c.sim_times[0:8])
    print('sim_values[0:20] =', c.sim_values[0:8])
    print()
    #--------------------------------------------------
    cost = c.compute_cost( SILENT=False )

#   test()
#------------------------------------------------------------------------
class calibrator:

    def __init__(self, cfg_file=None):

        self.version       = '0.5'
        self.cfg_file      = cfg_file
        self.set_file_info()

    #   __init__()
    #---------------------------------------------------------------------
    def read_cfg_file(self, delim='='):

        if not(os.path.exists( self.cfg_file ) ):
            print('Sorry, Could not find CFG file:')
            print('   ' + self.cfg_file )
            print('Provide full path or place in current directory.')
            print()
            return

        #---------------------------------------------
        # Save variables into self using key & value
        #---------------------------------------------
        cfg_unit = open(self.cfg_file, 'r')
        while (True):
            line = cfg_unit.readline()
            if (line == ''):
                break
            ## print( line )
            if (line[0] != '#') and (delim in line):
                parts = line.split( delim )
                key   = parts[0].strip()   # string
                value = parts[1].strip()   # string
                try:
                    value = int(value)   #######
                except:
                    pass   
                #-------------------------------------
                if (value == 'True'):  value = True
                if (value == 'False'): value = False
                if (value == 'None'):  value = None
                #-------------------------------------
                exec( "self." + key + " = value", {}, locals() )

        #--------------------------------------------        
        # Add a path separator at end of dir names?
        #--------------------------------------------
        if (self.basin_dir[-1] != os.sep):  
            self.basin_dir += os.sep
        if (self.topo_dir[-1] != os.sep):  
            self.topo_dir += os.sep
        if (self.cfg_dir[-1] != os.sep):  
            self.cfg_dir += os.sep
        if (self.output_dir[-1] != os.sep):  
            self.output_dir += os.sep
        if (self.obs_dir[-1] != os.sep):  
            self.obs_dir += os.sep
        
        #----------------------------------                   
        # These are built from CFG values
        #----------------------------------
        self.sim_dir       = self.output_dir
        self.sim_data_file = (self.sim_dir + self.cfg_prefix + '_0D-Q.txt')

        #----------------------------------
        # Standardize the datetime string
        #----------------------------------
        datetime_str = self.obs_start_datetime
        self.obs_start_datetime = tu.standardize_datetime_str( datetime_str )
        
    #   read_cfg_file()
    #---------------------------------------------------------------------
    def set_file_info(self):
 
        self.home_dir     = os.path.expanduser("~")
        self.examples_dir = emeli.paths['examples']
        
        if (self.cfg_file is not None):
            self.read_cfg_file()
            self.obs_data_file = (self.obs_dir + self.obs_data_file)
            return

        #----------------------------------------------
        # The remaining values use test data defaults
        #----------------------------------------------
        self.cost_fcn_name = 'Lp_norm'  # (or 'spearman')
        self.p_for_Lp_norm = 2.0        # L2 norm is default
        
        #----------------------------------
        # Remaining settings are defaults
        # for the Treynor example dataset
        #----------------------------------
        # Set cfg_prefix and directories
        #----------------------------------
        self.site_prefix  = 'Treynor'        
        self.cfg_prefix   = 'June_20_67'
        self.basin_dir    = self.examples_dir  + 'Treynor_Iowa_30m/'
        self.topo_dir     = self.basin_dir + '__topo/'
        self.cfg_dir      = self.basin_dir + '__No_Infil_June_20_67_rain/'
        self.output_dir   = self.home_dir + '/TF_Output/Treynor/'

        #---------------------------------------
        # Attributes of the observed data file
        #---------------------------------------
        self.obs_time_format    = 'hhmm'
        self.obs_time_irregular = True
        self.obs_time_units     = 'minutes'   ####################
        self.obs_header_lines   = 6
        self.obs_time_column    = 0    # (HHMM)
        self.obs_value_column   = 5
        self.obs_interp_method  = 'linear'  # or 'nearest', etc.
        #---------------------------------------------------
        obs_dir = self.basin_dir + '__observations/'
        self.obs_dir       = obs_dir
        self.obs_data_file = obs_dir + 'June_20_1967_Observed_Discharge.txt'
    
        #----------------------------------------        
        # Attributes of the simulated data file
        #-------------------------------------------------
        # obs is multi-column, sim is single-column now.
        #-------------------------------------------------
        self.sim_time_format    = 'timesince'
        self.sim_time_irregular = False
        self.sim_time_units     = 'minutes'   ####################
        self.sim_header_lines = 2
        self.sim_time_column  = 0    # [min]
        self.sim_value_column = 1
        #---------------------------------------------------
        # self.sim_data_file will be replaced each time
        self.sim_dir       = self.output_dir
        self.sim_data_file = self.sim_dir + (self.cfg_prefix + '_0D-Q.txt')
              
    #   set_file_info()
    #---------------------------------------------------------------------
    def print_file_info(self):
 
        #-----------------------------------
        # Print cfg_prefix and directories
        #-----------------------------------
        print('In print_file_info():')
        print('  site_prefix  =', self.site_prefix )
        print('  cfg_prefix   =', self.cfg_prefix )
        print('  home_dir     =', self.home_dir )
        ## print('  examples_dir =', self.examples_dir )
        print('  basin_dir    =', self.basin_dir )
        print('  topo_dir     =', self.topo_dir )
        print('  cfg_dir      =', self.cfg_dir )
        print('  output_dir   =', self.output_dir )
        print('  obs_dir      =', self.obs_dir )
        print('  sim_dir      =', self.sim_dir )
        print()

    #   print_file_info()
    #---------------------------------------------------------------------
    def print_cfg_info(self):
 
        #-------------------------------------------
        # Print information read from cal CFG file
        #-------------------------------------------
        print('In print_cfg_info():')
        print('site_prefix  =', self.site_prefix )
        print('cfg_prefix   =', self.cfg_prefix )
        print('home_dir     =', self.home_dir )
        print('basin_dir    =', self.basin_dir )
        print('topo_dir     =', self.topo_dir )
        print('cfg_dir      =', self.cfg_dir )
        print('output_dir   =', self.output_dir )
        print('------------------------------------------------------')
        print('obs_dir            =', self.obs_dir )
        print('obs_data_file      =', self.obs_data_file )
        print('obs_data_delim     =', self.obs_data_delim )
        print('obs_time_format    =', self.obs_time_format )
        print('obs_time_irregular =', self.obs_time_irregular )
        print('obs_time_units     =', self.obs_time_units )
        print('obs_start_datetime =', self.obs_start_datetime )
        print('obs_header_lines   =', self.obs_header_lines )
        print('obs_time_column    =', self.obs_time_column )
        print('obs_value_column   =', self.obs_value_column )
        print('obs_interp_method  =', self.obs_interp_method )
        print('------------------------------------------------------')
        print('sim_dir            =', self.sim_dir )
        print('sim_data_file      =', self.sim_data_file )
        print('sim_data_delim     =', self.sim_data_delim )
        print('sim_time_format    =', self.sim_time_format )
        print('sim_time_irregular =', self.sim_time_irregular )
        print('sim_time_units     =', self.sim_time_units )
        print('sim_start_datetime =', self.sim_start_datetime )
        print('sim_header_lines   =', self.sim_header_lines )
        print('sim_time_column    =', self.sim_time_column )
        print('sim_value_column   =', self.sim_value_column )
        print('------------------------------------------------------')
        print('cost_fcn_name      =', self.cost_fcn_name )
        print('p_for_Lp_norm      =', self.p_for_Lp_norm )
        print()

    #   print_cfg_info()      
    #---------------------------------------------------------------------
    def read_time_series(self, source='observed', SILENT=True,
                         time_type='float'):

        #---------------------------------------------
        # Note:  source = 'observed' or 'simulated'.
        #        time_type = 'float' or 'string'
        #---------------------------------------------
        if (source == 'observed'):
            txt_file  = self.obs_data_file
            delim     = self.obs_data_delim
            hdr_lines = self.obs_header_lines
            time_col  = self.obs_time_column
            value_col = self.obs_value_column
        else:     ## 'simulated'
            txt_file  = self.sim_data_file
            delim     = self.sim_data_delim
            hdr_lines = self.sim_header_lines
            time_col  = self.sim_time_column
            value_col = self.sim_value_column
            time_type = 'float'  # vs. string

        file_unit = open(txt_file, 'r')

        #--------------------    
        # Skip header lines
        #--------------------
        for k in range(hdr_lines):
            line = file_unit.readline()
        
        #------------------------    
        # Read times and values
        #----------------------------------------------
        # Note:  If (delim is None), "split()" method
        #        splits on any whitespace.
        #----------------------------------------------
        times  = list()
        values = list()
        while (True):
            line  = file_unit.readline()
            if (line == ''):
                break
            words = line.split( delim )   # None implies any whitespace
            ncols = len( words )
            if (time_col > ncols-1) or (value_col > ncols-1):
                print('ERROR: Number of columns =', ncols)
                print('   but time_col  =', time_col )
                print('   and value_col =', value_col )
                break
            if (time_type == 'float'):
                time  = np.float32( words[ time_col ] )
            else:
                time = words[ time_col ]  ## (string)
            value = np.float32( words[ value_col ] )
            times.append( time )
            values.append( value )
        #-------------------------------
        file_unit.close()
        times  = np.array( times )
        values = np.array( values )
    
        if (not(SILENT) and (time_type == 'float')):
            print('In read_time_series():')
            print('Min(times)  =', times.min())
            print('Max(times)  =', times.max())
            print('Min(values) =', values.min())
            print('Max(values) =', values.max())
            print()

        if (source == 'observed'):
            self.obs_times  = times
            self.obs_values = values
            #-----------------------------------------------       
            # Make a copy of the original info, because we
            # may apply conversion or regularize later.
            #-----------------------------------------------
            self.obs_orig_times  = times.copy()
            self.obs_orig_values = values.copy()
        else:
            self.sim_times  = times
            self.sim_values = values
                 
    #   read_time_series()
    #---------------------------------------------------------------------
    def regularize_obs_time_series( self, SILENT=True ): 

        #----------------------------------------------------------
        # Notes:  This function is for the special case where the
        #         time column contains times in hhmm format that
        #         are not evenly spaced in time.  It returns a
        #         new time series (decimal hours) with values at
        #         self.sim_times (whether regular or irregular).
        #----------------------------------------------------------
        # e.g. interp_method = linear, nearest, etc.
        #----------------------------------------------------------
        # The read_time_series() method saves a copy of the
        # original times and values, as:
        #    self.obs_orig_times and self.obs_orig_values
        # since we may apply conversion or regularize.
        #----------------------------------------------------------       
        times  = self.obs_times
        values = self.obs_values

        #----------------------------------------------             
        # First, define an interpolation function, f1
        #----------------------------------------------
        # f1 = interp1d( times, values, kind='nearest')
        # f1 = interp1d( times, values, kind=self.obs_interp_method)
        f1 = interp1d( times, values, kind=self.obs_interp_method,
                       fill_value='extrapolate')
 
        #--------------------------------------------------        
        # Note:  An alternative to performing a model run
        #        before we can regularize observed values
        #        is to use sim_time_interval like this.
        #        Could put sim_time_interval in CFG file.
        #--------------------------------------------------
        # dt   = self.sim_time_interval            
        # tmin = times.min()
        # tmax = times.max()
        # nt   = 1 + (tmax - tmin) / dt
        # new_times = (dt * np.arange( nt )) + tmin
        #------------------------------------------------------
        # Notes:  We currently regularize the observed values
        #         to match the times of the simulated values.
        #         This requires running the first simulation
        #         before we retrieve the observed values.
        #------------------------------------------------------
        #         This may also result in discarding some of
        #         the observed values.
        #------------------------------------------------------        
        new_times = self.sim_times 
        #-----------------------------
        new_values = f1( new_times )
    
        if not(SILENT):
            print('In regularize_obs_time_series():')
            print('   Orig number of times =', len(times) )
            print('   New number of times  =', len(new_times) )
            print()

        self.obs_times  = new_times
        self.obs_values = new_values

    #   regularize_obs_time_series()
    #---------------------------------------------------------------------
    def get_observed_values(self, SILENT=True):

        #----------------------------------------------------
        # For example:
        # Get the observed discharge data for Treynor Creek
        # for the extreme rainfall event of June 20, 1967.
        #----------------------------------------------------
        ## txt_file = 'June_20_1967_Observed_Discharge.txt'

        #--------------------------------------        
        # Set self.obs_times, self.obs_values
        #--------------------------------------
        if not(SILENT):
            print('Reading time series data...')
        otf = self.obs_time_format
        if (otf == 'hhmm'):
            self.read_time_series( source='observed', time_type='str' )
            times_hhmm = self.obs_times
            self.obs_times = tu.convert_times_from_hhmm_to_minutes( times_hhmm )
            ## self.convert_times_from_hhmm_to_minutes()
        elif (otf == 'datetime'):
            self.read_time_series( source='observed', time_type='str' )
            times_datetime = self.obs_times
            origin_datetime_str = self.obs_start_datetime  ################
            origin_datetime_obj = tu.get_datetime_obj_from_one_str( origin_datetime_str )
            self.obs_times = tu.convert_times_from_datetime_to_minutes( times_datetime,
                                        origin_datetime_obj )
            ## self.convert_times_from_datetime_to_minutes()        
        else:
            self.read_time_series( source='observed' )
                         
        #------------------------------------------------------  
        # Use 1D interpolation to create even spacing in time
        #------------------------------------------------------
        if (self.obs_time_irregular):
            if not(SILENT):
                print('Regularizing time series data...')
                print()
            self.regularize_obs_time_series()
            #--------------------------------------------------
            # Returning to hhmm-format times isn't necessary.
            # It depends on what we want to do next with the
            # times, like using them for x-axis of a plot.
            # They are not used to compute the cost function.
            #--------------------------------------------------
            times_min = self.obs_times
            # For testing
            ## print('###### times_min[0:10] =', times_min[0:10])
            ## print('###### times_min[-10:] =', times_min[-10:])
            if (otf == 'hhmm'):
                self.obs_times = tu.convert_times_from_minutes_to_hhmm( times_min )
                ## self.convert_times_from_minutes_to_hhmm()
            elif (otf == 'datetime'):
                origin_datetime_str = self.obs_start_datetime  ################
                origin_datetime_obj = tu.get_datetime_obj_from_one_str( origin_datetime_str )
                self.obs_times = tu.convert_times_from_minutes_to_datetime( times_min,
                                            origin_datetime_obj )
                ## self.convert_times_from_minutes_to_datetime()
        if not(SILENT):
            print('In get_observed_values():')
            print('Obs time format    =', self.obs_time_format )
            print('Obs time irregular =', self.obs_time_irregular )
            if (otf not in ['hhmm', 'datetime']):
                print('Min(obs_times)     =', self.obs_times.min() )
                print('Max(obs_times)     =', self.obs_times.max() )
            print('Min(obs_values)    =', self.obs_values.min() )
            print('Max(obs_values)    =', self.obs_values.max() )
            print()

    #   get_observed_values()
    #---------------------------------------------------------------------
    def get_simulated_values(self, SILENT=True ):

        #----------------------------------------------------
        # For example:
        # Get the observed discharge data for Treynor Creek
        # for the extreme rainfall event of June 20, 1967.
        #----------------------------------------------------
        # e.g. sim_data_file = (case_prefix + '_0D-Q.txt')
        
        #--------------------------------------        
        # Set self.sim_times, self.sim_values
        #--------------------------------------
        self.read_time_series( source='simulated' )
 
        if not(SILENT):
            print('In get_simulated_values():')
            print('Sim time format    =', self.sim_time_format )
            print('Sim time irregular =', self.sim_time_irregular)
            print('Min(sim_times)     =', self.sim_times.min() )
            print('Max(sim_times)     =', self.sim_times.max() )
            print('Min(sim_values)    =', self.sim_values.min() )
            print('Max(sim_values)    =', self.sim_values.max() )
            print()
             
    #   get_simulated_values()
    #---------------------------------------------------------------------
    def modify_cfg_file(self, p):

        # Or maybe modify_input_file(), like "chan-n.rtg"
        pass
 
    #   modify_cfg_file()
    #---------------------------------------------------------------------
    def backup_input_file(self, file_type='manning'):  

        topo_dir     = self.topo_dir
        site_prefix  = self.site_prefix
        cfg_dir      = self.cfg_dir
        
        if (file_type == 'manning'):
            input_file  = topo_dir + site_prefix + '_chan-n.rtg'
            backup_file = topo_dir + site_prefix + '_chan-n.ORIG'
        elif (file_type == 'width'):
            input_file  = topo_dir + site_prefix + '_chan-w.rtg'   
            backup_file = topo_dir + site_prefix + '_chan-w.ORIG'
                             
        shutil.copyfile( input_file, backup_file)

    #   backup_input_file()
    #---------------------------------------------------------------------
    def modify_channel_grid_file(self, cal_var, param, SILENT=True ):    

        REPORT       = not(SILENT)
        topo_dir     = self.topo_dir
        site_prefix  = self.site_prefix
        ## d8_area_file = site_prefix + '_area.rtg'
        d8_area_file = site_prefix + '_d8-area.rtg'
        ### d8_area_file = topo_dir + site_prefix + '_area.rtg'

        #----------------------------------------------------       
        # Note:  If g1 and g2 are given, it is assumed that
        #        g1 = c * Amax^p and g2 = c * Amin^p 
        #----------------------------------------------------
        if (cal_var == 'channel_width_power'):
            width_file = site_prefix + '_chan-w.rtg'
            parameterize.get_grid_from_TCA(site_prefix=site_prefix,
                         topo_dir=topo_dir, area_file=d8_area_file,
                         out_file=width_file, REPORT=REPORT,
                         g1=self.channel_width_max, p=param)
        if (cal_var == 'channel_width_max'):
            width_file = site_prefix + '_chan-w.rtg'
            parameterize.get_grid_from_TCA(site_prefix=site_prefix,
                         topo_dir=topo_dir, area_file=d8_area_file,
                         out_file=width_file, REPORT=REPORT,
                         g1=param, p=self.channel_width_power)                                 
        if (cal_var == 'manning_n_min'):
            manning_file = site_prefix + '_chan-n.rtg' 
            parameterize.get_grid_from_TCA(site_prefix=site_prefix,
                         topo_dir=topo_dir, area_file=d8_area_file,
                         out_file=manning_file, REPORT=REPORT,
                         g1=param, g2=self.manning_n_max )
        if (cal_var == 'manning_n_max'):
            manning_file = site_prefix + '_chan-n.rtg' 
            parameterize.get_grid_from_TCA(site_prefix=site_prefix,
                         topo_dir=topo_dir, area_file=d8_area_file,
                         out_file=manning_file, REPORT=REPORT,
                         g1=self.manning_n_max, g2=param )

    #   modify_channel_grid_file()
    #---------------------------------------------------------------------
    def compute_cost(self, SILENT=True ):

        Y_obs = self.obs_values
        Y_sim = self.sim_values

        #-----------------------------------        
        # Check if arrays are equal length
        #-----------------------------------
        n_obs = Y_obs.size   #  ndarrays
        n_sim = Y_sim.size
        if (n_obs != n_sim):
            if not(SILENT):
                print('WARNING: Y_obs and Y_sim have different size.')
                print('n_obs =', n_obs, 'and', 'n_sim =', n_sim )
                print('Will compare up to size of smallest.')
                print()
            if (n_obs < n_sim):
                Y_sim = Y_sim[0:n_obs]
            else:
                Y_obs = Y_obs[0:n_sim]

        #-----------------------------------------
        # Compute cost with chosen cost function
        #-----------------------------------------
        if (self.cost_fcn_name == 'Lp_norm'):
            #--------------------------------
            # Lp norm with array operations
            #--------------------------------
            p    = self.p_for_Lp_norm
            arg  = np.abs( Y_obs - Y_sim )**p
            cost = (arg.sum())**(1.0/p)
        elif (self.cost_fcn_name == 'spearman'):
            cost = spearman_metric(Y_obs, Y_sim)
    
        if not(SILENT):
            print('cost =', cost)
            print()

        return cost
    
    #   compute_cost()
    #---------------------------------------------------------------------
    def get_parameter_values(self, range=[0.1, 1.0], n=10 ): 
        
        vmin = np.float32( range[0] )
        vmax = np.float32( range[1] )
        ramp = np.arange(n, dtype='float32')/(n - 1) # in [0,1]
        p_values = (vmax - vmin)*ramp + vmin
        
        return p_values

    #   get_parameter_values()
    #---------------------------------------------------------------------
    def calibrate(self, PLOT=True, SILENT=True,
                  cal_var='channel_width_power' ):
                  ### cal_var='channel_width_max' ):
                  ### cal_var='manning_n_min' ):
                  ### cal_var='manning_n_max' ):

        #---------------------------------------
        # Create backup of grid to be modified
        #---------------------------------------
        if ('width' in cal_var):
            file_type = 'width'
        elif ('manning' in cal_var):
            file_type = 'manning'
        self.backup_input_file( file_type=file_type )
 
        #---------------------------   
        # Default parameter values
        # (unless made to vary)
        #--------------------------- 
        self.channel_width_power = 0.5
        self.channel_width_max   = 3.0  # [meters]
        self.manning_n_min       = 0.03
        self.manning_n_max       = 0.3
        # self.log_law_z0_min      =
        # self.log_law_z0_max      = 
        
        #--------------------   
        # Initialize values
        #--------------------
        index    = 0
        min_cost = 1e6
        tf_time_interp_method = 'None'
        ## tf_time_interp_method = 'Linear'
        
        if (cal_var == 'channel_width_power'):
            range = [0.2, 1.0]
        if (cal_var == 'channel_width_max'):
            range = [1.0, 10.0]
        if (cal_var == 'manning_n_min'):
            ## range = [0.04, 0.05]
            range = [0.02, 0.16]
        if (cal_var == 'manning_n_max'):
            range = [0.1, 0.9]
        npv = 15
        p_values = self.get_parameter_values( range=range, n=npv)
        costs = np.zeros( npv )

        #-------------------------------------------      
        # Run the model repeatedly, vary parameter
        #-------------------------------------------
        print('Working...')
        for param in p_values:
            print('   param =', param )
            ## self.modify_cfg_file( p )
            self.modify_channel_grid_file( cal_var, param )
            time.sleep( 1.0 )  ###### IS THIS NEEDED?
            ############################################
            #------------------------------------------------
            main.run_model(cfg_prefix=self.cfg_prefix,
                           cfg_directory=self.cfg_dir,
                           SILENT=True,  ### (always for calibration)
                           time_interp_method=tf_time_interp_method)

            #---------------------------------
            # Update self.sim_data_file
            # Files are not overwritten but
            # are appended with a number.
            #---------------------------------
            self.update_sim_data_file()
            self.get_simulated_values( SILENT=SILENT )
            
            #--------------------------------
            # Read the observed values file
            #---------------------------------------------
            # Do this after first simulation, so we can
            # interpolate the observed values to match
            # the times of the simulated values.
            #---------------------------------------------
            if (index == 0):
                self.get_observed_values( SILENT=SILENT )
            # Y_obs and Y_sim were saved into self. 
            cost = self.compute_cost( SILENT=SILENT )
            costs[ index ] = cost
            print('   cost  =', cost )

            if (cost < min_cost):
                min_cost   = cost
                best_index = index
            index += 1

            if (PLOT):
                str1 = cal_var.replace('_',' ').title()
                title = (str1 + ' = ' + str(param))
                ## print('title =', title)
                ## self.plot_time_series( title=title )
                self.plot_time_series( title=title, normalize=True, marker=',')

        #--------------------------------------
        # Save value of p that minimizes cost
        #--------------------------------------
        self.params     = p_values
        self.costs      = costs
        self.best_param = p_values[ best_index ]

    #   calibrate()
    #---------------------------------------------------------------------
    def plot_time_series(self, title=None, sim_only=False,
                         obs_only=False, normalize=False, marker='+' ):
    
        x  = self.sim_times
        y1 = self.sim_values
        y2 = self.obs_values
        n1 = y1.size
        n2 = y2.size
        n  = min(n1, n2)

        if (normalize):
            y1 = ( y1 / y1.max() )
            y2 = ( y2 / y2.max() )

        if (sim_only):
            tfvis.plot_data(x[0:n], y1[0:n], marker=marker,
                            title=title,
                            x_name='Time', x_units='minutes',
                            y_name='Discharge', y_units='m3/s',
                            x_size=10, y_size=3)          
        elif (obs_only):
            tfvis.plot_data(x[0:n], y2[0:n], marker=marker,
                            title=title,
                            x_name='Time', x_units='minutes',
                            y_name='Discharge', y_units='m3/s',
                            x_size=10, y_size=3)         
        else:
            tfvis.plot_data(x[0:n], y1[0:n], y2[0:n], marker=marker,
                            title=title,
                            x_name='Time', x_units='minutes',
                            y_name='Discharge', y_units='m3/s',
                            x_size=10, y_size=3)   

    #   plot_time_series()
    #---------------------------------------------------------------------
    def update_sim_data_file(self):

        #--------------------------------------------------------    
        # Note: There is a check_overwrite() method in 
        #       topoflow/utils/file_utils.py that is used
        #       by rts_files.py, ncgs_files.py, etc.
        #       Numbers are appended just before the extension.
        #       This prevents accidental overwrite of existing
        #       output files that could be important.
        #----------------------------------------------------------
        # This method finds the most recent simulation file
        # with the right "signature".
        # self.sim_data_file  = self.sim_dir + (self.cfg_prefix + '_0D-Q.txt')
        #----------------------------------------------------------
        ## sim_files = ['Test_0D-Q.txt', 'Test_0D-Q_1.txt', 'Test_0D-Q_2.txt']    
        pattern   = (self.sim_dir + self.cfg_prefix + '_0D-Q*')
        sim_files = glob.glob( pattern )
        maxnum = -1
        for fname in sim_files:
            TEST1 = ('_0D-Q' in fname)
            TEST2 = fname.endswith('.txt')
            if (TEST1 and TEST2):     
                parts1 = fname.split('_0D-Q')
                ending = parts1[1]    # (e.g. '_25.txt' or '.txt')
                numstr = ending.replace('.txt', '')
                numstr = numstr.replace('_','')
                if (numstr != ''):
                    num = int( numstr )
                else:
                    num = 0
                if (num > maxnum):
                    maxnum = num
                    target = fname
         
        print('target file =', target )       
        self.sim_data_file = target

    #   update_sim_data_file()
    #---------------------------------------------------------------------
    def delete_output_files(self):

        ### os.remove( self.sim_data_file )

        #-----------------------------------      
        # Note:  This seems too dangerous.
        #-----------------------------------          
        os.chdir( self.output_dir )
        #------------------------------------
        txt_files = glob.glob( '*.txt' )
        for fname in txt_files: 
            print( fname )   
            ## os.remove( fname )
        #------------------------------------
        rts_files = glob.glob( '*.rts' )
        for fname in rts_files:
            print( fname )    
            ## os.remove( fname )        
        #------------------------------------
        nc_files = glob.glob( '*.nc' )
        for fname in nc_files:
            print( fname )    
            ## os.remove( fname )
             
    #   delete_output_files()
    #---------------------------------------------------------------------
                         
#---------------------------------------------------------------------
def read_any_time_series(txt_file, hdr_lines=0, time_col=0,
                         value_col=1, SILENT=False,
                         time_type='float', delim=None ):

    file_unit = open(txt_file, 'r')

    #--------------------    
    # Skip header lines
    #--------------------
    for k in range(hdr_lines):
        line = file_unit.readline()
    
    #------------------------    
    # Read times and values
    #-----------------------------------------------
    # Note:  If (delim is None), "split()" method
    #        splits on any whitespace.
    #-----------------------------------------------
    times  = list()
    values = list()
    while (True):
        line  = file_unit.readline()
        if (line == ''):
            break
        words = line.split( delim )
        ncols = len( words )
        if (time_col > ncols-1) or (value_col > ncols-1):
            print('ERROR: Number of columns =', ncols)
            print('   but time_col  =', time_col )
            print('   and value_col =', value_col )
            break
        if (time_type == 'float'):
            time  = np.float32( words[ time_col ] )
        else:
            time = words[ time_col ]  ## (string)
        value = np.float32( words[ value_col ] )
        times.append( time )
        values.append( value )
    #-------------------------------
    file_unit.close()
    times  = np.array( times )
    values = np.array( values )

    if (not(SILENT) and (time_type == 'float')):
        print('In read_time_series():')
        print('  Min(times)  =', times.min())
        print('  Max(times)  =', times.max())
        print('  Min(values) =', values.min())
        print('  Max(values) =', values.max())
        print()

    return (times, values)
             
#   read_any_time_series()
#---------------------------------------------------------------------
def compute_rain_rates_from_accumulations( SILENT=False ):

    #-------------------------------------------------
    # Note:  This function is not general as is, but
    #        could be adapted to be more general.
    #        It was written for a specific need.
    #-------------------------------------------------
    ## cfg_prefix   = 'June_20_67'
    examples_dir  = emeli.paths['examples']
    basin_dir     = examples_dir  + 'Treynor_Iowa_30m/'
    obs_dir       = basin_dir + '__observations/'
    obs_data_file = obs_dir + 'June_20_1967_Observed_Discharge.txt'

    txt_file     = obs_data_file
    hdr_lines    = 6
    time_col     = 1   # (decimal times, hours)
    value_col    = 3
    time_units   = 'hours'
    val_units    = 'inches'
    mm_per_inch  = 25.4
    (times, accum) = read_any_time_series(txt_file, hdr_lines=hdr_lines,
                          time_col=time_col, value_col=value_col,
                          SILENT=False, time_type='float' )
    n  = times.size
    dt = (np.roll( times, -1 ) - times)  # (roll forward)
    dt[n-1] = dt[n-2]
    da = (np.roll( accum, -1 ) - accum)
    da[n-1] = da[n-2] 
    da    = da * mm_per_inch      #  now [mm]
    rates = (da / dt)             #  now [mm / hour]
 
    ## return rates

    times_min = np.int16( 60.0 * times )
          
    # f1 = interp1d( times, values, kind='nearest')
    f1 = interp1d( times_min, rates, kind='linear')
    # f1 = interp1d( times_min, rates, kind='linear', bounds_error=False,
    #                fill_value="extrapolate" )
    tmin = times_min.min()
    tmax = times_min.max()
    trange = (tmax - tmin) + 1
    new_times = np.arange( trange) + tmin
    new_rates = f1( new_times )
    rmin = new_rates.min()
    rmax = new_rates.max()
    
    if not(SILENT):
        print('new_times[0:20] =', new_times[0:20] )
        print()
        print('new_rates[0:20] =', new_rates[0:20] ) 
        print()
        print('min(da), max(da) =', da.min(), ',', da.max(), '[mm]')
        print('min(dt), max(dt) =', dt.min(), ',', dt.max(), '[hours]')
        print('min(rates), max(rates) =',rmin, ',', rmax, '[mm h-1]' )
        print()
        
    return new_rates

#   compute_rain_rates_from_accumulations()  
#---------------------------------------------------------------------
#---------------------------------------------------------------------
def spearman_metric(x, y):

    #------------------------------------------------------------
    # Note:  This function computes a "distance metric" from
    #        the Spearman correlation coefficient, which itself
    #        indicates the degree to which y is a monotonic
    #        function of x.  Here, x and y are two 1D vectors.
    #        Any monotonic increasing function of x has d=1.
    #        Any monotonic decreasing function of x has d=0.
    #------------------------------------------------------------
    #        If we assume that river width, w, is a monotonic
    #        function of river discharge, Q, then d(w,Q)=0.
    #        Tests:  Y = c*X (c > 0), Y = exp(X)
    #------------------------------------------------------------
    # Note:  Instead of d=0, tests give:  d=5.551115123e-17.
    #------------------------------------------------------------     
    coeff, p = spearmanr( x, y )   # p-value not used here.
    d = (1.0 - coeff) / 2.0
    return d

#   spearman_metric                         
#---------------------------------------------------------------------


         
    
