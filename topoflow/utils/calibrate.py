
#------------------------------------------------------------------------
#  Copyright (c) 2020, Scott D. Peckham
#
#  May 2020.  Created for TopoFlow calibration notebook.
#
#---------------------------------------------------------------------
#
#  test()
#
#  class calibrator():
#      __init__()
#      set_file_info()
#      print_file_info()
#
#      read_time_series()
#      convert_times_from_hhmm_to_minutes()
#      convert_times_from_minutes_to_hhmm()
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
#      plot_hydrographs()
#      update_sim_data_file()
#      delete_output_files()    ## commented out; dangerous
#
#--------------------------------------------
#  These functions are outside of the class
#--------------------------------------------
#  read_any_time_series()
#  compute_rain_rates_from_accumulations()
#
#---------------------------------------------------------------------

import glob, os, os.path, shutil, time

import numpy as np
from scipy.interpolate import interp1d

from topoflow.framework import emeli  # (for examples dir)
from topoflow.utils import parameterize
from topoflow.utils import visualize as tfvis
from topoflow import main

#------------------------------------------------------------------------
def test():

    c = calibrator()
    c.print_file_info()
    c.get_observed_values( SILENT=False  )
    c.get_simulated_values( SILENT=False )
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

    def __init__(self):

        self.version       = '0.5'
        self.p_for_Lp_norm = 2.0       # L2 norm is default
        self.set_file_info()

    #   __init__()
    #---------------------------------------------------------------------
    def set_file_info(self):
 
        #---------------------------------
        # Set cfg_prefix and directories
        #---------------------------------
        self.cfg_prefix   = 'June_20_67'
        self.site_prefix  = 'Treynor'
        self.home_dir     = os.path.expanduser("~")
        self.examples_dir = emeli.paths['examples']
        self.basin_dir    = self.examples_dir  + 'Treynor_Iowa_30m/'
        self.topo_dir     = self.basin_dir + '__topo/'
        self.cfg_dir      = self.basin_dir + '__No_Infil_June_20_67_rain/'
        self.output_dir   = self.home_dir + '/TF_Output/Treynor/'

        #---------------------------------------
        # Attributes of the observed data file
        #---------------------------------------
        self.HHMM_TIMES       = True
        self.IRREGULAR_TIMES  = True
        self.interp_method    = 'linear'  # or 'nearest', etc.
        self.obs_header_lines = 6
        self.obs_time_column  = 0    # (HHMM)
        self.obs_value_column = 5
        #---------------------------------------------------
        obs_dir = self.basin_dir + '__observations/'
        self.obs_dir       = obs_dir
        self.obs_data_file = obs_dir + 'June_20_1967_Observed_Discharge.txt'

        #----------------------------------------        
        # Attributes of the simulated data file
        #-------------------------------------------------
        # obs is multi-column, sim is single-column now.
        #-------------------------------------------------
        self.sim_header_lines = 2
        self.sim_time_column  = 0    # [min]
        self.sim_value_column = 1
        #---------------------------------------------------
        # self.sim_data_file will be replaced each time
        self.sim_dir        = self.output_dir
        self.sim_data_file0 = self.sim_dir + (self.cfg_prefix + '_0D-Q.txt')
        self.sim_data_file  = self.sim_dir + (self.cfg_prefix + '_0D-Q.txt')
              
    #   set_file_info()
    #---------------------------------------------------------------------
    def print_file_info(self):
 
        #-----------------------------------
        # Print cfg_prefix and directories
        #-----------------------------------
        print('In print_file_info():')
        print('  cfg_prefix   =', self.cfg_prefix )
        print('  site_prefix  =', self.site_prefix )
        print('  examples_dir =', self.examples_dir )
        print('  basin_dir    =', self.basin_dir )
        print('  topo_dir     =', self.topo_dir )
        print('  output_dir   =', self.output_dir )
        print('  cfg_dir      =', self.cfg_dir )
        print('  home_dir     =', self.home_dir )
        print('  obs_dir      =', self.obs_dir )
        print('  sim_dir      =', self.sim_dir )
        print()

    #   print_file_info()     
    #---------------------------------------------------------------------
    def read_time_series(self, source='observed', SILENT=True,
                         time_type='float' ):

        if (source == 'observed'):
            txt_file  = self.obs_data_file
            hdr_lines = self.obs_header_lines
            time_col  = self.obs_time_column
            value_col = self.obs_value_column
        else:     ## 'simulated'
            txt_file  = self.sim_data_file
            hdr_lines = self.sim_header_lines
            time_col  = self.sim_time_column
            value_col = self.sim_value_column
            time_type = 'float'

        file_unit = open(txt_file, 'r')

        #--------------------    
        # Skip header lines
        #--------------------
        for k in range(hdr_lines):
            line = file_unit.readline()
        
        #------------------------    
        # Read times and values
        #-----------------------------------------------
        # Note:  "split()" method splits on whitespace
        #-----------------------------------------------
        times  = list()
        values = list()
        while (True):
            line  = file_unit.readline()
            if (line == ''):
                break
            words = line.split()
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
    def convert_times_from_hhmm_to_minutes(self):
 
        if not(self.HHMM_TIMES):
            print('WARNING: self.HHMM_TIMES = False.')

        times_hhmm = self.obs_times
        
        n_times   = len( times_hhmm )
        times_min = np.zeros( n_times )
        for k in range(n_times):
            hhmm = times_hhmm[k]
            hour = np.int16( hhmm[:2] )
            min  = np.int16( hhmm[2:] )
            times_min[k] = (hour * 60) + min
 
        #-------------------------------------------   
        # Replace HHMM times with time in minutes.
        #-------------------------------------------
        self.obs_times = times_min
    
    #   convert_hhmm_to_minutes()
    #---------------------------------------------------------------------
    def convert_times_from_minutes_to_hhmm(self):

        if not(self.HHMM_TIMES):
            print('WARNING: self.HHMM_TIMES = False.')
            
        times_min = self.obs_times
        
        n_times   = len( times_min )
        times_hhmm = np.zeros( n_times )
        for k in range(n_times):
            hour = int( times_min[k] / 60 )
            min  = int( times_min[k] % 60 )
            #-----------------------------------        
            hh   = str(hour)
            if (len(hh) == 1):  hh = ('0' + hh)
            #-----------------------------------
            mm   = str(min)
            if (len(mm) == 1):  mm = ('0' + mm)
            #-----------------------------------
            times_hhmm[k] = int(hh + mm)
 
        #-------------------------------------------   
        # Replace time in minutes with HHMM times.
        #-------------------------------------------   
        self.obs_times = times_hhmm
    
    #   convert_minutes_to_hhmm()
    #---------------------------------------------------------------------
    def regularize_obs_time_series( self, SILENT=True ): 

        #----------------------------------------------------------
        # Notes:  This function is for the special case where the
        #         time column contains times in HHMM format that
        #         are not evenly spaced in time.  It returns a
        #         new time series (decimal hours) with a regular
        #         spacing of 1 minute.
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
              
        # f1 = interp1d( times, values, kind='nearest')
        f1 = interp1d( times, values, kind=self.interp_method)
        tmin = times.min()
        tmax = times.max()
        trange = (tmax - tmin) + 1
        new_times  = np.arange( trange) + tmin
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
        if (self.HHMM_TIMES):
            self.read_time_series( source='observed',
                                   time_type='str' )
            self.convert_times_from_hhmm_to_minutes()
        else:
            self.read_time_series( source='observed' )
                         
        #------------------------------------------------------  
        # Use 1D interpolation to create even spacing in time
        #------------------------------------------------------
        if (self.IRREGULAR_TIMES):
            self.regularize_obs_time_series()
            #--------------------------------------------------
            # Returning to HHMM-format times isn't necessary.
            # It depends on what we want to do next with the
            # times, like using them for x-axis of a plot.
            # They are not used to compute the cost function.
            #--------------------------------------------------
            if (self.HHMM_TIMES):
                self.convert_times_from_minutes_to_hhmm()
    
        if not(SILENT):
            print('HHMM_TIMES      =', self.HHMM_TIMES )
            print('IRREGULAR_TIMES =', self.IRREGULAR_TIMES )
            print('Min(obs_times)  =', self.obs_times.min() )
            print('Max(obs_times)  =', self.obs_times.max() )
            print('Min(obs_values) =', self.obs_values.min() )
            print('Max(obs_values) =', self.obs_values.max() )
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
            print('Min(sim_times)  =', self.sim_times.min() )
            print('Max(sim_times)  =', self.sim_times.max() )
            print('Min(sim_values) =', self.sim_values.min() )
            print('Max(sim_values) =', self.sim_values.max() )
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
        d8_area_file = site_prefix + '_area.rtg'
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

        #----------------------------------------------------
        # Note: Lp norm cost function with array operations
        #----------------------------------------------------
        p     = self.p_for_Lp_norm
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
            ### return None

        arg = np.abs( Y_obs - Y_sim )**p
        cost = (arg.sum())**(1.0/p)
    
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

        #-----------------------------------
        # Read the observed discharge data
        #-----------------------------------
        self.get_observed_values( SILENT=SILENT )

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
            range = [0.01, 0.09]
        if (cal_var == 'manning_n_max'):
            range = [0.1, 0.9]
        p_values = self.get_parameter_values( range=range, n=15)
        costs = np.zeros( 15 )

        #-------------------------------------------      
        # Run the model repeatedly, vary parameter
        #-------------------------------------------
        print('Working...')
        for param in p_values:
            print('   param =', param )
            ## self.modify_cfg_file( p )
            self.modify_channel_grid_file( cal_var, param )
            time.sleep( 1.0 ) 
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
            # Y_obs and Y_sim were saved into self. 
            cost = self.compute_cost( SILENT=SILENT )
            costs[ index ] = cost
                  
            if (cost < min_cost):
                min_cost   = cost
                best_index = index
            index += 1

            if (PLOT):
                str1 = cal_var.replace('_',' ').title()
                title = (str1 + ' = ' + str(param))
                ## print('title =', title)
                self.plot_hydrographs( title=title )

        #--------------------------------------
        # Save value of p that minimizes cost
        #--------------------------------------
        self.params     = p_values
        self.costs      = costs
        self.best_param = p_values[ best_index ]

    #   calibrate()
    #---------------------------------------------------------------------
    def plot_hydrographs(self, title=None ):
    
        x  = self.sim_times
        y1 = self.sim_values
        y2 = self.obs_values
        n1 = y1.size
        n2 = y2.size
        n  = min(n1, n2)

        tfvis.plot_data(x[0:n], y1[0:n], y2[0:n], marker='+',
                        title=title,
                        x_name='Time', x_units='minutes',
                        y_name='Discharge', y_units='m3/s',
                        x_size=10, y_size=3)   

    #   plot_hydrographs()
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
                         time_type='float' ):

    file_unit = open(txt_file, 'r')

    #--------------------    
    # Skip header lines
    #--------------------
    for k in range(hdr_lines):
        line = file_unit.readline()
    
    #------------------------    
    # Read times and values
    #-----------------------------------------------
    # Note:  "split()" method splits on whitespace
    #-----------------------------------------------
    times  = list()
    values = list()
    while (True):
        line  = file_unit.readline()
        if (line == ''):
            break
        words = line.split()
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

             
                          
             
             
    
