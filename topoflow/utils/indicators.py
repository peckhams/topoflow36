#
#  Copyright (c) 2021, Scott D. Peckham
#
#  Oct. 2021. Continued improvements.
#  Sept 2021. Functions to create hydro indicators
#               from TopoFlow output files (netCDF).
#             Wrote create_stat_grid_stacks().
#             Wrote create_day_count_grid_stack().
#             Wrote create_affected_population_grid_stack().
#
#-------------------------------------------------------------------
#
#  test()
#  create_indicator_grid_stacks()
#  create_stat_grid_stacks()  # Do min, max, mean, sdev at once
#  get_condition_map()
#  create_day_count_grid_stack()
#  create_affected_population_grid_stack()
#
#-------------------------------------------------------------------

import os, os.path
import numpy as np
from topoflow.utils import regrid as rg  # read/regrid population data
from topoflow.utils import ncgs_files
from topoflow.utils import time_utils as tu

#-------------------------------------------------------------------
def test():

    create_indicator_grid_stacks()

#   test()
#-------------------------------------------------------------------
def create_indicator_grid_stacks(case_prefix=None, output_dir=None,
           compute_stat_grids=True):

    if (case_prefix is None):
        case_prefix = 'Test1'
    if (output_dir is None):
        output_dir = '/Users/peckhams/TF_Tests/Tana_120sec/'
    
    # var_file       = case_prefix + '_2D-Q.nc'
    # day_count_file = case_prefix + '_2D-days-per-month-lo-Q.nc'
    
    # hydro_var_names = ['discharge']
    # conditions_map  = {'discharge': ['low_discharge']}

    hydro_var_names = ['discharge', 'd_flood' ]    
    # hydro_var_names = ['discharge', 'd_flood', 'rainrate']
    conditions_list = {
    'discharge' : ['low_discharge', 'high_discharge',
                   'incr_discharge', 'decr_discharge'],
    'd_flood'   : ['flooding'],
    'rainrate'  : ['low_rainrate', 'high_rainrate'] }
    
    for hydro_name in hydro_var_names:
        print('-------------------------------------------')
        print('Working on hydro var name:', hydro_name )
        print('-------------------------------------------')
 
        #-------------------------------------------------       
        # This should not be done every time.  It should
        # be computed for a long "baseline" period and
        # then used to assess shorter duration runs.
        #-------------------------------------------------
        if (compute_state_grids):  
            create_stat_grid_stacks( output_dir=output_dir,
                case_prefix=case_prefix, hydro_var_name=hydro_name )

        for condition in conditions_list[ hydro_name ]:
            #-------------------------------
            # Uses a default var_file name
            #-------------------------------
            print('Hydro condition =', condition, '...' )
            print('-------------------------------------')
            create_day_count_grid_stack( output_dir=output_dir,
                case_prefix=case_prefix,
                hydro_var_name=hydro_name,
                condition=condition )
                # var_file=None )

            #-------------------------------------
            # Uses a default day_count_file name
            #-------------------------------------
            create_affected_population_grid_stack( output_dir=output_dir,
                case_prefix=case_prefix, condition=condition )
                # day_count_file=None )
          
#   create_indicator_grid_stacks()
#-------------------------------------------------------------------
def create_stat_grid_stacks( output_dir=None,
        case_prefix=None, hydro_var_name='discharge',
        var_file=None, min_file=None, max_file=None,
        mean_file=None, sdev_file=None):

    #---------------------------------------------------------------
    # Note:  For each year, there is a mean value for each month.
    #        There is also a "long-term" monthly mean for each
    #          month, which is the mean value for that month
    #          averaged over several (e.g. 10) years.
    #        We will need a separate function to compute these.
    #---------------------------------------------------------------
    # Note:  It may be even better to first consider the long-
    #        term (multi-year) mean and standard deviation for
    #        each day of the year.  Then count the number of
    #        days with unusually low or high values for that day.
    #        The long-term grid stacks would have 365 grids.
    #---------------------------------------------------------------    
 
    #-----------------------------------------------------
    # Note: This creates monthly grid stacks of several
    #       statistics (e.g. min, max, mean, sdev) from
    #       a grid stack of a hydrologic variable.
    #       Spatial resolution remains the same.
    #       All input & output files are netCDF format.
    #-------------------------------------------------------------
    # https://en.wikipedia.org/wiki/Variance#Population_variance
    #-------------------------------------------------------------
    if (output_dir is not None):
        os.chdir(output_dir)

    #-------------------------------------------------
    # Get input filename with the simulated variable
    #-------------------------------------------------
    if (var_file is None):
        if (hydro_var_name == 'discharge'):
            var_file = case_prefix + '_2D-Q.nc'  # (river discharge)
        elif (hydro_var_name == 'd_flood'):
            var_file = case_prefix + '_2D-d-flood.nc' # (flood depth)
        elif (hydro_var_name == 'rainrate'):
            var_file = case_prefix + '_2D-rainrate.nc' # (rainrate)
        else:
            print('Unknown hydro var name: ' + hydro_var_name)
            print('Choices are: discharge, d_flood, rainrate.')
            print()
            return

    s1_map = {'discharge':'Q', 'd_flood':'d-flood', 'rainrate':'rain'}
    s1 = s1_map[ hydro_var_name ]
    
    #-----------------------------------------------
    # Construct default filenames for output files
    #-----------------------------------------------
    if (min_file is None):
        min_file  = case_prefix + '_2D-month-min-'  + s1 + '.nc'
    if (max_file is None):
        max_file  = case_prefix + '_2D-month-max-'  + s1 + '.nc'
    if (mean_file is None):
        mean_file = case_prefix + '_2D-month-mean-' + s1 + '.nc'
    if (sdev_file is None):
        sdev_file = case_prefix + '_2D-month-sdev-' + s1 + '.nc'

    #-------------------------------------------------     
    # Prepare to read grid stack for var from netCDF
    #-------------------------------------------------
    ncgs_var = ncgs_files.ncgs_file()    # ncgs object
    OK = ncgs_var.open_file( var_file )
    if not(OK): return
                             
    #-------------------------------------------
    # Get some metadata from var_file (netCDF)
    #-------------------------------------------
    var_names = ncgs_var.get_var_names( no_dim_vars=True )
    var_name  = var_names[0]  # (e.g. 'Q')
    n_grids   = ncgs_var.ncgs_unit.variables[var_name].n_grids
    var_units = ncgs_var.ncgs_unit.variables[var_name].units  #######
    ncols     = ncgs_var.ncgs_unit.variables['X'].size
    nrows     = ncgs_var.ncgs_unit.variables['Y'].size    
    datetimes = ncgs_var.ncgs_unit.variables[ 'datetime' ]
    # times      = ncgs_var.ncgs_unit.variables[ 'time' ]
    # time_units = ncgs_var.ncgs_unit.variables[ 'time' ].units

    #---------------------------------------
    # Get grid_info & time_info structures
    #---------------------------------------
    grid_info = ncgs_var.get_grid_info(var_name, var_file)
    time_info = ncgs_var.get_time_info()
    #---------------------------------
    # Override some of the time_info
    #---------------------------------
    time_info.duration_units = 'months'  ######## ??
    time_info.duration = 1               ######## ??
    
    #-------------------------------
    # Get long statistic var names
    #-------------------------------
    min_long_name  = hydro_var_name + '_monthly_min'
    max_long_name  = hydro_var_name + '_monthly_max'     
    mean_long_name = hydro_var_name + '_monthly_mean'
    sdev_long_name = hydro_var_name + '_monthly_sdev'

    #--------------------------------
    # Get short statistic var names
    #--------------------------------
    min_short_name  = s1 + '_monthly_min'
    max_short_name  = s1 + '_monthly_max'
    mean_short_name = s1 + '_monthly_mean'
    sdev_short_name = s1 + '_monthly_sdev'
                                                                      
    #--------------------------------------------     
    # Prepare to write new grid stacks to netCDF
    # Add support for time_units='months' ??       ################
    #--------------------------------------------
    # The units for these 4 statistics are all
    # the same as the units of the variable.
    #--------------------------------------------    
    ncgs_min   = ncgs_files.ncgs_file()  # new ncgs object
    time_units = 'months'
    time_res   = 1
    dtype      = 'float32'
    #--------------------------------------------------------
    ncgs_min.open_new_file( min_file,
         grid_info=grid_info, time_info=time_info,
         var_name=min_short_name, long_name=min_long_name, 
         units_name=var_units, dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res )
    #--------------------------------------------------------
    ncgs_max = ncgs_files.ncgs_file()  # new ncgs object
    ncgs_max.open_new_file( max_file,
         grid_info=grid_info, time_info=time_info,
         var_name=max_short_name, long_name=max_long_name, 
         units_name=var_units, dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res )
    #--------------------------------------------------------
    ncgs_mean = ncgs_files.ncgs_file()  # new ncgs object
    ncgs_mean.open_new_file( mean_file,
         grid_info=grid_info, time_info=time_info,
         var_name=mean_short_name, long_name=mean_long_name, 
         units_name=var_units, dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res )
    #--------------------------------------------------------
    ncgs_sdev = ncgs_files.ncgs_file()  # new ncgs object
    ncgs_sdev.open_new_file( sdev_file,
         grid_info=grid_info, time_info=time_info,
         var_name=sdev_short_name, long_name=sdev_long_name, 
         units_name=var_units, dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res )
   
    max_num   = 9999999.0
    min_num   = -max_num
    shape     = (nrows, ncols)
    #---------------------------
    min_grid  = np.zeros( shape, dtype='float32')
    max_grid  = np.zeros( shape, dtype='float32')
    mean_grid = np.zeros( shape, dtype='float32')
    sdev_grid = np.zeros( shape, dtype='float32')
    #--------------------------------------------------
    mean_sqr_grid = np.zeros( shape, dtype='float32')

    #-----------------------------------------------
    # Read grids one-by-one and compute statistics
    #-----------------------------------------------
    print('Reading grids from netCDF file: ')
    print('   ' + var_file)
    print('& computing grid stacks for min, max, mean, sdev...')
    month_num = -1
    time = 0  # (number of months)
    for time_index in range(n_grids):
        #----------------------------
        # Parse the datetime string
        #----------------------------
        datetime_str = datetimes[ time_index ]
        (y,m1,d,h,m2,s) = tu.split_datetime_str(datetime_str, ALL=True)
        time += 1  # (number of months)
        
        #------------------------
        # Is this a new month ?
        #------------------------
        if (m1 != month_num):
            if (month_num != -1):
                #-------------------------------------
                # Variance = E(X^2) - E(X)^2
                # Std. deviation = sqrt( Variance )
                #-------------------------------------
                mean_grid     /= n_grids_per_month
                mean_sqr_grid /= n_grids_per_month
                sdev_grid[:] = np.sqrt( mean_sqr_grid - mean_grid**2)
                #----------------------------------------
                # Write statistic grids to netCDF files
                #----------------------------------------
                ncgs_min.add_grid(  min_grid,  min_short_name,
                                    time=time, time_units=time_units )
                ncgs_max.add_grid(  max_grid,  max_short_name,
                                    time=time, time_units=time_units )
                ncgs_mean.add_grid( mean_grid, mean_short_name,
                                    time=time, time_units=time_units )
                ncgs_sdev.add_grid( sdev_grid, sdev_short_name,
                                    time=time, time_units=time_units )

            #--------------------------------------
            # New month; re-initialize stat grids
            #-------------------------------------
            # All grid changes are done in-place
            #-------------------------------------            
            month_num = m1    # start new month
            n_grids_per_month = 0
            #----------------------
            np.maximum( min_grid, max_num, min_grid)  # set to max_num
            np.minimum( max_grid, min_num, max_grid)  # set to min_num
            mean_grid[:] = 0
            mean_sqr_grid[:] = 0
            sdev_grid[:] = 0

        #---------------------------------------
        # Get next grid from input netCDF file
        #---------------------------------------
        grid = ncgs_var.get_grid(var_name, time_index)
        n_grids_per_month += 1

        #-------------------------
        # Update statistic grids
        #-----------------------------------
        # Will divide mean and mean_sqr by
        # n_grids_per_month after summing.
        #-----------------------------------
        np.minimum( min_grid, grid, min_grid)  # (in-place)
        np.maximum( max_grid, grid, max_grid)  # (in-place)
        mean_grid     += grid       # will divide by number later
        mean_sqr_grid += grid**2.0  # will divide by number later
            
    #-----------------------------------
    # Close all input and output files
    #-----------------------------------
    ncgs_var.close_file()   # variable not variance
    ncgs_min.close_file()
    ncgs_max.close_file()
    ncgs_mean.close_file()
    ncgs_sdev.close_file()
    print('Finished.')
    print(' ')
    
#   create_stat_grid_stacks()
#-------------------------------------------------------------------
def get_condition_map():

    cmap = {
    'low_discharge':  'lo-Q',
    'high_discharge': 'hi-Q',
    'flooding':       'flood',
    'low_rainrate':   'lo-rain',
    'high_rainrate':  'hi-rain',
    'incr_discharge': 'incr-Q',
    'decr_discharge': 'decr-Q' }
    return cmap
    
#   get_condition_map()
#-------------------------------------------------------------------
def create_day_count_grid_stack( output_dir=None,
           case_prefix=None, hydro_var_name='discharge',
           condition='low_discharge',
           var_file=None, mean_file=None, sdev_file=None,
           days_file=None):  # days_file is output

    #----------------------------------------------------
    # Note: This creates a monthly grid stack for the
    #       number of days in each month and each grid
    #       cell for which some condition is satisfied.
    #       Conditions are defined in terms of a hydro
    #       variable and its long-term monthly mean and
    #       standard deviation.
    #       Spatial resolution remains the same.
    #----------------------------------------------------
    # Allowed conditions:
    #   low_discharge, high_discharge, flooding,
    #   low_rainrate, high_rainrate,
    #   decr_discharge, incr_discharge
    #---------------------------------------------------- 
    if (output_dir is not None):
        os.chdir(output_dir)

    #-------------------------------------------------
    # Get input filename with the simulated variable
    #-------------------------------------------------        
    if (var_file is None):
        if (hydro_var_name == 'discharge'):
            var_file = case_prefix + '_2D-Q.nc'  # (river discharge)
        elif (hydro_var_name == 'd_flood'):
            var_file = case_prefix + '_2D-d-flood.nc' # (flood depth)
        elif (hydro_var_name == 'rainrate'):
            var_file = case_prefix + '_2D-rainrate.nc' # (rainrate)
        else:
            print('Unknown hydro var name: ' + hydro_var_name)
            print('Choices are: discharge, d_flood, rainrate.')
            print()
            return

    s1_map = {'discharge':'Q', 'd_flood':'d-flood', 'rainrate':'rain'}
    s1 = s1_map[ hydro_var_name ]
    
    #------------------------------
    # Construct default var names
    #------------------------------
    cmap = get_condition_map()
    out_var_name  = 'days-per-month-' + cmap[ condition ]
    out_long_name = 'days_per_month_with_' + condition
             
    #------------------------------
    # Construct default filenames
    #------------------------------
    if (mean_file is None):
        mean_file = case_prefix + '_2D-month-mean-' + s1 + '.nc'
    if (sdev_file is None):
        sdev_file = case_prefix + '_2D-month-sdev-' + s1 + '.nc'
    if (days_file is None):
        days_file = case_prefix + '_2D-' + out_var_name + '.nc'

    #--------------------------------
    # Get short statistic var names
    #--------------------------------
    mean_short_name = s1 + '_monthly_mean'
    sdev_short_name = s1 + '_monthly_sdev'
    
    #------------------------------------------     
    # Prepare to read grid stacks from netCDF
    #------------------------------------------
    # The 1st one has full time-resolution;
    # other two have one-month time resolution
    #------------------------------------------
    ncgs_var = ncgs_files.ncgs_file()    # ncgs object
    OK = ncgs_var.open_file( var_file )
    if not(OK): return
    #------------------------------------------
    ncgs_mean = ncgs_files.ncgs_file()   # ncgs object
    OK = ncgs_mean.open_file( mean_file )
    if not(OK): return
    #------------------------------------------
    ncgs_sdev = ncgs_files.ncgs_file()   # ncgs object
    OK = ncgs_sdev.open_file( sdev_file )
    if not(OK): return

    #-----------------------------------------
    # Get some metadata from var netCDF file
    #-----------------------------------------
    var_names = ncgs_var.get_var_names( no_dim_vars=True )
    var_name  = var_names[0]  # (e.g. 'Q')
    n_grids   = ncgs_var.ncgs_unit.variables[var_name].n_grids
    var_units = ncgs_var.ncgs_unit.variables[var_name].units  #######
    ncols     = ncgs_var.ncgs_unit.variables['X'].size
    nrows     = ncgs_var.ncgs_unit.variables['Y'].size    
    datetimes = ncgs_var.ncgs_unit.variables[ 'datetime' ]
    # times      = ncgs_var.ncgs_unit.variables[ 'time' ]
    # time_units = ncgs_var.ncgs_unit.variables[ 'time' ].units
 
    #---------------------------------------
    # Get grid_info & time_info structures
    #---------------------------------------
    grid_info = ncgs_var.get_grid_info(var_name, var_file)
    time_info = ncgs_var.get_time_info()
    #---------------------------------
    # Override some of the time_info
    #---------------------------------
    time_info.duration_units = 'months'  ######## ??
    time_info.duration = 1               ######## ??
            
    #--------------------------------------------     
    # Prepare to write new grid stack to netCDF
    #--------------------------------------------  
    ncgs_days = ncgs_files.ncgs_file()  # new ncgs object
    ncgs_days.open_new_file( days_file,
         grid_info=grid_info, time_info=time_info,
         var_name=out_var_name, long_name=out_long_name, 
         units_name=var_units, dtype='float32',
         time_units='months', time_res='1', 
         MAKE_RTI=True)
    shape       = (nrows, ncols)
    detect_grid = np.zeros( shape, dtype='bool_')
    days_grid   = np.zeros( shape, dtype='int16')   #######
    last_grid   = np.zeros( shape, dtype='float32')

    #---------------------------------------------------    
    # Read grids one-by-one and compute day count grid
    #---------------------------------------------------
    print('Reading grids from netCDF file: ')
    print('   ' + var_file)
    print('and computing days-per-month grid stack')
    print('for condition: ' + condition + '...')
    #-------------------------------------------   
    day_num     = -1   # will be in {1,...,31}
    month_num   = -1   # will be in {1,...,12}
    day_count   = 0    # will count the days
    month_count = 0    # will count the months
    time_units  = 'months'
    #------------------------------------------------------------
    # NOTE!  time_index is an index in the ncgs_var grid stack
    #        which may correspond to a duration like 6 hours.
    #------------------------------------------------------------
    # NOTE!  Code assumes that the time between grids in  #######
    #        ncgs_var grid stack is <= one *day*.         #######
    #------------------------------------------------------------
    for time_index in range(n_grids):
        #----------------------------
        # Parse the datetime string
        #----------------------------
        datetime_str = datetimes[ time_index ]
        (y,m1,d,h,m2,s) = tu.split_datetime_str(datetime_str, ALL=True)

        #------------------------
        # Is this a new month ?
        #------------------------
        start_new_month = (m1 != month_num)
        if (start_new_month):
            finished_previous_month = (month_num != -1)
            if (finished_previous_month):
                #-------------------------------------------
                # Write days_grid for last month to netCDF
                #-------------------------------------------
                ncgs_days.add_grid( days_grid, out_var_name,
                                    time=month_count,
                                    time_units=time_units )
                month_count += 1    #############
 
            #-------------------------------------
            # New month; re-initialize days_grid
            #-------------------------------------
            month_num = m1    # start new month
            days_grid[:] = 0     
            n_grids_per_month = 0

            #-------------------------------------------
            # Get mean & sdev grids for this new month
            # These are long-term monthly stats.
            #-------------------------------------------
            try:
                mean_grid = ncgs_mean.get_grid( mean_short_name, month_count )
                sdev_grid = ncgs_sdev.get_grid( sdev_short_name, month_count )
                # long_mean_grid = ncgs_mean.get_grid( mean_short_name, month_count )
                # long_sdev_grid = ncgs_sdev.get_grid( sdev_short_name, month_count )
            except:
                #-------------------------------------------
                # No more grids in the monthly "mean_grid"
                #-------------------------------------------
                print('########################################')
                print(' WARNING: No more monthly mean grids.')
                print('          month_count =', month_count)
                print('########################################')
                print()
                break
      
            #---------------------------------------------
            # Compute threshold grid once per month here
            #---------------------------------------------
            if ('low' in condition):
                threshold = (mean_grid - sdev_grid)  # can be < 0
            elif ('high' in condition):
                threshold = (mean_grid + sdev_grid)
            else:
                pass
                ## threshold = 0
        
        #----------------------        
        # Is this a new day ?
        #--------------------------------------------------
        # NOTE!  Code assumes that the time between grids
        #        in ncgs_var grid stack is <= one day.
        #--------------------------------------------------
        start_new_day = (d != day_num)
        if (start_new_day):
            finished_previous_day = (day_num != -1)
            if (finished_previous_day):
                #-----------------------------------
                # Update the days_grid.
                # Was condition detected this day?
                #-----------------------------------
                w = (detect_grid == True)
                days_grid[w] += 1
                day_count += 1   #### not used yet
                
            #-------------------------------------
            # New day; re-initialize detect_grid
            #-------------------------------------
            day_num = d    # start new day
            detect_grid[:] = False
            n_grids_per_day = 0

        #---------------------------------------    
        # Get next grid from input netCDF file
        #---------------------------------------
        grid = ncgs_var.get_grid(var_name, time_index)
        n_grids_per_month += 1
        n_grids_per_day   += 1

        #--------------------------------------------------
        # Update detect_grid. Did condition occur in day?
        # w is where condition is true (boolean array)
        #--------------------------------------------------
        if ('low' in condition):
            w = (grid < threshold)
        elif ('high' in condition):
            w = (grid < threshold)
        elif (condition == 'flooding'):
            w = (grid > 0)
            # w = (grid > 0.1)  # [meters]
        elif ('incr_' in condition):
            w = (grid > last_grid)
            last_grid[:] = grid   # in-place
        elif ('decr_' in condition):
            w = (grid < last_grid)
            last_grid[:] = grid   # in-place
        else:
            print('ERROR: Unknown condition.')
            print()
        #----------------------------------
        detect_grid[w] = True 
    
    #-------------------------------
    # Close input and output files
    #-------------------------------
    ncgs_var.close_file()
    ncgs_mean.close_file()
    ncgs_sdev.close_file()
    ncgs_days.close_file()     
    print('Finished.')
    print(' ')
    
#   create_day_count_grid_stack()
#-------------------------------------------------------------------
def create_affected_population_grid_stack( output_dir=None,
           case_prefix=None,
           condition='low_discharge', day_count_file=None,
           pop_file_in=None, pop_file_out=None,
           out_file=None, out_var_name=None):

    #-------------------------------------------------------
    # Note: This creates an indicator grid based on the
    #       number of people that are subjected to a given
    #       condition (e.g. hydrologic) in a given month.
    #-------------------------------------------------------
    # Note: day_count_file has number of days each month,
    #       in each grid cell, that a condition was True.
    #       day_count_file is computed by the function:
    #           create_day_count_grid_stack()
    #       The condition is indicated in the filename.    
    #-------------------------------------------------------
    # Note: Given a day_count_file for some condition, X,
    #       and a population count grid (e.g. GPW-v4)
    #       resampled to the same grid cell size and
    #       clipped to same bounding box, we can create
    #       an "affected population grid" via:
    #          w = (days_per_month_with_X_grid > 0)
    #          out_grid = pop_grid[ w ] 
    #-------------------------------------------------------
    if (output_dir is not None):
        os.chdir(output_dir)

    cmap = get_condition_map()
                
    #------------------------   
    # Set default filenames
    #------------------------
    if (pop_file_in is None):
        pop_file_dir = '/Users/peckhams/Data/GPW-v4/'
        pop_file     = 'gpw_v4_population_count_rev11_2020_30_sec.tif'
        pop_file_in  = pop_file_dir + pop_file
    #-------------------------------------------------------       
    if (day_count_file is None):
        prefix = case_prefix + '_2D-days-per-month-'
        suffix = cmap[ condition ] + '.nc'
        day_count_file = prefix + suffix
    #-------------------------------------------------------
    if (out_file is None):
        prefix    = case_prefix + '_2D-pop-1-or-more-dpm-'
        suffix   = cmap[ condition ] + '.nc'
        out_file = prefix + suffix
    #-------------------------------------------------------
    if (out_var_name is None):
        prefix = 'pop-1-or-more-dpm-'
        suffix = cmap[ condition ]
        out_var_name = prefix + suffix
    #-------------------------------------------------------
    # prefix = 'people~near-river'  # for discharge and d_flood
    prefix = 'people~in-grid-cell'
    midstr = '-with-1-more-days-per-month-of-'
    suffix = condition + '__count'
    out_long_name = prefix + midstr + suffix

    #----------------------------------------------------------     
    # Prepare to read "days-per-month" grid stack from netCDF
    #----------------------------------------------------------
    ncgs_days = ncgs_files.ncgs_file()
    OK = ncgs_days.open_file( day_count_file )  # days-per-month
    if not(OK): return

    #--------------------------------------------
    # Get some metadata from the day_count_file
    #--------------------------------------------
    var_names = ncgs_days.get_var_names( no_dim_vars=True )
    var_name  = var_names[0]  # (e.g. 'Q')
    n_grids   = ncgs_days.ncgs_unit.variables[var_name].n_grids
    var_units = ncgs_days.ncgs_unit.variables[var_name].units  #######
    ncols     = ncgs_days.ncgs_unit.variables['X'].size
    nrows     = ncgs_days.ncgs_unit.variables['Y'].size    
    datetimes = ncgs_days.ncgs_unit.variables[ 'datetime' ]
    # times      = ncgs_days.ncgs_unit.variables[ 'time' ]
    # time_units = ncgs_days.ncgs_unit.variables[ 'time' ].units

    #---------------------------------------
    # Get grid_info & time_info structures
    #---------------------------------------
    grid_info = ncgs_days.get_grid_info(var_name, day_count_file)
    time_info = ncgs_days.get_time_info()
    # time_info.duration_units = 'months'  # not needed here ?
    # time_info.duration = 1               # not needed here ?

    #--------------------------------------------     
    # Prepare to write new grid stack to netCDF
    #--------------------------------------------  
    ncgs_out   = ncgs_files.ncgs_file()  # new ncgs object
    ## dtype      = 'int32'
    dtype      = 'float32'  # (actual dtype in GPW-v4)
    time_units = 'months'
    time_res   = 1
    #---------------------------------
    ncgs_out.open_new_file( out_file,
         grid_info=grid_info, time_info=time_info,
         var_name=out_var_name, long_name=out_long_name, 
         units_name='1', dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res ) 
                   
    #-------------------------------------------
    # Get out_xres_sec, out_yres_sec, & bounds
    #-------------------------------------------
    out_xres_sec = grid_info.xres  # [arcseconds]
    out_yres_sec = grid_info.yres  # [arcseconds]
    out_bounds = ncgs_days.ncgs_unit.geospatial_bounds     #### CHECK ORDERING
    ######################################################
    out_grid = np.zeros( (nrows, ncols), dtype=dtype)

    #---------------------------------------------------
    # Use "out_xres_sec" for default pop_file_out name
    #---------------------------------------------------    
    if (pop_file_out is None):
        prefix   = 'Test1_2D_gpw_v4_pop_count_rev11_2020_'
        suffix   = '_sec.tif'
        suffix2  = '_sec.nc'
        xres_str = str(int(out_xres_sec))
        pop_file_out    = prefix + xres_str + suffix  ######
        pop_file_out_nc = prefix + xres_str + suffix2
       
    #-------------------------------------------   
    # Read population grid (GPW-v4 as geotiff).
    # Clip and resample as necessary.
    #-------------------------------------------
    if not(os.path.exists( pop_file_out )):
        print('Clipping & resampling GPW-v4 pop. count grid...')
        rg.regrid_geotiff(in_file=pop_file_in,
            out_file=pop_file_out, out_bounds=out_bounds,
            out_xres_sec=out_xres_sec, out_yres_sec=out_yres_sec,
            RESAMPLE_ALGO='bilinear', REPORT=True)

    #---------------------------------------------------
    # GPW-v4 data type is "float32".
    # GPW-v4 nodata value = -3.402823e+38 (Water Mask)
    #---------------------------------------------------
    pop_grid = rg.read_geotiff( pop_file_out )
    w1 = (pop_grid < 0)   # boolean array
    pop_grid[ w1 ] = 0.0

    #--------------------------------------     
    # Write new population grid to netCDF
    #--------------------------------------
    if not(os.path.exists( pop_file_out_nc )):
        print('Writing new population grid to netCDF...')
        ncgs_pop   = ncgs_files.ncgs_file()  # new ncgs object
        dtype      = 'float32'  # (actual dtype in GPW-v4)
        # dtype    = 'int32'
        #------------------------------------------
        ncgs_pop.open_new_file( pop_file_out_nc,
             grid_info=grid_info, time_info=time_info,
             var_name='pop_count', long_name='population__count', 
             units_name='1', dtype=dtype, MAKE_RTI=True,
             time_units='years', time_res=1 )  ####################
        ncgs_pop.add_grid( pop_grid, 'pop_count',
                          time=0, time_units='years' )
        ncgs_pop.close()
        print('Finished writing population grid.')
        print()
                                                     
    #------------------------------------------------
    # Read dpm grids one-by-one & compute indicator
    #------------------------------------------------
    print('Reading grids from netCDF file: ')
    print('   ' + day_count_file)
    print('and computing new indicator...')

    #-------------------------------------------------
    # Note: time_index = number of months from start
    #-------------------------------------------------
    for time_index in range(n_grids):   
        #-------------------------------------
        # Read next grid from day_count_file
        #-------------------------------------
        days_per_month_grid = ncgs_days.get_grid(var_name, time_index)
        
        #-------------------------------------------            
        # Compute indicator grid & write to netCDF
        #-------------------------------------------    
        w1 = (days_per_month_grid > 0)  # i.e. "one or more"
        w2 = np.invert( w1 )
        out_grid[w1] = pop_grid[ w1 ]
        out_grid[w2] = 0 
        ncgs_out.add_grid( out_grid, out_var_name,
                           time=time_index, time_units=time_units )
               
    #--------------
    # Close files
    #--------------
    # pop_file was closed already.
    ncgs_days.close_file()
    ncgs_out.close_file()   
    print('Finished computing indicator grid stack.')
    print(' ')
    
#   create_affected_population_grid_stack()
#-------------------------------------------------------------------


