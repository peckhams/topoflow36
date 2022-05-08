#
#  Copyright (c) 2021-2022, Scott D. Peckham
#
#  Feb. 2022. Minor bug fixes.
#             Improved create_stat_grid_stacks().
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
#  create_population_grid_file    # separate fcn on 2021-10-26
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

    pop_dir = '~/basins/kenya/tana_120sec/__misc/'
    pop_dir = os.path.expanduser( pop_dir )
    
    create_indicator_grid_stacks(
           case_prefix='Test1',
           output_dir=os.path.expanduser('~/output/'),
           pop_dir=pop_dir,
           compute_stat_grids=True,
           OVERWRITE_OK=False)

#   test()
#-------------------------------------------------------------------
def create_indicator_grid_stacks(case_prefix=None, output_dir=None,
           compute_stat_grids=True, stats_dir=None,
           pop_dir=None, pop_file=None, FLOODING=True,
           OVERWRITE_OK=False):

    if (case_prefix is None):
        case_prefix = 'Test1'
    #-------------------------------------------------------
    if (output_dir is None):
        output_dir = os.path.expanduser('~/output/')
    elif (output_dir[-1] != os.sep):
        output_dir += os.sep
    #-------------------------------------------------------
    if (stats_dir is None):
        stats_dir = os.path.expanduser('~/stats1/')
    elif (stats_dir[-1] != os.sep):
        stats_dir += os.sep
    #-------------------------------------------------------
    if (pop_dir is None):
        print('#### ERROR: pop_dir argument is required.')  
        ### pop_dir = '/Users/peckhams/Data/GPW-v4/'
    elif (pop_dir[-1] != os.sep):
        pop_dir += os.sep
    #-------------------------------------------------------
    if (pop_file is None):
        pass  # Let it be auto-constructed 
        
    # var_file       = case_prefix + '_2D-Q.nc'
    # day_count_file = case_prefix + '_2D-days-per-month-lo-Q.nc'
    
    # hydro_var_names = ['discharge']
    # conditions_map  = {'discharge': ['low_discharge']}

    if (FLOODING):
        hydro_var_names = ['discharge', 'd_flood' ]
    else:
        hydro_var_names = ['discharge']
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
        if (compute_stat_grids):
            create_stat_grid_stacks( output_dir=output_dir,
                stats_dir=stats_dir,
                case_prefix=case_prefix,
                hydro_var_name=hydro_name,
                OVERWRITE_OK=OVERWRITE_OK )

        for condition in conditions_list[ hydro_name ]:
            #-------------------------------
            # Uses a default var_file name
            #-------------------------------
            print('Hydro condition =', condition, '...' )
            print('-------------------------------------')
            create_day_count_grid_stack( output_dir=output_dir,
                stats_dir=stats_dir,
                case_prefix=case_prefix,
                hydro_var_name=hydro_name,
                condition=condition,
                OVERWRITE_OK=OVERWRITE_OK )
                # var_file=None )

            #-------------------------------------
            # Uses a default day_count_file name
            #-------------------------------------
            create_affected_population_grid_stack( output_dir=output_dir,
                case_prefix=case_prefix, condition=condition,
                pop_dir=pop_dir, pop_file=pop_file,
                OVERWRITE_OK=OVERWRITE_OK )
                # day_count_file=None )
            #--------------------------------------------------------------
#             create_affected_population_grid_stack( output_dir=None,
#                 pop_dir=None, pop_file=None, case_prefix=None,
#                 # Next 2 not needed if pop_file exists.
#                 global_pop_dir=None, global_pop_file=None, 
#                 condition=condition, day_count_file=None,
#                 out_file=None, out_var_name=None
                     
#   create_indicator_grid_stacks()
#-------------------------------------------------------------------
def create_stat_grid_stacks( output_dir=None, stats_dir=None,
        case_prefix=None, hydro_var_name='discharge',
        var_file=None, min_file=None, max_file=None,
        mean_file=None, sdev_file=None,
        OVERWRITE_OK=False):

    #-------------------------------------------------------------
    # Note: This creates monthly grid stacks of several stats
    #       (e.g. min, max, mean, sdev) from a grid stack of a
    #       hydrologic variable.
    #       Spatial resolution remains the same.
    #       All input & output files are netCDF format.
    #-------------------------------------------------------------
    # https://en.wikipedia.org/wiki/Variance#Population_variance
    #-----------------------------------------------------------------
    # Option 1:  First, for each of 10 years, compute the 4 stats
    #            for each month, for each grid cell.  The resulting
    #            grid stack for each stat will then have 120 grids.
    #            Next, for each month, compute the 4 stats as an
    #            average over 10 years.  The resulting grid stack
    #            for each stat will then have 12 grids. 
    #-----------------------------------------------------------------
    # Option 2:  For each grid cell, compute the long-term (10-yr)
    #            stats (mean, sdev, min and max) using all values
    #            for that grid cell (regardless of month). There
    #            will then be 4 grids, 1 for each of the 4 statistics.
    #            Since discharge is naturally much lower (or higher)
    #            in certain months, this doesn't seem that helpful.
    #-----------------------------------------------------------------
    # Option 3:  For each grid cell, compute the long-term (10-yr)
    #            stats (mean, sdev, min, & max) for each day of the
    #            year.  Each stat will then be computed from just 10
    #            values, and the grid stack for each stat will have
    #            365 grids.
    #-----------------------------------------------------------------
    home_dir = os.path.expanduser('~') + os.sep    
    if (output_dir is None):
        output_dir = home_dir + 'output/'
    if (stats_dir is None):
        stats_dir = home_dir + 'stats1/'
    #---------------------------------------
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
    # Construct default filenames for "stat" files
    #-----------------------------------------------
    prefix = (stats_dir + case_prefix)  ######
    if (min_file is None):
        min_file  = prefix + '_2D-month-min-'  + s1 + '.nc'
    if (max_file is None):
        max_file  = prefix + '_2D-month-max-'  + s1 + '.nc'
    if (mean_file is None):
        mean_file = prefix + '_2D-month-mean-' + s1 + '.nc'
    if (sdev_file is None):
        sdev_file = prefix + '_2D-month-sdev-' + s1 + '.nc'

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
    min_long_name  = hydro_var_name + '_monthly_min'   ###### CHECK NAMES
    max_long_name  = hydro_var_name + '_monthly_max'     
    mean_long_name = hydro_var_name + '_monthly_mean'
    sdev_long_name = hydro_var_name + '_monthly_sdev'

    #--------------------------------
    # Get short statistic var names
    #--------------------------------
    min_short_name  = s1 + '_monthly_min'    ####### CHECK NAMES
    max_short_name  = s1 + '_monthly_max'
    mean_short_name = s1 + '_monthly_mean'
    sdev_short_name = s1 + '_monthly_sdev'
                                                                      
    #--------------------------------------------     
    # Prepare to write new grid stacks to netCDF
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
         time_units=time_units, time_res=time_res,
         OVERWRITE_OK=OVERWRITE_OK )
    #--------------------------------------------------------
    ncgs_max = ncgs_files.ncgs_file()  # new ncgs object
    ncgs_max.open_new_file( max_file,
         grid_info=grid_info, time_info=time_info,
         var_name=max_short_name, long_name=max_long_name, 
         units_name=var_units, dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res,
         OVERWRITE_OK=OVERWRITE_OK )
    #--------------------------------------------------------
    ncgs_mean = ncgs_files.ncgs_file()  # new ncgs object
    ncgs_mean.open_new_file( mean_file,
         grid_info=grid_info, time_info=time_info,
         var_name=mean_short_name, long_name=mean_long_name, 
         units_name=var_units, dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res,
         OVERWRITE_OK=OVERWRITE_OK )
    #--------------------------------------------------------
    ncgs_sdev = ncgs_files.ncgs_file()  # new ncgs object
    ncgs_sdev.open_new_file( sdev_file,
         grid_info=grid_info, time_info=time_info,
         var_name=sdev_short_name, long_name=sdev_long_name, 
         units_name=var_units, dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res,
         OVERWRITE_OK=OVERWRITE_OK )
   
    max_num = 9999999.0
    min_num = -max_num
    n_years = np.zeros(12, dtype='float32')   # Num times each month occurs 
    #-------------------------------------------------
    shape0     = (12, nrows, ncols)
    min_grids  = np.zeros( shape0, dtype='float32') + max_num
    max_grids  = np.zeros( shape0, dtype='float32') + min_num
    mean_grids = np.zeros( shape0, dtype='float32')
    sdev_grids = np.zeros( shape0, dtype='float32')
    #-------------------------------------------------
    shape     = (nrows, ncols)
    min_grid  = np.zeros( shape, dtype='float32') + max_num
    max_grid  = np.zeros( shape, dtype='float32') + min_num
    mean_grid = np.zeros( shape, dtype='float32')
    sdev_grid = np.zeros( shape, dtype='float32')
    #-------------------------------------------------
    mean_sqr_grid = np.zeros( shape, dtype='float32')

    #-----------------------------------------------
    # Read grids one-by-one and compute statistics
    #-----------------------------------------------
    print('Reading grids from netCDF file: ')
    print('   ' + var_file)
    print('& computing grid stacks for min, max, mean, sdev...')
    month_num  = -1
    ## month_count = 0

    for time_index in range(n_grids):
        #----------------------------
        # Parse the datetime string
        #----------------------------
        datetime_str = datetimes[ time_index ]
        (y,m1,d,h,m2,s) = tu.split_datetime_str(datetime_str, ALL=True)
        
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

                #-----------------------------------
                # Update monthly long-term average
                #-----------------------------------
                m = (month_num - 1)
                n_years[ m ] += 1    ##########
                #--------------------------------------
                cur_min = min_grids[m, :,:]
                min_grids[m, :,:] = np.minimum( min_grid, cur_min )
                cur_max = max_grids[m, :,:]
                max_grids[m, :,:] = np.maximum( max_grid, cur_max )
                #------------------------------------------------------              
                mean_grids[m, :,:] += mean_grid  # (divide at end)
                sdev_grids[m, :,:] += sdev_grid  # (divide at end)
                
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
        #----------------------------------------------
        # For Q, u, etc. first grid has all zeros.      ###############
        # Probably shouldn't write it to file at all.
        # Exception:  d starts with d0.
        #----------------------------------------------
        if (time_index > 0):
            np.minimum( min_grid, grid, min_grid)  # (in-place)
        np.maximum( max_grid, grid, max_grid)  # (in-place)
        mean_grid     += grid       # will divide by number later
        mean_sqr_grid += grid**2.0  # will divide by number later

    #--------------------------------
    # Save info for the final month
    # (Bug fix: 2022-02-16)
    #--------------------------------
    mean_grid     /= n_grids_per_month
    mean_sqr_grid /= n_grids_per_month
    sdev_grid[:] = np.sqrt( mean_sqr_grid - mean_grid**2)
    #-----------------------------------
    # Update monthly long-term average
    #-----------------------------------
    m = (month_num - 1)
    cur_min = min_grids[m, :,:]
    min_grids[m, :,:] = np.minimum( min_grid, cur_min )
    cur_max = max_grids[m, :,:]
    max_grids[m, :,:] = np.maximum( max_grid, cur_max )
    #------------------------------------------------------              
    mean_grids[m, :,:] += mean_grid  # (divide at end)
    sdev_grids[m, :,:] += sdev_grid  # (divide at end)
      
    #----------------------------------------
    # Write statistic grids to netCDF files
    # There will be 12 grids for each stat.
    #----------------------------------------
    for m in range(12):    # {0,...,11}
        #------------------------------------------------
        # Divide by number of times each month occurred
        #------------------------------------------------
        mean_grids[m, :,:] /= n_years[ m ]
        sdev_grids[m, :,:] /= n_years[ m ]
        #------------------------------------------------------------    
        ncgs_min.add_grid(  min_grids[m, :,:],  min_short_name,
                            time=m, time_units=time_units )
        ncgs_max.add_grid(  max_grids[m, :,:],  max_short_name,
                            time=m, time_units=time_units )
        ncgs_mean.add_grid( mean_grids[m, :,:], mean_short_name,
                            time=m, time_units=time_units )
        ncgs_sdev.add_grid( sdev_grids[m, :,:], sdev_short_name,
                            time=m, time_units=time_units )
                                           
    #-----------------------------------
    # Close all input and output files
    #-----------------------------------
    print('Wrote 12 grids for each statistic.')
    ncgs_var.close_file()   # variable not variance
    ncgs_min.close_file()
    ncgs_max.close_file()
    ncgs_mean.close_file()
    ncgs_sdev.close_file()
    print('Finished.')
    print(' ')
    
#   create_stat_grid_stacks()
#-------------------------------------------------------------------
def create_stat_grid_stacks_old( output_dir=None,
        case_prefix=None, hydro_var_name='discharge',
        var_file=None, min_file=None, max_file=None,
        mean_file=None, sdev_file=None,
        OVERWRITE_OK=False):

    #-------------------------------------------------------------
    # Note: This creates monthly grid stacks of several stats
    #       (e.g. min, max, mean, sdev) from a grid stack of a
    #       hydrologic variable.
    #       Spatial resolution remains the same.
    #       All input & output files are netCDF format.
    #-------------------------------------------------------------
    # https://en.wikipedia.org/wiki/Variance#Population_variance
    #-----------------------------------------------------------------
    # Option 1:  First, for each of 10 years, compute the 4 stats
    #            for each month, for each grid cell.  The resulting
    #            grid stack for each stat will then have 120 grids.
    #            Next, for each month, compute the 4 stats as an
    #            average over 10 years.  The resulting grid stack
    #            for each stat will then have 12 grids. 
    #-----------------------------------------------------------------
    # Option 2:  For each grid cell, compute the long-term (10-yr)
    #            stats (mean, sdev, min and max) using all values
    #            for that grid cell (regardless of month). There
    #            will then be 4 grids, 1 for each of the 4 statistics.
    #            Since discharge is naturally much lower (or higher)
    #            in certain months, this doesn't seem that helpful.
    #-----------------------------------------------------------------
    # Option 3:  For each grid cell, compute the long-term (10-yr)
    #            stats (mean, sdev, min, & max) for each day of the
    #            year.  Each stat will then be computed from just 10
    #            values, and the grid stack for each stat will have
    #            365 grids.
    #-----------------------------------------------------------------    
 

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
         time_units=time_units, time_res=time_res,
         OVERWRITE_OK=OVERWRITE_OK )
    #--------------------------------------------------------
    ncgs_max = ncgs_files.ncgs_file()  # new ncgs object
    ncgs_max.open_new_file( max_file,
         grid_info=grid_info, time_info=time_info,
         var_name=max_short_name, long_name=max_long_name, 
         units_name=var_units, dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res,
         OVERWRITE_OK=OVERWRITE_OK )
    #--------------------------------------------------------
    ncgs_mean = ncgs_files.ncgs_file()  # new ncgs object
    ncgs_mean.open_new_file( mean_file,
         grid_info=grid_info, time_info=time_info,
         var_name=mean_short_name, long_name=mean_long_name, 
         units_name=var_units, dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res,
         OVERWRITE_OK=OVERWRITE_OK )
    #--------------------------------------------------------
    ncgs_sdev = ncgs_files.ncgs_file()  # new ncgs object
    ncgs_sdev.open_new_file( sdev_file,
         grid_info=grid_info, time_info=time_info,
         var_name=sdev_short_name, long_name=sdev_long_name, 
         units_name=var_units, dtype=dtype, MAKE_RTI=True,
         time_units=time_units, time_res=time_res,
         OVERWRITE_OK=OVERWRITE_OK )
   
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
    month_num  = -1
    month_count = 0

    for time_index in range(n_grids):
        #----------------------------
        # Parse the datetime string
        #----------------------------
        datetime_str = datetimes[ time_index ]
        (y,m1,d,h,m2,s) = tu.split_datetime_str(datetime_str, ALL=True)
        
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
                                    time=month_count, time_units=time_units )
                ncgs_max.add_grid(  max_grid,  max_short_name,
                                    time=month_count, time_units=time_units )
                ncgs_mean.add_grid( mean_grid, mean_short_name,
                                    time=month_count, time_units=time_units )
                ncgs_sdev.add_grid( sdev_grid, sdev_short_name,
                                    time=month_count, time_units=time_units )
                month_count += 1    ####
                
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

    #--------------------------------
    # Save info for the final month
    # (Bug fix: 2022-02-16)
    #--------------------------------
    mean_grid     /= n_grids_per_month
    mean_sqr_grid /= n_grids_per_month
    sdev_grid[:] = np.sqrt( mean_sqr_grid - mean_grid**2)
    #----------------------------------------
    # Write statistic grids to netCDF files
    #----------------------------------------
    ncgs_min.add_grid(  min_grid,  min_short_name,
                        time=month_count, time_units=time_units )
    ncgs_max.add_grid(  max_grid,  max_short_name,
                        time=month_count, time_units=time_units )
    ncgs_mean.add_grid( mean_grid, mean_short_name,
                        time=month_count, time_units=time_units )
    ncgs_sdev.add_grid( sdev_grid, sdev_short_name,
                        time=month_count, time_units=time_units )
    month_count += 1
    print('Wrote', month_count, 'grids for each statistic.')
                                                         
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
    
#   create_stat_grid_stacks_old()
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
def create_day_count_grid_stack( output_dir=None, stats_dir=None,
           case_prefix=None, hydro_var_name='discharge',
           condition='low_discharge',
           var_file=None, mean_file=None, sdev_file=None,
           days_file=None,   # days_file is output
           OVERWRITE_OK=False):

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
    if (stats_dir is None):
        stats_dir = '~/stats1/'
        ## stats_dir = '~/stats/2005-01_to_2015-01_stats/'
        stats_dir = os.path.expanduser( stats_dir )
    
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
#     if (mean_file is None):
#         mean_file = case_prefix + '_2D-month-mean-' + s1 + '.nc'
#     if (sdev_file is None):
#         sdev_file = case_prefix + '_2D-month-sdev-' + s1 + '.nc'
#     if (days_file is None):
#         days_file = case_prefix + '_2D-' + out_var_name + '.nc'

    #------------------------------------------------
    # Construct default filenames using "stats" dir
    #------------------------------------------------
    if (mean_file is None):
        mean_file = case_prefix + '_2D-month-mean-' + s1 + '.nc'
        mean_file = stats_dir + mean_file
    if (sdev_file is None):
        sdev_file = case_prefix + '_2D-month-sdev-' + s1 + '.nc'
        sdev_file = stats_dir + sdev_file
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
         units_name=var_units, dtype='int16',   ################# TEST THIS !!
         ### units_name=var_units, dtype='float32',
         time_units='months', time_res='1', 
         MAKE_RTI=True, OVERWRITE_OK=OVERWRITE_OK)
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
                # Write days_grid for prev month to netCDF
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

            #-----------------------------------------------
            # Get mean & sdev grids for this new month
            # These should be long-term monthly stats.
            # 2nd arg. to get_grid() should be month num.
            # get_grid() will jump to the requested month.
            # Note: month_num is an integer in {1,...,12}
            # Assume that these grid stacks have 12 grids.
            #-----------------------------------------------
            try:
                mean_grid = ncgs_mean.get_grid( mean_short_name, month_num - 1 )
                sdev_grid = ncgs_sdev.get_grid( sdev_short_name, month_num - 1 )
                #-----------------------------------------------------------------
                # mean_grid = ncgs_mean.get_grid( mean_short_name, month_count )
                # sdev_grid = ncgs_sdev.get_grid( sdev_short_name, month_count )
            except:
                # This should never occur now.  ################
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
            w = (grid > threshold)
        elif (condition == 'flooding'):
            ### w = (grid > 0)   ################# FLOOD THRESHOLD
            w = (grid > 0.1)  # [meters] FLOOD THRESHOLD
        elif ('incr_' in condition):
            w = (grid > last_grid)
            np.copyto(last_grid, grid)  # in-place
            ## last_grid[:] = grid   # in-place
        elif ('decr_' in condition):
            w = (grid < last_grid)
            np.copyto(last_grid, grid)  # in-place
            ## last_grid[:] = grid  # in-place
        else:
            print('ERROR: Unknown condition.')
            print()
        #----------------------------------
        detect_grid[w] = True

    #-------------------------------------------
    # Write days_grid for final month to netCDF
    # (Bug fix: 2022-02-16)
    #-------------------------------------------
    ncgs_days.add_grid( days_grid, out_var_name,
                        time=month_count,
                        time_units=time_units )
    month_count += 1
    print('Wrote', month_count, 'grids to file.')
                                                    
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
def create_population_grid_file( pop_dir=None, pop_file=None,
           out_xres_sec=None, out_yres_sec=None,
           out_bounds=None, case_prefix=None,
           grid_info=None, time_info=None,
           global_pop_dir=None, global_pop_file=None):
           
    #---------------------------   
    # Check required arguments
    #---------------------------
    if (out_xres_sec is None):
        print('ERROR: out_xres_sec argument is required.')
        return
    if (out_yres_sec is None):
        print('ERROR: out_yres_sec argument is required.')
        return
    if (out_bounds is None):
        print('ERROR: out_bounds argument is required.')
        return
    if (grid_info is None):
        print('ERROR: grid_info argument is required.')
        return
    if (time_info is None):
        print('ERROR: time_info argument is required.')
        return
                                
    #------------------------   
    # Set default filenames
    #------------------------
    if (case_prefix is None):
        case_prefix = 'Test1'
    if (pop_dir is None):
        pop_dir = ''
    if (pop_file is None):
        prefix   = case_prefix + '_2D_gpw_v4_pop_count_rev11_2020_'
        xres_str = str(int(out_xres_sec))
        suffix   = '_sec.tif'
        pop_file = prefix + xres_str + suffix
    #----------------------------------------------
    if (global_pop_dir is None):
        global_pop_dir = '/Users/peckhams/Data/GPW-v4/'
    if (global_pop_file is None):
        global_pop_file = 'gpw_v4_population_count_rev11_2020_30_sec.tif'
    #---------------------------------------------------
    pop_file_nc     = pop_file.replace('.tif', '.nc')
    global_pop_path = global_pop_dir + global_pop_file           
    pop_path_tif    = pop_dir + pop_file      # GeoTIFF
    pop_path_nc     = pop_dir + pop_file_nc   # netCDF

    #---------------------------------------------------   
    # Read global population grid (GPW-v4 as geotiff).
    # Clip and resample as necessary.
    # Save new population grid in GeoTIFF format.
    #---------------------------------------------------
    if (os.path.exists( pop_path_tif )):
        print('ERROR: A file already exists with the name:')
        print('   ' + pop_path_tif )
        print('   Skipping...')
    else:
        print('Clipping & resampling global GPW-v4 pop. count grid...')
        rg.regrid_geotiff(in_file=global_pop_path,
            out_file=pop_path_tif, out_bounds=out_bounds,
            out_xres_sec=out_xres_sec, out_yres_sec=out_yres_sec,
            RESAMPLE_ALGO='bilinear', REPORT=True)
        print('Finished writing population grid as GeoTIFF.')

    #--------------------------------------     
    # Write new population grid to netCDF
    #--------------------------------------
    if (os.path.exists( pop_path_nc )):
        print('ERROR: A file already exists with the name:')
        print('   ' + pop_path_nc )
        print('   Skipping...')
    else:
        print('Writing new population grid to netCDF...')
        ncgs_pop = ncgs_files.ncgs_file()  # new ncgs object
        dtype    = 'float32'  # (actual dtype in GPW-v4)
        # dtype  = 'int32'
        #--------------------------------------------------
        ncgs_pop.open_new_file( pop_path_nc,
             grid_info=grid_info, time_info=time_info,
             var_name='pop_count', long_name='population__count', 
             units_name='1', dtype=dtype, MAKE_RTI=True,
             time_units='years', time_res=1 )  #################
        #---------------------------------------------------
        # GPW-v4 data type is "float32".
        # GPW-v4 nodata value = -3.402823e+38 (Water Mask)
        #---------------------------------------------------
        pop_grid = rg.read_geotiff( pop_path_tif )
        w1 = (pop_grid < 0)   # boolean array
        pop_grid[ w1 ] = 0.0
        #----------------------------------------
        ncgs_pop.add_grid( pop_grid, 'pop_count',
                           time=0, time_units='years' )
        ncgs_pop.close()
        print('Finished writing population grid as netCDF.')
        print()
           
#   create_population_grid_file()
#-------------------------------------------------------------------
def create_affected_population_grid_stack( output_dir=None,
           pop_dir=None, pop_file=None, case_prefix=None,
           # Next 2 not needed if pop_file exists.
           global_pop_dir=None, global_pop_file=None, 
           condition='low_discharge', day_count_file=None,
           out_file=None, out_var_name=None,
           OVERWRITE_OK=False):

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
    if (global_pop_dir is None):
        global_pop_dir = '/Users/peckhams/Data/GPW-v4/'
    elif (global_pop_dir[-1] != os.sep):
        global_pop_dir += os.sep
    #-------------------------------------------------------
    if (global_pop_file is None):
        global_pop_file = 'gpw_v4_population_count_rev11_2020_30_sec.tif'
    #-------------------------------------------------------
    if (case_prefix is None):
        case_prefix = 'Test1'
    #------------------------------------------------------- 
    if (pop_dir is None):
        pop_dir = ''     # usually will be bmi.misc_dir
    elif (pop_dir[-1] != os.sep):
        pop_dir += os.sep
    #------------------------------------------------------- 
    if (pop_file is None):
        pass  # A default will be constructed below
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
         time_units=time_units, time_res=time_res,
         OVERWRITE_OK=OVERWRITE_OK ) 
                   
    #-------------------------------------------
    # Get out_xres_sec, out_yres_sec, & bounds
    #-------------------------------------------
    out_xres_sec = grid_info.xres  # [arcseconds]
    out_yres_sec = grid_info.yres  # [arcseconds]
    out_bounds = ncgs_days.ncgs_unit.geospatial_bounds     #### CHECK ORDERING
    ######################################################
    out_grid = np.zeros( (nrows, ncols), dtype=dtype)

    #------------------------------------------------
    # Construct default name for pop_file (GeoTIFF)
    # File may already exist
    #------------------------------------------------
    if (pop_file is None):
        prefix   = case_prefix + '_2D_gpw_v4_pop_count_rev11_2020_'
        xres_str = str(int(out_xres_sec))
        suffix   = '_sec.tif'
        pop_file = prefix + xres_str + suffix

    #--------------------------------------------------
    # Create a population grid file for this region ?
    #--------------------------------------------------
    # Set pop_file=None to use the default filename.
    #--------------------------------------------------
    pop_path_tif = (pop_dir + pop_file)
    if not(os.path.exists( pop_path_tif )):
        create_population_grid_file( pop_dir=pop_dir,
            pop_file=pop_file, case_prefix=case_prefix,
            global_pop_dir=global_pop_dir,  
            global_pop_file=global_pop_file,  
            out_xres_sec=out_xres_sec, out_bounds=out_bounds,
            out_yres_sec=out_yres_sec,
            grid_info=grid_info, time_info=time_info )

    #---------------------------------------------------
    # GPW-v4 data type is "float32".
    # GPW-v4 nodata value = -3.402823e+38 (Water Mask)
    #---------------------------------------------------
    pop_grid = rg.read_geotiff( pop_path_tif )
    w1 = (pop_grid < 0)   # boolean array
    pop_grid[ w1 ] = 0.0
                                                    
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
    print('Wrote', n_grids, 'grids to file.')   
    print('Finished computing indicator grid stack.')
    print(' ')
    
#   create_affected_population_grid_stack()
#-------------------------------------------------------------------


