#   
#  Copyright (c) 2020, Scott D. Peckham
#
#  Note: This file contains a set of functions for visualizing the
#        contents of TopoFlow output files in netCDF format.
#
#  May 2020.  Moved all routines from Jupyter notebook called
#             TopoFlow_Visualization.ipynb to here.
#             Tested all of them in a new Jupyter notebook called
#             TopoFlow_Visualization2.ipynb.
#
#--------------------------------------------------------------------
#
#  Define some stretch functions for 2D color images:
#  histogram_equalize()
#  power_stretch0()
#  power_stretch1()
#  power_stretch2()
#  power_stretch3()
#  log_stretch()
#  linear_stretch()
#  stretch_grid()
#
#  Define functions to show grids as color images:
#  read_grid_from_nc_file()
#  show_grid_as_image()
#  save_grid_stack_as_images()
#
#  Define some plotting functions:
#  plot_time_series()
#  plot_z_profile()
#  save_profile_series_as_images()
#
#  Create movies from set of images:
#     (works for grid images, profile images, etc.)
#  create_movie_from_images()
#
#  From Richards 1D Equation Jupyter notebook
#  plot_data()
#  plot_soil_profile()
#
#--------------------------------------------------------------------
# import os.path
# import shutil

import glob, os
import numpy as np
import matplotlib.pyplot as plt
import imageio

from topoflow.utils import ncgs_files
from topoflow.utils import ncts_files
from topoflow.utils import ncps_files

from topoflow.utils import rtg_files
from topoflow.utils import rts_files

#--------------------------------------------------------------------
def histogram_equalize( grid, PLOT_NCS=False):

    #  https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html
    (hist, bin_edges) = np.histogram( grid, bins=256)
    # hmin = hist.min()
    # hmax = hist.max()

    cs  = hist.cumsum()
    ncs = (cs - cs.min()) / (cs.max() - cs.min())
    ncs.astype('uint8');
    ############## ncs.astype('uint8') # no semi-colon at end ??????????
    if (PLOT_NCS):
        plt.plot( ncs )

    flat = grid.flatten()
    flat2 = np.uint8( 255 * (flat - flat.min()) / (flat.max() - flat.min()) )
    grid2 = ncs[ flat2 ].reshape( grid.shape )
    return grid2

#   histogram_equalize()
#--------------------------------------------------------------------
def power_stretch0( grid, p ):

    gmin = grid.min()
    gmax = grid.max()
    norm = (grid - gmin) / (gmax - gmin)
    
    return norm**p
    
#   power_stretch0()
#--------------------------------------------------------------------
def power_stretch1( grid, p ):
    return grid**p
    
#   power_stretch1()
#--------------------------------------------------------------------
def power_stretch2( grid, a=1000, b=0.5):
    # Note: Try a=1000 and b=0.5
    gmin = grid.min()
    gmax = grid.max()
    norm = (grid - gmin) / (gmax - gmin)
    return (1 - (1 + a * norm)**(-b))
    
#   power_stretch2()
#--------------------------------------------------------------------
def power_stretch3( grid, a=1, b=2):
    # Note:  Try a=1, b=2 (shape of a quarter circle)
    gmin = grid.min()
    gmax = grid.max()
    norm = (grid - gmin) / (gmax - gmin)
    return (1 - (1 - norm**a)**b)**(1/b)
    
#   power_stretch3()
#--------------------------------------------------------------------
def log_stretch( grid, a=1 ):
    return np.log( (a * grid) + 1 )
    
#   log_stretch()
#--------------------------------------------------------------------
def linear_stretch( grid ):

    gmin = grid.min()
    gmax = grid.max()
    norm = (grid - gmin) / (gmax - gmin)
    return norm
   
#   linear_stretch()
#--------------------------------------------------------------------
def stretch_grid( grid, stretch, a=1, b=2, p=0.5 ):

    name = stretch
    if   (name == 'hist_equal'):
        grid2 = histogram_equalize( grid, PLOT_NCS=False)    
    elif (name == 'linear'):
        grid2 = linear_stretch(grid)
    elif (name == 'log'):
        grid2 = log_stretch( grid, a=a )
    elif (name == 'power'):
        grid2 = power_stretch0( grid, p=p )
    elif (name == 'power1'): 
        # Try:  p = 0.3   
        grid2 = power_stretch1( grid, p)
    elif (name == 'power2'):
        # Try:  a=1000, b=0.5.
        grid2 = power_stretch2( grid, a=a, b=b )
    elif (name == 'power3'):        
        # Try:  a=1, b=2.
        grid2 = power_stretch3( grid, a=a, b=b)
    else:
        print('### SORRY, Unknown stretch =', name)
        return grid

    return grid2
 
#   stretch_grid()
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def read_grid_from_nc_file( nc_file, time_index=1, REPORT=True ):

    # Typical 2D nc files
    # nc_file = case_prefix + '_2D-Q.nc'
    # nc_file = case_prefix + '_2D-d-flood.nc'

    if ('_2D' not in nc_file):
        print('ERROR: This function is only for TopoFlow "2D" files.') 
        return
            
    ncgs = ncgs_files.ncgs_file()
    ncgs.open_file( nc_file )
    var_name_list = ncgs.get_var_names()
    if (REPORT):
        print('var_names in netCDF file =' )
        print( var_name_list )    
    ## var_index = 1   # 0=time, 1=X, 2=Y, 3=V  ###############
    var_index = 3   # 0=time, 1=X, 2=Y, 3=V  ###############
    var_name  = var_name_list[ var_index ]
    long_name = ncgs.get_var_long_name( var_name )
    var_units = ncgs.get_var_units( var_name )
    n_grids   = ncgs.ncgs_unit.variables[ var_name ].n_grids

    if (REPORT):
        print('long_name =', long_name)
        print('var_name  =', var_name)
        print('var_units =', var_units)
        print('n_grids   =', n_grids)

    #--------------------------------------------
    # Use these to set "extent" in plt.imshow()
    #--------------------------------------------
    minlon = ncgs.ncgs_unit.variables['X'].geospatial_lon_min
    maxlon = ncgs.ncgs_unit.variables['X'].geospatial_lon_max
    minlat = ncgs.ncgs_unit.variables['Y'].geospatial_lat_min
    maxlat = ncgs.ncgs_unit.variables['Y'].geospatial_lat_max
    extent = [minlon, maxlon, minlat, maxlat]
    
    #----------------------------------------------
    # Read grid from nc_file for given time_index
    #----------------------------------------------
    grid = ncgs.get_grid( var_name, time_index )
    
    if (REPORT):
        print( 'extent = ')
        print( extent )
        print( 'grid shape =', grid.shape )
        print( 'min(grid)  =', grid.min() )
        print( 'max(grid)  =', grid.max() )

    ncgs.close_file()
    return (grid, long_name, extent)
    
#   read_grid_from_nc_file()
#--------------------------------------------------------------------
def show_grid_as_image( grid, long_name, extent=None,
                        cmap='rainbow',
                        stretch='power3',
                        a=1, b=2, p=0.5,
                        NO_SHOW=False, im_file=None,
                        xsize=8, ysize=8, dpi=None): 

    # Note:  extent = [minlon, maxlon, minlat, maxlat]
    
    #-------------------------
    # Other color map names
    #--------------------------------------------
    # hsv, jet, gist_rainbow (reverse rainbow),
    # gist_ncar, gist_stern
    #--------------------------------------------    

    #--------------------------
    # Other stretch functions
    #--------------------------
    grid2 = stretch_grid( grid, stretch, a=a, b=b, p=p )
#     if (stretch == 'power_stretch3'):
#         grid2 = power_stretch3( grid, a=0.5 )
#     elif (stretch == 'power_stretch1a'):   
#         grid2 = power_stretch1( grid, 0.5)
#     elif (stretch == 'power_stretch1b'):
#         grid2 = power_stretch1( grid, 0.2)
#     elif (stretch == 'power_stretch2'):
#         grid2 = power_stretch2( grid )
#     elif (stretch == 'log_stretch'):
#         grid2 = log_stretch( grid )
#     elif (stretch == 'hist_equal'):
#         grid2 = histogram_equalize( grid, PLOT_NCS=True)
#     else:
#         print('SORRY, Unknown stretch =', stretch)
#         return 

    #----------------------------
    # Set up and show the image
    #----------------------------
    # figure = plt.figure(1, figsize=(xsize, ysize))
    fig, ax = plt.subplots( figsize=(xsize, ysize), dpi=dpi)
    im_title = long_name.replace('_', ' ').title()
    ax.set_title( im_title )
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')

    gmin = grid2.min()
    gmax = grid2.max()

    im = ax.imshow(grid2, interpolation='nearest', cmap=cmap,
                   vmin=gmin, vmax=gmax, extent=extent)

    #--------------------------------------------------------        
    # NOTE!  Must save before "showing" or get blank image.
    #        File format is inferred from extension.
    #        e.g. TMP_Image.png, TMP_Image.jpg.
    #--------------------------------------------------------
    if (im_file is not None):  
        plt.savefig( im_file )
    if not(NO_SHOW):
        plt.show()
 
    plt.close()
        
#   Information on matplotlib color maps
#   https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
# 
#   Information on matplotlib.pyplot.imshow
#   https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.imshow.html
# 
#   Information on matplotlib.pyplot.savefig
#   https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.savefig.html
# 
#   plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w',
#               orientation='portrait', papertype=None, format=None,
#               transparent=False, bbox_inches=None, pad_inches=0.1,
#               frameon=None, metadata=None)
  
#   show_grid_as_image()
#--------------------------------------------------------------------
def save_grid_stack_as_images( nc_file, png_dir, extent=None,
                               xsize=6, ysize=6, dpi=192,
                               REPORT=True):

    # Example nc_files:
    # nc_file = case_prefix + '_2D-Q.nc'
    # nc_file = case_prefix + '_2D-d-flood.nc'

    if ('_2D' not in nc_file):
        print('ERROR: This function is only for TopoFlow "2D" files.') 
        return

    ncgs = ncgs_files.ncgs_file()        
    ncgs.open_file( nc_file )
    var_name_list = ncgs.get_var_names()
    var_index = 3   # 0=time, 1=X, 2=Y, 3=V  ###############
    var_name  = var_name_list[ var_index ]
    long_name = ncgs.get_var_long_name( var_name )
    n_grids = ncgs.ncgs_unit.variables[var_name].n_grids

    im_title = long_name.replace('_', ' ').title()
    im_file_prefix = 'TF_Grid_Movie_Frame_'
    time_pad_map = {1:'0000', 2:'000', 3:'00', 4:'0', 5:''}
    cmap = 'rainbow'

    if (REPORT):
        print('Creating images from grid stack in nc_file:')
        print('  ' + nc_file )
        print('  ' + 'var name  =', var_name)
        print('  ' + 'long name =', long_name)
        print('  ' + 'n_grids   =', n_grids)
        print('This may take a few minutes.')
        print('Working...')
        
    time_index = 0
    while (True):
        # print('time index =', time_index )
        try:
            grid = ncgs.get_grid( var_name, time_index )
        except:
            break
        time_index += 1

        #----------------------------------------    
        # Build a filename for this image/frame
        #----------------------------------------
        tstr = str(time_index)
        pad = time_pad_map[ len(tstr) ]
        time_str = (pad + tstr)
        im_file = im_file_prefix + time_str + '.png' 
        im_file = (png_dir + '/' + im_file)
                
        show_grid_as_image( grid, long_name, cmap='rainbow',
                            stretch='power_stretch3',
                            extent=extent,
                            NO_SHOW=True, im_file=im_file,
                            xsize=xsize, ysize=ysize, dpi=dpi )

    ncgs.close_file()
    tstr = str(time_index)
    print('Finished saving ' + tstr + ' images to PNG files.')

#   save_grid_stack_as_images()
#--------------------------------------------------------------------
def plot_time_series(nc_file, output_dir=None, marker=',',
                     REPORT=True, xsize=11, ysize=6):

    # Example nc_files:
    # nc_file = case_prefix + '_0D-Q.nc'
    # nc_file = case_prefix + '_0D-d-flood.nc'

    #----------------------------------------------------
    # Possible plot point markers:
    # https://matplotlib.org/3.1.1/api/markers_api.html
    #----------------------------------------------------
    # marker = ','  # pixel
    # marker = '.'  # point (small circle)
    # marker = 'o'  # circle
    # marker = '+'
    # marker = 'x'

    if ('_0D' not in nc_file):
        print('ERROR: This function is only for TopoFlow "OD" files.') 
        return
       
    if (output_dir is not None):
        os.chdir( output_dir )

    ncts = ncts_files.ncts_file()
    ncts.open_file( nc_file )
    var_name_list = ncts.get_var_names()
    lon_list      = ncts.get_var_lons()
    lat_list      = ncts.get_var_lats()

    var_index = 1
    var_name = var_name_list[ var_index ]

    if (REPORT):
        print( 'var_names in netCDF file =' )
        print( var_name_list )
        print( 'var longitudes =')
        print( lon_list )
        print( 'var latitudes =')
        print( lat_list )
        print()

    (series, times) = ncts.get_series( var_name )
    long_name = series.long_name
    v_units   = series.units
    t_units   = times.units
    values    = np.array( series )
    times     = np.array( times )
    # values = series[:]   # also works
    # times  = times[:]    # also works

    if (t_units == 'minutes'):
        # times = times / 60.0
        # t_units = 'hours'
        times = times / (60.0 * 24)
        t_units = 'days'

    # ymin = values.min()
    ymin = 0.0
    ymax = values.max()

    figure = plt.figure(1, figsize=(xsize, ysize))
    # fig, ax = plt.subplots( figsize=(xsize, ysize))

    y_name = long_name.replace('_', ' ').title()

    plt.plot( times, values, marker=marker)
    plt.xlabel( 'Time' + ' [' + t_units + ']' )
    plt.ylabel( y_name + ' [' + v_units + ']' )
    plt.ylim( np.array(ymin, ymax) )

    plt.show()

    ncts.close_file()
    
#   plot_time_series()
#--------------------------------------------------------------------
def plot_z_profile(nc_file, time_index=50,
                   output_dir=None, marker=',',
                   REPORT=True, xsize=11, ysize=6):

    # Example nc_files:
    # nc_file = case_prefix + '_1D-q.nc'
    # nc_file = case_prefix + '_1D-p.nc'
    # nc_file = case_prefix + '_1D-K.nc'
    # nc_file = case_prefix + '_1D-v.nc'

    #----------------------------------------------------
    # Possible plot point markers:
    # https://matplotlib.org/3.1.1/api/markers_api.html
    #----------------------------------------------------
    # marker = ','  # pixel
    # marker = '.'  # point (small circle)
    # marker = 'o'  # circle
    # marker = '+'
    # marker = 'x'

    if ('_1D' not in nc_file):
        print('ERROR: This function is only for TopoFlow "1D" files.') 
        return
       
    if (output_dir is not None):
        os.chdir( output_dir )
 
    ncps = ncps_files.ncps_file()
    ncps.open_file( nc_file )
    var_name_list = ncps.get_var_names()
    lon_list      = ncps.get_var_lons()
    lat_list      = ncps.get_var_lats()

    var_index = 2  #### 0 = time, 1 = z, 2 = var
    var_name = var_name_list[ var_index ]

    if (REPORT):
        print( 'var_names in netCDF file =' )
        print( var_name_list )
        print( 'var longitudes =')
        print( lon_list )
        print( 'var latitudes =')
        print( lat_list )
        print()

    # (profile, z, time) = ncps.get_profile(var_name, time_index)

    (profiles, z, times) = ncps.get_profiles(var_name)
    long_name = profiles.long_name
    v_units   = profiles.units
    z_units   = z.units
    # t_units   = times.units

    profile    = profiles[ time_index]
    time       = times[ time_index ]
    values     = np.array( profile )
    z_values   = np.array( z )
    # times     = np.array( times )
    # values = profile[:]   # also works
    # times  = times[:]    # also works

    # if (t_units == 'minutes'):
    #     # times = times / 60.0
    #     # t_units = 'hours'
    #     times = times / (60.0 * 24)
    #     t_units = 'days'

    # ymin = 0.0
    ymin = values.min()
    ymax = values.max()

    figure = plt.figure(1, figsize=(11,6))
    # fig, ax = plt.subplots( figsize=(11,6))

    y_name = long_name.replace('_', ' ').title()

    plt.plot( z_values, values, marker=marker)
    plt.xlabel( 'Depth' + ' [' + z_units + ']' )
    plt.ylabel( y_name + ' [' + v_units + ']' )
    plt.ylim( np.array(ymin, ymax) )

    plt.show()

    ncps.close_file()
         
#   plot_z_profile()
#--------------------------------------------------------------------
def save_profile_series_as_images(nc_file, png_dir=None,
                 ymin=None, ymax=None, marker=',', REPORT=True,
                 xsize=11, ysize=6, dpi=192):

    # Examples of nc_files:
    # nc_file = case_prefix + '_1D-q.nc'
    # nc_file = case_prefix + '_1D-p.nc'
    # nc_file = case_prefix + '_1D-v.nc'
    # nc_file = case_prefix + '_1D-K.nc'

    if (png_dir is None):
        print('ERROR, PNG directory is not set.')
        return
    
    # If ymin or ymax is set, it is used for all plots
    ALL_SAME_YMIN = (ymin is not None)
    ALL_SAME_YMAX = (ymax is not None)
         
    ncps = ncps_files.ncps_file()
    ncps.open_file( nc_file )
    var_name_list = ncps.get_var_names()
    var_index = 2   #### 0 = time, 1 = z, 2 = var
    var_name  = var_name_list[ var_index ]
    long_name = ncps.get_var_long_name( var_name )
    ##  svo_name = ncps.get_var_svo_name( var_name )

    y_name = long_name.replace('_', ' ').title()

    im_title = long_name.replace('_', ' ').title()
    im_file_prefix = 'TF_Profile_Movie_Frame_'
    time_pad_map = {1:'0000', 2:'000', 3:'00', 4:'0', 5:''}
    ## cmap = 'rainbow'

    #---------------------------------
    # Read all of the profiles, etc.
    #---------------------------------
    try:
        (profiles, z, times) = ncps.get_profiles( var_name )
        n_profiles = len( profiles )  ########
        long_name = profiles.long_name
        v_units   = profiles.units
        z_values  = np.array( z )
        z_units   = z.units
        # t_units   = times.units
    except:
        print('ERROR, Could not read profiles from file:')
        print( nc_file )
        return
                   
    if (REPORT):
        print('Creating images from z profiles in nc_file:')
        print('  ' + nc_file )
        print('  ' + 'var name  =', var_name)
        print('  ' + 'long_name =', long_name )
        ## print('  ' + 'svo_name  =', svo_name )
        print('  ' + 'Number of profiles =', n_profiles )
        print('This may take a few minutes.')
        print('Working...')

    for time_index in range(n_profiles):
        # print('time index =', time_index )
        profile   = profiles[ time_index ]
        time      = times[ time_index ]
        time_index += 1
 
        values = np.array( profile )
        if not(ALL_SAME_YMIN):
            ymin = values.min()
        if not(ALL_SAME_YMAX):
            ymax = values.max()

        fig, ax = plt.subplots( figsize=(xsize, ysize), dpi=dpi) 
        ax.set_title( im_title )   
        plt.plot( z_values, values, marker=marker)
        plt.xlabel( 'Depth' + ' [' + z_units + ']' )
        plt.ylabel( y_name + ' [' + v_units + ']' )
        plt.ylim( ymin, ymax )

        #-------------------------------------------        
        # Using np.array() like this doesn't work.
        #-------------------------------------------
        ### plt.ylim( np.array(ymin, ymax) )

        #------------------------------------------        
        # Don't show each plot as it is generated
        #------------------------------------------
        ## plt.show()

        #------------------------------------------
        # Build a filename for this image/frame
        #------------------------------------------
        tstr = str(time_index)
        pad = time_pad_map[ len(tstr) ]
        time_str = (pad + tstr)
        im_file = im_file_prefix + time_str + '.png' 
        im_file = (png_dir + '/' + im_file)

        plt.savefig( im_file )
        plt.close()

    ncps.close_file()
    tstr = str(time_index)
    print('Finished saving ' + tstr + ' images to PNG files.')

#   save_profile_series_as_images()
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def create_movie_from_images( mp4_file, png_dir, fps=10, REPORT=True):

    #----------------------------------------
    # png_dir  = directory with PNG files
    # mp4_file = case_prefix + '_Movie.mp4'
    # fps      = frames per second
    #----------------------------------------
    im_file_list = sorted( glob.glob( png_dir + '/*.png' ) )
    n_frames = len( im_file_list )

    if (REPORT):
        print('Creating movie from', n_frames, 'PNG files.')
        ## print('This may take a few minutes.')
        print('Working...')
    writer = imageio.get_writer( mp4_file, fps=fps )

    for im_file in im_file_list:
        writer.append_data(imageio.imread( im_file ))
    writer.close()
   
    if (REPORT): 
        print('Finished creating movie, MP4 format.')
        print('  ' + mp4_file)
        print()

#   create_movie_from_images()
#--------------------------------------------------------------------
def plot_data( x, y, xmin=None, xmax=None, ymin=None, ymax=None,
               x_name='x', x_units='', marker=',', 
               y_name='y', y_units='',
               x_size=8,   y_size=4):

    figure = plt.figure(1, figsize=(x_size, y_size))
    # fig, ax = plt.subplots( figsize=(x_size, y_size))

    # Set the plot point marker
    # https://matplotlib.org/3.1.1/api/markers_api.html
    # marker = ','  # pixel
    # marker = '.'  # point (small circle)
    # marker = 'o'  # circle
    # marker = '+'
    # marker = 'x'

    #if (ymin is None):
    #    ymin = y.min()
    #if (ymax is None):
    #    ymax = y.max()
    #if (ymax - ymin < 0.1):
    #    ymin = ymin - 0.5
    #    ymax = ymin + 0.5

    # x_name2 = x_name.replace('_', ' ').title()
    # y_name2 = y_name.replace('_', ' ').title()
        
    plt.plot( x, y, marker=marker)
    plt.xlabel( x_name + ' [' + x_units + ']' )
    plt.ylabel( y_name + ' [' + y_units + ']' )
    
    plt.ylim( ymin, ymax )
    plt.xlim( xmin, xmax )
    #-------------------------------------
    # This may be necessary depending on
    # the data type of ymin, ymax
    #-------------------------------------
    ## plt.ylim( np.array([ymin, ymax]) )
    ## plt.xlim( np.array([xmin, xmax]) )
    plt.show()

#   plot_data()
#----------------------------------------------------------------------------   
def plot_soil_profile( z, var, var_name='theta', qs=None,
                       MMPH=True, SWAP_AXES=False ):

    if (var_name == 'K') or (var_name == 'v'):
        if (MMPH):
            y = var * (1000.0 * 3600.0)
            y_units = 'mmph'
        else:
            y = var
            y_units = 'm/s'
    else:
        y = var
    #----------------------------------------------
    x       = z
    ymin    = y.min()   # default
    ymax    = y.max()   # default
    x_name  = 'z'
    x_units = 'm'
    y_name  = var_name
    #----------------------------------------------
    if (var_name == 'theta'):
        y_units = '1'
        ymin = 0.0
        ymax = qs + 0.01
    elif (var_name == 'psi'):
        y_units = 'm'
    elif (var_name == 'K'):
        pass
    elif (var_name == 'v'):
        pass
    else:
        y_units = '1'
    #-----------------------------------------
    if not(SWAP_AXES):
        plot_data( x, y, ymin=ymin, ymax=ymax,
                   x_name=x_name, x_units=x_units,
                   y_name=y_name, y_units=y_units)
    else:
        plot_data( y, -x, xmin=ymin, xmax=ymax,
                   x_name=y_name, x_units=y_units,
                   y_name=x_name, y_units=x_units)

#   plot_soil_profile()
#----------------------------------------------------------------------------
   

   