#   
#  Copyright (c) 2020-2023, Scott D. Peckham
#
#  Note: This file contains a set of functions for visualizing the
#        contents of output files in netCDF format
#        (e.g. TopoFlow or Stochastic Conflict Model)
#
#  Sep 2023.  Added xticks and yticks to plot_data().
#  Aug 2023.  Search for 2023-08-24 to see fixed bugs.
#             Added rtg_type keyword to read_and_show_rtg().
#  Feb 2022.  create_visualization_files -> create_media_files.
#  Oct 2021.  create_visualization_files().
#  Sep 2021.  Added LAND_SEA_BACKDROP option: show_grid_as_image() 
#  May 2020.  Moved all routines from Jupyter notebook called
#             TopoFlow_Visualization.ipynb to here.
#             Tested all of them in the Jupyter notebook.
#
#--------------------------------------------------------------------
#
#  Define some stretch functions for 2D color images:
#  normalize_grid()
#  histogram_equalize()
#  power_stretch0()
#  power_stretch1()
#  power_stretch2()
#  power_stretch3()
#  log_stretch()
#  linear_stretch()
#  tanh_stretch()     # (2022-05-02, for slope grids)
#  stretch_grid()
#
#  Define functions to show grids as color images:
#  read_grid_from_nc_file()
#  read_and_show_rtg()
#  show_grid_as_image()
#  save_grid_stack_as_images()
#  save_rts_as_images()
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
#  Next function will be called from a Dojo Docker container
#      after topoflow_driver.finalize() to create images
#      and movies from netCDF output files.
#
#  get_image_dimensions()   2022-05
#  create_media_files()     2021-10
#  create_output_movies()   2022-02-23  (separate function)
#  create_stat_movies()     2022-02-23  (separate function)
#  create_met_movies()      2022-02-17  (separate function)
#  delete_png_files()       2021-10
#
#--------------------------------------------------------------------
# import os.path
# import shutil

import glob, os, os.path, time
import datetime as dt  #########
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates   ##### 2022-04-30
from matplotlib import cm
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable  #########
import imageio

from topoflow.utils import ncgs_files
from topoflow.utils import ncts_files
from topoflow.utils import ncps_files
from topoflow.utils import rtg_files
from topoflow.utils import rts_files
from topoflow.utils import time_utils as tu

#--------------------------------------------------------------------
def normalize_grid( grid ): 

    gmin = grid.min()
    gmax = grid.max()

    if (gmin != gmax):
        norm = (grid - gmin) / (gmax - gmin)
    else:
        # Avoid divide by zero
        norm = np.zeros( grid.shape, dtype=grid.dtype )
    return norm

#   normalize_grid()
#--------------------------------------------------------------------
def histogram_equalize( grid, PLOT_NCS=False):

    #  https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html
    (hist, bin_edges) = np.histogram( grid, bins=256)
    # hmin = hist.min()
    # hmax = hist.max()

    cs  = hist.cumsum()
    ncs = (cs - cs.min()) / (cs.max() - cs.min())
    ncs.astype('uint8')

    if (PLOT_NCS):
        plt.plot( ncs )

    flat = grid.flatten()
    if (flat.max() != flat.min()):
        flat2 = np.uint8( 255 * (flat - flat.min()) / (flat.max() - flat.min()) )
        grid2 = ncs[ flat2 ].reshape( grid.shape )
    else:
        flat2 = np.zeros( flat.size, dtype='uint8' )
        grid2 = ncs[ flat2 ].reshape( grid.shape )

    return grid2

#   histogram_equalize()
#--------------------------------------------------------------------
def power_stretch0( grid, p ):

    norm = normalize_grid( grid )
    
    return norm**p
    
#   power_stretch0()
#--------------------------------------------------------------------
def power_stretch1( grid, p ):
    return grid**p
    
#   power_stretch1()
#--------------------------------------------------------------------
def power_stretch2( grid, a=1000, b=0.5):

    # Note: Try a=1000 and b=0.5
    norm = normalize_grid( grid )
    return (1 - (1 + a * norm)**(-b))
    
#   power_stretch2()
#--------------------------------------------------------------------
def power_stretch3( grid, a=1, b=2):

    # Note:  Try a=1, b=2 (shape of a quarter circle)
    norm = normalize_grid( grid )
    return (1 - (1 - norm**a)**b)**(1/b)
    
#   power_stretch3()
#--------------------------------------------------------------------
def log_stretch( grid, a=1 ):
    return np.log( (a * grid) + 1 )
    
#   log_stretch()
#--------------------------------------------------------------------
def linear_stretch( grid ):

    norm = normalize_grid( grid )
    return norm
   
#   linear_stretch()
#--------------------------------------------------------------------
def tanh_stretch( grid, a=2 ):

    # Seems good for slope grids
    #-----------------------------------
    # This stretch maps [0,1] to [0,1]
    # and has y[0]=0, y[1]=1.
    #-----------------------------------
    norm = normalize_grid( grid )
    denom = np.tanh( a )
    return np.tanh( a * norm ) / denom
       
#   tanh_stretch()
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
    elif (name == 'tanh'):
        grid2 = tanh_stretch( grid, a=a)
    else:
        print('### SORRY, Unknown stretch =', name)
        return grid

    return grid2
 
#   stretch_grid()
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def read_grid_from_nc_file( nc_file, time_index=0, REPORT=True ):

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

    var_name_list = ncgs.get_var_names( no_dim_vars=True )  ####
    var_index = 0   # (dim vars are now excluded)
    var_name  = var_name_list[ var_index ]
    
    #----------------------------         
    # Determine valid var_index
    #-----------------------------------------
    # Old: 0=time, 1=X, 2=Y, 3=V
    # New: 0=time, 1=datetime, 2=X, 3=Y, 4=V
    #-----------------------------------------
    # var_index = 1
    # other_vars = ['time','datetime','X','Y','Z']
    # while (True):
    #     var_name = var_name_list[ var_index ]
    #     if (var_name not in other_vars):
    #         break
    #     var_index += 1    
    ### var_index = 3   # 0=time, 1=X, 2=Y, 3=V  ###############
    ### var_name  = var_name_list[ var_index ]
    
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
def read_and_show_rtg( rtg_filename, long_name, VERBOSE=True,
                       rtg_type='FLOAT',
                       cmap='jet', BLACK_ZERO=False,
                       stretch='hist_equal',
                       a=1, b=2, p=0.5, im_file=None,
                       xsize=8, ysize=8, dpi=None ):

    #-----------------------------------------------------------    
    # Note:  Added "rtg_type" keyword on 2023-08-24 to support
    #        data types other than FLOAT, like INTEGER.
    #        rtg_type is then passed to rtg.read_grid().
    #-----------------------------------------------------------
    rtg = rtg_files.rtg_file()
    OK  = rtg.open_file( rtg_filename )
    if not(OK):
        print('Sorry, Could not open RTG file:')
        print( rtg_filename )
        return
    
    grid   = rtg.read_grid( VERBOSE=VERBOSE, rtg_type=rtg_type )
    extent = rtg.get_bounds()
    rtg.close_file()

    if (VERBOSE):
        print('Byte swap needed =', rtg.byte_swap_needed())
        print('Reading grid from RTG file...')
        print('extent =', extent)
        print('min(grid), max(grid) =', grid.min(), grid.max())
        print('Finished.')
        print()

    show_grid_as_image( grid, long_name, extent=extent, cmap=cmap,
                        BLACK_ZERO=BLACK_ZERO, stretch=stretch,
                        a=a, b=b, p=p, im_file=im_file,
                        xsize=xsize, ysize=ysize, dpi=dpi)
                              
#   read_and_show_rtg()
#--------------------------------------------------------------------
def show_grid_as_image( grid, long_name, extent=None,
                        cmap='rainbow', BLACK_ZERO=False,
                        LAND_SEA_BACKDROP=False,
                        stretch='power3', a=1, b=2, p=0.5,
                        subtitle=None, nodata=-9999.0,  ####
                        ## units='unknown units',
                        NO_SHOW=False, im_file=None,
                        xsize=4, ysize=4, dpi=None): 
                        #### xsize=8, ysize=8, dpi=None): 

    #-------------------------------------------------------------
    # Note:  Seem to have a memory issue while creating movies.
    #        Crashes with message: "zsh: killed python ..."
    #        Changed xsize & ysize from 8 to 4. (2021-11-18).
    #        If problem is with ffmpeg, may be able to use
    #        indexmem and rtbufsize flags to increase memory.
    #        Actually seems due to issue just below.
    #-------------------------------------------------------------
    # Note:  When an image is used to make a movie, then size of
    #        image in pixels MUST be a multiple of the so-called
    #        macro_block_size = 16. Image_xsize = dpi * xsize2D.
    #        Image_ysize = dpi * ysize2D.  Note:  192 = 16 * 12.
    #        If dpi=192, xsize=4 and ysize=3.3333=(10/3), then:
    #        Image_xsize = 192 * 4 = 768 = 16 * 48,
    #        Image_ysize = 192 * (10/3) = 640 = 16 * 40, 
    #        4/(10/3) = 1.2 and 768/640 = 1.2.
    #        See: show_grid_as_image().
    #-------------------------------------------------------------                        
    # Note:  xsize=4, ysize=3.3333 works pretty well for:
    #        Baro River (Eth.):  ncols = 134, nrows = 129.
    #        Tana River (Kenya): ncols = 131, nrows = 115.
    #        Awash River (Eth.): ncols = 163, nrows = 141.
    #        Must accommodate, axes labels and colorbar. 
    #-------------------------------------------------------------
    # Note:  extent = [minlon, maxlon, minlat, maxlat]
    #-------------------------------------------------------------
        
    #-------------------------
    # Other color map names
    #--------------------------------------------
    # hsv, jet, gist_rainbow (reverse rainbow),
    # gist_ncar, gist_stern
    #--------------------------------------------

    #------------------------------------   
    # Record locations of nodata values
    #---------------------------------------------------------
    # Remove nodata values so they don't affect the stretch.
    # Will later change nodata values to NaNs in grid2,
    # and they will be colored white.
    #---------------------------------------------------------
    w_nodata = (grid <= nodata)  # boolean array
    emin = grid[ grid > nodata ].min()
    grid[ w_nodata ] = emin

    #---------------------------------------    
    # This seems to be needed (2023-08-24)
    #---------------------------------------
    if ('int' in str(grid.dtype)):
        grid = np.float32(grid)
 
    #-----------------------------
    # Apply a "stretch function"
    #-----------------------------
    grid2 = stretch_grid( grid, stretch, a=a, b=b, p=p )
    if ('float' in str(grid2.dtype)):
        grid2[ w_nodata ] = np.nan  # won't work for ints

    #---------------------------------------
    # Modify the colormap (0 = black) ?
    # cmap is name of colormap, a string
    #--------------------------------------------------------
    # cmap arg to imshow can be name (as str) or cmap array
    # 4th entry is opacity, or alpha channel (I think)
    #--------------------------------------------------------
    # See: "Creating listed colormaps" section at:
    # https://matplotlib.org/3.1.0/tutorials/colors/
    #         colormap-manipulation.html
    #--------------------------------------------------------
    # "Land green" = #c6e5bc = (198, 229, 188)
    # "Sea blue"   = #aad3df = (170, 211, 223)
    #--------------------------------------------------------
    if (BLACK_ZERO):
        n_colors = 256
        color_map  = cm.get_cmap( cmap, n_colors )
        new_colors = color_map( np.linspace(0, 1, n_colors) )
        black = np.array([0.0, 0.0, 0.0, 1.0])
        new_colors[0,:] = black
        new_cmap = ListedColormap( new_colors )
    elif (LAND_SEA_BACKDROP):
        n_colors = 256
        color_map  = cm.get_cmap( cmap, n_colors )
        new_colors = color_map( np.linspace(0, 1, n_colors) )
        land_green = np.array([198, 229, 188, 256]) / 256.0
        sea_blue   = np.array([170, 211, 223, 256]) / 256.0
        new_colors[0,:]   = land_green
        new_colors[255,:] = sea_blue
        new_cmap = ListedColormap( new_colors )
    else:
        # Note that cmap keyword is a string.
        n_colors = 256
        color_map  = cm.get_cmap( cmap, n_colors )
        new_colors = color_map( np.linspace(0, 1, n_colors) )
        new_cmap = ListedColormap( new_colors )
    #------------------------------------------------------
    # Also see "set_under" and "set_over":
    # https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.colors.Colormap.html
    #------------------------------------------------------
    new_cmap.set_bad( color = 'white')  # to color Nans white

    #--------------------
    # Set the fontsizes
    #------------------------------------------------------------------
    # fontsizes = ['xx-small', 'x-small', 'small', 'medium', 'large', 
    #              'x-large', 'xx-large', 'larger', 'smaller']
    #------------------------------------------------------------------
#     fontsize1 = 'large'
#     fontsize2 = 'medium'
#     fontsize3 = 'small'
    #---------------------
#     fontsize1 = 12
#     fontsize2 = 10
#     fontsize3 = 8
    #------------
    # Use these
    #------------    
    ## fontsize1 = 10
    if (xsize >= 4.0):
        fontsize1 = 8
        fontsize2 = 8
        fontsize3 = 6
    else:
        fontsize1 = 7
        fontsize2 = 5
        fontsize3 = 5

    #-------------------------------------------------    
    # For Tana River basin, use xsize=4, ysize=3.333
    #-------------------------------------------------
    if (ysize > 3.3) and (ysize < 3.334):
        ysize = 10/3.0
            
    #----------------------------
    # Set up and show the image
    #----------------------------
    # figure = plt.figure(1, figsize=(xsize, ysize), dpi=dpi)
    fig, ax = plt.subplots( figsize=(xsize, ysize), dpi=dpi)
    ax.set_xlabel('Longitude [deg]', fontsize=fontsize2)
    ax.set_ylabel('Latitude [deg]', fontsize=fontsize2)
    ax.tick_params(axis='both', labelsize=fontsize2)

    #---------------------------------------
    # Construct image title from long_name
    #---------------------------------------
    im_title = long_name.replace('_', ' ').title()
    acronyms = ['Chirps', 'Gpm', 'Gldas', 'Dem']
    for word in acronyms:
        if (word in im_title):
            im_title = im_title.replace(word, word.upper())
    #---------------------------------
    # Add measurement units to title
    #---------------------------------
    if (im_title.endswith('Area')):
        im_title += ' [km2]'
    elif (im_title.endswith('Elevation')):
        im_title += ' [m]'
    elif (im_title.endswith('Slope')):
        im_title += ' [m/m]'
    elif (im_title.endswith('Discharge')):
        im_title += ' [m3/s]'
    elif (im_title.endswith('Depth')):
        im_title += ' [m]'
    elif (im_title.endswith('Velocity')):
        im_title += ' [m/s]'
                
    #---------------------    
    # Adjust the margins
    #--------------------------------------
    # Right margin is wider for color bar
    #--------------------------------------
    # plt.margins(x=0.25, y=0.05)  # doesn't work??
    # This works well for xsize=4, ysize=3.1
    plt.subplots_adjust(left=0.16, bottom=0.05, right=0.80, top=0.97)
    ##plt.subplots_adjust(left=0.15, bottom=0.05, right=0.85, top=0.95)
        
    #-------------------------
    # Add title and subtitle
    #-------------------------
    ## subtitle = '2021-01-10'  # for testing
    if (subtitle is not None):
        # Spacing is affected by image ysize
        # plt.suptitle(im_title, fontsize=fontsize1)
        # plt.title( subtitle, fontsize=fontsize2 )
        ## plt.suptitle(im_title, fontsize=fontsize1, y=0.925)
        #---------------------------------------
        # This method does better with margins
        #---------------------------------------
        new_title = im_title + '\n' + subtitle
        plt.title( new_title, fontsize=fontsize2 )
    else:
        plt.title( im_title, fontsize=fontsize1 )
        # plt.title( im_title )
        # plt.title( im_title, fontsize=14)
        # ax.set_title( im_title )     
    gmin = np.nanmin( grid2 )  # ignore NaNs
    gmax = np.nanmax( grid2 )
    im = ax.imshow(grid2, interpolation='nearest', cmap=new_cmap,
                   vmin=gmin, vmax=gmax, extent=extent)

    #----------------------------------------------
    # Compute tick labels (in unstretched values)
    #---------------------------------------------
    # Change NaN back to nodata for this to work
    #---------------------------------------------
    grid2[ w_nodata ] = -9999.0   ## Need this!
    ## tickLocs = np.linspace(gmin, gmax, 8)
    tickLocs = np.linspace(gmin, gmax, 10)
    cbarLabels = []
    unStretched = grid
    for val in tickLocs:
        idx = np.abs(grid2.flatten() - val).argmin()
        # print(idx, unStretched.flatten()[idx])
        if ('int' in str(grid.dtype)):
            cbarLabels.append("%d"%(unStretched.flatten()[idx]))    ####### TEST THIS!
        else:
            cbarLabels.append("%.3f"%(unStretched.flatten()[idx]))
        ## cbarLabels.append("%.2f"%(unStretched.flatten()[idx]))
   
    #----------------
    # Add a colorbar
    #----------------
    # https://matplotlib.org/stable/tutorials/colors/colorbar_only.html
    # https://stackoverflow.com/questions/32462881/add-colorbar-to-existing-axis
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(im, cax=cax, ticks=tickLocs)
    cax.set_yticklabels(cbarLabels, fontsize=fontsize3)
    #------------------------------------
    # Similar, but can't add cbarLabels
    #------------------------------------
    # cbar = fig.colorbar(im, ax=ax, ticks=tickLocs, pad=0.03,
    #                     shrink=0.8, label='units')
                        # yticklabels=cbarLabels)  # doesn't work
    ## cbar.cax.set_yticklabels(cbarLabels)  # doesn't work  
             
    #--------------------------------------------------------        
    # NOTE!  Must save before "showing" or get blank image.
    #        File format is inferred from extension.
    #        e.g. TMP_Image.png, TMP_Image.jpg.
    #--------------------------------------------------------
    if (im_file is not None):
        plt.savefig( im_file, dpi=dpi )
        ## plt.savefig( im_file )  # dpi defaults to "figure"
    else:
        plt.show()   # Ignore NO_SHOW arg for now.   
    #-----------------------------------------------
#     if (im_file is not None):  
#         plt.savefig( im_file )
#     if not(NO_SHOW):
#         plt.show()
 
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
              stretch='power3', a=1, b=2, p=0.5,
              cmap='rainbow', REPORT=True,
              BLACK_ZERO=False, LAND_SEA_BACKDROP=False,
              xsize=4, ysize=4, dpi=192 ):

    #------------------------------------------------------------
    # Note:  Seem to have a memory issue while creating movies.
    #        Crashes with message: "zsh: killed python ..."
    #        Changed xsize & ysize from 6 to 4. (2021-11-18).
    #        If problem is with ffmpeg, may be able to use
    #        indexmem and rtbufsize flags to increase memory.
    #------------------------------------------------------------
    
    # Example nc_files:
    # nc_file = case_prefix + '_2D-Q.nc'
    # nc_file = case_prefix + '_2D-d-flood.nc'

    if ('_2D' not in nc_file):
        print('ERROR: This function is only for TopoFlow "2D" files.') 
        return

    ncgs = ncgs_files.ncgs_file()        
    ncgs.open_file( nc_file )
    var_name_list = ncgs.get_var_names( no_dim_vars=True )  ####
    var_index  = 0   # (dim vars are now excluded)
    var_name   = var_name_list[ var_index ]
    long_name  = ncgs.get_var_long_name( var_name )
    n_grids    = ncgs.ncgs_unit.variables[var_name].n_grids
    datetimes  = ncgs.ncgs_unit.variables[ 'datetime' ]  ########
    #---------------------------------------
    start_year = str(datetimes[0])[:4]
    # end_year   = str(datetimes[-1])[:4]     # Is also 2005 for stats...
    end_year   = '2015'                       ##### TEMP FIX #####
    years_str  = start_year + '-' + end_year
    month_map  = {1:'Jan', 2:'Feb', 3:'Mar', 4:'Apr', 5:'May', 6:'Jun',
                  7:'Jul', 8:'Aug', 9:'Sep', 10:'Oct', 11:'Nov', 12:'Dec'}

    #-----------------------------------------------
    # MINT netCDF conventions:
    #    bounds = [minlon, minlat, maxlon, maxlat]
    # For matplotlib.pyplot.imshow():
    #    extent = [minlon, maxlon, minlat, maxlat]
    #-----------------------------------------------
    if (extent is None):
        bounds = ncgs.ncgs_unit.geospatial_bounds 
        extent = [ bounds[0], bounds[2], bounds[1], bounds[3]]

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

        #-------------------------------------            
        # Prepare subtitle w/ date for image
        #-------------------------------------
        datetime = datetimes[ time_index ]
        if ('2D-month' in nc_file):
            month_num = int(str(datetime)[5:7])
            month_str = month_map[ month_num ] + ': '
            subtitle = 'Average for ' + month_str + years_str
        elif ('days-per-month' in nc_file):
            datetime = datetime[:10]   # (just the date)
            subtitle = str(datetime)
        elif ('pop-1-or-more' in nc_file):
            datetime = datetime[:10]   # (just the date)
            subtitle = str(datetime)
        else:
            subtitle = str(datetime)
        #-----------------------------
        time_index += 1

        #----------------------------------------    
        # Build a filename for this image/frame
        #----------------------------------------
        tstr = str(time_index)
        pad = time_pad_map[ len(tstr) ]
        time_str = (pad + tstr)
        im_file = im_file_prefix + time_str + '.png' 
        im_file = (png_dir + '/' + im_file)
                
        show_grid_as_image( grid, long_name, cmap=cmap,
                            stretch=stretch, a=a, b=b, p=p, 
                            extent=extent,
                            BLACK_ZERO=BLACK_ZERO,
                            LAND_SEA_BACKDROP=LAND_SEA_BACKDROP,
                            NO_SHOW=True, im_file=im_file,
                            subtitle=subtitle,
                            xsize=xsize, ysize=ysize,
                            dpi=dpi )  ## dpi=None
                            
    ncgs.close_file()
    tstr = str(time_index)
    print('Finished saving grid images to PNG files.')
    print('   Number of files = ' + tstr)
    print()

#   save_grid_stack_as_images()
#--------------------------------------------------------------------
def save_rts_as_images( rts_file, png_dir, extent=None,
                        long_name='River Discharge',
                        stretch='power3', a=1, b=2, p=0.5,
                        cmap='rainbow', BLACK_ZERO=False,
                        REPORT=True,
                        start_datetime='2015-10-01 00:00:00',
                        time_interval_hours=6,
                        xsize=4, ysize=4, dpi=192 ):
                        ### xsize=6, ysize=6, dpi=192):

    # Example rts_files:
    # rts_file = case_prefix + '_2D-Q.rts'
    # rts_file = case_prefix + '_2D-d-flood.rts'

    if ('.rts' not in rts_file):
        print('ERROR: This function is only for RTS files.') 
        return

    rts = rts_files.rts_file()
    OK  = rts.open_file( rts_file )
    if not(OK):
        print('Could not open RTS file.')
        return
    n_grids = rts.number_of_grids()
    print('Byte swap needed =', rts.byte_swap_needed())

    if (extent is None):
        extent = rts.get_bounds()

    im_title = long_name.replace('_', ' ').title()
    im_file_prefix = 'TF_RTS_Movie_Frame_'
    time_pad_map = {1:'0000', 2:'000', 3:'00', 4:'0', 5:''}

    if (REPORT):
        print('Creating images from grid stack in rts_file:')
        print('  ' + rts_file )
        print('  ' + 'long name =', long_name)
        print('  ' + 'n_grids   =', n_grids)
        print('  ' + 'extent    =', extent)
        print('This may take a few minutes.')
        print('Working...')
        
    time_index = 0
    time_units = 'hours'
    rts_min = 1e12
    rts_max = 1e-12

    while (True):
        # print('time index =', time_index )
        try:
            grid = rts.read_grid( time_index )   # alias to get_grid()
            gmin = grid.min()
            gmax = grid.max()
            rts_min = min( rts_min, gmin )
            rts_max = max( rts_max, gmax )
        except:
            break

        #----------------------------------------    
        # Build a filename for this image/frame
        #----------------------------------------
        tstr = str(time_index + 1)
        pad = time_pad_map[ len(tstr) ]
        time_str = (pad + tstr)
        im_file = im_file_prefix + time_str + '.png' 
        im_file = (png_dir + '/' + im_file)

        #---------------------------------------------
        # Create a timestamp subtitle for this image
        #---------------------------------------------
        time = (time_index * time_interval_hours)
        current_datetime = tu.get_current_datetime( start_datetime, time, time_units )
        subtitle = str(current_datetime)
        time_index += 1

        show_grid_as_image( grid, long_name, cmap=cmap,
                            stretch=stretch, a=a, b=b, p=p,
                            BLACK_ZERO=BLACK_ZERO, extent=extent,
                            NO_SHOW=True, im_file=im_file,
                            subtitle=subtitle,
                            xsize=xsize, ysize=ysize, dpi=dpi )

    rts.close_file()
    print('min(rts), max(rts) =', rts_min, rts_max)
    tstr = str(time_index)
    print('Finished saving grid images to PNG files.')
    print('   Number of files = ' + tstr)  
    print()

#   save_rts_as_images()
#--------------------------------------------------------------------
def plot_time_series(nc_file, output_dir=None, var_index=0,
                     marker=',', REPORT=True, xsize=11, ysize=6,
                     im_file=None, start_date=None, end_date=None,
                     time_interval_hours=6,
                     start_index=0, end_index=-1):

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
    var_name_list = ncts.get_var_names( no_dim_vars=True )
    var_index     = 0  # now dim vars are excluded
    var_name      = var_name_list[ var_index ]
    lon_list      = ncts.get_var_lons()
    lat_list      = ncts.get_var_lats()

    if (REPORT):
        print( 'var_names in netCDF file =' )
        print( var_name_list )
        print( 'var longitudes =')
        print( lon_list )
        print( 'var latitudes =')
        print( lat_list )
        print()

    # (series, times) = ncts.get_series( var_name )
    # long_name = series.long_name
    # v_units   = series.units
    # values    = np.array( series )
    #-----------------------------------------
    ts_values = ncts.get_values( var_name )
    ts_times  = ncts.get_times()
    long_name = ts_values.long_name
    v_units   = ts_values.units
    t_units   = ts_times.units
#     values    = ts_values[:]
#     times     = ts_times[:]
    values    = ts_values[start_index:end_index]  # 12/05/22
    times     = ts_times[start_index:end_index]   # 12/05/22
    # values    = np.array( ts_values )
    # times     = np.array( ts_times )

    if (t_units == 'minutes'):
        # times = times / 60.0
        # t_units = 'hours'
        times = times / (60.0 * 24)
        t_units = 'days'

    # For testing
    ####################
    # print(' max value =', values.max())
    # print(' min value =', values.min())
    # print(' values[-50:-1] =', values[-50:-1])
    # print(' times[-50:-1]  =', times[-50:-1])
    # print()
    
    # ymin = values.min()
    ymin = 0.0
    ymax = values.max()

    figure = plt.figure(1, figsize=(xsize, ysize))
    # fig, ax = plt.subplots( figsize=(xsize, ysize))

    y_name = long_name.replace('_', ' ').title()

    ### start_date = '2005-01-01'  # for testing
    ### end_date   = '2015-01-01'  # for testing
    if (start_date is not None):
        #----------------------------------
        # Tick labels on x-axis are dates
        #----------------------------------
        times = np.arange(np.datetime64( start_date ),
                          np.datetime64( end_date),
                          np.timedelta64( time_interval_hours, 'h'))
        n_times = len(times)
        #-------------------------------------------------------------------
        start_date_obj = dt.datetime.strptime( start_date, '%Y-%m-%d' )
        end_date_obj   = dt.datetime.strptime( end_date,   '%Y-%m-%d' )
        values = values[:n_times]  ########
        #---------------------------------------------
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d')) 
        plt.plot( times, values, marker=marker)
        plt.gcf().autofmt_xdate()  # Rotate the dates  ########
        ## plt.xlabel( 'Date', fontsize=16 )
        plt.ylabel( y_name + ' [' + v_units + ']', fontsize=16 )
        plt.ylim( np.array(ymin, ymax) )
        plt.xticks( fontsize=14 ) ####### 12/05/22

#     if (start_date is not None):
#         #--------------------------------------------
#         # Tick labels on x-axis are dates
#         # Plot a tick every 2 months w/ MonthLocator
#         #   if time range is 3 years.
#         #---------------------------------------------
#         times = np.arange(np.datetime64( start_date ),
#                           np.datetime64( end_date),
#                           np.timedelta64( time_interval_hours, 'h'))
#         n_times = len(times)
#         #-------------------------------------------------------------------
#         start_date_obj = dt.datetime.strptime( start_date, '%Y-%m-%d' )
#         end_date_obj   = dt.datetime.strptime( end_date,   '%Y-%m-%d' )
#         n_months = tu.get_month_difference( start_date_obj, end_date_obj )
#         n_ticks  = 20
#         interval = np.int16( np.ceil( n_months / n_ticks ) )
#         ## print('interval =', interval)
#         #-------------------------------------------------------------------
#         # print('len(times)  =', n_times)
#         # print('len(values) =', len(values))
#         values = values[:n_times]  ###############
#         #---------------------------------------------
#         plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
#         plt.gca().xaxis.set_major_locator(mdates.MonthLocator( interval=interval ))
#         ## plt.gca().xaxis.set_major_locator(mdates.MonthLocator())  ## too crowded  
#         plt.plot( times, values, marker=marker)
#         plt.gcf().autofmt_xdate()  # Rotate the dates  ########
#         plt.xlabel( 'Date' )
#         plt.ylabel( y_name + ' [' + v_units + ']' )
#         plt.ylim( np.array(ymin, ymax) )
    
    #--------------------------------------------------       
    # Tick labels on x-axis are elapsed time, in days
    #--------------------------------------------------
    else:    
        plt.plot( times, values, marker=marker)
        plt.xlabel( 'Time' + ' [' + t_units + ']' )
        plt.ylabel( y_name + ' [' + v_units + ']' )
        plt.ylim( np.array(ymin, ymax) )

    #--------------------------------------------------------        
    # NOTE!  Must save before "showing" or get blank image.
    #        File format is inferred from extension.
    #        e.g. TMP_Image.png, TMP_Image.jpg.
    #--------------------------------------------------------
    if (im_file is not None):  
        plt.savefig( im_file )
    else:
        plt.show()
    plt.close()        

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
    var_name_list = ncgs.get_var_names( no_dim_vars=True )  ####
    var_index = 0   # (dim vars are now excluded)
    var_name  = var_name_list[ var_index ]
    lon_list  = ncps.get_var_lons()
    lat_list  = ncps.get_var_lats()

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
        print('ERROR: PNG directory is not set.')
        return
    
    # If ymin or ymax is set, it is used for all plots
    ALL_SAME_YMIN = (ymin is not None)
    ALL_SAME_YMAX = (ymax is not None)
         
    ncps = ncps_files.ncps_file()
    ncps.open_file( nc_file )

    var_name_list = ncgs.get_var_names( no_dim_vars=True )  ####
    var_index = 0   # (dim vars are now excluded)
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
        print('ERROR: Could not read profiles from file:')
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
    print('Finished saving profile images to PNG files.')
    print('   Number of files = ' + tstr)
    print()
    
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
        print('Creating movie from PNG files.')
        print('   Number of PNG files =', n_frames)
        ## print('This may take a few minutes.')
        print('   Working...')

    #---------------------------------------------------------------        
    # Note:  The default codec for imageio is "libx264", which
    #        is widely used and supports mp4; H.264 codec.
    #        You can request  "mpeg4" also, and it works.
    #        If you use Get Info on an MP4, the More Info section
    #        give codecs for these as:  "H.264" or "MPEG-4 Video".
    #        If I copy an MP4 to: topoflow36/movies on GitHub,
    #        I can't get them to play in Jupyter notebook or lab.
    #        See the notebook:  MP4_Video_Test.ipynb for more info.
    # https://imageio.readthedocs.io/en/stable/format_ffmpeg.html
    #----------------------------------------------------------------
    writer = imageio.get_writer( mp4_file, fps=fps )
    ## writer = imageio.get_writer( mp4_file, fps=fps, codec='mpeg4' )
        
    for im_file in im_file_list:
        writer.append_data(imageio.imread( im_file ))
    writer.close()
   
    if (REPORT): 
        print('Finished creating movie, MP4 format.')
        print('  ' + mp4_file)
        print()

#   create_movie_from_images()
#--------------------------------------------------------------------
def plot_data( x, y, y2=None, y3=None, y4=None, y5=None,
               y6=None, y7=None, y8=None,
               xmin=None, xmax=None, ymin=None, ymax=None,
               x_name='x', x_units='', marker=',', title=None,
               y_name='y', y_units='',
               x_size=8,   y_size=4, xticks=None, yticks=None):

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
    if (y2 is not None):
        plt.plot(x, y2, marker=marker)
    if (y3 is not None):
        plt.plot(x, y3, marker=marker)
    if (y4 is not None):
        plt.plot(x, y4, marker=marker)
    if (y5 is not None):
        plt.plot(x, y5, marker=marker)
    if (y6 is not None):
        plt.plot(x, y6, marker=marker)
    if (y7 is not None):
        plt.plot(x, y7, marker=marker)
    if (y8 is not None):
        plt.plot(x, y8, marker=marker)
                        
    plt.xlabel( x_name + ' [' + x_units + ']' )
    plt.ylabel( y_name + ' [' + y_units + ']' )
    if (title is not None):
        plt.title( title )

    plt.ylim( ymin, ymax )
    plt.xlim( xmin, xmax )
    plt.xticks( xticks )
    plt.yticks( yticks )
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
def get_image_dimensions( ncols, nrows, dpi=192 ):

    #-------------------------------------------------------------
    # Note:  When an image is used to make a movie, then size of
    #        image in pixels must be a multiple of the so-called
    #        macro_block_size = 16. Image_xsize = dpi * xsize2D.
    #        Image_ysize = dpi * ysize2D.  Note:  192 = 16 * 12.
    #        If dpi=192, xsize=4 and ysize=3.3333=(10/3), then:
    #        Image_xsize = 192 * 4 = 768 = 16 * 48,
    #        Image_ysize = 192 * (10/3) = 640 = 16 * 40, 
    #        4/(10/3) = 1.2 and 768/640 = 1.2.
    #        See: show_grid_as_image().
    #-------------------------------------------------------------                        
    # Note:  xsize=4, ysize=3.3333 works pretty well for:
    #        Baro River (Eth.):  ncols = 134, nrows = 129.
    #        Tana River (Kenya): ncols = 131, nrows = 115.
    #        Awash River (Eth.): ncols = 163, nrows = 141.
    #        NOT for Omo River (Eth.): nc = 99, nr = 151.
    #        Must accommodate, axes labels and colorbar. 
    #-------------------------------------------------------------
    # xsize / ysize = (ncols / nrows)    # roughly
    # ysize = xsize * (nrows / ncols)    # roughly
    # im_xsize = xsize * dpi   # (must be multiple of 16)
    # im_ysize = ysize * dpi   # (must be multiple of 16)
    # If xsize = 4.0, and dpi=192, then
    #    im_xsize = 768 = 48 * 16
    # Possible im_ysizes = n * 16, and 48/n ~ ncols/nrows
    # n = floor( float(nrows) / ncols) * 48 )
    # im_ysize = n * 16
    # ysize    = im_ysize / 192
    #--------------------------------------------------------------
    # Example: xsize = 4, n1 = 48, ncols = 131, nrows = 115.
    #          n2 = 42, im_ysize = 42 * 16 = 672, ysize = 3.5
    #--------------------------------------------------------------
    # Example: xsize = 4, n1 = 48, ncols = 134, nrows = 129.
    #          n2 = 46, im_ysize = 46 * 16 = 736, ysize = 3.83333
    #--------------------------------------------------------------
    # Example: xsize = 4, n1 = 48, ncols = 99, nrows = 151.
    #          n2 = 73, im_ysize = 73 * 16 = 1168, ysize = 6.083333
    #--------------------------------------------------------------
    # Note:  fac is to adjust for width of color bar, etc.
    #        ncols_pad is, too, and seems better.
    #-------------------------------------------------------------- 
    # ncols_pad = ncols
    # fac = 0.87  # (0.85 is too low, 0.9 is too high)
    #---------------------------------------------------
    # xsize = 4.0
    # ncols2 = (ncols + 30)
    # nrows2 = nrows
    #---------------------------------------------------
#     xsize = 3.0  # Reduce size to create movies faster
#     xpad = 150  # (depends on font sizes)
#     ypad = 120
#     ncols2 = (ncols + xpad)
#     nrows2 = (nrows + ypad)
#     ## xsize = 2.0
#     # fac = 0.9
#     fac = 1.0
#     #----------------------------------
#     im_xsize = (xsize * dpi)
#     n1 = np.int16( im_xsize / 16 )
#     n2 = np.int16( np.floor( fac * n1 * nrows2 / ncols2) )
#     im_ysize = (n2 * 16)
#     ysize = (im_ysize / dpi)

    #------------------------------------------
    # This worked well for Omo & Awash Rivers
    # (Not tested for movies, though.)
    #------------------------------------------
    # Note: 1.5 * 192 = 288 = 18 * 16
    #       2.0 * 192 = 384 = 24 * 16
    #       2.5 * 192 = 480 = 30 * 16
    #------------------------------------------    
    xsize_im = 1.5   # [inches]
    grid_xsize = (xsize_im * dpi)
    n1 = np.int16( grid_xsize / 16 )
    n2 = np.int16( np.floor(n1 * nrows / ncols) )
    grid_ysize = (n2 * 16)

    #---------------------------------------
    # Give all margins as multiples of 16.
    # Margin sizes are given in pixels.
    #-------------------------------------------
    # 5 * 16 works for left_margin if lats > 0
    # but use 6 * 16 to support lats < 0.
    #--------------------------------------------------------------------
    # We're also using these settings for margins (2022-05-09)
    # plt.subplots_adjust(left=0.16, bottom=0.05, right=0.80, top=0.97)
    #--------------------------------------------------------------------
    left_margin   = 16 * 6   # (96)
    right_margin  = 16 * 9   # (144)  # bigger for colorbar
    bottom_margin = 16 * 5   # (80)
    top_margin    = 16 * 8    # title & subtitle
    im_xsize = grid_xsize + left_margin + right_margin
    im_ysize = grid_ysize + bottom_margin + top_margin
    xsize    = (im_xsize / dpi)   # [inches]
    ysize    = (im_ysize / dpi)   # [inches]
        
    return (xsize, ysize)

#   get_image_dimensions()
#----------------------------------------------------------------------------
def create_media_files( output_dir=None, media_dir=None,
                        topo_dir=None, met_dir=None, misc_dir=None,
                        site_prefix=None, case_prefix='Test1',
                        movie_fps=10, dpi=192,
                        DEM_ncols=None, DEM_nrows=None, ####
                        xsize2D=4, ysize2D=3.333,
                        start_date=None, end_date=None,  #####
                        time_interval_hours=6,
                        STAT_MOVIES=False):          #####

    #-------------------------------------------------------------
    # Note:  When an image is used to make a movie, then size of
    #        image in pixels must be a multiple of the so-called
    #        macro_block_size = 16. Image_xsize = dpi * xsize2D.
    #        Image_ysize = dpi * ysize2D.  Note:  192 = 16 * 12.
    #        If dpi=192, xsize=4 and ysize=3.3333=(10/3), then:
    #        Image_xsize = 192 * 4 = 768 = 16 * 48,
    #        Image_ysize = 192 * (10/3) = 640 = 16 * 40, 
    #        4/(10/3) = 1.2 and 768/640 = 1.2.
    #        See: show_grid_as_image().
    #-------------------------------------------------------------                        
    # Note:  xsize=4, ysize=3.3333 works pretty well for:
    #        Baro River (Eth.):  ncols = 134, nrows = 129.
    #        Tana River (Kenya): ncols = 131, nrows = 115.
    #        Awash River (Eth.): ncols = 163, nrows = 141.
    #        NOT for Omo River (Eth.): nc = 99, nr = 151.
    #        Must accommodate, axes labels and colorbar. 
    #-------------------------------------------------------------
    if (output_dir is None):
        print('SORRY, output_dir argument is required.')
        print()
        return

    if (media_dir is None):
        print('SORRY, media_dir argument is required.')
        print()
        return
        
    if (topo_dir is None):
        print('SORRY, topo_dir argument is required.')
        print()
        return

    if (met_dir is None):
        print('SORRY, met_dir argument is required.')
        print()
        return

    if (misc_dir is None):
        print('SORRY, misc_dir argument is required.')
        print()
        return
              
    if (site_prefix is None):
        print('SORRY, site_prefix argument is required.')
        print()
        return
               
    #------------------------------------
    # Setup required output directories
    #------------------------------------
    temp_png_dir = media_dir + 'temp_png/'
    if not(os.path.exists( temp_png_dir )):
        os.mkdir( temp_png_dir )
    else:
        delete_png_files( temp_png_dir )  ##### 11/28/22
    #-----------------------------------------
    movie_dir = media_dir + 'movies/'
    if not(os.path.exists( movie_dir )):
        os.mkdir( movie_dir )
    #-----------------------------------------
    image_dir = media_dir + 'images/'
    if not(os.path.exists( image_dir )):
        os.mkdir( image_dir )

    #---------------------------------
    # Compute best xsize and ysize ?
    #---------------------------------
    if (DEM_ncols is not None) and (DEM_nrows is not None):
        (xsize2D, ysize2D) = get_image_dimensions( DEM_ncols, DEM_nrows )
              
    #---------------------------------------------
    # Create time series plot for all "0D" files
    # e.g. Discharge, Flood Depth, etc.
    # marker=',' means to use pixel as marker
    #---------------------------------------------
    os.chdir( output_dir )    
    nc0D_file_list = glob.glob('*0D*nc')
    for nc_file in nc0D_file_list:
        im_file = nc_file.replace('.nc', '.png')
        im_path = (image_dir + im_file)
        var_index = 0   # (dimension vars are now excluded)
        plot_time_series(nc_file, output_dir=output_dir,
                         var_index=var_index, marker=',',
                         REPORT=False, xsize=11, ysize=6,
                         im_file=im_path,
                         start_date=start_date, end_date=end_date,  ####
                         time_interval_hours=time_interval_hours)   ####
                         
    #-----------------------------------------    
    # Create images for several single grids
    #-----------------------------------------
    rtg_extensions = ['_DEM.rtg', '_slope.rtg', '_d8-area.rtg']
    rtg_file_list  = [site_prefix + ext for ext in rtg_extensions]
    long_name_list = ['land_surface_elevation', 'land_surface_slope',
                      'total_contributing_area']
    k = 0
    for rtg_file in rtg_file_list:
        if (rtg_file.endswith('_d8-area.rtg')):
            stretch = 'power3'
        elif (rtg_file.endswith('_slope.rtg')):
            stretch = 'tanh'
        else:
            stretch = 'hist_equal'
        #--------------------------------------------------
        im_file   = rtg_file.replace('.rtg', '.png')
        im_path   = (image_dir + im_file)
        rtg_path  = (topo_dir + rtg_file)
        long_name = long_name_list[k]
        k += 1
        read_and_show_rtg( rtg_path, long_name, cmap='jet',
                           ### BLACK_ZERO=False,
                           stretch=stretch, VERBOSE=True,
                           xsize=xsize2D, ysize=ysize2D, dpi=dpi,
                           im_file=im_path)
        print('Finished saving grid as image.')
        print('  ' + im_path )
        print()

    #-----------------------------------------    
    # Create image for GPW-v4 pop count grid
    #-----------------------------------------
    os.chdir( misc_dir ) 
    nc_file_list = glob.glob('*gpw*nc')
    for nc_file in nc_file_list:
        (grid, long_name, extent) = read_grid_from_nc_file( nc_file )
                                       ### time_index=1, REPORT=True )
        if ('gpw_v4' in nc_file):
            im_file = site_prefix + '_gpw-v4-popcount.png'
        else:
            im_file = nc_file.replace('.nc', '.png')
        
        im_path = (image_dir + im_file)
        show_grid_as_image( grid, long_name, extent=extent,
                    cmap='rainbow',
                    ## BLACK_ZERO=False,
                    ## stretch='log', a=10,
                    ## stretch='log', a=0.2,  # (best yet)
                    stretch='log', a=0.1,  # (best yet)
                    ## stretch='power2', a=1000, b=0.5,  # (not bad)
                    ## stretch='power3', a=1, b=2, p=0.5,
                    ## stretch='hist_equal', a=1, b=2, p=0.5,
                    im_file=im_path,
                    xsize=xsize2D, ysize=ysize2D, dpi=dpi)
        print('Finished saving grid as image.')
        print('  ' + im_path )
        print()

    ## return #######################################################
      
    #----------------------------------------------
    # Create set of images and movie for all "2D"
    # output files which contain grid stacks
    # e.g. *_2D-Q.nc, *_2D-d-flood.nc'
    #----------------------------------------------
    create_output_movies(output_dir=output_dir, movie_dir=movie_dir,
                         temp_png_dir=temp_png_dir,
                         movie_fps=movie_fps, dpi=dpi, 
                         xsize2D=xsize2D, ysize2D=ysize2D )

    #----------------------------------------------
    # Create set of images and movie for all "2D"
    # stat files which contain grid stacks
    #----------------------------------------------
    if (STAT_MOVIES):
        create_stat_movies(stat_dir=None, movie_dir=movie_dir,
                    temp_png_dir=temp_png_dir,
                    movie_fps=movie_fps, dpi=dpi, 
                    xsize2D=xsize2D, ysize2D=ysize2D )
                                                    
#     os.chdir( output_dir ) 
#     nc2D_file_list = glob.glob('*2D*nc')
#     for nc_file in nc2D_file_list:
#         #------------------------------------------
#         # Change the stretch for specific files ?
#         #------------------------------------------
#         # Use 'linear' for:
#         # ('days-per-month' in nc_file)
#         # nc_file.endswith('-d-flood.nc')
#         # nc_file.endswith('-u.nc')
#         #------------------------------------------
#         stretch = 'linear'
#         a=1;   b=2;   p=0.5
#         if (nc_file.endswith('-Q.nc')):
#             stretch = 'power3'
#             a=1;  b=2
#         elif (nc_file.endswith('-d.nc')):
#             stretch = 'power3'
#             a=1;  b=2
#         elif ('pop-1-or-more' in nc_file):
#             stretch = 'log'
#             a = 0.2
#   
#         #------------------------------------
#         # First, create a set of PNG images
#         #------------------------------------
#         save_grid_stack_as_images( nc_file, temp_png_dir,
#                                    ##### extent=None,  # auto-computed
#                                    stretch=stretch, a=a, b=b, p=p,
#                                    cmap='rainbow', REPORT=True,
#                                    xsize=xsize2D, ysize=ysize2D, dpi=dpi)
# 
#         #----------------------------------------------
#         # Create movie from set of images in temp_png
#         #----------------------------------------------
#         # movie_fps = "frames per second"
#         mp4_file = nc_file.replace('.nc', '.mp4')
#         mp4_path = (movie_dir + mp4_file)
#         create_movie_from_images( mp4_path, temp_png_dir,
#                                   fps=movie_fps, REPORT=True)
# 
#         #-----------------------------------
#         # Delete all PNG files in temp_png
#         #-----------------------------------
#         delete_png_files( temp_png_dir )

    #-----------------------------------------------------
    # Create set of images and movie for CHIRPS rainfall
    # (or other meteorology) grid stacks in RTS format
    #------------------------------------------------------
    # Note: The problem with doing this here is that the
    # met_dir may have many other RTS files for different
    # time periods.  Calling create_met_movies() from
    # met_base.finalize() is an option, but takes a long
    # time and will affect model runtime, etc.
    # Also, it doesn't need to be done after each run
    # since those files don't change.
    #------------------------------------------------------    
    RAIN_MOVIES = False   # Movies can be very large
    if (RAIN_MOVIES):
        create_met_movies( met_dir=met_dir, media_dir=media_dir,
                           movie_fps=10,
                           xsize=xsize2D, ysize=ysize2D, dpi=dpi)
                           ## xsize2D=4, ysize2D=4, dpi=192)

#   create_media_files()
#----------------------------------------------------------------------------
def create_output_movies(output_dir=None, movie_dir=None,
                         temp_png_dir=None,
                         movie_fps=10,
                         xsize2D=4, ysize2D=3.333, dpi=192):

    if (output_dir is None):
        print('SORRY, output_dir argument is required.')
        print()
        return
    if (movie_dir is None):
        print('SORRY, output_dir argument is required.')
        print()
        return
    if (temp_png_dir is None):
        print('SORRY, output_dir argument is required.')
        print()
        return
        
    #----------------------------------------------
    # Create set of images and movie for all "2D"
    # files which contain grid stacks
    # e.g. *_2D-Q.nc, *_2D-d-flood.nc'
    #----------------------------------------------
    os.chdir( output_dir ) 
    nc2D_file_list = glob.glob('*2D*nc')
    for nc_file in nc2D_file_list:
        #------------------------------------------
        # Change the stretch for specific files ?
        #------------------------------------------
        # Use 'linear' for:
        # ('days-per-month' in nc_file)
        # nc_file.endswith('-d-flood.nc')
        # nc_file.endswith('-u.nc')
        #------------------------------------------
        stretch = 'linear'
        a=1;   b=2;   p=0.5
        if (nc_file.endswith('-Q.nc')):
            stretch = 'power3'
            a=1;  b=2
        elif (nc_file.endswith('-d.nc')):
            stretch = 'power3'
            a=1;  b=2
        elif ('pop-1-or-more' in nc_file):
            stretch = 'log'
            a = 0.2
  
        #------------------------------------
        # First, create a set of PNG images
        #------------------------------------
        save_grid_stack_as_images( nc_file, temp_png_dir,
                                   ##### extent=None,  # auto-computed
                                   stretch=stretch, a=a, b=b, p=p,
                                   cmap='rainbow', REPORT=True,
                                   xsize=xsize2D, ysize=ysize2D, dpi=dpi)

        #----------------------------------------------
        # Create movie from set of images in temp_png
        #----------------------------------------------
        # movie_fps = "frames per second"
        mp4_file = nc_file.replace('.nc', '.mp4')
        mp4_path = (movie_dir + mp4_file)
        create_movie_from_images( mp4_path, temp_png_dir,
                                  fps=movie_fps, REPORT=True)

        #-----------------------------------
        # Delete all PNG files in temp_png
        #-----------------------------------
        delete_png_files( temp_png_dir )
        
#   create_output_movies()
#----------------------------------------------------------------------------
def create_stat_movies(stat_dir=None, movie_dir=None,
                       temp_png_dir=None,
                       movie_fps=10, dpi=192, 
                       xsize2D=4.0, ysize2D=3.333 ):

    if (stat_dir is None):
        stat_dir = '~/stats1/'
        stat_dir = os.path.expanduser( stat_dir )

    #------------------------------------------      
    # Just need to set output_dir to stat_dir
    #------------------------------------------    
    create_output_movies(output_dir=stat_dir, movie_dir=movie_dir,
                         temp_png_dir=temp_png_dir,
                         movie_fps=movie_fps, dpi=dpi,
                         xsize2D=xsize2D, ysize2D=ysize2D )
        
#   create_stat_movies() 
#----------------------------------------------------------------------------
def create_met_movies( met_dir=None, media_dir=None,
                       movie_fps=10, dpi=192,
                       DEM_ncols=None, DEM_nrows=None,
                       xsize2D=4, ysize2D=3.333,
                       start_datetime='2015-10-01 00:00:00',  ######
                       time_interval_hours=6,  #######                          
                       rts_files=None):  # list of RTS filenames in met_dir

    if (media_dir is None):
        print('SORRY, media_dir argument is required.')
        print()
        return

    if (met_dir is None):
        print('SORRY, met_dir argument is required.')
        print()
        return

    #------------------------------------
    # Setup required output directories
    #------------------------------------
    temp_png_dir = media_dir + 'temp_png/'
    if not(os.path.exists( temp_png_dir )):
        os.mkdir( temp_png_dir )
    #----------------------------------------
    movie_dir = media_dir + 'movies/'
    if not(os.path.exists( movie_dir )):
        os.mkdir( movie_dir )

    #---------------------------------
    # Compute best xsize and ysize ?
    #---------------------------------
    if (DEM_ncols is not None) and (DEM_nrows is not None):
        (xsize2D, ysize2D) = get_image_dimensions( DEM_ncols, DEM_nrows )
                              
    #-----------------------------------------------------
    # Create set of images and movie for CHIRPS rainfall
    # (or other meteorology) grid stacks in RTS format
    #-----------------------------------------------------
    print('Creating movies for rainfall grid stacks...')
    os.chdir( met_dir )
    if (rts_files is None):
        rts_files = glob.glob('*.rts')
    for rts_file in rts_files:
        #------------------------------------
        # Generate long_name from filename
        #------------------------------------
        if ('CHIRPS' in rts_file):
            long_name='CHIRPS Rainfall Rate'
            time_interval_hours = 6   ########
        ## elif ('GLDAS' in rts_file) and ('Rain' in rts_file):
        elif ('GLDAS' in rts_file):
            long_name='GLDAS Rainfall Rate'
            time_interval_hours = 3   ########
        ## elif ('GPM' in rts_file) and ('Rain' in rts_file):
        elif ('GPM' in rts_file):
            long_name='GPM Rainfall Rate'
            time_interval_hours = 0.5   ######
        else:
            # Could be air temperature, RH, etc.
            long_name='Unknown Variable'  ##########
                
        #------------------------------------
        # First, create a set of PNG images
        #-------------------------------------------------------
        # Timestamps will be added as a subtitle to each image
        # using start_datetime and time_interval_hours
        #------------------------------------------------------- 
        save_rts_as_images( rts_file, temp_png_dir, extent=None,
                 long_name=long_name,
                 start_datetime=start_datetime,
                 time_interval_hours=time_interval_hours,
                 stretch='hist_equal', a=1, b=2, p=0.5,    ################
                 ## stretch='power3', a=1, b=2, p=0.5,
                 cmap='rainbow', BLACK_ZERO=False,
                 REPORT=True, xsize=xsize2D, ysize=ysize2D, dpi=dpi )

        #----------------------------------------------
        # Create movie from set of images in temp_png
        #----------------------------------------------
        # movie_fps = "frames per second"
        mp4_file = rts_file.replace('.rts', '.mp4')
        mp4_path = (movie_dir + mp4_file)
        create_movie_from_images( mp4_path, temp_png_dir,
                                  fps=movie_fps, REPORT=True)

        #-----------------------------------
        # Delete all PNG files in temp_png
        #-----------------------------------
        delete_png_files( temp_png_dir )

#   create_met_movies()
#----------------------------------------------------------------------------
def delete_png_files( temp_png_dir ):

    png_files = glob.glob( temp_png_dir + '*.png' )
    for file in png_files:
        try:
           os.remove( file )
        except OSError as e:
            print("Error: %s : %s" % (file, e.strerror)) 

    #-------------------
    # Is this needed ?
    #-------------------
    time.sleep( 1.0 )

    print('Finished deleting PNG files in:')
    print('  ' + temp_png_dir)
    print()

#   delete_png_files()
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# FOR AN IDEA IN THE VISUALIZATION NOTEBOOK.
#----------------------------------------------------------------------------
# import glob, os, os.path, shutil
# 
# class dataset_info():
#     
#     def __init__(self, name='Baro'):  
#         self.home_dir = os.path.expanduser("~")
#         if (name.lower() == 'baro'):
#             #-------------------------------------------------
#             # Baro River, with mouth near Gambella, Ethiopia
#             #-------------------------------------------------
#             self.site_prefix = 'Baro-Gam_60sec'
#             self.case_prefix = 'Test1'
#             ## self.case_prefix = 'Test2'
#             self.test_dir    = self.home_dir + '/TF_Tests3/'           ########
#         elif (name.lower() == 'treynor'):
#             #-------------------------------------------------
#             # Treynor River, in Iowa (part of Nishnabotna R.)
#             #-------------------------------------------------
#             self.site_prefix = 'Treynor'
#             self.case_prefix = 'June_20_67'
#             # self.case_prefix = 'June_07_67'
#             self.test_dir    = self.home_dir + '/TF_Output/'
#     
#         self.basin_dir  = self.test_dir + self.site_prefix + '/'    
#         self.output_dir = self.basin_dir    ######
#         self.topo_dir   = self.basin_dir + '__topo/'
#         self.met_dir    = self.basin_dir + '__met/'
#         self.soil_dir   = self.basin_dir + '__soil/'
# 
#         print('Home directory   =', self.home_dir)
#         print('Basin directory  =', self.basin_dir)
#         print('Output directory =', self.output_dir)
# 
#         os.chdir( self.output_dir )  ##########
# 
# ds = dataset_info( name='Baro' )

#----------------------------------------------------------------------------
   