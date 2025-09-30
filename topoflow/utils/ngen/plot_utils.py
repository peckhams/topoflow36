
# Copyright (c) 2023-2024, Scott D. Peckham
#
# Aug 2024. Moved create_swb_scatter_plot() from shape_utils.py
#           to here.
# Nov 2023. Wrote original version of create_swb_scatter_plot().
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import plo_utils as pu
#  >>> pu.create_swb_scatter_plot()
#  >>> pu.create_swb_scatter_plot( SWB1=True )
#
#---------------------------------------------------------------------
#
#  create_swb_scatter_plot()   # 2023-11-28 (original)
#
#---------------------------------------------------------------------
import numpy as np
from osgeo import ogr  ##, osr
import time

from topoflow.utils import projections as proj
from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import shape_utils as su
from topoflow.utils.ngen import gages2_utils as g2

from matplotlib import pyplot as plt   # for plot_shapefile()
from matplotlib import patches

# from mpl_toolkits.basemap import Basemap  ## not installed
# from osgeo import gdal

#---------------------------------------------------------------------
def create_swb_scatter_plot( data_dir=None, source='GAGES2_ALL',
                    ## or can set source='GAGES2_SB3',
                    shape_file=None,
                    fig_xsize=10, fig_ysize=6,
                    marker_size=3, font_size=12,
                    SWB1=False, SWAP_XY=False,
                    TIGER_2012=True, CENSUS_2014=False,
                    CONVERT=False,
                    SAVE_PNG=False, save_dir=None):

    #-----------------------------------------------------------
    # Note: This function currently plots US state boundaries
    #       and then plots points using SWB classes.
    #-----------------------------------------------------------
    # Note: This function uses only Python packages that are
    #       already used by TopoFlow and that are available
    #       in its "tf36" conda environment.  Gdal is used to
    #       read the shapefile and matplotlib is used to plot
    #       the shape boundaries (basins here).
    #       Other packages such as shapely, geopandas, and
    #       descartes could be used but are unnecessary and
    #       add extra dependencies.
    #-----------------------------------------------------------
    # Set SWB1 to True to use original SWB class ranges,
    #    which results in many unclassified basins.
    #    These will plot as white dots.
    #-----------------------------------------------------------
    # marker size for matplotlib.pyplot.plot is in points.
    # dot size (s) for matplotlib.pyplot.scatter is points^2.
    #-----------------------------------------------------------
    # If you increase fig xsize and ysize, you should also
    # increase marker_size and font_size.
    #-----------------------------------------------------------        
    start_time = time.time()   
    print('Running...')
    #----------------------------------------------
    if (shape_file is None):
        TIGER_2012 = True
    #----------------------------------------------  
    if (TIGER_2012):
        data_dir = dtu.get_data_dir( 'TIGER_2012' )
        shape_file = 'tl_2012_us_state.shp'
        #------------------------------------------
        # Use these bounds to plot x1,y1 directly
        # Mercator projection is EPSG:3857  ##### CHECK
        #------------------------------------------
        if not(CONVERT):
            xmin = -14500000.0  # Mercator bounding box
            xmax = -7000000.0
            ymin = 2500000.0
            ymax = 6500000.0
            #------------------------------------------
            # xmin = -15165054.230  # Mercator bounding box for USA
            # xmax = -6778397.8760
            # ymin = 2294608.0410
            # ymax = 6947666.110    
    #----------------------------------------------
    if (CENSUS_2014):
        data_dir = dtu.get_data_dir( 'CENSUS_2014' )
        shape_file = 'usa-states-census-2014.shp'
        if not(CONVERT):
            xmin = -126.0  # -124.736342
            xmax = -66.0   # -66.945392
            ymin = 24.0    # 24.521208
            ymax = 51.0    # 49.382808
    #----------------------------------------------                 
    prj_file   = shape_file.replace('.shp', '.prj')  
    shape_path = data_dir + shape_file
    prj_path   = data_dir + prj_file
    shape_unit = ogr.Open( shape_path )
    ## shape_unit = open_shapefile( shape_path )
    layer      = shape_unit.GetLayer()
    n_features = np.int32(0)
   
    #-------------------------
    # Create a pyplot figure
    #-------------------------
    fig_config = {'figsize': (fig_xsize, fig_ysize)}
    fig = plt.figure(**fig_config)
    ax = plt.gca()
    ax.set_aspect('equal')  ########

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
       
    #------------------------------------------------
    # Create a base plot of the US state boundaries
    #------------------------------------------------ 
    # Iterate over the features in the layer
    # GeometryType = 3 is POLYGON.
    #-----------------------------------------
    for feature in layer:
        geometry = feature.GetGeometryRef()
        n_features += 1
        ### n_rings = geometry.GetGeometryCount()
 
        #-------------------------------------               
        # Don't need attributes for plotting
        #-------------------------------------
        # attributes = feature.items()  # This is a Python dictionary

        #---------------------------------- 
        # Get all points in this geometry
        #----------------------------------
        x1, y1 = su.get_polygon_points( feature )
#         print('x1 =', x1)
#         print('y1 =', y1)
        # break

        #---------------------------------------------------       
        # Convert coords to Lambert Azimuthal Equal Area ?
        #---------------------------------------------------
        # EPSG = 2163 is a Lambert Azimuthal Equal Area
        # projection w/ almost identical info to our PRJ
        # https://spatialreference.org/ref/epsg/2163/
        # https://epsg.io/2163
        # Download PRJ file and compare.
        # Also very similar to EPSG = 9311, which is a
        #  "non-deprecated replacement" for EPSG 2163?
        #-------------------------------------------------
#         if (CONVERT): 
#             x2,y2 = convert_coords(x1,y1, inPRJfile=prj_path,
#                             outEPSG=2163, SWAP_XY=SWAP_XY,
#                             PRINT=False)
                                           
        #------------------------------------------       
        # Convert coords to Geo. lon/lat (WGS-84)
        #------------------------------------------
        if (CONVERT):
            x2,y2 = convert_coords(x1,y1, inPRJfile=prj_path,
                                   SWAP_XY=SWAP_XY, PRINT=False)
#         print('x2 =', x2)
#         print('y2 =', y2)
#         break       

        #-----------------------------------    
        # Plot a shapefile (basin) polygon
        #------------------------------------------
        # A random color is assigned, by default.
        #------------------------------------------
        if (CONVERT):
            plt.plot(x2, y2, clip_on=True, color='black')
        else:       
            plt.plot(x1, y1, clip_on=True, color='black')
        #------------------------------------------        
        ## plt.plot(x2, y2, clip_on=False)

    #-----------------------------------------------
    # Read LONG_CENT, LAT_CENT, SWB2_CLASS for all
    # GAGES-II Selected Basins v3.
    #--------------------------------------------------
    # See: https://matplotlib.org/stable/api/_as_gen/
    #              matplotlib.pyplot.plot.html
    #--------------------------------------------------        
    x_g2, y_g2, swb_g2 = g2.get_swb_points( ORIGINAL=SWB1, source=source )
    x_g2 = np.float64(x_g2)
    y_g2 = np.float64(y_g2)  # Needed for convert_coords?   #################
    print('min(x_g2), max(x_g2) =', min(x_g2), ',', max(x_g2))
    print('min(y_g2), max(y_g2) =', min(y_g2), ',', max(y_g2))
        
    #--------------------------------------------------
    # Convert coordinates from Geographic to Mercator
    #--------------------------------------------------
    # EPSG = 4326 is WGS 84 (Geographic coords)
    # EPSG = 3395 is WGS 84 / World Mercator
    #--------------------------------------------------
    x_g2m, y_g2m = proj.Mercator_XY( x_g2, y_g2 )    
#     x_g2m, y_g2m = convert_coords(x_g2,y_g2, inEPSG=4326,
#                                   ### inPRJfile=prj_path,
#                                   ### outPRJfile=prj_path,
#                                   outEPSG=3395,
#                                   SWAP_XY=True, PRINT=False)
    print('min(x_g2m), max(x_g2m) =', min(x_g2m), ',', max(x_g2m))
    print('min(y_g2m), max(y_g2m) =', min(y_g2m), ',', max(y_g2m))
    
    #--------------------------------------
    # Assign colors to the 10 SWB classes
    #--------------------------------------
    color_map = {
    'A1':'red',     'A2':'orange', 'A3':'yellow',
    'B1':'lime',    'B2':'green',  'B3':'olive',  # B3 doesn't occur?
    ## 'C1':'purple',  'C2':'magenta',
    'C1':'pink',    'C2':'magenta',
    'D1':'cyan',    'D2':'cornflowerblue', 'D3':'blue',
    ##'D1':'cyan',    'D2':'blue',   'D3':'darkblue',
    'None':'white'}  # "invisible" on a white background

    #-------------------------------------------------------
    # If swb_g2 is an ndarray, you can use np.vectorize()
    # to create a "callable" that can be used to apply the
    # color_map to every value in the array swb_g2.
    #-------------------------------------------------------
    # For an even faster method see:
    # https://stackoverflow.com/questions/16992713/
    # translate-every-element-in-numpy-array-according-to-key
    #----------------------------------------------------------
    swb_colors = np.vectorize(color_map.get)(swb_g2)

    #------------------------------------    
    # If swb_g2 is a list, you can use:
    #------------------------------------
    # swb_colors = list(map(color_dict.get, swb_g2))
    
    #-------------------------------------------
    # Now add a circle for each GAGES2 site.
    # Convert xy to another CRS, if necessary.
    #-----------------------------------------------------------
    # For matplotlib.pyplot.plot, markersize units are points.
    # For matplotlib.pyplot.scatter, s units are points^2.    
    # 1 point = 1/72 inch
    #-----------------------------------------------------------
    dot_size = marker_size**2  # See Notes above.
    ax.scatter(x_g2m, y_g2m, s=dot_size, color=swb_colors, marker='o',
               linewidths=0.5, edgecolors='black')
               ## linewidths=1.5, edgecolors='face')  # defaults

    #-----------------------------
    # TEST: Plot a single circle
    #-----------------------------
#     if (CONVERT):
#         x0 =  -106.0
#         y0 =  39.0
#     else:
#         x0 = -11800000.0  # Mercator bounding box
#         y0 = 4700000.0
#     plt.plot(x0,y0, marker='o', markersize=5, color='red')
    #------------------------------------------------------------------
    # for each in site_list:
    #     circle = plt.Circle((x, y), radius, color='r')
    #     ax.add_patch(circle)
    #     OR: plt.plot(x,y, marker='o', markersize=7, color='red')

    #---------------------------------------
    # Add a legend. Can change x_L and y_T
    # and everything else will adjust.
    #---------------------------------------
    fs = font_size   # set by keyword
    ## ms = marker_size + 3
    ms  = 8
    mkr = 'o'  # circle
    # mkr = 's'  # square
    mec = 'black'  # mec = markeredgecolor for plot command
    mew = 0.5      # mew = markeredgewidth for plot command
    dx_text = 150000.0
    dx = 500000.0
    dy = 180000.0
    #---------------------
#     x_L  = -14000000.0  # left bottom
#     y_T  = 3500000.0
    #---------------------
    x_L  = -8500000.0   # right bottom
    y_T  = 3700000.0
    #----------------------
    x_L2 = x_L + dx_text
    x_M  = x_L + dx
    x_M2 = x_M + dx_text
    #--------------------------------------------------------------------------
    plt.plot(x_L,  y_T,    marker=mkr, markersize=ms, color=color_map['A1'],
             mec=mec, mew=mew)
    plt.text(x_L2, y_T,    'A1', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x_L,  y_T-dy, marker=mkr, markersize=ms, color=color_map['A2'],
             mec=mec, mew=mew)           
    plt.text(x_L2, y_T-dy, 'A2', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x_L,  y_T-2*dy, marker=mkr, markersize=ms, color=color_map['A3'],
             mec=mec, mew=mew)           
    plt.text(x_L2, y_T-2*dy, 'A3', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x_L,  y_T-3*dy, marker=mkr, markersize=ms, color=color_map['B1'],
             mec=mec, mew=mew)           
    plt.text(x_L2, y_T-3*dy, 'B1', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x_L,  y_T-4*dy, marker=mkr, markersize=ms, color=color_map['B2'],
             mec=mec, mew=mew)           
    plt.text(x_L2, y_T-4*dy, 'B2', fontsize=fs, ha='left', va='center')
    #-----------------------
    # Right half of legend
    #-----------------------
    plt.plot(x_M,  y_T,    marker=mkr, markersize=ms, color=color_map['C1'],
             mec=mec, mew=mew)
    plt.text(x_M2, y_T,    'C1', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x_M,  y_T-dy, marker=mkr, markersize=ms, color=color_map['C2'],
             mec=mec, mew=mew)           
    plt.text(x_M2, y_T-dy, 'C2', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x_M,  y_T-2*dy, marker=mkr, markersize=ms, color=color_map['D1'],
             mec=mec, mew=mew)           
    plt.text(x_M2, y_T-2*dy, 'D1', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x_M,  y_T-3*dy, marker=mkr, markersize=ms, color=color_map['D2'],
             mec=mec, mew=mew)           
    plt.text(x_M2, y_T-3*dy, 'D2', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x_M,  y_T-4*dy, marker=mkr, markersize=ms, color=color_map['D3'],
             mec=mec, mew=mew)           
    plt.text(x_M2, y_T-4*dy, 'D3', fontsize=fs, ha='left', va='center')

    #------------------------------------- 
    # Draw a rectangle around the legend
    #-------------------------------------
    width  = 1100000.0
    height = 1050000.0
    x0     = (x_L - 120000.0)
    y0     = y_T - 5*dy 
    rect = patches.Rectangle((x0,y0), width, height,
                   linewidth=1, edgecolor='black', facecolor='none')
    ax.add_patch(rect)

    #----------------------------    
    # Save figure to PNG file ?
    #----------------------------
    if (SAVE_PNG):
        if (source == 'GAGES2_SB3'):
            png_file = 'SWB_GAGES2_SB3_CONUS.png'
            if (save_dir is None):
                save_dir = dtu.get_new_data_dir('USGS_GAGES2_SB3')
        else:
            png_file = 'SWB_GAGES2_ALL_CONUS.png'
            if (save_dir is None):
                save_dir = dtu.get_new_data_dir('USGS_GAGES2_all')
        fig.savefig(save_dir + png_file)

    #--------------------------------
    # Print info, close files, etc.
    #--------------------------------
    shape_unit = None ######## CHECK METHOD TO CLOSE
    run_time   = (time.time() - start_time)
    print('Plotted', n_features, 'features in US state layer.') 
    print('Run time =', run_time, '[sec]')
    print('Finished.')
                   
    #--------------------------
    # Show the completed plot
    #--------------------------
    plt.show()

#   create_swb_scatter_plot()
#---------------------------------------------------------------------


