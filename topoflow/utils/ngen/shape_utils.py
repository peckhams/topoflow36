
# Copyright (c) 2023-2024, Scott D. Peckham
#
# Jan 2024. Renamed all "ngen/utils" files to end in "_utils.py"
#           instead of "_tools.py"
#           Modified to use new data_utils.py.
# Nov 2023. Added: get_basin_repo_dir, plot_shapefile.
# Jun 2023. Started from hlr_utils.py to write shape_utils.py.
#           Wrote create_tsv_from_shapefile,
#           write_tsv_header, write_tsv_line, print_attributes,
#           get_polygon_points, convert_coords, get_bounding_box.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import shape_utils as su
#  >>> su.create_tsv_from_shapefile()
#
#---------------------------------------------------------------------
#
#  open_shapefile()
#  plot_shapefile()   ## 2023-11-28
#
#  create_tsv_from_shapefile()    # main function
#  write_tsv_header()
#  write_tsv_line()
#  print_attributes()
#  get_polygon_points()
#  get_polygon_points1()  # alt. version
#  convert_coords()
#  get_bounding_box()
#  check_bounding_box()
#
#---------------------------------------------------------------------
import numpy as np
from osgeo import ogr, osr
import json, time

from topoflow.utils import projections as proj
from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import gages2_utils as g2

from matplotlib import pyplot as plt   # for plot_shapefile()
from matplotlib import patches

# from mpl_toolkits.basemap import Basemap  ## not installed

# from osgeo import gdal

#---------------------------------------------------------------------
def open_shapefile( shape_path ):

    file_unit = ogr.Open( shape_path )
    return file_unit
    
#   open_shapefile()
#---------------------------------------------------------------------
def plot_shapefile( data_dir=None, shape_file=None,
                    fig_xsize=10, fig_ysize=6,
                    marker_size=7, font_size=12,
                    SWB1=False, SWAP_XY=False,
                    TIGER_2012=True, CENSUS_2014=False,
                    CONVERT=False, SAVE_PNG=False ):

    #-----------------------------------------------------------
    # Note: This function currently plots US state boundaries
    #       and then plots points using SWB classes.
    #       A more general version will be given later.
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
        x1, y1 = get_polygon_points( feature )
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
    x_g2, y_g2, swb_g2 = g2.get_swb_points( ORIGINAL=SWB1 )
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

    #----------------------------------------
    # If swb_g2 is an ndarray, you can use:
    #----------------------------------------------------------
    # For an even faster method see:
    # https://stackoverflow.com/questions/16992713/
    # translate-every-element-in-numpy-array-according-to-key
    #----------------------------------------------------------
    swb_colors = np.vectorize(color_map.get)(swb_g2)
    print('swb_colors[0:10]', swb_colors[0:10])

    #------------------------------------    
    # If swb_g2 is a list, you can use:
    #------------------------------------
    # swb_colors = list(map(color_dict.get, swb_g2))
    
    #-------------------------------------------
    # Now add a circle for each GAGES2 station
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
    # for each in station_list:
    #     circle = plt.Circle((x, y), radius, color='r')
    #     ax.add_patch(circle)
    #     OR: plt.plot(x,y, marker='o', markersize=7, color='red')

    #---------------------------------------
    # Add a legend. Can change x_L and y_T
    # and everything else will adjust.
    #---------------------------------------
    fs = font_size   # set by keyword
    ms = marker_size
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
       
    #--------------------------
    # Show the completed plot
    #--------------------------
    run_time = (time.time() - start_time)
    plt.show()

    #----------------------------    
    # Save figure to image file
    #----------------------------
    if (SAVE_PNG):
        fig.savefig(data_dir + 'SWB_GAGES2_CONUS.png')

    #--------------------------------
    # Print info, close files, etc.
    #--------------------------------
    shape_unit = None ######## CHECK METHOD TO CLOSE
    print('Plotted', n_features, 'features.') 
    print('Run time =', run_time, '[sec]')
    print('Finished.')
    
#   plot_shapefile()
#---------------------------------------------------------------------
def create_tsv_from_shapefile( data_dir=None, new_dir=None,
                    dem_dir=None,
                    shape_file='gagesII_9322_sept30_2011.shp',                    
                    prj_file  ='gagesII_9322_sept30_2011.prj',
                    tsv_file='new_gages2_all.tsv',
                    nf_max=50, REPORT=True, insert_key=None,
                    ADD_BOUNDING_BOX=True, SWAP_XY=True,
                    filter_key=None, filter_value=None,
                    out_lon_heading='Longitude',
                    out_lat_heading='Latitude'):

    start_time = time.time()
    print('Running...')
    #--------------------------------------------
    # You will need to modify these directories
    # and get all the necessary files.
    #--------------------------------------------
    if (data_dir is None):
        data_dir = dtu.get_data_dir( 'USGS_GAGES2_all' )
    if (new_dir is None):
        data_dir = dtu.get_new_data_dir( 'USGS_GAGES2_all' )
        
    shape_path = data_dir + shape_file
    prj_path   = data_dir + prj_file
    shape_unit = ogr.Open( shape_path )
    layer      = shape_unit.GetLayer()
    n_features = np.int32(0)
    
    tsv_path   = new_dir + tsv_file 
    tsv_unit   = open( tsv_path, 'w') 
    if (ADD_BOUNDING_BOX):
        new_headings=['MINLON','MAXLON','MINLAT','MAXLAT']

    #-----------------------------------------    
    # Iterate over the features in the layer
    # GeometryType = 3 is POLYGON.
    #-----------------------------------------
    for feature in layer:
        geometry   = feature.GetGeometryRef()
        attributes = feature.items()  # This is a Python dictionary
        # n_rings    = geometry.GetGeometryCount()

        #--------------------------------
        # Write header for new TSV file
        #--------------------------------
        if (n_features == 0):
            write_tsv_header( tsv_unit, attributes, insert_key=insert_key,
                  new_headings=new_headings )            
              
        #----------------------------------------
        # Print attributes in attribute table ?
        #----------------------------------------
        if (REPORT):
            print_attributes( attributes, n_features )

        #---------------------------------- 
        # Get all points in this geometry
        #----------------------------------
        x1, y1 = get_polygon_points( feature )
#         print('x1 =', x1)
#         print('y1 =', y1)
        # break

        #------------------------------------------       
        # Convert coords to Geo. lon/lat (WGS-84)
        #------------------------------------------
        x2,y2 = convert_coords(x1,y1, inPRJfile=prj_path,
                               SWAP_XY=SWAP_XY, PRINT=False)
#         print('x2 =', x2)
#         print('y2 =', y2)
#         break        

        #------------------------------ 
        # Get geographic bounding box
        #------------------------------
        minlon, maxlon, minlat, maxlat = get_bounding_box(x2, y2)
#         print()
#         print('minlon, maxlon =', minlon, maxlon)
#         print('minlat, maxlat =', minlat, maxlat)
 
        #------------------------------------- 
        # Is outlet inside of bounding box ?
        #-------------------------------------
        IN_BOX = check_bounding_box(attributes,
                       minlon, maxlon, minlat, maxlat,
                       NOTIFY=REPORT,
                       out_lon_heading=out_lon_heading,
                       out_lat_heading=out_lat_heading)

        #--------------------------------------------              
        # Apply some test or filter to this feature
        # to decide whether to write it to TSV file
        #--------------------------------------------
        if (filter_key is None):
            ALL_GOOD = True
        else:
            value = attributes[ filter_key ]
            if (isinstance(value, str)):
                ALL_GOOD = (value.lower() == filter_value.lower())
            else:
                ALL_GOOD = (value == filter_value)
                        
        #-------------------------------------              
        # Write line of info to new TSV file
        #-------------------------------------
        if (ALL_GOOD):
            new_values = []  # empty list
            if (ADD_BOUNDING_BOX):
                new_values = [minlon,maxlon,minlat,maxlat]
            write_tsv_line( tsv_unit, attributes, insert_key=insert_key,
                            new_values=new_values )

        n_features += 1
        if (n_features == nf_max):
            break
            
    #--------------------------------
    # Print info, close files, etc.
    #--------------------------------
    print('Processed', n_features, 'features.') 
    tsv_unit.close()
    run_time = (time.time() - start_time)
    print('Run time =', run_time, '[sec]')
    print('Finished.')
         
#   create_tsv_from_shapefile()
#---------------------------------------------------------------------
def write_tsv_header( tsv_unit, attributes, insert_key=None,
                      new_headings=[] ):   # default is empty list

    #-------------------------------------------------------- 
    # Note: This is for writing a shapefile attribute table
    #       to a new TSV file and inserting new columns.
    #       Insert nothing if new_headings is empty list.
    #--------------------------------------------------------
    #       TSV works better than CSV since basin names
    #       often contain a comma.
    #--------------------------------------------------------    
    key_list = list( attributes.keys() )
    if (insert_key is None):
        insert_key = key_list[0]  # insert after 1st key

    tsv_header = ''
    delim = '\t'  # tab delimited
    
    for key in key_list:
        tsv_header += (key + delim)
        if (key == insert_key):
            for heading in new_headings:
                tsv_header += (heading + delim)
    tsv_header = tsv_header[:-1] + '\n' # remove delim, add newline
    tsv_unit.write( tsv_header )

#   write_tsv_header()          
#---------------------------------------------------------------------
def write_tsv_line( tsv_unit, attributes, insert_key=None,
                    new_values=[]):  # default is empty list

    #-------------------------------------------------------- 
    # Note: This is for writing a shapefile attribute table
    #       to a new TSV file and inserting new columns.
    #       Insert nothing if new_headings is empty list.
    #--------------------------------------------------------   
    key_list = list( attributes.keys() )
    if (insert_key is None):
        insert_key = key_list[0]  # insert after 1st key

    k = 0
    line = ''
    delim = '\t'  # tab delimited
    for val in attributes.values():
#         if (isinstance(val, str)):
#             if (',' in val):
#                 val = val.replace(',',';')  ############# BAD IDEA

        line += str(val) + delim   ###### CHECK FOR LOSS
        key = key_list[ k ]
        if (key == insert_key):
            for value in new_values:
                val_str = str(value)  # works on strings, too
                ## val_str = '{x:.4f}'.format(x=value)
                line += val_str + delim
        k += 1
    line = line[:-1] + '\n'   # remove delim, add newline
    tsv_unit.write( line )

#   write_tsv_line()          
#---------------------------------------------------------------------
def print_attributes( attributes, n_features ):

    print('##### FOR FEATURE =', n_features)
    for key in attributes.keys():
        pad = ' '*( 12 - len(key))  # assume max key len = 12
        print(key + pad + '=', attributes[key] )
    print()
         
#   print_attributes()          
#---------------------------------------------------------------------
def get_polygon_points( feature, WARNING=True ):

    #---------------------------------------------------- 
    # Get all points in this geometry by exporting
    # feature to GeoJSON and converting json to dict().
    #----------------------------------------------------
    # Printing jstr shows that triple square brackets
    # are used for coordinates, so need subscript [0].
    # This is likely because polygons can have more
    # than one "ring" for holes (in general).
    #----------------------------------------------------
    jstr = feature.ExportToJson()
    ## print(jstr)
    jdict = json.loads( jstr ) # python dictionary
    gtype = jdict['geometry']['type'].lower()
    ## print('gtype =', gtype)
 
    #------------------------------------------------------
    #  There appears to be some error in the CAMELS:
    #  shapefile: 'shapefiles/merge/basinset_gf_nhru.shp'
    #------------------------------------------------------
    if ('07067000' in jstr):
        print('## NOTE: ID = 07067000.')
        print('## gtype =', gtype)
        print()

    if (gtype == 'point'):
        if (WARNING):
            print('## WARNING: Geometry type is point not polygon.')
        point = np.array( jdict['geometry']['coordinates'] )
        x1 = point[0]
        y1 = point[1]  
    elif (gtype == 'polygon'):
        points = np.array( jdict['geometry']['coordinates'][0] )
        x1 = points[:,0]
        y1 = points[:,1]
    elif (gtype == 'multipolygon'):
        #-----------------------------------------------
        # This was needed for a CAMELS shapefile:
        # 'shapefiles/merge/basinset_gf_nhru.shp' that
        # had features with quadruple square brackets.
        #-----------------------------------------------
#         print('## NOTE: gtype is multipolygon.')
#         print('## jstr =')
#         print(jstr)
#         print()
        points = np.array( jdict['geometry']['coordinates'][0][0] )
        x1 = points[:,0]
        y1 = points[:,1]    
    else:
        if (WARNING):
            print('## WARNING: Unknown gtype =', gtype)
            print('## Skipping to next feature.')

    return x1, y1

#     # Before adding "multipolygon" support
#     if (gtype == 'point'):
#         if (WARNING):
#             print('## WARNING: Geometry type is point not polygon.')
#         point = np.array( jdict['geometry']['coordinates'] )
#         x1 = point[0]
#         y1 = point[1]  
#     else:
#         points = np.array( jdict['geometry']['coordinates'][0] )
#         x1 = points[:,0]
#         y1 = points[:,1]
#     return x1, y1
        
#   get_polygon_points()
#---------------------------------------------------------------------
def get_polygon_points1( geometry ):

    #--------------------------------------------------
    # Note: About same speed as get_polygon_points().
    #--------------------------------------------------
    ring = geometry.GetGeometryRef(0)  # outer ring
    p  = ring.GetPoints()    # This is a list of xy tuples
    p2 = np.array(p)
    x1 = p2[:,0]
    y1 = p2[:,1]
    return x1, y1

#   get_polygon_points1()
#---------------------------------------------------------------------
def convert_coords(x1, y1, inEPSG=None, inPRJfile=None,
                   outEPSG=4326, outPRJfile=None,
                   PRINT=False, SWAP_XY=True):

    #--------------------------------------------------------
    # Note: EPSG=4326 is for WGS-84 (Geographic lon/lat)
    #--------------------------------------------------------
    # Read WKT (Well Known Text) from shapefile's PRJ file.
    #--------------------------------------------------------
    # Note: Must set SWAP_XY to True for most datasets.
    #--------------------------------------------------------

    #------------------------------------------------    
    # Get Spatial Reference System for input coords
    #------------------------------------------------
    srs_in = osr.SpatialReference() 
    if (inPRJfile is not None):
        prj_unit = open(inPRJfile, 'r')
        prj_wkt  = prj_unit.read()
        prj_unit.close()
        srs_in.ImportFromWkt( prj_wkt )
    elif (inEPSG is not None):
        srs_in.ImportFromEPSG( inEPSG )
    else:
        print('### SORRY: An input SRS must be specified')
        print('### using the inEPSG or inPRJfile keyword.')
        print()
        return

    #-------------------------------------------------    
    # Get Spatial Reference System for output coords
    #-------------------------------------------------  
    srs_out = osr.SpatialReference()
    if (outPRJfile is not None):
        prj_unit = open(outPRJfile, 'r')
        prj_wkt  = prj_unit.read()
        prj_unit.close()
        srs_out.ImportFromWkt( prj_wkt )
    elif (outEPSG is not None):
        srs_out.ImportFromEPSG( outEPSG )
    else:
        print('### SORRY: An output SRS must be specified')
        print('### using the outEPSG or outPRJfile keyword.')
        print()
        return
        
    #-----------------------------------
    # Create coordinate transformation
    #-----------------------------------
    coordTransform = osr.CoordinateTransformation(srs_in, srs_out)

    #------------------------------------------
    # Is (x1,y1) a point or a basin polygon ?
    #------------------------------------------
    if (x1.size == 1):
        #----------------------
        # Transform the point
        #----------------------
        point = ogr.Geometry(ogr.wkbPoint)  # WKB = Well Known Binary
        point.AddPoint(x1, y1)
        point.Transform( coordTransform )
        if (SWAP_XY):
            x2 = point.GetY()
            y2 = point.GetX()
        else:
            x2 = point.GetX()
            y2 = point.GetY()        
    else:
        #----------------------------------
        # Transform each point in polygon
        #----------------------------------
        x2 = np.zeros(x1.size, dtype='float64')
        y2 = np.zeros(x1.size, dtype='float64')
        for k in range(x1.size): 
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(x1[k], y1[k])
            point.Transform( coordTransform )
            if (SWAP_XY):
                x2[k] = point.GetY()
                y2[k] = point.GetX()
            else:
                x2[k] = point.GetX()
                y2[k] = point.GetY()            
            
        #----------------------------------
        # Transform each point in polygon
        #----------------------------------
#         ring = ogr.Geometry(ogr.wkbLinearRing)
#         for k in range(x1.size): 
#             point = ogr.Geometry(ogr.wkbPoint)
#             point.AddPoint(x1[k], y1[k])
#             point.Transform( coordTransform )
#             x2k = point.GetY()
#             y2k = point.GetX()  #### Need to swap.
#             x2[k] = x2k
#             y2[k] = y2k    
#             ring.AddPoint(x2k, y2k)
#         p  = ring.GetPoints()  # list of tuples
#         p2 = np.array(p)
#         x2 = p2[:,0]
#         y2 = p2[:,1]
        # x2 = ring.GetX()  # only gets x2[0]
        #------------------------------------------   
#         basin = ogr.Geometry(ogr.wkbPolygon)
#         basin.AddGeometry(ring)
#         basin.CloseRings()
#         p  = basin.GetPoints() # DOESN'T WORK
#       # This doesn't work
#       basin.Transform( coordTransform )
#       x2 = basin.GetX()
#       y2 = basin.GetY()
              
    #----------------------------------
    # Print input and output points ?
    #----------------------------------
    if (PRINT):
        print('Input:  x1, y1 =', x1, ',', y1 )
        print('Output: x2, y2 =', x2, ',', y2 )  # lon,lat
        print()

    return x2, y2

#   convert_coords()
#---------------------------------------------------------------------
def get_bounding_box(x2, y2, buffer=0.02, decimals=5):

    #---------------------------------------------------------
    # Note: Buffer is in decimal degrees, and 0.02 is
    #       roughly 2.2 km.  See Wikipedia: Decimal degrees.
    #       Recall that HYDRO1K DEM has 1 km grid cells.
    #---------------------------------------------------------
    if (not(isinstance(x2, float))):
        minlon = x2.min() - buffer
        maxlon = x2.max() + buffer
        minlat = y2.min() - buffer
        maxlat = y2.max() + buffer
    else:
        minlon = x2 - buffer
        maxlon = x2 + buffer
        minlat = y2 - buffer
        maxlat = y2 + buffer
            
    #--------------------------------------- 
    # Accuracy is roughly 11 meters if we
    # round to 4 places after decimal.
    #---------------------------------------
    # Accuracy is roughly 1.1 meters if we
    # round to 5 places after decimal.
    #------------------------------------------- 
    # Can only do in-place w/out for ndarrays.
    #-------------------------------------------
    minlon = np.around(minlon, decimals=decimals)
    maxlon = np.around(maxlon, decimals=decimals)
    minlat = np.around(minlat, decimals=decimals)
    maxlat = np.around(maxlat, decimals=decimals)
 
    return minlon, maxlon, minlat, maxlat
    
#   get_bounding_box()
#---------------------------------------------------------------------
def check_bounding_box( atts, minlon, maxlon, minlat, maxlat,
          out_lon_heading='Longitude', out_lat_heading='Latitude',
          NOTIFY=True):

    #-------------------------------------------------
    # Is the outlet lon/lat inside the bounding box?
    #-------------------------------------------------
    # Note: For GAGES-II shapefile, we have:
    #    out_lon_heading = 'LNG_GAGE'
    #    out_lat_heading = 'LAT_GAGE'
    #-------------------------------------------------
    # Note: For MOPEX shapefile, we have:
    #    out_lon_heading = 'Longitude'
    #    out_lat_heading = 'Latitude'
    #-------------------------------------------------
    # Note: For HLR shapefile, we have:
    #    out_lon_heading = 
    #    out_lat_heading = 
    #-------------------------------------------------     
    IN_RANGE = True
    outlet_lon = atts[ out_lon_heading ]
    outlet_lat = atts[ out_lat_heading ]                  
    if ((outlet_lon < minlon) or (outlet_lon > maxlon) ):
        IN_RANGE = False          
    if ((outlet_lat < minlat) or (outlet_lat > maxlat) ):
        IN_RANGE = False
        
    if (NOTIFY):
        if (not(IN_RANGE)):
            print('### WARNING: Outlet is not in bounding box.')
            print('Lon, minlon, maxlon =', outlet_lon, minlon, maxlon)
            print('Lat, minlat, maxlat =', outlet_lat, minlat, maxlat)
            print()
    
    return IN_RANGE

#   check_bounding_box()
#---------------------------------------------------------------------


