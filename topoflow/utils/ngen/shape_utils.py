
# Copyright (c) 2023, Scott D. Peckham
#
# Jun 2023. Started from hlr_tools.py to write shape_utils.py.
#           Wrote create_tsv_from_shapefile,
#           write_tsv_header, write_tsv_line, print_attributes,
#           get_polygon_points, convert_coords, get_bounding_box.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import gages2_tools as g2
#  >>> g2.create_tsv_from_shapefile()
#
#---------------------------------------------------------------------
#
#  open_shapefile()
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

# from osgeo import gdal

#---------------------------------------------------------------------
def open_shapefile( shape_path ):

    file_unit = ogr.Open( shape_path )
    return file_unit
    
#   open_shapefile()
#---------------------------------------------------------------------
def create_tsv_from_shapefile( data_dir=None, dem_dir=None,
                    shape_file='gagesII_9322_sept30_2011.shp',                    
                    prj_file  ='gagesII_9322_sept30_2011.prj',
#                     shape_file='bas_ref_all.shp',                    
#                     prj_file  ='bas_ref_all.prj',
                    tsv_file='new_gages2_all.tsv',
#                     tsv_file='new_gages2_conus.tsv',
                    nf_max=50, REPORT=True, insert_key=None,
                    ADD_BOUNDING_BOX=True, SWAP_XY=True,
                    filter_key=None, filter_value=None ):

    start_time = time.time()
    print('Running...')
    #--------------------------------------------
    # You will need to modify these directories
    # and get all the necessary files.
    #--------------------------------------------
    if (data_dir is None):
        data_dir  = '/Users/peckhams/Dropbox/NOAA_NextGen/'
        data_dir += '__NextGen_Example_Basin_Repo/'
        data_dir += 'GAGES-II/gagesII_9322_point_shapefile/'
        # data_dir += 'GAGES-II/boundaries-shapefiles-by-aggeco/'
        
    shape_path = data_dir + shape_file
    prj_path   = data_dir + prj_file
    shape_unit = ogr.Open( shape_path )
    layer      = shape_unit.GetLayer()
    n_features = np.int32(0)
    
    tsv_path   = data_dir + tsv_file 
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
                       NOTIFY=REPORT)
                   ## out_lon_heading='Longitude',
                   ## out_lat_heading='Latitude')

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
                line += str(value) + delim
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
                   outEPSG=4326, PRINT=False, SWAP_XY=True):

    #--------------------------------------------------------
    # Note: ESPG=4326 is for WGS-84 (Geographic lon/lat)
    #--------------------------------------------------------
    # Read WKT (Well Known Text) from shapefile's PRJ file.
    #--------------------------------------------------------
    # Note: Must set SWAP_XY to True for most datasets.
    #--------------------------------------------------------    
    if (inPRJfile is not None):
        prj_unit = open(inPRJfile, 'r')
        prj_wkt  = prj_unit.read()
        prj_unit.close()

    #-----------------------------------
    # Create coordinate transformation
    #-----------------------------------
    srs_in = osr.SpatialReference()
    if (inEPSG is None):
        srs_in.ImportFromWkt( prj_wkt )
    else:
        srs_in.ImportFromEPSG( inESPG )
    
    srs_out = osr.SpatialReference()
    srs_out.ImportFromEPSG( outEPSG )

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


