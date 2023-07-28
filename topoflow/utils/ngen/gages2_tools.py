
# Copyright (c) 2023, Scott D. Peckham
#
# Jun 2023. Started from hlr_tools.py and wrote a more general
#              set of shapefile functions in shape_utils.py.
#           Then wrote create_tsv0, create_tsv_from_shapefile,
#              get_poly_feature_v0, get_poly_feature,
#              get_gage_ID, check_gage_ID, check_overshoot,
#              check_bounding_box, check_basin_poly_file.
# Jul 2023. Added get_basin_repo_dir, get_gages2_data_dir.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import gages2_tools as g2
#  >>> g2.create_tsv_from_shapefile( nf_max=100, tsv_file="TEST.tsv" )
#  >>> g2.create_tsv_from_shapefile( nf_max=100, tsv_file="TEST2.tsv",
#         basin_option='ALL')
#  >>> g2.create_tsv_from_shapefile( nf_max=10000 )
#
#---------------------------------------------------------------------
#
#  create_tsv0()
#      
#  create_tsv_from_shapefile()   # main function
#  get_poly_feature_v0()
#  get_poly_feature()
#  get_gage_ID()
#  check_gage_ID()
#  check_overshoot()
#  check_bounding_box()
#  check_basin_poly_file()
#
#---------------------------------------------------------------------

import numpy as np
from osgeo import ogr, osr
import json, sys, time

from topoflow.utils import shape_utils as su

#---------------------------------------------------------------------
def get_basin_repo_dir():

    #-----------------------------------
    # Modify this directory as needed.
    #-----------------------------------
    repo_dir  = '/Users/peckhams/Dropbox/NOAA_NextGen/'
    repo_dir += '__NextGen_Example_Basin_Repo/'
    return repo_dir

#   get_basin_repo_dir()
#---------------------------------------------------------------------
def get_gages2_data_dir( gtype='root' ):

    #-----------------------------------
    # Modify this directory as needed.
    #-----------------------------------
    repo_dir  = get_basin_repo_dir()
    data_dir  = repo_dir + 'GAGES-II/'
    if (gtype == 'root'):
        pass
    elif (gtype == 'points'):
       data_dir += 'gagesII_9322_point_shapefile/'
    elif (gtype == 'polygons'):
       data_dir += 'boundaries-shapefiles-by-aggeco/'
       
    return data_dir
       
#   get_gages2_data_dir()
#---------------------------------------------------------------------
def create_tsv0( data_dir=None,
                 shape_file='gagesII_9322_sept30_2011.shp',                    
                 prj_file  ='gagesII_9322_sept30_2011.prj',
#                shape_file='bas_ref_all.shp',                    
#                prj_file  ='bas_ref_all.prj',
                 tsv_file='new_gages2_all.tsv',
#                 tsv_file='new_gages2_conus.csv',
                 nf_max=50, REPORT=True, basin_option='ALL'):

    #--------------------------------------------------------
    # Note:  Attributes in shapefile attribute table for
    #        gagesII_9322_sept30_2011.shp are:
    # 
    # g2_staid      = attributes['STAID']    # station ID  
    # g2_staname    = attributes['STANAME']  # station name
    # g2_class      = attributes['CLASS']    # Ref or Non-ref
    # g2_aggecoregi = attributes['AGGECOREGI'] 
    # g2_drain_sqkm = attributes['DRAIN_SQKM'] 
    # g2_huc02      = attributes['HUC02'] 
    # g2_lat_gage   = attributes['LAT_GAGE'] 
    # g2_lng_gage   = attributes['LNG_GAGE'] 
    # g2_state      = attributes['STATE'] 
    # g2_hcdn_2009  = attributes['HCDN_2009']
    # g2_active09   = attributes['ACTIVE09'] 
    # g2_flyrs1900  = attributes['FLYRS1900']                        
    # g2_flyrs1950  = attributes['FLYRS1950'] 
    # g2_flyrs1990  = attributes['FLYRS1990'] 
    #--------------------------------------------------------
    if (data_dir is None):
        data_dir = get_gages2_data_dir()

    insert_key = 'LNG_GAGE'
    if (basin_option == 'REF'):
        filter_key   = 'CLASS'
        filter_value = 'Ref'
    elif (basin_option == 'NONREF'):
        filter_key   = 'CLASS'
        filter_value = 'Non-ref'
    elif (basin_option == 'ALL'):
        filter_key   = None
        filter_value = None

    su.create_tsv_from_shapefile(data_dir=data_dir, shape_file=shape_file,
                      prj_file=prj_file, nf_max=nf_max,
                      insert_key=insert_key,
                      ADD_BOUNDING_BOX=True,
                      filter_key=filter_key, filter_value=filter_value)

#   create_tsv0()
#---------------------------------------------------------------------
def create_tsv_from_shapefile( data_dir=None, tsv_file=None,
        point_dir=None, polygon_dir=None,
        point_shp_file='gagesII_9322_sept30_2011.shp',
        nf_max=50, SWAP_XY=True,  # Must swap for GAGES-II
        REPORT=False, DEBUG=False,
        insert_key=None, ADD_BOUNDING_BOX=True, 
        filter_key=None, filter_value=None,
        basin_option='REF' ):

    #------------------------------------------------------------------
    # Note:  This routine creates a TSV file with attributes from the
    #        attribute table of point_shp_file and new attributes
    #        computed from a shapefile of basin polygons.
    #        New attributes:  minlon, maxlon, minlat, maxlat.
    #        Can also get AREA & PERIMETER from polygon shapefile.
    #        By default, it does this only for the "Reference" basins.
    #        There are 9322 total basins and 2057 reference basins.
    #        Could remove "CLASS" from TSV, since always "Ref".
    #        Could remove US state after comma in STANAME, since
    #        there is a separate STATE field.
    #------------------------------------------------------------------
    # Note:  This routine will not work unless the features in
    #        point_shp_file are sorted by station_ID and the features
    #        in reference basin polygon file are sorted by gage_ID.
    #        It turned out that this was the case for point_shp_file
    #        but needed to sort features in the ref basin poly file
    #        by gage_ID using QGIS as follows:

    #        To sort features in a shapefile by an attribute in QGIS:
    #        Select the layer: bas_ref_all.shp, in Layers panel.
    #        Choose Processing > Toolbox > Vector General >
    #           Order by Expression.
    #        Choose GAGE_ID from Expression droplist.
    #        Enter name of new shapefile, then click Run.   
    #------------------------------------------------------------------
    # See:   https://stackoverflow.com/questions/68693938/
    #          qgis-how-to-permanently-change-the-order-of-records-
    #          sort-records-features-perm
    #------------------------------------------------------------------
    # Note:  To get a TSV file for both Ref and Non-ref basins,
    #        we need to use a set of basin polygon shapefiles,
    #        one for each "region", and they each need to be
    #        sorted by GAGE_ID.
    #------------------------------------------------------------------ 
    start_time = time.time()
    print('Running...')
    
    if (data_dir is None):
        data_dir = get_gages2_data_dir('root')   
    if (point_dir is None):
        point_dir = get_gages2_data_dir('points')
    if (polygon_dir is None):
        polygon_dir = get_gages2_data_dir('polygons')

    basin_option = basin_option.upper()
    if (tsv_file is None):
        if (basin_option == 'REF'):
            tsv_file = 'new_gages2_ref.tsv'
        elif (basin_option == 'ALL'):
            tsv_file = 'new_gages2_all.tsv' 
    tsv_path = data_dir + tsv_file
    
    #---------------------------------------------    
    # Open shapefile for the gauge/outlet points
    # This doesn't have the basin polygon.
    #---------------------------------------------    
    point_shp_path = point_dir + point_shp_file
    point_prj_path = point_shp_path.replace('.shp', '.prj')
    point_unit     = su.open_shapefile( point_shp_path )
    ## point_unit  = ogr.Open( point_shp_path )
    point_layer = point_unit.GetLayer()

    #------------------------------------------    
    # Info for the set of regional shapefiles
    #------------------------------------------
    if (basin_option == 'REF'):
        poly_shp_files = [ 'bas_ref_all_sorted.shp' ]
    elif (basin_option == 'ALL'):
        poly_shp_files = [
            'bas_ref_all_sorted.shp',
            'nonref_sorted/bas_nonref_AKHIPR_sorted.shp',
            'nonref_sorted/bas_nonref_AKHIPR_sorted.shp',
            'nonref_sorted/bas_nonref_AKHIPR_sorted.shp',
            'nonref_sorted/bas_nonref_CntlPlains_sorted.shp',
            'nonref_sorted/bas_nonref_EastHghlnds_sorted.shp',
            'nonref_sorted/bas_nonref_MxWdShld_sorted.shp',
            'nonref_sorted/bas_nonref_NorthEast_sorted.shp',
            'nonref_sorted/bas_nonref_SECstPlain_sorted.shp',
            'nonref_sorted/bas_nonref_SEPlains_sorted.shp',
            'nonref_sorted/bas_nonref_WestMnts_sorted.shp',
            'nonref_sorted/bas_nonref_WestPlains_sorted.shp',
            'nonref_sorted/bas_nonref_WestXeric_sorted.shp' ]
    poly_shp_keys = [
        'Ref', 'Alaska', 'Hawaii', 'PuertoRico', 'CntlPlains',
        'EastHghlnds', 'MxWdShld', 'NorthEast', 'SECstPlain',
        'SEPlains', 'WestMnts', 'WestPlains', 'WestXeric' ]
    #---------------------------------------------
    # Create two Python dictionaries for regions
    # Open shapefiles for the basin boundaries.
    # These have the basin polygons.
    #---------------------------------------------
    poly_layers    = dict()
    poly_prj_paths = dict()
    poly_shp_units = list()
    for k in range(len(poly_shp_files)):
        # print('poly_shp_files[k] =', poly_shp_files[k])
        poly_shp_path = polygon_dir + poly_shp_files[k]    
        poly_prj_path = poly_shp_path.replace('.shp', '.prj')
        #-----------------------------------------------
        # Units must have unique names and remain open
        # so that layer.GetNextFeature method works.
        #-----------------------------------------------   
        ## poly_shp_unit = ogr.Open( poly_shp_path )
        ## poly_layer = poly_shp_unit.GetLayer()
        #----------------------------------------------- 
        ## poly_shp_units.append( ogr.Open( poly_shp_path ))
        poly_shp_units.append( su.open_shapefile( poly_shp_path ))
        poly_layer = poly_shp_units[k].GetLayer()
        ## poly_shp_unit = None  # don't close the files here
        poly_key      = poly_shp_keys[k]
        poly_layers[ poly_key ]    = poly_layer
        poly_prj_paths[ poly_key ] = poly_prj_path

    if (DEBUG):
        print('## basin_option =', basin_option)
        print('## nf_max =', nf_max)
             
    #--------------------    
    # Open new TSV file
    #--------------------
    # FIRST CHECK IF FILE EXISTS AND OK TO OVERWRITE?  #######
    tsv_unit = open( tsv_path, 'w') 
    if (ADD_BOUNDING_BOX):
        ## new_headings=['MINLON','MAXLON','MINLAT','MAXLAT']
        new_headings=['MINLON','MAXLON','MINLAT','MAXLAT','AREA','PERIMETER']
    n_features_tsv = np.int32(0)

    insert_key   = 'LNG_GAGE'
    filter_key   = None
    filter_value = None

    #---------------------------------------------    
    # Iterate over point features in point layer
    #---------------------------------------------
    bad_box_count = np.int32(0)
    n_features    = np.int32(0)
    ## poly_unit_num = np.int32(0)
    
    bad_station_IDs = list()
    long_gage_IDs   = list()
    station_ID      = '00000000'
    gage_ID         = '00000000'
        
    for point_feature in point_layer:
        ## print('n_features =', n_features)
        point_geom = point_feature.GetGeometryRef()
        point_atts = point_feature.items()  # This is a Python dictionary
        # n_rings    = point_geom.GetGeometryCount()

        #--------------------------------
        # Write header for new TSV file
        #--------------------------------
        if (n_features == 0):
            su.write_tsv_header( tsv_unit, point_atts, insert_key=insert_key,
                  new_headings=new_headings )            
              
        #----------------------------------------
        # Print attributes in attribute table ?
        #----------------------------------------
        if (REPORT):
            su.print_attributes( point_atts, n_features )

        #-------------------------------------- 
        # Skip to polygon for this station ID
        #--------------------------------------
        last_station_ID = station_ID
        station_ID      = point_atts['STAID']
        station_class   = point_atts['CLASS']
        REF_ONLY = (basin_option == 'REF')
        if (REF_ONLY or (station_class == 'Ref')):
            poly_file_key = 'Ref'
        else:
            poly_file_key = point_atts['AGGECOREGI']

        #------------------------------------------------------
        # Maybe write a separate function that first scans
        # point_shp_file to see if it is sorted by station_ID
        #------------------------------------------------------
        if (station_ID < last_station_ID):
            print('#### WARNING:  station IDs are not sorted.')
            print('####    last_station_ID =', last_station_ID)
            print('####    station_ID      =', station_ID)
            print()

        if (REF_ONLY and (station_class.lower() == 'non-ref')):
            n_features += 1
            if (DEBUG):
                print('Skipping non-ref station ID =', station_ID)
            continue

        if (DEBUG):
            print()
            print('Finding match for station ID =', station_ID)

        #-------------------------------------------------------
        # Note:  Both gage_IDs and station_IDs must be sorted.
        #-------------------------------------------------------
        FOUND_POLY_MATCH = False  # the default for any break
        ## while (gage_ID != station_ID):
        while (gage_ID < station_ID):
            ## poly_feature = get_poly_feature_v0( poly_unit_num, poly_layers)
            if (DEBUG):
                print('## Calling get_poly_feature...')
            poly_feature = get_poly_feature( poly_file_key, poly_layers)
            last_gage_ID = gage_ID
            if (DEBUG):
                print('## Calling get_gate_ID...')
            gage_ID      = get_gage_ID( poly_feature, DEBUG )
            check_gage_ID( gage_ID, last_gage_ID, long_gage_IDs, DEBUG)       

            #----------------------------------------------
            # Get basin area and perimeter from poly_file
            #----------------------------------------------
            poly_atts = poly_feature.items()  # a Python dictionary      
            area      = poly_atts['AREA']
            perimeter = poly_atts['PERIMETER']
    
            # This was needed before sorting features by GAGE_ID.
            OVERSHOOT = check_overshoot( gage_ID, station_ID, bad_station_IDs)
            if (OVERSHOOT): break

            if (gage_ID == station_ID):
                FOUND_POLY_MATCH = True
                if (DEBUG):
                    print('   Found match.') 

        #------------------------------------------ 
        # Get all points in basin polygon
        # Convert coords to Geo. lon/lat (WGS-84)
        # Get geographic bounding box
        #------------------------------------------
        if (FOUND_POLY_MATCH):
            ## poly_prj_path = poly_prj_paths[ poly_unit_num ]
            poly_prj_path = poly_prj_paths[ poly_file_key ]
            if (DEBUG):
                print('   Computing bounding box.')
                print()
            x1, y1 = su.get_polygon_points( poly_feature )
            x2, y2 = su.convert_coords(x1,y1, inPRJfile=poly_prj_path,
                                       SWAP_XY=SWAP_XY, PRINT=False)
            minlon, maxlon, minlat, maxlat = su.get_bounding_box(x2, y2)

            IN_RANGE = check_bounding_box( point_atts, minlon,
                             maxlon, minlat, maxlat)
            if not(IN_RANGE):
               bad_box_count += 1

        #--------------------------------------------              
        # Apply some test or filter to this feature
        # to decide whether to write it to TSV file
        #--------------------------------------------
        if (filter_key is None):
            ALL_GOOD = True
        else:
            value = point_atts[ filter_key ]
            if (isinstance(value, str)):
                ALL_GOOD = (value.lower() == filter_value.lower())
            else:
                ALL_GOOD = (value == filter_value)
                        
        #-------------------------------------              
        # Write line of info to new TSV file
        #-------------------------------------
        if (ALL_GOOD and FOUND_POLY_MATCH):
            new_values = []  # empty list;  default
            if (ADD_BOUNDING_BOX):
                ## new_values = [minlon,maxlon,minlat,maxlat]
                new_values = [minlon,maxlon,minlat,maxlat,area,perimeter]
            su.write_tsv_line( tsv_unit, point_atts, insert_key=insert_key,
                               new_values=new_values )
            n_features_tsv += 1

        n_features += 1
        if (n_features >= nf_max):   # May need >= vs. == now.
            break
            
    #--------------------------------
    # Print info, close files, etc.
    #--------------------------------
    print('Processed', n_features, 'features.')
    print('Wrote', n_features_tsv, 'rows to tsv file.')
    print('Bad station IDs =')
    print( bad_station_IDs )
    print('Missing station count  =', len(bad_station_IDs))
    print('Bad bounding box count =', bad_box_count)
    print('   (outlet outside bounding box)')
    print('Long gage IDs =')
    print( long_gage_IDs )
    tsv_unit.close()
    run_time = (time.time() - start_time)
    print('Run time =', run_time, '[sec]')
    print('Finished.')
         
#   create_tsv_from_shapefile()
#---------------------------------------------------------------------
def get_poly_feature_v0( poly_unit_num, poly_layers ):

    n_poly_units  = len( poly_layers )
    poly_layer    = poly_layers[ poly_unit_num ]
    ## n_polys       = poly_layer.GetFeatureCount()
    poly_feature  = poly_layer.GetNextFeature()
          
    if (poly_feature is None):
        #----------------------------------------------
        # Reached EOF; prepare to read from next file
        #----------------------------------------------
        # print('## WARNING: poly_feature is None.') 
        poly_unit_num += 1
        if (poly_unit_num >= n_poly_units):
            print('## ERROR: poly_feature is None.') 
            return None
        poly_layer    = poly_layers[ poly_unit_num ]
        poly_feature  = poly_layer.GetNextFeature()

    return poly_feature

#   get_poly_feature_v0()
#---------------------------------------------------------------------
def get_poly_feature( poly_file_key, poly_layers ):

    ## print('poly_file_key =', poly_file_key)
    poly_layer   = poly_layers[ poly_file_key ]
    ## print('poly_layer =', poly_layer)
    poly_feature = poly_layer.GetNextFeature()
    ## print('poly_feature =', poly_feature)
          
    if (poly_feature is None):
        #----------------------------------------------
        # Reached EOF; prepare to read from next file
        #----------------------------------------------
        print('## WARNING: No more features in file.')
        print('##          poly_feature is None.')

    return poly_feature  # May be None

#   get_poly_feature()
#---------------------------------------------------------------------
def get_gage_ID( poly_feature, DEBUG ):
 
    if (poly_feature is None):
        return '99999999'

    #-----------------------------------------------
    # Only poly atts are: AREA, PERIMETER, GAGE_ID
    #-----------------------------------------------
    ## poly_geom = poly_feature.GetGeometryRef()
    poly_atts = poly_feature.items()  # a Python dictionary      
    gage_ID   = poly_atts['GAGE_ID']
    ## gage_ID   = str(gage_ID).strip()   ##### Not needed

    if (DEBUG):
        print('   gage_ID =', gage_ID)

    return gage_ID
    
#   get_gage_ID()
#---------------------------------------------------------------------
def check_gage_ID( gage_ID, last_gage_ID, long_gage_IDs, DEBUG ):

    if (gage_ID == '99999999'):
        print('### ERROR: Station ID not found in')
        print('###        any of the polygon files.')
        ### print('### station_ID =', station_ID)
        print()
                
    if (gage_ID < last_gage_ID):
        print('### WARNING:  gage IDs are not sorted.')
        print('###    last_gage_ID =', last_gage_ID)
        print('###    gage_ID      =', gage_ID)
        print()

    if (len(gage_ID) > 8):
        if (DEBUG):
            print('### WARNING: Long gage ID.')
            print('###    gage_ID =', gage_ID, '(', len(gage_ID), 'chars)')
            print()
        long_gage_IDs.append( gage_ID )
        
#   check_gage_ID()
#---------------------------------------------------------------------
def check_overshoot( gage_ID, station_ID, bad_station_IDs):

    OVERSHOOT = (gage_ID > station_ID)  # comparing strings        
    if (OVERSHOOT):
        #-------------------------------------------------
        # This should only happen if a ref ID is missing
        # from the basin boundary shapefile or gage_IDs
        # are not sorted.  But it occurs for 18 stations:
        # e.g. station_ID = 01031500, 01047000
        #-------------------------------------------------
        # Note: Most IDs are 8 chars; some are 9+ chars.
        # e.g. 402114105350101
        #-------------------------------------------------
        ## if (REPORT or DEBUG):
        if (True):
            print('##   WARNING: gage_ID > station_ID')
            print('##            gage_ID    =', gage_ID)
            print('##            station_ID =', station_ID)
        bad_station_IDs.append( station_ID )

    return OVERSHOOT
    
#   check_overshoot()
#---------------------------------------------------------------------
def check_bounding_box( point_atts, minlon, maxlon, minlat, maxlat):

    IN_RANGE = True
    gage_lon = point_atts['LNG_GAGE']
    gage_lat = point_atts['LAT_GAGE']                  
    if ((gage_lon < minlon) or (gage_lon > maxlon) ):
        IN_RANGE = False          
    if ((gage_lat < minlat) or (gage_lat > maxlat) ):
        IN_RANGE = False    
    return IN_RANGE

#   check_bounding_box()
#---------------------------------------------------------------------
def check_basin_poly_file():

    ############################################
    # Check for problems reading polygon file
    # File is not sorted by gage_ID.
    #   e.g. 01031500, 01047000, 01187800.
    # But QGIS shows these w/ View Att. Table
    ############################################
    
    #--------------------------------------------
    # You will need to modify these directories
    # and get all the necessary files.
    #--------------------------------------------
    polygon_dir = get_gages2_data_dir('polygons')
    polygon_file = 'bas_ref_all.shp'
    ## polygon_file = 'bas_ref_all_sorted.shp'
    polygon_path = polygon_dir + polygon_file
      
    n_features = np.int32(0)
    poly_unit  = ogr.Open( polygon_path )
    poly_layer = poly_unit.GetLayer()
    gage_ID    = '00000000'
    while (True):    
        poly_feature = poly_layer.GetNextFeature()
        if (poly_feature is None):
            break        
        #-----------------------------------------------
        # Only poly atts are: AREA, PERIMETER, GAGE_ID
        #-----------------------------------------------
        ## poly_geom = poly_feature.GetGeometryRef()
        poly_atts = poly_feature.items()  # a Python dictionary
        last_gage_ID = gage_ID      
        gage_ID   = poly_atts['GAGE_ID']
        if (last_gage_ID > gage_ID):
            print('WARNING: gage_IDs are not sorted.')
            print('   last_gage_ID =', last_gage_ID)
            print('   gage_ID      =', gage_ID)
        print('gage_ID =', gage_ID)
        n_features += 1
    
    print('Number of features =', n_features)

#   check_basin_poly_file()   
#---------------------------------------------------------------------



