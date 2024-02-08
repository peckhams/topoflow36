
#---------------------------------------------------------------------
# Copyright (c) 2023-2024, Scott D. Peckham
#
# Jan 2024. Renamed all "ngen/utils" files to end in "_utils.py"
#           instead of "_tools.py"
#           Modified to use new data_utils.py.
# Nov 2023. Wrote get_swb_points to plot points on a map.
# Oct 2023. Moved get_swb_class_names() and get_swb_class() to
#           swb_utils.py.      
#           Added merge_SB_info_for_regions(),
#           get_csv_headings(), compute_swb_classes(),
#           get_snow_precip_fraction(), get_precip_timing_index(),
#           get_aridity_index(), get_swb_class_names(),
#           get_swb_class().      
# Jun 2023. Started from hlr_utils.py and wrote a more general
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
#  >>> from topoflow.utils.ngen import gages2_utils as g2
#  >>> g2.create_tsv_from_shapefile( nf_max=100, tsv_file="TEST.tsv" )
#  >>> g2.create_tsv_from_shapefile( nf_max=100, tsv_file="TEST2.tsv",
#         basin_option='ALL')
#  >>> g2.create_tsv_from_shapefile( nf_max=10000 )
#  >>> g2.compute_swb_classes(ORIGINAL=False)
#
#---------------------------------------------------------------------
#
#  get_gages2_data_dir()
#  create_tsv0()          # OBSOLETE?
#      
#  create_tsv_from_shapefile()
#  get_poly_feature()
#  get_gage_ID()
#  check_gage_ID()
#  check_overshoot()
#  check_bounding_box()
#  check_basin_poly_file()
#
#  The following functions are for GAGES-II "Selected Basins":
#  https://www.sciencebase.gov/catalog/item/5a4e9bcbe4b0d05ee8c6643f
#  which are a subset of the Reference basins.
#
#  merge_SB_info_for_regions()
#  get_csv_headings()
#  compute_swb_classes()
#     get_snow_precip_fraction()
#     get_precip_timing_index()
#     get_aridity_index()
#  get_swb_points()
#
#---------------------------------------------------------------------

import numpy as np
from osgeo import ogr, osr
import json, sys, time

from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import shape_utils as su
from topoflow.utils.ngen import swb_utils as swb

#---------------------------------------------------------------------
def get_gages2_data_dir( gtype='root' ):

    data_dir  = dtu.get_data_dir('USGS_GAGES2_all')
    if (gtype == 'root'):
        pass
    elif (gtype == 'points'):
       data_dir += 'gagesII_9322_point_shapefile/'
    elif (gtype == 'polygons'):
       data_dir += 'boundaries-shapefiles-by-aggeco/'
    elif (gtype == 'SB3'):
       data_dir = dtu.get_data_dir('USGS_GAGES2_SB3')
       
    return data_dir
       
#   get_gages2_data_dir()
#---------------------------------------------------------------------
# def create_tsv0( data_dir=None,
#                  shape_file='gagesII_9322_sept30_2011.shp',                    
#                  prj_file  ='gagesII_9322_sept30_2011.prj',
# #                shape_file='bas_ref_all.shp',                    
# #                prj_file  ='bas_ref_all.prj',
#                  tsv_file='new_gages2_all.tsv',
# #                 tsv_file='new_gages2_conus.csv',
#                  nf_max=50, REPORT=True, basin_option='ALL'):
# 
#     #--------------------------------------------------------
#     # Note:  Attributes in shapefile attribute table for
#     #        gagesII_9322_sept30_2011.shp are:
#     # 
#     # g2_staid      = attributes['STAID']    # station ID  
#     # g2_staname    = attributes['STANAME']  # station name
#     # g2_class      = attributes['CLASS']    # Ref or Non-ref
#     # g2_aggecoregi = attributes['AGGECOREGI'] 
#     # g2_drain_sqkm = attributes['DRAIN_SQKM'] 
#     # g2_huc02      = attributes['HUC02'] 
#     # g2_lat_gage   = attributes['LAT_GAGE'] 
#     # g2_lng_gage   = attributes['LNG_GAGE'] 
#     # g2_state      = attributes['STATE'] 
#     # g2_hcdn_2009  = attributes['HCDN_2009']
#     # g2_active09   = attributes['ACTIVE09'] 
#     # g2_flyrs1900  = attributes['FLYRS1900']                        
#     # g2_flyrs1950  = attributes['FLYRS1950'] 
#     # g2_flyrs1990  = attributes['FLYRS1990'] 
#     #--------------------------------------------------------
#     if (data_dir is None):
#         data_dir = get_gages2_data_dir()
# 
#     insert_key = 'LNG_GAGE'
#     if (basin_option == 'REF'):
#         filter_key   = 'CLASS'
#         filter_value = 'Ref'
#     elif (basin_option == 'NONREF'):
#         filter_key   = 'CLASS'
#         filter_value = 'Non-ref'
#     elif (basin_option == 'ALL'):
#         filter_key   = None
#         filter_value = None
# 
#     su.create_tsv_from_shapefile(data_dir=data_dir, shape_file=shape_file,
#                       prj_file=prj_file, nf_max=nf_max,
#                       insert_key=insert_key,
#                       ADD_BOUNDING_BOX=True,
#                       filter_key=filter_key, filter_value=filter_value)
# 
# #   create_tsv0()
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
    # Note:  In the "point" shapefile, which includes all 9322 basins,
    #        the USGS Site ID is called "STAID" (station ID).
    #        In all the "polygon" shapefiles, the USGS Site ID is
    #        called "GAGE_ID".  There is one polygon shapefile for
    #        all the Reference basins, and for NonReference basins
    #        there is one polygon shapefile for each ecoregion.
    #------------------------------------------------------------------
    # Note:  This routine will not work unless the features in
    #        point_shp_file are sorted by station_ID and the features
    #        in basin polygon files are sorted by gage_ID.
    #        It turned out that this was the case for point_shp_file
    #        but needed to sort features in the basin polygon files
    #        by gage_ID using QGIS as follows:
    #
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

    if (point_dir is None):
        point_dir = get_gages2_data_dir('points')
    #---------------------------------------------------------
    # Use the new, ID-sorted polygon shapefiles we created.
    # It appears that the original G2 "point shapefile" has
    # already been sorted by USGS site ID.
    #--------------------------------------------------------- 
#     if (polygon_dir is None):
#         polygon_dir = get_gages2_data_dir('polygons')

    basin_option = basin_option.upper()
    if (basin_option == 'ALL'):
        key = 'USGS_GAGES2_all'
    elif (basin_option == 'NONREF'):
        key = 'USGS_GAGES2_nonref'
    elif (basin_option == 'REF'):
        key = 'USGS_GAGES2_ref'
    else:
        key = 'USGS_GAGES2_SB3' 
    new_dir = dtu.get_new_data_dir( key )

    if (tsv_file is None):
        bop = basin_option.lower()
        tsv_file = 'new_gages2_' + bop + '.tsv'
    tsv_path = new_dir + tsv_file
    
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
    poly_ref_dir = dtu.get_new_data_dir( 'USGS_GAGES2_ref')
    if (basin_option == 'REF'):
        poly_shp_files = [ 'bas_ref_all_sorted.shp' ]
    elif (basin_option == 'ALL'):
        poly_nonref_dir = dtu.get_new_data_dir( 'USGS_GAGES2_nonref')
        poly_shp_files = [
            'bas_ref_all_sorted.shp',
            #------------------------------------
            'bas_nonref_AKHIPR_sorted.shp',
            'bas_nonref_CntlPlains_sorted.shp',
            'bas_nonref_EastHghlnds_sorted.shp',
            'bas_nonref_MxWdShld_sorted.shp',
            'bas_nonref_NorthEast_sorted.shp',
            'bas_nonref_SECstPlain_sorted.shp',
            'bas_nonref_SEPlains_sorted.shp',
            'bas_nonref_WestMnts_sorted.shp',
            'bas_nonref_WestPlains_sorted.shp',
            'bas_nonref_WestXeric_sorted.shp' ]
    poly_shp_keys = [
        'Ref', 'AKHIPR', 'CntlPlains',
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
        shp_filename = poly_shp_files[k]
        if ('_nonref_' in shp_filename):   ######
            polygon_dir = poly_nonref_dir
        else:
            polygon_dir = poly_ref_dir
        poly_shp_path = polygon_dir + shp_filename   
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
    bad_box_IDs   = list()
    n_features    = np.int32(0)
    
    bad_station_IDs = list()
    long_gage_IDs   = list()
    station_ID      = '00000000'
    gage_ID         = '00000000'
    not_conus_list  = ['Alaska', 'Hawaii', 'PuertoRico']

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
            #------------------------------------------------------
            # See "poly_shp_keys" above for a list of these keys.
            # Note that Alaska, Hawaii, and Puerto Rico are all
            # in one file, despite appearing as region names.
            #------------------------------------------------------
            region_name = point_atts['AGGECOREGI']
            if (region_name in not_conus_list):
                poly_file_key = 'AKHIPR'
            else:
                poly_file_key = region_name
#             poly_file_key = point_atts['AGGECOREGI']
#             if (poly_file_key in not_conus_list):
#                 poly_file_key = 'AKHIPR'

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
        FOUND_POLY_MATCH = False  # the default
        if (DEBUG):
            print('## Calling get_poly_feature...')
        #-----------------------------------------------------------
        # Note: Here we're using the poly_file_key we read from
        #       the "point shapefile" that has all station IDs
        #       in sorted order to get info from the corresponding
        #       basin "polygon shapefile".
        #-----------------------------------------------------------
        # Note: If EOF is reached, get_poly_feature returns None
        #       and then get_gage_ID returns '99999999', which
        #       is greater than any other gage ID.
        #-----------------------------------------------------------
        poly_feature = get_poly_feature( poly_file_key, poly_layers)
        last_gage_ID = gage_ID
        if (DEBUG):
            print('## Calling get_gage_ID...')
        gage_ID = get_gage_ID( poly_feature, DEBUG )
        check_gage_ID( gage_ID, last_gage_ID, long_gage_IDs, DEBUG)       

        #----------------------------------------------
        # Get basin area and perimeter from poly_file
        #----------------------------------------------
        poly_atts = poly_feature.items()  # a Python dictionary      
        area      = poly_atts['AREA']
        perimeter = poly_atts['PERIMETER']

        #------------------------------------------------------
        # This was needed before sorting features by GAGE_ID
        # and we used a while loop.
        #------------------------------------------------------
#         OVERSHOOT = check_overshoot( gage_ID, station_ID, bad_station_IDs)
#         if (OVERSHOOT): break

        if (gage_ID == station_ID):
            #------------------------------------------
            # Found matching polygon. 
            # Get all points in basin polygon
            # Convert coords to Geo. lon/lat (WGS-84)
            # Get geographic bounding box
            #------------------------------------------
            FOUND_POLY_MATCH = True
            poly_prj_path = poly_prj_paths[ poly_file_key ]
            if (DEBUG):
                print('   Found match.') 
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
                bad_box_IDs.append( station_id )

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
            if ((n_features_tsv % 100) == 0):
                print('n_features written =', n_features_tsv)

        n_features += 1
        if (n_features >= nf_max):   # May need >= vs. == now.
            break
            
    #--------------------------------
    # Print info, close files, etc.
    #--------------------------------
    print('Processed', n_features, 'features.')
    print('Wrote', n_features_tsv, 'rows to tsv file:')
    print('   ' + tsv_file)
#     print('Bad station IDs =')
#     print( bad_station_IDs )
#     print('Missing station count  =', len(bad_station_IDs))
    print('Bad bounding box count =', bad_box_count)
    print('   (outlet outside bounding box)')
    print('Bad bounding box IDs =')
    print( bad_box_IDs )
    print('Long gage IDs =')
    print( long_gage_IDs )
    tsv_unit.close()
    run_time = (time.time() - start_time)
    print('Run time =', run_time, '[sec]')
    print('Finished.')
         
#   create_tsv_from_shapefile()
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

    #-------------------------------------------------------------- 
    # Note:  If EOF was reached for a polygon shapefile, then
    #        poly_feature will be None and this function returns
    #        '99999999', which is greater than any other gage ID.
    #--------------------------------------------------------------
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

    #-----------------------------------------------------------------     
    # Note: This was getting triggered twice for basin_option='ALL',
    #       which includes the NonRef basins.  The issue had to do
    #       with reading the AKHIPR file 3 times vs. 1, since
    #       Alaska, Hawaii, and PuertoRico appear as distinct
    #       regions.  The issue is now fixed.
    #       But found no problem in the final TSV file.
    #
    #       In the first case,
    #       last_gage_ID = 16010000 and gage_ID = 15019990.
    #       The gage ID: 16010000 is a Ref basin in Hawaii, 
    #       The gage ID: 15019990 is a NonRef basin in Alaska,
    #       and is the first one in the file:
    #           bas_nonref_AKHIPR_sorted.shp.
    #       The ID of the last one in that file is: 50148890.
    
    #       In the second case,
    #       last_gage_ID = 402913084285400 and gage_ID = 15019990.
    #       The gage ID: 402913.." is a NonRef basin in CntlPlains,
    #       and is the last one in the file:
    #           bas_nonref_CntlPlains_sorted.shp.
    #       The ID of the first one in that file is: 03144816.     
    #-----------------------------------------------------------------           
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
def merge_SB_info_for_regions( out_csv_file=None, SORT_BY_ID=True,
                               REPORT=False):

    #----------------------------------------------------------------
    # Note: SB stands for "Selected Basins".  Info for these basins
    #       is saved in 19 different CSV files, for "regions".
    #       Region numbers are 1 to 9, 10L, 10U, and 11 to 18.
    #       Info in the new, merged CSV file will be used to
    #       compute the Seasonal Water Balance classes.
    #       Here we are using the "untransformed" values of vars,
    #       stored in CSV files with "untransfBC" in the name.
    #       These CSV files were exported from an Excel file that
    #       had multiple sheets: SB_untransfBCs.new.xlsx.
    #----------------------------------------------------------------
    #       Maybe modify so output file is tab-delimited (TSV).
    #       The input files for regions are CSV files.
    #----------------------------------------------------------------
    delim    = ','  # CSV files
    new_dir  = dtu.get_new_data_dir( 'USGS_GAGES2_SB3' )
    data_dir = new_dir + 'by_region/'

    if (out_csv_file is None):
        out_csv_file = 'SB3_all_regions_untransfBCs.csv'
    out_csv_path = new_dir + out_csv_file
    out_csv_unit = open( out_csv_path, 'w')

    region_list = ['1','2','3','4','5','6','7','8','9', '10L','10U',
                  '11','12','13','14','15','16','17','18']
    n_basins   = 0
  
    print('Working...')  
    for k in range(19):  
        region_str = region_list[k]   
        csv_file   = 'SB3_region' + region_str + '_untransfBCs.new.csv'
        csv_path   = data_dir + csv_file
        csv_unit   = open( csv_path, 'r' )
        if (REPORT):
            print('Merging info from file:')
            print('   ' + csv_file)

        #-----------------------------------------------
        # Note: Not all files have same header !!
        # For missing columns, see the file:
        #    SB3_BCs_byRegion.xlsx.
        #-----------------------------------------------       
        headings = get_csv_headings( csv_unit )
        headings.insert(2, 'REGION')
        #-----------------------------------------------
        # Create header for new file w/ extra headings
        # The first region/sheet is only missing HGBD.
        #-----------------------------------------------
        if (k == 0):
            all_headings = headings.copy()
            all_headings.insert(36, 'HGBD')
            ## print('### k=0:  all_headings =', all_headings)
            header = delim.join( all_headings ) + '\n'
            out_csv_unit.write( header )

        while (True):
            line = csv_unit.readline()
            if (line == ''):
                break   # (reached end of file)
            #----------------------------------------
            # Need to remove the newline char first
            #----------------------------------------
            line = line[:-1]
            values = line.split( delim )
            values.insert(2, region_str)
            #------------------------------------------------
            # Create a dictionary to map headings to values
            #------------------------------------------------
            val_dict = dict()
            for j in range(len(headings)):
                val_dict[ headings[j] ] = values[j]
            #----------------------------------------
            # For each heading in all_headings,
            # write corresponding value to new file
            #----------------------------------------
            out_line = ''
            for h in all_headings:
                try:
                   value = val_dict[h]
                except:
                   #--------------------------------------
                   # Heading & value missing in this CSV
                   #--------------------------------------
                   # Regions 8 and 12 don't provide the
                   # variable:  SNOW_PCT_PRECIP
                   #--------------------------------------
                   value = '-9999'
                out_line += value + ','
            out_line = out_line[:-1] + '\n'
            out_csv_unit.write( out_line ) 
            n_basins += 1
            
        csv_unit.close()
       
    out_csv_unit.close()
    print('Wrote info for', n_basins, 'basins to new CSV file:')
    print('  ' + out_csv_path)
    print('Finished merging info for SB regions.')
    print()

    if (SORT_BY_ID):
        #-------------------------------------------
        # Read all lines from new, merged CSV file
        #-------------------------------------------
        csv_unit1 = open( out_csv_path, 'r' )
        lines = csv_unit1.readlines()  # read all lines
        csv_unit1.close()
        #--------------------------------------------------
        # Sort lines by USGS station ID (in first column)
        #--------------------------------------------------
        header = lines[0]   # header line  
        lines  = lines[1:]    
        lines.sort()
        #------------------------------------
        # Save sorted lines to new CSV file
        #------------------------------------
        sorted_csv_path = out_csv_path.replace('.csv', '_sorted.csv')
        print('Sorting lines by USGS station ID to create new file:')
        print('  ' + sorted_csv_path)    
        csv_unit2 = open( sorted_csv_path, 'w' )
        csv_unit2.writelines( [header] + lines )
        csv_unit2.close()
        print('Finished.')
        print()
         
#   merge_SB_info_for_regions()
#---------------------------------------------------------------------
def get_csv_headings( csv_unit, delim=',' ):

    #------------------------------------------------------------
    # Note: Info for "selected basins" (sb) is saved in 19
    #       CSV files, one for each of 19 "regions".
    #       Region numbers are 1 to 9, 10L, 10U, & 11 to 18.
    #       However, not all files have the same col. headings.
    #       See: SB3_BCs_byRegion.xlsx. for missing columns.  
    #------------------------------------------------------------ 
    header   = csv_unit.readline()
    headings = header.split( delim )
    for k in range( len(headings) ):
        headings[k] = headings[k].strip()
    return headings
   
#   get_csv_headings()
#---------------------------------------------------------------------
def compute_swb_classes( in_csv_file=None, ORIGINAL=False,
                         NEW_TSV=False, USE_B3=False, REPORT=False ):

    #---------------------------------------------------------------
    # Note: In GAGES-II, extra information is provided for a set
    #       1947 "Selected Basins".  Recall that the GAGES-II data
    #       set includes 9322 basins, 2057 of which are classified
    #       as Reference basins (vs. Non-Reference).
    #       This function uses the info for the Selected Basins
    #       to classify these basins using the Seasonal Water
    #       Balance method described in Berghuijs et al. (2014).
    #       It uses 3 parameters: precip_timing_index (delta_p), 
    #       snow_precip_fraction (f_s), and aridity_index (phi).
    #       See the function:  get_swb_class() in swb_utils.py.
    #--------------------------------------------------------------
    #  In filenames, "SB" stands for "Selected Basins", and
    #  "untransf" is short for "untransformed" values of vars.
    #--------------------------------------------------------------
    #  If USE_B3, a new class, B3, is used so that B1, B2, and
    #  B3 "stack up" along phi-axis just as A1, A2, and A3 do.
    #--------------------------------------------------------------  
    if (in_csv_file is None):
        in_csv_file = 'SB3_all_regions_untransfBCs_sorted.csv'
    delim       = ','   # CSV file
    new_dir     = dtu.get_new_data_dir( 'USGS_GAGES2_SB3' )
    in_csv_path = new_dir + in_csv_file
    in_csv_unit = open( in_csv_path, 'r' )
    headings    = get_csv_headings( in_csv_unit, delim=delim )
            
    if (NEW_TSV):
        delim2 = '\t'  # TSV file
        out_tsv_file = in_csv_file.replace('.csv', '_SWB.tsv')
        out_tsv_path = new_dir + out_tsv_file
        out_tsv_unit = open( out_tsv_path, 'w' )
        #------------------------------------------------------------------
        # Headings to keep or rename (could remove others for now):
        # StnID, REGION, DRAIN_SQKM, LAT_CENT, LONG_CENT, PPTAVG_BASIN,
        # T_AVG_BASIN, RH_BASIN, PET (mm), SNOW_PCT_PRECIP, 
        # PRECIP_SEAS_IND ????, JUL_JAN_TMP_DIFF_DEGC, ELEV_RANGE_M_BASIN
        #------------------------------------------------------------------
        # Headings to add:
        # PRECIP_TIMING_INDEX, ARIDITY_INDEX, SWB_CLASS, SWB2_CLASS
        #------------------------------------------------------------------
        headings2 = headings + ['PRECIP_TIMING_INDEX', 'ARIDITY_INDEX',
                        'SNOW_PRECIP_FRACTION', 'SWB1_CLASS', 'SWB2_CLASS']
        #---------------------------------
        # Write header for new TSV file
        #---------------------------------
        tsv_header = ''
        for heading in headings2:
            tsv_header += (heading + delim2)
        tsv_header = tsv_header[:-1] + '\n' # remove delim, add newline
        out_tsv_unit.write( tsv_header )
        
    fs_min  = 1000.0
    fs_max  = -1000.0
    phi_min = 1000.0
    phi_max = -1000.0
    dp_min  = 1000.0
    dp_max  = -1000.0
    #-------------------
    n_basins = 0
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0

    class_count1 = dict()
    class_count2 = dict()    
    names = swb.get_swb_class_names( USE_B3=USE_B3 )
    for name in names:
        class_count1[ name ] = 0
        class_count2[ name ] = 0

    while (True):
        line = in_csv_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        #----------------------------------------
        # Need to remove the newline char first
        #----------------------------------------
        line = line[:-1]
        values = line.split( delim )
        #------------------------------------------------
        # Create a dictionary to map headings to values
        #------------------------------------------------
        val_dict = dict()
        for j in range(len(headings)):
            val_dict[ headings[j] ] = values[j]

        #---------------------------------------------------
        # Note: Using results from Berghuijs et al. paper:
        #       delta_p = precip_timing_index
        #       f_s     = snow_precip_fraction
        #       phi     = aridity_index
        #---------------------------------------------------
        f_s     = get_snow_precip_fraction( val_dict, UNITS='none',
                                            REPORT=REPORT )
        delta_p = get_precip_timing_index( val_dict, f_s, REPORT=REPORT )
        phi     = get_aridity_index( val_dict, REPORT=REPORT )
 
        #---------------------------------  
        # Count extreme cases of delta_p
        #---------------------------------
        if (delta_p == -0.999):
            count1 += 1
        if (delta_p == 0.999):
            count2 += 1
        if (delta_p == -1.0):
            count3 += 1
        if (delta_p == 1.0):
            count4 += 1
                        
        #---------------------------------------------
        # Save both SWB versions in the new TSV file
        #---------------------------------------------       
        swb1_class = swb.get_swb_class( delta_p, f_s, phi,
                                   ORIGINAL=True, USE_B3=USE_B3)                                   
        swb2_class = swb.get_swb_class( delta_p, f_s, phi,
                                   ORIGINAL=False, USE_B3=USE_B3)
        class_count1[ swb1_class ] += 1   # (swb_class may be 'None')
        class_count2[ swb2_class ] += 1   # (swb_class should never be 'None')               
        n_basins += 1
        
        dp_min  = np.minimum( dp_min, delta_p )
        dp_max  = np.maximum( dp_max, delta_p )
        fs_min  = np.minimum( fs_min, f_s )
        fs_max  = np.maximum( fs_max, f_s )
        phi_min = np.minimum( phi_min, phi )
        phi_max = np.maximum( phi_max, phi )

        #--------------------------------------------
        # Add new heading/value pairs to dictionary
        #---------------------------------------------------
        # f_s now has units = "none" vs. "percent".
        # Use: "SNOW_PRECIP_FRACTION" vs "SNOW_PCT_PRECIP"
        #---------------------------------------------------
        if (NEW_TSV):
            dp_str  = '{x:.6f}'.format(x=delta_p)  # 6 decimal accuracy
            fs_str  = '{x:.6f}'.format(x=f_s)
            phi_str = '{x:.6f}'.format(x=phi)
            #----------------------------------------------       
            val_dict['PRECIP_TIMING_INDEX']  = dp_str
            val_dict['SNOW_PRECIP_FRACTION'] = fs_str
            val_dict['ARIDITY_INDEX']        = phi_str
            val_dict['SWB1_CLASS']           = swb1_class
            val_dict['SWB2_CLASS']           = swb2_class
        
            #----------------------------------------
            # For each heading in headings2,
            # write corresponding value to new file
            #----------------------------------------
            out_line = ''
            for h in headings2:
                value = val_dict[h]
                out_line += value + delim2
            out_line = out_line[:-1] + '\n'
            out_tsv_unit.write( out_line )

    #-------------------
    # Close both files
    #-------------------
    in_csv_unit.close()
    if (NEW_TSV):
        out_tsv_unit.close()
    #---------------------------
    print('del_p_min, del_p_max =', dp_min, ',', dp_max)
    print('fs_min, fs_max       =', fs_min, ',', fs_max)
    print('phi_min, phi_max     =', phi_min, ',', phi_max)
    print()
    print('Number of basins     =', n_basins)
    c1_str = '{x:5d}'.format(x=class_count1['None'])
    c2_str = '{x:5d}'.format(x=class_count2['None'])
    print('Number unclassified  = ' + c1_str + ' (SWB1)' )
    print('Number unclassified  = ' + c2_str + ' (SWB2)' )
    print()
    print('Limiting values of delta_p:')
    print('# basins with delta_p = -0.999 =', count1)
    print('# basins with delta_p = +0.999 =', count2)
    print('# basins with delta_p = -1.000 =', count3, ' (T_bs <= 1)')
    print('# basins with delta_p = +1.000 =', count4, ' (T_bs >= 1)')
    print()
    for name in names:
        cc1_str = '{x:4d}'.format(x=class_count1[name])
        cc2_str = '{x:4d}'.format(x=class_count2[name])
        if (len(name) == 2):
            name = ' ' + name + ' '      
        print('   count[ ' + name + ' ] = v1:' + cc1_str + ',   v2:' + cc2_str)
    if (NEW_TSV):
        print('Wrote extra info to new TSV file: ')
        print('   ' + out_tsv_file)
    print('Finished.')
    print()

#   compute_swb_classes()
#---------------------------------------------------------------------
def get_snow_precip_fraction( val_dict, UNITS='none', REPORT=False ):

    #----------------------------------------
    # If f_s = 0, all precip falls as rain.
    # If f_s = 1, all precip falls as snow.
    #-----------------------------------------------
    # Regions 8 and 12 don't provide the variable: 
    # SNOW_PCT_PRECIP, so it will now be -9999.
    # See Over et al. (2018, Figure 1.)
    # Region 8 has parts of LA, AR, MS and TN.
    # Region 12 is in Texas.
    #-----------------------------------------------
    val = val_dict[ 'SNOW_PCT_PRECIP' ]  # [percent]
    val = np.float32( val )
    if (val == -9999.0):
        if (REPORT):
            print('### WARNING: No data for f_s; setting to 0.')
            print()
        val = 0.0
    #--------------------------------------------
    # Initial units are "percent", in [0,100].
    # Set UNITS keyword to anything else to
    # convert to unitless values in [0,1].
    #--------------------------------------------
    if not(UNITS.lower() == 'percent'):
        f_s = (val / 100.0)  # convert percent to [0,1]
    return f_s

#   get_snow_precip_fraction()
#---------------------------------------------------------------------
def get_precip_timing_index( val_dict, f_s, REPORT=False ):

    #------------------------------------------------
    # Uses equation (6a) in Berghuijs et al. (2014)
    #  from Woods (2009).
    # delta_p should range from -1 to 1.
    #------------------------------------------------

    #----------------------------------------------------
    # Compute T_bar_star = T_bs = (T_avg - T0) / T_amp
    # T0 = rain-to-snow temperature threshold (1 deg_C)
    #-----------------------------------------------------
    # T_amp is "seasonal amplitude, degrees C", under
    # the assumption that temp. varies as sine function.
    # This code block ASSUMES July is hottest month and
    # January is coldest month, everywhere in the USA,
    # so that difference can be used to get amplitude.
    # Need to check this assumption.  ###########
    #-----------------------------------------------------
    # From Woods (2009), near equation (5):
    # T_bs >  1 => climate is above freezing all year.
    # T_bs =  0 => climate is below freezing 1/2 year.
    # T_bs < -1 => climate is below freezing all year.
    #-----------------------------------------------------    
    T0    = 1.0  # [deg C]
    T_avg = val_dict[ 'T_AVG_BASIN' ]  # [deg C]
    T_avg = np.float32( T_avg )
    T_rng = np.float32( val_dict[ 'JUL_JAN_TMP_DIFF_DEGC'] )
    T_amp = np.abs( T_rng ) / 2.0
    T_bs  = (T_avg - T0) / T_amp
     
    #----------------------------------------------------------
    # Solve equation (13) in Woods (2009) for delta_p to
    # get a function of f_s and T_bar_star (T_bs).
    #----------------------------------------------------------
    # NOTE!!.  There is a typo in Berghuijs et al. (2014), in
    # equation (6a), where it should be T**2 vs. T inside the
    # sqrt().  The correct equation is in Woods (2009),
    # equation (13). The plot of f_s for T_bs in [-1,1] makes
    # more sense for the correct equation.
    # See Mathematica notebook: Seasonal_Water_Balance.nb.
    #----------------------------------------------------------
    # Arcsin function requires argument in [-1,1].
    #----------------------------------------------------------    
    if ((-1 < T_bs) and (T_bs < 1)):
        numer = np.pi - (2 * f_s * np.pi) - (2 * np.arcsin(T_bs))
        ## denom = 2 * np.sqrt(1 - T_bs)     # Berghuis et al. (2014)
        denom = 2 * np.sqrt(1 - T_bs**2)  # Woods (2009)
        delta_p = (numer / denom)
        #---------------------------------------------
        # The range of the delta_p function
        # Using -0.999 and 0.999 here to distinguish
        # from -1 and 1 used in the "else" part.
        #---------------------------------------------
        delta_p = np.maximum(delta_p, -0.999)   #############
        # This gets used twice for GAGES-II.
        delta_p = np.minimum(delta_p, 0.999)      # max now is 1.029
        ### delta_p = np.minimum(delta_p, 1.0)      # max now is 1.029
    else:
        #----------------------------------------------------------
        # When f_s = 0, all precip falls as rain.
        # When f_s = 1, all precip falls as snow.
        #
        # If (T_bs = 1), then for ANY delta_p:
        #   (1) Eq (12) in Woods (2009) simplifies to  f_s = 0.0. 
        #   (2) Eq (6a) in Berghuijs (2014) also gives f_s = 0.0.
        # If (T_bs = -1), then:
        #   (1) Eq (12) in Woods (2009) simplifies to f_s = 1.0.
        #   (2) Eq (6a) in Berghuijs (2014) simplifies to:
        #          f_s = 1 - delta_p*sqrt(2)/pi
        # In both cases for Woods (2009, eq 12), the delta_p term
        # drops out and we can't solve for delta_p.
        #------------------------------------------------------------
        # However, when we can solve for delta_p, we find that
        # for any fixed value of f_s in (0,1):
        # (1)  There is a range of T-values, [T1,T2],
        #      with T1(f_s) > -1, and T2(f_s) < 1,
        #      so that delta_p is a monotonically decreasing fcn.
        #      of T, with delta_p(T1(f_s), f_s) = 1 and
        #      delta_p(T2(f_s), f_s) = -1.
        # (2)  As f_s approaches 0, both T1(f_s) and T2(f_s) ->  1.
        # (3)  As f_s approaches 1, both T1(f_s) and T2(f_s) -> -1.
        # (4)  For continuity, we need delta_p = 1 for T < T1, and
        #       delta_p = -1 for T > T2.
        # (5)  We can't analytically solve Woods eq 12 for
        #      T(f_s, delta_p) - it is transcendental.
        # (6)  If T=0, delta_p = (pi/2)[1 - 2*f_s], and then
        #      delta_p(f_s=0) = pi/2, and delta_p(f_s=1) = -pi/2.
        # (7)  delta_p(0, T) = arccos(T)/sqrt(1 - T^2)
        # (8)  delta_p(1, T) = [arcos(T) - pi]/sqrt(1 - T^2)
        # (9)  Limit[ delta_p(0, T), T ->  1 ] =  1.
        # (10) Limit[ delta_p(1, T), T -> -1 ] = -1.
        #
        # In comments near eqn (5) in Berghuijs (2014): 
        # (1) T is a dimensionless measure of mean temperature.
        # (2) delta_p = +1, for strongly summer-dominant precip;
        #     this also seems to match up with f_s = 0.
        # (3) delta_p = -1, for strongly winter-dominant precip
        #     this also seems to match up with f_s = 1.
        #------------------------------------------------------------
        if (T_bs <= -1):
            # Woods (2009, eq 12) implies f_s(-1, delta_p) = 1.
            # Limit[ delta_p(f_s=1, T), T -> -1 ] = -1.
            delta_p = -1.0
        if (T_bs >= 1):
            #----------------------------------------------------
            # Woods (2009), near eq (5), says:
            # "Climates with values of T􏰘 > 1 are above freezing
            # all year round, while climates with T􏰘 < 􏰛1 are
            # below freezing all year round.""
            #----------------------------------------------------
            # Note: This could also be happening because the
            # denominator (amplitude) is too small.  This would
            # be the case if January and June are not the
            # coldest and hottest months (on avg).
            #----------------------------------------------------
            # Woods (2009, eq 12) implies f_s(1, delta_p) = 0.
            # Limit[ delta_p(f_s=0, T), T -> 1 ] =  1.
            #----------------------------------------------------
            # This gets used a lot for GAGES-II.
            delta_p = 1.0
 
    if (REPORT):
        print('    T_avg =', T_avg, '[deg C]')
        print('    T_rng =', T_rng, '[deg C]')
        print('    T_bs  =', T_bs)
        print('    f_s   =', f_s)
        print('    d_p   =', delta_p)
               
    if (delta_p < -1) or (delta_p > 1):
        print('### ERROR:  delta_p is not in [-1,1].')

    return delta_p
    
#   get_precip_timing_index()
#---------------------------------------------------------------------
def get_aridity_index( val_dict, REPORT=False ):

    #---------------------------------------------------------
    # Equation (7) in Berghuijs et al. (2014) paper is used.
    # In that paper, units of P and PET are both [mm/day].
    # Only need for PET and P to have same units.
    # phi ranges from 0 to Infinity (theoretical), but in
    # the Berghuijs paper the max observed for 471 MOPEX
    # basins was about 5.3.
    #---------------------------------------------------------
    # In the GAGES-II spreadsheet: SB3_BCdescriptions.xlsx,
    # the units of mean annual precip (PPTAVG_BASIN) are
    # given as [cm].  However, the units for mean annual
    # potential evaporation (PET) are not given.  However,
    # based on information given in Over et al. (2018),
    # Table 1 (search doc for "aridity"), the max value they
    # observe, in region 15, is 5.531.  Since this agrees
    # with the max value value of PET after dividing by 10,    
    # the units of PET in GAGES-II must be [mm].
    #---------------------------------------------------------
    PET = val_dict[ 'PET' ]  # [cm?, not given]
    P   = val_dict[ 'PPTAVG_BASIN' ]  # [cm]
    PET = np.float32( PET ) # [mm]
    PET = PET / 10.0        # [mm -> cm] See Notes above.
    P   = np.float32( P )   # [cm]
    phi = (PET / P)         # [unitless]
    if (phi < 0):
        print('### ERROR: phi must be > 0.')
    return phi

#   get_aridity_index()
#---------------------------------------------------------------------
def get_swb_points( in_tsv_file=None, ORIGINAL=False ):

    #--------------------------------------------------------
    # Note: This is used to plot points on a CONUS map as
    #       circles where color is determined by SWB class.
    #--------------------------------------------------------
    if (in_tsv_file is None):
        in_tsv_file = 'SB3_all_regions_untransfBCs_sorted_SWB.tsv'

    delim       = '\t'  # TSV file
    new_dir     = dtu.get_new_data_dir( 'USGS_GAGES2_SB3')
    in_tsv_path = new_dir + in_tsv_file
    in_tsv_unit = open( in_tsv_path, 'r' )

    #--------------------------------------------------
    # Get headings & column numbers for some headings
    #--------------------------------------------------
    headings  = get_csv_headings( in_tsv_unit, delim='\t' )
    x_index   = headings.index('LONG_CENT')
    y_index   = headings.index('LAT_CENT')
    if (ORIGINAL):
        swb_index = headings.index('SWB1_CLASS')  # Can be "None"
    else:
        swb_index = headings.index('SWB2_CLASS')

    #------------------------------------  
    # Numpy arrays to store values into
    #------------------------------------
    n_basins = dtu.get_n_records( 'USGS_GAGES2_SB3' )
    x_g2 = np.zeros(n_basins, dtype='float32') 
    y_g2 = np.zeros(n_basins, dtype='float32')
    c_g2 = np.zeros(n_basins, dtype='<U4')
    
    k = 0            
    while (True):
        line = in_tsv_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        #----------------------------------------
        # Need to remove the newline char first
        #----------------------------------------
        line = line[:-1]
        values = line.split( delim )

        #-----------------------------------------
        # Save values for some columns to arrays
        #-----------------------------------------
        x_g2[k] = np.float32( values[ x_index ])
        y_g2[k] = np.float32( values[ y_index ])
        c_g2[k] = values[ swb_index ]
        k += 1
  
    in_tsv_unit.close()
    # print('Finished.') 
    return x_g2, y_g2, c_g2

#   get_swb_points()
#---------------------------------------------------------------------


    


