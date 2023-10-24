
# Copyright (c) 2023, Scott D. Peckham
#
# Oct 2023. Added merge_SB_info_for_regions(),
#           get_sb_region_headings(), compute_swb_classes(),
#           get_snow_precip_fraction(), get_precip_timing_index(),
#           get_aridity_index(), get_seasonal_water_balance_class(),
#           get_swb_class().     
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
#  >>> g2.compute_swb_classes(ORIGINAL=False)
#
#---------------------------------------------------------------------
#
#  get_basin_repo_dir()
#  get_gages2_data_dir()
#  create_tsv0()
#      
#  create_tsv_from_shapefile()
#  get_poly_feature_v0()
#  get_poly_feature()
#  get_gage_ID()
#  check_gage_ID()
#  check_overshoot()
#  check_bounding_box()
#  check_basin_poly_file()
#
#  merge_SB_info_for_regions()
#  get_sb_region_headings()
#  compute_swb_classes()
#     get_snow_precip_fraction()
#     get_precip_timing_index()
#     get_aridity_index()
#     get_seasonal_water_balance_class()
#     get_swb_class()         # See detailed Notes.
#
#---------------------------------------------------------------------

import numpy as np
from osgeo import ogr, osr
import json, sys, time

from topoflow.utils.ngen import shape_utils as su

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
    delim = ','  # CSV files
    gages2_dir = get_gages2_data_dir()
    gages2_SB_dir = gages2_dir + '__Selected_Basin_Info/'

    if (out_csv_file is None):
        out_csv_file = 'SB3_all_regions_untransfBCs.csv'
    out_csv_path = gages2_SB_dir + out_csv_file
    out_csv_unit = open( out_csv_path, 'w')

    region_list = ['1','2','3','4','5','6','7','8','9', '10L','10U',
                  '11','12','13','14','15','16','17','18']
    n_basins   = 0
  
    print('Working...')  
    for k in range(19):  
        region_str = region_list[k]   
        csv_file   = 'SB3_region' + region_str + '_untransfBCs.new.csv'
        csv_path   = gages2_SB_dir + csv_file
        csv_unit   = open( csv_path, 'r' )
        if (REPORT):
            print('Merging info from file:')
            print('   ' + csv_file)

        #-----------------------------------------------
        # Note: Not all files have same header !!
        # For missing columns, see the file:
        #    SB3_BCs_byRegion.xlsx.
        #-----------------------------------------------       
        headings = get_sb_region_headings( csv_unit )
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
def get_sb_region_headings( csv_unit, delim=',' ):

    header   = csv_unit.readline()
    headings = header.split( delim )
    for k in range( len(headings) ):
        headings[k] = headings[k].strip()
    return headings
   
#   get_sb_region_headings()
#---------------------------------------------------------------------
def compute_swb_classes( in_csv_file=None, ORIGINAL=False,
                         USE_B3=False, REPORT=False ):

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
    #       See the function:  get_swb_class() below.
    #--------------------------------------------------------------
    #  In filenames, "SB" stands for "Selected Basins", and
    #  "untransf" is short for "untransformed" values of vars.
    #--------------------------------------------------------------
    #  If USE_B3, a new class is used so that B1, B2, and B3
    #  "stack up" along phi-axis just as A1, A2, and A3 do.
    #--------------------------------------------------------------  
    if (in_csv_file is None):
        in_csv_file = 'SB3_all_regions_untransfBCs_sorted.csv'
    out_csv_file = in_csv_file.replace('.csv', '_SWB.csv')

    delim = ','  # CSV files
    gages2_dir = get_gages2_data_dir()
    gages2_SB_dir = gages2_dir + '__Selected_Basin_Info/'
    in_csv_path   = gages2_SB_dir + in_csv_file
    in_csv_unit   = open( in_csv_path, 'r' )
    #----------------------------------------------
#     out_csv_path  = gages2_SB_dir + out_csv_file
#     out_csv_unit  = open( out_csv_path, 'w' )
    
    headings = get_sb_region_headings( in_csv_unit, delim=',' )
    fs_min  = 1000.0
    fs_max  = -1000.0
    phi_min = 1000.0
    phi_max = -1000.0
    dp_min  = 1000.0
    dp_max  = -1000.0
    #-----------------------
    n_basins   = 0
    class_count = dict()
    if (USE_B3):
        classes = ['A1','A2','A3','B1','B2','B3',
                   'C1','C2','D1','D2','D3','None']
    else:
        classes = ['A1','A2','A3','B1','B2','C1','C2','D1','D2','D3','None']
    for name in classes:
        class_count[ name ] = 0

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
        f_s     = get_snow_precip_fraction( val_dict, REPORT=REPORT )
        delta_p = get_precip_timing_index( val_dict, f_s, REPORT=REPORT )
        phi     = get_aridity_index( val_dict, REPORT=REPORT )
        swb_class = get_swb_class( delta_p, f_s, phi,
                                   ORIGINAL=ORIGINAL, USE_B3=USE_B3)
        class_count[ swb_class ] += 1  # (swb_class may be 'None')
        n_basins += 1
        
        dp_min  = np.minimum( dp_min, delta_p )
        dp_max  = np.maximum( dp_max, delta_p )
        fs_min  = np.minimum( fs_min, f_s )
        fs_max  = np.maximum( fs_max, f_s )
        phi_min = np.minimum( phi_min, phi )
        phi_max = np.maximum( phi_max, phi )

        #----------------------------------------
        # For each heading in all_headings,
        # write corresponding value to new file
        #----------------------------------------
#         out_line = ''
#         for h in headings:
#             value = val_dict[h]
#             out_line += value + ','
#         out_line = out_line[:-1] + '\n'
#         out_csv_unit.write( out_line )

    #-------------------
    # Close both files
    #-------------------
    in_csv_unit.close()
#     out_csv_unit.close()
    #---------------------------
    print('del_p_min, del_p_max =', dp_min, ',', dp_max)
    print('fs_min, fs_max       =', fs_min, ',', fs_max)
    print('phi_min, phi_max     =', phi_min, ',', phi_max)
    print()
    print('Number of basins     =', n_basins)
    print('Number unclassified  =', class_count['None'] )
    print()
    for name in classes:
        print('   count[ ' + name + ' ] =', class_count[name] )
    print('Finished.')
    print()

#   compute_swb_classes()
#---------------------------------------------------------------------
def get_snow_precip_fraction( val_dict, REPORT=False ):

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
    f_s = (val / 100.0)  # in [0,1]
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
    # This code block assumes July is hottest month and
    # January is coldest month, everywhere in the USA,
    # so that difference can be used to get amplitude.
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
    #----------------------------------------------------------
    # Arcsin function requires argument in [-1,1].
    #----------------------------------------------------------    
    if ((-1 < T_bs) and (T_bs < 1)):
        numer = np.pi - (2 * f_s * np.pi) - (2 * np.arcsin(T_bs))
        ## denom = 2 * np.sqrt(1 - T_bs)     # Berghuis et al. (2014)
        denom = 2 * np.sqrt(1 - T_bs**2)  # Woods (2009)
        delta_p = (numer / denom)
        delta_p = np.maximum(delta_p, -0.999)   #############
        delta_p = np.minimum(delta_p, 1.0)      # max now is 1.029
    else:
        #---------------------------------------------------------
        # If (T_bs = 1), then eq (12) in Woods (2009) simplifies
        # to f_s = 0.0.  If (T_bs = -1), then it simplifies to
        # f_s = 1.0.  In both cases, the delta_p term drops out
        # and can't solve for delta_p.  But see comments near
        # equation (5) in Berghuijs et al. (2014).
        #---------------------------------------------------------
        if (T_bs <= -1):
            delta_p = -1.0
        if (T_bs >= 1):
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
def get_seasonal_water_balance_class( delta_p, f_s, phi):

    return get_swb_class( delta_p, f_s, phi)
  
#   get_seasonal_water_balance_class()
#---------------------------------------------------------------------
def get_swb_class( delta_p, f_s, phi, ORIGINAL=False,
                   USE_B3=True, REPORT=False):

    #-------------------------------------------------------------
    # Note: Set ORIGINAL=True and USE_B3=False to use the
    #       original 10 class boundaries of the Seasonal Water
    #       Balance (SWB) classification system as defined by
    #       Berghuijs et al. (2014).
    #       This leaves 1008 of the 1947 basins unclassified.
    #-------------------------------------------------------------
    #       Set ORIGINAL=False and USE_B3=False to use expanded
    #       class boundaries that span a much larger portion of
    #       the 3-parameter space, without any overlap.
    #       This results in all 1947 basins being classified.
    #       The 3D parameter space has the following parameters. 
    # 
    #       delta_p = precip_timing_index   in [-1,1]
    #       f_s     = snow_precip_fraction  in [0,1]
    #       phi     = aridity_index         in [0, infinity]
    #-------------------------------------------------------------
    #       Set ORIGINAL=False and USE_B3=True to use expanded
    #       class boundaries and to introduce a new class, B3,
    #       so that classes span the full 3-parameter space,
    #       without any overlap.  The class B3 "stacks above"
    #       the classes B1 and B2, just as A3 does for A1 & A2.
    #       So far, this extra class/region has not been needed.
    #-------------------------------------------------------------
    #       A Mathematica program has been written to visualize
    #       the set of "cuboid" regions in the parameter space.
    #-------------------------------------------------------------   
    d_p = delta_p  # synonym
    swb_class = 'None'

    #-------------------------------------------------------------
    # All A classes have -0.4 as the same max value of delta_p.
    # All C classes have 0.0 as the same min value of delta_p.
    # Classes D2 and D3 have -0.1 as the same min value of
    # delta_p, while the min value for D1 is -0.4.
    #-------------------------------------------------------------
    # This leaves a big void region in the 3D parameter space of
    # (delta_p, f_s, phi) that can be closed if we set these 5
    # delta_p values to be the same value in [-0.4, -0.1].
    #-------------------------------------------------------------
    if (ORIGINAL):
        A_dp_max  = -0.4
        # C_dp_min = 0.3   # (from fig. 7) 
        C_dp_min  = 0.0    # (from table 3)
        D1_dp_min = -0.4
        D2_dp_min = -0.1
        D3_dp_min = -0.1
    else:
        #---------------------------------------------
        # Helps: 1008 goes down to 953 unclassified.
        #---------------------------------------------
        mid_delta_p = -0.2
        A_dp_max    = mid_delta_p 
        C_dp_min    = mid_delta_p
        D1_dp_min   = -0.4
        # D1_dp_min   = mid_delta_p
        D2_dp_min   = mid_delta_p
        D3_dp_min   = mid_delta_p    
    
    #----------------------------------------------------------
    # The max value of delta_p for all B classes is the same.
    # In Berghuijs et al. (2014) figure 7, it is -0.4.
    # In Berghuijs et al. (2014) table 3, it is 0.0.       
    # This max value could be raised further to any value
    # <= 1.0 without intersecting any other classes.
    #-----------------------------------------------------------
    if (ORIGINAL):     
        # B_dp_max = -0.4   # (figure 7)
        B_dp_max = 0.0      # (table 3)
    else:
        #--------------------------------------------
        # Helps: 953 goes down to 951 unclassified.
        #-------------------------------------------- 
        B_dp_max = 1.0    # (max allowed value)
    
    #----------------------------------------------------------
    # The max value of delta_p for all the D classes could be
    # raised to 1.0 without intersecting any other classes.
    #----------------------------------------------------------
    if (ORIGINAL):
        D1_dp_max = 0.3
        D2_dp_max = 0.3
        D3_dp_max = 0.4
    else:
        #--------------------------------------------
        # Helps: 951 goes down to 333 unclassified.
        #-------------------------------------------- 
        D1_dp_max = 1.0
        D2_dp_max = 1.0
        D3_dp_max = 1.0
    
    #--------------------------------------------------------------
    # The min value of phi for classes A1, B1, D1, D2, & D3 could
    # all be lowered to 0 without intersecting any other classes.
    #--------------------------------------------------------------
    # A1_phi_min = 0.35   # Berghuijs et al. (2014), table 3
    #--------------------------------------------------------------
    if (ORIGINAL):
        A1_phi_min = 0.0   # Berghuijs et al. (2014), figure 7.
        B1_phi_min = 0.4
        D1_phi_min = 0.5
        D2_phi_min = 0.5
        D3_phi_min = 0.4
    else:
        #-------------------------------------------
        # Helps: 333 goes down to 18 unclassified.
        #------------------------------------------- 
        A1_phi_min = 0.0
        B1_phi_min = 0.0
        D1_phi_min = 0.0
        D2_phi_min = 0.0
        D3_phi_min = 0.0
    
    #-------------------------------------------------------------
    # In Berghuijs et al. (2014), the min value of phi for class
    # C1 is 0.9, in both figure 6 or table 3, and the max value
    # of phi for all D classes is 0.9.  There is a void in the
    # parameter space for phi values < 0.9 below the C1 box.
    # But if we were to lower C1_phi_min, the box for class C1
    # would overlap all the D-class boxes.
    #-------------------------------------------------------------
    # A better alternative to filling this void seems to be to
    # increase the maximum value of delta_p for all of the D
    # classes to the max possible value of 1.0.   
    #-------------------------------------------------------------
    C1_phi_min = 0.9
    D1_phi_max = 0.9
    D2_phi_max = 0.9
    D3_phi_max = 0.9
  
    #-----------------------------------------------------------     
    # In Berghuijs, the max value of f_s for class A3 is 0.45.
    # The max value of f_s for classes C1 and C2 is 0.25.
    # All of these max values of f_s could be raised to 1.0
    # without intersecting any other class boxes.
    #-----------------------------------------------------------
    # However, if we introduce a class B3 analogous to A3,
    # as done now, then we cannot raise A3_fs_max.
    # If we USE_B3 and don't raise A3_fs_max, previously
    #   unclassified points will be assigned to class B3.
    # If we don't USE_B3 and raise A3_fs_max, previously
    #   unclassified points wil be assigned to class A3. 
    #-----------------------------------------------------------
    # As classification code is written now, we assume that
    # C1 and C2 both use the same C_fs_max.
    #-----------------------------------------------------------
    if (ORIGINAL or USE_B3):
        A3_fs_max = 0.45
    else:
        A3_fs_max = 1.0
    #-----------------------
    if (ORIGINAL):
        C_fs_max = 0.25
    else:
        #-----------------------------------------
        # Helps: 
        #--------------------------------------------------
        # Can't set this higher than 0.45 or we intersect
        # the new bounds of B1, B2, etc.
        #--------------------------------------------------
        C_fs_max = 0.45

    #----------------------------------------------------------------    
    # In Berghuijs et al. (2014), the maximum value of phi for the
    # classes A3 and C2 are slightly different, namely, 5 and 5.3.
    # They can be safely raised to any higher value without
    # intersecting any other class box.  The max value of phi for
    # class B2 is 1.75, but it could also be raised to any higher
    # value.  Alternately, we could introduce a class B3 that has
    # phi values ranging from 1.75 to 5.3 (or some other phi_max),
    # and that is done here when USE_B3=True.
    #----------------------------------------------------------------
    if (ORIGINAL):
        # There is no B3 in original SWB method.
        A3_phi_max = 5.0
        C2_phi_max = 5.3
        B3_phi_max = 5.3
    else:    
        #----------------------------------------------------------
        # In the GAGES-II Selected Basins (SB3), after converting
        # units from mm to cm, the max value of phi is 5.5314.
        # The theoretical range is 0 to Infinity.
        # In the Berghuijs (2014) paper, phi max is 5.3, so here
        # we increase the max allowed to 5.6.
        #---------------------------------------------------------
        # Helps: 1 goes down to 0 unclassified.
        #----------------------------------------
        A3_phi_max = 5.6
        C2_phi_max = 5.6
        B3_phi_max = 5.6
   
    #--------------------------------
    # Check "A" classes for a match
    # "Precipitation out of phase"
    #--------------------------------
    if ((-1 < d_p) and (d_p <= A_dp_max)):
        if ((A1_phi_min < phi) and (phi <= 0.75)):
            if ((0 < f_s) and (f_s <= 0.45)):
                swb_class = 'A1'   
        elif ((0.75 < phi) and (phi <= 1.75)):
            if ((0 < f_s) and (f_s <= 0.45)):
                swb_class = 'A2' 
        elif ((1.75 < phi) and (phi <= A3_phi_max)):
            #-------------------------------------------
            # We can either introduce a class B3, or
            # we can raise A3_fs_max from 0.45 to 1.0
            #-------------------------------------------
            if ((0 < f_s) and (f_s <= A3_fs_max)):                
                swb_class = 'A3'

    #--------------------------------
    # Check "B" classes for a match
    # These are snow-dominated.
    #--------------------------------
    if ((-1 < d_p) and (d_p <= B_dp_max)):    # (from table 3)
        if ((0.45 < f_s) and (f_s <= 1)):
           if ((B1_phi_min < phi) and (phi <= 0.75)):
               swb_class = 'B1'
           elif ((0.75 < phi) and (phi <= 1.75)):
               swb_class = 'B2'
           elif ((1.75 < phi) and (phi <= B3_phi_max)):
               #-------------------------------------------
               # B3 is not one of the original 10 classes
               # but it is analogous to A3.
               #-------------------------------------------
               if (USE_B3):
                   swb_class = 'B3'
                
    #--------------------------------
    # Check "C" classes for a match
    # "Precipitation in phase"
    #--------------------------------
    if ((C_dp_min < d_p) and (d_p <= 1)):
        if ((0 <= f_s) and (f_s <= C_fs_max)):
            if ((C1_phi_min < phi) and (phi <= 1.5)):
                swb_class = 'C1'
            elif ((1.5 < phi) and (phi <= C2_phi_max)):
                swb_class = 'C2'

    #--------------------------------
    # Check "D" classes for a match
    # "Mild seasonality and humid"
    #--------------------------------
    if ((D1_dp_min < d_p) and (d_p <= D1_dp_max)):
        if (f_s <= 0):
            if ((D1_phi_min < phi) and (phi <= 0.9)):
                swb_class = 'D1'
    #--------------------------------------------------
    if ((D2_dp_min < d_p) and (d_p <= D2_dp_max)):
        if ((0 < f_s) and (f_s <= 0.2)):
            if ((D2_phi_min < phi) and (phi <= 0.9)):
                swb_class = 'D2'
    #--------------------------------------------------
    if ((D3_dp_min < d_p) and (d_p <= D3_dp_max)):
        if ((0.2 < f_s) and (f_s <= 0.45)):
            if ((D3_phi_min < phi) and (phi <= 0.9)):
                swb_class = 'D3'
    
    if (REPORT):
        if (swb_class == 'None'):
            print('### SORRY, No matching SWB class for:')
            print('    d_p =', d_p)
            print('    f_s =', f_s)
            print('    phi =', phi)
            print()
        else:
            print('SWB class =', swb_class)
            print()

    return swb_class

#   get_swb_class() 
#---------------------------------------------------------------------

   
    


