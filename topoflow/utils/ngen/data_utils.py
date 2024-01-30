
#---------------------------------------------------------------------
# Copyright (c) 2023-2024, Scott D. Peckham
# 
# Jan 2024. Moved general utility functions from other files
#           to here and modified the "db" key strings.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import data_utils as dtu
#  >>> dtu.compare_basin_ids( db_str2='CAMELS' )
#  >>> dtu.compare_basin_ids( db_str2='MOPEX' )
#  >>> dtu.compare_basin_ids( db_str2='USGS_GAGES2_ref' )
#  >>> dtu.compare_basin_ids( db_str2='USGS_GAGES2_SB3' )
#
#  >>> dtu.compare_basin_names( db_str2='CAMELS ')
#  >>> dtu.compare_basin_names( db_str2='MOPEX ')
#  >>> dtu.compare_basin_names( db_str2='USGS_GAGES2_ref ')
#
#---------------------------------------------------------------------
#
#  get_repo_dir()
#  get_data_dir()
#  get_new_data_dir()
#  get_new_tsv_filepath()
#  get_n_records()
#
#  get_heading_map()   ### New approach
#  get_id_column2()
#  get_name_column2()

#  get_id_column()
#  get_name_column()
#  get_lon_column()
#  get_lat_column()
#  get_minlon_column()
#
#  get_header_lines()
#  skip_header_lines()
#
#  get_basin_ids()
#  compare_basin_ids()
#  get_basin_names()
#  compare_basin_names()
#  
#  get_dataset_key_names()
#
#  get_state_code_map()
#  get_state_code()
#
#  export_csv_to_tsv()
#  convert_dms_to_dec_deg()
#
#---------------------------------------------------------------------
import numpy as np
import os, os.path

from topoflow.utils.ngen import collate_basins as cb

# import time

# import re  # for regular expressions, not used now
# from osgeo import ogr, osr
# import json, sys, time

#---------------------------------------------------------------------
def get_repo_dir():

    home_dir = os.path.expanduser('~/')
    repo_dir = home_dir + 'Dropbox/NOAA_NextGen/'
    repo_dir += '__Combined_River_Basin_Repo/'
    ### repo_dir = home_dir + 'Combined_River_Basin_Repo/'
    return repo_dir

#   get_repo_dir()
#---------------------------------------------------------------------
def get_data_dir( db_str='USGS_NWIS_Web' ):

    repo_dir = get_repo_dir()

    #---------------------------------------------    
    # Directories for each of the basin datasets
    #---------------------------------------------
    if (db_str == 'CAMELS'):
        data_dir = repo_dir + 'CAMELS/Data/'
    elif (db_str == 'MOPEX'):
        data_dir = repo_dir + 'MOPEX/Data/MOPEX431/'
    #----------------------------------------------------------
    elif (db_str == 'NOAA_RFC'):
        data_dir = repo_dir + 'NOAA_RFCs/Data/'
    elif (db_str == 'NOAA_HADS'):
        data_dir = repo_dir + 'NOAA_HADS/Data/'
    elif (db_str == 'NOAA_via_API'):
        data_dir = repo_dir + 'NOAA_via_API/Data/'   
    #----------------------------------------------------------
    elif (db_str == 'NSF_CZO'):
        data_dir = repo_dir + 'NSF_CZO/Data/'
    elif (db_str == 'NSF_LTER'):
        data_dir = repo_dir + 'NSF_LTER/Data/'
    elif (db_str == 'NSF_NEON'):
        data_dir = repo_dir + 'NSF_NEON/Data/'
    #----------------------------------------------------------
    elif (db_str == 'US_EPA'):
        data_dir = repo_dir + 'US_EPA/Data/'
    #----------------------------------------------------------
    elif (db_str == 'USDA_ARS'):
        data_dir = repo_dir + 'USDA_ARS/Data/'
    #----------------------------------------------------------            
    elif (db_str == 'USGS_FPS'):
        data_dir = repo_dir + 'USGS_FPS/Data/'
    elif (db_str == 'USGS_GAGES2_all'):
        data_dir = repo_dir + 'USGS_GAGES2/Data/All/'
    elif (db_str == 'USGS_GAGES2_nonref'):
        data_dir = repo_dir + 'USGS_GAGES2/Data/NonRef/'
    elif (db_str == 'USGS_GAGES2_ref'):
        data_dir = repo_dir + 'USGS_GAGES2/Data/Ref/'
    elif (db_str == 'USGS_GAGES2_SB3'):
        data_dir = repo_dir + 'USGS_GAGES2/Data/SB3/'
    #----------------------------------------------------------                
    elif (db_str == 'USGS_HLR'):
        # For vector data, add "hlrshape/".
        # For raster data, add "arctar00000/hlrus/""
        data_dir = repo_dir + 'USGS_HLR/Data/_USA/'
    elif (db_str == 'USGS_HCDN'):
        data_dir = repo_dir + 'USGS_HCDN/Data/CDROM_via_FTP/hcdn/'
    #---------------------------------------------------------------
    elif (db_str == 'USGS_NWIS_Web'):
        data_dir = repo_dir + 'USGS_NWIS_Web/Data/'
    elif (db_str == 'USGS_NWIS_Web_Old'):
        data_dir = repo_dir + 'USGS_NWIS_Web_Old/Data/'
    elif (db_str == 'USGS_NWIS_WQP1'):
        data_dir = repo_dir + 'USGS_NWIS_WQP/Data/'  #####
    elif (db_str == 'USGS_NWIS_WQP2'):
        data_dir = repo_dir + 'USGS_NWIS_WQP/Data/'  #####
    elif (db_str == 'USGS_NWIS_WQP3'):
        data_dir = repo_dir + 'USGS_NWIS_WQP/Data/'  #####
    #---------------------------------------------------------------
    elif (db_str == 'CENSUS_2014'):
        data_dir = repo_dir + '__Extras/USA_Shapefiles/CENSUS_2014/'
    elif (db_str == 'TIGER_2012'):
        data_dir = repo_dir + '__Extras/USA_Shapefiles/TIGER_2012/'
    else:
        print('SORRY: No match for:', db_str)
        data_dir = None

    return data_dir

#   get_data_dir()
#---------------------------------------------------------------------
def get_new_data_dir( db_str='USGS_NWIS_Web' ):

    repo_dir = get_repo_dir()

    #---------------------------------------------    
    # Directories for each of the basin datasets
    #---------------------------------------------
    if (db_str == 'CAMELS'):
        data_dir = repo_dir + 'CAMELS/_New/'
    elif (db_str == 'MOPEX'):
        data_dir = repo_dir + 'MOPEX/_New/'
    #----------------------------------------------------------
    elif (db_str == 'NOAA_RFC'):    
        data_dir = repo_dir + 'NOAA_RFCs/_New/'
    elif (db_str == 'NOAA_HADS'):    
        data_dir = repo_dir + 'NOAA_HADS/_New/'
    elif (db_str == 'NOAA_via_API'):    
        data_dir = repo_dir + 'NOAA_via_API/_New/'
    #----------------------------------------------------------
    elif (db_str == 'NSF_CZO'):
        data_dir = repo_dir + 'NSF_CZO/_New/'
    elif (db_str == 'NSF_LTER'):
        data_dir = repo_dir + 'NSF_LTER/_New/'
    elif (db_str == 'NSF_NEON'):
        data_dir = repo_dir + 'NSF_NEON/_New/'
    #----------------------------------------------------------
    elif (db_str == 'US_EPA'):
        data_dir = repo_dir + 'US_EPA/_New/'
    #----------------------------------------------------------
    elif (db_str == 'USDA_ARS'):
        data_dir = repo_dir + 'USDA_ARS/_New/'
    #----------------------------------------------------------            
    elif (db_str == 'USGS_FPS'):
        data_dir = repo_dir + 'USGS_FPS/_New/'
    #---------------------------------------------------------- 
    elif (db_str == 'USGS_GAGES2_all'):
        data_dir = repo_dir + 'USGS_GAGES2/_New/All/'
    elif (db_str == 'USGS_GAGES2_nonref'):
        data_dir = repo_dir + 'USGS_GAGES2/_New/NonRef/'
    elif (db_str == 'USGS_GAGES2_ref'):
        data_dir = repo_dir + 'USGS_GAGES2/_New/Ref/'
    elif (db_str == 'USGS_GAGES2_SB3'):
        data_dir = repo_dir + 'USGS_GAGES2/_New/SB3/'
    #---------------------------------------------------------- 
    elif (db_str == 'USGS_HLR'):
        data_dir = repo_dir + 'USGS_HLR/_New/'
    elif (db_str == 'USGS_HCDN'):
        data_dir = repo_dir + 'USGS_HCDN/_New/'
    #----------------------------------------------------------  
    elif (db_str == 'USGS_NWIS_Web'):
        data_dir = repo_dir + 'USGS_NWIS_Web/_New/'
    elif (db_str == 'USGS_NWIS_Web_Old'):
        data_dir = repo_dir + 'USGS_NWIS_Web_Old/_New/'
    elif (db_str == 'USGS_NWIS_WQP1'):
        data_dir = repo_dir + 'USGS_NWIS_WQP/_New/'
    elif (db_str == 'USGS_NWIS_WQP2'):
        data_dir = repo_dir + 'USGS_NWIS_WQP/_New/'
    elif (db_str == 'USGS_NWIS_WQP3'):
        data_dir = repo_dir + 'USGS_NWIS_WQP/_New/'
    else:
        print('SORRY: No match for:', db_str)
        data_dir = None

    return data_dir

#   get_new_data_dir()
#---------------------------------------------------------------------
def get_new_tsv_filepath( db_str='USGS_NWIS_Web' ):

    data_dir = get_data_dir( db_str )
    new_dir  = get_new_data_dir( db_str )

    if (db_str == 'CAMELS'):
        file_path = new_dir + 'new_camels.tsv'
    elif (db_str == 'MOPEX'):
        file_path = new_dir + 'new_mopex431_sorted.tsv'
    #-----------------------------------------------------------------
    elif (db_str == 'NOAA_RFC'):
        ##########################################
        # This currently doesn't have USGS IDs.
        ##########################################
        file_path = new_dir + 'new_rfc_info.tsv'
        ## file_path = new_dir + 'new_rfc_info2.tsv'
    elif (db_str == 'NOAA_HADS'):
        file_path = new_dir + 'new_hads_info_sorted.tsv'
    elif (db_str == 'NOAA_via_API'):
        file_path = new_dir + 'USGS_NWIS_Web_Site_Info_via_API.tsv'
        ### file_path = new_dir + 'USGS_NWIS_WQP3_Site_Info_via_API.tsv'
    #-----------------------------------------------------------------
    elif (db_str == 'NSF_CZO'):
        file_path = new_dir + 'czo_basin_info_sheet1.tsv'  #####
    elif (db_str == 'NSF_LTER'):
        file_path = new_dir + 'new_lter_info.tsv'
    elif (db_str == 'NSF_NEON'):
        file_path = new_dir + 'NEON_Field_Site_Metadata_20231026.tsv'
    #-----------------------------------------------------------------
    elif (db_str == 'USDA_ARS'):
        file_path = new_dir + 'ars_basin_info.tsv'
    elif (db_str == 'US_EPA'):
        file_path = data_dir + 'EPA_WQX_sites.tsv'      ######### data_dir   
    #-----------------------------------------------------------------
    elif (db_str == 'USGS_FPS'):
        file_path = new_dir + 'FPS_basin_info.tsv'
    elif (db_str == 'USGS_GAGES2_all'):
        file_path = new_dir + 'new_gages2_all.tsv'
    elif (db_str == 'USGS_GAGES2_ref'):
        file_path = new_dir + 'new_gages2_ref.tsv'
    elif (db_str == 'USGS_GAGES2_SB3'):
        file_path = new_dir + 'SB3_all_regions_untransfBCs_sorted_SWB.tsv'   #####
    elif (db_str == 'USGS_HCDN'):
        file_path = new_dir + 'hcdn_stations.tsv'
    elif (db_str == 'USGS_HLR'):
        #------------------------------------------------------
        # Note: Only use ones where an outlet was determined.
        #------------------------------------------------------
        file_path = new_dir + 'new_hlr_na_with_outlets.tsv'
        ## file_path = new_dir + 'new_hlr_na_no_outlets.tsv'
    #-----------------------------------------------------------------
    # These currently use "data_dir" vs. "new_dir".
    # but apparently most don't measure discharge.
    #----------------------------------------------------------------- 
    elif (db_str == 'USGS_NWIS_Web'):
        file_path = data_dir + 'NWIS_Stream_Site_Daily_Data.tsv'
    elif (db_str == 'USGS_NWIS_Web_Old'):
        file_path = data_dir + 'USGS_stream_gauge_data.tsv'     ############
    #-----------------------------------------------------------------
    # The one for 'USGS_NWIS3' has 145375 "stream-type" records
    # but apparently most don't measure discharge.
    #-----------------------------------------------------------------          
    elif (db_str == 'USGS_NWIS_WQP1'):
        file_path = data_dir + 'WQP_USGS_Stream_Site_Daily_Q_Info.tsv'       
    elif (db_str == 'USGS_NWIS_WQP2'):
        file_path = data_dir + 'WQP_USGS_Stream_Site_Inst_Q_Info.tsv'
    elif (db_str == 'USGS_NWIS_WQP3'):
        file_path = data_dir + 'WQP_USGS_Stream_Site_Info.tsv'
    else:
        print('SORRY: No match for:', db_str)
        file_path = None
       
    return file_path
 
#   get_new_tsv_filepath()
#---------------------------------------------------------------------
def get_n_records( db_str ):

    if (db_str == 'CAMELS'):
        n_records = 671
    elif (db_str == 'MOPEX'):
        n_records = 431
    #--------------------------------------------------------
    elif (db_str == 'NOAA_RFC'):
        n_records = 7759
    elif (db_str == 'NOAA_HADS'):
        n_records = 10512   # not just streams
    elif (db_str == 'NOAA_via_API'):
        n_records = 7136
    #--------------------------------------------------------
    elif (db_str == 'NSF_CZO'):
        n_records = 15  ### w/ currently available data. 19 in TSV file.        
    elif (db_str == 'NSF_LTER'):
        n_records = 30
    elif (db_str == 'NSF_NEON'):
        n_records = 81
    #--------------------------------------------------------
    elif (db_str == 'USDA_ARS'):
        n_records = 771
    #--------------------------------------------------------
    elif (db_str == 'USGS_FPS'):
        n_records = 4756
    elif (db_str == 'USGS_GAGES2_all'):
        n_records = 9322
    elif (db_str == 'USGS_GAGES2_ref'):
        n_records = 2057
    elif (db_str == 'USGS_GAGES2_SB3'):
        n_records = 1947
    elif (db_str == 'USGS_HCDN'):
        n_records = 1703    ### should be 1639 ??
    elif (db_str == 'USGS_HLR'):
        n_records = 43391    # Total, correct number
        # n_records = 47479  # w/ repeated VALUEs
    elif (db_str == 'USGS_HLR_outlets'):
        n_records = 9539    # Have outlet info in new_hlr_na_all.csv
    elif (db_str == 'USGS_NHDplus_v1'):   # Mostly EPA
        n_records = 19031
    elif (db_str == 'USGS_NHDplus_v2'):   # Mostly EPA
        n_records = 28163
    #-----------------------------------------------
    # Not sure why this had 24520 vs. 27915.
    #-----------------------------------------------
#     elif (db_str == 'USGS_NWIS_Web_Old'):
#         n_records = 24520
    #-----------------------------------------------
    elif (db_str == 'USGS_NWIS_Web'):
        n_records = 27915             # Daily Q (via NWIS Web)
    elif (db_str == 'USGS_NWIS_WQP1'):
        n_records = 14707             # Daily Q (via WQP)
    elif (db_str == 'USGS_NWIS_WQP2'):
        n_records = 46760             # Inst. Q (via WQP)
    elif (db_str == 'USGS_NWIS_WQP3'):
        n_records = 145375            # Any stream-type msmt.
    else:
        print('SORRY: No match for:', db_str)
        n_records = None

    return n_records

#   get_n_records()
#---------------------------------------------------------------------
def get_heading_map():

    #------------------------------------------------------------
    # Note: Here "id" refers to a USGS Site ID (8 to 15 chars)
    #       Also, lon and lat are for the basin outlet or site.
    #------------------------------------------------------------
    h_map = dict()
    h_map['CAMELS'] = {
        'id':'gauge_id', 'name':'gauge_name',
        'lon':'gauge_lon', 'lat':'gauge_lat', 'minlon':'MINLON'}
    h_map['MOPEX'] = {
        'id':'SiteCode', 'name':'SiteName',
        'lon':'Longitude', 'lat':'Latitude', 'minlon':'MINLON'}
    #-------------------------------------------------------------
    h_map['NOAA_HADS'] = {
        'id':'USGS_ID', 'name':'USGS_name',
        'lon':'Longitude', 'lat':'Latitude', 'minlon':None}
    h_map['NOAA_RFC'] = {  ### for ba12my15.csv
        'id':None, 'name':'NAME',
        'lon':'LON', 'lat':'LAT', 'minlon':'MINLON'}  ########get
#     h_map['NOAA_RFC'] = {  ### for new_rfc_info.tsv
#         'id':None, 'name':None,
#         'lon':None, 'lat':None, 'minlon':'MINLON'} 
    h_map['NOAA_via_API'] = {
        'id':'usgs_id', 'name':'usgs_name',
        'lon':'longitude', 'lat':'latitude', 'minlon':None}         
    #-------------------------------------------------------------
    h_map['USDA_ARS'] = {
        'id':None, 'name':'Watershed Name',
        'lon':'Longitude (dec deg)', 'lat':'Latitude (dec deg)',
        'minlon':None}
    #-------------------------------------------------------------
    h_map['USGS_FPS'] = {
        'id':'SiteNumber', 'name':'SiteName',
        'lon':None, 'lat':None, 'minlon':None}
    h_map['USGS_HCDN'] = {
        'id':'usgs_id', 'name':None,
        'lon':None, 'lat':None, 'minlon':None}
    h_map['USGS_HLR'] = {
        'id':None, 'name':None,  # but have Closest_USGS_ID
        'lon':'OUTLON', 'lat':'OUTLAT', 'minlon':'MINLON'}
    #-------------------------------------------------------------
    h_map['USGS_GAGES2_all'] = {
        'id':'STAID', 'name':'STANAME',
        'lon':'LNG_GAGE', 'lat':'LAT_GAGE', 'minlon':'MINLON'}
    h_map['USGS_GAGES2_ref'] = {
        'id':'STAID', 'name':'STANAME',
        'lon':'LNG_GAGE', 'lat':'LAT_GAGE', 'minlon':'MINLON'}  
    h_map['USGS_GAGES2_SB3'] = {
        'id':'StnID', 'name':None,
        'lon':None, 'lat':None, 'minlon':None}      
    #-------------------------------------------------------------
    h_map['USGS_NWIS_Web'] = {
        'id':'site_no', 'name':'site_nm',
        'lon':'dec_long_va', 'lat':'dec_lat_va', 'minlon':None}      
    #-------------------------------------------------------------
    h_map['USGS_NWIS_WQP1'] = {
        'id':  'MonitoringLocationIdentifier',
        'name':'MonitoringLocationName',
        'lon': 'LongitudeMeasure',
        'lat': 'LatitudeMeasure', 'minlon':None}
    h_map['USGS_NWIS_WQP2'] = {
        'id':  'MonitoringLocationIdentifier',
        'name':'MonitoringLocationName',
        'lon': 'LongitudeMeasure',
        'lat': 'LatitudeMeasure', 'minlon':None}
    h_map['USGS_NWIS_WQP3'] = {
        'id':  'MonitoringLocationIdentifier',
        'name':'MonitoringLocationName',
        'lon': 'LongitudeMeasure',
        'lat': 'LatitudeMeasure', 'minlon':None}            
    #-------------------------------------------------------------
    h_map['NSF_CZO'] = {
        ##### Will be adding USGS id soon
        'id':None, 'name':'Catchment Name',
        'lon':'Longitude (dd)', 'lat':'Latitude (dd)',
        'minlon':None}
    h_map['NSF_LTER'] = {
        'id':None, 'name':'Site Name',
        'lon':'Longitude', 'lat':'Latitude', 'minlon':None}
    h_map['NSF_NEON'] = {
        # There is a "field_site_id" but its not a USGS ID.
        'id':None, 'name':'field_site_name',
        'lon':'field_longitude', 'lat':'field_latitude',
        'minlon':None}
    #-------------------------------------------------------------
    return h_map
    
#   get_heading_map()
#---------------------------------------------------------------------
def get_id_column2( db_str ):

    h_map        = get_heading_map()
    all_headings = cb.get_headings( key=db_str )
    id_heading   = h_map[ db_str ]['id']
    if (id_heading is not None):
        col = all_headings.index( id_heading )
    else:
        print('Sorry, this dataset does not have USGS ID.')
        col = -1

    return col

#   get_id_column2()
#---------------------------------------------------------------------
def get_name_column2( db_str ):

    h_map        = get_heading_map()
    all_headings = cb.get_headings( key=db_str )
    id_heading   = h_map[ db_str ]['name']
    if (id_heading is not None):
        col = all_headings.index( id_heading )
    else:
        print('Sorry, this dataset does not have site name.')
        col = -1
        
    return col

#   get_name_column2()
#---------------------------------------------------------------------
def get_id_column( db_str ):

    #--------------------------------------------------------
    # ID Column Heading for each dataset:
    #
    # CAMELS          = "gauge_id"
    # MOPEX           = "SiteCode"
    #--------------------------------------------------------
    # NOAA_RFC        = File doesn't have USGS site ID
    # NOAA_HADS       = "USGS_ID"
    # NOAA_via_API    = "usgs_id"
    #--------------------------------------------------------
    # NSF_CZO         = ?????????????
    # NSF_LTER        = "Site Acronym" (e.g. "AND")
    # NSF_NEON        = "field_site_id" (e.g. "ABBY")
    #--------------------------------------------------------
    # USDA_ARS        = "Watershed ID"
    #--------------------------------------------------------
    # USGS_FPS        = "SiteNumber" (prepend 0 if 7 digits)
    # USGS_GAGES2_all = "STAID"
    # USGS_GAGES2_ref = "STAID"
    # USGS_GAGES2_SB3 = "StnID"
    # USGS_HCDN       = "usgs_id" (in TSV file)
    # USGS_HLR        = "VALUE"
    # USGS_NWIS_Web   = "site_no"
    # USGS_NWIS_WQP1  = "MonitoringLocationIdentifier"
    # USGS_NWIS_WQP2  = "MonitoringLocationIdentifier"
    # USGS_NWIS_WQP3  = "MonitoringLocationIdentifier"
    #--------------------------------------------------------
    id_col = None
    list1 = ['CAMELS', 'MOPEX', 'NOAA_HADS', 'NOAA_via_API',
             'NSF_LTER', 'USDA_ARS',
             'USGS_FPS', 'USGS_GAGES2_all', 'USGS_GAGES2_ref',
             'USGS_GAGES2_SB3', 'USGS_HCDN', 'USGS_HLR']
    list2 = ['NSF_NEON', 'USGS_NWIS_Web']
    list3 = ['USGS_NWIS_WQP1', 'USGS_NWIS_WQP2', 'USGS_NWIS_WQP3']

    if (db_str in list1):
        id_col = 0
    elif (db_str in list2):
        id_col = 1
    elif (db_str in list3):
        id_col = 2
    else:
        #### CZO ####
        print('SORRY: No match for:', db_str)
        print('       Returning id_col = 0.')
        id_col = 0

    return id_col

#   get_id_column()
#---------------------------------------------------------------------
def get_name_column( db_str ):

    #--------------------------------------------------------
    # Name Column Heading for each dataset:
    #
    # CAMELS          = "gauge_name"
    # MOPEX           = "SiteName"
    #--------------------------------------------------------
    # NOAA_RFC        = 
    # NOAA_HADS       = "USGS_name"
    # NOAA_via_API    = "usgs_name"
    #--------------------------------------------------------
    # NSF_LTER        = "Site Name" (but not a basin name)
    # NSF_NEON        = "field_site_name"
    #--------------------------------------------------------
    # USDA_ARS        = "Watershed name" (vs. Location, etc.)
    #--------------------------------------------------------
    # USGS_FPS        = "SiteName"
    # USGS_GAGES2_all = "STANAME"
    # USGS_GAGES2_ref = "STANAME"
    # USGS_GAGES2_SB3 = No name given #####################
    # USGS_HCDN       = No header in stations.dat
    # USGS_HLR        = No name given
    # USGS_NWIS_Web   = "station_nm"
    # USGS_NWIS_WQP1  = "MonitoringLocationName"
    # USGS_NWIS_WQP2  = "MonitoringLocationName"
    # USGS_NWIS_WQP3  = "MonitoringLocationName"
    #--------------------------------------------------------
    name_col = None
    list1 = ['NOAA_HADS', 'NSF_LTER', 'USGS_FPS',
             'USGS_GAGES2_all', 'USGS_GAGES2_ref', 'USGS_HCDN'  ]
             ########## 'USGS_GAGES2_SB3' ]  # No name given
    list2 = ['CAMELS', 'MOPEX', 'NSF_NEON', 'USGS_NWIS_Web']
    list3 = ['USGS_NWIS_WQP1', 'USGS_NWIS_WQP2', 'USGS_NWIS_WQP3']
    list4 = ['USDA_ARS']
    list5 = ['NOAA_via_API']
    if (db_str in list1):
        name_col = 1
    elif (db_str in list2):
        name_col = 2
    elif (db_str in list3):
        name_col = 3
    elif (db_str in list4):
        name_col = 4  # Watershed name
    elif (db_str in list5):
        name_col = 7
    else:
        #### NOAA_HADS, NSF_CZO, and USGS_HLR ####
        print('SORRY: No match for:', db_str)
        print('       Returning name_col = 0.')
        name_col = 0
       
    return name_col

#   get_name_column()
#---------------------------------------------------------------------
def get_lon_column( db_str ):

    #-------------------------------------------------------
    # Name Longitude Heading for each dataset:
    #
    # CAMELS          = 
    # MOPEX           =
    #-------------------------------------------------------
    # NOAA_RFC        =
    # NOAA_HADS       = Longitude, 3
    # NOAA_via_API    = 
    #-------------------------------------------------------
    # NSF_CZO         =        
    # NSF_LTER        = 
    # NSF_NEON        = 
    #-------------------------------------------------------
    # USDA_ARS        =
    # US_EPA          =
    #-------------------------------------------------------
    # USGS_FPS        = 
    # USGS_GAGES2_all = "LNG_GAGE" 7
    # USGS_GAGES2_ref = "LNG_GAGE" 7
    # USGS_HCDN       =   
    # USGS_HLR        =   
    # USGS_NWIS_Web   = "dec_long_va" 5
    # USGS_NWIS_WQP1  = "LongitudeMeasure" 12
    # USGS_NWIS_WQP2  = "LongitudeMeasure" 12
    # USGS_NWIS_WQP3  = "LongitudeMeasure" 12
    #-------------------------------------------------------
    # NOTE:  For GAGES2_SB3, outlet lat/lon are not given,
    #        but we can find them using site ID.
    #        However, the centroid lat/lon are provided.
    #-------------------------------------------------------
    list1 = ['USGS_GAGES2_all', 'USGS_GAGES2_ref']
    list2 = ['USGS_NWIS_WQP1', 'USGS_NWIS_WQP2', 'USGS_NWIS_WQP3']   
    if (db_str in list1):
        lon_col = 7
    elif (db_str == 'USGS_NWIS_Web'):
        lon_col = 5
    elif (db_str == 'NOAA_HADS'):
        lon_col = 3
    elif (db_str == 'USGS_HCDN'):
        lon_col = 10
    elif (db_str in list2):
        lon_col = 12
    else:
        print('SORRY: No match for:', db_str)
        lon_col = None

    return lon_col

#   get_lon_column()
#---------------------------------------------------------------------
def get_lat_column( db_str ):

    #-------------------------------------------------------
    # Name Longitude Heading for each dataset:
    #
    # CAMELS          = 
    # MOPEX           = 
    #-------------------------------------------------------
    # NOAA_RFC        =
    # NOAA_HADS       = Latitude, 4
    # NOAA_via_API    = 
    #-------------------------------------------------------
    # NSF_CZO         =        
    # NSF_LTER        = 
    # NSF_NEON        = 
    #-------------------------------------------------------
    # USDA_ARS        =
    # US_EPA          =
    #-------------------------------------------------------
    # USGS_FPS        = 
    # USGS_GAGES2_all = "LAT_GAGE" 6
    # USGS_GAGES2_ref = "LAT_GAGE" 6
    # USGS_HCDN       =   
    # USGS_HLR        =   
    # USGS_NWIS_Web   = "dec_lat_va" 4
    # USGS_NWIS_WQP1  = "LatitudeMeasure" 11
    # USGS_NWIS_WQP2  = "LatitudeMeasure" 11
    # USGS_NWIS_WQP3  = "LatitudeMeasure" 11
    #-------------------------------------------------------
    # NOTE:  For GAGES2_SB3, outlet lat/lon are not given,
    #        but we can find them using site ID.
    #        However, the centroid lat/lon are provided.
    #-------------------------------------------------------
    list1 = ['USGS_GAGES2_all', 'USGS_GAGES2_ref']
    list2 = ['USGS_NWIS_WQP1', 'USGS_NWIS_WQP2', 'USGS_NWIS_WQP3']
    list3 = ['USGS_NWIS_Web', 'NOAA_HADS']   
    if (db_str in list1):
        lat_col = 6
    elif (db_str in list3):
        lat_col = 4
    elif (db_str == 'USGS_HCDN'):
        lat_col = 9
    elif (db_str in list2):
        lat_col = 11
    else:
        print('SORRY: No match for:', db_str)
        lat_col = None

    return lat_col

#   get_lat_column()
#---------------------------------------------------------------------
def get_minlon_column( db_str, SILENT=False ):

    list1 = ['USGS_GAGES2_all', 'USGS_GAGES2_ref']
    if (db_str == 'CAMELS'):
        col = 3
    elif (db_str == 'MOPEX'):
        col = 5
    elif (db_str in list1):
        col = 8     
    elif (db_str == 'USGS_HLR'):
        col = 5
    else:
        if not(SILENT):
            print('SORRY: No match for:', db_str)
            print('       Returning col = None.')
        col = None

    return col

#   get_minlon_column()
#---------------------------------------------------------------------
# def get_elev_column( db_str, SILENT=False ):
# 
#     # CAMELS only has elev_mean
#     # MOPEX has None in the "Elevation_" column.
#     
# #   get_elev_column()
# #---------------------------------------------------------------------
# def get_area_column( db_str, SILENT=False ):
# 
#     # CAMELS has "area_gages2" in col 11
#     # MOPEX has "DrainageAr" in col 20
#     
# #   get_area_column()
#---------------------------------------------------------------------
def get_header_lines( db_str='USGS_NWIS_Web'):

    if (db_str == 'USGS_NWIS_Web'):
        #-------------------------------------------------
        # This depends on how many attributes are saved.
        # Get 77 for all attributes on NWIS Web page.
        #-------------------------------------------------
        nh_lines = 77   # all attributes
    elif (db_str == 'USGS_NWIS_Web_Old'):
        nh_lines = 37   # some attributes
    elif (db_str == 'USGS_FPS'):
        nh_lines = 2
    else:
        nh_lines = 1
   
    return nh_lines
    
#   get_header_lines()
#---------------------------------------------------------------------
def skip_header_lines( file_unit, key='USGS_NWIS_Web'):

    nh_lines = get_header_lines( key )
    for k in range( nh_lines ):
        line = file_unit.readline()
    
#   skip_header_lines()
#---------------------------------------------------------------------
def get_basin_ids( db_str='USGS_NWIS_Web' ):

    file_path = get_new_tsv_filepath( db_str )
    file_unit = open( file_path, 'r' )
 
    #--------------------
    # Get the delimiter
    #--------------------
    if (file_path.endswith('.csv')):
        delim = ','
    elif (file_path.endswith('.tsv')):
        delim = '\t'
    else:
        delim = ';'

    #-----------------------------
    # Skip over the header lines
    #-----------------------------
    nh_lines = get_header_lines( db_str ) 
    for k in range(nh_lines):
        file_unit.readline()
        
    id_col  = get_id_column( db_str )    
    id_list = list()
    wqp_list = ['USGS_NWIS_WQP1', 'USGS_NWIS_WQP2', 'USGS_NWIS_WQP3']

    while (True):
        line = file_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        values = line.split( delim )     
        id = values[ id_col ].strip()
        if (db_str in wqp_list):
            id = id.replace('USGS-', '')
        if (len(id) == 7):
            id = '0' + id   ########## 2023-11-22. Is this ever needed?
        id_list.append( id )     

    file_unit.close()
    return id_list
      
#   get_basin_ids()
#---------------------------------------------------------------------
def compare_basin_ids( db_str1='USGS_NWIS_Web',
                       db_str2='USGS_GAGES2_all',
                       PRINT_DIFF1_IDs=False, PRINT_DIFF2_IDs=False):

    #-------------------------------------------------------------
    # Notes:
    # db_str1='USGS_GAGES2_SB3', db_str2='MOPEX' shows that:
    # 343 MOPEX basins are not in USGS_GAGES2_SB3.
    #-------------------------------------------------------------
    # db_str1='USGS_GAGES2_ref', db_str2='MOPEX' shows that:
    # 343 MOPEX basins are not in USGS_GAGES2_ref.
    #-------------------------------------------------------------
    # db_str1='USGS_GAGES2_all', db_str2='MOPEX' shows that:
    # only 7 MOPEX basins are not in USGS_GAGES2_all.
    #-------------------------------------------------------------
    # db_str1='USGS_NWIS_Web', db_str2='MOPEX' shows that:
    # 1 MOPEX basins is not in USGS_NWIS_Web.
    #-------------------------------------------------------------
    # db_str1='USGS_GAGES2_ref', db_str2='CAMELS' shows that:
    # 0 CAMELS basins are not in USGS_GAGES2_ref.
    #-------------------------------------------------------------
    # db_str1='USGS_NWIS_WQP3', db_str2='USGS_GAGES2_all' shows:
    # 18 USGS_GAGES2_all basins are not in USGS_NWIS_WQP3.
    #-------------------------------------------------------------      
    id_list1 = get_basin_ids( db_str1 )
    id_list2 = get_basin_ids( db_str2 )
    # Apply "set()" here removes possible duplicates (ARS has them)
    print(db_str1 + ' basins: n =', len(set(id_list1)))
    print(db_str2 + ' basins: n =', len(set(id_list2)))
    id_list3 = list( set(id_list1) & set(id_list2) )
    print('Basins in both', db_str1, 'and', db_str2 + ': ', 'n =', len(id_list3))
    #print()
#     print('Dataset 1 =', db_str1, '(', len(id_list1), 'basins)' )
#     print('Dataset 2 =', db_str2, '(', len(id_list2), 'basins)' )
        
    #-----------------------------
    # Compare the two name lists
    #----------------------------- 
    diff1 = list( set(id_list1) - set(id_list2))
    diff1.sort()
    print('Basins in', db_str1, 'but not', db_str2 + ': ', 'n =', len(diff1))
    if (PRINT_DIFF1_IDs):
        print('Basin IDs in dataset 1 not in dataset 2:')
        print(diff1)
        print()

    diff2 = list( set(id_list2) - set(id_list1))
    diff2.sort()
    print('Basins in', db_str2, 'but not', db_str1 + ': ', 'n =', len(diff2))
    if (PRINT_DIFF2_IDs):
        print('Basin IDs in dataset 2 not in dataset 1:')
        print(diff2)
        print('len(diff2) =', len(diff2))
        print()
    ## return diff2

#   compare_basin_ids()
#---------------------------------------------------------------------
def get_basin_names( db_str='USGS_NWIS_Web', UPPER=True,
                     REPLACE_PUNCTUATION=True,
                     REPLACE_ABBREVIATIONS=True,
                     FIX_SPELLINGS=True,
                     REPLACE_STATE_LETTERS=True,
                     REPLACE_CITY_LETTERS=True,
                     FIX_MISSING_STATE=True,
                     FIX_DIRECTIONS=True):

    #---------------------------------------
    # What about:
    # SUGAR C TR A LK FRSME NR ARCADIA 
    #---------------------------------------
    file_path = get_new_tsv_filepath( db_str )
    file_unit = open( file_path, 'r' )
 
    #--------------------
    # Get the delimiter
    #--------------------
    if (file_path.endswith('.csv')):
        delim = ','
    elif (file_path.endswith('.tsv')):
        delim = '\t'
    else:
        delim = ';'

    #-----------------------------
    # Skip over the header lines
    #-----------------------------
    nh_lines = get_header_lines( db_str ) 
    for k in range(nh_lines):
        file_unit.readline()
        
    name_col  = get_name_column( db_str )    
    name_list = list()

    while (True):
        line = file_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        values = line.split( delim )
        name = values[ name_col ].strip()
        name = usgs.expand_usgs_site_name( name )   #####       
        name_list.append( name )     

    file_unit.close()
    return name_list
      
#   get_basin_names()
#---------------------------------------------------------------------
def compare_basin_names( db_str1='USGS_NWIS_Web',
                         db_str2='USGS_GAGES2_all'):

    name_list1 = get_basin_names( db_str1 )
    name_list2 = get_basin_names( db_str2 )
    print('Dataset 1 =', db_str1)
    print('Dataset 2 =', db_str2)
        
    #-----------------------------
    # Compare the two name lists
    #-----------------------------
#     print()   
#     diff1 = list( set(name_list1) - set(name_list2))
#     print('Basin names in dataset 1 not in dataset 2:')
#     print(diff1)
#     print()

    diff2 = list( set(name_list2) - set(name_list1))
    print('Basin names in dataset 2 not in dataset 1:')
    print(diff2)
    print('len(diff2) =', len(diff2))
    print()
    ## return diff2

#   compare_basin_names()
#---------------------------------------------------------------------
def get_dataset_key_names():
  
    names = [
    'CAMELS', 'MOPEX', 'NOAA_HADS',
    'NSF_CZO', 'NSF_LTER', 'NSF_NEON',
    'USDA_ARS', 'USGS_FPS', 'USGS_GAGES2_all',
    'USGS_GAGES2_ref', 'USGS_GAGES2_SB3', 
    'USGS_HCDN',  ###########################################
    'USGS_HLR', 'USGS_NWIS_Web', 'USGS_NWIS_WQP1',
    'USGS_NWIS_WQP2', 'USGS_NWIS_WQP3']
    return names

#   get_dataset_key_names()
#---------------------------------------------------------------------
def get_state_code_map():

    code_map = {
    'Alabama':        'AL',
    'Alaska':         'AK',
    'Arizona':        'AZ', 
    'Arkansas':       'AR',
    'California':     'CA',
    'Colorado':       'CO',
    'Connecticut':    'CT',
    'Delaware':       'DE',
    'Florida':        'FL',
    'Georgia':        'GA',
    'Hawaii':         'HI',
    'Idaho':          'ID',
    'Illinois':       'IL',
    'Indiana':        'IN',
    'Iowa':           'IA',
    'Kansas':         'KS',
    'Kentucky':       'KY',
    'Louisiana':      'LA',
    'Maine':          'ME',
    'Maryland':       'MD',
    'Massachusetts':  'MA',
    'Michigan':       'MI',
    'Minnesota':      'MN',
    'Mississippi':    'MS',
    'Missouri':       'MO',
    'Montana':        'MT',
    'Nebraska':       'NE',
    'Nevada':         'NV',
    'New Hampshire':  'NH',
    'New Jersey':     'NJ',
    'New Mexico':     'NM',
    'New York':       'NY',
    'North Carolina': 'NC',
    'North Dakota':   'ND',
    'Ohio':           'OH',
    'Oklahoma':       'OK',
    'Oregon':         'OR',
    'Pennsylvania':   'PA',
    'Puerto Rico':    'PR',
    'Rhode Island':   'RI',
    'South Carolina': 'SC',
    'South Dakota':   'SD',
    'Tennessee':      'TN',
    'Texas':          'TX',
    'Utah':           'UT',
    'Vermont':        'VT',
    'Virginia':       'VA',
    'Washington':     'WA',
    'West Virginia':  'WV',
    'Wisconsin':      'WI',
    'Wyoming':        'WY',
    #-----------------------------------
    # Canadian provinces & territories
    # that border the USA
    #-----------------------------------
    'Alberta':          'AB',
    'British Columbia': 'BC',
    'Manitoba':         'MB',
    'New Brunswick':    'NB',
    'Ontario':          'ON',
    'Quebec':           'QC',
    'Saskatchewan':     'SK',
    'Yukon Territory':  'YT',    
    #-----------------------
    # US territories, etc.
    #-----------------------
    'American Samoa':    'AS',
    'CNMI' :             'MP',
    'FSM'  :             'FM',
    'Guam' :             'GU',
    'Japan':             'JA',
    'Okinawa':           'JA',
    'Palau':             'PS',
    'Samoa':             'AS',
    'Saint Croix':       'VI',
    'Saint John':        'VI',
    'Saint Thomas':      'VI', 
    'USVI' :             'VI',
    'US Virgin Islands': 'VI',
    'Yellowsone National Park': 'YNP' }   # (spans: WY, MT, ID) 

    return code_map
        
#   get_state_code_map()
#---------------------------------------------------------------------
def get_state_code( usgs_id, long_name ):

    #-------------------------------------------------------
    # Note: It is assumed that long_name arg was created
    #       with the usgs.expand_usgs_site_name() function,
    #       which calls "replace_state_letters" and
    #       "fix_spellings".
    #-------------------------------------------------------
    # Many names contain "state line" and 2 states.
    #-------------------------------------------------------    
    state_code = ''
    if (long_name[-3] == ' '):
        state_code = long_name[-2:].upper()

    #----------------------------------------------
    # There are 2-letter FIPS codes and ISO codes
    # https://www.census.gov/library/reference/
    #      code-lists/ansi/ansi-codes-for-states.html
    #----------------------------------------------
    # AS = American Samoa (FIPS)
    # DC = District of Columbia (FIPS)
    # FM = Federal States of Micronesia (FIPS, appears as FSM)
    # GU = Guam (FIPS)
    # JA = Japan, Okinawa (FIPS)
    # MH = Marshall Islands (FIPS)
    # MP = Commonwealth of the Northern Mariana Islands (FIPS)
    # PR = Puerto Rico (FIPS)
    # PS = Palau (FIPS)
    # UM = U.S. Minor Outlying Islands (FIPS)
    # VI = Virgin Islands (FIPS, USVI)
    # YNP = Yellowstone National Park (spans: WY, MT, ID)
    # YT  = Yukon Territory (Canada)
     
    #----------------------------
    # Borders between US states
    #----------------------------
    if (long_name.lower().endswith('al-ga')):
        state_code = 'AL/GA'   
    #--------------------------------------------
    if (long_name.lower().endswith('az-ca')):
        state_code = 'CA/AZ'
    if (long_name.lower().endswith('az-cal')):
        state_code = 'CA/AZ'    
    #--------------------------------------------       
    if (long_name.lower().endswith('ca-az')):
        state_code = 'CA/AZ'
    #-------------------------------------------- 
    if (long_name.lower().endswith('colorado-nebraska')):
        state_code = 'CO/NE'
    #--------------------------------------------
    if (long_name.lower().endswith('il-ky')):
        state_code = 'IL/KY'
    #--------------------------------------------
    if (long_name.lower().endswith('ky-tn')):
        state_code = 'KY/TN' 
                    
    #------------------------------        
    # For places in Canada
    # Alberta, Ontario, Manitoba,
    # Quebec, Saskatchewan, etc.
    #------------------------------
    if (long_name.lower().endswith('alberta')):
        state_code = 'IB-ALB-CAN'
    if (long_name.lower().endswith('view alta')):
        state_code = 'IB-ALB-CAN'
    if (long_name.lower().endswith('manitoba')):
        state_code = 'IB-MAN-CAN'
    if (long_name.lower().endswith('manitoba canada')):
        state_code = 'IB-MAN-CAN'
    if (long_name.lower().endswith('ontario')):
        state_code = 'IB-ONT-CAN'
    if (long_name.lower().endswith('ontario can')):
        state_code = 'IB-ONT-CAN'
    if (long_name.lower().endswith('ontario canada')):
        state_code = 'IB-ONT-CAN'
    if (long_name.lower().endswith('quebec')):
        state_code = 'IB-QBC-CAN'
    if (long_name.lower().endswith('sask')):
        state_code = 'IB-SASK-CAN'
    if (long_name.lower().endswith('intnl bndary sask')):
        state_code = 'IB-SASK-CAN'
    if (long_name.lower().endswith(' yt')):
        state_code = 'IB-YUK-CAN'
                                                                                
    #------------------------------      
    # For CNMI (2 in USGS_gauged)
    #------------------------------
    if (long_name.lower().endswith('cnmi')):
        state_code = 'MP'

    #------------------------------      
    # For FSM (27 in USGS_gauged)
    #------------------------------
    if (long_name.lower().endswith('fsm')):
        state_code = 'FM'
        
    #----------------------------------------       
    # For Guam (19 in USGS_gauged)
    # Same string occurs in other names as:
    # GUAMANI (PR), STILLAGUAMISH (WA),
    #----------------------------------------
    if (long_name.lower().endswith('guam')):
        state_code = 'GU'

    #----------------------------------
    # "International Boundary"
    # Could be US/MEXICO or US/CANADA
    #----------------------------------
    if (long_name.lower().endswith("international boundary")):
        state_code = 'IB'
    if (long_name.lower().endswith("int'l boundary")):
        state_code = 'IB'
    if (long_name.lower().endswith("int boundary")):
        state_code = 'IB'
                        
    #-----------------------------------------        
    # For Okinawa (27 in USGS_gauged)
    # "Okinawa, Japan" occurs just once.
    #-----------------------------------------
    if (long_name.lower().endswith('okinawa')):
        state_code = 'JA'  # FIPS code
    if (long_name.lower().endswith('japan')):
        state_code = 'JA'
            
    #------------------------------------------        
    # For Palau (7 in USGS_gauged)
    # Same string occurs in:
    # PAPALAUA STREAM NEAR PUKOO, MOLOKAI, HI
    #------------------------------------------
    if (long_name.lower().endswith('palau')):
        state_code = 'PS'

    #-----------------------------------------
    # For American Samoa (17 in USGS_gauged)
    # "Am. Samoa" occurs 15 times
    # "American Samoa" occurs 2 times.
    # Recall that "." was removed already.
    #-----------------------------------------
    if (long_name.lower().endswith('samoa')):
        state_code = 'AS'  # FIPS code
#     if (long_name.lower().endswith('am samoa')):
#         state_code = 'AS'  # FIPS code
#     if (long_name.lower().endswith('american samoa')):
#         state_code = 'AS'  # FIPS code
                              
    #------------------------        
    # For US Virgin Islands
    #------------------------
    if (long_name[-4:].lower() == 'usvi'):
        state_code = 'VI'
    if (long_name.lower().endswith('saint croix')):
        state_code = 'VI'
    if (long_name.lower().endswith('saint john')):
        state_code = 'VI'
    if (long_name.lower().endswith('saint thomas')):
        state_code = 'VI'

    #----------------------------------        
    # For Yellowstone National Park
    #  (spans 3 states: WY, MT, ID)
    # 32 matches in USGS_gauged
    # 31 end in "YNP"
    # The remaining one is:
    # MIDDLE CR AT E ENTRANCE YNP  WY
    #----------------------------------
    if (long_name[-3:].lower() == 'ynp'):
        state_code = 'YNP'

    #-------------------------------------------
    # Option to try to find missing state code
    #-------------------------------------------
#     if (state_code == ''):
#         try:
#             state_name, state_code = usgs.get_usgs_missing_state( usgs_id, state_code_map)
#         except:
#             print('WARNING: Could not fix state code for:')
#             print('         site_id =', usgs_id)
#             print()
                
    return state_code
   
#   get_state_code()
#---------------------------------------------------------------------
def export_csv_to_tsv( csv_file='new_gages2_all.csv',
                       tsv_file='new_gages2_all.tsv',
                       data_dir=None):

    if (data_dir is None):
        data_dir = get_gages2_data_dir( gtype='root' )
    csv_path   = data_dir + csv_file
    tsv_path   = data_dir + tsv_file
    
    csv_unit   = open(csv_path, 'r')
    tsv_unit   = open(tsv_path, 'w')

    print('Exporting CSV to TSV...')
    while (True):
        line = csv_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        line2 = line.replace(',', '\t')
        line2 = line2.replace(';', ',')
        tsv_unit.write( line2 )
        
    csv_unit.close()
    tsv_unit.close()
    print('Finished.')
    print()

#   export_csv_to_tsv()
#---------------------------------------------------------------------
def convert_dms_to_dec_deg( dms_str, delim=None,
                            REMOVE_SYMBOLS=False ):

    #-----------------------------------------------------
    # Note: Set delim=None in split for any white space.
    #       Set REMOVE_SYMBOLS=True to remove the degree
    #       symbol, and the minute (') and second (")
    #       symbols from dms_str.
    #-----------------------------------------------------
    # Note: Could use this in usgs_utils.py.
    #-----------------------------------------------------    
    if (REMOVE_SYMBOLS):
        dms_str = dms_str.replace('&#176;',' ')  # degree symbol
        dms_str = dms_str.replace("'", ' ')
        dms_str = dms_str.replace('"', '')
        dms_str = dms_str.replace(',', '')
    
    dms_list = dms_str.split( delim )
    D = np.float64( dms_list[0] )
    M = np.float64( dms_list[1] )
    S = np.float64( dms_list[2] )
    if (D < 0):
        dec_deg_val = D - (M/60.0) - (S/3600.0)
    else:
        dec_deg_val = D + (M/60.0) + (S/3600.0)    

    return dec_deg_val

#   convert_dms_to_dec_deg()
#---------------------------------------------------------------------


