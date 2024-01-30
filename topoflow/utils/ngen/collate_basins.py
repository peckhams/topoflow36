
# Next Steps:
#
# * GAGES-II Selected Basins (there are 1947 of them):
#    - Can be assigned an HLR code using basin centroid lat/lon.
#    - Can be assigned an SWB/SWB2 class.
#    - Have basin shapefiles and have geographic bounding box.
#    - Have lots of additional metadata.
#    - Probably include many CAMELS, FPS, MOPEX, and RFC basins.
#    - Some of the included RFC basins likely have hydrograph type.
#
#---------------------------------------------------------------------
# Copyright (c) 2023-2024, Scott D. Peckham
#
# Jan 2024. Wrote write_new_values() due to multiple use.
#           Wrote more general version of get_new_usgs_line().
#           Renamed all "ngen/utils" files to end in "_utils.py"
#           instead of "_tools.py"
#           Modified to use new data_utils.py.
# Oct 2023. Modified collate() to get HLR code for any lon/lat
#           using new routines in hlr_utils.py.
#           Added get_lon_col, get_lat_col.
#           Moved some functions to usgs_utils.py:  haversine,
#           distance_on_sphere, get_closest_site.
#           Added RFC basins to collate() function.
# Sep 2023. Finished add_hlr_basins (new tools in usgs_utils.py).
#           Updated get_closest_site() for nodata lon or lat.
#           Started work on adding NOAA RFC basins.
# Jul 2023. Added haversine, distance_on_sphere,
#           get_closest_site, add_ars_basins, add_czo_basins,
#           add_hlr_basins, add_lter_basins, add_neon_basins,
#           and their subfunctions.
#           Moved all USGS-specific utils to usgs_utils.py.
# Jun 2023. Utils to collate several river basin datasets.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import collate_basins as cb
#  >>> cb.collate()
#
#---------------------------------------------------------------------
#
#  get_next_line()
#  get_next_id()
#  get_bounding_box()
#  skip_header_lines()
#  get_headings_line()
#  get_headings()
#  write_new_header()
#  write_new_header_OLD()   #########
#  write_new_values()
#
#  get_inclusion_keys()
#  get_inclusion_dict()
#  get_inclusion_list()
#  get_all_site_ids()
#  get_att_key_list()
#  get_heading_list()
#  get_default_att_dict()
#  update_attributes()
#
#  collate()                          ##### In progress #####
#  get_all_nwis_web_field_names()
#  get_new_usgs_line()
#
#  add_ars_basins()
#  add_czo_basins()
#  add_hlr_basins()
#  add_lter_basins()
#  add_neon_basins()
#
#---------------------------------------------------------------------
from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import usgs_utils as usgs
from topoflow.utils.ngen import usda_utils as usda
from topoflow.utils.ngen import hlr_utils  as hlr
from topoflow.utils.ngen import rfc_utils  as rfc
from topoflow.utils.ngen import czo_utils as cz

from topoflow.utils import regrid  # (to read hlr code geotiff)

import numpy as np
import time

# import re  # for regular expressions, not used now
# from osgeo import ogr, osr
# import json, sys, time

#---------------------------------------------------------------------
def get_next_line( key, file_units, NO_NEWLINE=True ):

    line = file_units[ key ].readline()
    if (NO_NEWLINE):
        line = line[:-1]
    return line
    
#   get_next_line()
#---------------------------------------------------------------------
def get_next_id( line, key, delim ):

    # print('### key =', key)
    id_col = dtu.get_id_column( key )
    vals   = line.split( delim )
    id     = vals[ id_col ].strip()
    # print('### id_col =', id_col)

    #---------------------------------------
    # Prepend '0' if only 7 digits for FPS
    #---------------------------------------    
    if (key == 'USGS_FPS') and (len(id) == 7):
        id = '0' + id
    return id
  
#   get_next_id()
#---------------------------------------------------------------------
def get_bounding_box( line, key, delim ):

    try:
        minlon_col = dtu.get_minlon_column( key, SILENT=True )
        if (minlon_col is None):
           bbox = ['-999','-999','-999','-999']  # strings
        else:
            vals   = line.split( delim )
            minlon = vals[ minlon_col     ].strip()
            maxlon = vals[ minlon_col + 1 ].strip()
            minlat = vals[ minlon_col + 2 ].strip()
            maxlat = vals[ minlon_col + 3 ].strip()
            bbox = [minlon, maxlon, minlat, maxlat]
    except:
        print('ERROR in get_bounding_box:')
        print('key =', key)
        print('line =')
        print(line)
        print()
    return bbox

#   get_bounding_box()
#---------------------------------------------------------------------
def get_headings_line( key='USGS_NWIS_Web'):

    file_path = dtu.get_new_tsv_filepath( key )
    file_unit = open( file_path, 'r' )
    nh_lines  = dtu.get_header_lines( key )
    nh_skip   = nh_lines

    #-----------------------------------------------
    # This one is unusual. It has another "format"
    # line after the main column headings line.
    #-----------------------------------------------    
    if (key == 'USGS_NWIS_Web'):
        nh_skip = (nh_lines - 1)
    
    for k in range( nh_skip ):
        line = file_unit.readline()        
    file_unit.close()
    
    headings_line = line  # last one read   
    return headings_line
     
#   get_headings_line()
#---------------------------------------------------------------------
def get_headings( key='USGS_NWIS_Web', delim='\t'):

    headings_line = get_headings_line( key=key)
    headings_line = headings_line[:-1]  # remove newline char
    headings      = headings_line.split( delim )
    return headings

#   get_headings()
#---------------------------------------------------------------------
def write_new_header( out_tsv_unit, delim ):

    att_headings = get_heading_list()  ## just main atts

    #----------------------------    
    # See: get_inclusion_keys()
    #----------------------------------
    # Last 4 don't have USGS site IDs
    #----------------------------------    
    incl_headings = [
    'Is_USGS_NWIS_Web', 'Is_GAGES2_Any',
    'Is_GAGES2_Ref', 'Is_GAGES2_SB3',
    'Is_FPS', 'Is_HCDN',
    'Is_RFC', 'Is_CAMELS', 'Is_MOPEX',
    'Is_CZO', 'Is_LTER', 'Is_NEON', 'Is_ARS' ]
    
    new_headings = att_headings + incl_headings

    header = ''
    for h in new_headings:
        header += (h + delim)
    header = header[:-1]  # remove last delim
    out_tsv_unit.write( header + '\n')
            
#   write_new_header()
#---------------------------------------------------------------------
def write_new_header_OLD( out_tsv_unit, delim ):
       
    #------------------------------------------------
    # Use a subset of the USGS headings.
    # These must agree with "fields_to_keep" in
    # get_new_usgs_line() and those written to TSV.
    #------------------------------------------------

    #-----------------------------------------------------   
    # These are for case where base_key = 'USGS_NWIS_Web'
    #-----------------------------------------------------    
    base_headings = [
    'Agency', 'SiteID', 'SiteName', 'SiteType',
    'SiteLat', 'SiteLon', 
    'CoordAcyCode', 'CoordDatum',
    ## 'DistrictCode', 'StateCode', 'CountyCode'
    'CountryCode', 
    'SiteElev', 'ElevAcy', 'ElevDatum',
    'HUC', 'BasinCode',
    ## 'DataTypesCode'
    'DrainArea', 'ReliabilityCode' ]
 
    new_headings  = [
    'LongName', 'State', 'SiteURL',
    'Closest_Site_ID', 'Closest_Site_Dist',
    'StartDate', 'EndDate',
    'MinLon', 'MaxLon', 'MinLat', 'MaxLat',
    #-------------------------------------------------------------
    'Is_USGS_NWIS_Web', 'Is_GAGES2_Ref',
    'Is_GAGES2_NonRef', 'Is_GAGES2_SB3',
    ###### 'Is_HCDN',
    'Is_CAMELS', 'Is_FPS', 'Is_MOPEX',
    'Is_ARS', 'Is_CZO', 'Is_HLR', 'Is_LTER', 'Is_NEON', 'Is_RFC',
    #-------------------------------------------------------------    
    'HLR_Code', 'SWB_Class', 'Hgraph_Type']
    new_headings = base_headings + new_headings
           
    header = ''
    for h in new_headings:
        header += (h + delim)
    header = header[:-1]  # remove last delim
    out_tsv_unit.write( header + '\n')
    
#   write_new_header_OLD()
#---------------------------------------------------------------------
def write_new_values( out_tsv_unit, val_list, delim ):

    new_line = ''
    for val in val_list:
        new_line += (val + delim)
    new_line = new_line[:-1]    # remove last delim
    out_tsv_unit.write( new_line + '\n' )
        
#   write_new_values()
#---------------------------------------------------------------------
def get_inclusion_keys():

    ikeys = [
    'USGS_NWIS_Web', 'USGS_GAGES2_all', 'USGS_GAGES2_ref',
    'USGS_GAGES2_SB3', 'USGS_FPS', 'USGS_HCDN',
    'NOAA_HADS', 'CAMELS', 'MOPEX' ]
    ### 'NSF_CZO', 'NSF_LTER', 'NSF_NEON', 'USDA_ARS']
 
    return ikeys
    
#   get_inclusion_keys()
#---------------------------------------------------------------------
def get_inclusion_dict( key_names=None ):
 
     #------------------------------------------------------  
    # Note: USGS_FPS has many strange (non-USGS) site IDs.
    #       Most have fewer than 8 characters.
    #-------------------------------------------------------
    # Note: USGS_HLR does not have USGS site IDs but has
    #       a numerical "VALUE".
    #-------------------------------------------------------    
    if (key_names is None):
        key_names = get_inclusion_keys()

    #---------------------------------------------------------
    # Note: This function returns a dictionary such that an
    # "inclusion" list is returned given a USGS site ID key.
    # This is very fast.
    #---------------------------------------------------------
    inclusion_dict   = dict()
    default_sub_dict = dict()
    for key in key_names:
        default_sub_dict[ key ] = 'N'

    for key in key_names:
        site_ids = dtu.get_basin_ids( key )
        for sid in site_ids:
            #-------------------------------------------------
            # Use next 2 lines to exclude non-USGS-style IDs
            #-------------------------------------------------
            if (len(sid) < 7) or (sid[0].isalpha()):
                continue
            if (sid in inclusion_dict):  
                inclusion_dict[ sid ][ key ] = 'Y'
            else:
                sub_dict = default_sub_dict.copy()   ## Need .copy!
                sub_dict[ key ] = 'Y'
                inclusion_dict[ sid ] = sub_dict
            #---------------------
            # Alternate approach
            #---------------------
#             if (sid in inclusion_dict):
#                 sub_dict = inclusion_dict[ sid ]
#             else:
#                 sub_dict = default_sub_dict.copy()
#             sub_dict[ key ] = 'Y'
#             inclusion_dict[ sid ] = sub_dict

    return inclusion_dict

#   get_inclusion_dict()
#---------------------------------------------------------------------
def get_inclusion_list( site_id, inclusion_dict=None ):

    if (inclusion_dict is None):
        inclusion_dict = get_inclusion_dict()
        
    sub_dict = inclusion_dict[ site_id ]
    inclusion_keys = get_inclusion_keys()
    inclusion_list = list()
    for key in inclusion_keys:
        inclusion_list.append( sub_dict[key] )

    ## Won't be in a particular order.
    ## inclusion_list = list( sub_dict.values() )

    return inclusion_list

#   get_inclusion_list()
#---------------------------------------------------------------------
def get_all_site_ids():

    d = get_inclusion_dict()
    ids = list( d.keys() )
    ids.sort()
    return ids
    
#   get_all_site_ids()
#---------------------------------------------------------------------
def get_att_key_list():

    att_key_list = [
    'site_id', 'nws_loc_id', 'rfc', 'huc', 'site_name',
    'site_type',
    'stage_data',
    'state_code', 'country_code',
    'lon', 'lat',
    'elev', 'elev_units',
    'area', 'area_units',
    'horiz_datum', 'vert_datum',
    'minlon', 'maxlon', 'minlat', 'maxlat',      
    'long_name', 'closest_site_id', 'closest_site_dist',
    'site_url', 'huc_url',
    'status', 'start_date', 'end_date',
    'hlr_code', 'swb_class', 'hgraph_type' ]

    return att_key_list
    
#   get_att_key_list()
#---------------------------------------------------------------------
def get_heading_list():

    #-------------------
    # Do this for now.
    #-------------------
    heading_list = list()
    att_key_list = get_att_key_list()
    for key in att_key_list:
        heading_list.append( key.title() )

    return heading_list

#   get_heading_list()
#---------------------------------------------------------------------
def get_default_att_dict( site_id ):

    att_dict = dict()
    #---------------------------------------  
    att_dict[ 'site_id' ]      = site_id
    att_dict[ 'site_name' ]    = '-'
    att_dict[ 'site_type' ]    = '-'
    #---------------------------------------
    att_dict[ 'nws_loc_id' ]   = '-'
    att_dict[ 'rfc' ]          = '-'
    att_dict[ 'huc' ]          = '-'
    #---------------------------------------
    #### att_dict[ 'has_discharge'] = '-'
    att_dict[ 'stage_data' ]    = '-'
    #---------------------------------------
    att_dict[ 'state_code' ]   = '-'
    att_dict[ 'country_code' ] = '-'
    #---------------------------------------
    att_dict[ 'lon' ]          = '-9999'
    att_dict[ 'lat' ]          = '-9999'
    att_dict[ 'elev' ]         = '-9999'
    att_dict[ 'area' ]         = '-9999'
    att_dict[ 'elev_units' ]   = '-'
    att_dict[ 'area_units' ]   = '-'
    att_dict[ 'horiz_datum' ]  = '-'
    att_dict[ 'vert_datum' ]   = '-'
    #---------------------------------------
    att_dict[ 'minlon' ]       = '-9999'
    att_dict[ 'maxlon' ]       = '-9999'
    att_dict[ 'minlat' ]       = '-9999'
    att_dict[ 'maxlat' ]       = '-9999'       
    #---------------------------------------
    att_dict[ 'long_name' ]    = '-'
    att_dict[ 'closest_site_id' ]   = '-'
    att_dict[ 'closest_site_dist' ] = '-'
    att_dict[ 'site_url' ]     = usgs.get_usgs_site_url( site_id )
    att_dict[ 'huc_url' ]      = '-'
    ##       usgs.get_usgs_huc_url( huc_num )
    #---------------------------------------
    att_dict[ 'status' ]       = '-'
    att_dict[ 'start_date' ]   = '-'
    att_dict[ 'end_date' ]     = '-'
    #---------------------------------------
    att_dict[ 'eco_region' ]   = '-'
    att_dict[ 'hlr_code' ]     = '0'  ######## OR use '-' ????
    att_dict[ 'swb_class' ]    = '-'
    att_dict[ 'hgraph_type' ]  = '-'
    #---------------------------------------    
    return att_dict

#   get_default_att_dict()
#---------------------------------------------------------------------
def update_attributes( val_list, key, att_dict ):

    UPDATE_ATTS = False
    all_headings = get_headings( key=key )
        
    if (key == 'USGS_NWIS_Web'):
        headings = ['site_no', 'station_nm', 'site_tp_cd',
        'dec_long_va', 'dec_lat_va', 'dec_coord_datum_cd',
        'state_cd', 'country_cd',
        'alt_va', 'alt_datum_cd', 'huc_cd',
        'drain_area_va',  #### 'contrib_drain_area_va',
        'sv_begin_date', 'sv_end_date']  ### 'sv_count'
        #----------------------------------------------------
        att_keys = ['site_id', 'site_name', 'site_type',
        'lon', 'lat', 'horiz_datum',
        'state_code', 'country_code',
        'elev', 'vert_datum', 'huc',
        'area', 'start_date', 'end_date']
        #------------------------------------
        # Spot-checked several site ID URLs
        #------------------------------------
        att_dict[ 'elev_units' ] = 'ft'
        att_dict[ 'area_units' ] = 'mi2'
        UPDATE_ATTS = True

        #--------------------------------------------
        # Now get "stage_data" from "data_types_cd"
        # This code has 30 characters, which may
        # be: A, I, or N
        # A = Active data collection site
        # I = Inactive or discontd data coll. site
        # N = Never collected ???
        #--------------------------------------------
        # See:
        # https://help.waterdata.usgs.gov/
        #   codes-and-parameters/type-of-data-collected
        # There is no code for "discharge"
        #------------------------------------------------
        dcol = all_headings.index( 'data_types_cd' )
        data_types = val_list[ dcol ].strip()
        if (len(data_types) >= 2):
            d1_pos = 0  # Stage or water levels--continuous
            d2_pos = 1  # Stage or water levels--intermittent
            d1_code = 'c' + data_types[ d1_pos ]
            d2_code = 'i' + data_types[ d2_pos ]
            att_dict[ 'stage_data' ] = d1_code + d2_code
    #------------------------------------------------------
    if (key == 'USGS_GAGES2_all'):
        headings = ['STAID', 'STANAME',
        ## 'CLASS', 
        'AGGECOREGI',
        'DRAIN_SQKM', ## 'HUC02',
        'LAT_GAGE', 'LNG_GAGE',
        'MINLON', 'MAXLON', 'MINLAT', 'MAXLAT',
        ## 'AREA', 'PERIMETER'   # Both in meters
        'STATE'] 
        ### 'HCDN_2009', 'ACTIVE09']
        ### 'FLYRS1900', 'FLYRS_1950', 'FLYRS_1990']
        #-----------------------------------------
        att_keys = ['site_id', 'site_name',
        'eco_region', 'area', 'lat', 'lon',
        'minlon', 'maxlon', 'minlat', 'maxlat',
        'state_code'] ###, 'status']
        #----------------------------------------------------------
        # NOTE: elev, elev_units, area, area_units would already
        #       have been set if this is an NWIS_Web site, and
        #       those units are "ft" and "mi2".
        #----------------------------------------------------------
        #       Even though GAGES2 column heading is 'DRAIN_SQKM'
        #       going to the site ID URL for several shows that
        #       area units are "mi2" and elev units are "ft".
        #----------------------------------------------------------
        area_units = att_dict[ 'area_units' ]
        if (area_units == '-'):
            att_dict[ 'area_units' ] = 'km2'
        #----------------------------------------
        elev_units = att_dict[ 'elev_units' ]
        if (elev_units == '-'):
            att_dict[ 'elev_units' ] = 'm'
        #-----------------------------------------------
        vert_datum = att_dict[ 'vert_datum' ]
        if (vert_datum == '-'):
            att_dict[ 'vert_datum' ] = 'CHECK_GAGES2'
        #-----------------------------------------------
        vert_datum = att_dict[ 'horiz_datum' ]
        if (vert_datum == '-'):
            att_dict[ 'horiz_datum' ] = 'NAD83'  # checked
        #-----------------------------------------------
        site_type = att_dict[ 'site_type' ]
        if (site_type == '-'): 
            att_dict[ 'site_type' ] = 'ST'  ####### (for all GAGES2?)
        #-----------------------------------------------
        country_code = att_dict[ 'country_code' ]
        if (country_code == '-'):          
            att_dict[ 'country_code'] = 'US'
        #----------------------------------------------------
        UPDATE_ATTS = True
    #------------------------------------------------------
    if (key == 'USGS_GAGES2_REF'):
        #---------------------------------------
        # Fully contained in 'USGS_GAGES2_all'
        #---------------------------------------
        pass 
    #------------------------------------------------------
    if (key == 'USGS_GAGES2_SB3'):
        #------------------------------------------------
        # Fully contained in 'USGS_GAGES2_REF', but
        # many more attributes available if needed,
        # such as:
        # LAT_CENT, LONG_CENT, met vars, and soil vars.
        #------------------------------------------------
        # All CAMELS basins are in this dataset.
        #-----------------------------------------
        ## att_dict[ 'swb_class' ] = val_list[68] # original SWB  
        att_dict[ 'swb_class' ] = val_list[69] # extended SWB
    #------------------------------------------------------
    if (key == 'NOAA_HADS'):
        headings = ['USGS_ID', 'USGS_name', 'Station_type',
        'Longitude', 'Latitude', 'RFC_ID', 'NWS_loc_ID']
        ### 'NWS_HSA_ID', 'GOES_ID', 'HUC8', 'Hgraph_type']
        UPDATE_ATTS = True
        #-----------------------------------------
        att_keys = ['site_id', 'site_name', 'site_type',
        'lon', 'lat', 'rfc', 'nws_loc_id'] 
        ### 'nws_hsa_id', 'goes_id', 'huc', 'hgraph_type']
        UPDATE_ATTS = True
    #------------------------------------------------------
    if (key == 'NOAA_via_API'):
        headings = ['usgs_id', 'nws_loc_id', 'rfc_abbrev',
        'usgs_name', 'state_code', 'longitude', 'latitude',
        'elevation' ]
        att_keys = ['site_id', 'nws_loc_id', 'rfc',
        'site_name', 'state_code', 'lon', 'lat', 'elev' ]
        ## 'in_service', 'pedts_obs', 'pedts_pred']
        UPDATE_ATTS = True    
    #------------------------------------------------------
    if (key == 'USGS_FPS'):
        #------------------------------------
        # Check/include the associated date     ############
        #------------------------------------
        headings = ['Status']   # Active or Inactive
        att_keys = ['status']
        UPDATE_ATTS = True
    #------------------------------------------------------
    if (key == 'USGS_HCDN'):
        headings = ['usgs_stn_id', 'usgs_stn_name', 'huc',
        'area_sqmi', 'state_code', 'gage_elev_ft']
        ## 'gage_lat_dms', 'gage_lon_dms'
        att_keys = ['site_id', 'site_name', 'huc',
        'area', 'state_code', 'elev']
        #----------------------------------------
        area_units = att_dict[ 'area_units' ]
        if (area_units == '-'):
            att_dict[ 'area_units' ] = 'mi2'
        #----------------------------------------
        elev_units = att_dict[ 'elev_units' ]
        if (elev_units == '-'):
            att_dict[ 'elev_units' ] = 'ft'
        UPDATE_ATTS = True
    #------------------------------------------------------
    if (key == 'MOPEX'):
        #-------------------------------------------
        # Try to get bounding box, if not set. 
        # Only a few missing (like 7) before this.
        #-------------------------------------------
        headings = ['SiteCode', 'SiteName', 'Latitude', 'Longitude',
        'MINLON', 'MAXLON', 'MINLAT', 'MAXLAT', 'State', 'HUC8',
        'DrainageAr']
        att_keys = [ 'site_id', 'site_name', 'lat', 'lon',
        'minlon', 'maxlon', 'minlat', 'maxlat', 'state_code',
        'huc', 'area']
        #----------------------------------------
        area_units = att_dict[ 'area_units' ]
        if (area_units == '-'):
            att_dict[ 'area_units' ] = 'mi2'  ######## CHECK
        #----------------------------------------
        elev_units = att_dict[ 'elev_units' ]
        if (elev_units == '-'):
            att_dict[ 'elev_units' ] = 'ft'   ######## CHECK
        UPDATE_ATTS = True
    #------------------------------------------------------
    if (UPDATE_ATTS):
        k = 0
        unset_list = ['-', 'unknown', '-9999']
        for h in headings:
            col = all_headings.index( h )
            att_key = att_keys[k]
            att_val = att_dict[ att_key ]     
            if (att_val in unset_list):
                val_str = val_list[ col ].strip()
                if (att_key == 'horiz_datum'):
                    # Remove the space character
                    # But Numbers app puts it back? ##########
                    val_str = val_str.replace(' ', '')
                att_dict[ att_key ] = val_str
            k += 1
        
#   update_attributes()
#---------------------------------------------------------------------
def get_primary_keys():

    #-------------------------------------------------------------
    # Note: We can only use this if all of the site IDs that
    #       were submitted to the API are in get_all_site_ids().
    #       Now using results from doing this for all site IDs
    #       in the USGS_NWIS_Web dataset.
    #-------------------------------------------------------------
    key_names = get_inclusion_keys()
    key_names = key_names + ['NOAA_via_API']
    return key_names

#   get_primary_keys()
#---------------------------------------------------------------------
def collate( out_tsv_file=None, max_count=1000, DEBUG=False ):

    #-------------------------------------------------------------
    # Note that the "USGS_NWIS_Web" dataset includes many gauges
    # that are no longer active. The "USGS_GAGES2_all" dataset
    # is newer and contains both "ref" and "non-ref" basins.
    #-------------------------------------------------------------
    #    key = 'USGS_NWIS_WQP3'     # 145375 site IDs
    #    key = 'USGS_NWIS_Web'      #  27915 site IDs
    #    key = 'USGS_GAGES2_all'    #   9322 site IDs
    #-------------------------------------------------------------
    start_time = time.time()

    if (out_tsv_file is None):
        out_tsv_file  = 'collated_basins_'
        out_tsv_file += str(max_count) + '.tsv'

    repo_dir = dtu.get_repo_dir()
    out_dir  = repo_dir + '__Collated/'   ##################
    out_tsv_path = out_dir + out_tsv_file
    out_tsv_unit = open( out_tsv_path, 'w')      #########

    #---------------------------------------  
    # Get key names for "primary" datasets
    #---------------------------------------
    ## key_names = get_inclusion_keys()
    key_names = get_primary_keys()  # See notes for this function.
    state_fips_code_map = usgs.get_state_fips_code_map()  # num to (2-letter, name)   
    state_code_map = dtu.get_state_code_map()  # full name to 2-letter
    n_bad_state_codes   = 0
    n_fixed_state_codes = 0
    
 
    #---------------------------------------------------
    # Get dictionary to map USGS ID to hydrograph type
    #---------------------------------------------------
    hydrograph_type_dict = rfc.get_hydrograph_type_dict()

    #-------------------------------------------------------
    # Get info to map HLR code to lon-lat pair
    # This info is passed to: hlr.get_hlr_code_for_point()
    #-------------------------------------------------------
    hlr_grid_info = hlr.get_hlr_grid_info()
    hlr_dir = dtu.get_new_data_dir( 'USGS_HLR' )
    hlr_raster_dir = hlr_dir + 'Raster/'  
    hlr_code_tif_file = hlr_grid_info['code_file']
    hlr_code_tif_path = hlr_raster_dir + hlr_code_tif_file
    hlr_code_grid = regrid.read_geotiff( in_file=hlr_code_tif_path,
                           GEO=False, REPORT=False )

    #-------------------------------------------
    # Get list of USGS site IDs for CZO basins
    #-------------------------------------------
    czo_site_ids = cz.get_usgs_site_ids()
    
    #------------------------------------------------
    # Open TSV files for each of the basin datasets
    #------------------------------------------------
    delim = '\t'  # tab
    file_units = dict()
    current_id   = dict()
    current_line = dict()
    for key in key_names:
        file_path = dtu.get_new_tsv_filepath( key )
        file_unit = open(file_path, 'r')
        file_units[ key ] = file_unit
        #--------------------        
        # Skip header lines
        #--------------------
        dtu.skip_header_lines( file_unit, key=key)
        #--------------------------------------
        # Read first site ID for each dataset
        #--------------------------------------
        line    = get_next_line( key, file_units )
        next_id = get_next_id( line, key, delim )
        current_id[ key ]   = next_id
        current_line[ key ] = line

    #-----------------------------------------
    # Write column headings for new TSV file
    #-----------------------------------------
    write_new_header( out_tsv_unit, delim)
    
    inclusion_dict  = get_inclusion_dict()
    att_key_list    = get_att_key_list()
    site_type_codes = ['ST', 'ST-CA', 'ST-DCH', 'ST-TS']
    site_type_map   = {'ST':'Stream', 'ST-CA':'Stream-Canal',
        'ST-DCH':'Stream-Ditch', 'ST-TS':'Stream-Tidal'}      
    #----------------------------------------------
    # Loop through all site_IDs, and use the fact
    # that each TSV is sorted by USGS site ID.
    #----------------------------------------------
    site_id_list = get_all_site_ids()
    counter = 0
    #### UPDATE next line
    print('Working on USGS, GAGES2, CAMELS, FPS, MOPEX, & RFC basins...')
        
    for site_id in site_id_list:
        #--------------------------------------------
        # Get the "inclusion" list that tells which
        # datasets include the current site_id
        #--------------------------------------------
        # CZO, LTER, NEON, and USDA-ARS don't have
        # USGS site IDs and are added at the end.
        #--------------------------------------------
        inclusion_list = get_inclusion_list( site_id, inclusion_dict=inclusion_dict )
        if (site_id in czo_site_ids):
            extra_list = ['Y', 'N', 'N', 'N']  # CZO, LTER, NEON, ARS
        else:
            extra_list = ['N', 'N', 'N', 'N']  # CZO, LTER, NEON, ARS
        inclusion_list = inclusion_list + extra_list 
      
        #------------------------------
        # Set defaults for attributes
        #------------------------------
        att_dict = get_default_att_dict( site_id )
      
        #---------------
        # For testing
        #---------------
        if (DEBUG):
#             for key in key_names:
#                 print(key + ' ID =', current_id[ key ])
            print('NWIS ID       =', current_id['USGS_NWIS_Web'])
            print('GAGES2_all ID =', current_id['USGS_GAGES2_all'])
            print('GAGES2_SB3 ID =', current_id['USGS_GAGES2_SB3'])
            print('CAMELS ID     =', current_id['CAMELS']) 
            print('FPS ID        =', current_id['USGS_FPS']) 
            print('MOPEX ID      =', current_id['MOPEX']) 
            print('NOAA HADS ID  =', current_id['NOAA_HADS']) 

        #--------------------------------------               
        # Does this dataset contain site_id ?
        #------------------------------------------------
        # Note that loop is over all possible site IDs.
        # Order of key_names determines att preference.
        #------------------------------------------------
        for key in key_names:
            db_site_id = current_id[ key ]                
            if (db_site_id == site_id):
                if (DEBUG):
                    print('  Found match for ' + key + '.')
                # val_list = get_line_values( current_line, key, delim)
                line     = current_line[ key ]
                val_list = line.split( delim )
                #-----------------------------------------
                # Update the attributes for this site_id
                #-----------------------------------------
                update_attributes( val_list, key, att_dict ) 

                #------------------------------------------------------
                # If bounding box is still not set, try to get it
                # Note that ordering of key_names matters here.
                #------------------------------------------------------
#                 bbox = get_bounding_box( line, key, delim )
#                 if (bbox[0] != '-9999'):
#                     att_dict[ 'minlon' ] = bbox[0]
#                     att_dict[ 'maxlon' ] = bbox[1]
#                     att_dict[ 'minlat' ] = bbox[2]
#                     att_dict[ 'maxlat' ] = bbox[3]
                #------------------------------------------------------                               
                line = get_next_line( key, file_units )
                if (line != ''):
                    next_id = get_next_id( line, key, delim )
                    current_id[ key ]   = next_id
                    current_line[ key ] = line

        #------------------------------------
        # Expand abbreviations in site_name
        #-------------------------------------------
        # Fix issue where names end like "Site 14"
        #   or "(Sw Site 3)" or "(S2)" or "(Sf8)"
        #-------------------------------------------
        # Expand "Irr Dist" and "Gv Irr Dist"
        # Expand "Xing" to "Crossing" ??
        #-------------------------------------------
        site_name = att_dict[ 'site_name' ]
        if (site_name != '-'):
            long_name  = usgs.get_usgs_long_name( site_name )
            att_dict[ 'long_name' ] = long_name

        #---------------------------------
        # Try to get 2-letter state code
        #---------------------------------
        state_code = att_dict[ 'state_code' ]
        if (state_code.isnumeric()):
            code, name = state_fips_code_map[ int(state_code) ]
            att_dict[ 'state_code' ] = code  # (2 letters)
        if (state_code == '-'):
            long_name  = att_dict[ 'long_name' ]
            if (long_name != '-'):
                site_id = att_dict[ 'site_id' ]
                code    = dtu.get_state_code( site_id, long_name )
                att_dict[ 'state_code' ] = code  # (2 letters)

        #------------------------------------------
        # Try to get site_type string from abbrev
        #------------------------------------------
        site_type = att_dict[ 'site_type' ]
        if (site_type in site_type_codes):
            stype = site_type_map[ site_type ]
            att_dict[ 'site_type' ] = stype
  
        #----------------------------------------------------
        # Get the USGS HLR code for this site via lon & lat
        #----------------------------------------------------
        lon = att_dict[ 'lon' ]
        lat = att_dict[ 'lat' ]
        if (lon.isnumeric() and lat.isnumeric()):
            if (lon != '-9999') and (lat != '-9999'):
                hlr_code = hlr.get_hlr_code_for_point( lon, lat,
                    hlr_grid=hlr_code_grid, grid_info=hlr_grid_info )
                att_dict[ 'hlr_code' ] = hlr_code  # as a string

        #-----------------------------------------------
        # Try to get the hydrograph type for this site
        #-----------------------------------------------
        if (site_id in hydrograph_type_dict):
            hgraph_info = hydrograph_type_dict[ site_id ]
            att_dict[ 'hgraph_type' ] = hgraph_info[ 'htype' ]
            att_dict[ 'rfc_name' ]    = hgraph_info[ 'rfc_name' ]
            att_dict[ 'huc8' ]        = hgraph_info[ 'huc8' ]
#             att_dict[ 'lon' ]         = hgraph_info[ 'lon' ]
#             att_dict[ 'lat' ]         = hgraph_info[ 'lat' ]

        #----------------------------------------
        # Set closest USGS site ID and distance
        #----------------------------------------
        site_id = att_dict[ 'site_id' ]
        att_dict[ 'closest_site_id' ]   = site_id
        att_dict[ 'closest_site_dist' ] = '0.0'
        #-------------------------------------
        # This is used for non-USGS site IDs
        #-------------------------------------
#         closest_id, clon, clat, dmin = \
#             usgs.get_closest_site(usgs_site_coords, 
#                 lon_str, lat_str, REPORT=False)
#         closest_site_id   = closest_id
#         closest_site_dist = '{x:.2f}'.format(x=dmin)  # string
         
        #--------------------------------------------
        # Build ordered att_list from att_dict here
        #--------------------------------------------
        att_list = list()
        for att_key in att_key_list:
            att_list.append( att_dict[ att_key ])
        
        out_val_list = (att_list + inclusion_list)
        write_new_values( out_tsv_unit, out_val_list, delim)
        counter += 1
        if ((counter % 100) == 0):
            print('  processed', counter, 'basins')
        if (counter == max_count):
            break

    #------------------------------------
    # Close the input files used so far
    #------------------------------------
    for key in key_names:
        file_units[ key ].close()
            
    #---------------------------------------------
    # Get coords of all USGS gauged basins so we
    # can find the closest site and distance
    #---------------------------------------------
    ## usgs_site_coords = usgs.get_usgs_site_coords( NWIS_ALL=False )
    usgs_site_coords = usgs.get_usgs_site_coords( NWIS_ALL=True )
        
    #-----------------------------------    
    # Add rows for the USDA ARS basins
    #-----------------------------------
    key = 'USDA_ARS'
    file_path = dtu.get_new_tsv_filepath( key )
    file_unit = open(file_path, 'r')
    dtu.skip_header_lines( file_unit, key=key)
    add_ars_basins( file_unit, out_tsv_unit,
                    usgs_site_coords,
                    hlr_grid_info, hlr_code_grid,
                    att_key_list, delim )
    file_unit.close()
    
    #------------------------------    
    # Add rows for the CZO basins
    #------------------------------
    key = 'NSF_CZO'
    file_path = dtu.get_new_tsv_filepath( key )
    file_unit = open(file_path, 'r')
    dtu.skip_header_lines( file_unit, key=key)
    add_czo_basins( file_unit, out_tsv_unit,
                    usgs_site_coords,
                    hlr_grid_info, hlr_code_grid,
                    att_key_list, delim )
    file_unit.close()

    #------------------------------
    # Add rows for the HLR basins
    #------------------------------
#     key = 'USGS_HLR'
#     file_path = dtu.get_new_tsv_filepath( key )
#     file_unit = open(file_path, 'r')
#     dtu.skip_header_lines( file_unit, key=key)
#     add_hlr_basins( file_unit, out_tsv_unit,
#                     usgs_site_coords,
#                     hlr_grid_info, hlr_code_grid,
#                     att_key_list, delim )
#     file_unit.close()
    
    #-------------------------------    
    # Add rows for the LTER basins
    #-------------------------------
    key = 'NSF_LTER'
    file_path = dtu.get_new_tsv_filepath( key )
    file_unit = open(file_path, 'r')
    dtu.skip_header_lines( file_unit, key=key)
    add_lter_basins( file_unit, out_tsv_unit,
                    usgs_site_coords,
                    hlr_grid_info, hlr_code_grid,
                    att_key_list, delim )
    file_unit.close()
    
    #-------------------------------    
    # Add rows for the NEON basins
    #-------------------------------
    key = 'NSF_NEON'
    file_path = dtu.get_new_tsv_filepath( key )
    file_unit = open(file_path, 'r')
    dtu.skip_header_lines( file_unit, key=key)
    add_neon_basins( file_unit, out_tsv_unit,
                    usgs_site_coords,
                    hlr_grid_info, hlr_code_grid,
                    att_key_list, delim )
    file_unit.close()
            
    #------------------  
    # Close all files
    #------------------
    out_tsv_unit.close()
    run_time = (time.time() - start_time)
    print('run_time =', run_time, '[secs]')
    ## print('n_bad_state_codes   =', n_bad_state_codes)
    ## print('n_fixed_state_codes =', n_fixed_state_codes)
    print('Finished.')
    print()
                     
#   collate()
#---------------------------------------------------------------------
def get_all_nwis_web_field_names():

    #------------------------------------------------------------
    # For USGS_NWIS_Web, here is the complete list of headings,
    # even though a user may have chosen a subset of these when
    # downloading data from the USGS NWIS Web site.
    #
    #  agency_cd       -- Agency
    #  site_no         -- Site identification number
    #  station_nm      -- Site name
    #  site_tp_cd      -- Site type
    #  lat_va          -- DMS latitude
    #  long_va         -- DMS longitude
    #  dec_lat_va      -- Decimal latitude
    #  dec_long_va     -- Decimal longitude
    #  coord_meth_cd   -- Latitude-longitude method
    #  coord_acy_cd    -- Latitude-longitude accuracy
    #  coord_datum_cd  -- Latitude-longitude datum
    #  dec_coord_datum_cd -- Decimal Latitude-longitude datum
    #  district_cd     -- District code
    #  state_cd        -- State code
    #  county_cd       -- County code
    #  country_cd      -- Country code
    #  land_net_ds     -- Land net location description
    #  map_nm          -- Name of location map
    #  map_scale_fc    -- Scale of location map
    #  alt_va          -- Altitude of Gage/land surface
    #  alt_meth_cd     -- Method altitude determined
    #  alt_acy_va      -- Altitude accuracy
    #  alt_datum_cd    -- Altitude datum
    #  huc_cd          -- Hydrologic unit code
    #  basin_cd        -- Drainage basin code
    #  topo_cd         -- Topographic setting code
    #  data_types_cd   -- Flags for the type of data collected
    #  instruments_cd  -- Flags for instruments at site
    #  construction_dt -- Date of first construction
    #  inventory_dt    -- Date site established or inventoried
    #  drain_area_va   -- Drainage area
    #  contrib_drain_area_va -- Contributing drainage area
    #  tz_cd           -- Mean Greenwich time offset
    #  local_time_fg   -- Local standard time flag
    #  reliability_cd  -- Data reliability code
    #  gw_file_cd      -- Data-other GW files
    #  nat_aqfr_cd     -- National aquifer code
    #  aqfr_cd         -- Local aquifer code
    #  aqfr_type_cd    -- Local aquifer type code
    #  well_depth_va   -- Well depth
    #  hole_depth_va   -- Hole depth
    #  depth_src_cd    -- Source of depth data
    #  project_no      -- Project number
    #  rt_bol          -- Real-time data flag
    #  peak_begin_date -- Peak-streamflow data begin date
    #  peak_end_date   -- Peak-streamflow data end date
    #  peak_count_nu   -- Peak-streamflow data count
    #  qw_begin_date   -- Water-quality data begin date
    #  qw_end_date     -- Water-quality data end date
    #  qw_count_nu     -- Water-quality data count
    #  gw_begin_date   -- Field water-level measurements begin date
    #  gw_end_date     -- Field water-level measurements end date
    #  gw_count_nu     -- Field water-level measurements count
    #  sv_begin_date   -- Site-visit data begin date
    #  sv_end_date     -- Site-visit data end date
    #  sv_count_nu     -- Site-visit data count
    #-----------------------------------------------------------------
    field_names = [
    'agency_cd', 'site_no', 'station_nm', 'site_tp_cd',
    'lat_va', 'long_va', 'dec_lat_va', 'dec_long_va',
    'coord_meth_cd', 'coord_acy_cd', 'coord_datum_cd',
    'dec_coord_datum_cd',
    'district_cd', 'state_cd', 'county_cd', 'country_cd',
    'land_net_ds', 'map_nm', 'map_scale_fc',
    'alt_va', 'alt_meth_cd', 'alt_acy_va', 'alt_datum_cd',
    'huc_cd', 'basin_cd', 'topo_cd',
    'data_types_cd', 'instruments_cd',
    'construction_dt', 'inventory_dt',
    'drain_area_va', 'contrib_drain_area_va',
    'tz_cd', 'local_time_fg',
    'reliability_cd',
    'gw_file_cd', 'nat_aqfr_cd', 'aqfr_cd', 'aqfr_type_cd',
    'well_depth_va', 'hole_depth_va', 'depth_src_cd',
    'project_no', 'rt_bol',
    'peak_begin_date', 'peak_end_date', 'peak_count_nu',
    'qw_begin_date', 'qw_end_date', 'qw_count_nu',
    'gw_begin_date', 'gw_end_date', 'gw_count_nu',
    'sv_begin_date', 'sv_end_date', 'sv_count_nu' ]   
    
    return field_names

#   get_all_nwis_web_field_names()
#---------------------------------------------------------------------
def get_new_usgs_line( usgs_line, base_key_headings, delim,
                       base_key='USGS_NWIS_Web' ):

    #------------------------------------------------------------        
    # Note: The "fields_to_keep" must match those in:
    #       write_new_header().
    #------------------------------------------------------------
    # For USGS_NWIS_Web, here is the complete list of headings,
    # even though a user may have chosen a subset of these when
    # downloading data from the USGS NWIS Web site.
    #
    #  agency_cd       -- Agency
    #  site_no         -- Site identification number
    #  station_nm      -- Site name
    #  site_tp_cd      -- Site type
    #  lat_va          -- DMS latitude
    #  long_va         -- DMS longitude
    #  dec_lat_va      -- Decimal latitude
    #  dec_long_va     -- Decimal longitude
    #  coord_meth_cd   -- Latitude-longitude method
    #  coord_acy_cd    -- Latitude-longitude accuracy
    #  coord_datum_cd  -- Latitude-longitude datum
    #  dec_coord_datum_cd -- Decimal Latitude-longitude datum
    #  district_cd     -- District code
    #  state_cd        -- State code
    #  county_cd       -- County code
    #  country_cd      -- Country code
    #  land_net_ds     -- Land net location description
    #  map_nm          -- Name of location map
    #  map_scale_fc    -- Scale of location map
    #  alt_va          -- Altitude of Gage/land surface
    #  alt_meth_cd     -- Method altitude determined
    #  alt_acy_va      -- Altitude accuracy
    #  alt_datum_cd    -- Altitude datum
    #  huc_cd          -- Hydrologic unit code
    #  basin_cd        -- Drainage basin code
    #  topo_cd         -- Topographic setting code
    #  data_types_cd   -- Flags for the type of data collected
    #  instruments_cd  -- Flags for instruments at site
    #  construction_dt -- Date of first construction
    #  inventory_dt    -- Date site established or inventoried
    #  drain_area_va   -- Drainage area
    #  contrib_drain_area_va -- Contributing drainage area
    #  tz_cd           -- Mean Greenwich time offset
    #  local_time_fg   -- Local standard time flag
    #  reliability_cd  -- Data reliability code
    #  gw_file_cd      -- Data-other GW files
    #  nat_aqfr_cd     -- National aquifer code
    #  aqfr_cd         -- Local aquifer code
    #  aqfr_type_cd    -- Local aquifer type code
    #  well_depth_va   -- Well depth
    #  hole_depth_va   -- Hole depth
    #  depth_src_cd    -- Source of depth data
    #  project_no      -- Project number
    #  rt_bol          -- Real-time data flag
    #  peak_begin_date -- Peak-streamflow data begin date
    #  peak_end_date   -- Peak-streamflow data end date
    #  peak_count_nu   -- Peak-streamflow data count
    #  qw_begin_date   -- Water-quality data begin date
    #  qw_end_date     -- Water-quality data end date
    #  qw_count_nu     -- Water-quality data count
    #  gw_begin_date   -- Field water-level measurements begin date
    #  gw_end_date     -- Field water-level measurements end date
    #  gw_count_nu     -- Field water-level measurements count
    #  sv_begin_date   -- Site-visit data begin date
    #  sv_end_date     -- Site-visit data end date
    #  sv_count_nu     -- Site-visit data count
    #-----------------------------------------------------------------
    # For USGS_NWIS_WQP,  here is the complete list of headings:
    #
    # OrganizationIdentifier
    # OrganizationFormalName
    # MonitoringLocationIdentifier
    # MonitoringLocationName
    # MonitoringLocationTypeName
    # MonitoringLocationDescriptionText
    # HUCEightDigitCode
    # DrainageAreaMeasure/MeasureValue
    # DrainageAreaMeasure/MeasureUnitCode
    # ContributingDrainageAreaMeasure/MeasureValue
    # ContributingDrainageAreaMeasure/MeasureUnitCode
    # LatitudeMeasure
    # LongitudeMeasure
    # SourceMapScaleNumeric
    # HorizontalAccuracyMeasure/MeasureValue
    # HorizontalAccuracyMeasure/MeasureUnitCode
    # HorizontalCollectionMethodName
    # HorizontalCoordinateReferenceSystemDatumName
    # VerticalMeasure/MeasureValue
    # VerticalMeasure/MeasureUnitCode
    # VerticalAccuracyMeasure/MeasureValue
    # VerticalAccuracyMeasure/MeasureUnitCode
    # VerticalCollectionMethodName
    # VerticalCoordinateReferenceSystemDatumName
    # CountryCode
    # StateCode
    # CountyCode
    # AquiferName
    # LocalAqfrName
    # FormationTypeText
    # AquiferTypeName
    # ConstructionDateText
    # WellDepthMeasure/MeasureValue
    # WellDepthMeasure/MeasureUnitCode
    # WellHoleDepthMeasure/MeasureValue
    # WellHoleDepthMeasure/MeasureUnitCode
    # ProviderName
    #-----------------------------------------------------------------    
    # Note: 'NAD27' shows as "NAD 27.00" in Numbers
    #       but string value is just 'NAD27'.
    #-----------------------------------------------------------------
    usgs_line = usgs_line[:-1]  # remove newline
    usgs_fields = usgs_line.split( delim )
    ## n_fields = len(usgs_fields)

    nwis_web_keys = ['USGS_NWIS_Web', 'USGS_NWIS_Web_Old']
    wqp_keys = ['USGS_NWIS_WQP1', 'USGS_NWIS_WQP2',
                'USGS_NWIS_WQP3']

    if (base_key in nwis_web_keys):
        fields_to_keep = [
        'agency_cd', 'site_no', 'station_nm', 'site_tp_cd',
        'dec_lat_va', 'dec_long_va',
        'coord_acy_cd', 'coord_datum_cd',
        ## 'district_cd', 'state_cd', 'county_cd',
        'country_cd',
        'alt_va', 'alt_acy_va', 'alt_datum_cd',
        'huc_cd', 
        ## 'basin_cd',   # usually missing
        ## 'data_types_cd',
        'drain_area_va' ]
        ## 'reliability_cd' ]  # usually missing
        ## 'sv_begin_date', 'sv_end_date', 'sv_count_nu' ] 
    elif (base_key in wqp_keys):
        fields_to_keep = [
        'OrganizationIdentifier', 
        ## OrganizationFormalName,
        'MonitoringLocationIdentifier',
        'MonitoringLocationName',
        'MonitoringLocationTypeName',
        ## MonitoringLocationDescriptionText,
        'HUCEightDigitCode',
        'DrainageAreaMeasure/MeasureValue',
        'DrainageAreaMeasure/MeasureUnitCode',
        'LatitudeMeasure',
        'LongitudeMeasure',
        'HorizontalCoordinateReferenceSystemDatumName',
        'VerticalMeasure/MeasureValue',
        'VerticalMeasure/MeasureUnitCode',
        'VerticalCoordinateReferenceSystemDatumName' ]
        ## 'CountryCode', 'StateCode', 'CountyCode' ] 
    else:
        fields_to_keep = []
        print('ERROR in get_new_usgs_line:')
        print('  base_key not supported:', base_key)

    #---------------------------------------------------    
    # Create new USGS line with a subset of all fields
    #---------------------------------------------------
    new_usgs_line = ''
    for field in fields_to_keep:
        col = base_key_headings.index( field )
        val = usgs_fields[ col ]
        new_usgs_line += (val + delim) 
    new_usgs_line = new_usgs_line[:-1]  # remove last delim 
   
    return new_usgs_line
    
#   get_new_usgs_line()
#---------------------------------------------------------------------
def add_ars_basins( ars_unit, out_tsv_unit,
                    usgs_site_coords,
                    hlr_grid_info, hlr_code_grid,
                    att_key_list, delim ):

    ars_id_col   = dtu.get_id_column( 'USDA_ARS' )
    ars_name_col = dtu.get_name_column( 'USDA_ARS' )

    print('Working on USDA/ARS basins...')
    while (True):
        ars_line = ars_unit.readline()
        if (ars_line == ''):
            break  # (reached end of file)
        vals = ars_line.split( delim )
        #-------------------------------------
        # Start with defaults for attributes
        #-------------------------------------
        att_dict = get_default_att_dict( site_id='-' )
        att_dict['agency']       = 'USDA-ARS'
        site_id = 'USDA-ARS-' + vals[ ars_id_col ].strip()
        att_dict['site_id']      = site_id
        att_dict['site_name']    = usda.get_ars_name( vals )
        att_dict['site_type']    = 'Stream' ### DOUBLE CHECK        
        att_dict['state_code']   = vals[2].strip()  # 2-letter state code
        att_dict['country_code'] = 'US'
        att_dict['long_name']    = usda.get_ars_long_name( vals )
        att_dict['lat']          = vals[7].strip()   # outlet lat
        att_dict['lon']          = vals[8].strip()   # outlet lon
        att_dict['area']         = vals[10].strip()  # basin area
        att_dict['area_units']   = 'km2'
        att_dict['start_date']   = vals[11].strip()  ###########
        att_dict['end_date']     = vals[12].strip()
        # This URL no longer works.
        state_code = vals[2].strip()
        att_dict['site_url']     = usda.get_ars_url( state_code )
        att_dict['notes']        = vals[15].strip()
        att_dict['horiz_datum']   = 'NAD83?' ### FIND OUT FOR SURE
        #---------------------------------------------------
        if (att_dict['lon'] == ''):
            att_dict['lon'] = '-9999'
        if (att_dict['lat'] == ''):
            att_dict['lat'] = '-9999'
        #----------------------------------------
        # Get closest USGS gauge from lon & lat
        #----------------------------------------
        lon = att_dict['lon']
        lat = att_dict['lat']           
        closest_id, clon, clat, dmin = \
            usgs.get_closest_site(usgs_site_coords, 
                lon, lat, REPORT=False)
        att_dict['closest_site_id']   = closest_id
        att_dict['closest_site_dist'] = '{x:.2f}'.format(x=dmin) 
        #---------------------------------------------
        # Don't know bounding box, elev, elev_units,
        # vert_datum, huc, etc. 
        #---------------------------------------------
        
        #----------------------------------------------------
        # Get the USGS HLR code for this site via lon & lat
        #----------------------------------------------------
        if (lon.isnumeric() and lat.isnumeric()):
            if (lon != '-9999') and (lat != '-9999'):
                hlr_code = hlr.get_hlr_code_for_point( lon, lat,
                    hlr_grid=hlr_code_grid, grid_info=hlr_grid_info )
                att_dict[ 'hlr_code' ] = hlr_code  # as a string         
 
        #--------------------------------------------
        # Build ordered att_list from att_dict here
        #--------------------------------------------
        att_list = list()
        for att_key in att_key_list:
            att_list.append( att_dict[ att_key ])

        # 'USGS_NWIS_Web', 'USGS_GAGES2_all', 'USGS_GAGES2_ref',
        # 'USGS_GAGES2_SB3', 'USGS_FPS', 'USGS_HCDN',
        # 'NOAA_HADS', 'CAMELS', 'MOPEX' ]
        ### 'NSF_CZO', 'NSF_LTER', 'NSF_NEON', 'USDA_ARS']
        #-------------------------
        # Get the inclusion list
        #-------------------------
        inclusion_list = list()
        for k in range(13):
            inclusion_list.append('N')
        inclusion_list[12] = 'Y'

        #------------------------------------------
        # Write out values for USDA ARS watershed
        #------------------------------------------        
        out_val_list = (att_list + inclusion_list)
        write_new_values( out_tsv_unit, out_val_list, delim)
        
#   add_ars_basins()
#---------------------------------------------------------------------
def add_czo_basins( czo_unit, out_tsv_unit,
                    usgs_site_coords,
                    hlr_grid_info, hlr_code_grid,
                    att_key_list, delim ):

    #------------------------------------------------
    # Note: CZO Lidar DEMs are available from:
    # https://portal.opentopography.org/dataSearch?
    #    search=Critical%20Zone%20Observatory
    #------------------------------------------------------------------
    # Johnston Draw in Reynolds Creek  ###########################
    # https://data.nal.usda.gov/search?query=%22Johnston%20Draw%22
    # https://data.nal.usda.gov/search?query=%22Reynolds%20Creek%22
    #---------------------------------------------------------------------------
    # https://data.nal.usda.gov/dataset/data-hydrological-modeling-dataset-
    # johnston-draw-catchment-reynolds-creek-experimental-watershed-idaho-usa
    #---------------------------------------------------------------------------
    # https://data.nal.usda.gov/dataset/data-eleven-years-mountain-weather-
    # snow-soil-moisture-and-stream-flow-data-rain-snow-transition-zone
    # -johnston-draw-catchment-reynolds-creek-experimental-watershed-
    # and-critical-zone-observatory-usa-v11
    #------------------------------------------------------------------
    print('Working on CZO basins...')
    while (True):
        # Don't add .strip() to czo_line;
        # will remove tabs for blank entries
        czo_line = czo_unit.readline()
        if (czo_line == ''):
            break  # (reached end of file)
        vals = czo_line.split( delim )
#         czo_oname  = vals[0].strip()
#         wshed_name = vals[1].strip()
#         location   = vals[3].strip()
        #-------------------------------------
        # Start with defaults for attributes
        #-------------------------------------  
        att_dict = get_default_att_dict( site_id='-' )
        site_id  = cz.get_czo_site_id( vals )
        att_dict['agency']       = 'NSF-CZO'
        att_dict['site_id']      = site_id
        # This function doesn't exist
        ## att_dict['site_name']    = cz.get_czo_site_name( vals )
        att_dict['site_name']    = cz.get_czo_long_name( vals )
        att_dict['site_type']    = 'Stream' ### DOUBLE CHECK        
        att_dict['state_code']   = vals[2].strip()  # 2-letter state code
        att_dict['country_code'] = 'US'
        att_dict['long_name']    = cz.get_czo_long_name( vals )
        att_dict['lat']          = vals[4].strip()   # outlet lat
        att_dict['lon']          = vals[5].strip()   # outlet lon
        att_dict['elev']         = vals[6].strip()   # elevation
        att_dict['elev_units']   = 'm'  # (masl)
        att_dict['area']         = vals[10].strip()  # basin area
        att_dict['area_units']   = 'km2'
        att_dict['start_date']   = vals[7].strip()  ###########
        att_dict['end_date']     = vals[8].strip()
        att_dict['site_url']     = cz.get_czo_site_url( site_id )
        att_dict['notes']        = vals[12].strip()   # credits
        att_dict['horiz_datum']  = 'NAD83?' ### FIND OUT FOR SURE
        #---------------------------------------------------
        if (att_dict['lon'] == ''):
            att_dict['lon'] = '-9999'
        if (att_dict['lat'] == ''):
            att_dict['lat'] = '-9999'
        #----------------------------------------
        # Get closest USGS gauge from lon & lat
        #------------------------------------------------
        # Most CZO sites are close to USGS NWIS gauges,
        # so distance should not be very big.
        #------------------------------------------------
        lon = att_dict['lon']
        lat = att_dict['lat']           
        closest_id, clon, clat, dmin = \
            usgs.get_closest_site(usgs_site_coords, 
                lon, lat, REPORT=False)
        att_dict['closest_site_id']   = closest_id
        att_dict['closest_site_dist'] = '{x:.2f}'.format(x=dmin) 
        #-------------------------------------------------
        # Don't know bounding box, vert_datum, huc, etc. 
        #-------------------------------------------------
        
        #----------------------------------------------------
        # Get the USGS HLR code for this site via lon & lat
        #----------------------------------------------------
        if (lon.isnumeric() and lat.isnumeric()):
            if (lon != '-9999') and (lat != '-9999'):
                hlr_code = hlr.get_hlr_code_for_point( lon, lat,
                    hlr_grid=hlr_code_grid, grid_info=hlr_grid_info )
                att_dict[ 'hlr_code' ] = hlr_code  # as a string         
 
        #--------------------------------------------
        # Build ordered att_list from att_dict here
        #--------------------------------------------
        att_list = list()
        for att_key in att_key_list:
            att_list.append( att_dict[ att_key ])

        # 'USGS_NWIS_Web', 'USGS_GAGES2_all', 'USGS_GAGES2_ref',
        # 'USGS_GAGES2_SB3', 'USGS_FPS', 'USGS_HCDN',
        # 'NOAA_HADS', 'CAMELS', 'MOPEX' ]
        ### 'NSF_CZO', 'NSF_LTER', 'NSF_NEON', 'USDA_ARS']
        #-------------------------
        # Get the inclusion list
        #-------------------------
        inclusion_list = list()
        for k in range(13):
            inclusion_list.append('N')
        inclusion_list[9] = 'Y'

        #------------------------------
        # Extra info for a few basins
        #-----------------------------------------------
        # Note that info from CUAHSI about the nearest
        # USGS site ID for 22 CZO basins was already
        # included before calling this function.
        #-----------------------------------------------
        site_id = att_dict['site_id']
        if (site_id == 'CZO-RC-JD'):
            att_dict['site_name'] = 'Johnston Draw in Reynolds Creek, ID'
            inclusion_list[12] = 'Y'  # Also a USDA-ARS basin
        elif (site_id == 'CZO-LUQ-IC'):
            # att_dict['site_id'] = '50075000'  ######
            att_dict['site_name'] = 'Rio Icacos near Naguabo, PR'
            inclusion_list[0] = 'Y'  # also a USGS NWIS ID
        if (site_id == 'CZO-LUQ-MAM'):
            ## att_dict['site_id'] = '50065500'
            att_dict['site_name'] = 'Rio Mameyes near Sabana, PR'
            inclusion_list[0] = 'Y'  # also a USGS NWIS ID
                      
        #------------------------------------------
        # Write out values for USDA ARS watershed
        #------------------------------------------        
        out_val_list = (att_list + inclusion_list)
        write_new_values( out_tsv_unit, out_val_list, delim)
    
#   add_czo_basins()
#---------------------------------------------------------------------
def add_hlr_basins( hlr_unit, out_tsv_unit,
                    usgs_site_coords,
                    hlr_grid_info, hlr_code_grid,
                    att_key_list, delim ):
                    
    hlr_id_col   = get_id_column( 'USGS_HLR' )
    ## hlr_name_col = get_name_column( 'USGS_HLR' )
    site_info = usgs.get_usgs_site_info_dict()
    mapped_to_gauge_count = 0

    print('Working on USGS HLR basins...')
    while (True):
        hlr_line = hlr_unit.readline()
        if (hlr_line == ''):
            break  # (reached end of file)
        vals = hlr_line.split( delim )
        #-------------------------------------
        # Start with defaults for attributes
        #-------------------------------------
        att_dict = get_default_att_dict( site_id='-' )
        att_dict['agency']       = 'USGS-HLR'
        att_dict['site_id']      = 'HLR-' + vals[ hlr_id_col ].strip()
#         att_dict['site_name']    = 
        att_dict['site_type']    = 'Stream' ### DOUBLE CHECK        
#         att_dict['state_code']   = # 2-letter state code
        att_dict['country_code'] = 'US'
#         att_dict['long_name']    = 
        att_dict['lat']          = vals[2].strip()   # outlet lat
        att_dict['lon']          = vals[1].strip()   # outlet lon
        att_dict['minlon']       = vals[5].strip()
        att_dict['maxlon']       = vals[6].strip()
        att_dict['minlat']       = vals[7].strip()
        att_dict['maxlat']       = vals[8].strip()
        # MINELE = outlet elevation in m
        att_dict['elev']         = vals[17].strip()   # elevation
        att_dict['elev_units']   = 'm'  # (masl)
        # COUNT = basin area in km2, since grid cell is 1 km2
        att_dict['area']         = vals[9].strip()  # basin area
        att_dict['area_units']   = 'km2'
#         att_dict['start_date']   = 
#         att_dict['end_date']     = 
#         att_dict['site_url']     = 
#         att_dict['notes']        = 
        att_dict['horiz_datum']  = 'NAD83?' ### FIND OUT FOR SURE
        #---------------------------------------------------
        if (att_dict['lon'] == ''):
            att_dict['lon'] = '-9999'
        if (att_dict['lat'] == ''):
            att_dict['lat'] = '-9999'

        #----------------------------------------
        # Get closest USGS gauge from lon & lat
        #------------------------------------------------
        # Most CZO sites are close to USGS NWIS gauges,
        # so distance should not be very big.
        #------------------------------------------------
        lon = att_dict['lon']
        lat = att_dict['lat']           
        closest_id, clon, clat, dmin = \
            usgs.get_closest_site(usgs_site_coords, 
                lon, lat, REPORT=False)
        att_dict['closest_site_id']   = closest_id
        att_dict['closest_site_dist'] = '{x:.2f}'.format(x=dmin) 
        #-------------------------------------------------
        # Don't know bounding box, vert_datum, huc, etc. 
        #-------------------------------------------------

        #-----------------------------------------------     
        # Try to get these from closest USGS site info
        #-----------------------------------------------
        IS_USGS_NWIS_ID = (dmin < 1.0)  # within 1 km (found 1506)
        ## IS_USGS_NWIS_ID = (dmin < 2.0)  # within 2 km (found 3046)
        if (IS_USGS_NWIS_ID):
            usgs_id = closest_id
            # att_dict['site_id'] = usgs_id
            # Use usgs_id to get more attributes.
            try:
                site_rec  = site_info[ usgs_id ]  # site_info is dict.
                site_name = site_rec['name']
                att_dict['site_name']  = site_name
                att_dict['long_name']  = usgs.get_usgs_long_name( site_name )  
                att_dict['state_code'] = site_rec['state']   # 2-letter state code
                att_dict['huc_code']   = site_rec['huc']       
                att_dict['site_url']   = site_rec['url']
                # Almost all of them & only one given
                # Note: NAD83 uses meters;  NAD93 uses feet?
                att_dict['horiz_daum'] = 'NAD83'
                mapped_to_gauge_count += 1
            except:
                pass
               
        #----------------------------------------------------
        # Get the USGS HLR code for this site via lon & lat
        #----------------------------------------------------
        att_dict['hlr_code'] = vals[22].strip() # code in {1,...,20
#         if (lon.isnumeric() and lat.isnumeric()):
#             if (lon != '-9999') and (lat != '-9999'):
#                 hlr_code = hlr.get_hlr_code_for_point( lon, lat,
#                     hlr_grid=hlr_code_grid, grid_info=hlr_grid_info )
#                 att_dict[ 'hlr_code' ] = hlr_code  # as a string         
 
        #--------------------------------------------
        # Build ordered att_list from att_dict here
        #--------------------------------------------
        att_list = list()
        for att_key in att_key_list:
            att_list.append( att_dict[ att_key ])

        # 'USGS_NWIS_Web', 'USGS_GAGES2_all', 'USGS_GAGES2_ref',
        # 'USGS_GAGES2_SB3', 'USGS_FPS', 'USGS_HCDN',
        # 'NOAA_HADS', 'CAMELS', 'MOPEX' ]
        ### 'NSF_CZO', 'NSF_LTER', 'NSF_NEON', 'USDA_ARS']
        #-------------------------
        # Get the inclusion list
        #-------------------------
        inclusion_list = list()
        for k in range(13):
            inclusion_list.append('N')
        ### inclusion_list[??] = 'Y'
        if (IS_USGS_NWIS_ID):
            inclusion_list[0] = 'Y'
                 
        #------------------------------------------
        # Write out values for USDA ARS watershed
        #------------------------------------------        
        out_val_list = (att_list + inclusion_list)
        write_new_values( out_tsv_unit, out_val_list, delim)
      
    print('# HLR basins mapped to USGS gauge =', mapped_to_gauge_count)

#   add_hlr_basins()
#---------------------------------------------------------------------
def add_lter_basins( lter_unit, out_tsv_unit,
                     usgs_site_coords,
                     hlr_grid_info, hlr_code_grid,
                     att_key_list, delim ):
                    
    #--------------------------------------------
    # Notes:  Metadata is from:
    # https://lternet.edu/site/
    # https://lternet.edu/site-characteristics/
    # https://lternet.edu/using-lter-data/
    #
    # LTER uses EDI data repository at:
    # https://portal.edirepository.org/nis/advancedSearch.jsp
    #
    # Connected to USDA somehow:
    # https://www.fs.usda.gov/rm/boise/AWAE/labs/
    #   awae_flagstaff/watersheds/fs-sites.html
    #
    # LTER "core" watersheds and areas
    # http://lter.konza.ksu.edu/lter-core-data
    #-------------------------------------------------------
    # See next URL for table of info on some small basins,
    #   and scroll down for data links.
    # https://www.fs.usda.gov/rm/boise/research/watershed/
    #   harvester/silvercreekdata.htm
    #
    # Benton Creek, Priest River Experimental Forest, ID
    # Bonanza Creek,
    # USGS: 15209750 BONANZA C NR MC CARTHY AK
    # Caspar Creek, Mendicino County, CA
    # Cedar Creek, 
    # Horse Creek, ID
    # Hubbard Brook,
    # Jornada Basin,
    # Konza Prairie LTER watersheds
    # Parker, Rowley, Ipswich (Plum Island Ecosystems LTER)
    # Silver Creek, ID
    #-------------------------------------------------------         
    print('Working on LTER basins...')
    while (True):
        lter_line = lter_unit.readline()
        if (lter_line == ''):
            break  # (reached end of file)
        vals = lter_line.split( delim )
        #-------------------------------------
        # Start with defaults for attributes
        #-------------------------------------
        att_dict = get_default_att_dict( site_id='-' )
        att_dict['agency']       = 'NSF-LTER'
        att_dict['site_id']      = 'LTER-' + vals[0].strip()
        att_dict['site_name']    = vals[1].strip()
#        att_dict['site_type']    = 'Stream' ### DOUBLE CHECK
        #-------------------------------------------
        # Could extract state codes from Site Name
        #-------------------------------------------        
#         att_dict['state_code']   = # 2-letter state code
        att_dict['country_code'] = 'US'
#         att_dict['long_name']    = 
        att_dict['lat']          = vals[4].strip()   # outlet lat
        att_dict['lon']          = vals[5].strip()   # outlet lon
        att_dict['elev']         = vals[6].strip()   # elevation
        att_dict['elev_units']   = '-'
#         att_dict['area']         =  # basin area
#         att_dict['area_units']   = '-'
#         att_dict['start_date']   = 
#         att_dict['end_date']     = 
#         att_dict['site_url']     = 
#         att_dict['notes']        = 
#         att_dict['horiz_datum']  = 'NAD83?' ### FIND OUT FOR SURE
        #---------------------------------------------------
        if (att_dict['lon'] == ''):
            att_dict['lon'] = '-9999'
        if (att_dict['lat'] == ''):
            att_dict['lat'] = '-9999'
        #----------------------------------------
        # Get closest USGS gauge from lon & lat
        #------------------------------------------------
        # Most LTER sites are close to USGS NWIS gauges,
        # so distance should not be very big.
        #------------------------------------------------
        lon = att_dict['lon']
        lat = att_dict['lat']          
        closest_id, clon, clat, dmin = \
            usgs.get_closest_site(usgs_site_coords, 
                lon, lat, REPORT=False)
        att_dict['closest_site_id']   = closest_id
        att_dict['closest_site_dist'] = '{x:.2f}'.format(x=dmin)
        #----------------------------------
        # Update the site URL if missing?
        #----------------------------------
        site_url = att_dict[ 'site_url' ]
        if (site_url[-1] == '-'):
            att_dict[ 'site_url' ] = usgs.get_usgs_site_url( closest_id )
 
        #-------------------------------------------------
        # Don't know bounding box, vert_datum, huc, etc. 
        #-------------------------------------------------
           
        #----------------------------------------------------
        # Get the USGS HLR code for this site via lon & lat
        #----------------------------------------------------
        if (lon.isnumeric() and lat.isnumeric()):
            if (lon != '-9999') and (lat != '-9999'):
                hlr_code = hlr.get_hlr_code_for_point( lon, lat,
                    hlr_grid=hlr_code_grid, grid_info=hlr_grid_info )
                att_dict[ 'hlr_code' ] = hlr_code  # as a string         
 
        #--------------------------------------------
        # Build ordered att_list from att_dict here
        #--------------------------------------------
        att_list = list()
        for att_key in att_key_list:
            att_list.append( att_dict[ att_key ])

        # 'USGS_NWIS_Web', 'USGS_GAGES2_all', 'USGS_GAGES2_ref',
        # 'USGS_GAGES2_SB3', 'USGS_FPS', 'USGS_HCDN',
        # 'NOAA_HADS', 'CAMELS', 'MOPEX' ]
        ### 'NSF_CZO', 'NSF_LTER', 'NSF_NEON', 'USDA_ARS']
        #-------------------------
        # Get the inclusion list
        #-------------------------
        inclusion_list = list()
        for k in range(13):
            inclusion_list.append('N')
        inclusion_list[10] = 'Y'
                 
        #------------------------------------------
        # Write out values for USDA ARS watershed
        #------------------------------------------        
        out_val_list = (att_list + inclusion_list)
        write_new_values( out_tsv_unit, out_val_list, delim)
      
#   add_lter_basins()
#---------------------------------------------------------------------
def add_neon_basins( neon_unit, out_tsv_unit,
                     usgs_site_coords,
                     hlr_grid_info, hlr_code_grid,
                     att_key_list, delim ):

    #--------------------------------------------------------------
    # Notes:  Metadata is from:
    # https://www.neonscience.org/field-sites/explore-field-sites
    #--------------------------------------------------------------       
    print('Working on NEON basins...')
    while (True):
        neon_line = neon_unit.readline()
        if (neon_line == ''):
            break  # (reached end of file)
        vals = neon_line.split( delim )
        #-------------------------------------------------------------
        field_domain_id        = vals[0].strip()
        field_site_id          = vals[1].strip()
        field_site_name        = vals[2].strip()
#         field_site_type        = vals[3].strip()
#         field_site_subtype     = vals[4].strip()
#         field_collocated_site  = vals[5].strip()
#         field_site_host        = vals[6].strip()
        field_site_url         = vals[7].strip()
        field_nonneon_res_alwd = vals[8].strip()
        field_access_details   = vals[9].strip()
        field_neon_op_office   = vals[10].strip()
        field_latitude         = vals[11].strip()
        field_longitude        = vals[12].strip()
        field_geodetic_datum   = vals[13].strip()
        field_utm_north        = vals[14].strip()
        field_utm_easting      = vals[15].strip()
        field_utm_zone         = vals[16].strip()
#         field_site_county      = vals[17].strip()
        field_site_state       = vals[18].strip()  # 2-letter code
        field_site_country     = vals[19].strip()
#         field_site_mean_elev_m = vals[20].strip()
        field_site_min_elev_m  = vals[21].strip()  ####
#         field_site_max_elev_m  = vals[22].strip()
#         field_mean_annu_temp_c = vals[23].strip()
        #-----------------------------------------------
#         field_mean_ann_precip_mm = vals[24].strip()
#         field_dom_wind_direction = vals[25].strip()
#         field_mean_canopy_height = vals[26].strip()
#         field_dom_nlcd_classes   = vals[27].strip()
#         field_dom_plant_species  = vals[28].strip()
        # e.g. [h10250001](https://water.usgs.gov/lookup/getwatershed?10250001)
        field_usgs_huc           = vals[29].strip()                       
        field_watershed_name     = vals[30].strip()
        field_watershed_size_km2 = vals[31].strip()
#         field_lake_depth_mean_m  = vals[32].strip()
#         field_lake_depth_max_m   = vals[33].strip()
#         field_tower_height_m     = vals[34].strip()
#         field_usgs_geology_unit  = vals[35].strip()
#         field_megapit_soil_famly = vals[36].strip()
#         field_soil_subgroup      = vals[37].strip()
#         field_avg_num_green_days = vals[38].strip()
#         field_avg_green_incr_doy = vals[39].strip()
#         field_avg_green_max_doy  = vals[40].strip()
#         field_avg_green_decr_doy = vals[41].strip()
#         field_avg_green_min_doy  = vals[42].strip()
#         field_phenocams          = vals[43].strip()
#         field_num_tower_levels   = vals[44].strip()       
        #-----------------------------------------------
        p1 = field_usgs_huc.find('[')
        p2 = field_usgs_huc.find(']')
        huc_str = field_usgs_huc[p1+2: p2]
        huc_url = 'https://water.usgs.gov/lookup/getwatershed?' + huc_str 
        #---------------------------------------------------------            
        site_id   = 'NEON-' + field_domain_id + '-' + field_site_id
        long_name = (field_site_name + ', ' + 
                     field_watershed_name + ', ' +
                     field_site_state)
        #---------------------------------------------------------
        country_code = field_site_country
        if (field_site_country == 'USA'):
            country_code = 'US'
                    
        #-------------------------------------
        # Start with defaults for attributes
        #-------------------------------------  
        att_dict = get_default_att_dict( site_id='-' )
        att_dict['agency']       = 'NSF-NEON'
        att_dict['site_id']      = site_id
        att_dict['site_name']    = field_watershed_name
        att_dict['site_type']    = 'Stream' ### DOUBLE CHECK        
        att_dict['state_code']   = field_site_state  # 2-letter state code
        att_dict['country_code'] = country_code
        att_dict['long_name']    = long_name
        att_dict['lat']          = field_latitude    # outlet lat
        att_dict['lon']          = field_longitude   # outlet lon
        att_dict['elev']         = field_site_min_elev_m  # elevation
        att_dict['elev_units']   = 'm'
        att_dict['area']         = field_watershed_size_km2  # basin area
        att_dict['area_units']   = 'km2'
#         att_dict['start_date']   = 
#         att_dict['end_date']     = 
        att_dict['site_url']     = field_site_url
#         att_dict['notes']        =
        att_dict['horiz_datum']  = field_geodetic_datum
        att_dict['huc']          = huc_str
        att_dict['huc_url']      = huc_url
        #---------------------------------------------------
        if (att_dict['lon'] == ''):
            att_dict['lon'] = '-9999'
        if (att_dict['lat'] == ''):
            att_dict['lat'] = '-9999'
        #----------------------------------------
        # Get closest USGS gauge from lon & lat
        #------------------------------------------------
        # Most NEON sites are close to USGS NWIS gauges,
        # so distance should not be very big.
        #------------------------------------------------
        lon = att_dict['lon']
        lat = att_dict['lat']           
        closest_id, clon, clat, dmin = \
            usgs.get_closest_site(usgs_site_coords, 
                lon, lat, REPORT=False)
        att_dict['closest_site_id']   = closest_id
        att_dict['closest_site_dist'] = '{x:.2f}'.format(x=dmin)
        #-----------------------------------
        # Update the site_url if missing ?
        #-----------------------------------
        if (att_dict[ 'site_url' ] == ''):
            att_dict[ 'site_url' ] = usgs.get_usgs_site_url( closest_id )
 
        #-------------------------------------------------
        # Don't know bounding box, vert_datum, huc, etc. 
        #-------------------------------------------------
        
        #----------------------------------------------------
        # Get the USGS HLR code for this site via lon & lat
        #----------------------------------------------------
        if (lon.isnumeric() and lat.isnumeric()):
            if (lon != '-9999') and (lat != '-9999'):
                hlr_code = hlr.get_hlr_code_for_point( lon, lat,
                    hlr_grid=hlr_code_grid, grid_info=hlr_grid_info )
                att_dict[ 'hlr_code' ] = hlr_code  # as a string         
 
        #--------------------------------------------
        # Build ordered att_list from att_dict here
        #--------------------------------------------
        att_list = list()
        for att_key in att_key_list:
            att_list.append( att_dict[ att_key ])

        # 'USGS_NWIS_Web', 'USGS_GAGES2_all', 'USGS_GAGES2_ref',
        # 'USGS_GAGES2_SB3', 'USGS_FPS', 'USGS_HCDN',
        # 'NOAA_HADS', 'CAMELS', 'MOPEX' ]
        ### 'NSF_CZO', 'NSF_LTER', 'NSF_NEON', 'USDA_ARS']
        #-------------------------
        # Get the inclusion list
        #-------------------------
        inclusion_list = list()
        for k in range(13):
            inclusion_list.append('N')
        inclusion_list[11] = 'Y'

        #--------------------------------------
        # Write out values for NEON watershed
        # Skip any other NEON sites
        #--------------------------------------
        if (field_watershed_name != ''):      
            out_val_list = (att_list + inclusion_list)
            write_new_values( out_tsv_unit, out_val_list, delim)
        
#   add_neon_basins()
#---------------------------------------------------------------------
     





