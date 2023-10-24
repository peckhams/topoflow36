
# Copyright (c) 2023, Scott D. Peckham
#
# Oct 2023. create_tsv() and all supporting functions working. 
# Sep 2023. Started.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36
#  % python
#  >>> from topoflow.utils.ngen import rfc_tools as rfc
#  >>> rfc.create_tsv()
#  >>> rfc.sort_by_site_code()
#
#---------------------------------------------------------------------
#
#  get_basin_repo_dir()
#  get_rfc_data_dir()
#  get_usgs_rfc_dict()
#  get_hydrograph_type_dict()
#  create_tsv()
#  sort_by_site_code()
#
#---------------------------------------------------------------------

import numpy as np
import os, os.path

import pickle    # to save usgs station info map
from topoflow.utils.ngen import usgs_utils as usgs
# import time

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
def get_rfc_data_dir():

    #-----------------------------------
    # Modify this directory as needed.
    #-----------------------------------
    repo_dir  = get_basin_repo_dir()
    data_dir  = repo_dir + 'NOAA_RFCs/'
    data_dir += 'USGS-HADS_basins/'
    return data_dir
       
#   get_rfc_data_dir()
#---------------------------------------------------------------------
def get_usgs_rfc_dict():

    rfc_dir   = get_rfc_data_dir()
    dict_file = 'USGS_RFC_station_info_dict.pkl'
    file_path = rfc_dir + dict_file
    if (os.path.exists( file_path )):
        print('Reading saved RFC station info dictionary...')
        file_unit = open( file_path, 'rb')
        station_info = pickle.load( file_unit )
        file_unit.close()
        print('Finished.')
        print()
        return station_info
    
    print('Getting USGS/RFC station info as dictionary...')
    info_file = 'usgsMasterMeta.txt'
    info_path = rfc_dir + info_file
    info_unit = open( info_path, 'r' )
    delim     = ','

    #----------------------------------
    # Skip over the header (one line)
    #----------------------------------
    line = info_unit.readline()

    #-------------------------------------------------------------    
    # Construct the station_info dictionary, where USGS gauge id
    # is the key used to get a dictionary with info for a gauge
    # The file: usgsMasterMeta.txt has info for 8170 basins.
    # Note that this file is sorted by USGS site ID.
    #-------------------------------------------------------------
    station_info = dict()
    k = 0

    while (True):
        info_line = info_unit.readline()
        if (info_line == ''):
            break  # (reached end of file)
        vals = info_line.split( delim )
        sid  = vals[0].strip()
        lat  = vals[1].strip()
        lon  = vals[2].strip()
        huc8 = vals[6].strip()
        rfc  = vals[7].strip()
        if (rfc.strip() == ''):
            rfc = 'unknown'
        ## url  = usgs.get_usgs_site_url( sid )
        station_info[ sid ] = \
            {'rfc_name':rfc, 'huc8':huc8, 'lat':lat, 'lon':lon }
        k += 1

    print('Processed', k, 'USGS/RFC stations.')
    
    SAVE_TO_FILE = True
    if (SAVE_TO_FILE):
        rfc_dir   = get_rfc_data_dir()
        file_path = rfc_dir + dict_file
        file_unit = open( file_path, 'wb' )
        pickle.dump( station_info, file_unit, protocol=pickle.HIGHEST_PROTOCOL)
        file_unit.close()
        print('Saved USGS/RFC station info dictionary to file:')
        print('  ' + file_path)
        print()

    return station_info
                    
#   get_usgs_rfc_dict()
#---------------------------------------------------------------------
def get_hydrograph_type_dict():

    rfc_dir   = get_rfc_data_dir()
    dict_file = 'RFC_hydrograph_type_dict.pkl'
    file_path = rfc_dir + dict_file
    if (os.path.exists( file_path )):
        print('Reading saved RFC hydrograph type dictionary...')
        file_unit = open( file_path, 'rb')
        hydro_type = pickle.load( file_unit )
        file_unit.close()
        print('Finished.')
        print()
        return hydro_type
    
    print('Getting RFC hydrograph type as dictionary...')

    #-------------------------------------------------    
    # Open each of the hydrograph type files to read
    # Each file is sorted by USGS station ID
    #-------------------------------------------------
    delim = ','
    htype_dict = dict()

    #--------------------------------    
    # Make a list of CSV file names
    #--------------------------------
    file_names = list()
    file_names.append( 'NWMv3_basin_type_flashy.txt' ) 
    file_names.append( 'NWMv3_basin_type_slow.txt' )
    file_names.append( 'NWMv3_basin_type_regulation.txt' )
    file_names.append( 'NWMv3_basin_type_snow.txt' )

    #----------------------------------    
    # Make a list of hydrograph types
    #----------------------------------
    htypes = ['flashy', 'slow', 'regulated', 'snow-dom']
    
    #------------------------------------------------   
    # Create a dictionary with hydrograph types
    # and other info for each USGS station ID found
    #------------------------------------------------
    for k in range(len(file_names)):
        htype = htypes[k]
        file_path = rfc_dir + file_names[k]
        file_unit = open( file_path, 'r' )
        line = file_unit.readline()  # skip 1 header line
        while (True):
            line = file_unit.readline()
            if (line == ''):
                break  # (reached end of file)
            vals = line.split( delim )
            sid  = vals[0].strip()
            lat  = vals[1].strip()
            lon  = vals[2].strip()
            huc8 = vals[6].strip()
            rfc  = vals[7].strip()
            ## url  = usgs.get_usgs_site_url( sid )
            if (lat != 'NA'):  ################### ADD MISSING INFO LATER #######
                htype_dict[ sid ] = \
                    {'htype':htype, 'rfc_name':rfc, 'huc8':huc8,
                     'lat':lat, 'lon':lon }       
        file_unit.close() 

    return htype_dict 
    
#   get_hydrograph_type_dict()
#---------------------------------------------------------------------
def create_tsv( data_dir=None, tsv_file='new_rfc_hads.tsv',
                REPORT=True):

    if (data_dir is None):
        data_dir = get_rfc_data_dir()
    hads_file = 'hads.ncep.noaa.gov_USGS_ALL_USGS-HADS_SITES.txt'
    hads_path  = data_dir + hads_file
    hads_unit  = open( hads_path, 'r' )
    hads_delim = '|'

    #---------------------------------------
    # Skip over HADS file header (4 lines)
    #---------------------------------------
    for i in range(4):
       line = hads_unit.readline()

    #-----------------------------
    # Open new TSV file to write
    #-----------------------------
    tsv_path = data_dir + tsv_file
    tsv_unit = open( tsv_path, 'w')

    #--------------------------------
    # Write header for new TSV file
    #--------------------------------
    delim  = '\t'  # tab character  
    header = ''
    header += 'USGS_ID'      + delim
    header += 'USGS_name'    + delim
    header += 'Station_type' + delim
    header += 'Longitude'    + delim
    header += 'Latitude'     + delim
    header += 'RFC_ID'       + delim
    header += 'NWS_loc_ID'   + delim
    header += 'NWS_HSA_ID'   + delim
    header += 'GOES_ID'      + delim
    header += 'HUC8'         + delim
    header += 'Hgraph_type'  + '\n'  # newline at end  
    tsv_unit.write( header )   

    rfc_map   = get_usgs_rfc_dict()
    htype_map = get_hydrograph_type_dict()
    null_info = {'rfc_name': 'unknown', 'huc8': 'unknown',
                 'lat':'-999', 'lon':'-999'}
    station_type_dict = usgs.get_station_type_dict()

    k = 0
    lat_diff_count = 0
    lon_diff_count = 0

    while (True):
        hads_line = hads_unit.readline()
        if (hads_line == ''):
            break  # (reached end of file)
        vals = hads_line.split( hads_delim )
        nws_loc_id = vals[0].strip()
        usgs_id    = vals[1].strip()
        goes_id    = vals[2].strip()
        nws_hsa_id = vals[3].strip()
        lat_dms    = vals[4].strip()
        lon_dms    = vals[5].strip()
        if (lon_dms[0] != '-'):
            #--------------------------
            # Add missing minus sign.
            # One entry does have one.
            #--------------------------
            lon_dms = '-' + lon_dms
        loc_name   = vals[6].strip()
        #-----------------------------------
        try:
            rfc_info = rfc_map[ usgs_id ]
        except:
            rfc_info = null_info
        rfc_name   = rfc_info['rfc_name']
        huc8_code  = rfc_info['huc8']
        rfc_lat    = rfc_info['lat']
        rfc_lon    = rfc_info['lon']
        if (rfc_lat == 'NA'):
            rfc_lat = '-999'
        if (rfc_lon == 'NA'):
            rfc_lon = '-999'
        #-----------------------------------
        dms_list = lat_dms.split()  # split on white space
        D = np.float64( dms_list[0] )
        M = np.float64( dms_list[1] )
        S = np.float64( dms_list[2] )
        if (D < 0):
            lat = D - (M/60.0) - (S/3600.0)
        else:
            lat = D + (M/60.0) + (S/3600.0)
        #-----------------------------------
        dms_list = lon_dms.split()  # split on white space
        D = np.float64( dms_list[0] )
        M = np.float64( dms_list[1] )
        S = np.float64( dms_list[2] )
        if (D < 0):
            lon = D - (M/60.0) - (S/3600.0)
        else:
            lon = D + (M/60.0) + (S/3600.0)

        #---------------------------------
        # Format the lat and lon strings
        #---------------------------------
        lon_str = '{x:.5f}'.format(x=lon)
        lat_str = '{x:.5f}'.format(x=lat)
        
        #-------------------------------------
        # Compare lon/lat to rfc_lon/rfc_lat
        #-----------------------------------------------------------
        # Note:  rfc_lon & rfc_lat are only available for the 8170
        # basins listed in usgsMasterMeta.txt, but lon & lat from
        # DMS are available for the 10512 basins in the HADS file.
        #-----------------------------------------------------------
        # 1601 USGS IDs in the HADS file do not occur in the
        # USGS_NWIS_all (streams only), and appear to be other
        # site types like atmosphere, canal, lake, and well.
        # See: get_station_type_dict() in usgs_utils.py. 
        #-----------------------------------------------------------
        # See:  https://en.wikipedia.org/wiki/Decimal_degrees
        #-----------------------------------------------------------
        ## tol = 0.0001  # difference over 10 meters
        tol = 0.001   # difference over 100 meters
        if (rfc_lat != '-999'):
            lat_diff = np.abs(lat - np.float64(rfc_lat))
            if (lat_diff > tol):
                lat_diff_count += 1
        if (rfc_lon != '-999'):
            lon_diff = np.abs(lon - np.float64(rfc_lon))
            if (lon_diff > tol):
                lon_diff_count += 1        

        #-------------------------------------
        # Get hydrograph type classification
        # (not available for all stations)
        #-------------------------------------
        try:
            htype_info = htype_map[ usgs_id ]
            htype = htype_info[ 'htype' ]
        except:
            htype = 'unknown'

        #------------------------------          
        # Get the station type string
        #------------------------------
        try:
            station_type = station_type_dict[ usgs_id ] 
        except:
            station_type = 'unknown'
          
        #-----------------------------
        # Write line to new TSV file
        #-----------------------------
        line = ''
        line += usgs_id      + delim
        line += loc_name     + delim
        line += station_type + delim
        line += lon_str      + delim
        line += lat_str      + delim
        line += rfc_name     + delim
        line += nws_loc_id   + delim
        line += nws_hsa_id   + delim
        line += goes_id      + delim
        line += huc8_code    + delim
        line += htype        + '\n'   # add newline at end
        tsv_unit.write( line )       
        k += 1

    #------------------ 
    # Close the files
    #------------------
    hads_unit.close()
    tsv_unit.close()
    print('Finished creating new TSV file for RFC basins.')
    print('  lon_diff_count =', lon_diff_count)
    print('  lat_diff_count =', lat_diff_count)
    print()
                
#   create_tsv()
#---------------------------------------------------------------------
def sort_by_site_code(tsv_file1='new_rfc_hads.tsv',
                      tsv_file2='new_rfc_hads_sorted.tsv'):

    rfc_dir   = get_rfc_data_dir()
    tsv_path1 = rfc_dir + tsv_file1
    tsv_path2 = rfc_dir + tsv_file2
    
    tsv_unit1 = open(tsv_path1, 'r')
    lines = tsv_unit1.readlines()
    tsv_unit1.close()
        
    #--------------------------------
    # Sort all lines (by "USGS_ID")
    #--------------------------------
    header = lines[0]   # header line  
    lines2 = lines[1:]    
    lines2.sort()
    tsv_unit2 = open(tsv_path2, 'w')
    tsv_unit2.writelines( [header] + lines2 )
    tsv_unit2.close()

#   sort_by_site_code()
#---------------------------------------------------------------------

