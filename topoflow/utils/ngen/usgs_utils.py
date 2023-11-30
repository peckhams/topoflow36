
# Copyright (c) 2023, Scott D. Peckham
#
# Oct 2023. Added get_station_type_dict().
#           Moved more functions here from collate_basins.py:
#           haversine, distance_on_sphere, get_closest_station.
#           Added get_hlr_outlet_info & get_closest_hlr_outlet.
# Sep 2023. Added get_usgs_station_info_dict() and
#           usgs_state_code_map() to map HLR basins to stations.
# Jul 2023. Utils to work with USGS station info.
#           Moved several functions here from collate_basins.py.
#           These are now called from collate_basins.py.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import usgs_utils as ut
#  >>> ut.get_station_subset()
#  >>> coords = ut.get_usgs_station_coords( NWIS_ALL=False )
#  >>> coords = ut.get_usgs_station_coords( NWIS_ALL=True )
#  >>> coords = ut.compare_usgs_station_coords( tol=1e-6 )
#
#---------------------------------------------------------------------
#
#  get_repo_dir()
#  get_usgs_dir()
#  get_usgs_filepath()
#  get_usgs_header_lines()
#  get_n_records()
#
#  get_station_subset()       ###
#  get_station_type_dict()    ###
#
#  get_usgs_station_coords()
#  compare_usgs_station_coords()
#
#  usgs_state_code_map()
#  get_usgs_station_info_dict()
#
#  replace_periods()
#  replace_punctuation()
#  replace_st()
#  replace_ft()
#  replace_abbreviations()
#  replace_state_letters()
#  fix_spellings()
#  replace_city_letters()
#  expand_usgs_site_name()   ######## or use "station_name" ??
#  get_usgs_long_name()
#
#  get_usgs_site_url()
#  get_usgs_huc_url()
#
#  get_usgs_missing_lat_lon()
#  get_usgs_missing_state()
#
#  haversine()            # moved here from collate_basins.py
#  distance_on_sphere()
#  get_closest_station()
#  get_hlr_outlet_info()
#  get_closest_hlr_outlet()
#
#---------------------------------------------------------------------

import numpy as np
import os, os.path
from topoflow.utils.ngen import hlr_tools as hlr

import pickle    # to save usgs station info map
import requests  # for get_usgs_missing_lat_lon, get_usgs_missing_state
import time

# import re  # for regular expressions. not used now

# from osgeo import ogr, osr
# import json, sys, time

#---------------------------------------------------------------------
def get_repo_dir():

    repo_dir = '/Users/peckhams/Dropbox/NOAA_NextGen/'
    repo_dir += '__NextGen_Example_Basin_Repo/'
    return repo_dir

#   get_repo_dir()
#---------------------------------------------------------------------
def get_usgs_dir( NWIS_ALL=False ):

    repo_dir = get_repo_dir()
    if (NWIS_ALL):
        usgs_dir = repo_dir + 'USGS_NWIS/'
    else:
        usgs_dir = repo_dir + 'USGS_Gauged_Basins/'   
    return usgs_dir

#   get_usgs_dir()
#---------------------------------------------------------------------
def get_usgs_filepath( NWIS_ALL=False ):

    usgs_dir = get_usgs_dir( NWIS_ALL=NWIS_ALL)
    if (NWIS_ALL):
        usgs_file = 'USGS_NWIS_stream_stations_usgs.tsv'
        ## usgs_file = 'USGS_NWIS_stream_stations_all.tsv'
    else:
        usgs_file = 'USGS_stream_gauge_data.tsv'
    return (usgs_dir + usgs_file)

#   get_usgs_filepath()
#---------------------------------------------------------------------
def get_usgs_header_lines( NWIS_ALL=False ):

    if (NWIS_ALL):
        nh_lines = 1
    else:
        nh_lines = 37
    return nh_lines

#   get_usgs_header_lines()
#---------------------------------------------------------------------
def get_n_records(type='usgs_streams'):

    if (type == 'usgs_streams'):
        #-------------------------------------------------
        # Any w/ monitoring location type of 'stream'
        # & station_id starts with 'USGS-'
        #-------------------------------------------------
        # Others in Alaska start with 'CAX' or 'USFWS'.
        # Some start with USAID (e.g. Afghanistan, Iraq)
        # Some start with RQ (e.g. Honduras)
        # Others start with 2-letter state code:
        #  AR, FL, IN, OR, TX
        # Others start with USARS, USBLM,
        #   USCE (USACE, USCOE?), USEPA, USFS, USFWS.
        #-------------------------------------------------
        return 145375
    elif (type == 'streams'):
        # Any w/ monitoring location type of 'stream'
        return 152926
    elif (type == 'all'):
        # All NWIS records: streams, wells, springs, etc.   
        return 1739528
    
#   get_n_records()
#---------------------------------------------------------------------
def get_station_subset(loc_type='usgs_stream',
                       in_tsv_file='USGS_NWIS_stations.tsv', 
                       out_tsv_file='USGS_NWIS_stream_stations.tsv'):

    #------------------------------------------------------------
    # Note: USGS_NWIS_stations.tsv contains info for 1,739,528
    #       basins, which have different "station types", such
    #       as Stream, Atmosphere, Lake and Well and different
    #       responsible agencies.
    #------------------------------------------------------------
    #       USGS_NWIS_stream_stations.tsv contains info for
    #       152,926 USGS stations of type "Stream".
    #------------------------------------------------------------    
    delim = '\t'  # tab character
    print('Working...')
    in_count  = 0
    out_count = 0
    
    #------------------------------
    # Open input and output files
    #------------------------------
    usgs_dir     = get_usgs_dir( NWIS_ALL=True )
    in_tsv_path  = usgs_dir + in_tsv_file
    out_tsv_path = usgs_dir + out_tsv_file
    in_tsv_unit  = open(in_tsv_path,  'r')
    out_tsv_unit = open(out_tsv_path, 'w')

    #--------------------------    
    # Process the header line
    #--------------------------
    header = in_tsv_unit.readline()
    out_tsv_unit.write( header )
    
    while (True):
        line = in_tsv_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        values = line.split( delim )
        in_count  += 1

        station_id   = values[2].strip()
        station_type = values[4].strip()
        
        #---------------------------------------        
        # Subset by "monitoring location type"
        # See notes in get_n_records().
        #---------------------------------------
        if (loc_type == 'usgs_stream'):
            if (station_type.lower() == 'stream'):
                if (station_id.startswith('USGS-')):
                    out_tsv_unit.write( line )
                    out_count += 1        
        elif (loc_type == 'any_stream'):
            if (station_type.lower() == 'stream'):
                out_tsv_unit.write( line )
                out_count += 1
        elif (loc_type == 'usaid_stream'):
            if (station_type.lower() == 'stream'):
                if (station_id.startswith('USAID-')):
                    out_tsv_unit.write( line )
                    out_count += 1         
#         elif (loc_type == 'any_well'):
#             if (station_type.lower() == 'well'):
#                 out_tsv_unit.write( line )
#                 out_count += 1
#         elif (loc_type == 'any_spring'):
#             if (station_type.lower() == 'spring'):
#                 out_tsv_unit.write( line )
#                 out_count += 1
                                
        if ((in_count % 10000) == 0):
            print('  in_count =', in_count)

    #-------------------------------   
    # Close input and output files
    #-------------------------------
    in_tsv_unit.close()
    out_tsv_unit.close()
    print('Total in  count =', in_count)
    print('Total out count =', out_count)
    print('Finished.')
    
#   get_station_subset()
#---------------------------------------------------------------------
def get_station_type_dict( in_tsv_file='USGS_NWIS_stations.tsv' ): 

    #------------------------------------------------------------
    # Note: USGS_NWIS_stations.tsv contains info for 1,739,528
    #       basins, which have different "station types", such
    #       as Stream, Atmosphere, Lake and Well.
    #------------------------------------------------------------
    usgs_dir  = get_usgs_dir( NWIS_ALL=True )
    dict_file = 'USGS_NWIS_station_type_dict.pkl'
    file_path = usgs_dir + dict_file
    if (os.path.exists( file_path )):
        print('Reading saved USGS-NWIS station type dictionary...')
        file_unit = open( file_path, 'rb')
        station_type_dict = pickle.load( file_unit )
        file_unit.close()
        print('Finished.')
        print()
        return station_type_dict

    delim = '\t'  # tab character
    print('Working...')
    in_count = 0
    usgs_station_count = 0
 
    #------------------
    # Open input file
    #------------------
    in_tsv_path = usgs_dir + in_tsv_file
    in_tsv_unit = open(in_tsv_path,  'r')

    #-----------------------    
    # Skip the header line
    #-----------------------
    header = in_tsv_unit.readline()
    station_type_dict = dict()
    station_type_list = list()

    while (True):
        line = in_tsv_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        values = line.split( delim )
        in_count  += 1

        station_id   = values[2].strip()
        station_type = values[4].strip()
        
        if (station_id.startswith('USGS-')):
            sid = station_id.replace('USGS-', '')
            station_type_dict[ sid ] = station_type
            usgs_station_count += 1

        #----------------------------------------------------         
        # E.g. atmosphere, canal, lake, spring stream, well
        #---------------------------------------------------- 
        if (station_type not in station_type_list):
            station_type_list.append( station_type )
                                
        if ((in_count % 10000) == 0):
            print('  in_count =', in_count)

    print('Unique station types =')
    print( station_type_list )

    SAVE_TO_FILE = True
    if (SAVE_TO_FILE):
        file_unit = open( file_path, 'wb' )
        pickle.dump( station_type_dict, file_unit, protocol=pickle.HIGHEST_PROTOCOL)
        file_unit.close()
        print('Saved USGS-NWIS station type dictionary to file:')
        print('  ' + file_path)
        print()
        
    #-------------------------------   
    # Close input and output files
    #-------------------------------
    in_tsv_unit.close()
    print('Total in count =', in_count)
    print('Number of USGS stations =', usgs_station_count )
    print('Finished.')
    
    return station_type_dict

#   get_station_type_dict()
#---------------------------------------------------------------------
def get_usgs_station_coords( delim='\t', NWIS_ALL=False,
                             SAVE_TO_FILE=False):

    #--------------------------------------------------------------
    # Note: Even though the NWIS_ALL spreadsheet has many more
    #       stations, the "USGS Gauged Basins" spreadsheet
    #       contains 16 station IDs that are not in the
    #       NWIS_ALL spreadsheet for station type "Stream"
    #       because they have some other station type.
    #       These IDs are:
    # 02203536, 04095090, 06146500, 06647000, 06780500, 08366400,
    # 08385503, 09196960, 09442958, 09527660, 10010500, 11204680,
    # 12485500, 13080000, 14246900, 14313600
    #--------------------------------------------------------------
    # These "missing" station IDs all seem to have valid URLs:
    # https://waterdata.usgs.gov/nwis/inventory/?site_no=02203536
    #    Tidal stream site
    # https://waterdata.usgs.gov/nwis/inventory/?site_no=04095090
    #    Ditch site
    # https://waterdata.usgs.gov/nwis/inventory/?site_no=06146500
    #    Ditch site
    # https://waterdata.usgs.gov/nwis/inventory/?site_no=06647000
    #    Canal site
    # https://waterdata.usgs.gov/nwis/inventory/?site_no=06780500
    #    Canal site
    # https://waterdata.usgs.gov/nwis/inventory/?site_no=08366400
    #    Canal site
    # https://waterdata.usgs.gov/nwis/inventory/?site_no=08385503
    #    Canal site
    #--------------------------------------------------------------
    usgs_dir = get_usgs_dir( NWIS_ALL=NWIS_ALL )
    # Note: file_name is same, but directories are different.
    file_path = usgs_dir + 'USGS_station_coords.npy'
    if (os.path.exists( file_path )):
        print('Reading saved USGS station coords...')
        station_coords = np.load( file_path )
        print('Finished.')
        print()
        return station_coords
      
    #----------------------------------------------------
    # Note: Use set of all USGS gauged basins for this.
    #----------------------------------------------------
    print('Getting USGS station coords...')
    usgs_path = get_usgs_filepath( NWIS_ALL=NWIS_ALL )
    usgs_unit = open( usgs_path, 'r' )

    #-----------------------
    # Skip over the header
    #-----------------------
    nh_lines = get_usgs_header_lines( NWIS_ALL=NWIS_ALL )
    for j in range( nh_lines ):
        line = usgs_unit.readline()

    #-------------------------------------    
    # Construct the station_coords array
    #-----------------------------------------------
    # Numpy Unicode string array w/ up to 16 chars
    # Every element is initialized as null string
    # station_id, lon, lat, huc_id (for NWIS_ALL)
    #-----------------------------------------------
    k = 0
    if (NWIS_ALL):
        n_stations = 145375
        station_coords = np.zeros((n_stations,4), dtype='<U16')
    else:
        n_stations = 24520
        station_coords = np.zeros((n_stations,3), dtype='<U16')
    n_bad_lats = 0
    n_bad_lons = 0

    while (True):
        usgs_line = usgs_unit.readline()
        if (usgs_line == ''):
            break  # (reached end of file)
        vals = usgs_line.split( delim )
        if (NWIS_ALL):
            station_id  = vals[2].strip()
            station_id  = station_id.replace('USGS-', '')
            station_huc = vals[6].strip()
            station_lat = vals[11].strip()  # lat, then lon in file
            station_lon = vals[12].strip()        
        else:
            station_id  = vals[1].strip()
            station_lat = vals[4].strip()  # lat, then lon in file
            station_lon = vals[5].strip()

        #------------------------------------------------
        # This was triggered 236 times for NWIS_ALL and
        # 58 times for the other "USGS_gauged" option.
        #------------------------------------------------        
        if (station_lat == '') or (station_lon == ''):
            # This occurs for Okinawa, Guam, Palau, FSM, etc.
            ########## datum is not used yet ##########
            lon, lat, datum = get_usgs_missing_lat_lon( station_id )
            station_lat = lat
            station_lon = lon
            n_bad_lats += 1
            n_bad_lons += 1

        station_coords[k,0] = station_id
        station_coords[k,1] = station_lon
        station_coords[k,2] = station_lat
        if (NWIS_ALL):
            station_coords[k,3] = station_huc
        k += 1

    if (n_bad_lats > 0) or (n_bad_lons > 0):
        print('Number of missing latitudes  =', n_bad_lats)
        print('Number of missing longitudes =', n_bad_lons)
        print('These were fixed by scraping values from the')
        print('USGS station URL.  This occurs for:')
        print('  Federal States of Micronesia (FSM), Guam,')
        print('  Okinawa, Palau, etc.')
        print()
   
    if (NWIS_ALL):
        #-----------------------------
        # Need to sort by station id
        #-----------------------------
        station_ids = station_coords[:,0]
        w = np.argsort( station_ids )
        station_coords = station_coords[w,:]
 
    if (SAVE_TO_FILE):
        usgs_dir  = get_usgs_dir( NWIS_ALL=NWIS_ALL )
        file_path = usgs_dir + 'USGS_station_coords.npy'
        np.save( file_path, station_coords)
        print('Saved USGS station coords to file:')
        print('  ' + file_path)
        print()

    TEST_IF_NUMERIC = False
    if (TEST_IF_NUMERIC):
        try:
            lons = np.float64( station_coords[:,1] )
        except:
            print('ERROR: Non-numeric longitudes found.')
        try:
            lats = np.float64( station_coords[:,2] )
        except:
            print('ERROR: Non-numeric latitudes found.')
         
    return station_coords

#   get_usgs_station_coords()
#---------------------------------------------------------------------
def compare_usgs_station_coords( tol=0.001 ):

    #-------------------------------------------------
    # Note: tol=0.01 corresponds to 1111 m or less
    #          n_lon_diffs = 3, n_lat_diffs = 2
    #       tol=0.001 corresponds to 111 m or less
    #          n_lon_diffs = 13, n_lat_diffs = 10
    #       tol=1e-4 corresponds to 11.1 m or less
    #          n_lon_diffs = 84, n_lat_diffs = 72  
    #       tol=1e-5 corresponds to 1.11 m or less
    #          n_lon_diffs = 105, n_lat_diffs = 104
    #       tol=1e-6 corresponds to 0.11 m or less
    #          n_lon_diffs = 106, n_lat_diffs = 106
    #       tol=1e-8 corresponds to 0.0011 m or less
    #          n_lon_diffs = 106, n_lat_diffs = 106
    #-------------------------------------------------
    # See:  Wikipedia: Decimal degrees
    #-------------------------------------------------
    coords1 = get_usgs_station_coords( NWIS_ALL=True )
    ids1 = coords1[:,0]
    w1 = np.argsort( ids1 )
    coords1 = coords1[w1,:]  # sort coord1 by station id string
    #----------------------------------------    
    coords2 = get_usgs_station_coords( NWIS_ALL=False )
    ids2 = coords2[:,0]
    w2 = np.argsort( ids2 )
    coords2 = coords2[w2,:]  # sort coord2 by station id string 
    #-----------------------------------------
    print('Comparing lons & lats for 2 USGS datasets...')
    print('   coords1.shape =', coords1.shape)  
    print('   coords2.shape =', coords2.shape)
    n_stations1 = 145375  # for NWIS_ALL option
    n_stations2 = 24520
    k = 0
    j = 0
    n_id_matches = 0
    n_lon_diffs  = 0
    n_lat_diffs  = 0

    while ((k < n_stations1) and (j < n_stations2)):
        usgs_id1 = coords1[k,0]
        usgs_id2 = coords2[j,0]
#         print('usgs_id1 =', usgs_id1)
#         print('usgs_id2 =', usgs_id2)
        if (usgs_id1 == usgs_id2):
            n_id_matches += 1
            usgs_lon1 = np.float64( coords1[k,1] )
            usgs_lat1 = np.float64( coords1[k,2] )
            usgs_lon2 = np.float64( coords2[j,1] )
            usgs_lat2 = np.float64( coords2[j,2] )
            if (np.abs(usgs_lon1 - usgs_lon2) > tol):
                print('WARNING for USGS ID =', usgs_id1)
                print('  lon1 =', usgs_lon1, ', lon2 =', usgs_lon2)
                n_lon_diffs += 1
            if (np.abs(usgs_lat1 - usgs_lat2) > tol):
                print('WARNING for USGS ID =', usgs_id1)
                print('  lat1 =', usgs_lat1, ', lat2 =', usgs_lat2)
                n_lat_diffs += 1
            k += 1
            j += 1
        elif (usgs_id1 < usgs_id2):
            k += 1
        elif (usgs_id1 > usgs_id2):
            j += 1

    print('id1 count    =', k)
    print('id2 count    =', j)
    print('tolerance    =', tol)
    print('n_id_matches =', n_id_matches)
    print('n_lon_diffs  =', n_lon_diffs)
    print('n_lat_diffs  =', n_lat_diffs)
    print('Finished.')
    print()

#   compare_usgs_station_coords()
#---------------------------------------------------------------------
def usgs_state_code_map():

    #-------------------------------------------------------    
    # Note: These are FIPS codes.
    #       FIPS = Federal Information Processing Standard
    # https://en.wikipedia.org/wiki/Federal_Information
    #         Processing_Standard_state_code
    #-------------------------------------------------------
    map = {
    0:  ['--', 'Not in the USA.'],
    1:  ['AL', 'Alabama'],
    2:  ['AK', 'Alaska'],
    60: ['AS', 'American Samoa'],
    3:  ['AS', 'American Samoa (FIPS 5-1)'],
    4:  ['AZ', 'Arizona'],
    5:  ['AR', 'Arkansas'],
    81: ['BI', 'Baker Island'],
    6:  ['CA', 'California'],
    7:  ['CZ', 'Canal Zone (FIPS 5-1)'],
    8:  ['CO', 'Colorado'],
    9:  ['CT', 'Connecticut'],
    10: ['DE', 'Delaware'],
    11: ['DC', 'District of Columbia'],
    12: ['FL', 'Florida'],
    64: ['FM', 'Federal States of Micronesia'],
    13: ['GA', 'Georgia'],
    14: ['GU', 'Guam (FIPS 5-1)'],
    66: ['GU', 'Guam'],
    15: ['HI', 'Hawaii'],
    84: ['HI', 'Howard Island'],
    16: ['ID', 'Idaho'],
    17: ['IL', 'Illinois'],
    18: ['IN', 'Indiana'],
    19: ['IA', 'Iowa'],
    86: ['JI', 'Jarvis Island'],
    67: ['JA', 'Johnston Atoll'],
    20: ['KS', 'Kansas'],
    21: ['KY', 'Kentucky'],
    89: ['KR', 'Kingman Reef'],
    22: ['LA', 'Louisiana'],
    23: ['ME', 'Maine'],
    68: ['MH', 'Marshall Islands'],
    24: ['MD', 'Maryland'],
    25: ['MA', 'Massachusetts'],
    26: ['MI', 'Michigan'],
    71: ['MI', 'Midway Islands'],
    27: ['MN', 'Minnesota'],
    28: ['MS', 'Mississippi'],
    29: ['MO', 'Missouri'],
    30: ['MT', 'Montana'],
    76: ['NI', 'Navassa Island'],
    31: ['NE', 'Nebraska'],
    32: ['NV', 'Nevada'],
    33: ['NH', 'New Hampshire'],
    34: ['NJ', 'New Jersey'],
    35: ['NM', 'New Mexico'],
    36: ['NY', 'New York'],
    37: ['NC', 'North Carolina'],
    38: ['ND', 'North Dakota'],
    69: ['MP', 'Northern Mariana Islands'],
    39: ['OH', 'Ohio'],
    40: ['OK', 'Oklahoma'],
    41: ['OR', 'Oregon'],
    70: ['PW', 'Palau'],
    95: ['PA', 'Palmyra Atoll'],
    42: ['PA', 'Pennsylvania'],
    43: ['PR', 'Puerto Rico (FIPS 5-1)'],
    72: ['PR', 'Puerto Rico'],
    44: ['RI', 'Rhode Island'],
    45: ['SC', 'South Carolina'],
    46: ['SD', 'South Dakota'],
    47: ['TN', 'Tennessee'],
    48: ['TX', 'Texas'],
    74: ['UM', 'U.S. Minor Outlying Islands'],
    49: ['UT', 'Utah'],
    50: ['VT', 'Vermont'],
    51: ['VA', 'Virginia'],
    52: ['VI', 'Virgin Islands (FIPS 5-1)'],
    78: ['VI', 'Virgin Islands'],
    79: ['WI', 'Wake Island'],
    53: ['WA', 'Washington'],
    54: ['WV', 'West Virginia'],
    55: ['WI', 'Wisconsin'],
    56: ['WY', 'Wyoming'],
    83: ['MX', 'In Mexico'],
    84: ['MX', 'In Mexico'],
    86: ['MX', 'In Mexico'],
    90: ['NB', 'New Brunswick, Canada'],
    91: ['QC', 'Quebec, Canada'],
    92: ['ON', 'Ontario, Canada'],
    93: ['MB', 'Manitoba, Canada'],
    94: ['SK', 'Saskatchewan, Canada'],
    95: ['AB', 'Alberta, Canada'],
    96: ['BC', 'British Columbia, Canada'],
    97: ['YT', 'Yukon Territory, Canada'] }
    
    return map

#   usgs_state_code_map()
#---------------------------------------------------------------------
def get_usgs_station_info_dict():

    usgs_dir = get_usgs_dir( NWIS_ALL=True )
    station_info_file = 'USGS_station_info_dict.pkl'
    file_path = usgs_dir + station_info_file
    if (os.path.exists( file_path )):
        print('Reading saved USGS station info dictionary...')
        file_unit = open( file_path, 'rb')
        station_info = pickle.load( file_unit )
        file_unit.close()
        print('Finished.')
        print()
        return station_info
      
    #----------------------------------------------------
    # Note: Use set of all USGS gauged basins for this.
    #----------------------------------------------------
    print('Getting USGS station info as dictionary...')
    usgs_path = get_usgs_filepath( NWIS_ALL=True )
    usgs_unit = open( usgs_path, 'r' )
    delim = '\t'

    #-----------------------
    # Skip over the header
    #-----------------------
    nh_lines = get_usgs_header_lines( NWIS_ALL=True )
    for j in range( nh_lines ):
        line = usgs_unit.readline()

    #-------------------------------------------------------------    
    # Construct the station_info dictionary, where USGS gauge id
    # is the key used to get a dictionary with info for a gauge
    #-------------------------------------------------------------
    # Using NWIS_ALL, so there will be 145,375 basins.
    #-------------------------------------------------------------
    sc_map = usgs_state_code_map()   
    station_info = dict()
    k = 0
    n_bad_lats = 0
    n_bad_lons = 0

    while (True):
        usgs_line = usgs_unit.readline()
        if (usgs_line == ''):
            break  # (reached end of file)
        vals = usgs_line.split( delim )

        sid   = vals[2].strip()
        sid   = sid.replace('USGS-', '')
        name  = vals[3].strip()
        # stype = vals[4].strip()  # always "Stream"
        huc   = vals[6].strip()
        lat   = vals[11].strip()  # lat, then lon in file
        lon   = vals[12].strip()
        fips  = vals[25].strip()  # state FIPS code        
        sc    = sc_map[ int(fips) ][0]  # 2-letter code
        url   = get_usgs_site_url( sid )

        #------------------------------------------------
        # This was triggered 236 times for NWIS_ALL and
        # 58 times for the other "USGS_gauged" option.
        #------------------------------------------------        
        if (lat == '') or (lon == ''):
            # This occurs for Okinawa, Guam, Palau, FSM, etc.
            ########## datum is not used yet ##########
            lon, lat, datum = get_usgs_missing_lat_lon( sid )
            n_bad_lats += 1
            n_bad_lons += 1

        station_info[ sid ] = \
            {'name':name, 'huc':huc, 'lat':lat, 'lon':lon,
             'state':sc, 'url':url }
        k += 1

    if (n_bad_lats > 0) or (n_bad_lons > 0):
        print('Number of missing latitudes  =', n_bad_lats)
        print('Number of missing longitudes =', n_bad_lons)
        print('These were fixed by scraping values from the')
        print('USGS station URL.  This occurs for:')
        print('  Federal States of Micronesia (FSM), Guam,')
        print('  Okinawa, Palau, etc.')
        print()
   
    SAVE_TO_FILE = True
    if (SAVE_TO_FILE):
        usgs_dir  = get_usgs_dir( NWIS_ALL=True )
        file_path = usgs_dir + station_info_file
        file_unit = open( file_path, 'wb' )
        pickle.dump( station_info, file_unit, protocol=pickle.HIGHEST_PROTOCOL)
        file_unit.close()
        print('Saved USGS station info dictionary to file:')
        print('  ' + file_path)
        print()

#     TEST_IF_NUMERIC = False
#     if (TEST_IF_NUMERIC):
#         try:
#             lons = np.float64( station_coords[:,1] )
#         except:
#             print('ERROR: Non-numeric longitudes found.')
#         try:
#             lats = np.float64( station_coords[:,2] )
#         except:
#             print('ERROR: Non-numeric latitudes found.')
         
    return station_info

#   get_usgs_station_info_dict()
#---------------------------------------------------------------------
#---------------------------------------------------------------------
def replace_periods( name, new_char='' ):

    #----------------------------------------------------------
    # Note: This function replaces periods "." in a string
    #       except when they are within a decimal number.
    #       It uses "finditer" in the re module.
    #----------------------------------------------------------
    # e.g. SALT RIVER NEAR CHRYSOTILE, ARIZ. MILE 34.8 -> 348
    # if we use: name = name.replace('.','')
    #----------------------------------------------------------
    # Something like this could also work:
    # re.sub('\.(?!\s|\d|$)', '', text)
    # return text
    #----------------------------------------------------------           
    char_list  = list( name )
    char_list2 = list()
    n = len(char_list)

    for k in range(n):
        char = char_list[k]
        if (char == '.'):
            char_L = char_list[k-1] if (k > 0)   else ''
            char_R = char_list[k+1] if (k < n-1) else ''
            if not(char_L.isnumeric() and char_R.isnumeric()):
               char = new_char
        char_list2.append( char )
                     
    name = "".join(char_list2)
    return name
    
#   replace_periods()
#---------------------------------------------------------------------
def replace_punctuation( name ):

    name = name.replace(', ', ',')
    name = name.replace(',', ' ')
    name = name.replace('; ', ';')
    name = name.replace(';', ' ')
    name = replace_periods( name, '')  # preserve decimals
    name = name.replace('@', 'AT')
    # SAUK RIVER AB WHITE CHUCK RIVER NR  DARRINGTON, WA
    name = name.replace('  ', ' ')
    
    return name
    
#   replace_punctuation()     
#---------------------------------------------------------------------
def replace_st( name ):

    #----------------------------------------------------------
    # Note: The abbreviation "ST" can be 'SAINT', 'STREET',
    #       or 'STATE'.  All occurrences of "ST" for "STATE"
    #       are: "ST HW", "ST HWY", "ST HIGHWAY" or "ST_LINE".
    #       Cases w/ 'SAINT' are more rare that 'STREET',
    #       so we use a list gleaned from USGS gauged basins.
    #----------------------------------------------------------
    name = name.replace(' ST LINE ', ' STATE LINE ')
    name = name.replace(' ST HWY ',  ' STATE HIGHWAY ')
    name = name.replace(' ST HW ',   ' STATE HIGHWAY ')
            
    saint_list = [
    'ALBAN', 'ALBANS',
    'ANTHONY', 'ANTHONYS',
    'AUGUST', 'AUGUSTINE',
    'BENEDICT', 'BENEDICTS',
    'CLEMENT', 'CLEMENTS',
    'CROIX',
    'DAVID', 'DAVIDS',
    'FRANCIS',
    'GEORGE', 'GEORGES',
    'HELENA', 'HELENS',
    'JOE',
    'JOHN', 'JOHNS',
    'JOHNSBURY',
    'JOSEPH', 'JOSEPHS',
    'LANDRY',
    'LEONARD', 'LEONARDS',
    'LOUIS',
    'MARIES',
    'MARY',    'MARYS',
    'MATTHEW', 'MATTHEWS',
    'MICHAEL', 'MICHAELS',
    'PARIS',
    'PAUL',    'PAULS',
    'SEBASTIAN', 'SEBASTIANS',
    'STEPHEN', 'STEPHENS',
    'THOMAS' ]

    #------------------------------------
    # Works also if string starts w/ ST
    #------------------------------------
    for saint in saint_list:
        saint_abbr = 'ST ' + saint
        if (saint_abbr in name):
            saint_name = 'SAINT ' + saint
#             if ('MICHAEL' in name):
#                 print('saint_abbr =', saint_abbr)
#                 print('saint_name =', saint_name)
            name = name.replace(saint_abbr, saint_name )

    #-------------------------------    
    # Anything else must be street
    #-------------------------------
    if (' ST ' in name):
        name = name.replace(' ST ', ' STREET ')

    return name

#   replace_st()
#---------------------------------------------------------------------
def replace_ft( name ):

    #--------------------------------------------------------
    # For these, FT is always preceeded by numeric string.
    #--------------------------------------------------------
    # ASSABET RIVER 200 FT BELOW RT 85 AT HUDSON, MA
    # HUNT R 250 FT DS FRY BRK AT FRENCHTOWN, RI
    # QUEEN R 1400 FT UPSTR WM REYNOLDS RD AT EXETER RI
    # TEMECULA C 500 FT DS OF VAIL DAM NR TEMECULA CA
    # MOHIHI STREAM AT ALT 3420 FT NR WAIMEA,KAUAI,HI
    # KOAIE STREAM AT ALT 3770 FT NEAR WAIMEA,KAUAI,HI
    # Waialae Str at alt 3,820 ft nr Waimea, Kauai, HI
    # WAIALAE STREAM AT ALT 800 FT NR WAIMEA,KAUAI,HI
    # Halaulani Str at alt 400 ft nr Kilauea, Kauai, HI
    # HANALEI RIVER AT ALT 625 FT NR HANALEI,KAUAI,HI
    # Kaupuni Str at alt 374 ft nr Waianae, Oahu, HI
    # Luluku Str at alt 220 ft nr Kaneohe, Oahu, HI
    # WAIHEE STREAM AT ALT 260 FT NEAR HEEIA, OAHU, HI
    # WAIAHOLE STR AT ALT 250 FT WAIAHOLE,OAHU,HI
    # Waikane Str at alt 75 ft at Waikane, Oahu, HI
    #--------------------------------------------------------
    # SPIRIT CK .35 mi DS OF McCOYS CRK AT FT GORDON,GA
    # FT DRUM CREEK AT SUNSHINE ST PKY NEAR FT DRUM, FL
    # MIDDLE RIVER CANAL AT S-36 NR FT LAUDERDALE, FL
    # NORTH NEW RIVER CA AT 20 MI BEND NR FT LAUD, FL
    # LAKE OUTFALL, HENDRY CREEK AT GLADIOLUS DR FT MYER
    # LAKE OUTFALL TO HENDRY CREEK SUMMERLIN RD FT MYERS
    # JACKS BRANCH AT CR 78 NR FT DENAUD FL
    # SHELL CREEK AT CIRCLE K GROVES NEAR FT OGDEN FL
    # LT MANATEE R AT TAYLOR-GILL RD NEAR FT LONESOME FL
    # CHATTAHOOCHEE R .36 MI DS WFG DAM NR FT GAINES, GA
    # TURTLE CREEK NR FT WALTON BEACH, FL
    # Smith River near Ft Logan MT
    # SF L WIND RIV AB WASHAKIE RES, NR FT WASHAKIE, W
    # PORCUPINE CREEK NR FT YATES, N
    # Tyler Falls at Ft Niobrara NWR nr Valentine, Nebr.
    # CACHE LA POUDRE RIV AT MO OF CN, NR FT COLLINS, CO
    # QUARRY C AT MISSOURI R, FT LEAVENWORTH, KS
    # Neosho River blw Ft Gibson Lake nr Ft Gibson, OK
    # W Fk Trinity Rv at Lk Worth Dam abv Ft Worth, TX
    # Clear Fk Brazos Rv at Ft Griffin, TX
    # TRINCHERA C AB TURNERS RANCH, NR FT GARLAND, CO.
    # SAND GATE DIV FROM FT SUMNER CANAL AT FT SUMNER,NM
    # Limpia Ck abv Ft Davis, TX
    # Coyanosa Draw nr Ft Stockton, TX
    # HUACHUCA CAN NR FT HUACHUCA ARIZ
    # WARNER C NR FT BRAGG CA
    #--------------------------------------------------------
    name  = name.upper()
    if not(' FT ' in name):
        return name  # for speed

    parts = name.split(' FT ')
    words = parts[0].split(' ')
    last_word = words[-1]
    if (last_word.isnumeric()):
        name = name.replace(' FT ',   ' FEET ')
    else:
        name = name.replace(' FT ',   ' FORT ')

    return name

#   replace_ft()
#---------------------------------------------------------------------
def replace_abbreviations( name ):

    #----------------------------------------------
    # Call replace_periods() before replace_st().
    #----------------------------------------------
    name = replace_st( name ) # STATE, SAINT, OR STREET
    name = replace_ft( name ) # FEET, FORT

    #---------------------------------------------------    
    # Investigate these:
    # ' FM ' (28 matches)
    # ' GL ' (5 matches), Maybe Glacier ??
    # BREWERY CREEK-UPSTREAM SITE-AT CROSS PLAINS, WI
    # Red Bluff Ck at FM 1283 nr Pipe Creek, TX
    #---------------------------------------------------
    
    #---------------------------------------------
    # e.g. COMPTON C A 120TH ST NR COMPTON CA
    #      ARCADIA WASH A GRAND AVE A ARCADIA CA
    #---------------------------------------------
    name = name.replace(' A ',    ' AT ')
    name = name.replace(' AB ',   ' ABOVE ')
    name = name.replace(' ABV ',  ' ABOVE ')
    name = name.replace(' ALT ',  ' ALTITUDE ') # e.g. ALT 250 FT
    #--------------------------------------------
    name = name.replace(' BL ',   ' BELOW ')
    ## No match for 'BLO' in NWIS, but many in RFC lists.
    ## e.g. "St. Vrain Ck. blo Longmont"
    name = name.replace(' BLO ',  ' BELOW ')   
    name = name.replace(' BLW ',  ' BELOW ')
    name = name.replace(' DS ',   ' DOWNSTREAM ')  ## CHECK ALL
    # name = name.replace(' MI ',   ' MILES ')  # Or Michigan ???
    name = name.replace(' NR ',   ' NEAR ')
    # https://waterdata.usgs.gov/wa/nwis/current?type=rivermi
    # RM followed by decimal number
    # M Branch Shepherd Creek at RM 0.04 near Mt Airy OH
    name = name.replace(' RM ',   ' RIVER MILE ')  # occurs 5 times
    name = name.replace(' UPSTR ',' UPSTREAM ')
            
    #---------------------------------------------
    # NOTE: These all start w/o a space.
    #--------------------------------------------- 
    name = name.replace('EF ',   'EAST FORK ')
    name = name.replace('MF ',   'MIDDLE FORK ')
    name = name.replace('NF ',   'NORTH FORK ')
    name = name.replace('SF ',   'SOUTH FORK ')
    name = name.replace('WF ',   'WEST FORK ')
    #---------------------------------------------
    name = name.replace('E F ',  'EAST FORK ')
    name = name.replace('M F ',  'MIDDLE FORK ')
    name = name.replace('N F ',  'NORTH FORK ')
    name = name.replace('S F ',  'SOUTH FORK ')
    name = name.replace('W F ',  'WEST FORK ')
    #---------------------------------------------
    name = name.replace('E FK ', 'EAST FORK ')
    name = name.replace('M FK ', 'MIDDLE FORK ')
    name = name.replace('N FK ', 'NORTH FORK ')
    name = name.replace('S FK ', 'SOUTH FORK ')
    name = name.replace('W FK ', 'WEST FORK ')
    #---------------------------------------------
    name = name.replace('EFK ',  'EAST FORK ')
    name = name.replace('NFK ',  'NORTH FORK ')
    name = name.replace('SFK ',  'SOUTH FORK ')
    name = name.replace('WFK ',  'WEST FORK ')
    #---------------------------------------------
    name = name.replace('EFORK ', ' EAST FORK ')
    name = name.replace('NFORK ', ' NORTH FORK ')
    name = name.replace('SFORK ', ' SOUTH FORK ')
    name = name.replace('WFORK ', ' WEST FORK ')  # starts as "W.FORK"
    #---------------------------------------------
    name = name.replace(' FK ',   ' FORK ')

    #-----------------------------------   
    # More water feature abbreviations   
    #---------------------------------------------                    
    # 'NB POTOMAC RIVER AT STEYER MD'
    # But next 2 also contain "NB" vs. "NB ":
    # 'GREENBRIER RIVER AT DURBIN, WV'
    # 'DRIFTWOOD RIVER NEAR EDINBURGH IND'
    #---------------------------------------------
    name = name.replace(' B ',   ' BRANCH ')
    name = name.replace(' BR ',  ' BRANCH ')
    name = name.replace('NB ', 'NORTH BRANCH ')
    name = name.replace('SB ', 'SOUTH BRANCH ')             
    #---------------------------------------------
    name = name.replace(' BYU ',   ' BAYOU ')      ## CHECK MORE
    #---------------------------------------------    
    name = name.replace(' CANL ',  ' CANAL ')
    name = name.replace(' CONF ',  ' CONFLUENCE ')
    name = name.replace(' CONFL ', ' CONFLUENCE ') 
    #---------------------------------------------
    name = name.replace(' C ',   ' CREEK ')
    name = name.replace(' CK ',  ' CREEK ')
    name = name.replace(' CR ',  ' CREEK ')
    name = name.replace(' CRK ', ' CREEK ')
    #---------------------------------------------
    # HUACHUCA CAN NR FT HUACHUCA ARIZ
    # HUACHUCA CANYON NEAR FORT HUACHUCA, AZ
    # BEAR RIVER DUCK CLUB CAN NR BEAR RIVER CITY, UT   # CAN = CANAL, checked
    #---------------------------------------------
    name = name.replace(' CAN ',  ' CANYON ')   #### (9 matches; CAN BE CANAL)
    name = name.replace(' CN ',   ' CANYON ')   #### (90 matches; CAN BE CANAL ??)
    name = name.replace(' CNYN ', ' CANYON ') 
    name = name.replace(' CYN ',  ' CANYON ')
    #---------------------------------------------
    # BIG SAND CL AB ST DITCH NR BADGER BASIN, WY
    # MALAD R BL BRN DUCK CL CAN NR BEAR RIVER CITY, UT
    # JW SOUTHERN PP A SND BAR DIV DAM NR LNG BRN CA

    # MALAD R BL BRN DUCK CL CAN NR BEAR RIVER CITY, UT
    #---------------------------------------------
    # name = name.replace(' BRN ',  ' BROWN ')   # 2 matches
    # name = name.replace(' CL ',   ' COULEE ')  # 2 matches
    # name = name.replace(' SND ',  ' SAND ')    # 2 matches
    # name = name.replace(' SND BAR ',  ' SAND BAR ')    # 2 matches
    #---------------------------------------------
    # US could be "upstream" or for US highway
    #--------------------------------------------- 
    name = name.replace(' DIV ',  ' DIVERSION ')
    name = name.replace(' DTCH ', ' DITCH ')     
    #--------------------------------------------- 
    name = name.replace(' INTK ',  ' INTAKE ')
    name = name.replace(' INTKE ', ' INTAKE ')     
    #---------------------------------------------
    name = name.replace(' L ',     ' LITTLE ') 
    name = name.replace('LTL ',    'LITTLE ')  # 16 matches, no start space
    #---------------------------------------------
    name = name.replace(' LK ',     ' LAKE ')
    name = name.replace(' MTH ',    ' MOUTH ')  # 5 matches
    name = name.replace(' PD ',     ' POND ')    
    #---------------------------------------------                  
    name = name.replace(' R ',     ' RIVER ')
    name = name.replace(' RV ',    ' RIVER ')
    name = name.replace(' RIV ',   ' RIVER ')
    name = name.replace(' RVR ',   ' RIVER ')   # In RFC names
    #---------------------------------------------
    name = name.replace(' RES ',   ' RESERVOIR ')
    name = name.replace(' RESV ',  ' RESERVOIR ')
    #--------------------------------------------- 
    # NF SKOKOMISH R BL STAIRCASE RPDS NR HOODSPORT, WA
    name = name.replace(' RPDS ',   ' RAPIDS ')
    #--------------------------------------------- 
    name = name.replace(' SPGS ',   ' SPRINGS ')
    name = name.replace(' SPRGS ',  ' SPRINGS ')
    name = name.replace(' SPRNGS ', ' SPRINGS ')
    name = name.replace(' SPS ',    ' SPRINGS ')  ##### CHECK
    #---------------------------------------------     
    # 'HAIKU STR NEAR HEEIA OAHU HI'
    name = name.replace(' STR ', ' STREAM ')
    name = name.replace(' TRIB ', ' TRIBUTARY ')
  
    #-----------------------------
    # Longer place abbreviations
    #------------------------------------------------------------------------
    # name = name.replace('CRIR ', 'Colorado River Indian Reservation')
    # name = name.replace('C R I R ', 'Colorado River Indian Reservation')
    name = name.replace(' O R N L ', ' ORNL ')
    name = name.replace(' OAK RIDGE NATL LAB ', ' ORNL ')
    # name = name.replace('PVID ', 'Palo Verde Irrigation District)
    # name = name.replace('P V I D ', 'Palo Verde Irrigation District)
    # 18 matches for "(TVA)", all at end of string.
    # name = name.replace(' (TVA) ', ' (Tennessee Valley Authority) ')
    # YNP = Yellowstone National Park.
    # YNP spans 3 states: WY, MT, ID
    #------------------------------------------------------------------------

    #-------------------------------
    # Abbreviations at end of name
    #-----------------------------------
    # TVA = Tennessee Valley Authority
    #-----------------------------------
    if (name.endswith('TENN (TVA)')):
        name = name[:-10] + '(TVA) TN'
    if (name.endswith('TN (TVA)')):
        name = name[:-8] + '(TVA) TN'
    #-----------------------------------
    if (name.endswith('TENN (CE)')):
        name = name[:-9] + '(CE) TN'
    if (name.endswith('TN (CE)')):
        name = name[:-7] + '(CE) TN'
    #-----------------------------------
                       
    #---------------------------------
    # Abbreviations at start of name
    #---------------------------------    
    # 'G MIAMI RIVER AT HAMILTON OH'
    if (name[:2] == 'G '):
        name = 'GREAT ' + name[2:]
    # 'L BEAVER CREEK NEAR EAST LIVERPOOL OH' 
    if (name[:2] == 'L '):
        name = 'LITTLE ' + name[2:]
    if (name[:3] == 'ST '):
        name = 'SAINT ' + name[3:]    ### see replace_st() function
    #---------------------------------
    if (name[:2] == 'E '):
        name = 'EAST ' + name[2:]
    if (name[:2] == 'N '):
        name = 'NORTH ' + name[2:]
    if (name[:2] == 'S '):
        name = 'SOUTH ' + name[2:]
    if (name[:2] == 'W '):
        name = 'WEST ' + name[2:]
    #--------------------------------
    if (name[:3] == 'SO '):
        name = 'SOUTH ' + name[3:]
    if (name[:3] == 'NO '):
        name = 'NORTH ' + name[3:] 
    #--------------------------------
    if (name[:5] == 'TRIB '):
        name = 'TRIBUTARY ' + name[5:]

    #------------------------------       
    # Other feature abbreviations           
    #------------------------------
    name = name.replace(' BUTE ',   ' BUTTE ')
    name = name.replace(' ISL ',    ' ISLAND ')   # e.g. Twitchell Isl      
    name = name.replace(' MT ',     ' MOUNT ')    # e.g. MT CLEMENS
    name = name.replace(' MTN ',    ' MOUNTAIN ')
    name = name.replace(' PT ',   ' POINT ')
    name = name.replace(' RH ',   ' RANCH ')
    name = name.replace(' RCH ',  ' RANCH ')
    name = name.replace(' RNCH ', ' RANCH ')
    name = name.replace(' RNGR ', ' RANGER ')
    # name = name.replace()' SR ', ' STATE ROAD ')  ## ??????????????
    ########  See replace_st() function   ########
    name = name.replace(' STA ',  ' STATION ')  ## "RANGER STA"
    name = name.replace(' STRM SWR ',  ' STORM SEWER ')
    # RAPID CREEK ABOVE WRF NR RAPID CITY, SD
    name = name.replace(' WRF ',  ' WHARF ')  ####### CHECK
    ## SEWAGE TREATMENT PL  (PL = plant or place)    

    #-----------------------------    
    # Road-related abbreviations
    #-----------------------------
    name = name.replace(' AVE ',  ' AVENUE ')
    name = name.replace(' BLVD ',  ' BOULEVARD ')
    name = name.replace(' BRG ',  ' BRIDGE ')
    # In USGS_gauged as "Rd Crsg"; In CAMELS as "Rd"
    name = name.replace(' CRSG ', ' CROSSING ')  # 3 matches
    name = name.replace(' CTR ',  ' CENTER ')
    name = name.replace(' DR ',   ' DRIVE ')  ###### CHECK
    name = name.replace(' HWY ',  ' HIGHWAY ')
    name = name.replace(' HIGHWY ', ' HIGHWAY ')
    ## name = name.replace(' INTK ',  ??????
    name = name.replace(' JCT ',  ' JUNCTION ')
    name = name.replace(' LN ',   ' LANE ')
    name = name.replace(' PKWY ', ' PARKWAY ')
    name = name.replace(' RD ',   ' ROAD ')
    name = name.replace(' RT ',   ' ROUTE ')
    name = name.replace(' RTE ',  ' ROUTE ')
    #------------------------------------------------------
    # These are now in replace_st() function.
    ## name = name.replace(' ST LINE ', ' STATE LINE ')
    ## name = name.replace(' ST HWY ',  ' STATE HIGHWAY ')
    ## name = name.replace(' ST HW ',   ' STATE HIGHWAY ')
    #------------------------------------------------------
    # CA AQUEDUCT A N PORTAL TEHACHAPI TNL NR GORMAN CA
    name = name.replace(' TNL ',  ' TUNNEL ')  # JUST ONE
    name = name.replace(' TUN ',  ' TUNNEL ')  ### CHECK ALL

    return name

#   replace_abbreviations()
#---------------------------------------------------------------------
def replace_state_letters( name ):

    #-------------------------------------------
    # What about these endings:
    # Az-Ca, Az-Cal, Az-Ca, -Ariz-Calif, Ca-Az
    #-------------------------------------------
#     if (name[-2:] == ' A'):
#         name = name[:-2] + ' AL'   
    #--------------------------------
    if (name[-4:] == ' ALA'):
        name = name[:-4] + ' AL'
    if (name[-4:] == ' ARK'):     # don't match "PARK"
        name = name[:-4] + ' AR'
    if (name[-4:] == ' CAL'):
        name = name[:-4] + ' CA'
    if (name[-4:] == ' FLA'):
        name = name[:-4] + ' FL'
    if (name[-4:] == ' ILL'):
        name = name[:-4] + ' IL'
    if (name[-4:] == ' IND'):
        name = name[:-4] + ' IN'
    if (name[-4:] == ' N C'):
        name = name[:-4] + ' NC'
    if (name[-4:] == ' N D'):
        name = name[:-4] + ' ND'
    if (name[-4:] == ' N M'):
        name = name[:-4] + ' NM'
    if (name[-4:] == ' N Y'):
        name = name[:-4] + ' NY'
    if (name[-4:] == ' S C'):
        name = name[:-4] + ' SC'
    if (name[-4:] == ' S D'):
        name = name[:-4] + ' SD'
    if (name[-4:] == ' TEN'):
        name = name[:-4] + ' TN'
    if (name[-4:] == ' TEX'):
        name = name[:-4] + ' TX'
    if (name[-4:] == ' WYO'):
        name = name[:-4] + ' WY'
    #--------------------------------
    if (name[-5:] == ' ARIZ'):
        name = name[:-5] + ' AZ'
    if (name[-5:] == ' COLO'):
        name = name[:-5] + ' CO'
    if (name[-5:] == ' CONN'):
        name = name[:-5] + ' CT'
    if (name[-5:] == ' MASS'):
        name = name[:-5] + ' MA'
    if (name[-5:] == ' MICH'):
        name = name[:-5] + ' MI'
    if (name[-5:] == ' MINN'):
        name = name[:-5] + ' MN'
    if (name[-5:] == ' MISS'):
        name = name[:-5] + ' MS'
    if (name[-5:] == ' NEBR'):
        name = name[:-5] + ' NE'
    if (name[-5:] == ' OREG'):
        name = name[:-5] + ' OR'
    if (name[-5:] == ' TENN'):
        name = name[:-5] + ' TN'
    if (name[-5:] == ' UTAH'):
        name = name[:-5] + ' UT'
    if (name[-5:] == ' WASH'):
        name = name[:-5] + ' WA'
    if (name[-5:] == ' W VA'):
        name = name[:-5] + ' WV'
    #--------------------------------
    if (name[-6:] == ' AL-GA'):
        name = name[:-6] + ' AL/GA'   ##########
    if (name[-6:] == ' AZ-CA'):
        name = name[:-6] + ' AZ/CA'   ##########
    if (name[-6:] == ' CA-AZ'):
        name = name[:-6] + ' AZ/CA'   ##########
    if (name[-6:] == ' IL-KY'):
        name = name[:-6] + ' IL/KY'   ##########
    if (name[-6:] == ' KY-TN'):
        name = name[:-6] + ' KY/TN'   ##########    
    #--------------------------------    
    if (name[-6:] == ' CALIF'):
        name = name[:-6] + ' CA'
    if (name[-6:] == ' IDAHO'):
        name = name[:-6] + ' ID'
    if (name[-6:] == ' MAINE'):
        name = name[:-6] + ' ME'
    if (name[-6:] == ' N MEX'):
        name = name[:-6] + ' NM'
    if (name[-6:] == ' TEXAS'):
        name = name[:-6] + ' TX'
    #--------------------------------
    if (name[-7:] == ' NEVADA'):
        name = name[:-7] + ' NV'
    #--------------------------------
    if (name[-8:] == ' ALABAMA'):
        name = name[:-8] + ' AL'
    if (name[-8:] == ' FLORIDA'):
        name = name[:-8] + ' FL'
    if (name[-8:] == ' GEORGIA'):
        name = name[:-8] + ' GA'
    if (name[-8:] == ' MIAMI F'):
        name = name[:-8] + ' MIAMI FL'
    #--------------------------------
    if (name[-9:] == ' ARKANSAS'):
        name = name[:-9] + ' AR'
    if (name[-9:] == ' KENTUCKY'):
        name = name[:-9] + ' KY'
    if (name[-9:] == ' MARYLAND'):
        name = name[:-9] + ' MD'
    if (name[-9:] == ' MISSOURI'):
        name = name[:-9] + ' MO'
    #--------------------------------
    if (name[-10:] == ' LOUISIANA'):
        name = name[:-10] + ' LA'
    if (name[-10:] == ' TENNESSEE'):
        name = name[:-10] + ' TN'
    #--------------------------------
    if (name[-14:] == ' NEW HAMPSHIRE'):
        name = name[:-14] + ' NH'
    #------------------------------
    # These may only be caused by
    # replacing "N" w/ NORTH
    # Moved that to the end.
    #------------------------------ 
#             if (name[-8:] == ' NORTH C'):
#                 name = name[:-8] + ' NC'
#             if (name[-8:] == ' SOUTH C'):
#                 name = name[:-8] + ' SC'
#             # Comes from "N MEX"
#             if (name[-10:] == ' NORTH MEX'):
#                 name = name[:-10] + 'NM'

    return name
    
#   replace_state_letters()
#---------------------------------------------------------------------
def fix_spellings( name ):

    #-------------------------------------------------
    # Note: Call this after "replace_state_letters".
    #-------------------------------------------------
    if (name == 'ARROYO VALLE BELOW LANG CANAL NEAR LIVERMORE CA'):
        name = name.replace(' CANAL ', ' CANYON ')
        ###################### CHECK ABOVE ###############
    if ('BOGUECHITTO' in name):
        name = name.replace('BOGUECHITTO', 'BOGUE CHITTO')
    if (name == 'BOGUE CHITTO NEAR BUSH LA'):
        name = 'BOGUE CHITTO RIVER NEAR BUSH LA'
    if ('CHYSOTILE' in name):
        # This occurs in CAMELS data set.
        name = name.replace('CHYSOTILE', 'CHRYSOTILE')
    if ('D ALENE' in name):
        # Note: "D ALENE" is used in USGS_gauged.
        name = name.replace('D ALENE', "D'ALENE")
    if ('ELEVENPOINT' in name):
        name = name.replace('ELEVENPOINT', 'ELEVEN POINT')
    if ('HAZLEGREEN' in name):
        name = name.replace('HAZLEGREEN', 'HAZELGREEN')
    if ('HOBOLOCHITTO' in name):
        name = name.replace('HOBOLOCHITTO', 'HOBOLO CHITTO')
    if (name == 'HORSE CREEK NEAR ARCADIA FL'):
        name = 'HORSE CREEK AT SR 72 NEAR ARCADIA FL'
    if (name == 'KINCHAFOONEE CREEK NEAR DAWSON GA'):
        name = 'KINCHAFOONEE CREEK AT PINEWOOD ROAD NEAR DAWSON GA'
    if (name == 'LAMAR RIVER NEAR TOWER RANGER STATION YNP'):
        name = 'LAMAR RIVER NEAR TOWER FALLS RANGER STATION YNP'
    if (name == 'LITTLE CYPRESS CREEK NEAR JEFFERSON TX'):
        name == 'LITTLE CYPRESS BAYOU NEAR JEFFERSON TX' # (in USGS gauged)
    if (name == 'MANATEE RIVER NEAR MYAKKA HEAD FL'):
        name = 'MANATEE RIVER AT SR 64 NEAR MYAKKA HEAD FL'
    if ('MARKLEEVILLECA' in name):
        name = name.replace('MARKLEEVILLECA', 'MARKLEEVILLE CA')
    if (name == 'MCDONALDS BRANCH IN LEBANON STATE FOREST NJ'):
        name = 'MCDONALDS BRANCH IN BYRNE STATE FOREST NJ'
        ###################### CHECK ABOVE ###############
    if (name == 'MINAM RIVER NEAR MINAM OR'):
        name = 'MINAM RIVER AT MINAM OR'
    if (name == 'MUDDY CREEK BELOW CLEAR CREEK NEAR COUGAR WA'):
        name = 'MUDDY RIVER BELOW CLEAR CREEK NEAR COUGAR WA'
    if ('NORTH EAST' in name):
        name = name.replace('NORTH EAST', 'NORTHEAST')
    if ('NORTH WEST' in name):
        name = name.replace('NORTH WEST', 'NORTHWEST')
        # careful for: "South West City, MO"
        ###################### CHECK ABOVE ###############
    if (name == 'NORTH FORK CLEARWATER RIVER NEAR CANYON RANGER STATION'):
        name += ' ID'
    if (name == 'PEACE RIVER AT ARCADIA FL'):
        name = 'PEACE RIVER AT SR 70 AT ARCADIA FL'
    if (name.startswith('RD CHUTE BYU')):
        name = 'RED CHUTE BAYOU TRIBUTARY AT SKNL P2 NEAR SLIGO LA'
        ###################### CHECK ABOVE ###############
    if (name == 'SAUK RIVER ABOVE WHITECHUCK RIVER NEAR DARRINGTON WA'):
        name = 'SAUK RIVER ABOVE WHITE CHUCK RIVER NEAR DARRINGTON WA'
    if (name == 'SHIAWASSEE RIVER AT BYRON MI'):
        name = 'SHIAWASSEE RIVER AT BATH ROAD AT BYRON MI'
    # Skannal Ponds in Skannal Plantation, near Sligo, LA
    # ESTATE OF JOHN C. SKANNAL (also Sligo Plantation)
    # Also misspelled as: Skannel
    if (' SKNL P2' in name):
        name = name.replace(' SKNL P2', ' SKANNEL POND #2')
    if (' SKNL PD' in name):
        name = name.replace(' SKNL PD', ' SKANNEL POND')
    if (' SKNL POND' in name):
        name = name.replace(' SKNL POND', ' SKANNEL POND')
    if ('TCHEFUNCTA' in name):
        name = name.replace('TCHEFUNCTA', 'TCHEFUNCTE')
    if ('WHISKEY CHITTO' in name):
        name = name.replace('WHISKEY', 'OUISKA')
     
    return name
       
#   fix_spellings()
#---------------------------------------------------------------------
def replace_city_letters( name ):

    if (name[-4:] == ' SLC'):
        name = name[:-4] + ' SALT LAKE CITY UT'
    if (name[-8:] == ' S L CITY'):
        name = name[:-8] + ' SALT LAKE CITY UT'
    #---------------------------------------------
    if (' SLGO ' in name):
        # Sligo, Louisiana
        name = name.replace(' SLGO ', ' SLIGO ')
    if (' SNTA BRB ' in name):
        name = name.replace(' SNTA BRB ', ' SANTA BARBARA ')
    return name
 
#   replace_city_letters()
#---------------------------------------------------------------------
def expand_usgs_site_name( name, UPPER=True,
                       REPLACE_PUNCTUATION=True,
                       REPLACE_ABBREVIATIONS=True,
                       FIX_SPELLINGS=True,
                       REPLACE_STATE_LETTERS=True,
                       REPLACE_CITY_LETTERS=True,
                       FIX_MISSING_STATE=True,
                       FIX_DIRECTIONS=True):

    #---------------------------------------------------------
    # Note: USGS station names contain many abbreviations,
    #       some standard and some unusual (e.g. SKNL).
    #       This function (and get_usgs_long_name) expands
    #       most of the abbrevations, corrects misspellings,
    #       adds missing state code, etc.
    #       Abbreviations like "ST" and "FT" can mean more
    #       than one thing: see replace_st, and replace_ft.
    #---------------------------------------------------------
    # What about:
    # SUGAR C TR A LK FRSME NR ARCADIA 
    #---------------------------------------     
    if (UPPER):
        name = name.upper()
    #-------------------------------------------
    # e.g. "CLINTON RIVER NEAR FRASER, MICH.               ."
    if (name.endswith('.')):
        name = name[:-1].strip()
    # "N.F. CLEARWATER RIVER NR CANYON RANGER STATION," (MOPEX)
    if (name.endswith(',')):
        name = name[:-1].strip()
    #-------------------------------------------
    if (REPLACE_PUNCTUATION):  ##### MAYBE MOVE THIS TOWARD END?
        name = replace_punctuation( name )
    #-------------------------------------------
    if (REPLACE_ABBREVIATIONS):
        name = replace_abbreviations( name )

    #---------------------------------------
    # Fix non-2-letter state abbreviations
    #---------------------------------------
    if (REPLACE_STATE_LETTERS):
        name = replace_state_letters( name )
 
    #-------------------
    # Fix misspellings
    #-------------------
    if (FIX_SPELLINGS):
        name = fix_spellings( name )
                                                                    
    if (REPLACE_CITY_LETTERS):
        name = replace_city_letters( name )
  
    if (FIX_MISSING_STATE):
        if (name == 'GALENA CREEK AT GALENA CREEK STATE PARK'):
            name += ' NV'
        #-------------------------------------------
        # There are over 600 more w/ missing state
        #-------------------------------------------
        # YNP = Yellowstone National Park.
        # YNP spans 3 states: WY, MT, ID
#             if (name[-4:] == ' YNP')
#                 name = name + ' WY'

    if (FIX_DIRECTIONS):
        # Moved these last because otherwise get:
        # 'MORA RIVER NEAR SHOEMAKER NORTH MEX'
        # 'VALLEY RIVER AT TOMOTLA NORTH C'
        name = name.replace(' N ',   ' NORTH ')
        name = name.replace(' S ',   ' SOUTH ')
        name = name.replace(' E ',   ' EAST ')
        name = name.replace(' W ',   ' WEST ')
        #---------------------------------------------
        name = name.replace(' NE ',  ' NORTHEAST ')
        name = name.replace(' NW ',  ' NORTHWEST ')
        name = name.replace(' SE ',  ' SOUTHEAST ')
        name = name.replace(' SW ',  ' SOUTHWEST ')
        if (name[:3] == 'NE '):
            name = 'NORTHEAST ' + name[3:]
        if (name[:3] == 'NW '):
            name = 'NORTHWEST ' + name[3:]
        if (name[:3] == 'SE '):
            name = 'SOUTHEAST ' + name[3:]
        if (name[:3] == 'SW '):
            name = 'SOUTHWEST ' + name[3:]
        
    return name

#   expand_usgs_site_name()
#---------------------------------------------------------------------
def get_usgs_long_name( usgs_name ):

    # Note: YT = Yukon Territory, Int'l Boundary, AK-Canada
    state_list = [
    ' AL', ' AK', ' AZ', ' AR', ' CA', ' CO', ' CT', ' DE', ' FL',
    ' GA', ' HI', ' ID', ' IL', ' IN', ' IA', ' KS', ' KY', ' LA',
    ' ME', ' MD', ' MA', ' MI', ' MN', ' MS', ' MO', ' MT', ' NE',
    ' NV', ' NH', ' NJ', ' NM', ' NY', ' NC', ' ND', ' OH', ' OK',
    ' OR', ' PA', ' PR', ' RI', ' SC', ' SD', ' TN', ' TX', ' UT',
    ' VT', ' VA', ' WA', ' WV', ' WI', ' WY', 'YT' ]
    
    long_name  = expand_usgs_site_name( usgs_name )
    long_name  = long_name.title()
    if (long_name[-3:].upper() in state_list):
        #----------------------------------
        # Capitalize 2-letter state codes
        #----------------------------------
        state_code = long_name[-2:].upper()
        long_name  = long_name[:-2] + state_code

    if (long_name[-4:].upper() == ' YNP'):
        long_name = long_name[:-3] + 'YNP'  # (Yellowstone spans 3 states)

    return long_name
    
#   get_usgs_long_name()
#---------------------------------------------------------------------
def get_usgs_site_url( site_no ):

    #--------------------------------------------------
    # Note: site_no is an 8-digit USGS station code.
    # e.g., 01206900.  Each has its own web page.
    #--------------------------------------------------
    # Note: The ID: 01095505 occurs in GAGES2_all but
    #       not in USGS_gauged, along w/ 119 others.
    #       But their URLs seem OK.
    #       See compare_basin_ids().
    #--------------------------------------------------    
    url = 'http://waterdata.usgs.gov/nwis/nwisman/?site_no='
    url += str(site_no)
    return url
 
#   get_usgs_site_url()
#---------------------------------------------------------------------
def get_usgs_huc_url( huc_no ):

    #-------------------------------------------------------------
    # Legacy HUC numbers have 8 digits, new ones have 12.
    # e.g. https://water.usgs.gov/lookup/getwatershed?10250001
    #-------------------------------------------------------------
    # HUCs subdivide water resource regions of the United States 
    # from Regions (2-digit) to Subregions (4-digit),
    # Basins (6-digit) to Subbasins (8-digit),
    # Watersheds (10-digit) to Subwatersheds (12-digit) -
    # through the use of the hydrologic unit codes (HUCS).
    # This applies to all of the 22 U.S. regions, which includes
    # Alaska (19), Hawaii (20), Puerto Rico (21) and Guam (22).
    #-------------------------------------------------------------
    url = 'https://water.usgs.gov/lookup/getwatershed?'
    url += str(huc_no)
    return url
    
#   get_usgs_huc_url()
#---------------------------------------------------------------------
def get_usgs_missing_lat_lon( site_id, REPORT=False ):

    #---------------------------------------------
    # Scrape lat and lon from the USGS site page
    #---------------------------------------------
    url = get_usgs_site_url( site_id )
    page_bytes = requests.get(url).content 
    page_str   = page_bytes.decode()
    lines      = page_str.split('\n')

    # Example line:
    # <dd>Latitude  26&#176;21'49", &nbsp; Longitude -127&#176;45'23" &nbsp; NAD27<br /></dd>    
    for line in lines:
        if ('Latitude' in line):
           break
    words    = line.split()  # on any whitespace (even 2 spaces)
    #--------------------------------------------
    lat_word = words[1]
    lat_word = lat_word.replace('&#176;',' ')  # degree symbol
    lat_word = lat_word.replace("'", ' ')
    lat_word = lat_word.replace('"', '')
    lat_word = lat_word.replace(',', '')
    dms_list = lat_word.split()
    D = np.float64( dms_list[0] )
    M = np.float64( dms_list[1] )
    S = np.float64( dms_list[2] )
    if (D < 0):
        lat = D - (M/60.0) - (S/3600.0)
    else:
        lat = D + (M/60.0) + (S/3600.0)    
    #---------------------------------------------
    lon_word = words[4]
    lon_word = lon_word.replace('&#176;',' ')
    lon_word = lon_word.replace("'", ' ')
    lon_word = lon_word.replace('"', '')
    dms_list = lon_word.split()
    D = np.float64( dms_list[0] )
    M = np.float64( dms_list[1] )
    S = np.float64( dms_list[2] )
    if (D < 0):
        lon = D - (M/60.0) - (S/3600.0)
    else:
        lon = D + (M/60.0) + (S/3600.0)      
    #---------------------------------------------
    datum = words[6]
    datum = datum[0:5]    
    if (REPORT):
        print('lon, lat, datum =', lon, ', ', lat, ', ', datum)
    return lon, lat, datum
        
#   get_usgs_missing_lat_lon()
#---------------------------------------------------------------------
def get_usgs_missing_state( site_id, code_map, REPORT=False ):

    #--------------------------------------------
    # Scrape state name from the USGS site page
    #--------------------------------------------
    # Doesn't work for these 6 station ids:
    # 01011500, 01129300, 06144250,
    # 06144350, 12300000, 12355000
    # These are on "International Boundary" and
    # in Canada.
    #--------------------------------------------
    if (site_id == '01011500'):
        return 'New Brunswick', 'NB'
    if (site_id == '01129300'):
        return 'Quebec', 'QC'
    if (site_id == '06144250'):
        return 'Alberta', 'AB'
    if (site_id == '06144350'):
        return 'Saskatchewan', 'SK'
    # These 2 may work w/ code_map now.
    if (site_id == '12300000'):
        return 'British Columbia', 'BC'
    if (site_id == '12355000'):
        return 'British Columbia', 'BC'     

    url = get_usgs_site_url( site_id )
    page_bytes = requests.get(url).content 
    page_str   = page_bytes.decode()
    lines      = page_str.split('\n')

    # Example line:
    # <dd>Yavapai County, Arizona,  Hydrologic Unit 15070102</dd>
    for line in lines:
        if ('hydrologic unit' in line.lower()):
           break
    words = line.split(',')
    state_name = words[1].strip()
    state_code = code_map[ state_name ] 
    if (REPORT):
        print('state_name, state_code =', state_name, ', ', state_code)
    return state_name, state_code
        
#   get_usgs_missing_state()
#---------------------------------------------------------------------
def haversine( lon1, lat1, lon2, lat2):

    # p = 0.017453292519943295
    p     = (np.pi / 180)
    term1 = np.cos((lat2 - lat1)*p)/2
    term2 = np.cos(lat1*p) * np.cos(lat2*p)
    term3 = (1-np.cos((lon2 - lon1)*p)) / 2
    haver = 0.5 - term1 + (term2 * term3) 
    return haver
    
#   haversine()
#---------------------------------------------------------------------
def distance_on_sphere(lon1, lat1, lon2, lat2):

    #----------------------------------
    # 12742 = diameter of Earth in km
    #--------------------------------------------------
    # https://stackoverflow.com/questions/41336756/
    #         find-the-closest-latitude-and-longitude
    # https://en.wikipedia.org/wiki/Haversine_formula
    #--------------------------------------------------
    haver = haversine( lon1, lat1, lon2, lat2)
    return 12742 * np.arcsin(np.sqrt( haver ))  # [km]
    
#   distance_on_sphere()
#---------------------------------------------------------------------
def get_closest_station(station_coords, lon, lat, REPORT=False):

    #------------------------------------------------------
    # Call usgs.get_usgs_station_coords() once to extract
    # USGS gauging station data from USGS TSV file.
    # station_data = [[id1,x1,y1],[id2,x2,y2],...]
    #------------------------------------------------------
    station_ids  = station_coords[:,0]
    station_lons = np.float64( station_coords[:,1] )
    station_lats = np.float64( station_coords[:,2] )

    lon2 = np.float64( lon )  # could be string
    lat2 = np.float64( lat )
    if ((lon2 <= -999) or (lat2 < -90.0) or (lat2 > 90.0)):
        closest_id  = 'unknown'
        closest_lon = 'unknown'
        closest_lat = 'unknown'
        dmin        = -999.0
        return closest_id, closest_lon, closest_lat, dmin
        
    dist_vals = distance_on_sphere(station_lons, station_lats, lon2, lat2)
    dmin = dist_vals.min()
    imin = np.argmin( dist_vals )
    
    closest_id  = station_ids[ imin ]
    closest_lon = station_lons[ imin ]
    closest_lat = station_lats[ imin ]

    if (REPORT):
        print('Closest USGS station ID =', closest_id)
        print('   Station longitude    =', closest_lon)
        print('   Station latitude     =', closest_lat)
        print('   Station distance     =', dmin, '[km]')
        print()

    return closest_id, closest_lon, closest_lat, dmin

#   get_closest_station()
#---------------------------------------------------------------------
def get_hlr_outlet_info():

    usgs_dir  = get_usgs_dir( NWIS_ALL=True )
    info_path = usgs_dir + 'USGS_HLR_outlet_info.npy'
    if (os.path.exists( info_path )):
        print('Reading saved HLR basin outlet info...')
        hlr_outlet_info = np.load( info_path )
        print('Finished.')
        print()
        return hlr_outlet_info
      
    #----------------------------------------------------
    # Note: Use the set of USGS HLR basins for which
    # we have been able to determine basin outlet info.
    # This is 9539 of the 43931 basins in HLR set.
    #----------------------------------------------------
    print('Getting HLR basin outlet info...')
    hlr_dir  = hlr.get_hlr_data_dir()
    hlr_file = 'new_hlr_na_all.tsv'
    hlr_path = hlr_dir + hlr_file
    hlr_unit = open( hlr_path, 'r' )
    hlr_line = hlr_unit.readline()   # skip one-line header
    n_stations = 9539
    
    #--------------------------------------    
    # Construct the hlr_outlet_info array
    #-----------------------------------------------
    # Numpy Unicode string array w/ up to 16 chars
    # Every element is initialized as null string
    # basin_id, lon, lat, hlr_code
    #-----------------------------------------------
    k = 0
    delim = '\t'  # tab character

    hlr_outlet_info = np.zeros((n_stations,4), dtype='<U16')
    n_bad_lats = 0
    n_bad_lons = 0

    while (True):
        hlr_line = hlr_unit.readline()
        if (hlr_line == ''):
            break  # (reached end of file)
        vals = hlr_line.split( delim )
        #--------------------------------
        hlr_id   = vals[0].strip()
        hlr_lon  = vals[3].strip()
        hlr_lat  = vals[4].strip()        
        hlr_code = vals[24].strip()

        #----------------------        
        # Check lon and lat ?
        #----------------------
        if (hlr_lat == '-999'):
            n_bad_lats += 1
        if (hlr_lon == '-999'):
            n_bad_lons += 1

        #-------------------------------------    
        # Save info in hlr_outlet_info array
        #-------------------------------------
        hlr_outlet_info[k,0] = hlr_id
        hlr_outlet_info[k,1] = hlr_lon
        hlr_outlet_info[k,2] = hlr_lat
        hlr_outlet_info[k,3] = hlr_code
        k += 1
   
    #----------------------------------
    # Sort by HLR outlet id ? (VALUE)
    #----------------------------------
#     hlr_ids = hlr_outlet_info[:,0]
#     w = np.argsort( hlr_ids )
#     hlr_outlet_info = hlr_outlet_info[w,:]

    #------------------------- 
    # Close the HLR TSV file
    #-------------------------
    hlr_unit.close()
    
    SAVE_TO_FILE = True
    if (SAVE_TO_FILE):
        np.save( info_path, hlr_outlet_info )
        print('Saved HLR outlet info to file:')
        print('  ' + info_path)
        print()
         
    return hlr_outlet_info
    
#   get_hlr_outlet_info()
#---------------------------------------------------------------------
def get_closest_hlr_outlet(usgs_id, usgs_station_info,
                           hlr_outlet_info, REPORT=False):

    #------------------------------------------------------------
    # This function uses usgs_station_info and hlr_outlet_info.
    # This function's caller should obtain these via functions:
    # get_usgs_station_info_dict() & get_hlr_outlet_info().
    #------------------------------------------------------------
    station_info = usgs_station_info[ usgs_id ]
    usgs_lon = station_info[ 'lon' ]
    usgs_lat = station_info[ 'lat' ]
    
    hlr_ids   = hlr_outlet_info[:,0]
    hlr_lons  = np.float64( hlr_outlet_info[:,1] )
    hlr_lats  = np.float64( hlr_outlet_info[:,2] )
    hlr_codes = np.float64( hlr_outlet_info[:,3] )
    
    lon2 = np.float64( usgs_lon )  # could be string
    lat2 = np.float64( usgs_lat )
    if ((lon2 <= -999) or (lat2 < -90.0) or (lat2 > 90.0)):
        closest_id  = 'unknown'
        closest_lon = 'unknown'
        closest_lat = 'unknown'
        dmin        = -999.0
        hlr_code    = -1
        return closest_id, closest_lon, closest_lat, dmin, hlr_code
        
    dist_vals = distance_on_sphere(hlr_lons, hlr_lats, lon2, lat2)
    dmin = dist_vals.min()
    imin = np.argmin( dist_vals )
    
    closest_id  = hlr_ids[ imin ]   # an HLR outlet ID ("value")
    closest_lon = hlr_lons[ imin ]
    closest_lat = hlr_lats[ imin ]
    hlr_code    = hlr_codes[ imin ]

    if (REPORT):
        print('Closest HLR outlet ID  =', closest_id)
        print('   Outlet longitude    =', closest_lon)
        print('   Outlet latitude     =', closest_lat)
        print('   Outlet distance     =', dmin, '[km]')
        print('   Basin HLR code      =', hlr_code)    # in {1,...,20}
        print()

    return closest_id, closest_lon, closest_lat, dmin, hlr_code

#   get_closest_hlr_outlet()
#---------------------------------------------------------------------






   