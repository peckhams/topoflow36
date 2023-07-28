
# Copyright (c) 2023, Scott D. Peckham
#
# Jun 2023. Utils to collate several river basin datasets.
# Jul 2023. Added haversine, distance_on_sphere,
#           get_closest_station, add_ars_basins, add_czo_basins,
#           add_hlr_basins, add_lter_basins, add_neon_basins,
#           and their subfunctions.
#           Moved all USGS-specific utils to usgs_utils.py.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import collate_basins as cb
#  >>> cb.compare_basin_ids( db_str2='CAMELS' )
#  >>> cb.compare_basin_ids( db_str2='MOPEX' )
#  >>> cb.compare_basin_ids( db_str2='GAGES2_ref' )
#  >>> cb.compare_basin_names( db_str2='CAMELS ')
#  >>> cb.compare_basin_names( db_str2='MOPEX ')
#  >>> cb.compare_basin_names( db_str2='GAGES2_ref ')
#  >>> cb.collate()
#
#---------------------------------------------------------------------
#
#  get_repo_dir()
#  get_data_directory()
#  get_data_filepath()
#  get_n_records()
#  get_id_column()
#  get_name_column()
#  get_minlon_column()
#  get_header_lines()
#  get_basin_ids()
#  compare_basin_ids()
#
#  get_basin_names()
#  compare_basin_names()
#
#  haversine()
#  distance_on_sphere()
#  get_closest_station()
#  
#  get_dataset_key_names()
#  get_next_line()
#  get_next_id()
#  get_bounding_box()
#
#  get_state_code_map()
#  get_state_code()
#
#  skip_header_lines()
#  get_headings_line()
#  get_headings()
#  write_new_header()
#
#  collate()            ##### In progress #####
#  get_new_usgs_line()
#
#  add_ars_basins()
#  get_ars_name()
#  get_ars_long_name()
#  get_ars_site_url()
#
#  add_czo_basins()
#  get_czo_id()
#  get_czo_long_name()
#  get_czo_site_url()
#
#  add_hlr_basins()
#  add_lter_basins()
#  add_neon_basins()
#
#  export_csv_to_tsv()
#
#---------------------------------------------------------------------
from topoflow.utils import usgs_utils as usgs

import numpy as np
import time

# import re  # for regular expressions, not used now
# from osgeo import ogr, osr
# import json, sys, time

#---------------------------------------------------------------------
def get_repo_dir():

    repo_dir = '/Users/peckhams/Dropbox/NOAA_NextGen/'
    repo_dir += '__NextGen_Example_Basin_Repo/'
    return repo_dir

#   get_repo_dir()
#---------------------------------------------------------------------
def get_data_directory( db_str='USGS_gauged' ):

    repo_dir = get_repo_dir()

    #---------------------------------------------    
    # Directories for each of the basin datasets
    #---------------------------------------------
    if (db_str == 'ARS'):
        data_dir = repo_dir + 'USDA_ARS/'
    elif (db_str == 'CAMELS'):
        data_dir = repo_dir + 'CAMELS/'
    elif (db_str == 'CZO'):
        data_dir = repo_dir + 'CZO/'
    elif (db_str == 'EPA'):
        data_dir = repo_dir + 'EPA/'
    elif (db_str == 'FPS'):
        data_dir = repo_dir + 'USGS_FPS/'
    elif (db_str == 'GAGES2_all'):
        data_dir = repo_dir + 'GAGES-II/'
    elif (db_str == 'GAGES2_ref'):
        data_dir = repo_dir + 'GAGES-II/'
    elif (db_str == 'HLR'):
        data_dir = repo_dir + 'HLR/hlrshape/'
    elif (db_str == 'LTER'):
        data_dir = repo_dir + 'LTER/'
    elif (db_str == 'MOPEX'):
        data_dir = repo_dir + 'MOPEX/'
        data_dir += 'Hydrologic Synthesis Project 2009_18GB/'
        data_dir += 'spatialdata/MOPEX431_basins/'
    elif (db_str == 'NEON'):
        data_dir = repo_dir + 'NEON/'
    elif (db_str == 'RFC'):    
        data_dir = repo_dir + 'NOAA_RFCs/'
    elif (db_str == 'USGS_gauged'):
        data_dir = repo_dir + 'USGS_Gauged_Basins/'
    else:
        print('SORRY: No match for:', db_str)
        data_dir = None

    return data_dir

#   get_data_directory()
#---------------------------------------------------------------------
def get_data_filepath( db_str='USGS_gauged' ):

    data_dir = get_data_directory( db_str )
    
    if (db_str == 'ARS'):
        file_path = data_dir + 'ars_basin_info.tsv'
    elif (db_str == 'CAMELS'):
        file_path = data_dir + 'new_camels.tsv'
    elif (db_str == 'CZO'):
        file_path = data_dir + 'czo_basin_info_sheet1.tsv'  #####
    elif (db_str == 'EPA'):
        file_path = data_dir + 'EPA_WQX_stations.tsv'
    elif (db_str == 'FPS'):
        file_path = data_dir + 'FPS_basin_info.tsv'
    elif (db_str == 'GAGES2_all'):
        file_path = data_dir + 'new_gages2_all.tsv'
    elif (db_str == 'GAGES2_ref'):
        file_path = data_dir + 'new_gages2_ref.tsv'
    elif (db_str == 'HLR'):
        file_path = data_dir + 'new_hlr_na_all_v2.csv'  ##########
    elif (db_str == 'LTER'):
        file_path = data_dir + 'new_lter_info.tsv'
    elif (db_str == 'MOPEX'):
        file_path = data_dir + 'new_mopex431_sorted.tsv'
    elif (db_str == 'NEON'):
        file_path = data_dir + 'NEON_Field_Site_Metadata_20230309.tsv'
    elif (db_str == 'RFC'):
        print('SORRY: TSV file is not ready yet.', db_str)
        file_path = None
        # file_path = data_dir + 'new_rfc_info.csv'
    elif (db_str == 'USGS_gauged'):
        file_path = data_dir + 'USGS_stream_gauge_data.tsv'
    else:
        print('SORRY: No match for:', db_str)
        file_path = None
       
    return file_path
 
#   get_data_filepath()
#---------------------------------------------------------------------
def get_n_records( db_str ):

    if (db_str == 'ARS'):
        n_records = 771
    elif (db_str == 'CAMELS'):
        n_records = 671
    elif (db_str == 'CZO'):
        n_records = 15  ### w/ currently available data
    elif (db_str == 'FPS'):
        n_records = 4756
    elif (db_str == 'GAGES2_all'):
        n_records = 9322
    elif (db_str == 'GAGES2_ref'):
        n_records = 2057
    elif (db_str == 'HLR'):
        n_records = 43391    # Correct number
        # n_records = 47479  # w/ repeated VALUEs
        # n_records = 9796   # in new_hlr_na_all_v2.csv
    elif (db_str == 'LTER'):
        n_records = 30
    elif (db_str == 'MOPEX'):
        n_records = 431
    elif (db_str == 'NEON'):
        n_records = 81
    elif (db_str == 'RFC'):
        print('SORRY: TSV file is not ready yet.', db_str)
        n_records = None
    elif (db_str == 'USGS_gauged'):
        n_records = 25420
    else:
        print('SORRY: No match for:', db_str)
        n_records = None

    return n_records

#   get_n_records()
#---------------------------------------------------------------------
def get_id_column( db_str ):

    #--------------------------------------------------------
    # ID Column Heading for each dataset:
    # ARS ID      = "Watershed ID"
    # CAMELS ID   = "gauge_id"
    # FPS ID      = "SiteNumber" (prepend 0 if 7 digits)
    # GAGES2_all  = "STAID"
    # GAGES2_ref  = "STAID"
    # HLR         = "VALUE"
    # LTER        = "Site Acronym" (e.g. "AND")
    # MOPEX       = "SiteCode"
    # NEON        = "field_site_id" (e.g. "ABBY")
    # RFC         = "id" (for Gages.csv and Catchments.csv)
    # USGS_gauged = "site_no"
    #--------------------------------------------------------
    id_col = None
    list1 = ['ARS', 'CAMELS', 'FPS', 'GAGES2_all',
             'GAGES2_ref', 'HLR', 'LTER', 'MOPEX']
    list2 = ['NEON', 'USGS_gauged']
    
    if (db_str in list1):
        id_col = 0
    elif (db_str in list2):
        id_col = 1
    else:
        #### CZO and RFC ####
        print('SORRY: No match for:', db_str)
        print('       Returning id_col = 0.')
        id_col = 0

    return id_col

#   get_id_column()
#---------------------------------------------------------------------
def get_name_column( db_str ):

    #--------------------------------------------------------
    # Name Column Heading for each dataset:
    # ARS ID      = "Watershed name" (vs. Location, etc.)
    # CAMELS ID   = "gauge_name"
    # FPS ID      = "SiteName"
    # GAGES2_all  = "STANAME"
    # GAGES2_ref  = "STANAME"
    # HLR         =  No name given
    # LTER        = "Site Name" (but not a basin name)
    # MOPEX       = "SiteName"
    # NEON        = "field_site_name"
    # RFC         = "description"
    # USGS_gauged = "station_nm"
    #--------------------------------------------------------
    name_col = None
    list1 = ['FPS', 'GAGES2_all', 'GAGES2_ref', 'LTER'  ]
    list2 = ['CAMELS', 'MOPEX', 'NEON', 'USGS_gauged']
    list3 = ['ARS']

    if (db_str in list1):
        name_col = 1
    elif (db_str in list2):
        name_col = 2
    elif (db_str in list3):
        name_col = 4  # Watershed name
    else:
        #### CZO, HLR and RFC ####
        print('SORRY: No match for:', db_str)
        print('       Returning name_col = 0.')
        name_col = 0
       
    return name_col

#   get_name_column()
#---------------------------------------------------------------------
def get_minlon_column( db_str, SILENT=False ):

    if (db_str == 'CAMELS'):
        col = 3
    elif (db_str == 'GAGES2_all'):
        col = 8
    elif (db_str == 'GAGES2_ref'):
        col = 8      
    elif (db_str == 'HLR'):
        col = 5
    elif (db_str == 'MOPEX'):
        col = 5
    else:
        if not(SILENT):
            print('SORRY: No match for:', db_str)
            print('       Returning col = None.')
        col = None

    return col

#   get_minlon_column()
#---------------------------------------------------------------------
def get_header_lines( db_str='USGS_gauged'):

    if (db_str == 'USGS_gauged'):
        nh_lines = 37
    elif (db_str == 'FPS'):
        nh_lines = 2
    else:
        nh_lines = 1
   
    return nh_lines
    
#   get_header_lines()
#---------------------------------------------------------------------
def get_basin_ids( db_str='USGS_gauged' ):

    file_path = get_data_filepath( db_str )
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

    while (True):
        line = file_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        values = line.split( delim )
        id = values[ id_col ].strip()  
        id_list.append( id )     

    file_unit.close()
    return id_list
      
#   get_basin_ids()
#---------------------------------------------------------------------
def compare_basin_ids( db_str1='USGS_gauged', db_str2='GAGES2_all'):

    id_list1 = get_basin_ids( db_str1 )
    id_list2 = get_basin_ids( db_str2 )
    print('Dataset 1 =', db_str1)
    print('Dataset 2 =', db_str2)
        
    #-----------------------------
    # Compare the two name lists
    #-----------------------------
#     print()   
#     diff1 = list( set(id_list1) - set(id_list2))
#     diff1.sort()
#     print('Basin IDs in dataset 1 not in dataset 2:')
#     print(diff1)
#     print()

    diff2 = list( set(id_list2) - set(id_list1))
    diff2.sort()
    print('Basin IDs in dataset 2 not in dataset 1:')
    print(diff2)
    print('len(diff2) =', len(diff2))
    print()
    ## return diff2

#   compare_basin_ids()

#---------------------------------------------------------------------
def get_basin_names( db_str='USGS_gauged', UPPER=True,
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
    file_path = get_data_filepath( db_str )
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
def compare_basin_names( db_str1='USGS_gauged', db_str2='GAGES2'):

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
def get_dataset_key_names():
  
    names = [
    'ARS', 'CAMELS', 'CZO', 'FPS', 'GAGES2_all', 'HLR',
    'LTER', 'MOPEX', 'NEON', 'USGS_gauged']
    return names

#   get_dataset_key_names()
#---------------------------------------------------------------------
def get_next_line( key, file_units ):

    line = file_units[ key ].readline()
    return line
    
#   get_next_line()
#---------------------------------------------------------------------
def get_next_id( line, key, delim ):

    # print('### key =', key)
    id_col = get_id_column( key )
    vals   = line.split( delim )
    id     = vals[ id_col ].strip()
    # print('### id_col =', id_col)

    #---------------------------------------
    # Prepend '0' if only 7 digits for FPS
    #---------------------------------------    
    if (key == 'FPS') and (len(id) == 7):
        id = '0' + id
    return id
  
#   get_next_id()
#---------------------------------------------------------------------
def get_bounding_box( line, key, delim ):

    try:
        minlon_col = get_minlon_column( key, SILENT=True )
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
def get_state_code( long_name ):

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

    return state_code
   
#   get_state_code()
#---------------------------------------------------------------------
def skip_header_lines( file_unit, key='USGS_gauged'):

    nh_lines = get_header_lines( key )
    for k in range( nh_lines ):
        line = file_unit.readline()
    
#   skip_header_lines()
#---------------------------------------------------------------------
def get_headings_line( key='USGS_gauged'):

    file_path = get_data_filepath( key )
    file_unit = open( file_path, 'r' )
    nh_lines  = get_header_lines( key )
    nh_skip   = nh_lines

    #-----------------------------------------------
    # This one is unusual. It has another "format"
    # line after the main column headings line.
    #-----------------------------------------------    
    if (key == 'USGS_gauged'):
        nh_skip = (nh_lines - 1)
    
    for k in range( nh_skip ):
        line = file_unit.readline()        
    file_unit.close()
    
    headings_line = line  # last one read   
    return headings_line
     
#   get_headings_line()
#---------------------------------------------------------------------
def get_headings( key='USGS_gauged', delim='\t'):

    headings_line = get_headings_line( key=key)
    headings_line = headings_line[:-1]  # remove newline char
    headings      = headings_line.split( delim )
    return headings

#   get_headings()
#---------------------------------------------------------------------
def write_new_header( out_tsv_unit, base_header_line, delim ):
    
#     base_header_line2 = base_header_line[:-1]  # remove new_line
#     base_headings     = base_header_line2.split( delim )
    
    #---------------------------------------
    # Use a subset of the USGS headings
    # Must match values written to new TSV
    #---------------------------------------
    base_headings = [
    'Agency', 'SiteID', 'SiteName', 'SiteType',
    'SiteLat', 'SiteLon', 'Coord Acy. Code', 'CoordDatum',
    ## 'State Code',
    'SiteElev',
    ## 'Elev Acy', 'Elev Datum', 'Basin Code',
    'DrainArea' ]
 
    new_headings  = [
    'LongName', 'State', 'SiteURL',
    'Closest_Station_ID', 'Closest_Station_Dist',
    'StartDate', 'EndDate',
    'MinLon', 'MaxLon', 'MinLat', 'MaxLat',
    #-------------------------------------------------------------
    'Is_USGS_Gauged', 'Is_GAGES2_Ref', 'Is_GAGES2_NonRef',
    'Is_CAMELS', 'Is_FPS', 'Is_MOPEX',
    'Is_ARS', 'Is_CZO', 'IS_HLR', 'IS_LTER', 'Is_NEON' ]   
    new_headings = base_headings + new_headings
    
    header = ''
    for h in new_headings:
        header += (h + delim)
    header = header[:-1]  # remove last delim
    out_tsv_unit.write( header + '\n')
    
#   write_new_header()
#---------------------------------------------------------------------
def collate( out_tsv_file='all_basins.tsv',
             DEBUG=False, max_count=1000):

    #-----------------------------------------------------------
    # Note that the "USGS_gauged" dataset includes many gauges
    # that are no longer active.  The "GAGES2_all" dataset is
    # newer and contains both "ref" and "non-ref" basins.
    #-----------------------------------------------------------
    start_time = time.time() 
    base_key = 'USGS_gauged'   # 25420 data rows
    # base_key = 'GAGES2_all'  # 9322 data rows
 
    repo_dir = get_repo_dir()
    out_tsv_path = repo_dir + out_tsv_file
    out_tsv_unit = open( out_tsv_path, 'w')      #########

    #-----------------------------------------   
    # Get all of the basin dataset key names
    #-----------------------------------------
    key_names = get_dataset_key_names()
    base_header_line = get_headings_line( key=base_key )
    state_code_map = get_state_code_map()

    #------------------------------------------------
    # Open TSV files for each of the basin datasets
    #------------------------------------------------
    delim = '\t'  # tab
    file_units   = dict()
    for key in key_names:
        file_name = get_data_filepath( key )
        file_unit = open(file_name, 'r')
        file_units[ key ] = file_unit
        #--------------------        
        # Skip header lines
        #--------------------
        skip_header_lines( file_unit, key=key)

    #-----------------------------------------
    # Write column headings for new TSV file
    #-----------------------------------------
    write_new_header( out_tsv_unit, base_header_line, delim) 
    
    #---------------------------------------------------
    # Read first station ID from some of the data sets
    #---------------------------------------------------
    key_list1  = ['GAGES2_all', 'CAMELS', 'FPS', 'MOPEX']
    current_id   = dict()
    current_line = dict()
    for key in key_list1:
        line = get_next_line( key, file_units )
        id   = get_next_id( line, key, delim )
        current_id[ key ]   = id
        current_line[ key ] = line

    #-----------------------------------------------------
    # Loop through all USGS gauged basins and use ID
    # name or lon/lat to determine the new fields:
    #    IS_USGS, IS_GAGES2_REF, IS_GAGES2_NONREF,
    #    IS_CAMELS, IS_MOPEX, IS_USGS_FPS.
    # Use the fact that each is sorted by USGS site no.
    #-----------------------------------------------------
#     usgs_id_col   = get_id_column( base_key )
#     usgs_name_col = get_name_column( base_key )
#     usgs_unit     = file_units[ base_key ]
# 
#     IS_USGS_GAUGE_ID = (base_key == 'USGS_gauged')
#     IS_GAGES2        = (base_key == 'GAGES2_all')  # also USGS
#     IS_ARS_ID        = False
#     IS_CZO_ID        = False
#     IS_HLR_ID        = False
#     IS_LTER_ID       = False
#     IS_NEON_ID       = False
# 
#     n_bad_state_codes   = 0
#     n_fixed_state_codes = 0
#     counter = 0
#     bmap = {True:'Y', False:'N'}
#     print('Working on USGS, GAGES2, CAMELS, FPS, & MOPEX basins...')
#     while (True):
#         usgs_line = usgs_unit.readline()
#         if (usgs_line == ''):
#             break  # (reached end of file)
#         values = usgs_line.split( delim )
#         
#         usgs_id    = values[ usgs_id_col ].strip()
#         usgs_name  = values[ usgs_name_col ].strip()
#         long_name  = usgs.get_usgs_long_name( usgs_name )
#         state_code = get_state_code( long_name )
#         if (state_code == ''):
#             n_bad_state_codes += 1
#             try:
#                 state_name, state_code = usgs.get_usgs_missing_state( usgs_id, state_code_map)
#                 n_fixed_state_codes += 1
#             except:
#                 print('WARNING: Could not fix state code for:')
#                 print('         station_id =', usgs_id)
#                 print()
#         site_url   = usgs.get_usgs_site_url( usgs_id )
#         # huc_url  = usgs.get_usgs_huc_url( usgs_id )
#         #-----------------------------------------------
#         closest_station_id   = usgs_id
#         closest_station_dist = '0.0'    # as a string
# 
#         #---------------
#         # For testing
#         #---------------
#         if (DEBUG):
#             counter +=1
#             print('USGS ID       =', usgs_id)
#             print('GAGES2_all ID =', current_id['GAGES2_all'])
#             print('CAMELS ID     =', current_id['CAMELS']) 
#             print('FPS ID        =', current_id['FPS']) 
#             print('MOPEX ID      =', current_id['MOPEX']) 
#        
#         #-----------------------------
#         # Is this also a GAGES2 ID ?
#         #-----------------------------
#         IS_GAGES2_ref     = False
#         IS_GAGES2_non_ref = False
#         key = 'GAGES2_all'
#         gages2_id = current_id[ key ]
#         if (usgs_id == gages2_id):
#             if (DEBUG):
#                 print('  Found match for GAGES2_all.')
#             line = current_line[ key ]
#             vals = line.split( delim )
#             gages2_class = vals[2].lower()  # 3rd column
#             IS_GAGES2_ref     = (gages2_class == 'ref')
#             IS_GAGES2_non_ref = (gages2_class == 'non-ref')
#             bbox = get_bounding_box( line, key, delim )
#         #------------------------------------------------------
#         # Use >=. Dataset may contain IDs not in USGS_gauged.
#         #------------------------------------------------------                               
#         if (usgs_id >= gages2_id):
#             line = get_next_line( key, file_units )
#             if (line != ''):
#                 id = get_next_id( line, key, delim )
#                 current_id[ key ]   = id
#                 current_line[ key ] = line
#  
#         #-----------------------------               
#         # Is this also a CAMELS ID ?
#         #-----------------------------
#         key = 'CAMELS'
#         camels_id = current_id[ key ]
#         IS_CAMELS_ID = (usgs_id == camels_id)
#         if (IS_CAMELS_ID):
#             if (DEBUG):
#                 print('  Found match for CAMELS.')
#             line = current_line[ key ]
#             vals = line.split( delim )
#             bbox = get_bounding_box( line, key, delim )
#         #------------------------------------------------------
#         # Use >=. Dataset may contain IDs not in USGS_gauged.
#         #------------------------------------------------------                             
#         if (usgs_id >= camels_id):
#             line = get_next_line( key, file_units )
#             if (line != ''):
#                 id = get_next_id( line, key, delim )
#                 current_id[ key ]   = id
#                 current_line[ key ] = line
#  
#         #---------------------------               
#         # Is this also an FPS ID ?
#         #---------------------------
#         key = 'FPS'
#         fps_id = current_id[ key ]
#         IS_FPS_ID = (usgs_id == fps_id)
#         if (IS_FPS_ID):
#             if (DEBUG):
#                 print('  Found match for FPS.')
#             line = current_line[ key ]
#             vals = line.split( delim )
#             bbox = get_bounding_box( line, key, delim )
#         #------------------------------------------------------
#         # Use >=. Dataset may contain IDs not in USGS_gauged.
#         #------------------------------------------------------                             
#         if (usgs_id >= fps_id):
#             line = get_next_line( key, file_units )
#             if (line != ''):
#                 id = get_next_id( line, key, delim )
#                 current_id[ key ]   = id
#                 current_line[ key ] = line
# 
#         #----------------------------               
#         # Is this also a MOPEX ID ?
#         #----------------------------
#         key = 'MOPEX'
#         mopex_id = current_id[ key ]
#         IS_MOPEX_ID = (usgs_id == mopex_id)
#         if (IS_MOPEX_ID):
#             if (DEBUG):
#                 print('  Found match for MOPEX.')
#             line = current_line[ key ]
#             vals = line.split( delim )
#             bbox = get_bounding_box( line, key, delim )
#         #------------------------------------------------------
#         # Use >=. Dataset may contain IDs not in USGS_gauged.
#         #------------------------------------------------------                             
#         if (usgs_id >= mopex_id):
#             line = get_next_line( key, file_units )
#             if (line != ''):
#                 id = get_next_id( line, key, delim )
#                 current_id[ key ]   = id
#                 current_line[ key ] = line
# 
#         if (DEBUG):
#             print('==========================================')
#             if (counter == max_count):
#                 break
# 
#         #----------------------------------------
#         # Extract desired fields from usgs_line
#         #----------------------------------------
#         new_usgs_line = get_new_usgs_line( usgs_line, delim, base_key )
#              
#         #-----------------------------------------------------
#         # Write out USGS_gauged line info plus the values
#         # of these boolean vars and site num url & state &
#         # long basin name, & bounding box, if available.
#         #-----------------------------------------------------
#         start_date = 'start_date'  ##### placeholder
#         end_date   = 'end_date'    ##### placeholder
#         
#         val_list = [
#         new_usgs_line,
#         # agency_cd, site_no, station_nm, site_tp_cd,
#         # dec_lat_va, dec_lon_va, coord_acy_cd, coord_datum_cd,
#         # alt_va, drain_area_va,
#         #----------------------------------------------
#         long_name, state_code, site_url,
#         closest_station_id, closest_station_dist,
#         start_date, end_date,
#         bbox[0], bbox[1], bbox[2], bbox[3],
#         # minlon, maxlon, minlat, maxlat
#         #----------------------------------------------
#         bmap[IS_USGS_GAUGE_ID],  bmap[IS_GAGES2_ref],
#         bmap[IS_GAGES2_non_ref], bmap[IS_CAMELS_ID],
#         bmap[IS_FPS_ID], bmap[IS_MOPEX_ID],
#         bmap[IS_ARS_ID], bmap[IS_CZO_ID], bmap[IS_HLR_ID],
#         bmap[IS_LTER_ID], bmap[IS_NEON_ID] ]
#         new_line = ''
#         for val in val_list:
#             new_line += (val + delim)
#         out_tsv_unit.write( new_line + '\n' )

    #---------------------------------------------
    # Get coords of all USGS gauged basins so we
    # can find the closest station and distance
    #---------------------------------------------
    usgs_station_coords = usgs.get_usgs_station_coords()
    
    #-----------------------------------    
    # Add rows for the USDA ARS basins
    #-----------------------------------
    ars_unit = file_units['ARS']
    add_ars_basins( ars_unit, out_tsv_unit,
                    usgs_station_coords, delim)

    #------------------------------    
    # Add rows for the CZO basins
    #------------------------------
    czo_unit = file_units['CZO']
    add_czo_basins( czo_unit, out_tsv_unit,
                    usgs_station_coords, delim)

    #------------------------------
    # Add rows for the HLR basins
    #------------------------------
#     hlr_unit = file_units['HLR']
#     add_hlr_basins( hlr_unit, out_tsv_unit,
#                     usgs_station_coords, delim)

    #-------------------------------    
    # Add rows for the LTER basins
    #-------------------------------
#     lter_unit = file_units['LTER']
#     add_lter_basins( lter_unit, out_tsv_unit,
#                      usgs_station_coords, delim)

    #-------------------------------    
    # Add rows for the NEON basins
    #-------------------------------
    neon_unit = file_units['NEON']
    add_neon_basins( neon_unit, out_tsv_unit,
                     usgs_station_coords, delim)

    #-----------------------------------    
    # Add rows for the NOAA RFC basins
    #-----------------------------------
#     rfc_unit = file_units['RFC']
#     add_rfc_basins( rfc_unit, out_tsv_unit,
#                     usgs_station_coords, delim)
        
    #------------------  
    # Close all files
    #------------------
    out_tsv_unit.close()
    for key in key_names:
        file_units[ key ].close()
    run_time = (time.time() - start_time)
    print('run_time =', run_time, '[secs]')
    ## print('n_bad_state_codes   =', n_bad_state_codes)
    ## print('n_fixed_state_codes =', n_fixed_state_codes)
    print('Finished.')
    print()
                     
#   collate()
#---------------------------------------------------------------------
def get_new_usgs_line( usgs_line, delim, base_key='USGS_gauged' ):

    #-----------------------------------------------
    # 0:  agency_code, e.g. 'USGS'
    # 1:  station/site id (8 digits)
    # 2:  station name w/ abbreviations
    # 3:  site_type_code; canal, ditch, etc.
    # 4:  decimal lat.
    # 5:  decimal lon. (dec_long_va)
    # 6:  coord accuracy code
    # 7:  coord datum code
    # 8:  state_code (2-digit numeric)
    # 9:  altitude of gauge or land
    # 10: altitude accuracy value (often missing)
    # 11: altitude datum code (always missing)
    # 12: drainage basin code (often missing)
    # 13: drainage area value
    # 14: contrib. drain. area (usu. missing)
    # 15: reliability code (usually missing)
    #-------------------------------------------------
    # Note: 'NAD27' shows as "NAD 27.00" in Numbers
    #       but string value is just 'NAD27'.
    #-------------------------------------------------
    usgs_line = usgs_line[:-1]  # remove newline
    usgs_fields = usgs_line.split( delim )
    n_fields = len(usgs_fields)
    new_usgs_line = ''

    if (base_key == 'USGS_gauged'):
        skip_cols = [8, 10, 11, 12, 14, 15]
        for col in range( n_fields ):
            if (col not in skip_cols):
                field = usgs_fields[ col ]
                # print('field =', field)                  
                new_usgs_line += (field + delim)
        new_usgs_line = new_usgs_line[:-1]  # remove last delim 
    else:
        print('ERROR in get_new_usgs_line:')
        print('  base_key not supported:', base_key)
           
    return new_usgs_line
    
#   get_new_usgs_line()
#---------------------------------------------------------------------
def add_ars_basins( ars_unit, out_tsv_unit,
                    usgs_station_coords, delim ):

    ars_id_col   = get_id_column( 'ARS' )
    ars_name_col = get_name_column( 'ARS' )
    agency_code  = 'USDA/ARS'

    print('Working on USDA/ARS basins...')
    while (True):
        ars_line = ars_unit.readline()
        if (ars_line == ''):
            break  # (reached end of file)
        values = ars_line.split( delim )
        #-------------------------------------------------
        ars_id     = values[ ars_id_col ].strip()
        ## ars_name   = values[ ars_name_col ].strip()
        ars_name   = get_ars_name( values )
        state_code = values[2].strip()   # 2-letter state code
        long_name  = get_ars_long_name( values )
        lat_str    = values[7].strip()   # outlet latitude
        lon_str    = values[8].strip()   # outlet longitude
        area       = values[10].strip()  # basin area
        start_date = values[11].strip()  ###########################
        end_date   = values[12].strip()
        site_url   = get_ars_site_url( state_code )
        notes      = values[15].strip()
        #--------------------------------------------
        if (lat_str == ''):
            lat_str = '-999'
        if (lon_str == ''):
            lon_str = '-999'
        #--------------------------------------------           
        closest_id, clon, clat, dmin = \
            get_closest_station(usgs_station_coords, 
                lon_str, lat_str, REPORT=False)
        closest_station_id   = closest_id
        closest_station_dist = str(dmin)    # as a string
        #--------------------------------------------
        # Don't know these
        #--------------------
        bbox = ['-999', '-999', '-999', '-999']
        site_type      = '--'
        coord_acy_cd   = '--'
        coord_datum_cd = 'NAD93??'   ####### FIND OUT
        elev           = '-9999'
    
        #------------------------------------------
        # Write out values for USDA ARS watershed
        #------------------------------------------
        val_list = [
        agency_code, ars_id, ars_name, site_type, lat_str, lon_str,
        coord_acy_cd, coord_datum_cd, elev, area,
        #----------------------------------------------
        long_name, state_code, site_url,
        closest_station_id, closest_station_dist,
        start_date, end_date,
        bbox[0], bbox[1], bbox[2], bbox[3],
            # minlon, maxlon, minlat, maxlat
        #------------------------------------------------------
        # IS_USGS_GAUGE_ID, IS_GAGES2_ref, IS_GAGES2_non_ref,
        # IS_CAMELS_ID, IS_FPS_ID, IS_MOPEX_ID, IS_ARS_ID, ###########
        # IS_CZO_ID, IS_HLR_ID, IS_LTER_ID, IS_NEON_ID
        #------------------------------------------------------
        'N','N','N','N','N','N','Y','N','N','N','N' ]

        new_line = ''
        for val in val_list:
            new_line += (val + delim)
        out_tsv_unit.write( new_line + '\n' )
        
#   add_ars_basins()
#---------------------------------------------------------------------
def get_ars_name( values ):

    #-----------------------------------------------------
    # Note:  Alternate name is not given in source data,
    #        but was found by other means in some cases.
    #-----------------------------------------------------
    location   = values[1].strip()  # often city name
    state_code = values[2].strip()  # 2-letter state code
    # alt_name   = values[3].strip()  # alternate river name
    wshed_name = values[4].strip()  # watershed name
    
    ars_name = wshed_name + ', ' + location + ', ' + state_code

    return ars_name
       
#   get_ars_name()
#---------------------------------------------------------------------
def get_ars_long_name( values ):

    location   = values[1].strip()  # often city name
    state_code = values[2].strip()  # 2-letter state code
    alt_name   = values[3].strip()  # alternate river name
    wshed_name = values[4].strip()  # watershed name
    
    if (alt_name != 'Alternate_Name'):
        long_name = alt_name + ', '
    else:
        long_name = ''
    long_name += wshed_name + ', '
    long_name += location + ', ' + state_code

    return long_name
       
#   get_ars_long_name()
#---------------------------------------------------------------------
def get_ars_site_url( state_code ):

    url  = 'https://hrsl.ba.ars.usda.gov/wdc/'
    url += state_code.lower() + '.htm'
    return url

#   get_ars_site_url()
#---------------------------------------------------------------------
def add_czo_basins( czo_unit, out_tsv_unit,
                    usgs_station_coords, delim ):

    #------------------------------------------------
    # Note: CZO Lidar DEMs are available from:
    # https://portal.opentopography.org/dataSearch?
    #    search=Critical%20Zone%20Observatory
    #------------------------------------------------
    agency_code  = 'NSF-CZO'

    print('Working on CZO basins...')
    while (True):
        # Don't add .strip() to czo_line;
        # will remove tabs for blank entries
        czo_line = czo_unit.readline()
        if (czo_line == ''):
            break  # (reached end of file)
        values = czo_line.split( delim )
        #-------------------------------------------------
        czo_oname  = values[0].strip()
        wshed_name = values[1].strip()
        state_code = values[2].strip()
        location   = values[3].strip()
        lat_str    = values[4].strip()
        lon_str    = values[5].strip()
        elev       = values[6].strip()
        start_date = values[7].strip()
        end_date   = values[8].strip()
        # timestep   = values[9].strip()
        area       = values[10].strip()
        # instrument = values[11].strip()
        credit     = values[12].strip()
        #-------------------------------------------------
        site_id    = get_czo_site_id( values )
        # czo_name   = get_czo_site_name( values )
        long_name  = get_czo_long_name( values )
        site_url   = get_czo_site_url( site_id )
        #-------------------------------------------------          
        closest_id, clon, clat, dmin = \
            get_closest_station(usgs_station_coords, 
                lon_str, lat_str, REPORT=False)
        closest_station_id   = closest_id
        closest_station_dist = str(dmin)    # as a string
        #--------------------------------------------
        # Don't know these
        #--------------------
        bbox = ['-999', '-999', '-999', '-999']
        site_type      = '--'
        coord_acy_cd   = '--'
        coord_datum_cd = 'NAD93??'   ####### FIND OUT

        #-------------------------------------------------------           
        # It looks like 2 of the CZO basins are, or are close
        # to, USGS gauged basins. And 1 is a USDA basin
        #-------------------------------------------------------
        IS_USGS_GAUGE_ID = False
        IS_ARS_ID = False
        bmap = {True:'Y', False:'N'}

        if (site_id == 'CZO-LUQ-IC'):
            usgs_id   = '50075000'
            usgs_name = 'Rio Icacos near Naguabo, PR'
            IS_USGS_GAUGE_ID = True
            # IS_USGS_GAUGE_ID = (closest_station_id == usgs_id)
            if (closest_station_id != usgs_id):
                print('WARNING: Possible issue with CZO-LUQ-IC basin')
                print('   and the closest USGS station ID.') 
        if (site_id == 'CZO-LUQ-MAM'):
            usgs_id   = '50065500'
            usgs_name = 'Rio Mameyes near Sabana, PR'
            IS_USGS_GAUGE_ID = True
            # IS_USGS_GAUGE_ID = (closest_station_id == usgs_id)
            if (closest_station_id != usgs_id):
                print('WARNING: Possible issue with CZO-LUQ-MAM basin')
                print('   and the closest USGS station ID.')                 
        if (site_id == 'CZO-RC-JD'):
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
            ars_id   = ''
            ars_name = 'Johnston Draw in Reynolds Creek, ID'
            IS_ARS_ID = True

        #------------------------------------------
        # Write out values for CZO watershed
        #------------------------------------------
        val_list = [
        agency_code, site_id, wshed_name, site_type, lat_str, lon_str,
        coord_acy_cd, coord_datum_cd, elev, area,
        #----------------------------------------------
        long_name, state_code, site_url,
        closest_station_id, closest_station_dist,
        start_date, end_date,
        bbox[0], bbox[1], bbox[2], bbox[3],
            # minlon, maxlon, minlat, maxlat              
        #------------------------------------------------------
        # IS_USGS_GAUGE_ID, IS_GAGES2_ref, IS_GAGES2_non_ref,
        # IS_CAMELS_ID, IS_FPS_ID, IS_MOPEX_ID, IS_ARS_ID, ###########
        # IS_CZO_ID, IS_HLR_ID, IS_LTER_ID, IS_NEON_ID
        #------------------------------------------------------
        bmap[IS_USGS_GAUGE_ID],'N','N','N','N','N',
        bmap[IS_ARS_ID], 'Y','N','N','N' ]

        new_line = ''
        for val in val_list:
            new_line += (val + delim)
        out_tsv_unit.write( new_line + '\n' )
        
#   add_czo_basins()
#---------------------------------------------------------------------
def get_czo_site_id( values ):

    #-------------------------------------------------
    # Note:  CZO stations don't have an official ID,
    #        so we make one from other info.
    #-------------------------------------------------
    obs_name = values[0].strip()
    p1 = obs_name.find('(')
    obs_abbr = obs_name[p1:]
    obs_abbr = obs_abbr[1:-1]
    #--------------------------------
    catch_name = values[1].strip()
    p2 = catch_name.find('(')
    catch_abbr = catch_name[p2:]
    catch_abbr = catch_abbr[1:-1]  
    #----------------------------------------------
    czo_id = 'CZO-' + obs_abbr + '-' + catch_abbr
    return czo_id

#   get_czo_site_id()
#---------------------------------------------------------------------
def get_czo_long_name( values ):

    catch_name = values[1].strip()
    state_code = values[2].strip()
        
    if ('(P30' in catch_name):
        long_name = catch_name
    else:
        p = catch_name.find('(')
        long_name = catch_name[:p].strip()
    long_name += ', ' + state_code

    return long_name
    
#   get_czo_long_name() 
#---------------------------------------------------------------------
def get_czo_site_url( site_id, REPORT=False,
                      ARCHIVE_SITES=True,ARCHIVE_BASE=False,
                      HYDROSHARE=False ):

    #----------------------------------------------------------
    # Note: Santa Catalina Mountains (SCM) and Jemez River
    # Basin (JRB) are now treated as one CZO for some reason,
    # called "Catalina-Jemez CZO".
    #----------------------------------------------------------
    # See CZO Catalina-Jemez User Profile in HydroShare:
    # https://www.hydroshare.org/user/5407/
    #----------------------------------------------------------
    czo_archive_base_url = 'https://czo-archive.criticalzone.org/'
    hydroshare_base_url  = 'https://www.hydroshare.org/group/'
    # earthchem_base_url   = 'https://www.earthchem.org/' # has chem data?
 
    #-----------------------------------------------------------------    
    # Notes:  (1) in CH CZO, there are 8 research areas; which one?
    #         (2) In Eel CZO, could not find URL for Elder Creek.
    #         (3) In JRB CZO, could not find URL for Upper Jaramillo.
    #         (4) In JRB CZO, could not find URL for La Jara.
    #         (5) In SH CZO,  could not determine watershed.
    #------------------------------------------------------------------           
    site_id_map = {
    'CZO-BC-GG':     ['boulder',         'gordon-gulch', 140],
    'CZO-CH-WS4':    ['calhoun',         'calhoun-czo-research-area-1', 148], # there are 8?
    'CZO-CRB-CR':    ['christina',       'christina-river-basin', 146],
    'CZO-ER-EC':     ['eel',             'eel-river-watershed', 141],  # not for elder creek
    'CZO-IML-A-US':  ['iml',             'sangamon-river-basin', 150],
    'CZO-IML-A2-US': ['iml',             'sangamon-river-basin', 150],
    'CZO-IML-B-CC':  ['iml',             'clear-creek-watershed', 150], 
    'CZO-IML-C-MR':  ['iml',             'minnesota-river-basin', 150],
    'CZO-JRB-UJ':    ['catalina-jemez',  'santa-catalina-mountains', 142], ## better site?
    'CZO-JRB-LJ':    ['catalina-jemez',  'santa-catalina-mountains', 142], ## better site?
    'CZO-LUQ-MAM':   ['luquillo',        'rio-icacos',  144],
    'CZO-LUQ-IC':    ['luquillo',        'rio-mameyes', 144],
    'CZO-RC-JD':     ['reynolds',        'johnston-draw', 143],
    'CZO-SCM-OR':    ['catalina-jemez',  'oracle-ridge-mid-elevation', 142],
    'CZO-SCM-MG':    ['catalina-jemez',  'bigelow-tower-marshall-gulch-high-elevation', 142],
    'CZO-SS-P301':   ['sierra',          'providence-creek-subcatchment-p301', 145],
    'CZO-SS-P303':   ['sierra',          'providence-creek-subcatchment-p303', 145],
    'CZO-SS-P304':   ['sierra',          'providence-creek-subcatchment-p304', 145],
    'CZO-SH-SH':     ['shale-hills',     'cole-farm-agricultural-site', 147] }  # what river?
    
    # Note:  hs_group_num = HydroShare Group number
    czo_abbr, site_abbr, hs_group_num = site_id_map[ site_id ]
    czo_field_areas_str = 'infrastructure/field-areas-'
    czo_field_area_str  = 'infrastructure/field-area/'  # specific area
    #----------------------------------------------------------------------
    czo_archive_url1  = czo_archive_base_url + czo_abbr + '/'
    czo_archive_url2  = czo_archive_url1
    czo_archive_url2 += czo_field_area_str   + site_abbr + '/'
    ## czo_archive_url2 += czo_field_areas_str  + czo_abbr + '/'
    #----------------------------------------------------------------------
    hydroshare_url = hydroshare_base_url + str(hs_group_num)
    
    if (REPORT):
        print('czo_url1 =', czo_archive_url1)
        print('czo_url2 =', czo_archive_url12)
        print('hydroshare_url =', hydroshare_url)
        print()

    if (ARCHIVE_SITES):
        return czo_archive_url2   # includes url1
    elif (ARCHIVE_BASE):
        return czo_archive_url1
    if (HYDROSHARE):
        return hydroshare_url

    #--------------------------------------------------
    # Note that CZEN also has some CZO websites like:
    #--------------------------------------------------
    # https://www.czen.org/content/calhoun-czo-1
        
    #------------------------------------------------------------
    # Note that individual field areas have URLs that are like
    # the one below.  Could scrape these to check data values.
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/catalina-jemez/infrastructure/
    #         field-area/oracle-ridge-mid-elevation/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/national/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/boulder/infrastructure/field-area/betasso/
    # https://czo-archive.criticalzone.org/boulder/infrastructure/field-area/gordon-gulch/
    # https://czo-archive.criticalzone.org/boulder/infrastructure/field-area/green-lakes-valley/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-1/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-2/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-3/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-4/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-5/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-6/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-7/
    # https://czo-archive.criticalzone.org/calhoun/infrastructure/field-area/calhoun-czo-research-area-8/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/catalina-jemez/infrastructure/field-area/santa-catalina-mountains/
    # https://czo-archive.criticalzone.org/catalina-jemez/infrastructure/field-area/oracle-ridge-mid-elevation/
    # https://czo-archive.criticalzone.org/catalina-jemez/infrastructure/field-area/bigelow-tower-marshall-gulch-high-elevation/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/christina/about/
    # https://czo-archive.criticalzone.org/christina/infrastructure/field-area/christina-river-basin/
    #------------------------------------------------------------
    # THERE SEEMS TO BE NO PAGE FOR ELDER CREEK.
    # https://czo-archive.criticalzone.org/eel/infrastructure/field-area/eel-river-watershed/
    # https://czo-archive.criticalzone.org/eel/infrastructure/field-area/rivendell/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/iml/infrastructure/field-area/sangamon-river-basin/
    # https://czo-archive.criticalzone.org/iml/infrastructure/field-area/clear-creek-watershed/
    # https://czo-archive.criticalzone.org/iml/infrastructure/field-area/minnesota-river-basin/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/luquillo/infrastructure/field-area/bisley/
    # https://czo-archive.criticalzone.org/luquillo/infrastructure/field-area/puente-roto/
    # Puente Roto USGS station:  USGS 50065500 RIO MAMEYES NR SABANA, PR
    # https://czo-archive.criticalzone.org/luquillo/infrastructure/field-area/rio-blanco-nr-florida/
    # https://czo-archive.criticalzone.org/luquillo/infrastructure/field-area/rio-icacos/
    # https://czo-archive.criticalzone.org/luquillo/infrastructure/field-area/rio-mameyes/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/reynolds-creek-experimental-watershed/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/murphy-creek/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/tollgate/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/upper-sheep-creek/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/johnston-draw/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/mountain-big-sage/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/lower-sage/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/wyoming-big-sage/
    # https://czo-archive.criticalzone.org/reynolds/infrastructure/field-area/flats/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/shale-hills/infrastructure/field-area/cole-farm-agricultural-site/
    #------------------------------------------------------------
    # https://czo-archive.criticalzone.org/sierra/infrastructure/field-area/providence-creek-subcatchment-p301/        
    # https://czo-archive.criticalzone.org/sierra/infrastructure/field-area/providence-creek-subcatchment-p303/
    # https://czo-archive.criticalzone.org/sierra/infrastructure/field-area/providence-creek-subcatchment-p304/

#   get_czo_site_url()
#---------------------------------------------------------------------
def add_hlr_basins( hlr_unit, out_tsv_unit,
                    usgs_station_coords, delim ):

    hlr_id_col   = get_id_column( 'HLR' )
    hlr_name_col = get_name_column( 'HLR' )
    agency_code  = 'USGS-HLR'

    print('Working on USGS HLR basins...')
    while (True):
        hlr_line = hlr_unit.readline()
        if (hlr_line == ''):
            break  # (reached end of file)
        values = hlr_line.split( delim )
        #-------------------------------------------------
#         ars_id     = values[ ars_id_col ].strip()
#         ## ars_name   = values[ ars_name_col ].strip()
#         ars_name   = get_ars_name( values )
#         state_code = values[2].strip()   # 2-letter state code
#         long_name  = get_ars_long_name( values )
#         lat_str    = values[7].strip()   # outlet latitude
#         lon_str    = values[8].strip()   # outlet longitude
#         area       = values[10].strip()  # basin area
#         start_date = values[11].strip()  ###########################
#         end_date   = values[12].strip()
#         site_url   = get_ars_site_url( state_code )
#         notes      = values[15].strip()
        #--------------------------------------------
        if (lat_str == ''):
            lat_str = '-999'
        if (lon_str == ''):
            lon_str = '-999'
        #--------------------------------------------           
        closest_id, clon, clat, dmin = \
            get_closest_station(usgs_station_coords, 
                lon_str, lat_str, REPORT=False)
        closest_station_id   = closest_id
        closest_station_dist = str(dmin)    # as a string
        #--------------------------------------------
        # Don't know these
        #--------------------
        # bbox = ['-999', '-999', '-999', '-999']
        site_type      = '--'
        coord_acy_cd   = '--'
        coord_datum_cd = 'NAD93??'   ####### FIND OUT
        elev           = '-9999'
    
        #------------------------------------------
        # Write out values for USDA ARS watershed
        #------------------------------------------
        val_list = [
        agency_code, ars_id, ars_name, site_type, lat_str, lon_str,
        coord_acy_cd, coord_datum_cd, elev, area,
        #----------------------------------------------
        long_name, state_code, site_url,
        closest_station_id, closest_station_dist,
        start_date, end_date,
        bbox[0], bbox[1], bbox[2], bbox[3],
            # minlon, maxlon, minlat, maxlat
        #------------------------------------------------------
        # IS_USGS_GAUGE_ID, IS_GAGES2_ref, IS_GAGES2_non_ref,
        # IS_CAMELS_ID, IS_FPS_ID, IS_MOPEX_ID, IS_ARS_ID, ###########
        # IS_CZO_ID, IS_HLR_ID, IS_LTER_ID, IS_NEON_ID
        #------------------------------------------------------
        'N','N','N','N','N','N','N','N','Y','N','N' ]

        new_line = ''
        for val in val_list:
            new_line += (val + delim)
        out_tsv_unit.write( new_line + '\n' )
        
#   add_hlr_basins()
#---------------------------------------------------------------------
def add_lter_basins( lter_unit, out_tsv_unit,
                     usgs_station_coords, delim ):

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
    agency_code  = 'NSF-LTER'

    print('Working on LTER basins...')
    while (True):
        lter_line = lter_unit.readline()
        if (lter_line == ''):
            break  # (reached end of file)
        values = lter_line.split( delim )
        #-------------------------------------------------
#         czo_oname  = values[0].strip()
#         wshed_name = values[1].strip()
#         state_code = values[2].strip()
#         location   = values[3].strip()
#         lat_str    = values[4].strip()
#         lon_str    = values[5].strip()
#         elev       = values[6].strip()
#         start_date = values[7].strip()
#         end_date   = values[8].strip()
#         # timestep   = values[9].strip()
#         area       = values[10].strip()
#         # instrument = values[11].strip()
#         credit     = values[12].strip()
#         #-------------------------------------------------
#         site_id    = get_lter_site_id( values )
#         # lter_name   = get_lter_site_name( values )
#         long_name  = get_lter_long_name( values )
#         site_url   = get_lter_site_url( site_id )
        #-------------------------------------------------          
        closest_id, clon, clat, dmin = \
            get_closest_station(usgs_station_coords, 
                lon_str, lat_str, REPORT=False)
        closest_station_id   = closest_id
        closest_station_dist = str(dmin)    # as a string
        #--------------------------------------------
        # Don't know these
        #--------------------
        bbox = ['-999', '-999', '-999', '-999']
        site_type      = '--'
        coord_acy_cd   = '--'
        coord_datum_cd = 'NAD93??'   ####### FIND OUT

        #-------------------------------------------------------           
        # It looks like 2 of the CZO basins are, or are close
        # to, USGS gauged basins. And 1 is a USDA basin
        #-------------------------------------------------------
        IS_USGS_GAUGE_ID = False
        IS_ARS_ID = False
        bmap = {True:'Y', False:'N'}

        #--------------------------------------
        # Write out values for LTER watershed
        #--------------------------------------
        val_list = [
        agency_code, site_id, wshed_name, site_type, lat_str, lon_str,
        coord_acy_cd, coord_datum_cd, elev, area,
        #----------------------------------------------
        long_name, state_code, site_url,
        closest_station_id, closest_station_dist,
        start_date, end_date,
        bbox[0], bbox[1], bbox[2], bbox[3],
            # minlon, maxlon, minlat, maxlat              
        #------------------------------------------------------
        # IS_USGS_GAUGE_ID, IS_GAGES2_ref, IS_GAGES2_non_ref,
        # IS_CAMELS_ID, IS_FPS_ID, IS_MOPEX_ID, IS_ARS_ID, ###########
        # IS_CZO_ID, IS_HLR_ID, IS_LTER_ID, IS_NEON_ID
        #------------------------------------------------------
        bmap[IS_USGS_GAUGE_ID],'N','N','N','N','N',
        bmap[IS_ARS_ID], 'Y','N','N','N' ]

        new_line = ''
        for val in val_list:
            new_line += (val + delim)
        out_tsv_unit.write( new_line + '\n' )
        
#   add_lter_basins()
#---------------------------------------------------------------------
def add_neon_basins( neon_unit, out_tsv_unit,
                     usgs_station_coords, delim ):

    #-------------------------------------------------------
    # Notes:  Metadata is from:
    #-------------------------------------------------------         
    agency_code  = 'NSF-NEON'

    print('Working on NEON basins...')
    while (True):
        neon_line = neon_unit.readline()
        if (neon_line == ''):
            break  # (reached end of file)
        values = neon_line.split( delim )
        #-------------------------------------------------
        field_domain_id        = values[0].strip()
        field_site_id          = values[1].strip()
        field_site_name        = values[2].strip()
#         field_site_type        = values[3].strip()
#         field_site_subtype     = values[4].strip()
#         field_collocated_site  = values[5].strip()
#         field_site_host        = values[6].strip()
        field_site_url         = values[7].strip()
        field_nonneon_res_alwd = values[8].strip()
        field_access_details   = values[9].strip()
        field_neon_op_office   = values[10].strip()
        field_latitude         = values[11].strip()
        field_longitude        = values[12].strip()
        field_geodetic_datum   = values[13].strip()
        field_utm_north        = values[14].strip()
        field_utm_easting      = values[15].strip()
        field_utm_zone         = values[16].strip()
#         field_site_county      = values[17].strip()
        field_site_state       = values[18].strip()  # 2-letter code
#         field_site_country     = values[19].strip()
#         field_site_mean_elev_m = values[20].strip()
        field_site_min_elev_m  = values[21].strip()  ####
#         field_site_max_elev_m  = values[22].strip()
#         field_mean_annu_temp_c = values[23].strip()
        #-----------------------------------------------
#         field_mean_ann_precip_mm = values[24].strip()
#         field_dom_wind_direction = values[25].strip()
#         field_mean_canopy_height = values[26].strip()
#         field_dom_nlcd_classes   = values[27].strip()
#         field_dom_plant_species  = values[28].strip()
        # e.g. [h10250001](https://water.usgs.gov/lookup/getwatershed?10250001)
        field_usgs_huc           = values[29].strip()                       
        field_watershed_name     = values[30].strip()
        field_watershed_size_km2 = values[31].strip()
#         field_lake_depth_mean_m  = values[32].strip()
#         field_lake_depth_max_m   = values[33].strip()
#         field_tower_height_m     = values[34].strip()
#         field_usgs_geology_unit  = values[35].strip()
#         field_megapit_soil_famly = values[36].strip()
#         field_soil_subgroup      = values[37].strip()
#         field_avg_num_green_days = values[38].strip()
#         field_avg_green_incr_doy = values[39].strip()
#         field_avg_green_max_doy  = values[40].strip()
#         field_avg_green_decr_doy = values[41].strip()
#         field_avg_green_min_doy  = values[42].strip()
#         field_phenocams          = values[43].strip()
#         field_num_tower_levels   = values[44].strip()       
        #-----------------------------------------------
        p1 = field_usgs_huc.find('[')
        p2 = field_usgs_huc.find(']')
        huc_str = field_usgs_huc[p1+1: p2]
        huc_url = 'https://water.usgs.gov/lookup/getwatershed?' + huc_str 
        #-----------------------------------------------
        # print('field_site_id =', field_site_id )              
        site_id = 'NEON-' + field_domain_id + '-' + field_site_id
        wshed_name     = field_watershed_name
        # print('Watershed name =', wshed_name)
        lat_str        = field_latitude
        lon_str        = field_longitude
        coord_datum_cd = field_geodetic_datum
        elev           = field_site_min_elev_m
        area           = field_watershed_size_km2
        state_code     = field_site_state
        long_name      = field_site_name + ', ' + wshed_name + ', ' + state_code
        site_url       = field_site_url
        #-------------------------------------------------          
        closest_id, clon, clat, dmin = \
            get_closest_station(usgs_station_coords, 
                lon_str, lat_str, REPORT=False)
        closest_station_id   = closest_id
        closest_station_dist = str(dmin)    # as a string
        #--------------------------------------------
        # Don't know these
        #--------------------
        start_date = '--'
        end_date   = '--'
        bbox = ['-999', '-999', '-999', '-999']
        site_type      = '--'
        coord_acy_cd   = '--'

        #############################
        # Need to check these 
        #############################
        IS_USGS_GAUGE_ID = False
        IS_ARS_ID = False
        bmap = {True:'Y', False:'N'}

        #--------------------------------------
        # Write out values for NEON watershed
        # Skip any other NEON sites
        #--------------------------------------
        if (field_watershed_name != ''):
            val_list = [
            agency_code, site_id, wshed_name, site_type, lat_str, lon_str, 
            coord_acy_cd, coord_datum_cd, elev, area,
            #----------------------------------------------
            long_name, state_code, site_url,
            closest_station_id, closest_station_dist,
            start_date, end_date,
            bbox[0], bbox[1], bbox[2], bbox[3],
                # minlon, maxlon, minlat, maxlat              
            #------------------------------------------------------
            # IS_USGS_GAUGE_ID, IS_GAGES2_ref, IS_GAGES2_non_ref,
            # IS_CAMELS_ID, IS_FPS_ID, IS_MOPEX_ID, IS_ARS_ID, ###########
            # IS_CZO_ID, IS_HLR_ID, IS_LTER_ID, IS_NEON_ID
            #------------------------------------------------------
            bmap[IS_USGS_GAUGE_ID],'N','N','N','N','N',
            bmap[IS_ARS_ID], 'N','N','N','Y' ]

            new_line = ''
            for val in val_list:
                new_line += (val + delim)
            out_tsv_unit.write( new_line + '\n' )
        
#   add_neon_basins()
#---------------------------------------------------------------------      
#---------------------------------------------------------------------
def export_csv_to_tsv( csv_file='new_gages2_all.csv',
                       tsv_file='new_gages2_all.tsv'):

    gages2_dir = get_gages2_data_dir( gtype='root' )
    csv_path   = gages2_dir + csv_file
    tsv_path   = gages2_dir + tsv_file
    
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





