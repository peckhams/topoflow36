
# Copyright (c) 2023, Scott D. Peckham
#
# Jun 2023. Utils to collate several river basin datasets.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils import collate_basins as cb
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
#  get_data_directory()
#  get_data_filepath()
#  get_id_column()
#  get_name_column()
#  get_header_lines()
#  get_basin_ids()
#  compare_basin_ids()
#
#  replace_periods()
#  replace_punctuation()
#  replace_st()
#  replace_state_letters()
#  fix_spellings()
#  replace_city_letters()
#  expand_basin_name()
#  get_basin_names()
#  compare_basin_names()
#
#  collate()   ##### Not finished yet.
#
#---------------------------------------------------------------------

import numpy as np
import re  # for regular expressions

# from osgeo import ogr, osr
# import json, sys, time

#---------------------------------------------------------------------
def get_data_directory( db_str='USGS_gauged' ):

    repo_dir = '/Users/peckhams/Dropbox/NOAA_NextGen/'
    repo_dir += '__NextGen_Example_Basin_Repo/'

    #---------------------------------------------    
    # Directories for each of the basin datasets
    #---------------------------------------------
    if (db_str == 'ARS'):
        data_dir = repo_dir + 'USDA_Watersheds/ARS_Water_Database/'
    elif (db_str == 'CAMELS'):
        data_dir = repo_dir + 'CAMELS/'
    elif (db_str == 'CZO'):
        data_dir = repo_dir + 'CZO/'
    elif (db_str == 'GAGES2_all'):
        data_dir = repo_dir + 'GAGES-II/'
    elif (db_str == 'GAGES2_ref'):
        data_dir = repo_dir + 'GAGES-II/'
    elif (db_str == 'HLR'):
        data_dir = repo_dir + 'HLR_Data_Sets/hlrshape/'
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
        file_path = data_dir + 'ars_db.csv'
    elif (db_str == 'CAMELS'):
        file_path = data_dir + 'new_camels.tsv'
    elif (db_str == 'CZO'):
        file_path = data_dir + 'czo_metadata.csv'  #### NOT READY YET
    elif (db_str == 'GAGES2_all'):
        file_path = data_dir + 'new_gages2_all.csv'
    elif (db_str == 'GAGES2_ref'):
        file_path = data_dir + 'new_gages2_ref.csv'
    elif (db_str == 'HLR'):
        file_path = data_dir + 'new_hlr_na_all_v2.csv'
    elif (db_str == 'LTER'):
        file_path = data_dir + 'new_lter_info.csv'
    elif (db_str == 'MOPEX'):
        file_path = data_dir + 'new_mopex431.tsv'
    elif (db_str == 'NEON'):
        file_path = data_dir + 'NEON_Field_Site_Metadata_20230309.csv'
#     elif (db_str == 'RFC'):
#         file_path = data_dir + 'new_rfc_info.csv'
    elif (db_str == 'USGS_gauged'):
        file_path = data_dir + 'USGS_stream_gauge_data.tsv'
    else:
        print('SORRY: No match for:', db_str)
        file_path = None
       
    return file_path
 
#   get_data_filepath()
#---------------------------------------------------------------------
def get_id_column( db_str ):

    if (db_str == 'ARS'):
        id_col = 0  # heading: "Watershed ID"
    elif (db_str == 'CAMELS'):
        id_col = 0    # heading: "gauge_id"
#     elif (db_str == 'CZO'):
#         id_col = 
    elif (db_str == 'GAGES2_all'):
        id_col = 0    # heading: "STAID" 
    elif (db_str == 'GAGES2_ref'):
        id_col = 0    # heading: "STAID"
    elif (db_str == 'HLR'):
        id_col = -1
        print('SORRY: No names available for HLR basins.')
        print('       Returning -1.')
    elif (db_str == 'LTER'):
        id_col = 0   # heading: "Site Acronym" (e.g. "AND")
    elif (db_str == 'MOPEX'):
        id_col = 1   # heading: "SiteCode"
    elif (db_str == 'NEON'):
        id_col = 1   # heading: "field_site_id" (e.g. "ABBY")
#     elif (db_str == 'RFC'):
#         id_col = 
    elif (db_str == 'USGS_gauged'):
        id_col = 1   # heading: "site_no"
    else:
        print('SORRY: No match for:', db_str)
        id_col = None
       
    return id_col

#   get_id_column()
#---------------------------------------------------------------------
def get_name_column( db_str ):

    if (db_str == 'ARS'):
        name_col = 3  # Use "Alternate Name" (vs. Location, etc.)
    elif (db_str == 'CAMELS'):
        name_col = 2    # heading: "gauge_name"
#     elif (db_str == 'CZO'):
#         name_col = 
    elif (db_str == 'GAGES2_all'):
        name_col = 1    # heading: "STANAME" 
    elif (db_str == 'GAGES2_ref'):
        name_col = 1    # heading: "STANAME"
    elif (db_str == 'HLR'):
        name_col = -1
        print('SORRY: No names available for HLR basins.')
        print('       Returning -1.')
    elif (db_str == 'LTER'):
        name_col = 1   # heading: "Site Name" (but not a basin name)
    elif (db_str == 'MOPEX'):
        name_col = 2   # heading: "SiteName"
    elif (db_str == 'NEON'):
        name_col = 2   # heading: "field_site_name"
#     elif (db_str == 'RFC'):
#         name_col = 
    elif (db_str == 'USGS_gauged'):
        name_col = 2   # heading: "station_nm"  #### BIG HEADER
    else:
        print('SORRY: No match for:', db_str)
        name_col = None
       
    return name_col

#   get_name_column()
#---------------------------------------------------------------------
def get_header_lines( db_str='USGS_gauged'):

    if (db_str == 'USGS_gauged'):
        nh_lines = 37
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
def compare_basin_ids( db_str1='USGS_gauged', db_str2='GAGES2'):

    id_list1 = get_basin_ids( db_str1 )
    id_list2 = get_basin_ids( db_str2 )
    print('Dataset 1 =', db_str1)
    print('Dataset 2 =', db_str2)
        
    #-----------------------------
    # Compare the two name lists
    #-----------------------------
#     print()   
#     diff1 = list( set(id_list1) - set(id_list2))
#     print('Basin IDs in dataset 1 not in dataset 2:')
#     print(diff1)
#     print()

    diff2 = list( set(id_list2) - set(id_list1))
    print('Basin IDs in dataset 2 not in dataset 1:')
    print(diff2)
    print('len(diff2) =', len(diff2))
    print()
    ## return diff2

#   compare_basin_ids()
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
    'BENEDICT', 'BENEDICTS',
    'AUGUST', 'AUGUSTINE',
    'CLEMENT', 'CLEMENTS',
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
def replace_state_letters( name ):

    if (name[-4:] == ' ARK'):     # don't match "PARK"
        name = name[:-4] + ' AR'
    if (name[-4:] == ' CAL'):
        name = name[:-4] + ' CA'
    if (name[-4:] == ' FLA'):
        name = name[:-4] + ' FL'
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
    if (name[-4:] == ' TEX'):
        name = name[:-4] + ' TX'
    if (name[-4:] == ' WYO'):
        name = name[:-4] + ' WY'
    #--------------------------------
    if (name[-5:] == ' ARIZ'):
        name = name[:-5] + ' AZ'
    if (name[-5:] == ' CONN'):
        name = name[:-5] + ' CT'
    if (name[-5:] == ' MICH'):
        name = name[:-5] + ' MI'
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
    if (name[-8:] == ' GEORGIA'):
        name = name[:-8] + ' GA'
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
    if (name[-10:] == ' TENNESSEE'):
        name = name[:-10] + ' TN'
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
    if (name == 'SAUK RIVER ABOVE WHITECHUCK RIVER NEAR DARRINGTON WA'):
        name = 'SAUK RIVER ABOVE WHITE CHUCK RIVER NEAR DARRINGTON WA'
    if (name == 'SHIAWASSEE RIVER AT BYRON MI'):
        name = 'SHIAWASSEE RIVER AT BATH ROAD AT BYRON MI'
    if ('SNTA BRB' in name):
        name = name.replace('SNTA BRB', 'SANTA BARBARA')
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
    return name
 
#   replace_city_letters()
#---------------------------------------------------------------------
def expand_basin_name( name, UPPER=True,
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
        # Call replace_periods() before replace_st().
        name = replace_st( name ) # STATE, SAINT, OR STREET
        #---------------------------------------------
        name = name.replace(' NR ',  ' NEAR ')
        name = name.replace(' BL ',  ' BELOW ')        
        name = name.replace(' AB ',  ' ABOVE ')
        name = name.replace(' BLW ', ' BELOW ')        
        name = name.replace(' ABV ', ' ABOVE ')
        #---------------------------------------------
        name = name.replace('EF ',   'EAST FORK ')
        name = name.replace('WF ',   'WEST FORK ')
        name = name.replace('SF ',   'SOUTH FORK ')
        name = name.replace('NF ',   'NORTH FORK ')
        name = name.replace('MF ',   'MIDDLE FORK ')
        name = name.replace('E F ',  'EAST FORK ')
        name = name.replace('W F ',  'WEST FORK ')
        name = name.replace('S F ',  'SOUTH FORK ')
        name = name.replace('N F ',  'NORTH FORK ')
        name = name.replace('M F ',  'MIDDLE FORK ')
        name = name.replace('S FK ', 'SOUTH FORK ')
        name = name.replace('N FK ', 'NORTH FORK ')
        name = name.replace('M FK ', 'MIDDLE FORK ')
        #---------------------------------------------
        # 'NB POTOMAC RIVER AT STEYER MD'
        # But next 2 also contain "NB" vs. "NB ":
        # 'GREENBRIER RIVER AT DURBIN, WV'
        # 'DRIFTWOOD RIVER NEAR EDINBURGH IND'
        #---------------------------------------------
        name = name.replace('NB ', 'NORTH BRANCH ')
        name = name.replace('SB ', 'SOUTH BRANCH ')             
        #---------------------------------------------            
        name = name.replace(' R ',   ' RIVER ')
        name = name.replace(' RV ',  ' RIVER ')
        name = name.replace(' RIV ', ' RIVER ')
        name = name.replace(' C ',   ' CREEK ')
        name = name.replace(' CK ',  ' CREEK ')
        name = name.replace(' CR ',  ' CREEK ')
        name = name.replace(' CRK ', ' CREEK ')
        # 'HAIKU STR NEAR HEEIA OAHU HI'
        name = name.replace(' STR ', ' STREAM ')
        name = name.replace(' B ',   ' BRANCH ')
        name = name.replace(' BR ',  ' BRANCH ')         
        name = name.replace(' CN ',    ' CANYON ')   ########### CHECK (OR CANAL)
        name = name.replace(' CNYN ',  ' CANYON ') 
        name = name.replace(' CYN ',   ' CANYON ')           
        name = name.replace(' TRIB ', ' TRIBUTARY ')
        name = name.replace(' FK ',   ' FORK ')
        name = name.replace(' WFORK ',' WEST FORK ')  # starts as "W.FORK"
        name = name.replace(' DIV ',  ' DIVERSION ')
        #---------------------------------------------
        name = name.replace(' O R N L ', ' ORNL ')
        name = name.replace(' OAK RIDGE NATL LAB ', ' ORNL ')
        #---------------------------------------------
        # 'G MIAMI RIVER AT HAMILTON OH'
        if (name[:2] == 'E '):
            name = 'EAST ' + name[2:]
        if (name[:2] == 'G '):
            name = 'GREAT ' + name[2:]
        # 'L BEAVER CREEK NEAR EAST LIVERPOOL OH'
        if (name[:2] == 'L '):
            name = 'LITTLE ' + name[2:]
        if (name[:2] == 'N '):
            name = 'NORTH ' + name[2:]
        if (name[:2] == 'S '):
            name = 'SOUTH ' + name[2:]
        if (name[:2] == 'W '):
            name = 'WEST ' + name[2:]
        #--------------------------------
        if (name[:3] == 'SO '):
            name = 'SOUTH ' + name[3:]
        if (name[:3] == 'ST '):
            name = 'SAINT ' + name[3:]   #########
        if (name[:3] == 'NO '):
            name = 'NORTH ' + name[3:]            
        #---------------------------------------------
        name = name.replace(' BYU ',    ' BAYOU ')   ######## CHECK       
        name = name.replace(' LK ',     ' LAKE ')
        name = name.replace(' MT ',     ' MOUNT ')    ##### MT CLEMENS
        name = name.replace(' RES ',    ' RESERVOIR ')
        name = name.replace(' RESV ',   ' RESERVOIR ')
        name = name.replace(' SPGS ',   ' SPRINGS ')
        name = name.replace(' SPRGS ',  ' SPRINGS ')
        name = name.replace(' SPRNGS ', ' SPRINGS ')
        name = name.replace(' SPS ',    ' SPRINGS ')  ##### CHECK
        #---------------------------------------------
        # e.g. COMPTON C A 120TH ST NR COMPTON CA
        name = name.replace(' A ',    ' AT ')    #############
        name = name.replace(' CTR ',  ' CENTER ')
        name = name.replace(' DR ',   ' DRIVE ')  ###### CHECK
        name = name.replace(' HWY ',  ' HIGHWAY ')
        name = name.replace(' HIGHWY ', ' HIGHWAY ')
        name = name.replace(' LN ',   ' LANE ')
        name = name.replace(' RD ',   ' ROAD ')
        name = name.replace(' RT ',   ' ROUTE ')
        name = name.replace(' RTE ',  ' ROUTE ')
        name = name.replace(' BRG ',  ' BRIDGE ')
        # ARCADIA WASH A GRAND AVE A ARCADIA CA
        name = name.replace(' AVE ',  ' AVENUE ')
        name = name.replace(' JCT ',  ' JUNCTION ')
        name = name.replace(' FT ',   ' FORT ')
        name = name.replace(' PT ',   ' POINT ')
        name = name.replace(' RNGR ', ' RANGER ')
        name = name.replace(' STA ',  ' STATION ')  ## "RANGER STA"
        name = name.replace(' STRM SWR ',  ' STORM SEWER ')
        # RAPID CREEK ABOVE WRF NR RAPID CITY, SD
        name = name.replace(' WRF ',  ' WHARF ')  ####### CHECK
    
        ## SEWAGE TREATMENT PL  (PL = plant or place)
          
        # name = name.replace()' SR ', ' STATE ROAD ')
    
        # In USGS_gauged as "Rd Crsg"; In CAMELS as "Rd"
        name = name.replace(' CRSG ', ' ')
        # name = name.replace(' CRSG ', ' CROSSING ')
        #---------------------------------------------
        name = name.replace(' RH ',   ' RANCH ')
        name = name.replace(' RCH ',  ' RANCH ')
        name = name.replace(' RNCH ', ' RANCH ')

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

#   expand_basin_name()
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
        name = expand_basin_name( name )   #####       
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
def collate( out_csv_file='all_basins.csv'):

    out_csv_path = repo_dir + out_csv_file
    out_csv_unit = open( out_csv_path, 'w')      #########
   
    #-------------------------------------------
    # CSV files for each of the basin datasets
    #-------------------------------------------
    ars_file    = get_data_filepath( 'ARS' )
    camels_file = get_data_filepath( 'CAMELS' )
    czo_file    = get_data_filepath( 'CZO' )
    gages2_file = get_data_filepath( 'GAGES2_all' )
    hlr_file    = get_data_filepath( 'HLR' )
    lter_file   = get_data_filepath( 'LTER' )
    mopex_file  = get_data_filepath( 'MOPEX' )
    neon_file   = get_data_filepath( 'NEON' )
    rfc_file    = get_data_filepath( 'RFC' )
    usgs_file   = get_data_filepath( 'USGS_gauged' )

#     csv_files = [
#       ars_file, camels_file, czo_file, gages2_file, hlr_file,
#       lter_file, mopex_file, neon_file, rfc_file, usgs_file ]   
#     for csv_file in csv_files:
#         csv_unit = open(csv_file, 'r')
               
#   collate()
#---------------------------------------------------------------------



