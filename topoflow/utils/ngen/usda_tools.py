
###### DOUBLE CHECK THE LATS & LONS  #######

# Copyright (c) 2023, Scott D. Peckham
#
# Jul 2023. Added: get_basin_repo_dir, get_ars_data_dir,
#           get_ars_basin_states, get_ars_url.
# May 2023. Wrote create_tsv, get_alternate_name,
#           get_watershed_notes.
#
#---------------------------------------------------------------------
#
#  get_basin_repo_dir()
#  get_ars_data_dir()
#  get_ars_state_codes()
#  get_ars_url()
#
#  write_tsv_header()
#  create_tsv()
#  get_alternate_name()
#  get_watershed_notes()
#
#---------------------------------------------------------------------
#
#  % python
#  >>> from topoflow.utils.ngen import usda_tools as ut
#  >>> ut.create_tsv( tsv_file='TEST.tsv' )
#
#---------------------------------------------------------------------
import numpy as np

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
def get_ars_data_dir():

    #-----------------------------------
    # Modify this directory as needed.
    #-----------------------------------
    repo_dir  = get_basin_repo_dir()
    data_dir  = repo_dir + 'USDA_ARS/'
    return data_dir
       
#   get_ars_data_dir()
#---------------------------------------------------------------------
def get_ars_state_codes():

    #---------------------------------------------
    # URLs may only work for these states.
    # Can get 2-letter state codes using:
    # get_state_code_map() in collate_basins.py.
    #---------------------------------------------
    
    #----------------------------------------------------  
    # These states have an associated URL for more info
    #----------------------------------------------------
    ars_state_codes = [
    'AL', 'AR', 'AZ', 'FL', 'GA', 'HI', 'IA', 'ID', 'IL',
    'IN', 'MD', 'MO', 'MS', 'NC', 'NE', 'NM', 'OH', 'OK',
    'PA', 'SD', 'TX', 'VA', 'VT', 'WA', 'WI', 'WV' ]

    return ars_state_codes

    #---------------------------------------
    # See: README_2.txt in USDA_ARS folder
    #---------------------------------------    
#     ars_states = [
#     'Arizona',                         
#     'Florida',                        
#     'Georgia',
#     'Hawaii',
#     'Idaho',
#     'Illinois',
#     'Iowa',
#     'Mississippi',
#     'Missouri',
#     'Nebraska',
#     'New Mexico',
#     'North Carolina',
#     'Ohio',
#     'Oklahoma',
#     'Pennsylvania',
#     'South Dakota',
#     'Texas',
#     'Vermont',
#     'Virginia',
#     'West Virginia',
#     'Wisconsin' ]
 
#   get_ars_state_codes()
#---------------------------------------------------------------------
def get_ars_url( state_code ):

    # Note that null string is also excluded.

    #--------------------------------------------------    
    # Note: These URLs no longer worked in July 2023.
    # e.g. https://hrsl.ba.ars.usda.gov/wdc/az.htm
    #--------------------------------------------------
    # Can use Wayback Machine by adding this to URL:
    # https://web.archive.org/web/20150923052139/
    #--------------------------------------------------
    # Some (e.g. Walnut Gulch) are now LTAR sites at:
    # https://ltar.ars.usda.gov/sites
    #--------------------------------------------------
    # HRSL = Hydrology and Remote Sensing Laboratory
    #--------------------------------------------------
    sc  = state_code.lower()
    url = 'https://hrsl.ba.ars.usda.gov/wdc/' + sc + '.htm'
    return url
        
#   get_ars_url()      
#---------------------------------------------------------------------
def write_tsv_header( out_unit, delim='\t' ):

    out_unit.write('Watershed ID'    + delim)
    out_unit.write('Location'        + delim)
    out_unit.write('State Code'      + delim)
    out_unit.write('Alternate Name'  + delim)
    out_unit.write('Watershed Name'  + delim)
    out_unit.write('Latitude (DMS)'  + delim)
    out_unit.write('Longitude (DMS)' + delim)
    out_unit.write('Latitude (dec deg)'  + delim)
    out_unit.write('Longitude (dec deg)' + delim)
    out_unit.write('Acres'        + delim)
    out_unit.write('Area (km2)'   + delim)
    out_unit.write('Start Date'   + delim)
    out_unit.write('End Date'     + delim)
    out_unit.write('ARS DB Years' + delim)
    out_unit.write('Data URL'     + delim)
    out_unit.write('Notes')
    out_unit.write('\n')
 
   
#   write_tsv_header()
#---------------------------------------------------------------------
def create_tsv( data_dir=None, tsv_file='ARS_db.tsv', delim='\t' ):

    #------------------------------------------------------
    # Sometimes there are multiple lines in INDEX.TXT for
    # a given Watershed ID.  Sometimes the lat/lon values
    # are different.  How best to handle this?
    #------------------------------------------------------
    if (data_dir is None):
        data_dir = get_ars_data_dir()   

    index_path = data_dir + 'INDEX.TXT'
    tsv_path   = data_dir + tsv_file
      
    in_unit  = open(index_path, 'r')
    out_unit = open(tsv_path,   'w')

    #----------------------------------------------------  
    # These states have an associated URL for more info
    #----------------------------------------------------
    state_codes = get_ars_state_codes()

    #----------------------------------------
    # Write column headers for new TSV file
    #----------------------------------------
    write_tsv_header( out_unit )
       
    #--------------------------------------
    # Skip over header lines in INDEX.TXT
    #--------------------------------------
    for k in range(8):
        line = in_unit.readline()
        
    while (True):
        line = in_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        if (line == '\n'):  # to avoid indent
            continue  # skip rest of loop
        #----------------------------------------------------------
        # Sometimes there are multiple lines in INDEX.TXT for a
        # given Watershed ID, but 1st few columns are left blank.
        # Sometimes lat/lon or acres values differ.
        # If left blank, sorting on columns will not work.
        #----------------------------------------------------------
        # Note: If the line length does not support the range
        #       of indices, then a null string is returned.
        #       See Watershed IDs = 54001 and 54002.
        #----------------------------------------------------------
        watershed_ID   = line[0:5]
        location       = line[5:20].strip()
        state_code     = line[20:22]
        watershed_name = line[23:36].strip()
        #------------------------------------------ 
        if (watershed_ID.strip() == ''):
            watershed_ID = last_watershed_ID
        if (location == ''):
            location = last_location 
        if (state_code.strip() == ''):
            state_code = last_state_code
        if (watershed_name == ''):
            watershed_name = last_watershed_name   
        #------------------------------------------                  
        out_unit.write( watershed_ID + delim)   # Watershed ID
        out_unit.write( location     + delim)   # Location
        out_unit.write( state_code   + delim)    # State Code
        #--------------------------------------------
        # Add the "Alt. Name" or specific location
        #--------------------------------------------
        wn = np.int32(watershed_ID)
        alt_name = get_alternate_name( wn )
        out_unit.write( alt_name + delim)   # Alt Name or placeholder
        #------------------------------------------
        if (wn == 9002):
            watershed_name = 'W-II'  # See Note column. 
        out_unit.write( watershed_name + delim)  # Watershed Name
        #--------------------------------------
        last_watershed_ID   = watershed_ID
        last_location       = location
        last_state_code     = state_code
        last_watershed_name = watershed_name          
        #-------------------------------------------
        # Write Latitude (DMS) and Longitude (DMS)
        #--------------------------------------------
        # A lat dms string is always 6 characters,
        #   since degrees can only have 2 digits.
        # A lon dms string may be 6 or 7 characters,
        #   since degrees may be 2 or 3 digits.
        #--------------------------------------------
        lat_dms_str = line[76:82].strip()     # always 6 characters #########
        out_unit.write( lat_dms_str + delim)  # Latitude (DMS)
        lon_dms_str = line[83:90].strip()     # 6 or 7 characters
        #--------------------------------------------------
        # Fix an error in the lon data.  See Note column.
        # lon_dms_str may be 7 chars, starting with "0",
        # so get last 4 as shown.
        #--------------------------------------------------
        if (wn >= 69051) and (wn <= 69053):
            lon_dms_str = '97' + lon_dms_str[-4:]  # 98 -> 97

        #----------------------------------------------------------- 
        # In CONUS, all lats>0, and all lons<0.  But in the file
        # INDEX.TXT, DMS lons don't have negative sign, so add it.
        #-----------------------------------------------------------       
        if (lon_dms_str != ''):
            # Longitudes west of prime meridian are negative.
            out_unit.write( '-' + lon_dms_str + delim)  # Longitude (DMS) 
        else:
            out_unit.write( lon_dms_str + delim)  # Longitude (DMS)
       
        #--------------------------------------------
        # Convert latitude (DMS) to decimal degrees
        # May be missing/null.
        # In at least 4 cases, the S part is "--"
        #--------------------------------------------
        lat_dec_deg_str = ''  # default
        if (lat_dms_str.strip() != ''):        
            D = np.float64( lat_dms_str[0:2] )
            M = np.float64( lat_dms_str[2:4] )
            ## S = np.float64( lat_dms_str[4:6] )
            S_str = lat_dms_str[4:6]
            S = 0.0  # default
            if (S_str != '--'):
                S = np.float64( S_str )
            if (D < 0):
               lat_dec_deg = D - (M/60) - (S/3600)  #### NEED THIS
            else:
               lat_dec_deg = D + (M/60) + (S/3600)
            #----------------------------------------------
            # Should only need 5 digits after the decimal
            # See: Wikipedia - Decimal Degrees
            #----------------------------------------------
            lat_dec_deg_str = str( round(lat_dec_deg, 5) )
        out_unit.write( lat_dec_deg_str + delim)
        
        #---------------------------------------------
        # Convert longitude (DMS) to decimal degrees
        # D may be 2 or 3 digits.
        # In at least 4 cases, the S part is "--"
        #---------------------------------------------
        lon_dec_deg_str = ''  # default
        if (lon_dms_str.strip() != ''):
            if (len(lon_dms_str) == 6):
                i1 = 2
            else:
                i1 = 3
            D = np.float64( lon_dms_str[0:i1] )
            M = np.float64( lon_dms_str[i1:i1+2] )
            ## S = np.float64( lon_dms_str[i1+2:i1+4] )
            S_str = lon_dms_str[i1+2:i1+4]
            S = 0.0  # default
            if (S_str != '--'):
                S = np.float64( S_str )
            if (D < 0):
                lon_dec_deg = D - (M/60) - (S/3600)   #### NEED THIS
            else:
                lon_dec_deg = D + (M/60) + (S/3600)
            #----------------------------------------------
            # Should only need 5 digits after the decimal
            # See: Wikipedia - Decimal Degrees
            #----------------------------------------------
            lon_dec_deg_str = str( round(lon_dec_deg, 5) )
            # Longitudes west of prime meridian are negative.
            lon_dec_deg_str = '-' + lon_dec_deg_str
        out_unit.write( lon_dec_deg_str + delim)           
        #-----------------------------------------------------
        # Note: Need strip() here to remove trailing newline
        # char for Watershed IDs: 54001 and 54002.
        #-----------------------------------------------------            
        acre_str = line[36:47].strip()
        out_unit.write( acre_str + delim)  # Acres
        
        #-----------------------------------------------
        # Convert acres to square km & add Area column
        #-----------------------------------------------
        area_str = ''
        if (acre_str != ''):
            area = np.float32(acre_str) * 0.00404686 # [sq km]
            area_str = str( round(area,4) )
        out_unit.write( area_str + delim)  # Area ()

        #-----------------------------------------------
        # Convert start date string to ISO date format
        # All dates in file are before 1998.
        #-----------------------------------------------
        s = line[47:55]
        start_date_iso = ''
        if (s != ''):
            start_date_iso = '19' + s[-2:] + '-' + s[:2] + '-' + s[3:5]
        out_unit.write( start_date_iso + delim)  # Start Date
    
        #---------------------------------------------
        # Convert end date string to ISO date format
        # End date string may be "Present"
        #---------------------------------------------
        e = line[56:64].strip()
        end_date_iso = ''
        if (e != ''):
            if (e.lower() != 'present'):
                end_date_iso = '19' + e[-2:] + '-' + e[:2] + '-' + e[3:5]
            else:
                end_date_iso = 'Present'
        out_unit.write( end_date_iso + delim)  # End Date  
        
        #---------------------------------------------------
        # Write date range available in ARS Water Database
        # May be missing/null
        #---------------------------------------------------
        out_unit.write( line[65:76].strip() + delim)  # ARS DB Years

        #--------------------------
        # Write a "more data" URL
        #--------------------------
        if (state_code in state_codes):
            url = get_ars_url( state_code )
            out_unit.write( url + delim )
        else:
            out_unit.write( '' + delim )

        #-----------------------------
        # Write any additional notes
        #-----------------------------
        note = get_watershed_notes( wn )
        out_unit.write( note )
        ## out_unit.write( note + delim )
                     
        #----------------------------
        # Write a newline character
        #----------------------------
        out_unit.write('\n')

    #-------------------------------
    # Close input and output files
    #-------------------------------
    in_unit.close()
    out_unit.close()
                                       
#   make_new_csv()
#---------------------------------------------------------------------
def get_alternate_name( wn ):

    #-------------------------------------------------------
    # Note: These names come from the file "inventor.txt".
    #       That file has watersheds w/ data collection as
    #       of the 1994 to 1995 timeframe.
    #-------------------------------------------------------
    if (wn == 16001):
        alt_name = 'Mahantango Cr. at Malta'
    elif (wn == 16003):
        alt_name = 'Mahantango Cr. at Klgtn.'
    elif (wn == 16006):
        alt_name = 'Stehr Bros. Weir'
    elif (wn >= 56001) and (wn <= 56004):
        alt_name = 'Moscow'
    elif (wn >= 62001) and (wn <= 62019):
        alt_name = 'Pigeon Roost'
    elif (wn >= 63001) and (wn <= 63015):
        alt_name = 'Walnut Gulch'
    elif (wn >= 63101) and (wn <= 63106):
        alt_name = 'Lucky Hills'
    elif (wn >= 63111) and (wn <= 63113):
        alt_name = 'Kendall'
    elif (wn == 63121):
        alt_name = 'Holiday Ranch'
    elif (wn >= 63122) and (wn <= 63124):
        alt_name = 'Cowen'
    elif (wn >= 63201) and (wn <= 63223):
        alt_name = 'Walnut Gulch Pond'
    elif (wn == 68001):
        alt_name = 'Reynolds Creek Outlet'
    elif (wn == 68002):
        alt_name = 'Salmon Creek'
    elif (wn == 68003):
        alt_name = 'Reynolds Creek'
    elif (wn == 68004):
        alt_name = 'Reynolds Creek at Tollgate'
    elif (wn >= 68011) and (wn <= 68012):
        alt_name = 'Reynolds Creek'
    elif (wn == 68013):
        alt_name = 'Reynolds Mtn. East'
    elif (wn == 68014):
        alt_name = 'Lower Sheep Creek'
    elif (wn >= 68015) and (wn <= 68020):
        alt_name = 'Reynolds Creek'
    elif (wn == 68021):
        alt_name = 'Upper Sheep Creek 1'
    elif (wn == 68022):
        alt_name = 'Upper Sheep Creek 2'
    elif (wn == 68023):
        alt_name = 'Flats Watershed'
    elif (wn == 68024):
        alt_name = "Nancy's Gulch"
    elif (wn >= 68033) and (wn <= 68034):
        alt_name = 'Reynolds Creek'                      
    elif (wn >= 74002) and (wn <= 74011):
        alt_name = 'Little River'
    elif (wn >=83001) and (wn <= 83014):
        alt_name = 'Goodwin Creek'
    else:
        alt_name = 'Alternate_Name'
        
    ## alt_name = 'Alternate_Name'
    
    return alt_name
   
#   get_alternate_name()
#---------------------------------------------------------------------
def get_watershed_notes( wn ):

    #---------------------------------------------
    # Note:  Don't use commas in notes.  Add one
    #        space at end of each note line.
    #---------------------------------------------
    note = ' '  # default
    if (wn == 9001):
        note  = 'Notice that first watershed name is W-1 '
        note += 'but others use roman numerals.'
    if (wn == 9002):
        note  = 'Watershed name in INDEX.TXT is W-I but according '
        note += 'to info at URL it should be W-II.'
    if (wn == 68021):
        note  = 'Lats and lons in INDEX.TXT and inventor.txt do '
        note += 'not match exactly but number of acres is an '
        note += 'exact match to Upper Sheep Creek 1.' 
    if (wn == 68022):
        note  = 'Lats and lons in INDEX.TXT and inventor.txt do '
        note += 'not match exactly but number of acres is an '
        note += 'exact match to Upper Sheep Creek 2.'     
    if (wn >= 69051) and (wn <= 69053):
        #---------------------------------------------------------
        # Note: In Google Maps go to: 34°56'41.0"N 97°57'08.0"W.
        #       Location is labeled "Ninnekah" and is a bridge
        #       on "Old Highway 81" that crosses a tributary
        #       that flows east to join the Washita River.
        #---------------------------------------------------------
        # Note: Other Little Washita lons may also be wrong.
        #---------------------------------------------------------        
        note  = 'Longitude DMS in INDEX.txt does not match value in inventor.txt. '
        note += 'From Google Maps degrees should be -97 not -98. '
        note += 'Outlet is at a bridge near the town of Ninnekah.'

    return note

#   get_watershed_notes()
#---------------------------------------------------------------------



            
