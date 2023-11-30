
# Copyright (c) 2023, Scott D. Peckham
#
# Nov 2023. Incorporated info in new RFC data sets.  Wrote
#           merge_rfc_basin_info().  Wrote convert_json_to_csv()
#           to work around issues with WGRFC shapefile.
#           Can now include RFC_name and basin bounding box
#           for 7759 basins in 13 NOAA RFCs.
#           Wrote:  fix_nws_ids_for_MBRFC().
#           Wrote:  cull_usa_dcp_data_rows().
#           Wrote:  get_station_data_from_api().
#           Still need to update create_tsv().
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
#  cull_usa_dcp_data_rows()
#  fix_nws_ids_for_MBRFC()        ####
#  convert_json_to_csv()
#  merge_rfc_basin_info()
#
#  create_tsv()
#  sort_by_site_code()
#
#---------------------------------------------------------------------

import numpy as np
import os, os.path
import json, requests, time

import pickle    # to save usgs station info map
from topoflow.utils.ngen import collate_basins as cb
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
def get_rfc_data_dir( ROOT=False ):

    #-----------------------------------
    # Modify this directory as needed.
    #-----------------------------------
    repo_dir  = get_basin_repo_dir()
    data_dir  = repo_dir + 'NOAA_RFCs/'
    if not(ROOT):
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
            if (lat != 'NA'):  ############## ADD MISSING INFO LATER #######
                htype_dict[ sid ] = \
                    {'htype':htype, 'rfc_name':rfc, 'huc8':huc8,
                     'lat':lat, 'lon':lon }       
        file_unit.close() 

    return htype_dict 
    
#   get_hydrograph_type_dict()
#---------------------------------------------------------------------
def cull_usa_dcp_data_rows( in_csv_file=None, out_csv_file=None):

    #--------------------------------------------------------------
    #  Note: This website at Iowa State Univ. provides DCP data
    #        for many networks that can be downloaded as CSV
    #        or as a shapefile:
    #  https://mesonet.agron.iastate.edu/sites/networks.php
    #--------------------------------------------------------------
    # Downloaded a CSV with DCP info for all All Networks.
    # This function culls just the ones for US states.
    #--------------------------------------------------------------
    # If viewing HTML Table for a state, it has an "Archive Ends"
    # column for inactive DCPs.  However, if saving as CSV, only
    # ACTIVE DCPs are included and this column is not provided.
    #--------------------------------------------------------------
    # It is likely that many of these DCPs are not of the type
    # "Stream", vs. "Atmosphere", "Well", etc.
    #--------------------------------------------------------------  
    # Wrote an email to Daryl Herzmann (akrherz@iastate.edu) 
    # on 2023-11-09 to request more info.
    # See: https://mesonet.agron.iastate.edu/info/contacts.php
    #--------------------------------------------------------------           
    if (in_csv_file is None):
        in_csv_file = 'All_Networks_Active_DCP_Data.csv'
    data_dir  = get_rfc_data_dir( ROOT=True )
    data_dir += 'DCP_Network_Data_from_IEM/'
    in_csv_path = data_dir + in_csv_file
    #-------------------------------------
    if (out_csv_file is None):
        out_csv_file = 'US_State_Active_DCP_Data.csv'
    out_csv_path = data_dir + out_csv_file
    #-------------------------------------    
    in_csv_unit  = open( in_csv_path,  'r' )
    out_csv_unit = open( out_csv_path, 'w' )  
    comma  = ','
    n_rows_in   = 0
    n_rows_out  = 0
    n_rows_skip = 0
    #------------------------------------------------
    # Skip these territories, etc.
    # Keep PR (Puerto Rico) and VI (Virgin Islands)
    # GU_DCP is Guam.
    #------------------------------------------------
    skip_list  = ['GU_DCP','P1_DCP','P2_DCP','P3_DCP','P4_DCP']
    state_dict = dict()
    #---------------------------------------
    # Copy header from in_file to out_file
    #---------------------------------------
    header = in_csv_unit.readline()
    out_csv_unit.write( header ) 
            
    #--------------------------------------
    # Read each data line in the CSV File
    # All US DCP networks end with "_DCP"
    #--------------------------------------
    # Canadian networks have names like:
    #    CA_AB_DCP (Canada, Alberta)
    # Mexican networks have names like:
    #    MX_BJ_DCP (Mexico, Baja)
    # Puerto Rico is included:  PR_DCP
    # Guam is included: GU_DCP
    # Barbuda is included: AG__DCP
    # Others are BD_DCP, BM_DCP (Bermuda)
    #--------------------------------------                                       
    while (True):
        line = in_csv_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        n_rows_in += 1
        #-------------------------------------
#         if (line.endswith('_DCP\n')):
#             out_csv_unit.write( line )
#             n_rows_out += 1
        #-------------------------------------        
        vals = line.split( comma )
        name = vals[-1].strip()  # last entry
        KEEP = not(name in skip_list)
        if (KEEP and name.endswith('_DCP') and (len(name) == 6)):
            #--------------------------------------
            # e.g. AK_DCP, CO_DCP, GU_PCP, PR_DCP
            # Could exclude Guam, others.
            #--------------------------------------
            out_csv_unit.write( line )
            n_rows_out += 1
            #--------------------------------------
            state_code = name[0:2]
            try:
                state_dict[ state_code ] += 1
            except:
                state_dict[ state_code ] = 0
        else:
            #--------------------------------------
            # Skip Canada, Mexico, & Barbuda DCPs
            #--------------------------------------
            n_rows_skip += 1  
    
    in_csv_unit.close()
    out_csv_unit.close()
 
    #-------------------------------------   
    # Print the count of each state code
    #-------------------------------------
    key_list = list( state_dict.keys())
    key_list.sort()
    for key in key_list:
        print('# ' + key, state_dict[key])
    print()

    print('Read   ', n_rows_in,  'rows.')
    print('Skipped', n_rows_skip,'rows.')
    print('Wrote  ', n_rows_out, 'rows.')
    print('Finished.')
    print()

#   cull_usa_dcp_data_rows()
#---------------------------------------------------------------------
def fix_nws_ids_for_MBRFC():

    #-------------------------------------------------------
    # Read: Station-ID_to_RFC_Crosswalk.csv
    # If RFC column is MBRFC and ID is numeric, then
    #    Read the string in the Name column
    #    Expand the string (USGS-style)
    #    Find same string in new_rfc_hads_sorted.tsv,
    #       and get corresponding NWS_loc_ID.
    #    Save new NWS IDs to: 
    #-------------------------------------------------------
    # This does not fix all of the numeric NWS IDs.
    #-------------------------------------------------------
    # Sorting the crosswalk on the ID column shows that
    # many numeric IDs are repeated w/ same row info
    # except for OBJECTID.
    # There are 1605 NWS IDs for MBRFC that are numeric.
    # IDs start at 101 and go up to 3297; many sequential.
    # Note that crosswalk also provides the CWA code.
    # In crosswalk, are lon/lat for outlet or centroid?
    #-------------------------------------------------------
    pass
    
#   fix_nws_ids_for_MBRFC()
#---------------------------------------------------------------------
def convert_json_to_csv( json_file=None, csv_file=None, 
                         rfc_name='WGRFC', REPORT=True ):

    #-------------------------------------------------------------
    # Note: The shapefile for the West Gulf RFC (WGRFC) had an
    #       empty attribute table.  However, was able to find
    #       a JSON file on an ESRI website with the attributes.
    #       The URL was:
    #       https://www.arcgis.com/home/item.html?
    #           id=fcbec367432d4b5a829a250906189fd6
    #       Then clicked on "Source: Feature Collection",
    #       downloaded to text and saved with JSON extension.
    #       This JSON file has info for 624 basins.
    #-------------------------------------------------------------
    #       It may be possible to access info for other RFCs in
    #       this same way, which gives additional attributes.
    #-------------------------------------------------------------    
    if (json_file is None):
        json_file = rfc_name + '_Data_from_ESRI_Online.json'
    data_dir  = get_rfc_data_dir( ROOT=True )
    data_dir += 'Basin_Shapefiles_by_RFC/'
    json_path = data_dir + json_file
    #-------------------------------------
    if (csv_file is None):
        csv_file = rfc_name + '_data.csv'
    csv_path = data_dir + csv_file
    #-------------------------------------    
    json_unit = open( json_path, 'r' )
    csv_unit  = open( csv_path,  'w' )     #########
    #-------------------------------------     
    data_dict    = json.load( json_unit )
    layer_list   = data_dict['layers']
    layer0_dict  = layer_list[0]
    featSet_dict = layer0_dict['featureSet']
    feature_list = featSet_dict['features']

    if (REPORT):
        print('data_dict keys    =', list(data_dict.keys()))
        print('len(layer_list)   =', len(layer_list))  # Is "1" here.
        print('layer0_dict keys  =', list(layer0_dict.keys()))
        print('featSet_dict keys =', list(featSet_dict.keys()))
        print('len(feature_list) =', len(feature_list))  # Is "624" here.        
        print('')

    #--------------------------------
    # Write header for new CSV file
    #--------------------------------
    delim  = ',' 
    header = ''
    header += 'RFC_ID'     + delim
    header += 'NWS_loc_ID' + delim
    header += 'MINLON'     + delim    # MIN_X_AXIS
    header += 'MAXLON'     + delim    # MAX_X_AXIS
    header += 'MINLAT'     + delim    # MIN_Y_AXIS
    header += 'MAXLAT'     + delim    # MAX_Y_AXIS
    header += 'NAME'       + '\n'     # newline at end
    csv_unit.write( header )   

    n_basins = 0
    for feature in feature_list:
        #--------------------------------------
        # Each feature is a Python dictionary
        # Keys are: 'geometry', 'attributes'.
        #--------------------------------------
        att_dict = feature['attributes']
        if (REPORT and (n_basins == 0)):
            print('feature_dict keys =', list(feature.keys()) )
            print('att_dict keys     =', list(att_dict.keys()) )
            print()

        name   = att_dict['NAME']        # string
        nws_id = att_dict['CH5_ID']      # string
        minlon = att_dict['MIN_X_AXIS']  # float
        maxlon = att_dict['MAX_X_AXIS']
        minlat = att_dict['MIN_Y_AXIS']
        maxlat = att_dict['MAX_Y_AXIS']
        n_basins += 1      

        #---------------------------------
        # Format the lat and lon strings
        #---------------------------------
        minlon_str = '{x:.5f}'.format(x=minlon)
        maxlon_str = '{x:.5f}'.format(x=maxlon)
        minlat_str = '{x:.5f}'.format(x=minlat)
        maxlat_str = '{x:.5f}'.format(x=maxlat)
            
        #-----------------------------       
        # Write info to new CSV file
        #-----------------------------
        csv_line = ''
        csv_line += rfc_name   + delim
        csv_line += nws_id     + delim
        csv_line += minlon_str + delim
        csv_line += maxlon_str + delim
        csv_line += minlat_str + delim
        csv_line += maxlat_str + delim
        csv_line += name + '\n'      # newline at end
        csv_unit.write( csv_line )
        
    json_unit.close()     
    csv_unit.close()
    print('Wrote info for', n_basins, 'basins.')
    print('Finished.')
    print()
   
#   convert_json_to_csv()
#---------------------------------------------------------------------
def merge_rfc_basin_info( tsv_file='new_rfc_info.tsv', REPORT=True):

    #----------------------------------------------------------------
    # Note: Found a site that provides a shapefile of basins for
    #       each of the 13 NOAA RFCs at:
    #          https://www.nohrsc.noaa.gov/gisdatasets/
    #       Exported shapefile attributes to CSV from QGIS.
    #       The 13 shapefiles do not provide the same set of
    #       attributes but they have several in common, such as
    #       the bounding box info for each basin.
    #       This also allows us to map 5-char NWS IDs to RFC names.
    #----------------------------------------------------------------
    #       The West Gulf RFC shapefile has no attributes and also
    #       seems to not have all of the basin boundaries.
    #----------------------------------------------------------------
    #       The Missouri Basin RFC has numbers instead of 5-char
    #       NWS location IDs for many basins.  On this page:
    #       https://water.weather.gov/ahps/region.php?rfc=mbrfc
    #       after clicking the radio button "River Observations",
    #       it says there are 1185 total gauges.  But the shapefile
    #       has info for 1412 basins.
    #----------------------------------------------------------------    
    data_dir  = get_rfc_data_dir( ROOT=True )
    data_dir += 'Basin_Shapefiles_by_RFC/'

    rfc_list = [
    'abrfc', 'aprfc', 'cbrfc', 'cnrfc', 'lmrfc',
    'marfc', 'mbrfc', 'ncrfc', 'nerfc', 'nwrfc',
    'ohrfc', 'serfc']
    #### 'ohrfc', 'serfc', 'wgrfc']   # wgrfc shapefile is corrupt
    prefix_list = ['b_' + s for s in rfc_list ]
    csv_paths = list()
    
    for prefix in prefix_list:
        csv_path = data_dir + prefix + '/' + prefix + '.csv'
        csv_paths.append( csv_path )

    #-----------------------------
    # Open new TSV file to write
    #-----------------------------
    tsv_path = data_dir + tsv_file
    tsv_unit = open( tsv_path, 'w')

    #--------------------------------
    # Write header for new TSV file
    #--------------------------------
    tab    = '\t'  # tab character  
    header = ''
    header += 'RFC_ID'     + tab
    header += 'NWS_loc_ID' + tab
    header += 'MINLON'     + tab    # MIN_X_AXIS_
    header += 'MAXLON'     + tab    # MAX_X_AXIS_
    header += 'MINLAT'     + tab    # MIN_Y_AXIS_
    header += 'MAXLAT'     + tab    # MAX_Y_AXIS_
    header += 'NAME'       + '\n'   # newline at end
    tsv_unit.write( header )   
       
    k = 0
    comma = ','
    n_total = 0
 
    for prefix in prefix_list:
        rfc_name = rfc_list[k].upper()
        csv_unit = open( csv_paths[k], 'r' )
        if (REPORT):
            print('For RFC = ' + rfc_name + ':')
            
        #--------------------------------------
        # Get info from the CSV header line
        #--------------------------------------
        # Find col nums of certain attributes
        # NCRFC has many other attributes.
        # Several have centroid lon, lat:
        #   e.g. NCRFC, NERFC, 
        #--------------------------------------
        n_basins = 0
        header   = csv_unit.readline()
        headings = header.split( comma )
        if (REPORT):
            print('   number of columns =', len(headings) )
        if (len(headings) <= 1):
            print('   ### ERROR: No data found for this RFC.')
            break
        name_col   = headings.index('NAME')
        minlon_col = headings.index('MIN_X_AXIS_')
        maxlon_col = headings.index('MAX_X_AXIS_')
        minlat_col = headings.index('MIN_Y_AXIS_')
        maxlat_col = headings.index('MAX_Y_AXIS_')
        if ('NWS_ID' in headings):
            nwsid_col = headings.index('NWS_ID')
        elif ('NWSID' in headings):
            nwsid_col = headings.index('NWSID')  # LMRFC
        elif ('CH5_ID' in headings):
            nwsid_col = headings.index('CH5_ID')
        elif ('BASIN' in headings):
            nwsid_col = headings.index('BASIN')  # CNRFC
        elif ('NAME' in headings):
            # For MBRFC, NAME is NWS_ID + "UPR" or "LWR"
            # It also has SHAPE_LENG, SHAPE_AREA,
            #  CENTROID_L (lat), CENTROID_1 (lon)
            nwsid_col = headings.index('NAME')  # MBRFC
        else:
            break  # WGRFC (West Gulf) has no attributes.

        #--------------------------------------
        # Read each data line in the CSV File
        #--------------------------------------                                   
        while (True):
            csv_line = csv_unit.readline()
            if (csv_line == ''):
                break  # (reached end of file)
            vals = csv_line.split( comma )
                                                              
            nws_loc_id = vals[ nwsid_col ].strip()
            name       = vals[ name_col ].strip()
            if (len(nws_loc_id) == 0):
                nws_loc_id = name  # this occurs for OHRFC
            #----------------------------------------
            # Get geographic bounding box for basin
            #----------------------------------------
            minlon = np.float64( vals[ minlon_col ].strip() )
            maxlon = np.float64( vals[ maxlon_col ].strip() )
            minlat = np.float64( vals[ minlat_col ].strip() )
            maxlat = np.float64( vals[ maxlat_col ].strip() )
            #---------------------------------
            # Format the lat and lon strings
            #---------------------------------
            minlon_str = '{x:.5f}'.format(x=minlon)
            maxlon_str = '{x:.5f}'.format(x=maxlon)
            minlat_str = '{x:.5f}'.format(x=minlat)
            maxlat_str = '{x:.5f}'.format(x=maxlat)
                
            #-----------------------------       
            # Write info to new TSV file
            #-----------------------------
            tsv_line = ''
            tsv_line += rfc_name   + tab
            tsv_line += nws_loc_id + tab
            tsv_line += minlon_str + tab
            tsv_line += maxlon_str + tab
            tsv_line += minlat_str + tab
            tsv_line += maxlat_str + tab
            tsv_line += name + '\n'      # newline at end
            tsv_unit.write( tsv_line )
            n_basins += 1
            n_total  += 1

        if (REPORT):
            print('   number of basins  =', n_basins)
        csv_unit.close()
        k += 1

    #-----------------------------------------------
    # Add info for the WGRFC by a different method
    #-----------------------------------------------
    n_basins = 0
    rfc_name       = 'WGRFC'
    wgrfc_csv_file = 'WGRFC_data.csv'
    wgrfc_csv_path = data_dir + wgrfc_csv_file
    csv_unit = open( wgrfc_csv_path, 'r' )
    header   = csv_unit.readline()
    if (REPORT):
        print('For RFC = ' + rfc_name + ':')
        print('   by reading info from JSON file')
    #--------------------------------------
    # Read each data line in the CSV File
    #--------------------------------------                                   
    while (True):
        csv_line = csv_unit.readline()
        if (csv_line == ''):
            break  # (reached end of file)
        vals = csv_line.split( comma )
                                                          
        nws_loc_id = vals[1].strip()
        minlon_str = vals[2].strip()  # already formatted
        maxlon_str = vals[3].strip()
        minlat_str = vals[4].strip()
        maxlat_str = vals[5].strip()        
        name       = vals[6].strip()
        #-----------------------------       
        # Write info to new TSV file
        #-----------------------------
        tsv_line = ''
        tsv_line += rfc_name   + tab
        tsv_line += nws_loc_id + tab
        tsv_line += minlon_str + tab
        tsv_line += maxlon_str + tab
        tsv_line += minlat_str + tab
        tsv_line += maxlat_str + tab
        tsv_line += name + '\n'      # newline at end
        tsv_unit.write( tsv_line )
        n_basins += 1
        n_total  += 1

    csv_unit.close()
    if (REPORT):
        print('   number of basins  =', n_basins)
    tsv_unit.close()
    print('Total number of RFC basins =', n_total)
    print('Finished merging RFC basin info to create:')
    print('   ' + tsv_path)
    print()
         
#   merge_rfc_basin_info()
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
def get_station_data_from_api( gauge_id='13334300', REPORT=True):

    #-----------------------------------------------------------
    # Note: This function uses a new (beta) API from NOAA
    #       to retrieve information for a given gauge.
    #       In the API, gauge_ID can be USGS ID or NWS Loc ID.
    #-----------------------------------------------------------
    # Q:  What is best way to check if this is stream gauge?
    # A:  The PEDTS Code seems to contain this information.
    #-----------------------------------------------------------
    # The "PEDTS Code" is explained in the NWS SHEF Manual.
    # The first 2 letters are "Physical Element" (PE) code.
    #    HG stands for "gauge height".
    #    QR stands for "discharge".
    #    PP stands for "precipitation".
    # The 3rd letter is the Duration Code and can be:
    #    I (instantaneous), H (hourly), D (daily), etc.
    # The 4th letter is the Type Code and can be:
    #    R (reading/observed), F (forecast), etc.
    # The 5th letter is the Source Code and is used to
    #    specify how the data was created or transmitted.
    # Note that "RG" = GOES observation.
    # Note that "FZ" = nonspecific forecast data (default)
    # See Tables 3 and 4 in the NWS SHEF Manual for details.
    #-----------------------------------------------------------    
    url = 'https://preview-api.water.noaa.gov/v1/gauges/' + gauge_id
    response  = requests.get( url )  # returns JSON as dict
    data_dict = response.json()

    #------------------------------------------------------------------------------
    # Note: If the API can't find a gauge ID, it returns this:
    # {'code': 5, 'message': '[01016500] could not find gauge ID', 'details': []}
    #------------------------------------------------------------------------------
    if ('code' in data_dict):   # Python 3 syntax
        if (data_dict['code'] == 5):
            if (REPORT):
                print('### SORRY, API could not find gauge_ID:', gauge_id)
            return None

    #------------------------------------------------
    # Extract some info and store in new dictionary
    #------------------------------------------------
    new_dict = dict()
    try:
        new_dict['usgs_id'] = data_dict['usgsId']
    except:
        #------------------------------------
        # Still need this for USGS_NWIS_all
        #------------------------------------
        new_dict['usgs_id'] = gauge_id   
    new_dict['nws_loc_id']   = data_dict['lid']
    new_dict['upstrm_lid']   = data_dict['upstreamLid']   
    new_dict['dnstrm_lid']   = data_dict['downstreamLid']
    new_dict['nwm_reach_id'] = data_dict['reachId']
    new_dict['rfc_abbrev']   = data_dict['rfc']['abbreviation']
    new_dict['wfo_abbrev']   = data_dict['wfo']['abbreviation']
    new_dict['usgs_name']    = data_dict['name']
    new_dict['description']  = data_dict['description']
    new_dict['state_code']   = data_dict['state']['abbreviation']
    new_dict['county_name']  = data_dict['county']
    new_dict['longitude']    = data_dict['longitude']
    new_dict['latitude']     = data_dict['latitude']
    new_dict['elevation']    = data_dict['datum']['elevation']
    new_dict['updated_time'] = data_dict['status']['updatedTime']
    if (new_dict['updated_time'] is None):
        new_dict['updated_time'] = ''
    new_dict['in_service']   = data_dict['inService']['enabled']
    new_dict['pedts_obs']    = data_dict['pedts']['observed'] 
    new_dict['pedts_pred']   = data_dict['pedts']['forecast']
                 
    if (REPORT):
        print('USGS ID        =', new_dict['usgs_id'])
        print('NWS LID        =', new_dict['nws_loc_id'])
        print('Upstream LID   =', new_dict['upstrm_lid'])
        print('Downstream LID =', new_dict['dnstrm_lid'])
        print('NWM Reach ID   =', new_dict['nwm_reach_id'])
        print('RFC name       =', new_dict['rfc_abbrev'])
        print('WFO name       =', new_dict['wfo_abbrev'])
        print('USGS name      =', new_dict['usgs_name'])
        print('Description    =', new_dict['description'])            
        print('US State Code  =', new_dict['state_code'])
        print('US County      =', new_dict['county_name'])
        print('Longitude      =', new_dict['longitude'])
        print('Latitude       =', new_dict['latitude'])
        print('Elevation      =', new_dict['elevation'])
        print('Last updated   =', new_dict['updated_time'])
        print('InService      =', new_dict['in_service'])
        print('PEDTS (obs)    =', new_dict['pedts_obs'])
        print('PEDTS (pred)   =', new_dict['pedts_pred'])
        print()
    
    return new_dict

#   get_station_data_from_api()
#---------------------------------------------------------------------
def create_station_info_tsv( max_lines=100 ):

    start_time = time.time()
    base_key = 'USGS_NWIS_all'  # 145375 data rows
    # base_key = 'USGS_gauged'    # 25420 data rows
    # base_key = 'GAGES2_all'     # 9322 data rows
    print('Base key =', base_key)
    print('Working...')
    usgs_id_col = cb.get_id_column( base_key )
    n_stations = 0
    n_skipped  = 0
    n_lines    = 0
    out_tsv_file = 'All_Station_Info_' + base_key + '.tsv'
    rfc_root_dir = get_rfc_data_dir( ROOT=True )
    out_tsv_path = rfc_root_dir + out_tsv_file
    out_tsv_unit = open( out_tsv_path, 'w')      #########    

    #-------------------------------------------
    # Open the TSV file with USGS IDs and info
    #-------------------------------------------
    delim = '\t'  # tab
    usgs_file = cb.get_data_filepath( base_key )
    usgs_unit = open(usgs_file, 'r')
    #-------------------        
    # Skip header line
    #-------------------
    cb.skip_header_lines( usgs_unit, key=base_key)

    #-----------------------------------------
    # Write column headings for new TSV file
    #-----------------------------------------
    headings = [
    'usgs_id', 'nws_loc_id', 'upstrm_lid', 'dnstrm_lid', 'nwm_reach_id',
    'rfc_abbrev', 'wfo_abbrev', 'usgs_name', 'description', 'state_code',
    'county_name', 'longitude', 'latitude', 'elevation', 'updated_time',
    'in_service', 'pedts_obs', 'pedts_pred' ]
    header = ''
    for h in headings:
        header += h + delim
    header = header[:-1] + '\n'
    out_tsv_unit.write( header ) 

    while (True):
        usgs_line = usgs_unit.readline()
        if (usgs_line == ''):
            break  # (reached end of file)
        n_lines += 1
#         if (n_lines > max_lines):
#             break

        values = usgs_line.split( delim )   
        usgs_id = values[ usgs_id_col ].strip()
        if (base_key == 'USGS_NWIS_all'):
            usgs_id = usgs_id.replace('USGS-','')  ######### NEED THIS!
        new_dict = get_station_data_from_api( gauge_id=usgs_id, REPORT=False)
        if (new_dict is None):
            # print('### Skipping USGS ID:', usgs_id)
            n_skipped += 1
            if ((n_skipped % 50) == 0):
                print('status: n_skipped  =', n_skipped)
            continue

        #------------------------------
        # Format some numeric strings
        #------------------------------
        lon      = new_dict['longitude']
        lat      = new_dict['latitude']
        elev      = new_dict['elevation']
        lon_str  = '{x:.5f}'.format(x=lon)
        lat_str  = '{x:.5f}'.format(x=lat)
        elev_str = '{x:.5f}'.format(x=elev)
                
        #-----------------------------------------
        # Write all station data to new TSV file
        #-----------------------------------------
        line = ''
        line += new_dict['usgs_id']         + delim
        line += new_dict['nws_loc_id']      + delim
        line += new_dict['upstrm_lid']      + delim 
        line += new_dict['dnstrm_lid']      + delim
        line += new_dict['nwm_reach_id']    + delim
        line += new_dict['rfc_abbrev']      + delim
        line += new_dict['wfo_abbrev']      + delim
        line += new_dict['usgs_name']       + delim
        line += new_dict['description']     + delim
        line += new_dict['state_code']      + delim
        line += new_dict['county_name']     + delim
        line += lon_str                     + delim
        line += lat_str                     + delim 
        line += elev_str                    + delim
        line += new_dict['updated_time']    + delim
        line += str(new_dict['in_service']) + delim
        line += new_dict['pedts_obs']       + delim 
        line += new_dict['pedts_pred']      + '\n'
        out_tsv_unit.write( line )
        n_stations += 1
        if ((n_stations % 10) == 0):
            print('status: n_stations =', n_stations)
  
    usgs_unit.close()
    out_tsv_unit.close()
    run_time = (time.time() - start_time)
    print('Run time =', run_time, '[secs].')
    print('Wrote info for', n_stations, 'stations.')
    print('Finished.')
    print()
                
#   create_station_info_tsv()
#---------------------------------------------------------------------




