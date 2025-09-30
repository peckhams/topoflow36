
#---------------------------------------------------------------------
# Copyright (c) 2025, Scott D. Peckham
#
#  Copyright (c) 2025, Scott D. Peckham
#
#  2025-06-03. Wrote add_class_nums_to_tsv_file() and
#              get_cluster_center_site_ids().
#  2025-05.    Changes to support RaFTS project.
#  2025-04.    Wrote initial version.  Extends gages2_utils.py
#              to compute indicators for all GAGES-II CONUS,
#              vs. SB3, which only has Reference basins in CONUS.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import gages2B_utils as g2B
#  >>> g2B.create_climate_data_tsv_file()
#  >>> g2B.add_class_nums_to_tsv_file()
#  >>> g2B.get_cluster_center_site_ids()
#
#---------------------------------------------------------------------
#
#  get_harbor_dir()
#  get_gages2B_data_dir()
#  get_csv_headings()
#  get_tsv_headings()
#  create_climate_data_tsv_file()
#  add_class_nums_to_tsv_file()
#  get_cluster_center_site_ids()
#  get_temperature_range()
#
#---------------------------------------------------------------------

import numpy as np
import sys

# from osgeo import ogr, osr
# import json, sys, time

from topoflow.utils.ngen import gages2_utils as g2
from topoflow.utils.ngen import usgs_utils as usgs
from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import swb_utils as swb
from topoflow.utils.ngen import classify as cfy

# from topoflow.utils.ngen import shape_utils as su

#---------------------------------------------------------------------
def get_gages2B_data_dir():

    data_dir  = dtu.get_data_dir('USGS_GAGES2_all')
    data_dir += 'basinchar_and_report_sept_2011/'
    data_dir += 'spreadsheets-in-csv-format/'
     
    return data_dir
       
#   get_gages2B_data_dir()
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
def get_tsv_headings( tsv_unit ):

    headings = get_csv_headings(tsv_unit, delim='\t')
    return headings

#   get_tsv_headings()
#---------------------------------------------------------------------
def create_climate_data_tsv_file( dataset='gages2_conus'):

    if (dataset == 'gages2_conus'):
        #----------------------------------------------
        # This GAGES-II file contains data for 9067
        # sites, apparently all in CONUS ("conterm").
        #----------------------------------------------
        # There are 9322 sites in GAGES-II.
        #   - 2057 Reference  (110 in AK/HI/PR)
        #   - 7265 Non-Reference
        #----------------------------------------------
        in_csv_file  ='conterm_climate.txt'
        in_dir       = get_gages2B_data_dir()
        in_csv_path  = in_dir + in_csv_file
        #------------------------------------------------
        out_tsv_file = 'gages2_climate_indicators.tsv'
        out_dir  = dtu.get_harbor_dir()
        out_dir += 'USGS_GAGES2/_New/All/'
        out_tsv_path = out_dir + out_tsv_file
    
    #----------------------------------------------------------
    # Read data from climate CSV file for all GAGES-II basins
    # then compute some hydroclimate indicators, and save
    # everything in a new TSV file.
    #----------------------------------------------------------
    # It looks like "conterm_climate.txt" only has data for
    # basins in CONUS.  But there is no file called:
    # AKHIPR_climate.txt in the same directory.
    #----------------------------------------------------------
    # The GAGES-II CONUS SB3 dataset has the attribute:
    #    'JUL_JAN_TMP_DIFF_DEGC'
    # but the in_csv_file here does not.  However, we have
    #    T_MAX_BASIN, T_MIN_BASIN, T_MAX_SITE, T_MIN_SITE,
    # which seems to be even better, because we don't need to
    # assume that January is always the coldest month and
    # July is always the hottest month for all locations.
    #----------------------------------------------------------
    # Note that 'PRECIP_SEAS_IND' is a seasonality indicator
    # that is defined in Dingman (2002, p. 143).
    # It seems comparable to the Knoben seasonality indicator.
    #----------------------------------------------------------
          
    #-------------------------------------
    # Open the input CSV file to read
    #--------------------------------------
    in_csv_unit = open( in_csv_path, 'r' )
    headings = get_csv_headings( in_csv_unit )
    comma    = ','
    print('in_dir =', in_dir)
    
    #------------------------------------
    # Open the output TSV file to write
    #------------------------------------
    out_tsv_unit = open( out_tsv_path, 'w')
    tab = '\t'
    
    #----------------------------------------------
    # Create header for new file with a subset of
    # original headings plus some new headings
    #----------------------------------------------
#     headings2 = ['PPT_AVG_BASIN', 'PET',
#     'SNOW_PCT_PRECIP', 'PRECIP_SEAS_IND',
#     'T_MIN_BASIN', 'T_MAX_BASIN',
#     #------------------------------------------------------
#     'JAN_PPT7100_CM', 'FEB_PPT7100_CM', 'MAR_PPT7100_CM',
#     'APR_PPT7100_CM', 'MAY_PPT7100_CM', 'JUN_PPT7100_CM',
#     'JUL_PPT7100_CM', 'AUG_PPT7100_CM', 'SEP_PPT7100_CM',
#     'OCT_PPT7100_CM', 'NOV_PPT7100_CM', 'DEC_PPT7100_CM',
#     #------------------------------------------------------------   
#     'JAN_TMP7100_DEGC', 'FEB_TMP7100_DEGC', 'MAR_TMP7100_DEGC',
#     'APR_TMP7100_DEGC', 'MAY_TMP7100_DEGC', 'JUN_TMP7100_DEGC',
#     'JUL_TMP7100_DEGC', 'AUG_TMP7100_DEGC', 'SEP_TMP7100_DEGC',
#     'OCT_TMP7100_DEGC', 'NOV_TMP7100_DEGC', 'DEC_TMP7100_DEGC']

    #-----------------------------------------------------
    # Note:  There aren't that many extras that we don't
    #        need now, so just include all of them.
    #-----------------------------------------------------
    new_headings = ['SNOW_FRACTION', 'TEMP_RANGE', 'SEASONALITY1',
       'ARIDITY', 'SWB_CLASS', 'LATITUDE', 'LONGITUDE',
       'ARIDITY2', 'ARIDITY_DIFF',
       'K_ARIDITY', 'K_SEASONALITY', 'K_SNOW_FRAC']
    headings2 = headings.copy()
    headings2 = headings2 + new_headings
    header = tab.join( headings2 ) + '\n'
    out_tsv_unit.write( header )

    problem_ids = g2.get_problem_ids()

    ####################################################################
    usgs_site_info_dict = usgs.get_usgs_site_info_dict( NWIS_ALL=True,
                               SAVE_TO_FILE=False )
    print('Working...')     # Due to last "Finished" message.
    ####################################################################
    n_sites = 0
    while (True):
        line = in_csv_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        #----------------------------------------
        # Need to remove the newline char first
        #----------------------------------------
        line = line[:-1]
        values = line.split( comma )

        #-----------------------------------------------------
        # Add a leading "0" that is missing from problem IDs
        ############# IS THIS NEEDED HERE ? ##########
        #-----------------------------------------------------
        site_id = values[0]
        if (site_id in problem_ids):
            print('### Fixed ID with leading 0 missing.')
            values[0] = '0' + values[0]
        
        #------------------------------------------------
        # Create a dictionary to map headings to values
        #------------------------------------------------
        val_dict = dict()
        for j in range(len(headings)):
            val_dict[ headings[j] ] = values[j]

        #----------------------------------------------        
        # Compute the Knoben et al. (2018) indicators
        #----------------------------------------------
        if (site_id not in usgs_site_info_dict):
            print('## site_id not in dictionary =', site_id)
            continue
        site_info = usgs_site_info_dict[ site_id ]
        longitude = np.float32(site_info['lon'])
        latitude  = np.float32(site_info['lat'])
        B_aridity, K_aridity, K_seasonality, K_snow_frac = \
            cfy.calc_knoben_indicators( latitude, val_dict )
        #####################################################
        #####################################################
                   
        #----------------------------------------
        # For each heading in all_headings,
        # write corresponding value to new file
        #----------------------------------------
        out_line = ''
        for h in headings2:
            if (h in val_dict):
                value = val_dict[h]
            elif (h == 'SNOW_FRACTION'):
                f_s = g2.get_snow_precip_fraction( val_dict )
                value = '{:.4f}'.format(f_s)
            elif (h == 'SEASONALITY1'):
                #----------------------------------------------------------
                # Note: "SB3=False" setting causes the use of T_MAX_BASIN
                #       and T_MIN_BASIN vs. JUL_JAN_TMP_DIFF_DEGC
                #----------------------------------------------------------
                f_s   = g2.get_snow_precip_fraction( val_dict )
                delta = g2.get_precip_timing_index( val_dict, f_s, SB3=False)
                value = '{:.4f}'.format(delta)
            elif (h == 'ARIDITY'):
                phi = g2.get_aridity_index( val_dict )
                value = '{:.4f}'.format(phi)
            elif (h == 'TEMP_RANGE'):
                T_range = get_temperature_range( val_dict )
                value = '{:.4f}'.format(T_range)
            elif (h == 'SWB_CLASS'):
                phi   = g2.get_aridity_index( val_dict )
                f_s   = g2.get_snow_precip_fraction( val_dict )
                delta = g2.get_precip_timing_index( val_dict, f_s, SB3=False)
                value = swb.get_swb_class( delta, f_s, phi) 
            ##############################################
            elif (h == 'LATITUDE'):
                value = '{:.4f}'.format(latitude)
            elif (h == 'LONGITUDE'):
                value = '{:.4f}'.format(longitude)
            elif (h == 'ARIDITY2'):
                value = '{:.4f}'.format(B_aridity)  # standard
            elif (h == 'ARIDITY_DIFF'):
                value = '{:.4f}'.format(phi - B_aridity)
            elif (h == 'K_ARIDITY'):
                value = '{:.4f}'.format(K_aridity)
            elif (h == 'K_SEASONALITY'):
                value = '{:.4f}'.format(K_seasonality)
            elif (h == 'K_SNOW_FRAC'):
                value = '{:.4f}'.format(K_snow_frac)
            ##############################################
            else:
                value = '-9999'
            #----------------------------
            # Add this value to outline
            #----------------------------
            out_line += (value + tab)

        #---------------------------
        # Compare annual aridities
        #---------------------------
#         a_diff = np.abs(phi - B_aridity)
#         if (a_diff > 0.05):
#             print('Annual aridity difference exceeds 0.05')
#             print('  for site ID =', site_id)
#             print('  phi from GAGES-II (PET/PPTAVG_BASIN) =', phi)
#             print('  phi from monthly PET (Hamon 1963)    =', B_aridity)
#             print('------------------------------------------------------')
            
        #------------------------------------            
        # Write out_line to output TSV file
        #------------------------------------
        out_line = out_line[:-1] + '\n'
        out_tsv_unit.write( out_line ) 
        n_sites += 1

    #------------------------------     
    # Close input and output file
    #------------------------------   
    in_csv_unit.close()  
    out_tsv_unit.close()            
    print('Wrote info for', n_sites, 'sites to new TSV file:')
    print('  ' + out_tsv_path)
    print('Finished.')
    print()

#   create_climate_data_tsv_file()
#---------------------------------------------------------------------
def add_class_nums_to_tsv_file( dataset='gages2_conus',
              indicators='knoben', n_classes_max=10):

    if (dataset == 'gages2_conus'):
        in_tsv_file='gages2_climate_indicators.tsv'
        in_dir = dtu.get_harbor_dir()
        in_dir += 'USGS_GAGES2/_New/All/'
        in_tsv_path = in_dir + in_tsv_file
        #---------------------------------------
        out_tsv_file='gages2_climate_indicators_FCM.tsv'
        out_dir = dtu.get_harbor_dir()
        out_dir += 'USGS_GAGES2/_New/All/'
        out_tsv_path = out_dir + out_tsv_file                
        
    #-------------------------------------------------------------
    # Note: In order to compute the FCM class numbers, we need
    #       need to first have computed and saved the required
    #       hydroclimate indicators, such as the Knoben
    #       indicators computed from GAGES-II data.
    #       This routine reads indicator data for all sites,
    #       then computes the FCM class number for all sites
    #       for the specified number of classes.  It adds
    #       columns to the original TSV file with class numbers.
    #-------------------------------------------------------------
    print('Computing FCM class numbers...')
#     if (indicators == 'knoben'):

    class_nums3, ctrs = cfy.get_FCM_class_nums( dataset='gages2_knoben',
                                                n_classes=3 )                                   
    n_sites = class_nums3.size 
    #-----------------------------------------------------------------                                        
    fcm_class_nums = np.zeros((n_classes_max, n_sites), dtype='uint8')
    fcm_class_nums[2,:] = class_nums3
    for k in range(3, n_classes_max):
         class_nums, ctrs = cfy.get_FCM_class_nums( dataset='gages2_knoben',
                                                    n_classes=(k+1) )                                                   
         fcm_class_nums[k,:] = class_nums
                
    #-----------------------------
    # Open the existing TSV file
    # and get the headings
    #-----------------------------
    in_tsv_unit = open( in_tsv_path, 'r')
    tab = '\t'
    headings = get_tsv_headings( in_tsv_unit )
    
    #-----------------------------------
    # Open the new, augmented TSV file
    #-----------------------------------
    out_tsv_unit = open( out_tsv_path, 'w')

    #-----------------------
    # Get the new headings
    #-----------------------
    new_headings = [
    'FCM_3CLASS_NUM', 'FCM_4CLASS_NUM', 'FCM_5CLASS_NUM',
    'FCM_6CLASS_NUM', 'FCM_7CLASS_NUM',
    'FCM_8CLASS_NUM',
    'FCM_9CLASS_NUM', 'FCM_10CLASS_NUM']
    headings2 = headings.copy()
    headings2 = headings2 + new_headings
    header = tab.join( headings2 ) + '\n'
    #-----------------------------
    out_tsv_unit.write( header )
    
    #-----------------------------------------------------     
    # For each line in original TSV file, copy all info
    # info and add some new info, then write to new file
    #-----------------------------------------------------
    print('Writing class numbers to new TSV file...')
    site_index = 0    
    while (True):
        line = in_tsv_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        #----------------------------------------
        # Need to remove the newline char first
        #----------------------------------------
        line = line[:-1]
        vals = line.split( tab )

        #------------------------------------------------
        # Create a dictionary to map headings to values
        #------------------------------------------------
#         val_dict = dict()
#         for j in range(len(headings)):
#             val_dict[ headings[j] ] = vals[j]
 
        #------------------------
        # Add values to outline
        #------------------------
        line += tab
        for k in range(2, n_classes_max):
            value = str(fcm_class_nums[k, site_index])
            line += (value + tab)
        site_index += 1

        #------------------------------------            
        # Write out_line to output TSV file
        #------------------------------------
        line = line[:-1] + '\n'
        out_tsv_unit.write( line ) 
                               
    #------------------ 
    # Close the files
    #------------------
    in_tsv_unit.close()
    out_tsv_unit.close()
    print('Finished.')
       
#   add_class_nums_to_tsv_file()
#---------------------------------------------------------------------
def get_cluster_center_site_ids( dataset='gages2_conus',
                       indicators='knoben',    ## OR 'berghuijs'
                       n_classes_max=10):

    if (dataset == 'gages2_conus'):
        in_tsv_file  = 'gages2_climate_indicators.tsv'
        out_tsv_file = 'G2_sites_closest_to_FCM_cluster_centers.tsv'
    elif (dataset == 'gages2_CHB'):
        in_tsv_file  = 'gages2_CHB_climate_indicators.tsv'
        out_tsv_file = 'CHB_sites_closest_to_FCM_cluster_centers.tsv'
 
    #----------------------       
    # Build the filepaths
    #----------------------
    in_dir = dtu.get_harbor_dir()
    in_dir += 'USGS_GAGES2/_New/All/'
    in_tsv_path = in_dir + in_tsv_file  
    #---------------------------------------------------------------
    out_dir = dtu.get_harbor_dir()
    out_dir += 'USGS_GAGES2/_New/All/'
    out_tsv_path = out_dir + out_tsv_file
                              
    #-------------------------------------------------------------
    # Note: These FCM class centers are computed using all 9050
    #       GAGES-II CONUS basins.  However, in_tsv_file can be
    #       for a subset of these, such as the Calibratable
    #       Headwater Basins (CHB).
    #-------------------------------------------------------------
    print('Computing FCM class centers...')
    centers = np.zeros((n_classes_max, n_classes_max, 3), dtype='float32')    
    for n in range(2, n_classes_max):
         n_classes = (n+1)
         class_nums, cntrs = cfy.get_FCM_class_nums( dataset='gages2_knoben',
                                                     n_classes=n_classes )
         #---------------------------------------
         # Note: cntrs has shape (n_classes, 3)
         #---------------------------------------
         for k in range(n_classes):
             centers[n,k,:] = cntrs[k,:]

    #--------------------------
    # Get the number if sites
    #--------------------------
    if (dataset == 'gages2_conus'):
        n_sites = class_nums.size
    elif (dataset == 'gages2_CHB'):
        n_sites = 1451
    
    print('Computing distances to class centers...')
    distances  = np.zeros((n_classes_max, n_classes_max, n_sites), dtype='float32')
    phi_vals   = np.zeros(n_sites, dtype='float32')
    delta_vals = np.zeros(n_sites, dtype='float32')
    fs_vals    = np.zeros(n_sites, dtype='float32')
    site_ids   = list()
    
    #-----------------------------
    # Open the existing TSV file
    # and get the headings
    #-----------------------------
    in_tsv_unit = open( in_tsv_path, 'r')
    tab = '\t'   
    headings = get_tsv_headings( in_tsv_unit )

    site_index = 0       
    while (True):
        line = in_tsv_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        #----------------------------------------
        # Need to remove the newline char first
        #----------------------------------------
        line = line[:-1]
        vals = line.split( tab )

        #------------------------------------------------
        # Create a dictionary to map headings to values
        #------------------------------------------------
        val_dict = dict()
        for j in range(len(headings)):
            val_dict[ headings[j] ] = vals[j]    

        #-----------------------------------------------    
        # Get 3 indicators for this site_id/site_index
        #-----------------------------------------------
        site_id = val_dict['STAID']
        #### if (indicators == 'knoben'):
        phi     = np.float32(val_dict['K_ARIDITY'])
        delta   = np.float32(val_dict['K_SEASONALITY'])
        fs      = np.float32(val_dict['K_SNOW_FRAC'])
        site_ids.append(site_id)
        #------------------------------
        phi_vals[site_index]   = phi
        delta_vals[site_index] = delta
        fs_vals[site_index]    = fs
        
        #------------------------------------------- 
        # Compute distances to all cluster centers
        #-------------------------------------------
#         d1   = (fs    - centers[:, :, 0])**2
#         d2   = (phi   - centers[:, :, 1])**2
#         d3   = (delta - centers[:, :, 2])**2
#         dist = np.sqrt(d1**2 + d2**2 + d3**2)
#         distances[:, :, site_index] = dist

        #------------------------------------------- 
        # Compute distances to all cluster centers
        # Order is (fs, phi, delta).
        #-------------------------------------------
        for n in range(2, n_classes_max):
            n_classes = (n+1)
            for k in range(n_classes):        
                d1   = (fs     - centers[n, k, 0])**2
                d2   = (phi    - centers[n, k, 1])**2
                d3   = (delta  - centers[n, k, 2])**2
                dist = np.sqrt(d1**2 + d2**2 + d3**2)
                distances[n, k, site_index] = dist       
        
        #---------------------------        
        # Increment the site_index
        #---------------------------
        site_index += 1

    #--------------------------------------
    # Open new TSV file with site_ids
    # that are closest to cluster centers
    #--------------------------------------
    out_tsv_unit = open( out_tsv_path, 'w')
    newline = '\n'
    
    #-------------------------------    
    # Write header to new TSV file
    #-------------------------------
    #### if (indicators == 'knoben'):
    headings = [
    'N_CLASSES', 'CLASS_NUM', 'DISTANCE', 'CLOSEST_SITE_ID',
    'K_ARIDITY',     'K_ARIDITY_CNTR',
    'K_SEASONALITY', 'K_SEASONALITY_CNTR',
    'K_SNOW_FRAC',   'K_SNOW_FRAC_CNTR']
    header = tab.join( headings ) + '\n'
    out_tsv_unit.write( header )
    
    #-----------------------------------------------
    # For each valid pair (n_classes, class_num)
    # find site_id/site_index with lowest distance            
    #-----------------------------------------------
    ## w = np.where(distances == distances
 
    for n in range(2, n_classes_max):
        n_classes = (n + 1)
        for k in range(n_classes):
            dmin = distances[n,k,:].min()
            dmin_index = np.argmin(distances[n,k,:])
            closest_site_id = site_ids[dmin_index]
            #--------------------------------------------
            line  = str(n_classes) + tab + str(k) + tab
            line += str(dmin) + tab
            line += closest_site_id + tab
            #--------------------------------------------
            # Order of centers is (fs, phi, delta).
            #----------------------------------------------          
            line += str(phi_vals[dmin_index]) + tab
            line += "{:.4}".format(centers[n,k,1]) + tab            
            #----------------------------------------------  
            line += str(delta_vals[dmin_index]) + tab
            line += "{:.4}".format(centers[n,k,2]) + tab
            #----------------------------------------------                         
            line += str(fs_vals[dmin_index]) + tab
            line += "{:.4}".format(centers[n,k,0]) + newline  #####            
            out_tsv_unit.write(line)

    #----------------------
    # Close the TSV files
    #----------------------
    in_tsv_unit.close()
    out_tsv_unit.close()
    print('Finished.')
    print()

#   get_cluster_center_site_ids()
#---------------------------------------------------------------------
def get_temperature_range( val_dict, SITE=False ):

    if not(SITE):
        T_min = val_dict[ 'T_MIN_BASIN' ]
        T_max = val_dict[ 'T_MAX_BASIN' ]
    else:
        T_min = val_dict[ 'T_MIN_SITE' ]
        T_max = val_dict[ 'T_MAX_SITE' ]

    T_range = (np.float32(T_max) - np.float32(T_min))
    return T_range  #  (float32 not string)
        
#   get_temperature_range()
#------------------------------------------------------------------------







