#
#  Copyright, 2025, Scott D. Peckham
#
#  This file is currently in the HARBOR code collection at:
#     github.com/peckhams/topoflow36/topoflow/utils/ngen.
#  Start by cloning the tf36 conda environment as "fuzzy" and then:
#     conda activate fuzzy
#     pip install scikit-fuzzy  # (not avail. for conda install)
#
#  See Berghuijs et al. (2014), Knoben et al. (2018), Mai et al. (2022).
# 
#  Original data for Knoben et al. (2018) was downloaded from:
#     https://data.bris.ac.uk/data/dataset/16ctquxqxk46h2v61gz7drcdz3
#
#  FCM = Fuzzy C-Means (classification method)
#
#  classify.py
#
#------------------------------------------------------------------------
#  To get started:
#
#  % conda activate fuzzy
#  % cd <path-to-classify.py>
#  % python
#  >>> import classify as c
#  >>> c.plot_knoben_indicators_3d()
#  >>> c.plot_gages2_sb3_indicators_3d()
#  >>> c.plot_knoben_indicators_3d()
#  >>> c.plot_FCM_clusters(dataset='sb3', n_classes=10)
#  >>> c.plot_FCM_clusters(dataset='knoben', n_classes=8)
#  >>> c.plot_two_indicators(d1='knoben', ds2='gages2', ind1='fs', ind2='fs')
#
#------------------------------------------------------------------------
#
#  get_ten_colors()
#
#  get_rafts_dir()
#  get_knoben_filepath()
#  get_gages2_sb3_indicators_filepath()
#  get_gages2_indicators_filepath()
#  get_gages2_all_filepath()
#  get_harbor_basin_filepath()
#  get_ngen_headwater_basin_filepath()    ## 2025-05
#  get_gages2_knoben_indicators_filepath()      #####
#  get_gages2_knoben_CHB_indicators_filepath()  #####
#
#  get_site_ids()
#  get_common_site_ids()
#  get_knoben_site_ids_in_gages2_sb3()
#
#  save_indicators_to_tsv()
#  read_indicators_from_tsv()
#  plot_indicators_from_tsv()
#
#  get_gages2_sb3_indicators()
#  plot_gages2_sb3_indicators_3d()
#
#  get_gages2_indicators()
#  plot_gages2_indicators_3d()
#
#  get_FCM_data()   # for SB3 or GAGES2 or KNOBEN dataset
#  get_FCM_model()
#  compare_FCM_models()
#  plot_each_FCM_cluster()
#  plot_FCM_clusters()
#  get_FCM_class_nums()
#
#  get_knoben_indicators()
#  get_gages2_knoben_indicators()  ######
#  plot_knoben_indicators_3d()
#
#  convert_knoben_aridity()
#  plot_two_indicators()
#
#  ----------------------------------------------
#  These functions allow Knoben indicators to
#  be computed from GAGES-II hydroclimate data
#  ----------------------------------------------
#  CHB = Calibratable Headwater Basins
#  This is a set of basins selected for testing
#  in NextGen.  A subset of about 300 will be
#  used for RaFTS.
#  ----------------------------------------------
#  get_day_angle()
#  get_declination()
#  get_total_PET_for_day()   ## via Hamon 1963 method
#
#  get_days_per_month()
#  get_month_list()
#  get_gages2_monthly_precip_array()  # read 12 values from GAGES-II
#  get_gages2_monthly_temp_array()    # read 12 values from GAGES-II
#  get_total_PET_for_month()          # sum daily pet for all days in month
#  get_harmon_monthly_PET_array()
#
#  get_knoben_monthly_moisture_index()
#  get_knoben_monthly_MI_array()
#  get_annual_aridity()          # standard definition, as a check
#  get_knoben_annual_aridity()
#  get_knoben_annual_seasonality()
#  get_knoben_annual_snow_fraction()
#  calc_knoben_indicators()
#
#  plot_monthly_precip_values()  # not ready
#  plot_monthly_temp_values()    # not ready
#
#  write_tsv_for_CHB_basins()
#
#------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz

from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import gages2B_utils as g2B

## from topoflow.utils.ngen import knoben

#------------------------------------------------------------------------
def get_ten_colors():

    colors = ['b', 'orange', 'g', 'r', 'c', 'm', 'y', 'k',
              'Brown', 'ForestGreen']
    return colors
    
#   get_ten_colors()
#------------------------------------------------------------------------
def get_rafts_dir():

    dir1 = '/Users/peckhams/Data/RaFTS_Work_2025/'
    return dir1

#   get_rafts_dir()
#------------------------------------------------------------------------
def get_knoben_filepath( raven=False ):

    dir1 = get_rafts_dir()
    dir1 += 'xSSA-North-America-master/data_in/'
    dir1 += 'basin_metadata/'
    if not(raven):
        file = dir1 + 'basin_climate_index_knoben_snow-knoben.txt'
    else:
        file = dir1 + 'basin_climate_index_knoben_snow-raven.txt'
    return file, dir1

#   get_knoben_filepath()
#------------------------------------------------------------------------
def get_gages2_sb3_indicators_filepath():

    dir1 = dtu.get_harbor_dir()
    dir1 += 'USGS_GAGES2/_New/SB3/'
    file = dir1 + 'SB3_all_regions_untransfBCs_sorted_SWB.tsv'
    return file, dir1

#   get_gages2_sb3_indicators__filepath()
#------------------------------------------------------------------------
def get_gages2_indicators_filepath():

    dir1 = dtu.get_harbor_dir()
    dir1 += 'USGS_GAGES2/_New/All/'
    file = dir1 + 'gages2_climate_indicators.tsv'
    return file, dir1

#   get_gages2_indicators__filepath()
#------------------------------------------------------------------------
def get_gages2_all_filepath():

    dir1 = dtu.get_harbor_dir()
    dir1 += 'USGS_GAGES2/_New/All/'
    file = dir1 + 'new_gages2_all.tsv'
    return file, dir1

#   get_gages2_all_filepath()
#------------------------------------------------------------------------
def get_noaa_rfc_filepath():

    dir1 = dtu.get_harbor_dir()
    dir1 += 'NOAA_RFCs/_New/'
    file = dir1 + 'new_combo_rfc_info.tsv'
    return file, dir1

#   get_noaa_rfc_filepath()
#------------------------------------------------------------------------
def get_harbor_basin_filepath():

    dir1 = dtu.get_harbor_dir()
    dir1 += '__Collated/'
    file = dir1 + 'collated_basins_all.tsv'
    return file, dir1

#   get_harbor_basin_filepath()
#------------------------------------------------------------------------
def get_ngen_headwater_basin_filepath():

    dir1 = get_rafts_dir()
    dir1 += 'Calibratable_Headwater_Basins/'
    file = dir1 + 'ngen_headwater_basin_gages.tsv'
    return file, dir1

#   get_ngen_headwater_basin_filepath()
#------------------------------------------------------------------------
def get_gages2_knoben_indicators_filepath():

    dir1 = dtu.get_harbor_dir()
    dir1 += 'USGS_GAGES2/_New/All/'
    file = dir1 + 'gages2_climate_indicators.tsv'
    return file, dir1    

#     dir1 = get_rafts_dir()
#     dir1 += 'Calibratable_Headwater_Basins/'
#     file = dir1 + 'gages2_climate_indicators.tsv'
#     return file, dir1
    
#   get_gages2_knoben_indicators_filepath()
#------------------------------------------------------------------------
def get_gages2_knoben_CHB_indicators_filepath():

    dir1 = dtu.get_harbor_dir()
    dir1 += 'USGS_GAGES2/_New/All/'
    file = dir1 + 'ngen_basin_climate_indicators.tsv'
    return file, dir1 
    
#   get_gages2_knoben_CHB_indicators_filepath()
#------------------------------------------------------------------------
def get_site_ids( file_path, col=0, delim='\t', numeric=False ):

    site_IDs = list()
    file_unit = open(file_path, 'r')
    line = file_unit.readline()   # skip header line

    while (True):
        #----------------------------------------
        # Read data from TSV file at filepath
        # Assume file is sorted by USGS site ID
        #----------------------------------------
        line = file_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        line = line[:-1]  # remove newline char at end
        vals = line.split(delim)
        site_ID = vals[ col ].strip()
        site_IDs.append( site_ID )

    file_unit.close()
    file_name = file_path.split('/')[-1]
    print('In the file:', file_name)
    
    #-----------------------------------------------    
    # Exclude non-numeric site IDs (e.g. Canada) ?
    #-----------------------------------------------
    if (numeric):
        site_IDs2 = list()
        for ID in site_IDs:
            if (ID.isdigit()):
                site_IDs2.append(ID)
        site_IDs = site_IDs2
        print('   Number of all-digit site_IDs =', len(site_IDs))
    else:
        print('   Number of site_IDs =', len(site_IDs))
    return sorted(site_IDs)
    
#   get_site_ids()
#------------------------------------------------------------------------
def get_common_site_ids( ds1='knoben', ds2='sb3', numeric=False,
                         raven=False, number_only=False,
                         option='common' ):

    ds1 = ds1.lower()
    ds2 = ds2.lower()
    delim1 = '\t'
    delim2 = '\t'
    id_col1 = 0
    id_col2 = 0
    if (ds1 == 'knoben'):
        file_path1, dir1 = get_knoben_filepath( raven=raven )
        delim1 = ';'
    elif (ds1 == 'sb3'):
        file_path1, dir1 = get_gages2_sb3_indicators_filepath()
    elif (ds1 == 'gages2'):
        file_path1, dir1 = get_gages2_indicators_filepath()
    elif (ds1 == 'gages2_psi'):
        file_path1, dir1 = get_gages2_indicators_filepath()
    elif (ds1 == 'gages2_all'):
        file_path1, dir1 = get_gages2_all_filepath()
    elif (ds1 == 'harbor'):
        file_path1, dir1 = get_harbor_basin_filepath()
    elif (ds1 == 'noaa_rfc'):
        file_path1, dir1 = get_noaa_rfc_filepath()
    elif (ds1 == 'ngen_gages'):
        file_path1, dir1 = get_ngen_headwater_basin_filepath()
        id_col1 = 2  ## site IDs start in 3rd column 
    else:
        print('ERROR: Unknown dataset =', ds1)                
    #--------------------------------------------------------------
    if (ds2 == 'knoben'):
        file_path2, dir2 = get_knoben_filepath( raven=raven )
        delim2 = ';'
    elif (ds2 == 'sb3'):
        file_path2, dir2 = get_gages2_sb3_indicators_filepath()
    elif (ds2 == 'gages2'):
        file_path2, dir2 = get_gages2_indicators_filepath()
    elif (ds2 == 'gages2_psi'):
        file_path2, dir2 = get_gages2_indicators_filepath()
    elif (ds2 == 'gages2_all'):
        file_path2, dir2 = get_gages2_all_filepath()
    elif (ds2 == 'harbor'):
        file_path2, dir2 = get_harbor_basin_filepath()
    elif (ds2 == 'noaa_rfc'):
        file_path2, dir2 = get_noaa_rfc_filepath()
    elif (ds2 == 'ngen_gages'):
        file_path2, dir2 = get_ngen_headwater_basin_filepath()
        id_col2 = 2  ## site IDs start in 3rd column    
    else:
        print('ERROR: Unknown dataset =', ds2)  
    #--------------------------------------------------------------
    site_IDs1  = get_site_ids( file_path1, col=id_col1,
                               delim=delim1, numeric=numeric )
    site_IDs2  = get_site_ids( file_path2, col=id_col2,
                               delim=delim2, numeric=numeric )
 
    #--------------------------                              
    # Process various options
    #--------------------------
    if (option == 'common'):                               
        common_IDs = list( set(site_IDs1) & set(site_IDs2) )
        print('Number of common site_IDs =', len(common_IDs))
        if (number_only):
            return
        return sorted(common_IDs)  # Is sorting necessary?
    elif (option == 'in1not2'):
        #--------------------------------------------------------------
        # Note: 'ngen_gages' has 2 USGS site IDs not in HARBOR,
        #       namely: '05316970', '05316992'.
        #       Both are near Leavenworth, Minnesota.
        #       Both have end dates in late 1997.
        # https://waterdata.usgs.gov/nwis/inventory/?site_no=05316970
        # https://waterdata.usgs.gov/nwis/inventory/?site_no=05316992
        #--------------------------------------------------------------
        IDs = list( set(site_IDs1) - set(site_IDs2) )
        print('Number of site_IDs in ds1 not in ds2 =', len(IDs))
        if (number_only):
            return
        return sorted(IDs)  # Is sorting necessary?
            
#   get_common_site_ids()
#------------------------------------------------------------------------
def get_knoben_site_ids_in_gages2( raven=False):

    file_path1, dir1 = get_knoben_filepath( raven=raven )
    file_path2, dir2 = get_gages2_all_filepath()
    site_IDs1  = get_site_ids( file_path1, delim=';' )
    site_IDs2  = get_site_ids( file_path2, delim='\t' )
    
    common_IDs = list( set(site_IDs1) & set(site_IDs2) )
    print('Number of Knoben site IDs in GAGES-II =', len(common_IDs))
    return sorted(common_IDs)  # Is sorting necessary?
    
#   get_knoben_site_ids_in_gages2
#------------------------------------------------------------------------
def save_indicators_to_tsv( raven=False ):

    #---------------------------------------------------------
    # Notes: Download the code and data from GitHub at:
    #   https://github.com/julemai/xSSA-North-America
    # That repo contains a folder:  data_in/basin_metadata
    # that includes data from Knoben et al. (2018) as
    # multi-column text files.
    #
    # This routine reads data for 3 climate indices from
    # the Knoben text files (aridity, seasonality, and
    # snow fraction) and also retrieves similar indices
    # from GAGES-II CONUS SB3 data set for the 5798 basins
    # in the Knoben dataset.  Some basins are in Canada
    # and have Canadian site IDs like:  04FC001.
    #---------------------------------------------------------
    # Two of the 3 climate indices, though similar in both,
    # are defined differently.
    # Aridity ranges from:
    #    -1 to 1 in Knoben et al. (2018), and from
    #     0 to 5.6 (inf) in Berghuijs et al. (2014).
    #        phi = (PET_avg / P_avg), Budyko (1974).
    # Seasonality ranges from:
    #    0 to 2 in Knoben et al. (2018), and from
    #   -1 to 1 in Berghuijs et al. (2014).
    # For both studies, snow fraction ranges from 0 to 1.
    #---------------------------------------------------------
    # Mai et al. (2022) used fuzzy c-means (FCM) clustering
    # to assign a subset of HYSETS to one of 8 classes.
    #---------------------------------------------------------
    # The 3 Berguijs (2014) indicators are in HARBOR:
    #    USGS_GAGES2/_New/SB3/
    #    SB3_all_regions_untransfBCs_sorted_SWB.tsv
    # The 3 Knoben (2018) indicators are in:
    #    xSSA-North-America-master/data_in/basin_metadata
    #    basin_climate_index_knoben_snow-knoben.txt and
    #    basin_climate_index_knoben_snow-raven.txt
    #-----------------------------------------------------------
    # B_aridity = standard, Budyko aridity (in Berghuijs)
    # K_aridity = aridity as defined in Knoben (2018)
    # Kc_aridity = K aridity converted to standard
    # B2K_aridity = standard aridity converted to Knoben
    #-----------------------------------------------------------
    # Note that K_aridity is in [-1,1], where most arid is
    # -1 and most humid is 1, so it is inversely proportional
    # to the standard, Budyko aridity, in [0, infinity].
    # Real Budyko values range are in [0, 5.6], though.
    #-----------------------------------------------------------
    
    #---------------------------------------------    
    # Open the Knoben climate indicator data set
    # and skip the one header line
    #---------------------------------------------
    file1_path, dir1 = get_knoben_filepath( raven=raven )
    file1_unit = open(file1_path, 'r')
    line1  = file1_unit.readline()  # skip header line
    delim1 = ';'

    #--------------------------------------------------------- 
    # Open the GAGES-II SB3 climate indicator data in HARBOR
    #---------------------------------------------------------
    file2_path, dir2 = get_gages2_sb3_indicators_filepath()
    file2_unit = open(file2_path, 'r')
    line2  = file2_unit.readline()  # skip header line
    delim2 = '\t'  # tab

    #----------------------------------------------
    # Open new TSV file to write indicators, etc.
    #----------------------------------------------
    new_file = dir2 + 'RaFTS_Hydro_Indicators.tsv'
    new_unit = open(new_file, 'w')

    #---------------------------------
    # Write header for new TSV file
    #---------------------------------
    headings = [
    'site_ID',
    'K_aridity', 'K_seasonality', 'K_snow_fraction',
    'B_aridity', 'B_seasonality', 'B_snow_fraction',
    'K2B_aridity', 'B2K_aridity',
    #--------------------------------------------------
    'G_ppt_avg_basin', 'G_t_avg_basin', 'G_PET',
    'G_snow_pct_precip', 'G_precip_seas_ind']
    new_header = ''
    for heading in headings:
        new_header += (heading + delim2)
    new_header = new_header[:-1] + '\n' # remove delim, add newline
    new_unit.write( new_header )
        
    #-----------------------------------------------    
    # Get common USGS site IDs in both input files
    #-----------------------------------------------
    common_site_IDs = get_common_site_ids( raven=raven )
    n_sites = len( common_site_IDs )
    n_K_bigger = 0
    n_a_diff_gt_half = 0
    #------------------------
    convert_aridity = True
    max_a_diff = 0.0
    max_a_diff_pct = 0.0
           
    for k in range(n_sites):
        site_ID = common_site_IDs[k]
        #----------------------------------------
        # Read data from the Knoben data file
        # Assume file is sorted by USGS site ID
        #----------------------------------------
        while (True):
            line1 = file1_unit.readline()
            if (line1 == ''):
                break   # (reached end of file)
            line1 = line1[:-1]  # remove newline char at end
            vals = line1.split(delim1)
            K_site_ID     = vals[0]
            K_aridity     = vals[1].strip()
            K_seasonality = vals[2].strip()
            K_snow_frac   = vals[3].strip()

            if (convert_aridity):
                #------------------------------------------------ 
                # Invert formula in Knoben (2018, p. 5093):
                #     If (phi < 1), M = (1 - phi)
                #     If (phi = 1), M = 0
                #     If (phi > 1), M = (1/phi) - 1                
                # But mean is over 12 months, so inverting this
                # formula may not be exactly right.
                #------------------------------------------------       
                K_aridity_val     = np.float32(vals[1].strip())
                K_seasonality_val = np.float32(vals[2].strip())
                if (K_aridity_val >= 0):
                    K2B_aridity_val = (1 - K_aridity_val)
                else:
                    K2B_aridity_val = (1 / (K_aridity_val + 1))
                K2B_aridity = '{:.5f}'.format(K2B_aridity_val)
                #-------------------------------------------------
                # Knoben (2018) define seasonality as the
                # range of aridity over the 12 monthly averages.
                # They say Berghuijs (2014) definition is good
                # for CONUS but perhaps not good globally.
                #-------------------------------------------------
                # This is not saved or used yet
                ##################################
                K2B_seasonality_val = (K_seasonality_val - 1)  # not this simple
                K2B_seasonality = '{:.5f}'.format(K2B_seasonality_val)
            if (K_site_ID == site_ID):
                break
                
        #-----------------------------------
        # Read data from a HARBOR TSV file
        #-----------------------------------
        while (True):         
            line2 = file2_unit.readline()
            if (line2 == ''):
                break   # (reached end of file)
            line2 = line2[:-1]  # remove newline char at end
            vals = line2.split(delim2)
            B_site_ID     = vals[0]
            B_aridity     = vals[66].strip()
            B_seasonality = vals[65].strip()
            B_snow_frac   = vals[67].strip()
            #---------------------------------------------
            # Convert standard, Budyko aridity, to the
            # normalized aridity index in Knoben paper.
            #---------------------------------------------            
            B_aridity_val = np.float32( B_aridity )
            if (B_aridity_val <= 1):
                B2K_aridity_val = (1 - B_aridity_val)
            else:
                B2K_aridity_val = (1/B_aridity_val) - 1
            B2K_aridity = '{:.5f}'.format(B2K_aridity_val)
            #---------------------------------------------
            # Next value is in: [22.26, 460.62]
            ppt_avg_basin   = vals[8].strip()
            t_avg_basin     = vals[9].strip()
            # Next value is in: [290.3, 1231.3]
            PET             = vals[16].strip()
            # Next val is in: [0, 100]; nodata=-9999.0
            snow_pct_precip = vals[17].strip()
            #-------------------------------------------------------------------
            # PRECIP_SEAS_IND = Precipitation seasonality index (Markham, 1970;
            # Dingman, 2002).  Index of how much annual precipitation falls
            # seasonally (high values) or spread out over the year (low values).
            # Based on monthly precip values from 30 year (1971-2000) PRISM
            # (PRISM Climate Group, 2008).  Range is 0 (precip spread out
            # exactly evenly in each month) to 1 (all precip falls in a
            # single month).
            #-------------------------------------------------------------------
            # Next value is in: [0,1]
            precip_seas_ind =  vals[18].strip()
            #----------------------------------------------
#             # Next value is in: [22.26, 460.62]
#             ppt_avg_basin   = np.float32( vals[8].strip())
#             t_avg_basin     = np.float32( vals[9].strip())
#             # Next value is in: [290.3, 1231.3]
#             PET             = np.float32( vals[16].strip())
#             # Next val is in: [0, 100]; nodata=-9999.0
#             snow_pct_precip = np.float32( vals[17].strip())
#             # Next value is in: [0,1]
#             precip_seas_ind = np.float32( vals[18].strip())
            #--------------------------
            if (B_site_ID == site_ID):
                break
                                  
        #--------------------------------       
        # Compare B and K aridity values
        #--------------------------------   
        a_diff = np.abs(B_aridity_val - K2B_aridity_val)
        max_a_diff = np.maximum( max_a_diff, a_diff)
        a_diff_pct = (a_diff / B_aridity_val) * 100
        max_a_diff_pct = np.maximum( max_a_diff_pct, a_diff_pct)
        if (a_diff > 0.5):
            print('For Site ID    =', site_ID)
            print('   B_aridity   =', B_aridity)
            print('   K2B_aridity =', K2B_aridity)
            print('   a_diff      =', a_diff)
            print('   a_diff_pct  =', a_diff_pct)
            print('   max_a_diff  =', max_a_diff)
            print('------------------------------------')
            n_a_diff_gt_half += 1
        if (K2B_aridity_val > B_aridity_val):
            n_K_bigger += 1
 
        #-------------------------------        
        # Write values to new TSV file
        #-------------------------------          
        vals = [site_ID,
        K_aridity, K_seasonality, K_snow_frac,
        B_aridity, B_seasonality, B_snow_frac,
        K2B_aridity, B2K_aridity,
        ppt_avg_basin, t_avg_basin, PET,
        snow_pct_precip, precip_seas_ind]
        line = ''
        for val in vals:
            line += (val + delim2)   # all strings already
        line = line[:-1] + '\n' # remove delim, add newline
        new_unit.write( line )       
                          
    #------------------   
    # Close the files
    #------------------
    file1_unit.close()
    file2_unit.close()
    new_unit.close()
                     
    #----------------------    
    # Print final message
    #----------------------
    print('Finished.')
    ## print('Number of common USGS site IDs =', n_sites)  # already printed
    print('Number of times K2B_aridity > B_aridity =', n_K_bigger)
    print('Number of times aridity diff > 0.5 =', n_a_diff_gt_half)
    print('Max aridity difference =', max_a_diff)    
    print('Max aridity diff pct   =', max_a_diff_pct)
    print()
    
#   save_indicators_to_tsv()
#------------------------------------------------------------------------
def read_indicators_from_tsv( X='K2B_aridity', Y='B_aridity'):

    n_sites = 507 ########
    x = np.zeros( n_sites, dtype='float32' )
    y = np.zeros( n_sites, dtype='float32' )

#     vals = [site_ID,
#     K_aridity, K_seasonality, K_snow_frac,
#     B_aridity, B_seasonality, B_snow_frac,
#     K2B_aridity, B2K_aridity,
#     ppt_avg_basin, t_avg_basin, PET,
#     snow_pct_precip, precip_seas_ind] 
       
    column_map = {
    'K_aridity':1, 'K_seasonality':2, 'K_snow_frac':3 ,
    'B_aridity':4, 'B_seasonality':5, 'B_snow_frac':6 ,
    #----------------------------------------------------
    'K2B_aridity':7, 'B2K_aridity':8}
    xcol = column_map[ X ]
    ycol = column_map[ Y ]
    
#     print('X =', X)
#     print('Y =', Y)
#     print('xcol =', xcol)
#     print('ycol =', ycol)

    #------------------------------------------
    # Open TSV file to read indicators, etc.
    #------------------------------------------
    dummy, tsv_dir = get_gages2_sb3_indicators_filepath()
    tsv_file = tsv_dir + 'RaFTS_Hydro_Indicators.tsv'
    tsv_unit = open(tsv_file, 'r')
    line = tsv_unit.readline()  # skip header line
    delim = '\t'  # tab
    k = 0

    while (True):         
        line = tsv_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        line = line[:-1]  # remove newline char at end
        vals = line.split(delim)
        site_ID = vals[0]
        x[k] = np.float32( vals[ xcol ].strip() )
        y[k] = np.float32( vals[ ycol ].strip() )
        k += 1        
    
    tsv_unit.close()
    return x,y
    
#   read_indicators_from_tsv()
#------------------------------------------------------------------------
def plot_indicators_from_tsv( X='B_aridity', Y='K_aridity',
                              Y2=None):

    #-----------------------------------
    # Read data for plot from TSV file
    #-------------------------------------------------------
    # Options for X and Y:
    #   K_aridity, K_seasonality, K_snow_frac, B2K_aridity
    #   B_aridity, B_seasonality, B_snow_frac, K2B_aridity
    #-------------------------------------------------------
    x1, y1 = read_indicators_from_tsv(X=X, Y=Y)
    if (Y2 is not None):
        # Y2 = 'B2K_aridity'
        x2, y2 = read_indicators_from_tsv(X=X, Y=Y2)    
        if (Y2 == 'B2K_aridity'):
            yoffset = 0.0
            ## yoffset = 0.1        
            y2 = y2 - yoffset
            
    #------------------
    # Get plot labels
    #------------------
    label_map = {
    'K_aridity':'Knoben Aridity', 'K_seasonality':'Knoben Seasonality',
    'K_snow_frac': 'Knoben Snow Fraction',
    'B_aridity':'Standard Aridity', 'B_seasonality':'GAGES-II Seasonality',
    'B_snow_frac': 'GAGES-II Snow Fraction',
    #-----------------------------------------------------------------
    'B2K_aridity':'B-to-K Aridity', 'K2B_aridity':'K-to-B Aridity'}

    #--------------------------
    # Create the X vs. Y plot
    #--------------------------
    fig, ax = plt.subplots()
    ax.plot( x1, y1, '.')
    if (Y2 is not None):
        ax.plot( x1, y2, 'r.')
    ax.set_xlabel( label_map[X] )
    ax.set_ylabel( label_map[Y] )
    title = label_map[Y] + ' vs. ' + label_map[X]
    ax.set_title( title )
    plt.show()
    
#   plot_indicators_from_tsv()
#------------------------------------------------------------------------
def get_gages2_sb3_indicators(IDs=[]):

    #--------------------------------------------------------- 
    # Open the GAGES-II SB3 climate indicator data in HARBOR
    #---------------------------------------------------------
    file_path, dir2 = get_gages2_sb3_indicators_filepath()
    file_unit = open(file_path, 'r')
    line  = file_unit.readline()  # skip header line
    delim = '\t'  # tab

    n_sites = 1947
    phi   = np.zeros( n_sites, dtype='float32'  )  # aridity
    fs    = np.zeros( n_sites, dtype='float32'  )  # snow fraction
    delta = np.zeros( n_sites, dtype='float32'  )  # seasonality
    cnum  = np.zeros( n_sites, dtype='int16'    )  # seasonality

    #-------------------------------------------------   
    # Map SWB2 classes to integers, as color indices
    #-------------------------------------------------
    cmap = {'A1':0, 'A2':1, 'A3':2, 'B1':3, 'B2':4, 'C1':5, 'C2':6,
            'D1':7, 'D2':8, 'D3':9}
                    
    #-----------------------------------
    # Read data from a HARBOR TSV file
    #-----------------------------------
    k = 0
    while (True):         
        line = file_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        line = line[:-1]  # remove newline char at end
        vals = line.split(delim)
        site_ID     = vals[0]
        seasonality = vals[65].strip()
        aridity     = vals[66].strip()
        snow_frac   = vals[67].strip()
        swb2_class  = vals[69].strip()   # string, like A1
        #----------------------------------
        if (IDs == []) or (site_ID in IDs):
            phi[k]   = np.float32( aridity )
            fs[k]    = np.float32( snow_frac )
            delta[k] = np.float32( seasonality )
            cnum[k]  = np.int16( cmap[ swb2_class ] )
            k += 1

    #-------------------------------
    # Remove unused slots in array
    #-------------------------------
    phi = phi[:k+1]
    fs  = fs[:k+1]
    delta = delta[:k+1]
    cnum  = cnum[:k+1]

    #-----------------------     
    # Close the input file
    #-----------------------
    file_unit.close()
    return fs, phi, delta, cnum
            
#   get_gages2_sb3_indicators()
#------------------------------------------------------------------------
def plot_gages2_sb3_indicators_3d():

    #---------------------------------------------------------
    # Note: cnum is a number in {0,...,9} that corresponds
    #       to one of the 10 Seasonal Water Balance classes.
    #---------------------------------------------------------
    fs, phi, delta, cnum = get_gages2_sb3_indicators()
    print('Number of GAGES-II CONUS SB3 Sites =', fs.size)

    w1 = (delta == 1.0)
    n1 = w1.sum()
    print('Number of sites with seasonality = 1.0:', n1)
    w2 = (delta <= -0.98)
    n2 = w2.sum()
    print('Number of sites with seasonality <= -9.8:', n2)
        
    colors = get_ten_colors()
    p_colors = list()
    for k in cnum:
        p_colors.append( colors[k] )
    ## p_colors = colors[ cnum ]  # doesn't work
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(fs, phi, delta, marker='.', color=p_colors)
    ax.set_xlabel('Snow fraction')
    ax.set_ylabel('Aridity')
    ax.set_zlabel('Seasonality')
    #--------------------------------
    # Reduce range for aridity axis
    #--------------------------------
    plt.ylim(0, 2.0)
    plt.show()
  
#   plot_gages2_sb3_indicators_3d()
#------------------------------------------------------------------------
def get_gages2_indicators( PRECIP_SEAS_IND=False, IDs=[] ):

    #--------------------------------------------------------- 
    # Open the GAGES-II SB3 climate indicator data in HARBOR
    #---------------------------------------------------------
    file_path, dir2 = get_gages2_indicators_filepath()
    file_unit = open(file_path, 'r')
    line  = file_unit.readline()  # skip header line
    delim = '\t'  # tab

    n_sites = 9067
    phi   = np.zeros( n_sites, dtype='float32'  )  # aridity
    fs    = np.zeros( n_sites, dtype='float32'  )  # snow fraction
    delta = np.zeros( n_sites, dtype='float32'  )  # seasonality
               
    #-----------------------------------
    # Read data from a HARBOR TSV file
    #-----------------------------------
    k = 0
    while (True):         
        line = file_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        line = line[:-1]  # remove newline char at end
        vals = line.split(delim)
        site_ID     = vals[0]
        if not(PRECIP_SEAS_IND):
            seasonality = vals[52].strip()
        else:
            seasonality = vals[25].strip()
        aridity     = vals[53].strip()
        snow_frac   = vals[50].strip()
        #----------------------------------
        if (IDs == []) or (site_ID in IDs):
            phi[k]   = np.float32( aridity )
            fs[k]    = np.float32( snow_frac )
            delta[k] = np.float32( seasonality )
            k += 1

    #-------------------------------
    # Remove unused slots in array
    #-------------------------------
    phi = phi[:k+1]
    fs  = fs[:k+1]
    delta = delta[:k+1]
    
    #-----------------------     
    # Close the input file
    #-----------------------
    file_unit.close()
    return fs, phi, delta
    
#   get_gages2_indicators()
#------------------------------------------------------------------------
def plot_gages2_indicators_3d():

    fs, phi, delta = get_gages2_indicators()
    print('Number of GAGES-II CONUS Sites (Ref + Non-Ref) =', fs.size)

    w1 = (delta == 1.0)
    n1 = w1.sum()
    print('Number of sites with seasonality = 1.0:', n1)
    w2 = (delta <= -0.98)
    n2 = w2.sum()
    print('Number of sites with seasonality <= -9.8:', n2)
        
    colors = get_ten_colors()
    p_colors = list()
    for k in cnum:
        p_colors.append( colors[k] )
    ## p_colors = colors[ cnum ]  # doesn't work
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(fs, phi, delta, marker='.', color=p_colors)
    ax.set_xlabel('Snow fraction')
    ax.set_ylabel('Aridity')
    ax.set_zlabel('Seasonality')
    #--------------------------------
    # Reduce range for aridity axis
    #--------------------------------
    plt.ylim(0, 2.0)
    plt.show()
    
#   plot_gages2_indicators_3d()
#------------------------------------------------------------------------
def get_FCM_data( dataset='sb3' ):

    #---------------
    # Get the data
    #---------------
    dataset = dataset.lower()
    if (dataset == 'sb3'):
        fs, phi, delta, cnum = get_gages2_sb3_indicators()
    elif (dataset == 'gages2'):
        fs, phi, delta = get_gages2_indicators() 
    elif (dataset == 'gages2_psi'):
        fs, phi, delta = get_gages2_indicators(PRECIP_SEAS_IND=True) 
    elif (dataset == 'knoben'):
        fs, phi, delta = get_knoben_indicators()
    elif (dataset == 'gages2_knoben'):
        fs, phi, delta = get_gages2_knoben_indicators()
    alldata = np.vstack((fs, phi, delta))
    
    return alldata
    
#   get_FCM_data()
#------------------------------------------------------------------------
def get_FCM_model( dataset='sb3', n_classes=3, show=True ):

    alldata = get_FCM_data( dataset )
    #------------------------------------------------------
    # Generate fuzzy model with ncenters cluster centers.
    # Note: Center ordering is random in this clustering
    # algorithm, so the centers may change places
    #------------------------------------------------------
    cntrs, u_orig, _, _, _, _, fpc = fuzz.cluster.cmeans(
        alldata, n_classes, 2, error=0.005, maxiter=1000)

    print('Number of classes =', n_classes)   
    print('Fuzzy partition coefficient =', fpc)
    if not(show):
        return fpc, cntrs

    #-------------------------------------
    # Show the model for given n_classes
    #-------------------------------------
    colors = get_ten_colors()
    colors = (colors + colors)  # cycle to get 20 colors
    for j in range( n_classes ):
        ## fig, ax = plt.subplots( projection='3d')
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        x  = alldata[0, u_orig.argmax(axis=0) == j]
        y  = alldata[1, u_orig.argmax(axis=0) == j]
        z  = alldata[2, u_orig.argmax(axis=0) == j]
        ax.plot(x, y, z, marker='.', linestyle='None',
                label='series ' + str(j), color=colors[j])
        ## OR use ax.scatter() here.

        ax.set_xlim(0, 1)
        ax.set_ylim(0,5.6)
        ax.set_zlim(-1,1)
        ax.set_xlabel('Snow fraction')
        ax.set_ylabel('Aridity')
        ax.set_zlabel('Seasonality')
        ax.set_title('Trained model')
        ax.legend()

    #-----------------
    # Show the plots
    #-----------------
    plt.show() 
    return fpc, cntrs
    
#   get_FCM_model()
#------------------------------------------------------------------------
def compare_FCM_models( dataset='sbe', cmin=3, cmax=13 ):

    #--------------------------------------------------------
    # Note: cmin and cmax are min & max number of classes.
    #       Each cluster represents a class.
    #--------------------------------------------------------
    dataset = dataset.lower()
    if (dataset == 'sb3'):
        print('Comparing FCM models for SB3 data...')
    elif (dataset == 'gages2'):
        print('Comparing FCM models for All GAGES2 data...')
    elif (dataset == 'gages2_psi'):
        print('Comparing FCM models for All GAGES2_PSI data...')      
    elif (dataset == 'knoben'):
        print('Comparing FCM models for Knoben data...')

    fpc_list = list()
    for nc in range(cmin, cmax+1):
        fpc, cntrs = get_FCM_model(dataset=dataset, n_classes=nc, show=False)
        fpc_list.append(fpc)
 
    fig, ax = plt.subplots()
    x = np.arange(cmax-cmin+1) + cmin
    y = fpc_list
    ax.plot(x, y, marker='.', ms=10)  # ms="marker size" in pts
    ax.set_xlabel("Number of clusters")
    ax.set_ylabel("Fuzzy partition coefficient") 
    ax.set_title('FPCs for ' + dataset + ' dataset')  
    plt.show()

#   compare_FCM_models()
#------------------------------------------------------------------------
def plot_each_FCM_cluster( dataset='sb3', n_classes=10 ):

    #------------------------------------------------------
    # Note: This plots each cluster in a separate window.
    #------------------------------------------------------
    fpc, cntrs = get_FCM_model( dataset=dataset, n_classes=n_classes, show=True)

#   plot_each_FCM_cluster()
#------------------------------------------------------------------------
def plot_FCM_clusters( dataset='sb3', n_classes=10 ):

    #---------------------------------------------------------
    # Two of the 3 climate indices, though similar in both,
    # are defined differently.
    # Aridity ranges from:
    #    -1 to 1 in Knoben et al. (2018), and from
    #     0 to 5.6 (inf) in Berghuijs et al. (2014).
    #        phi = (PET_avg / P_avg), Budyko (1974).
    # Seasonality ranges from:
    #    0 to 2 in Knoben et al. (2018), and from
    #   -1 to 1 in Berghuijs et al. (2014).
    # For both studies, snow fraction ranges from 0 to 1.
    #---------------------------------------------------------
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('Snow fraction')
    ax.set_ylabel('Aridity')
    ax.set_zlabel('Seasonality')
    ax.set_title('Trained model')
    ## ax.legend()
 
    #---------------
    # Get the data
    #---------------
    dataset = dataset.lower()
    if (dataset == 'sb3'):
        fs, phi, delta, cnum = get_gages2_sb3_indicators()
        ax.set_xlim(0, 1)
        ## ax.set_ylim(0,5.6)
        ax.set_ylim(0,3)
        ax.set_zlim(-1,1)
    elif (dataset == 'gages2'):
        fs, phi, delta = get_gages2_indicators()
        ax.set_xlim(0, 1)
        ## ax.set_ylim(0,5.6)
        ax.set_ylim(0,3)
        ax.set_zlim(-1,1)    
    elif (dataset == 'gages2_psi'):
        fs, phi, delta = get_gages2_indicators(PRECIP_SEAS_IND=True)
        ax.set_xlim(0, 1)
        ## ax.set_ylim(0,5.6)
        ax.set_ylim(0,3)
        ax.set_zlim(0,1)  ##### PRECIP_SEAS_IND is in [0,1].  
    elif (dataset == 'knoben'):
        fs, phi, delta = get_knoben_indicators()
        ## ax.set_xlim(0, 1)
        ax.set_xlim(0, 0.8)
        ax.set_ylim(-1,1)
        ax.set_zlim(0,2)
    elif (dataset == 'gages2_knoben'):
        fs, phi, delta = get_gages2_knoben_indicators()
        ## ax.set_xlim(0, 1)
        ax.set_xlim(0, 0.8)
        ax.set_ylim(-1,1)
        ax.set_zlim(0,2)
            
    alldata = np.vstack((fs, phi, delta))
    ### alldata = get_FCM_data( dataset=dataset )
    colors = get_ten_colors()
    colors = (colors + colors)

    cntrs, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(
        alldata, n_classes, 2, error=0.005, maxiter=1000, init=None)

    #-----------------------------------
    # Plot assigned clusters, for each
    # data point in training set
    #-----------------------------------
    cluster_membership = np.argmax(u, axis=0)
    for j in range(n_classes):
        x = fs[cluster_membership == j]
        y = phi[cluster_membership == j]
        z = delta[cluster_membership == j]
        ax.plot(x, y, z, marker='.', ms=4, linestyle='None',
                color=colors[j])

    #-----------------------------------------
    # Get info for the cluster centers array
    #-----------------------------------------
#     print('### cntrs =')
#     print(cntrs)
#     print('### cntrs.shape =', cntrs.shape) # (8,3)
#     print('### type(cntrs) =', type(cntrs)) # ndarray
#     print('### cntrs.size  =', cntrs.size)  # 24

    #---------------------------------------------------
    # Plot center of each fuzzy cluster as big red dot
    #---------------------------------------------------
    xc = cntrs[:,0]
    yc = cntrs[:,1]
    zc = cntrs[:,2]
    ax.plot(xc, yc, zc, marker='.', linestyle='None', ms=12,
            color='red')

    #----------------
    # Show the plot
    #----------------
    plt.show()

#   plot_FCM_clusters()
#------------------------------------------------------------------------
def get_FCM_class_nums( dataset='sb3', n_classes=10 ):

    #---------------------------------------------------------
    # Two of the 3 climate indices, though similar in both,
    # are defined differently.
    # Aridity ranges from:
    #    -1 to 1 in Knoben et al. (2018), and from
    #     0 to 5.6 (inf) in Berghuijs et al. (2014).
    #        phi = (PET_avg / P_avg), Budyko (1974).
    # Seasonality ranges from:
    #    0 to 2 in Knoben et al. (2018), and from
    #   -1 to 1 in Berghuijs et al. (2014).
    # For both studies, snow fraction ranges from 0 to 1.
    #---------------------------------------------------------

    #---------------
    # Get the data
    #---------------
    dataset = dataset.lower()
    if (dataset == 'sb3'):
        fs, phi, delta, cnum = get_gages2_sb3_indicators()
    elif (dataset == 'gages2'):
        fs, phi, delta = get_gages2_indicators()  
    elif (dataset == 'gages2_psi'):
        fs, phi, delta = get_gages2_indicators(PRECIP_SEAS_IND=True)
    elif (dataset == 'knoben'):
        fs, phi, delta = get_knoben_indicators()
    elif (dataset == 'gages2_knoben'):
        fs, phi, delta = get_gages2_knoben_indicators()
            
    alldata = np.vstack((fs, phi, delta))
    ### alldata = get_FCM_data( dataset=dataset )

    #-------------------------------------------
    # Apply FCM cluster analysis for n_classes
    #-------------------------------------------
    cntrs, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(
        alldata, n_classes, 2, error=0.005, maxiter=1000, init=None)

    #------------------------------------------
    # For each data point in alldata, get the
    # class with highest FCM membership value
    #------------------------------------------
    class_nums = np.argmax(u, axis=0)
    ## print('class_nums.dtype =', class_nums.dtype)  # 'int64'
    class_nums = np.uint8( class_nums )  #####
    
    #-----------------------------------------
    # Get info for the cluster centers array
    #-----------------------------------------
#     print('### cntrs =')
#     print(cntrs)
#     print('### cntrs.shape =', cntrs.shape) # (8,3)
#     print('### type(cntrs) =', type(cntrs)) # ndarray
#     print('### cntrs.size  =', cntrs.size)  # 24

    return class_nums, cntrs

#   get_FCM_class_nums()
#------------------------------------------------------------------------
def get_knoben_indicators( raven=False, IDs=[]):

    #-------------------------------------------------------
    # Open the Knoben et al. (2018) climate indicator data
    #------------------------------------------------------
    file_path, dir = get_knoben_filepath( raven=raven )
    file_unit = open(file_path, 'r')
    line  = file_unit.readline()  # skip header line
    delim = ';'

    n_sites = 5797
    phi   = np.zeros( n_sites, dtype='float32'  )  # aridity
    fs    = np.zeros( n_sites, dtype='float32'  )  # snow fraction
    delta = np.zeros( n_sites, dtype='float32'  )  # seasonality
                    
    #----------------------
    # Read data from file
    #----------------------
    k = 0
    while (True):         
        line = file_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        line = line[:-1]  # remove newline char at end
        vals = line.split(delim)
        site_ID     = vals[0]
        aridity     = vals[1].strip()
        seasonality = vals[2].strip()
        snow_frac   = vals[3].strip()
        #----------------------------------
        if (IDs == []) or (site_ID in IDs):
            phi[k]   = np.float32( aridity )
            fs[k]    = np.float32( snow_frac )
            delta[k] = np.float32( seasonality )
            k += 1

    #-------------------------------
    # Remove unused slots in array
    #-------------------------------
    phi = phi[:k+1]
    fs  = fs[:k+1]
    delta = delta[:k+1]

    #-----------------------     
    # Close the input file
    #-----------------------
    file_unit.close()
    return fs, phi, delta  ##, cnum
            
#   get_knoben_indicators()
#------------------------------------------------------------------------
def get_gages2_knoben_indicators(IDs=[]):

    #-----------------------------------------------------
    # Open the TSV file with Knoben indicators that were
    # computed from GAGES-II CONUS hydroclimate date
    # Now called:  'gages2_climate_indicators.tsv'
    #-----------------------------------------------------
    file_path, dir1 = get_gages2_knoben_indicators_filepath()
    file_unit = open(file_path, 'r')
    delim = '\t'

    #-------------------------------------
    # Get the headings from the TSV file
    #-------------------------------------
    headings = g2B.get_tsv_headings( file_unit )
    ### line  = file_unit.readline()  # skip header line
    
    n_sites = 9050
    phi   = np.zeros( n_sites, dtype='float32'  )  # aridity
    fs    = np.zeros( n_sites, dtype='float32'  )  # snow fraction
    delta = np.zeros( n_sites, dtype='float32'  )  # seasonality
                    
    #----------------------
    # Read data from file
    #----------------------
    k = 0
    while (True):         
        line = file_unit.readline()
        if (line == ''):
            break   # (reached end of file)
        line = line[:-1]  # remove newline char at end
        vals = line.split(delim)
        
        #------------------------------------------------
        # Create a dictionary to map headings to values
        #------------------------------------------------
        val_dict = dict()
        for j in range(len(headings)):
            val_dict[ headings[j] ] = vals[j]         
        
        site_ID       = val_dict['STAID']
        K_aridity     = val_dict['K_ARIDITY']
        K_seasonality = val_dict['K_SEASONALITY']
        K_snow_frac   = val_dict['K_SNOW_FRAC']
        #-----------------------------------------
        if (IDs == []) or (site_ID in IDs):
            phi[k]   = np.float32( K_aridity )
            fs[k]    = np.float32( K_snow_frac )
            delta[k] = np.float32( K_seasonality )
            k += 1

    #-------------------------------
    # Remove unused slots in array
    #-------------------------------
    phi = phi[:k+1]
    fs  = fs[:k+1]
    delta = delta[:k+1]

    #-----------------------     
    # Close the input file
    #-----------------------
    file_unit.close()
    return fs, phi, delta
            
#   get_gages2_knoben_indicators()
#------------------------------------------------------------------------
def plot_knoben_indicators_3d( convert_aridity=False ):

    #----------------------------------------------------------
    # Knoben hydroclimate indicators have these ranges:
    # snow fraction, fs:         range=[0,1]
    # normalized aridity, phi2:  range=[-1,1] (arid to humid)
    # seasonality, delta2        range=[0,2]
    #----------------------------------------------------------  
    fs, phi2, delta2 = get_knoben_indicators()
    n_sites = fs.size
    print('Number of Knoben Sites =', n_sites)

    #-----------------------------------------
    # Convert the aridity to standard?
    # Not perfect, due to monthly averaging.
    #-----------------------------------------
    if (convert_aridity):
        phi3 = np.zeros( n_sites, dtype='float32' )
        w1 = (phi2 >= 0)
        w2 = np.invert(w1)
        phi3[ w1 ] = (1 - phi2[w1])
        phi3[ w2 ] = (1 / (phi2[w2] + 1))
        phi2 = phi3   
    
    w1 = (delta2 < 0.1)
    n1 = w1.sum()
    print('Number of sites with seasonality < 0.4:', n1)
    w2 = (delta2 > 1.9)
    n2 = w2.sum()
    print('Number of sites with seasonality > 1.9:', n2)

    #--------------------- 
    # This doesn't work.
    #---------------------                   
#     colors = get_ten_colors()
#     p_colors = list()
#     for k in cnum:
#         p_colors.append( colors[k] )
#     p_colors = colors[ cnum ]
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(fs, phi2, delta2, marker='.')  ###, color=p_colors)
    ax.set_xlabel('Snow fraction')
    ax.set_ylabel('Aridity')
    ax.set_zlabel('Seasonality')
    #---------------------------------
    # Reduce range for aridity axis?
    #---------------------------------
    ### plt.ylim(0, 2.0)
    plt.show()
  
#   plot_knoben_indicators_3d()
#------------------------------------------------------------------------
def convert_knoben_aridity( K_aridity ):

    #--------------------------------------------
    # Let Knoben aridity = K_aridity = M. 
    # Invert formula in Knoben (2018, p. 5093):
    #     If (phi < 1), M = (1 - phi)
    #     If (phi = 1), M = 0
    #     If (phi > 1), M = (1/phi) - 1                
    # But mean is over 12 months, so this
    # inverse formula is not exactly right.
    #--------------------------------------------
    phi = np.zeros( K_aridity.size, dtype='float32')
    w1 = (K_aridity >= 0)
    w2 = np.invert( w1 )
    phi[ w1 ] = 1 - K_aridity[w1]
    phi[ w2 ] = (1 / (K_aridity[w2] + 1))
    return phi
               
#   convert_knoben_aridity()
#------------------------------------------------------------------------
def plot_two_indicators(ds1='knoben', ds2='gages2',
                        ind1='fs', ind2='fs',
                        convert_k_aridity=False):

    ds1 = ds1.lower()
    ds2 = ds2.lower()
    IDs = get_common_site_ids( ds1=ds1, ds2=ds2 )

    if (ds1 == 'knoben'):
        fs1, phi1, delta1 = get_knoben_indicators(IDs=IDs)
    elif (ds1 == 'sb3'):
        fs1, phi1, delta1, cnum = get_gages2_sb3_indicators(IDs=IDs) 
    elif (ds1 == 'gages2'):
        fs1, phi1, delta1 = get_gages2_indicators(IDs=IDs) 
    elif (ds1 == 'gages2_psi'):
        fs1, phi1, delta1 = get_gages2_indicators(IDs=IDs,
                                PRECIP_SEAS_IND=True) 
    else:
        print('ERROR: Unknown dataset =', ds1)
    #------------------------------------------------------------------           
    if (ds2 == 'knoben'):
        fs2, phi2, delta2 = get_knoben_indicators(IDs=IDs)
    elif (ds2 == 'sb3'):
        fs2, phi2, delta2, cnum = get_gages2_sb3_indicators(IDs=IDs) 
    elif (ds2 == 'gages2'):
        fs2, phi2, delta2 = get_gages2_indicators(IDs=IDs) 
    elif (ds2 == 'gages2_psi'):
        fs2, phi2, delta2 = get_gages2_indicators(IDs=IDs,
                                PRECIP_SEAS_IND=True) 
    else:
        print('ERROR: Unknown dataset =', ds2)
    #------------------------------------------------
    imap1 = {'fs':fs1, 'phi':phi1, 'delta':delta1}
    imap2 = {'fs':fs2, 'phi':phi2, 'delta':delta2}

    #------------------------------------
    # Attempt to convert Knoben aridity
    # to standard aridity?
    #------------------------------------
    if (convert_k_aridity):
        if (ds1 == 'knoben') and (ind1 == 'phi'):
            phi1 = convert_knoben_aridity( phi1 )
        if (ds2 == 'knoben') and (ind2 == 'phi'):
            phi2 = convert_knoben_aridity( phi2 )

    #--------------------------
    # Create the X vs. Y plot
    #--------------------------
    x = imap1[ ind1 ]
    y = imap2[ ind2 ]
    fig, ax = plt.subplots()
    ax.plot( x, y, '.')
#     if (Y2 is not None):
#         ax.plot( x1, y2, 'r.')
    x_label = ds1.title() + '_' + ind1
    y_label = ds2.title() + '_' + ind2
    ax.set_xlabel( x_label )
    ax.set_ylabel( y_label )
    title = (x_label + ' vs. ' + y_label)
    ax.set_title( title )
    plt.show()
   
#   plot_two_indicators()
#------------------------------------------------------------------------
#------------------------------------------------------------------------
def get_day_angle( Julian_day, DEGREES=False ):

    #---------------------------------------------------------
    # Notes:  The Julian day does not need to be an integer;
    #         decimal values can be used for more precision.
    #---------------------------------------------------------

    #-------------------------------------    
    # Use this if Julian Day starts at 1
    #-------------------------------------
    ## angle = (2 * np.pi) * (Julian_day - np.float64(1)) / np.float64(365)

    #-------------------------------------    
    # Use this if Julian Day starts at 0
    #-------------------------------------
    # Don't use Days_Per_Year() here.
    #-----------------------------------------------------------
    # We should be using 366 vs. 365 for leap years, but would
    # then need to pass year to every Day_Angle() call.
    #-----------------------------------------------------------
    angle = (2 * np.pi) * Julian_day / np.float64(365)
            
    if (DEGREES):    
        angle = angle * (np.float64(180) / np.pi)
    
    return angle
    
#   get_day_angle()
#------------------------------------------------------------------------
def get_declination( day_angle, DEGREES=False, DMS=False ):

    ########################################################
    # NB! Make sure that DEGREES and DMS default to False.
    ########################################################

    #-----------------------------------------------------------
    # Note:  The declination reaches its lowest value of -23.5
    #        degrees on the Winter Solstice (Dec. 21/22) and
    #        reaches its highest value of 23.5 degrees on the
    #        Summer Solstice (June 21/22).  It is zero for
    #        both the Vernal Equinox (Mar. 20/21) and the
    #        Autumnal Equinox (Sept. 22/23).  The value of
    #        23.4397 degrees is the fixed tilt angle of the
    #        Earth's axis from from the plane of the ecliptic.
    #-----------------------------------------------------------  
    delta = np.float64(0.006918) - \
            (np.float64(0.399912) * np.cos(day_angle)) + \
            (np.float64(0.070257) * np.sin(day_angle)) - \
            (np.float64(0.006758) * np.cos(np.float64(2) * day_angle)) + \
            (np.float64(0.000907) * np.sin(np.float64(2) * day_angle)) - \
            (np.float64(0.002697) * np.cos(np.float64(3) * day_angle)) + \
            (np.float64(0.001480) * np.sin(np.float64(3) * day_angle))
    
    #------------------------------------
    # Convert from radians to degrees ?
    #------------------------------------
    if (DEGREES):    
        delta = delta * (np.float64(180) / np.pi)
    
    #----------------------------------------
    # Convert from radians to "decimal DMS"
    #----------------------------------------
    if (DMS):    
        delta = delta * (np.float64(180) / np.pi)  # [decimal degrees]
        deg = np.int16(delta)
        min = np.int16((delta - deg) * np.float64(60))
        sec = np.int16(((delta - deg) * np.float64(60) - min) * np.float64(60))
        delta = deg + (min / np.float64(100)) + (sec / np.float64(10000))   # [decimal DMS, DD.MMSS]
    
    return delta
    
#   get_declination()
#---------------------------------------------------------------------
def get_total_PET_for_day( julian_day, T_avg, lat_deg ):

    #-----------------------------------------------------------------
    # Note: Compute the daily-average PET (potential evaporation)
    #       using the Hamon 1963 method.  
    #           PET = c * (n_daylight_hours/12) * P_sat
    # 
    #       See:
    #       https://www.hec.usace.army.mil/confluence/hmsdocs/
    #       hmstrm/evaporation-and-transpiration/hamon-method
    #
    #       e_sat = saturation vapor pressure at T_avg
    #       P_sat = saturation vapor density at T_avg
    #
    #       Source of formulas:
    #       PET (Hamon, 1963)
    #       solar_dec (solar declination)
    #       n_daylight_hours (Allen et al., 1998)
    #       sunset_hour_angle (Allen et al., 1998)
    #       e_sat (Allen et al., 1998)
    #       P_sat (Wiederhold, 1997)
    #-----------------------------------------------------------------
    # Allen, R.G., L.S. Pereira, D. Raes, M. Smith (1998)
    # Crop evapotranspiration-guidelines for computing crop water 
    # requirements - FAO irrigation and drainage
    #-----------------------------------------------------------------
    day_angle = get_day_angle(julian_day)
    solar_dec = get_declination(day_angle)
    lat_rad   = lat_deg * (np.pi / 180)
    tan_product = np.tan(lat_rad) * np.tan(solar_dec)
    sunset_hour_angle = np.arccos(-1 * tan_product)
    n_daylight_hours  = (24 / np.pi) * sunset_hour_angle
    #-----------------------------------------------------
    # This formula for e_sat is due to Brutsaert and the
    # units are kPa.  For more info, see the function:
    # Saturation_Vapor_Pressure() in solar_funcs.py
    #-----------------------------------------------------
    # Pa = Pascals = unit of pressure = N / m2.
    # N  = Newton  = unit of force    = kg * m / s2
    #-----------------------------------------------------
    # The denominator in P_sat formula is a conversion
    # from degrees Celsius to Kelvin.  So units of P_sat
    # are (kPa/K).
    #-----------------------------------------------------
    # A typical value of e_sat at 25 degC = 3.165 kPa
    # A typical value of P_sat at 20 degC = 17.2 g/m3
    
    #  kPa        Pa           kg * m
    #  ---- = ---------- =  --------------
    #   K     1000 * K      s2 * 1000 * K 
    #-----------------------------------------------------    
    e_sat = 0.6108 * np.exp( 17.27 * T_avg / (T_avg + 237.3)) # [kPa]

    #--------------------------------------------------    
    # Units of P_sat should be either g/m3 or kg/m3.
    # The factor of 216.7 must also cancel the Kelvin 
    # units in the denominator.
    #--------------------------------------------------
    P_sat = 216.7 * e_sat / (T_avg + 273.16)

    #-----------------------------------------------------
    # Hamon constant should be calibrated and validated.
    # The default value of the Hamon constant used in
    # HEC HMS (see URL above) are given as:
    #    0.0065 [in/(g/m3)]
    #    0.1651 [mm/(g/m3)] = 0.01651 [cm/(g/m3)]
    # Confirmed that:  0.0065 inches = 0.1651 mm.    
    #-----------------------------------------------------
    hamon_constant = 0.1651
    PET = hamon_constant * (n_daylight_hours / 12) * P_sat

    #-------------------------------------------------------
    # Note: Using 0.1651 gave results that are consistent
    #       with the annual total PET value from GAGES-II
    #       and the aridity values.
    #       So actual units of PET here must be cm.
    #-------------------------------------------------------
    return PET

#   get_total_PET_for_day()
#---------------------------------------------------------------------
def get_days_per_month( month_num ):

    day_map = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    day_map = np.array(day_map)
    n_days  = day_map[ month_num ]
    return n_days

#   get_days_per_month()
#---------------------------------------------------------------------
def get_month_list():

    month_list = [
    'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
    'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
    return month_list
    
#   get_month_list()
#---------------------------------------------------------------------
def get_gages2_monthly_precip_array( gages2_val_dict ):

    #-------------------------------------------
    # Note: Measurement units are centimeters.
    #       Make sure PET has the same units.
    #--------------------------------------------------------
    #       The mean annual PET given in GAGES-II has units
    #       of millimeters.  See notes in gages2_utils.py.
    #--------------------------------------------------------
    month_list = get_month_list()
    P_monthly_arr = np.zeros(12, dtype='float32')
    for k in range(12):
        P_key = month_list[k] + '_PPT7100_CM'
        P_month = np.float32(gages2_val_dict[P_key])
        P_monthly_arr[k] = P_month

    #-------------------------------------------------------------    
    # Note: In GAGES-II, PPTAVG_BASIN is the average (over many
    #       years) of the total amount of precip, in cm, for
    #       a given site.  So it is should equal the sum of
    #       all the monthly precip values, which are also
    #       averages over many years. 
    #-------------------------------------------------------------      
    TEST = False
    if (TEST):
        pptavg_basin = gages2_val_dict['PPTAVG_BASIN']
        P_sum = P_monthly_vals.sum()
        print('In GAGES-II, PPTAVG_BASIN =', pptavg_basin)
        print('Sum of monthly precip values =', P_sum)
    
    return P_monthly_arr

#   get_gages2_monthly_precip_array()
#---------------------------------------------------------------------
def get_gages2_monthly_temp_array( gages2_val_dict ):

    month_list = get_month_list()
    T_monthly_arr = np.zeros(12, dtype='float32')
    for k in range(12):
        T_key = month_list[k] + '_TMP7100_DEGC'
        T_month = np.float32(gages2_val_dict[T_key])
        T_monthly_arr[k] = T_month
    
    return T_monthly_arr
    
#   get_gages2_monthly_temp_array()
#---------------------------------------------------------------------
def get_total_PET_for_month( month_num, T_avg, lat_deg ):

    #------------------------------------------------
    # Note: Get estimate of average PET for a given
    #       month_num that ranges from 1 to 12.
    #       PET units here are centimeters.
    #------------------------------------------------
    day_map  = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    day_map  = np.array(day_map)
    start_days = np.zeros(13, dtype='int16')
    start_days[1:] = day_map.cumsum()
    
    PET_sum = np.float32(0)
    day1 = start_days[month_num - 1]
    dayz = start_days[month_num]
    for julian_day in range(day1, dayz):
        PET_sum += get_total_PET_for_day( julian_day, T_avg, lat_deg )
    PET_month = PET_sum  # [cm]

    #------------------------------    
    # This would give the average
    #------------------------------
    ### PET_month = (PET_sum / day_map[month_num-1])

    return PET_month

#   get_total_PET_for_month()
#----------------------------------------------------------------------
def get_harmon_monthly_PET_array( T_monthly_arr, lat_deg):

    PET_monthly_arr = np.zeros(12, dtype='float32')

    for month_num in range(1,13):
        T_mon   = T_monthly_arr[ month_num-1 ]
        PET_mon = get_total_PET_for_month( month_num, T_mon, lat_deg )
        PET_monthly_arr[ month_num-1] = PET_mon
        
    return PET_monthly_arr

#   get_harmon_monthly_PET_array()
 #---------------------------------------------------------------------
def get_knoben_monthly_moisture_index(P_month, PET_month):

    #------------------------------------------------------------
    # Note: MI is a version of Thornthwaite's "moisture index",
    #       computed for a given month and latitude.
    #       See Knoben et al. (2018)
    #------------------------------------------------------------
    ## print('## P_month, PET_month [cm] =', P_month, PET_month)
    ## print('## standard aridity = ', PET_month / P_month)
    if (PET_month < P_month):
        MI = 1 - (PET_month / P_month)
    elif (PET_month == P_month):
        MI = 0
    else:
        MI = (P_month / PET_month) - 1
    return MI

#   get_knoben_monthly_moisture_index()
#---------------------------------------------------------------------
def get_knoben_monthly_MI_array( P_monthly_arr, PET_monthly_arr ):

    #---------------------------------------------    
    # Compute and save average MI for each month
    #---------------------------------------------
    MI_monthly_arr = np.zeros(12, dtype='float32')
    for k in range(12):
        P_month   = P_monthly_arr[k]
        PET_month = PET_monthly_arr[k]
        #-------------------------------
        MI = get_knoben_monthly_moisture_index(P_month, PET_month)
        MI_monthly_arr[k] = MI
    
    return MI_monthly_arr

#   get_knoben_MI_monthly_array()
#---------------------------------------------------------------------
def get_annual_aridity(P_monthly_arr, PET_monthly_arr):

    #----------------------------------------------------------
    # Note: Both precip and PET must have same units, which
    #       will typically be mm or cm.  In GAGES-II, precip
    #       has units of cm and the annual PET it provides
    #       has units of mm.  Values are usually in [0, 5.7].
    #---------------------------------------------------------- 
    P_sum   = P_monthly_arr.sum()    # [cm]
    PET_sum = PET_monthly_arr.sum()  # [cm]
    aridity = (PET_sum / P_sum)
    return aridity

#   get_annual_aridity()
#---------------------------------------------------------------------
def get_knoben_annual_aridity( MI_monthly_arr ):

    #------------------------------------------------------------
    # Note: I_m is the annual average aridity as defined by
    #       Knoben et al. (2018), which is computed from a
    #       version of Thornthwaite's monthly "moisture index".
    #       See Knoben et al. (2018)
    #------------------------------------------------------------
    MI_sum = np.sum( MI_monthly_arr )
    I_m = (MI_sum / 12)
    
    return I_m

#   get_knoben_annual_aridity()
#---------------------------------------------------------------------
def get_knoben_annual_seasonality( MI_monthly_arr ):

    #------------------------------------------------------------
    # Note: I_mr is a measure of the "seasonality of aridity",
    #       as defined by Knoben et al. (2018). It is simply
    #       the range of the monthly aridity values.
    #       See Knoben et al. (2018)
    #------------------------------------------------------------
    I_mr = MI_monthly_arr.max() - MI_monthly_arr.min()
    return I_mr

#   get_knoben_annual_seasonality()
#---------------------------------------------------------------------
def get_knoben_annual_snow_fraction( P_monthly_arr, T_monthly_arr ):

    T0 = 0  # (rain-to-snow temp threshold from Knoben 2018)
    w = np.where(T_monthly_arr < T0)
    fs = (P_monthly_arr[w].sum() / P_monthly_arr.sum())
    return fs

#   get_knoben_annual_snow_fraction()
#---------------------------------------------------------------------
def calc_knoben_indicators( lat_deg, gages2_val_dict ):

    #---------------------------------------------
    # Note:  lat = latitude for the USGS site ID
    #---------------------------------------------
    P_monthly_arr = get_gages2_monthly_precip_array( gages2_val_dict)
    T_monthly_arr = get_gages2_monthly_temp_array( gages2_val_dict)

    PET_monthly_arr = get_harmon_monthly_PET_array( T_monthly_arr, lat_deg )
    #-----------------------------------------------------
    # Compute standard aridity as "ARIDITY2", as a check
    #-----------------------------------------------------
    aridity = get_annual_aridity( P_monthly_arr, PET_monthly_arr)
       
    MI_monthly_arr = get_knoben_monthly_MI_array( P_monthly_arr, PET_monthly_arr )
    ### get_monthly_MI_array( gages2_val_dict )

    K_aridity = get_knoben_annual_aridity( MI_monthly_arr )
    K_seasonality = get_knoben_annual_seasonality( MI_monthly_arr )
    K_snow_frac = get_knoben_annual_snow_fraction( P_monthly_arr, T_monthly_arr )
    
    return aridity, K_aridity, K_seasonality, K_snow_frac
   
#   get_knoben_indicators()
#---------------------------------------------------------------------
def plot_monthly_precip_values( gages2_val_dict, site_id ):

    # Plot and fit a sine curve.
    # P_monthly_arr = get_gages2_monthly_precip_array( gages2_val_dict)
    pass

#   plot_monthly_precip_values()
#---------------------------------------------------------------------
def plot_monthly_temp_values( gages2_val_dict, site_id ):

    # Plot and fit a sine curve.
    # T_monthly_arr = get_gages2_monthly_temp_array( gages2_val_dict)
    pass

#   plot_monthly_temp_values()
#---------------------------------------------------------------------
def write_tsv_for_CHB_basins(FCM=False):

    #--------------------------------------------
    # Note: CHB = Calibratable Headwater Basins
    #--------------------------------------------
    common_site_ids = get_common_site_ids( ds1='gages2_all',
                          ds2='ngen_gages', numeric=True,
                          option='common' )

    #----------------------------------
    # Open the input TSV file to read
    #----------------------------------
    if (FCM):
        in_tsv_file='gages2_climate_indicators_FCM.tsv'
    else:
        in_tsv_file='gages2_climate_indicators.tsv'
    in_dir = dtu.get_harbor_dir()
    in_dir += 'USGS_GAGES2/_New/All/'
    in_tsv_path = in_dir + in_tsv_file
    in_tsv_unit = open( in_tsv_path, 'r')
    tab = '\t'

    #------------------------------------
    # Open the output TSV file to write
    #------------------------------------
    if (FCM):
        out_tsv_file='gages2_CHB_climate_indicators_FCM.tsv'    
    else:
        out_tsv_file='gages2_CHB_climate_indicators.tsv'
    out_dir = dtu.get_harbor_dir()
    out_dir += 'USGS_GAGES2/_New/All/'
    out_tsv_path = out_dir + out_tsv_file
    out_tsv_unit = open( out_tsv_path, 'w')
    tab = '\t'

    #---------------------------------    
    # Read and write the header line
    #---------------------------------
    line = in_tsv_unit.readline()
    out_tsv_unit.write(line)
    
    print('Working...')   
    while (True):
        line = in_tsv_unit.readline()
        if (line == ''):
            break
        words = line.split(tab)
        site_id = words[0].strip()
        if (site_id in common_site_ids):
            out_tsv_unit.write(line)  # line should end with newline

    in_tsv_unit.close()
    out_tsv_unit.close()
    
    print('Finished creating TSV file with climate indicators')
    print('for ngen calibratable headwater basins in GAGES-II.')
    print()
     
#   write_tsv_for_CHB_basins()
#---------------------------------------------------------------------

## site_id not in dictionary = 01372058
## site_id not in dictionary = 01484085
## site_id not in dictionary = 02300021
## site_id not in dictionary = 02300042
## site_id not in dictionary = 02310747
## site_id not in dictionary = 02322698
## site_id not in dictionary = 02326550
## site_id not in dictionary = 02343801
## site_id not in dictionary = 02403500
## site_id not in dictionary = 07250085
## site_id not in dictionary = 08041749
## site_id not in dictionary = 08041780
## site_id not in dictionary = 08074000
## site_id not in dictionary = 08075500
## site_id not in dictionary = 10172860
## site_id not in dictionary = 10309101
## site_id not in dictionary = 14211720











