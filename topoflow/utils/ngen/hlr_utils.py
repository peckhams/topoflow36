#
# Copyright (c) 2023-2024, Scott D. Peckham
#
# Jan 2024. Renamed all "ngen/utils" files to end in "_utils.py"
#           instead of "_tools.py"
#           Modified to use new data_utils.py.
# Nov 2023. Wrote plot_color_legend().
# Oct 2023. Wrote get_hlr_info_dict(), create_hlr_grid(),
#           get_hlr_grid_info(), get_hlr_code_for_point(), &
#           Lambert_Azimuthal_Equal_Area_XY().
#           Wrote get_dem_info() as separate function.
#           Wrote try_to_get_outlet_info() as separate function
#           that calls several others to clean up code.
#           Modified create_tsv() to include the closest USGS
#           site ID and distance, by calling:
#           usgs.get_closest_site().
#           Changed string formatting for distance: 2-digit prec.
# Jul 2023. Modified to use newer shape_utils.py, which
#           contains many of the original functions from here.
#           Added get_hlr_data_dir, get_hlr_dem_dir.
#           Added get_repeated_value_list, to exclude them.
# Jun 2023. Wrote convert_vat_adf_to_csv, get_bounding_box.
#           Moved many functions into shape_utils.py.
# May 2023. Wrote: create_tsv, get_polygon_points (all versions),
#           get_cols_and_rows, get_polygon_zvals,
#           check_lon_lats, get_zmin_lon_and_lat, and
#           convert_coords.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import hlr_utils as hlr
#  >>> hlr.create_tsv(nb_max=100, tsv_file='new_hlr_na_100.tsv')
#  >>> hlr.create_tsv(nb_max=50000, tsv_file='new_hlr_na_all2.tsv',
#          SAVE_OUTLET_INFO=False, REPORT=False)
#  >>> hlr.create_tsv(nb_max=50000, tsv_file='new_hlr_na_all.tsv',
#          SAVE_OUTLET_INFO=True, REPORT=False)
#  Last line results in 9539 basins that have a single, unique
#  zmin (outlet elevation) value found from shapefile.
#
#---------------------------------------------------------------------
# Notes:
#
# The HLR (Hydrologic Landscape Regions) data set from USGS has
# some metadata, including HLR codes in {0,...,20}, for what is
# claimed to be 43931 basins.
# See:  https://www.sciencebase.gov/catalog/item/
# 4f4e4786e4b07f02db485758
# 
# However, the HLR shapefile in the "hlrshape" folder of this
# dataset may be a draft version because it has several problems,
# such as:
#
#    * 47479 features/basins (3548 more than 43931)
#
#    * 659 basins with invalid HLR code of 0.
#      These have COUNT=0, VALUE=-9999
#
#    * Multiple features/basins with the same VALUE field, which
#      is supposed to be a unique watershed ID.
#
#    * 3391 basins w/ area < 1e7 sqm = 10 km2, even though it is
#      claimed that the pruning threshold is 100 km2.
#
#    * Strange coordinates for the points on the basin boundaries.
#
#  To save the shapefile attribute table to a CSV file in QGIS:
#    - Right click on the layer in the Layers panel
#    - Choose Export > Save Features As...
#    - Choose Comma Separated Values (.csv) from Format droplist

#  There is a "raster" version of the HLR dataset that consists of
#  a set of ESRI ADF files in the arctar00000/hlrus folder.  This
#  dataset has fewer problems (e.g. 43931 basins) and seems more
#  consistent with a "final" version than the shapefile/vector
#  dataset.  While attributes for the HLR shapefile are included
#  in the shapefile attribute table, the attributes of the HLR
#  raster dataset are in the file: vat.adf.  The function called:
#  convert_vat_adf_to_csv() reads attributes from the vat.adf file
#  and writes them to a CSV file: vat_adf.csv.  The resulting CSV
#  file has 43931 rows, no HLR codes of 0, and (it appears) no
#  repeated values in the VALUE (watershed ID) column.  It seemed
#  to work but got this error message at the end:
#  ERROR 5: OSRCalcInvFlattening(): Wrong input values
##########################################################

#  A new HLR shapefile (hlrus4.shp) was created by reading the
#  raster version of the HLR data set into QGIS and then using
#  its raster to vector conversion tool as follows:
#    Raster > Conversion > Polygonize (Raster to Vector).
#    Change DN to VALUE.
#    The resulting shapefile attribute table has only one column:
#    VALUE (watershed ID). 
#  The new shapefile was saved into the HLR shapefiles folder. 
#  To add the (raster) attribute table (vat_adf.csv; see above)
#  to this new shapefile, QGIS was again used and the steps are 
#  described in the text file:
#  "__HLR_shapefile_notes.txt" in the hlrshape folder.

###############################################################
#  NOTE:  When the attribute table of this new shapefile was
#  saved to CSV (see above), it had many of the same problems
#  as the original shapefile.
###############################################################

#  Note that neither the vector or raster version of the HLR data
#  attribute tables provide names for the basins.  However, it is
#  likely that many of them are also USGS gauged basins, so it
#  may be possible to get names and USGS basin IDs from outlet
#  lons and lats.  Wolock is retiring soon but plans to generate
#  an updated version of the HLR data set in the near future.

#  Note that neither the vector or raster version of the HLR data
#  attribute tables provides the coordinates of the basin outlet
#  nor the geographic bounding box coordinates.  These are needed
#  for modeling these watersheds in NextGen.  The set of utilities
#  in this file attempt to rectify these omissions.  First, the
#  coordinates of each watershed polygon are read from the (new)
#  shapefile so that the min & max of lon & lat can be found for
#  each basin.  These bounds are padded slightly larger.  Second,
#  the HYDRO1K DEM for North America --- which was used to create
#  the HLR data set in the first place --- is used in an effort to
#  determine the outlet coordinates (lon/lat) for each basin.  The
#  elevation values for every grid cell on the basin boundary are
#  found, and the outlet is assumed to be the grid cell with the
#  lowest elevation.  However, this elevation is not always unique.

#  The source DEM is the GTOPO30 HYDRO1K DEM and uses a Lambert
#  Azimuthal Equal Area projection as described in its PRJ file.
        
#  In all of these utilities, we use GDAL and the PRJ files to
#  convert coordinates (e.g. Lambert Azimuthal Equal Area) to WGS84
#  lon/lat.  See the "convert_coords()" function in shape_utils.py.

#  Many code comments are for original (old) shapefile, that exhibit
#  that exhibits the numerous issues listed above.  The "new"
#  shapefile (hlrus4.shp) comes from QGIS work described above.
#
#---------------------------------------------------------------------
#
#  get_hlr_dem_dir()
#  get_repeated_values_list()
#  get_missing_basin_name()   # NOT READY YET
#  get_dem_info()
#
#  create_tsv()      # main function
#
#  basin_in_dem()
#  try_to_get_outlet_info()
#  get_cols_and_rows()
#  get_polygon_zvals()
#  check_lons_lats()
#  get_zmin_lon_and_lat()
#  get_zmin_lon_and_lat1()  # uses np.argmin()
#
#  convert_vat_adf_to_csv()
#  get_hlr_info_dict()
#  create_hlr_grid()
#  get_hlr_grid_info()
#  Lambert_Azimuthal_Equal_Area_XY()
#  get_hlr_code_for_point()
#
#  point_in_rectangle()
#  save_hlr_code_grid_to_rtg()
#  create_hlr_conus_grids
#
#  get_hlr_outlet_info()
#  get_closest_hlr_outlet()
#
#  plot_color_legend()
#
#---------------------------------------------------------------------
import numpy as np
import os, os.path, pickle  # to save dictionary

from osgeo import ogr, osr
import json, time

from matplotlib import pyplot as plt   # for plot_shapefile()
from matplotlib import patches

from osgeo import gdal
from topoflow.utils.ngen import shape_utils as su
from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import usgs_utils as usgs

from topoflow.utils import regrid     # to create grid of HLR codes
from topoflow.utils import rtg_files  # to save HLR grid to RTG
from topoflow.utils import rti_files  # to create RTI file for RTG

# Note: Lambert_Azimuthal_Equal_Area_XY is in both this file
#       and in topoflow/utils/projections.py 
## from topoflow.utils import projections as proj  ##########

# For convert_vat_adf_to_csv()
import csv, sys, gdal

#---------------------------------------------------------------------
def get_hlr_dem_dir( region='North_America' ):
   
    if (region == 'North_America'):
        repo_dir = dtu.get_repo_dir()
        dem_dir = repo_dir + 'USGS_HLR/Data/HYDRO1K_DEM_NA/'
    elif (region == 'CONUS'):
        data_dir = dtu.get_new_data_dir( 'USGS_HLR' )
        dem_dir  = data_dir + 'Raster/' 
    return dem_dir
    
#   get_hlr_dem_dir
#---------------------------------------------------------------------
def get_repeated_values_list():

    #------------------------------------------------------------
    # Note: The following watershed IDs (VALUE column) are
    #       repeated and will be removed or handled separately.
    #       There are only 103 of them out of 43931.
    #       They seem to be fragmented polygons.
    #------------------------------------------------------------
    # Note: The raster version of the HLR dataset is in the
    #       folder "arctar00000/hlrus" and it contains a "value
    #       attribute table" (vat.adf) which does not have
    #       repeated values. A function convert_vat_adf_to_csv
    #       herein converts this table to CSV: "vat_adf.csv".
    #------------------------------------------------------------    
    vals_list = [   
    25,    26,    30,    52,    64,    85,    111,   148,   217,  
    251,   311,   312,   334,   337,   375,   415,   433,   518,
    546,   880,   1009,  1452,  1516,  1775,  3826,  3972,  4226, 
    5020,  5139,  5206,  5269,  5387,  5439,  5559,  5580,  5699,
    5799,  5844,  5865,  5909,  5937,  5950,  5955,  6024,  6058,
    6062,  6106,  6141,  6202,  6252,  6259,  6348,  6434,  6450, 
    6482,  6551,  6552,  6592,  6604,  6619,  6628,  6638,  6645,
    6661,  7286,  7509,  8008,  8046,  8571,  9209,  10119, 11297,
    11563, 12197, 13489, 18233, 20580, 22586, 23175, 24662, 24971,
    27963, 28267, 28525, 30702, 31507, 33222, 33421, 34024, 35179,
    35309, 36331, 37980, 38642, 39686, 40056, 41545, 42160, 42783,
    43533, 43849, 43922, 44217 ]

    #--------------------------------- 
    # Corresponding histogram values
    #---------------------------------
#     hist = [
#     4 3 2 7 4 3 2 2 2 2 2 3 5 3 2 2 2 5 2 2 2 2 3 2 2 2 3 2 4
#     4 2 2 3 2 4 2 2 4 2 2 2 4 2 2 2 2 2 7 4 2 3 2 3 4 3 2 2 3
#     4 2 2 4 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#     2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3]
    return vals_list
    
#   get_repeated_values_list()
#---------------------------------------------------------------------
def get_missing_basin_name( lon, lat ):

    #------------------------------------------------------------
    # Note: The HLR data set does not provide basin names.
    #       Try to match lon/lat to a USGS site to get name.
    #------------------------------------------------------------
    pass
    
#   get_missing_basin_name()
#---------------------------------------------------------------------
def get_dem_info( dem_path ):

    #-------------------------------------------
    # HYDRO-1K DEM Info (Lambert Azimuthal EA)
    #-----------------------------------------------
    # The HYDRO-1K DEM used for the original HLR
    # work appears to have: ncols=9106, nrows=6855
    # based on the HLR metadata page.
    #-----------------------------------------------
    dem_info = dict()
    dem_info['path']  = dem_path
    dem_info['ncols'] = 9102
    dem_info['nrows'] = 8384
    dem_info['dx']    = 1000.0  # meters
    dem_info['dy']    = 1000.0  # meters
    #-------------------------------------------
    # Use half grid cell offsets: 500 = 1000/2
    #-------------------------------------------
    dem_info['xmin']  = -4462500.0
    dem_info['xmax']  =  4639500.0
    dem_info['ymin']  = -3999500.0
    dem_info['ymax']  =  4384500.0
    dem_info['zmin']  = -207
    dem_info['zmax']  = 6106  # Denali, AK
    dem_info['nodata'] = -9999  # used for ocean grid cells
    return dem_info

    #----------------------
    # From blw world file
    #----------------------
#     xmin  = -4462000.0
#     xmax  = xmin + (ncols * dx)
#     ymax  = 4384000.0
#     ymin  = ymax - (nrows * dy)
    
#   get_dem_info()
#---------------------------------------------------------------------
def create_tsv( shape_dir=None, dem_dir=None, new_dir=None,
                ## shape_file='hlrus4.shp',
                shape_file='hlrus.shp', prj_file='hlrus.prj',                     
                dem_file='na_dem.bil',
                old_shapefile = False,  # old = original shapefile
                tsv_file='new_hlr_na_with_outlets.tsv',  # only North America
                SAVE_OUTLET_INFO=True,
                nb_max=50, SWAP_XY=True, REPORT=True ):

    #---------------------------------------------------------
    # Note: To use original shapefile "as is", set:
    #       shape_file = 'hlrus.shp' and old_shapefile=True.
    #---------------------------------------------------------
    # Note: The projection details in hlrus.prj match those
    #       of the North America HYDRO1K DEM in na_dem.prj.
    #---------------------------------------------------------
    start_time = time.time() 
    new_dir = dtu.get_new_data_dir( 'USGS_HLR' )
    if (shape_dir is None):
        shape_dir = dtu.get_data_dir( 'USGS_HLR' )
        shape_dir += 'hlrshape/'
    if (dem_dir is None):
        dem_dir = get_hlr_dem_dir()

    ## prj_file = shape_file.replace('.shp', '.prj')
    
    shape_path = shape_dir + shape_file
    prj_path   = shape_dir + prj_file
    dem_path   = dem_dir   + dem_file
    tsv_path   = new_dir   + tsv_file
    shape_unit = su.open_shapefile( shape_path )
    layer      = shape_unit.GetLayer()

    dem_info = get_dem_info( dem_path )
    dem_unit = open( dem_path, 'rb')
    tsv_unit = open( tsv_path, 'w') 

    counters = dict()
    counters['total_features']    = np.int32(0)
    counters['unique_basins']     = np.int32(0)
    counters['excluded_features'] = np.int32(0)
    counters['saved_basins']      = np.int32(0)
    #------------------------------------------
    counters['one_cell_basin']  = np.int32(0)
    counters['code_is_zero']    = np.int32(0)
    #------------------------------------------
    if (SAVE_OUTLET_INFO):   
        counters['hlr_zmin_match'] = np.int32(0)   # hlr_zmin = zmin
        counters['one_zmin']       = np.int32(0)
        counters['one_zmin_match'] = np.int32(0)
        counters['two_tries']      = np.int32(0)
        counters['not_in_dem']     = np.int32(0)   # dem is only N. America
        counters['no_zvalues']     = np.int32(0)   # can't happen now?
        not_in_dem_vals = list()
    count_area_count = np.int32(0)  # for old_shapefile only

    ## This is outdated, or maybe applies only to old_shapefile.
    ## repeated_val_list = get_repeated_values_list()
    
    max_val = 48000 
    hlr_histogram = np.zeros(21, dtype='int32')
    val_histogram = np.zeros(max_val, dtype='int32')
    insert_key    = 'VALUE'  # insert new cols after this one
    hlr_zmin_match = np.zeros(max_val, dtype='bool')
    
    #-------------------------------------------
    # Get column headings for the new TSV file
    #-------------------------------------------
    if (SAVE_OUTLET_INFO):
        new_headings  = ['Closest_USGS_ID', 'dist_to_closest']
        new_headings += ['OUTLON', 'OUTLAT', 'OUTCOL', 'OUTROW']
        new_headings += ['MINLON', 'MAXLON', 'MINLAT', 'MAXLAT']
    else:
        new_headings = ['MINLON', 'MAXLON', 'MINLAT', 'MAXLAT']

    #---------------------------------------------
    # Get coords of all USGS gauged basins so we
    # can find the closest site and distance
    #---------------------------------------------
    if (SAVE_OUTLET_INFO):
        # NWIS_ALL = False  # 27,915 USGS stream gauges (NWIS_Web)
        # NWIS_ALL = False  # 25,420 USGS stream gauges (NWIS_Web_Old)
        NWIS_ALL = True   # 145,375 USGS stream gauges  (NWIS_WQP)
        usgs_site_coords = usgs.get_usgs_site_coords( NWIS_ALL=NWIS_ALL )
    
    #-----------------------------------------    
    # Iterate over the features in the layer
    # GeometryType = 3 is POLYGON.
    #---------------------------------------------------------
    # Names in the shapefile attribute table are:
    # AREA, PERIMETER, HLRUS_COV_, HLRUS_COV1, VALUE, COUNT,
    # AQPERMNEW, SLOPE, TAVE, PPT, PET, SAND, PMPE, MINELE,
    # RELIEF, PFLATTOT, PFLATLOW, PFLATUP, HLR
    #---------------------------------------------------------
    print('Working...')
    for feature in layer:
        geometry   = feature.GetGeometryRef()
        attributes = feature.items()  # This is a Python dictionary
        # n_rings    = geometry.GetGeometryCount()

        #--------------------------------
        # Write header for new TSV file
        #--------------------------------
        if (counters['total_features'] == 0):
            su.write_tsv_header( tsv_unit, attributes,
                insert_key=insert_key, new_headings=new_headings )
              
        #-----------------------------
        # Get some of the attributes
        #--------------------------------------------
        # In some cases, AREA = COUNT * 1e6, but
        # this is not consistent. For basin #6,
        # (just 1 cell), area=1e6 but count = 362.
        # In most, but not all cases,
        # PERIMETER = n_pts * 1000.  For basin #7,
        # n_pts = 124 but PERIMETER = 144000.0.
        # For basin #28 (2 cells), AREA=2e6,
        # COUNT=377, PERIMETER=6e3
        #--------------------------------------------
        hlr_value  = attributes['VALUE']  # watershed ID  
        hlr_count  = attributes['COUNT'] 
        hlr_zmin   = attributes['MINELE']
        hlr_relief = attributes['RELIEF']
        hlr_slope  = attributes['SLOPE']
        hlr_code   = attributes['HLR']
        #-----------------------------------
        EXCLUDE = False
        if (hlr_code == 0):
            counters['code_is_zero'] += 1
            EXCLUDE = True  # skip this basin
        #-----------------------------------            
        if (old_shapefile):
            hlr_area   = attributes['AREA']
            hlr_perim  = attributes['PERIMETER']
            if (hlr_area == 1e6 * hlr_count):
                count_area_count += 1
            if (hlr_area == 1e6):
                counters['one_cell_basin'] += 1
                EXCLUDE = True  ## skip this basin
        else:
            if (hlr_count == 1):
                #----------------------------------
                # In the new shapefile, hlr_count
                # has (min,max) = (31, 23476)
                #----------------------------------
                counters['one_cell_basin'] += 1
                EXCLUDE = True ## skip this basin

        #----------------------------------------
        # Print attributes in attribute table ?
        #----------------------------------------
        if (REPORT):
            n_features = counters['total_features']
            su.print_attributes( attributes, n_features )
 
        #----------------------------------------------
        # Skip all watersheds with "repeated VALUE"
        # that seem to have fragmented basin polygons
        #----------------------------------------------
        # Could return bounding box of the fragments
        # but then which outlet would we choose?
        #----------------------------------------------
        counters['total_features'] += 1   # even if repeated VALUE
        if (val_histogram[hlr_value] == 0):
            counters['unique_basins'] += 1
        val_histogram[ hlr_value ] += 1
#         if (val_histogram[ hlr_value ] > 1):
#             print('Skipping repeated HLR basin ID =', hlr_value)
#             EXCLUDE = True
        #----------------------------------------------
        # This list only has 103 repeated values and
        # doesn't seem to include all repeated values.
        # See get_repeated_values_list().
        #----------------------------------------------                
#         if (int(hlr_value) in repeated_val_list):
#             print('Skipping repeated HLR basin ID =', hlr_value)
#             EXCLUDE = True

        #---------------------------------------------
        # If EXCLUDE has been set, skip this feature
        #---------------------------------------------
        if (EXCLUDE):
            counters['excluded_features'] += 1
            continue
                    
        #------------------------------------------- 
        # Get all points/vertices in this geometry
        #-------------------------------------------
        poly_x, poly_y = su.get_polygon_points( feature )

        #------------------------------------ 
        # Exclude the last, repeated vertex
        #------------------------------------
        poly_x = poly_x[:-1]
        poly_y = poly_y[:-1]

        #---------------------------------------
        # Snap polygon coordinates to DEM grid
        #---------------------------------------
        poly_x, poly_y = snap_polygon_points( poly_x, poly_y, hlr_value, old_shapefile )
        # print('poly_x =', poly_x)
        # print('poly_y =', poly_y)
        # break

        #--------------------------------------------------        
        # Convert Lambert coords to Geo. lon/lat (WGS-84)
        #--------------------------------------------------
        poly_lons, poly_lats = su.convert_coords(poly_x, poly_y,
                                  inPRJfile=prj_path,
                                  SWAP_XY=SWAP_XY, PRINT=False)
#         print('poly_lons =', poly_lons)
#         print('poly_lats =', poly_lats)
#         break

        #------------------------------ 
        # Get geographic bounding box
        #------------------------------
        minlon, maxlon, minlat, maxlat = \
                su.get_bounding_box(poly_lons, poly_lats)

        #------------------------------------------------ 
        # Try to determine the basin outlet coordinates
        #------------------------------------------------
        if (SAVE_OUTLET_INFO):
            if not(basin_in_dem(poly_x, poly_y, dem_info, hlr_value, 
                                counters, not_in_dem_vals )):
                continue
            lon, lat, col, row, zmin, one_zmin, \
               counters, not_in_dem_vals = \
            try_to_get_outlet_info( dem_unit, dem_info,
                                    hlr_zmin, hlr_value,
                                    hlr_zmin_match,
                                    poly_x,poly_y, poly_lons,poly_lats,
                                    counters, not_in_dem_vals,
                                    REPORT=REPORT)
            print('###### For hlr_value =', hlr_value)
            print('   Outlet col, row =', col, ',', row)
            print('   count =', hlr_count)
            print('   zmin, hlr_zmin =', zmin, ',', hlr_zmin)

            #---------------------------------------------- 
            # Note: For hlr_value = 115
            #       Outlet col, row = 2260 , 698
            #            count = 115
            #            zmin, hlr_zmin = 149 , 177
            # But outlet col should be 2261, which would
            #    give: zmin = 177 and count/area = 117
            #----------------------------------------------         
        if (REPORT):
            print()
            n_features = counters['total_features']
            print('####### Basin number:', n_features, '########')
            if (SAVE_OUTLET_INFO):
                print('At lon,lat =', lon, ',', lat)
                print('At col,row =', col, ',', row)
                print('    zmin   =', zmin)       # meters
            print('hlr_zmin   =', hlr_zmin)   # meters
            print('hlr_relief =', hlr_relief) # meters
            print('hlr_value  =', hlr_value, '(watershed ID)')
            print('hlr_count  =', hlr_count)
            print('hlr_slope  =', hlr_slope)  # m/m
            print('hlr_code   =', hlr_code)   # 0 to 20
            if (old_shapefile):
                print('hlr_area   =', hlr_area)   # m^2
                print('hlr_perim  =', hlr_perim)  # meters
            print()

        #------------------------
        # Get closest USGS site
        #------------------------
        if (SAVE_OUTLET_INFO):
            lon_str = '{x:.4f}'.format(x=lon)  # not accurate to 4 decimals
            lat_str = '{x:.4f}'.format(x=lat)
            closest_id, clon, clat, dmin = \
                usgs.get_closest_site(usgs_site_coords, 
                     lon_str, lat_str, REPORT=False)
            closest_site_id   = closest_id
            closest_site_dist = '{x:.2f}'.format(x=dmin)  # string

        #-----------------------------------
        # Save info for which HLR basins ?
        #-----------------------------------
        GOOD_INFO = ((hlr_code > 0) and (hlr_count > 1))
#         if (old_shapefile):
#             GOOD_INFO = ((hlr_area == 1e6 * hlr_count) and
#                          (hlr_code > 0) and (hlr_area > 1e6))
#         else:
#             GOOD_INFO = ((hlr_code > 0) and (hlr_count > 1))
        #--------------------------------------------------------        
        if (SAVE_OUTLET_INFO):
            GOOD_OUTLET = ((zmin == hlr_zmin) and (one_zmin))
            ## GOOD_OUTLET = one_zmin
            SAVE_INFO = (GOOD_INFO and GOOD_OUTLET)
        else:
            SAVE_INFO = GOOD_INFO
    
        #------------------------------------------------               
        # Save info for some HLR basins to new TSV file
        #------------------------------------------------
        if (SAVE_INFO):
            counters['saved_basins'] += 1

            if (SAVE_OUTLET_INFO):
                su.write_tsv_line( tsv_unit, attributes,
                   insert_key=insert_key,
                   new_values=[closest_id,closest_site_dist,
                               lon,lat,col,row,minlon,maxlon,minlat,maxlat])
            else:
                su.write_tsv_line( tsv_unit, attributes,
                   insert_key=insert_key,
                   new_values=[minlon,maxlon,minlat,maxlat])
            #--------------------------------------     
            # How many basins with each HLR code?
            #--------------------------------------
            hlr_histogram[ hlr_code ] += 1

        if (counters['total_features'] >= nb_max):
            break
            
    #--------------------------------
    # Print info, close files, etc.
    #--------------------------------
    print('========================================================')
    print('# total features    =', counters['total_features'])
    print('# excluded features =', counters['excluded_features'])
    print('# unique basins     =', counters['unique_basins'])
    print('# saved basins      =', counters['saved_basins'], '(written to TSV)')
    print()
    if (SAVE_OUTLET_INFO):    
        ## print('# hlr zmin matches  =', counters['hlr_zmin_match'])
        print('# hlr zmin matches  =', hlr_zmin_match.sum(), ' (unique)')
        print('# one zmin in poly  =', counters['one_zmin'])
        print('# one zmin & match  =', counters['one_zmin_match'])
        print('# two try count     =', counters['two_tries'])
        print('# not in dem count  =', counters['not_in_dem'])
        # print('# no zvalues count  =', counters['no_zvalues'])
        print()
    print('# code zero count   =', counters['code_is_zero'])
    print('# one cell count    =', counters['one_cell_basin'])

    if (old_shapefile):
        print('# count-area count  =', count_area_count)
        
    print()
    print('hlr_histogram (for saved basins) =')
    print( hlr_histogram )
    print()
    #-------------------------
    w = (val_histogram > 1)
    i = np.arange( max_val )
    print('# Skipped repeated watershed values.')
    print('# Repeated watershed IDs (VALUE) =', w.sum())
    print('vals =')
    print( i[w] )
    print('hist =')
    print( val_histogram[w] )
    if (SAVE_OUTLET_INFO):
        print('not_in_dem_vals =')
        print( not_in_dem_vals )
    print()

    dem_unit.close()
    tsv_unit.close()
    shape_unit = None  # a way to close file
    run_time = (time.time() - start_time)
    print('Run time =', run_time, '[sec]')
    print('Finished.')
    
#   create_tsv()
#---------------------------------------------------------------------
def basin_in_dem( poly_x, poly_y, dem_info, hlr_value,
                  counters, not_in_dem_vals, REPORT=True):

    #--------------------------------------------------------
    # Need this because Hawaii, and maybe some other basins
    # are not in the HYDRO1K DEM we have for North America
    # See Notes at the top of this file for more info.
    #--------------------------------------------------------
    xmin = dem_info['xmin']
    ymin = dem_info['ymin']
    xmax = dem_info['xmax']
    ymax = dem_info['ymax']
    #------------------------
    w1 = np.logical_and( poly_x >= xmin, poly_x < xmax)
    w2 = np.logical_and( poly_y >= ymin, poly_y < ymax)
    w  = np.logical_and( w1, w2 )
    in_dem = (w.sum() > 0)

    if not(in_dem):
        counters['not_in_dem'] += 1
        not_in_dem_vals.append( hlr_value )
        if (REPORT):
            print('##### WARNING #####')
            print('Basin is not in the North America DEM.')
            print('Watershed ID (VALUE) =', hlr_value)
            print('It is likely in Hawaii.')
            print('Skipping to next basin.')
            print()
    return in_dem

#     x2 = poly_x[w]
#     y2 = poly_y[w]

#   basin_in_dem()
#---------------------------------------------------------------------
def snap_polygon_points( poly_x, poly_y, hlr_value, old_shapefile ):

    #----------------------------------
    # Round to nearest 1000 or 500 ??
    #--------------------------------------------
    # Most poly_x values end in 16.5 or 16.625,
    # and most poly_y values end in 190.5.
    #--------------------------------------------
#         poly_x = np.around(poly_x - 500, -3)
#         poly_y = np.around(poly_y - 500, -3)
#         poly_x = np.around(poly_x, -3) + 500
#         poly_y = np.around(poly_y, -3) - 500
    if (old_shapefile):
        poly_x += 16.5    # (old shapefile)
        w0 = ((poly_x % 1000) != 0)
        n0 = w0.sum()
        if (n0 > 0):
            poly_x[ w0 ] += 0.125
        poly_y -= 190.5
    else: 
        poly_x += 16.6    # (new shapefile)  
        poly_y -= 190.4
#         poly_x += (16.6  + 500)   # (new shapefile)  
#         poly_y -= (190.4 + 500)
    #----------------------------
    w1 = ((poly_x % 1000) != 0)   # boolean array
    n1 = w1.sum()
    w2 = ((poly_y % 1000) != 0)
    n2 = w2.sum()
    if (n1 > 0):
        print('### WARNING: poly_x is not a multiple of 1000.')
        print('    poly_x[w1] =', poly_x[w1] )
        print('    hlr_value  =', hlr_value)
    if (n2 > 0):
        print('### WARNING: poly_y is not a multiple of 1000.')
        print('    poly_y[w2] =', poly_y[w2] )
        print('    hlr_value  =', hlr_value)
        
    return poly_x, poly_y
    
#   snap_polygon_points()     
#---------------------------------------------------------------------
def try_to_get_outlet_info( dem_unit, dem_info, hlr_zmin,
                            hlr_value, hlr_zmin_match,
                            x1,y1, lons,lats,
                            counters, not_in_dem_vals,
                            REPORT=False):

    #----------------------------------------------------------
    # The source DEM is the GTOPO30 HYDRO1K DEM and uses a
    # Lambert Azimuthal EqArea projection as described above.
    # Can index this 1-km DEM with the Lambert xy coords
    # to get an elevation, then compare to MINELE attribute
    # as a way to get the outlet xy.  Or just find the xy
    # pair for each polygon with the lowest elevation.
    # But the DEM will be pretty big.  The HYDRO1K DEM for
    # North America includes Alaska, CONUS, PR, etc. but not
    # Hawaii. It has 9102 cols and 8384 rows.  Need to
    # convert Lambert xy coords to DEM column and row.
    #----------------------------------------------------------
    # See: http://www.ncgia.ucsb.edu/SANTA_FE_CD-ROM/
    #       sf_papers/verdin_kristine/santafe2.html
    # HYDRO1K DEM:  https://d9-wret.s3.us-west-2.amazonaws.com/
    #    assets/palladium/production/s3fs-public/atoms/files/
    #    HYDRO1kDocumentationReadMe.pdf
    #----------------------------------------------------------
    ncols  = dem_info['ncols']
    nrows  = dem_info['nrows']
    xmin   = dem_info['xmin']
    ymax   = dem_info['ymax']
    dx     = dem_info['dx']
    dy     = dem_info['dy']
    nodata = dem_info['nodata']

    #-------------------------------------------------------    
    # Get cols and rows of dem grid cells on basin polygon
    #-------------------------------------------------------
    cols, rows = get_cols_and_rows( x1, y1, xmin, ymax, dx, dy,
                     round_method='floor', print_c2_r2=False )
    ## break
           
    #---------------------------------------------------       
    # Get elevation (z) values for every polygon point
    # then get zmin and corresponding lon and lat.
    #---------------------------------------------------
    zvals = get_polygon_zvals( cols, rows, ncols, dem_unit,
                               REPORT=REPORT )

    #------------------------------------------      
    # Get lon and lat for grid cell with zmin
    # Recall grid cell size is 1 km (HYDRO1K)
    #------------------------------------------      
    zmin, zmin_arg, lon, lat, one_zmin = get_zmin_lon_and_lat( zvals,
          lons, lats, nodata, REPORT=REPORT)
#     if (zmin != hlr_zmin):
#         #----------------------------------------------------
#         # Note: This may only happen for "fragmented" basin
#         #       polygons where hlr_value is repeated.
#         #----------------------------------------------------
#         print('zmin != hlr_zmin for VALUE =', hlr_value)
#     else:
#         print('zmin == hlr_zmin for VALUE =', hlr_value)
                 
    #--------------------------------------------------------        
    # If zmin doesn't match hlr_zmin, then try again with
    # round_method = 'ceil'.  Still don't get all matches.
    # For old_shapefile,
    # Got 168 matches out of first 200 polygons (16% fails)
    # Also got 20/20, 28/30, 42/50.
    #--------------------------------------------------------
    # Note: Could have one_zmin=True and (hlr_zmin != zmin)
    #-------------------------------------------------------- 
#     if (zmin != hlr_zmin):
#         counters['two_tries'] += 1
#         cols, rows = get_cols_and_rows( x1, y1, xmin, ymax, dx, dy,
#                          round_method='ceil', print_c2_r2=False )
#         zvals = get_polygon_zvals( cols, rows, ncols, dem_unit,
#                                    REPORT=REPORT )    
#         zmin, zmin_arg, lon, lat, one_zmin = get_zmin_lon_and_lat( zvals,
#               lons, lats, nodata, REPORT=REPORT )

    #-----------------------------------        
    # Count when zmin matches hlr_zmin
    #--------------------------------------------------------
    # The polygon with HLR VALUE=25 is fragmented/repeated,
    # and zmin == hlr_zmin for some pieces and not others.
    #-------------------------------------------------------- 
    if (zmin == hlr_zmin):
        hlr_zmin_match[ hlr_value ] = True  ################
        counters['hlr_zmin_match'] += 1
        if (one_zmin):
            counters['one_zmin_match'] += 1
    if (one_zmin):
        counters['one_zmin'] += 1
    #------------------------------
    col = cols[ zmin_arg ]
    row = rows[ zmin_arg ]

    return lon, lat, col, row, zmin, one_zmin, \
           counters, not_in_dem_vals
 
#   try_to_get_outlet_info()
#---------------------------------------------------------------------
def get_cols_and_rows( x1, y1, xmin, ymax, dx, dy,
                       round_method=None, print_c2_r2=None):

    #-----------------------------------------------------------    
    # Note: The polygon vertices are in clockwise order and
    # they trace around the outer edges of the DEM grid cells
    # that are contained in a given basin.  In order to get
    # the cols and rows of grid cells in the basin, we note
    # the following.
    #
    # (1) If next vertex is to the right of the last one, then
    #     the grid cell below this edge is on the boundary.
    # (2) If next vertex is to the left of the last one, then
    #     the grid cell above this edge is on the boundary.
    # (3) If next vertex is below the last one, then the grid
    #     cell to the left of this edge is on the boundary.
    # (4) If next vertex is above the last one, then the grid
    #     cell to the right of this edge is on the boundary.
    #
    # A given grid cell may be included more than once, so
    # need to remove duplicates (non-unique pairs) at the end.
    #
    # This won't work as written for hlrus4.shp, because it
    # doesn't include all boundary vertices.      
    #------------------------------------------------------------
    # Formula 1 comes from taking the smallest x1
    # value of x1 = (xmin + dx/2) to get 0th col.
    # Formula 2 comes from taking the largest y1
    # value of y1 = (ymax - dy/2) to get 0th row.
    #------------------------------------------------
    DEBUG = False

    #-----------------------------------------------------
    # Note: The first hlr_values are for Alaska and 
    #       there are many fragmented basins along the
    #       northern coastline.  The following info is
    #       for the original shapefile: hlrus.shp,
    #       and with only one "try".
    #       get_cols_and_rows() does not yet support
    #       hlrus4.shp.
    #-----------------------------------------------------
    # For xoffset=500, yoffset=0:
    #   zmin = hlr_zmin matches are 56 of first 200.
    #   zmin = hlr_zmin matches are 119 of first 500.
    #      one_zmin_in_poly   = 146
    #      one_zmin_and_match = 119
    #   zmin = hlr_zmin matches are 249 of first 1000.
    #      hlr zmin matches   = 698 (incl. more than one)
    #      one_zmin_in_poly   = 318
    #      one_zmin_and_match = 249
    #   zmin = hlr_zmin matches are 863 of first 3000 "features".
    #      unique basins      = 2582
    #      excluded features  = 304
    #      hlr zmin matches   = 2310 (incl. more than one)
    #      one_zmin_in_poly   = 1124
    #      one_zmin_and_match = 863
    #   zmin = hlr_zmin matches are 2589 of first 10000 "features".
    #      unique basins      = 8495
    #      excluded features  = 1117
    #      hlr zmin matches   = 6888 (incl. more than one)
    #      one_zmin_in_poly   = 4014 
    #      one_zmin_and_match = 2589
    #   This includes hlr_values = 115 and 50.
    #   But not hlr_value = 64 or 148. (32767 vs. 0)
    #-----------------------------------------------------
    # For xoffset=500, yoffset=-500:
    #   zmin = hlr_zmin matches are 53 of first 200.
    #   This includes hlr_value = 115.
    #-----------------------------------------------------
    # For xoffset=1000, yoffset=0:
    #   zmin = hlr_zmin matches are 43 of first 200.
    #   This doesn't include hlr_value = 115.
    #-----------------------------------------------------
    # For xoffset=400, yoffset=0:
    #   zmin = hlr_zmin matches are 42 of first 200.
    #   This doesn't include hlr_value = 115.
    #-----------------------------------------------------
    # For xoffset=0, yoffset=0:
    #   zmin = hlr_zmin matches are 42 of first 200.
    #   This doesn't include hlr_value = 115.
    #-----------------------------------------------------
    # For xoffset=600, yoffset=0:
    #   zmin = hlr_zmin matches are 34 of first 200.
    #   This doesn't include hlr_value = 115.
    #----------------------------------------------------
    # For xoffset=500, yoffset=500:
    #   zmin = hlr_zmin matches are 34 of first 200.
    #   This doesn't include hlr_value = 115.
    #----------------------------------------------------
    ### CHECK hlr_value = 50
    ## xoffset = 0
    ## xoffset = 1000
    xoffset =  500  # grid cell centers
    # yoffset = -500
    yoffset = 0
    ## yoffset = 500
    cell_x = np.zeros( x1.size, dtype='int32')
    cell_y = np.zeros( y1.size, dtype='int32')
    x1_next = np.roll(x1, -1)  # last item rolls to first position
    y1_next = np.roll(y1, -1)
    #------------------------------------
    w1 = (x1_next > x1)
    w1[-1] = False  # exclude last item
    n1 = w1.sum()
    i1 = 0
    i2 = i1 + n1
    cell_x[i1:i2] = x1[ w1 ] + xoffset
    cell_y[i1:i2] = y1[ w1 ] - yoffset
    #------------------------------------
    w2 = (x1_next < x1)
    w2[-1] = False
    n2 = w2.sum()
    i1 = i2
    i2 = i1 + n2
    cell_x[i1:i2] = x1[ w2 ] - xoffset
    cell_y[i1:i2] = y1[ w2 ] + yoffset
    #------------------------------------
    w3 = (y1_next < y1)
    w3[-1] = False
    n3 = w3.sum()
    i1 = i2
    i2 = i1 + n3
    cell_x[i1:i2] = x1[ w3 ] - xoffset
    cell_y[i1:i2] = y1[ w3 ] - yoffset
    #------------------------------------
    w4 = (y1_next > y1)
    w4[-1] = False
    n4 = w4.sum()
    i1 = i2
    i2 = i1 + n4
    cell_x[i1:i2] = x1[ w4 ] + xoffset
    cell_y[i1:i2] = y1[ w4 ] + yoffset
    #------------------------------------
    cell_x = cell_x[:i2]  # remove extras
    cell_y = cell_y[:i2]
   
    #----------------------------------------  
    # Keep only unique cell_x, cell_y pairs
    #----------------------------------------    
    xstr = cell_x.astype('<U16')
    ystr = cell_y.astype('<U16')
    str_pairs = np.char.add(xstr, ystr)
    uniq_pairs, indices_of_uniq = np.unique(str_pairs, return_index=True) 
    cell_x = cell_x[ indices_of_uniq ]
    cell_y = cell_y[ indices_of_uniq ]
    #-----------------------------------------------
    # Recall that row 0 is top row; goes with ymax
    #-----------------------------------------------
    c1   = (cell_x - xmin) / dx  # formula1, float
    r1   = (ymax - cell_y) / dy  # formula2, float
    cols = np.int32(c1)
    rows = np.int32(r1)
    
    if (DEBUG):
        print('### x1.size     =', x1.size)
        print('### len(cell_x) =', len(cell_x))
        # print('cell_x =')
        # print( cell_x )
        print('### size(cols)  =', cols.size)
        print()

    return cols, rows
    
#   get_cols_and_rows()                 
#---------------------------------------------------------------------
def get_cols_and_rows2( x1, y1, xmin, ymax, dx, dy,
                       round_method='floor',
                       int_method='int32',
                       print_c2_r2=False ):

    #------------------------------------------------              
    # Notice that polygon coords always end in ".5"
    # But they should be multiples of either 1000,
    # or 500 if they are grid cell centers??
    #------------------------------------------------
    # Formula 1 comes from taking the smallest x1
    # value of x1 = (xmin + dx/2) to get 0th col.
    # Formula 2 comes from taking the largest y1
    # value of y1 = (ymax - dy/2) to get 0th row.
    #------------------------------------------------

    #-----------------------------------------    
    # Get 39/50 w/ floor. 30/50 w/ ceil.
    # 39/50 w/ trunc. 30/50 w/ around.
    # 39/50 w/ ceilfloor, 32/50 w/ floorceil 
    #-----------------------------------------
#     c1 = (x1 - (xmin - dx/2)) / dx  # float
#     r1 = ((ymax + dy/2) - y1) / dy  # float
    
    #-----------------------------------------    
    # Get 21/50 w/ floor. 39/50 w/ ceil.
    # 21/50 w/ trunc. 39/50 w/ around.
    # 25/50 w/ ceilfloor, 32/50 w/ floorceil 
    #-----------------------------------------
#     c1 = (x1 - (xmin + dx/2)) / dx  # float
#     r1 = ((ymax - dy/2) - y1) / dy  # float
    #-----------------------------------------
    # Get 39/50 w/ floor. 30/50 w/ ceil.
    # 39/50 w/ trunc. 39/50 w/ around.
    # 39/50 w/ ceilfloor, 32/50 w/ floorceil 
    #-----------------------------------------
    c1 = (x1 - xmin) / dx  # float
    r1 = (ymax - y1) / dy  # float
    #-----------------------------------------    
    # Get 39/50 w/ floor. 30/50 w/ ceil.
    # */50 w/ trunc. */50 w/ around.
    # */50 w/ ceilfloor, */50 w/ floorceil 
    #-----------------------------------------
#     c1 = (x1 - (xmin + 0.001)) / dx  # float
#     r1 = ((ymax - 0.001) - y1) / dy  # float
    
    #--------------------------    
    # These return float type
    #--------------------------
    if (round_method == 'floor'):
        #-----------------------------------------
        # zmin matches hlr_zmin for basins:
        #   1,2,3,#,5,6,7,8,9,#,11,12,13,14,15,
        #   #,17,18,19,20.  But not 4,10, or 16. 
        #-----------------------------------------
        c2 = np.floor( c1 )
        r2 = np.floor( r1 )
    elif (round_method == 'ceil'):
        #-----------------------------------------
        # zmin matches hlr_zmin for basins:
        #   #,#,3,4,5,6,7,8,9,10,11,12,13,14,15,
        #   16,17,18,19,#.  But not 1,2,or 20.
        #-----------------------------------------
        c2 = np.ceil( c1 )
        r2 = np.ceil( r1 )
    elif (round_method == 'trunc'):
        c2 = np.trunc( c1 )
        r2 = np.trunc( r1 )    
    elif (round_method == 'around'):
        #----------------------------------------
        # zmin matches hlr_zmin for basins:
        #   1,2,3,#,5,6,7,8,9,#,11,12,13,14,15,
        #   #,17,18,19,20. But not 4,10, or 16.
        #----------------------------------------
        c2 = np.around( c1 )
        r2 = np.around( r1 )
    elif (round_method == 'ceilfloor'):
        #----------------------------------------
        # zmin matches hlr_zmin for basins:
        #----------------------------------------
        c2 = np.ceil( c1 )
        r2 = np.floor( r1 )
    elif (round_method == 'floorceil'):
        #----------------------------------------
        # zmin matches hlr_zmin for basins:
        #----------------------------------------
        c2 = np.floor( c1 )
        r2 = np.ceil( r1 )
    else:
        c2 = np.around( c1 )
        r2 = np.around( r1 )
        
    if (print_c2_r2):
        print('c2 =', c2)
        print('r2 =', r2)

    #------------------------------    
    # Convert to int32 and return
    #------------------------------
    cols = np.int32( c2 )
    rows = np.int32( r2 )
    return cols, rows    

#   get_cols_and_rows()
#---------------------------------------------------------------------
def get_polygon_zvals( cols, rows, ncols, dem_unit,
                       REPORT=True ):

    #--------------------------------------------
    # Note: Both methods here avoid reading the
    #       very large DEM into memory.
    #--------------------------------------------
    zvals = np.zeros( cols.size, dtype='int16' )
    
    for k in range(cols.size):
        ID = (rows[k] * ncols) + cols[k]
        offset = np.int32(2) * ID  # (2-bytes per elev value)
        # print('k, offset =', k, offset)
        #----------------------
        # Use the seek method
        #----------------------
        dem_unit.seek( offset, 0 )
        b = dem_unit.read(2)  # read 2 bytes
        val = int.from_bytes(b, byteorder='big') # big=MSB
        zvals[k] = val
        # print('zval =', val)

        #---------------------------------------
        # Use the fromfile method (also works)
        #---------------------------------------
        # val = np.fromfile(dem_unit, dtype='int16', count=1, offset=offset)
        # val.byteswap(inplace=True)  # DEM has MSB byte order
        # dem_unit.seek( 0, 0 )  # reset to start of file
        # zvals[k] = val[0]
        ## print('type(zval) =', type(val))
        ## print('zval.dtype =', val.dtype)
        ## print('zval =', val[0])

    if (REPORT):
        print('zvals =', zvals)  #######
        print('n_pts =', zvals.size)
        # Note: Use this if first and last zvals are the same
        #       meaning last poly vertex is repeated.
        ## print('n_pts =', zvals.size - 1)
        print()
    return zvals
 
#   get_polygon_zvals()
#---------------------------------------------------------------------
def check_lon_lats( cols, rows, ncols, nrows, x2, y2):

    # Now obsolete?
    #----------------------------------------------
    # Need this because Hawaii is not in this DEM
    #----------------------------------------------
    w1 = np.logical_and( cols >= 0, cols < ncols)
    w2 = np.logical_and( rows >= 0, rows < nrows)
    w  = np.logical_and( w1, w2 )
    cols = cols[w]
    rows = rows[w]
    lons = x2[w]
    lats = y2[w]
    return lons, lats, cols, rows
    
#   check_lons_lats()
#---------------------------------------------------------------------
def get_zmin_lon_and_lat( zvals, lons, lats, nodata,
                          REPORT=True ):

    #------------------------------------------------------
    # This version does not use np.argmin(), since that
    # function only returns index of the first min value.
    #------------------------------------------------------
    zvals[ zvals == nodata ] = 32767  # Need this
    zmin = zvals.min()
    ivals = np.arange( zvals.size )
    w = (zvals == zmin)
    zmin_args = ivals[w]
    if (zmin_args.size == 1):
        one_zmin = True
    else:
        if (REPORT):
            print('###  WARNING  ###')
            print('Unsure about outlet lat and lon.')
            print('Multiple cells have min elevation')
            print()
        one_zmin = False
    zmin_arg = zmin_args[0]
    outlon   = lons[ zmin_arg ]  # for outlet
    outlat   = lats[ zmin_arg ]
    
    #--------------------------------------- 
    # Round to 4 places after decimal;
    # roughly 11 meters.  Can only do
    # in-place w/out keyword for ndarrays.
    #---------------------------------------
    outlon = np.around(outlon, decimals=4)
    outlat = np.around(outlat, decimals=4)

    return zmin, zmin_arg, outlon, outlat, one_zmin

#   get_zmin_lon_and_lat()           
#---------------------------------------------------------------------
# def get_zmin_lon_and_lat1( zvals, lons, lats, nodata ):
# 
#     zvals[ zvals == nodata ] = 32767  # Need this
#     zmin = zvals.min()
#     zmin_args = zvals.argmin()
#     if (zmin_args.size == 1):
#         zmin_arg = zmin_args
#     else:
#         #-----------------------------------------
#         # NOTE!! This will not necessarily give
#         #        the correct outlet lon and lat!
#         #-----------------------------------------
#         print('###  WARNING  ###')
#         print('Unsure about outlet lat and lon.')
#         print()
#         zmin_arg = zmin_args[0]
#     lon = lons[ zmin_arg ]
#     lat = lats[ zmin_arg ]
#     return zmin, lon, lat
# 
# #   get_zmin_lon_and_lat1()           
#---------------------------------------------------------------------
def convert_vat_adf_to_csv( vat_adf_file=None, csv_file=None):

    #----------------------------------------------------------
    # Note: This function was adapted from:
    #    https://gis.stackexchange.com/questions/84700/
    #    showing-only-raster-attribute-table-using-gdal/84729
    # VAT = "Value Attribute Table" (value vs. vector)
    #----------------------------------------------------------
    # Also see:
    # https://gdal.org/drivers/raster/arcinfo_grid_format.html
    #----------------------------------------------------------
    # This function seems to work but gives an error message:
    #
    # ERROR 5: OSRCalcInvFlattening(): Wrong input values
    # This error is triggered immediately when calling
    # gdal.Open() or ogr.Open().  Note that there is a file
    # called prj.adf in the same directory.  Can open this
    # file with QGIS, then click in bottom right on the
    # "Unknown CRS" button to see this:
    #     BASEGEOGCRS["unknown",
    #         DATUM["unknown",
    #             ELLIPSOID["unknown",6370997,0,
    # It seems to be using a sphere vs. ellipsoid, in which
    # case InverseFlattening may be undefined.
    # If we hide the file: prj.adf, there is no error msg,
    # and the resulting CSV file is identical (says diff),
    # so it seems we can safely ignore it.
    #
    # See "__README_RAT_TO_CSV.txt." about error message.
    #----------------------------------------------------------
    new_dir    = dtu.get_new_data_dir( 'USGS_HLR' )
    data_dir   = dtu.get_data_dir( 'USGS_HLR' )
    raster_dir = data_dir + 'arctar00000/hlrus/'
    
    if (vat_adf_file is None):
        vat_adf_file = raster_dir + 'vat.adf'
    if (csv_file is None):
        csv_file = new_dir + 'vat_adf.csv'

    #------------------------------------------
    # Note: ogr is a vector library;
    #       gdal is a raster library
    #------------------------------------------
    ds = gdal.Open( vat_adf_file )
    # print(ds)
    # return  # shows that error is triggered immediately.
    #------------------------------------
    # Not sure if ogr can open vat_adf.
    #------------------------------------
#     ds = ogr.Open( vat_adf_file )
#     print(ds)  # Gives "None"
#     return  # shows that error is triggered immediately.
    
    rat = ds.GetRasterBand(1).GetDefaultRAT()
    
    with open(csv_file, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)

        # Write out column headers
        icolcount=rat.GetColumnCount()
        cols=[]
        for icol in range(icolcount):
            cols.append(rat.GetNameOfCol(icol))
        csvwriter.writerow(cols)

        # Write out each row.
        irowcount = rat.GetRowCount()
        for irow in range(irowcount):
            cols=[]
            for icol in range(icolcount):
                itype=rat.GetTypeOfCol(icol)
                if itype==gdal.GFT_Integer:
                    value='%s'%rat.GetValueAsInt(irow,icol)
                elif itype==gdal.GFT_Real:
                    value='%.16g'%rat.GetValueAsDouble(irow,icol)
                else:
                    value='%s'%rat.GetValueAsString(irow,icol)
                cols.append(value)
            csvwriter.writerow(cols)
            
    # Note: "with" closes the new CSV file.
    # ds = None  # (Close vat.adf file)

#   convert_vat_adf_to_csv()
#---------------------------------------------------------------------
def get_hlr_info_dict( in_csv_file='vat_adf.csv', CODE_ONLY=False ): 

    #-------------------------------------------------------------
    # Note: The file vat_adf.csv contains attributes for the
    #       43931 HLR basins, obtained from the raster version
    #       of the HLR data set in the arctar00000/hlrus folder.
    #       See Notes for convert_vat_adf_to_csv, above.
    #----------------------------------------------------------------------------------
#     Shapefile attributes are defined here:
#     https://water.usgs.gov/GIS/metadata/usgswrd/XML/hlrus.xml
# 
#     Value     = Watershed ID, in range [1,44594]  ######## IN MULTIPLE ROWS
#     Count     = Watershed area in square kilometers, in range [31, 23476]
#                 (really a # of grid cells)
#     Aqpermnew = Aquifer permeability class (1-7, lowest-highest)
#     Slope     = Slope in percent rise
#     Tave      = Mean annual temperature in degrees Fahrenheit, in [18.7, 75.0]
#     Ppt       = Mean annual precipitation in inches per year, in [3.2, 136.8]
#     PET       = Mean annual potential evapotrans. in inches per year, in [12.8, 52,7]
#     Sand      = Percent sand in soil, in [2.5, 100]
#     PMPE      = Mean annual precipitation minus potential evapotranspiration
#                 in inches per year, in [-48.6, 99.7]
#     Minele    = Minimum elevation in watershed in meters, in [-79, 3169]
#     Relief    = Maximum elevation in watershed minus minimum elevation
#                 in watershed in meters, in range [0, 5618]
#     Pflattot  = Total percent flat land (slope less than 1 percent)
#                 in watershed, in range [0,100]
#     Pflatlow  = Percent flat land (slope less than 1 percent) in watershed
#                 lowland (elevation less than midpoint between minimum
#                 and maximum elevation), in range [0, 100]
#     Pflatup   = Percent flat land (slope less than 1 percent) in watershed
#                 upland (elevation greater than or equal to midpoint
#                 between minimum and maximum elevation), in range [0, 100]
#     Hlr       = Hydrologic landscape region identification number, in [1,20]
    #----------------------------------------------------------------------------------
    new_dir   = dtu.get_new_data_dir( 'USGS_HLR' )
    dict_file = 'HLR_basin_info_dict.pkl'
    file_path = new_dir + dict_file
    if (os.path.exists( file_path ) and not(CODE_ONLY)):
        print('Reading saved HLR basin info dictionary...')
        file_unit = open( file_path, 'rb')
        hlr_info_dict = pickle.load( file_unit )
        file_unit.close()
        print('Finished.')
        print()
        return hlr_info_dict

    delim = ','  # CSV file
    print('Working...')
    hlr_basin_count = 0
 
    #------------------
    # Open input file
    #------------------
    in_csv_path = new_dir + in_csv_file
    in_csv_unit = open(in_csv_path,  'r')

    #-----------------------    
    # Skip the header line
    #-----------------------
    header = in_csv_unit.readline()
    hlr_info_dict = dict()

    #------------------------------
    # Map zero to zero for nodata
    #------------------------------
    if (CODE_ONLY):
        hlr_info_dict[0] = 0
    else:
        hlr_info_dict[0] = {'area': '0.0', 'slope':'0.0', 'T_avg':'0.0',
            'P_avg':'0.0', 'PET_avg':'0.0', 'min_elev':'0', 'relief':'0',
            'hlr_code':0}  # integer, not string

    while (True):
        line = in_csv_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        values = line.split( delim )
        hlr_basin_count  += 1

        value    = values[0].strip()   # this is the basin ID
        area     = values[1].strip()   # this is COUNT; each cell is 1 km2.
        slope    = values[3].strip()   # [%]
        T_avg    = values[4].strip()   # [deg_F]
        P_avg    = values[5].strip()   # [in/yr]
        PET_avg  = values[6].strip()   # [in/yr]
        min_elev = values[9].strip()   # [meters]
        relief   = values[10].strip()  # [meters]
        hlr_code = values[14].strip()  # in {1,...,20}

        #------------------------------------------        
        # Save info in dictionary w/ value as key
        # Save each entry as original string or
        # as a re-formatted string.
        #------------------------------------------
#         hlr_info_dict[ value ] = { 'area':area, 'slope':slope, 'T_avg':T_avg,
#             'P_avg':P_avg, 'PET_avg':PET_avg, 'min_elev':min_elev,
#             'relief':relief, 'hlr_code':hlr_code }
        #------------------------------------------
        val_num = np.int32(value)
        if (CODE_ONLY):
            hlr_info_dict[ val_num ] = np.int16(hlr_code)
        else:      
            hlr_info_dict[ val_num ] = {
            'area':   '{x:.2f}'.format(x=np.float32(area)),  # int to float, [km2]
            'slope':  '{x:.4f}'.format(x=np.float32(slope)),
            'T_avg':  '{x:.4f}'.format(x=np.float32(T_avg)),
            'P_avg':  '{x:.4f}'.format(x=np.float32(P_avg)),
            'PET_avg':'{x:.4f}'.format(x=np.float32(PET_avg)),
            #-----------------------------------------------------
            'min_elev': np.int16(min_elev),    # an integer
            'relief':   np.int16(relief),      # an integer
            'hlr_code': np.int16(hlr_code) }   # an integer
                
    SAVE_TO_FILE = not(CODE_ONLY)
    if (SAVE_TO_FILE):
        file_unit = open( file_path, 'wb' )
        pickle.dump( hlr_info_dict, file_unit, protocol=pickle.HIGHEST_PROTOCOL)
        file_unit.close()
        print('Saved HLR basin info dictionary to file:')
        print('  ' + file_path)
        print()
        
    #-------------------------------   
    # Close input and output files
    #-------------------------------
    in_csv_unit.close()
    print('Number of HLR basins =', hlr_basin_count )
    print('Finished.')
    
    return hlr_info_dict

#   get_hlr_info_dict()
#---------------------------------------------------------------------   
def create_hlr_grid( value_tif_file='HLR_value_grid_9106x6855_1K.tif',
                     code_tif_file='HLR_code_grid_9106x6855_1K.tif',
                     code_rtg_file='HLR_code_grid_9106x6855_1K.rtg',
                     REPORT=True ):

    #---------------------------------------------------------
    # Note: We can export a grid of HLR values from QGIS.
    #       This routine creates a similar grid of HLR codes
    #       by mapping each HLR value to its code.
    #       Using this grid, we can assign an HLR code to
    #       any lon/lat pair that lies within the grid.
    #---------------------------------------------------------
    # max_hlr_value = 44594
    # value_nodata  = -2147483647.0 
    # code_nodata   = 0

    hlr_info_dict = get_hlr_info_dict( CODE_ONLY=True )
    
    #------------------------------    
    # Read the grid of HLR values
    #------------------------------
    new_dir = dtu.get_new_data_dir( 'USGS_HLR' )
    raster_dir = new_dir + 'Raster/'
    #------------------------------------------------
    value_tif_path = (raster_dir + value_tif_file)
    code_tif_path  = (raster_dir + code_tif_file)
    code_rtg_path  = (raster_dir + code_rtg_file)
    value_grid = regrid.read_geotiff( in_file=value_tif_path )
                        ##### MAKE_RTG=True, rtg_file=code_rtg_path )
    ## grid_shape = value_grid.shape

    #------------------------------    
    # Map HLR values to HLR codes
    #------------------------------ 
    print('Mapping HLR values to HLR codes...')
    w0 = (value_grid < 1)
    value_grid[ w0 ] = 0
    hlr_map  = np.vectorize(hlr_info_dict.__getitem__) 
    hlr_grid = hlr_map( value_grid )

    #-------------------------------------------------------------------
    # np.vectorize() has an "otypes" keyword to set output type
    #    but this doesn't work:
    # hlr_map = np.vectorize(hlr_info_dict.__getitem__, otypes='int16')
    # So convert to int16 afterwards.  Also see save_grid_to_geotiff().
    #--------------------------------------------------------------------    
    hlr_grid = hlr_grid.astype('int16')  ###############

    if (REPORT):
        print('hlr_grid.shape =', hlr_grid.shape )
        print('hlr_grid.dtype =', hlr_grid.dtype )
        print('hlr_grid.min() =', hlr_grid.min() )
        print('hlr_grid.max() =', hlr_grid.max() )

    #-------------------------------------
    # Save new grid of HLR codes to file
    #-------------------------------------
    # bounds = [ulx, lry, lrx, uly]  # [xmin, ymin, xmax, ymax]
    ulx  = -6175016.6
    uly  = 4324190.4
    xres = 1000.0
    yres = 1000.0
    regrid.save_grid_to_geotiff( code_tif_path, hlr_grid,
                        ulx, uly, xres, yres, 
                        dtype='int16', nodata=0)
    print('Finished.')

#   create_hlr_grid()
#---------------------------------------------------------------------
def get_hlr_grid_info():

    #------------------------------------
    # bounds = [ulx, lry, lrx, uly]
    #        = [xmin, ymin, xmax, ymax]
    #--------------------------------------------------------
    # data type may be LONG in GeoTIFF, but only need BYTE.
    #--------------------------------------------------------
    grid_info = dict()
    grid_info['ncols'] = 9106
    grid_info['nrows'] = 6855
    grid_info['dx']    = 1000.0
    grid_info['dy']    = 1000.0
    grid_info['xmin']  = -6175016.6
    grid_info['xmax']  =  2930983.4
    grid_info['ymin']  = -2530809.6
    grid_info['ymax']  =  4324190.4
    grid_info['pixel_geom'] = 1   # Code used in RT files (fixed-length)
    grid_info['zres']     = 1
    grid_info['z_units']  = 'none'
    grid_info['xy_units'] = 'meters'
    #------------------------------------------------------------
    grid_info['source']     = 'Wolock et al. (2004)'
    grid_info['code_file']  = 'HLR_code_grid_9106x6855_1K.tif'
    grid_info['value_file'] = 'HLR_value_grid_9106x6855_1K.tif'
    
    return grid_info

#   get_hlr_grid_info()
#---------------------------------------------------------------------
def Lambert_Azimuthal_Equal_Area_XY( lon, lat,
            lon0=-100.0, lat1=45.0):

    #--------------------------------------------------
    # Note: This function is also in projections.py.
    #--------------------------------------------------
    # Note: Formulas from Snyder's Book, p. 185.
    #--------------------------------------------------
    # R = radius of Clarke 1866 Authalic Sphere
    # Central_Meridian   = lon0 = -100.0
    # Latitude of Origin = lat0 = 45.0
    #--------------------------------------------------
    # Equivalent PROJ4 command string:
    # +proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0
    #    +a=6370997 +b=6370997 +units=m +no_defs
    #--------------------------------------------------
    # Was able to test using this PROJ4 string at:
    #    https://mygeodata.cloud/cs2cs/
    # and got perfect agreement.
    #--------------------------------------------------     
    R   = 6370997.0   # [meters]
    d2r = (np.pi / 180.0)
    lon_rad  = lon  * d2r   # psi in Snyder
    lon0_rad = lon0 * d2r
    lon_diff_rad = (lon_rad - lon0_rad)
    lat_rad  = lat  * d2r   # phi in Snyder
    lat1_rad = lat1 * d2r

    #----------------------------------------------------------
    # Note: Origin of map projection is at (lon0, lat1).
    #       x is positive to the east of the map origin
    #       and increases from west to east as it should.
    #       (-100, 45)   -> (0.0, 0.0)
    #       (-100.5, 45) -> (-39313.012957, 121.2945992)
    #       (-99.5, 45)  -> ( 39313.012957, 121.2945992)
    #----------------------------------------------------------
    # Note: y is positive to the north of the map origin
    #       and increases from south to north as it should.
    #       (-100, 45.5) -> (0.0,  55597.26072638)
    #       (-100, 44.5) -> (0.0, -55597.26072638)
    #----------------------------------------------------------   
    kterm1 = np.sin(lat1_rad) * np.sin(lat_rad)
    kterm2 = np.cos(lat1_rad) * np.cos(lat_rad) * np.cos(lon_diff_rad)
    kp     = np.sqrt(2 / (1 + kterm1 + kterm2))  # k_prime in Snyder
    #----------------------------------------------------------
    x      = R * kp * np.cos(lat_rad) * np.sin(lon_diff_rad)
    #----------------------------------------------------------
    yterm1 = np.cos(lat1_rad) * np.sin(lat_rad)
    yterm2 = np.sin(lat1_rad) * np.cos(lat_rad) * np.cos(lon_diff_rad)
    y      = R * kp * (yterm1 - yterm2)
    return x,y

#   Lambert_Azimuthal_Equal_Area_XY()
#---------------------------------------------------------------------
def get_hlr_code_for_point( lon, lat, hlr_grid=None, grid_info=None,
                            out_prj_file=None, USE_GDAL=False,
                            REPORT=False, WARNINGS=True):

    #--------------------------------------------
    # Note: float() can also handle '1e5', etc.
    #--------------------------------------------
    try:
        lon2 = float(lon)
        lat2 = float(lat)
    except:
        ## print('### ERROR: lon or lat is not a number.')
        return '0'
    #---------------------------------------------
    if (lon == '-9999') or (lat == '-9999'):
        ## print('### ERROR: lon or lat is -9999.')
        return '0'

    #----------------------------------------------------------
    # Note: isnumeric() returns False if there is anything
    #       other than a digit, such as minus sign or decimal
    #------------------------------------------------------
#     lon2 = lon.replace('-', '').replace('.', '')
#     lat2 = lat.replace('-', '').replace('.', '')
#     if not(lon2.isnumeric() and lat2.isnumeric()):
#         print('### ERROR: lon or lat is not numeric.')
#         return '0'

            
    if (grid_info is None):
        grid_info = get_hlr_grid_info()

    if (hlr_grid is None):
        new_dir = dtu.get_new_data_dir( 'USGS_HLR' )
        raster_dir = new_dir + 'Raster/'
        code_tif_file = grid_info['code_file']
        code_tif_path = raster_dir + code_tif_file
        hlr_grid = regrid.read_geotiff( in_file=code_tif_path, GEO=False )
        print('## hlr_grid shape =', hlr_grid.shape)
        print('## hlr_grid dtype =', hlr_grid.dtype)
        print('## hlr_grid min   =', hlr_grid.min())
        print('## hlr_grid max   =', hlr_grid.max())
        print()

    if (out_prj_file is None):
        new_dir = dtu.get_new_data_dir( 'USGS_HLR' )
        raster_dir = new_dir + 'Raster/'
        out_prj_file = 'hlrus.prj'
        out_prj_file = (raster_dir + out_prj_file)

    xmin = grid_info['xmin']
    xmax = grid_info['xmax']
    ymin = grid_info['ymin']
    ymax = grid_info['ymax']
    dx   = grid_info['dx']     # 1000.0 for HLR data set
    dy   = grid_info['dy']     # 1000.0 

    #----------------------------------------------
    # EPSG = 4326 is Geographic coords (WGS 84)
    # HLR data uses: Lambert_Azimuthal_Equal_Area
    #   with particular projection parameters.
    #-------------------------------------------------
    # EPSG = 2163 is a Lambert Azimuthal Equal Area
    # projection w/ almost identical info to our PRJ
    # https://spatialreference.org/ref/epsg/2163/
    # https://epsg.io/2163
    # Download PRJ file and compare.
    # Also very similar to EPSG = 9311, which is a
    #  "non-deprecated replacement" for EPSG 2163?
    #-------------------------------------------------
    # Bounds of conterminous US are:
    #     maxlat = 49.382808
    #     minlat = 24.521208
    #     maxlon = -66.945392
    #     minlon = -124.736342
    # Also see Notes in hydrofab_tools.py.
    #-------------------------------------------------
    if (USE_GDAL):
        #--------------------------------------------------------- 
        # Note: This also works, but only if we switch the order
        #       lon and lat args and set SWAP_XY=False and use
        #       out_proj_file. It fails if we keep lon,lat order
        #       and set SWAP_XY=True or False.
        #---------------------------------------------------------   
        lon_np = np.float64(lon)  # convert_coords needs this type
        lat_np = np.float64(lat)
        ## x,y = su.convert_coords(lon_np, lat_np,
        x,y = su.convert_coords(lat_np, lon_np,
                         inEPSG=4326,  inPRJfile=None,
                         ## outEPSG=9311, outPRJfile=None,
                         ## outEPSG=2163, outPRJfile=None,
                         outEPSG=None, outPRJfile=out_prj_file, 
                         PRINT=False, SWAP_XY=False)
    else:
        try:
            lon_np = np.float64(lon)  # convert_coords needs this type
            lat_np = np.float64(lat)
        except:
            print('### ERROR converting lon and lat to floats:')
            print('  lon =', lon)
            print('  lat =', lat)
        x,y = Lambert_Azimuthal_Equal_Area_XY( lon_np, lat_np )
        ## x,y = Lambert_Azimuthal_Equal_Area_XY( lon, lat )
        ### x,y = proj.Lambert_Azimuthal_Equal_Area_XY( lon, lat )

    #----------------------------------------
    # Check if x and y fall inside the grid
    #----------------------------------------
    in_xrange = (x > xmin) and (x < xmax)
    in_yrange = (y > ymin) and (y < ymax)
    if not( in_xrange and in_yrange ):
        if (WARNINGS):
            print('### ERROR: Computed x,y falls outside of grid.')
            print('###        Returning hlr_code = 0.')
            print('### lon =', lon)
            print('### lat =', lat)
            #--------------------------------------------
            # Note: These are Lambert_Azimuthal coords.
            #--------------------------------------------
#             if not( in_xrange ):
#                 print('### xmin, xmax =', xmin, ',', xmax)
#                 print('### x =', x)
#             if not (in_yrange ):
#                 print('### ymin, ymax =', ymin, ',', ymax)
#                 print('### y =', y) 
            print()
        return '0' 
            
    #------------------------------------------------
    # Recall that row 0 is top row & goes with ymax
    #------------------------------------------------
    col = np.int32((x - xmin) / dx)
    row = np.int32((ymax - y) / dy)
    hlr_code = hlr_grid[ row, col ] 

    if (REPORT):
        print('## HLR Information for Lon/Lat:')
        print('lon, lat =', lon, ',', lat)
        print('x =', x)
        print('y =', y)
        print('col  =', col)
        print('row  =', row)
        print('code =', hlr_code)
                       
    return str(hlr_code)   ### Return as string

#   get_hlr_code_for_point()
#---------------------------------------------------------------------   
# def point_in_rectangle(xP,yP, xA,yA, xB,yB, xC,yC, xD,yD):
def point_in_rectangle(xP,yP, xA,yA, xC,yC):

    #-------------------------------------------------------------
    # Note:  This is for a rectangle with arbitrary orientation,
    #        but in our case the top and bottom edges are
    #        parallel to the x-axis (east-west axis).
    #        This version allows the arguments to be arrays.
    #
    #        . (xA,yA)        . (xB,yB)
    #                                        . (xP,yP)
    #        . (xD,yD)        . (xC,yC)
    #    
    #        So in our case:
    #           xA = xD = minlon
    #           xB = xC = maxlon
    #           yC = yD = minlat
    #           yA = yB = maxlat
    #   
    # See:   https://martin-thoma.com/
    #        how-to-check-if-a-point-is-inside-a-rectangle/
    #-------------------------------------------------------------
    # Rectangle with arbitrary rotation or orientation.
    #----------------------------------------------------
#     ABCD = 0.5 * np.abs((yA - yC)*(xD - xB) + (yB - yD)*(xA - xC))
#     ABP  = 0.5 * np.abs(xA*(yB - yP) + xB*(yP - yA) + xP*(yA - yB))
#     BCP  = 0.5 * np.abs(xB*(yC - yP) + xC*(yP - yB) + xP*(yB - yC))
#     CDP  = 0.5 * np.abs(xC*(yD - yP) + xD*(yP - yC) + xP*(yC - yD))
#     DAP  = 0.5 * np.abs(xD*(yA - yP) + xA*(yP - yD) + xP*(yD - yA))
#     return np.invert(ABCD < (ABP + BCP + CDP + DAP))

    #---------------------------------------
    # Rectangle as geographic bounding box
    #-------------------------------------------------- 
    # This lets us eliminate xB,yB and xD,yD to get:
    #--------------------------------------------------  
#     ABCD = (xC - xA) * (yA - yC)
#     ABP  = 0.5 * np.abs(xA*(yA - yP) + xC*(yP - yA) )
#     BCP  = 0.5 * np.abs(xC*(yC - yP) + xC*(yP - yA) + xP*(yA - yC) )
#     CDP  = 0.5 * np.abs(xC*(yC - yP) + xA*(yP - yC) )
#     DAP  = 0.5 * np.abs(xA*(yA - yP) + xA*(yP - yC) + xP*(yC - yA) )
#     return np.invert(ABCD < (ABP + BCP + CDP + DAP))

    #-------------------------------------------------------
    # Rectangle as geographic bounding box: simpler method
    # Here again, xA,yA, xC,yC can be arrays.
    #-------------------------------------------------------
    in_xrange = np.logical_and( xP >= xA, xP <= xC )
    in_yrange = np.logical_and( yP >= yC, yP <= yA )
    return np.logical_and( in_xrange, in_yrange )

#   point_in_rectangle()
#--------------------------------------------------------------------- 
def save_hlr_code_grid_to_rtg( REPORT=True ):

    new_dir = dtu.get_new_data_dir( 'USGS_HLR' )
    raster_dir = new_dir + 'Raster/'
    
    hlr_grid_info     = get_hlr_grid_info()
    hlr_code_tif_file = hlr_grid_info['code_file']
    hlr_code_tif_path = raster_dir + hlr_code_tif_file
    rtg_path = hlr_code_tif_path.replace('.tif', '.rtg')

    #---------------------------------------     
    # Read HLR code grid from GeoTIFF file
    #---------------------------------------     
    hlr_code_grid = regrid.read_geotiff( in_file=hlr_code_tif_path, GEO=False)
    if (REPORT):
        print('Before conversion to BYTE type:')
        print('   min(hlr_code_grid) =', hlr_code_grid.min())
        print('   max(hlr_code_grid) =', hlr_code_grid.max())
    #--------------------------------------------------------    
    # This makes several incorrect assumptions re: RTI file
    #--------------------------------------------------------
#   hlr_code_grid = regrid.read_geotiff( in_file=hlr_code_tif_path, GEO=False,     
#                           MAKE_RTG=True, rtg_file=rtg_path )

    #------------------------------------------------
    # Convert grid to BYTE type (values in 1 to 20)
    # Use value "0" for anything out of range.
    #------------------------------------------------
    w = np.logical_or(hlr_code_grid > 20, hlr_code_grid < 1)
    hlr_code_grid[w] = 0 
    code_grid = np.int8( hlr_code_grid )
    data_type = 'BYTE'
    if (REPORT):
        print('After conversion to BYTE type:')
        print('   min(code_grid) =', code_grid.min())
        print('   max(code_grid) =', code_grid.max())
        print()

    #----------------------------------------
    # Extract grid info needed for RTI file
    #----------------------------------------
    ncols        = hlr_grid_info['ncols']
    nrows        = hlr_grid_info['nrows']
    xres         = hlr_grid_info['dx']
    yres         = hlr_grid_info['dy']
    zres         = hlr_grid_info['zres']
    x_west_edge  = hlr_grid_info['xmin']
    x_east_edge  = hlr_grid_info['xmax']
    y_south_edge = hlr_grid_info['ymin']
    y_north_edge = hlr_grid_info['ymax'] 
    pixel_geom   = hlr_grid_info['pixel_geom']
    z_units      = hlr_grid_info['z_units']
    xy_units     = hlr_grid_info['xy_units']
    data_source  = hlr_grid_info['source']
  
    #--------------------------------------
    # Create an RTI file for new RTG file
    #--------------------------------------
    rti = rti_files.make_info(grid_file=rtg_path,
              ncols=ncols, nrows=nrows,
              xres=xres, yres=yres,
              #--------------------------------------
              data_source=data_source,
              data_type=data_type,  ## Must agree with dtype  #######
              byte_order=rti_files.get_rti_byte_order(),
              pixel_geom=pixel_geom,
              zres=zres, z_units=z_units,
              y_south_edge=y_south_edge,
              y_north_edge=y_north_edge,
              x_west_edge=x_west_edge,
              x_east_edge=x_east_edge,
              box_units=xy_units, gmin=0, gmax=20)
    rti_files.write_info( rtg_path, rti )    
    rtg = rtg_files.rtg_file() 
    OK  = rtg.open_new_file( rtg_path, rti )
    if not(OK):
        print('ERROR during open_new_file().')
        return       
    rtg.write_grid( code_grid, VERBOSE=True )
    rtg.close_file()
              
#   save_hlr_code_grid_to_rtg()
#---------------------------------------------------------------------
def create_hlr_conus_grids( REPORT=True, MAKE_RTG=True ):

    #----------------------------------------------------------------
    # Note: All of the grids here have the same projection:
    #       Lambert Azimuthal Equal Area, with same parameters.
    #       See na_dem.prj.  They also have 1000 meter resolution.
    #----------------------------------------------------------------
    # This function clips larger North America (NA) grids to the
    # CONUS bounding box without changing the 1k spatial resolution
    # HYDRO1K.
    #----------------------------------------------------------------
        
    #----------------------------------
    # CONUS bounding box coordinates
    # (Lambert Azimuthal Equal Area)
    #----------------------------------
    # With xres = yres = 1000 meters,
    # ncols = 5020, nrows = 3230
    #----------------------------------
    xmin = -2285016.599999999627
    xmax = 2734983.400000000373
    ymin = -2260809.600000000093
    ymax = 969190.399999999907
    conus_bounds = [xmin, ymin, xmax, ymax]

    #-------------------------------------------
    # This NA DEM has 9102 cols and 8384 rows.
    #-------------------------------------------
    na_dem_dir = get_hlr_dem_dir( 'North_America' )    
    dem_file  = 'na_dem.tif'
    dem_path  = na_dem_dir + dem_file
    #----------------------------------------------
    new_dir = dtu.get_new_data_dir( 'USGS_HLR' )
    raster_dir = new_dir + 'Raster/'
    #----------------------------------------------
    code_file = 'HLR_code_grid_9106x6855_1k.tif'
    code_path = raster_dir + code_file

    #----------------------------------------------
    # Next 2 grids have 5020 cols and 3230 rows.
    #----------------------------------------------
    conus_dem_file  = 'HLR_DEM_CONUS_1k.tif'
    conus_dem_path  = raster_dir + conus_dem_file
    #----------------------------------------------
    conus_code_file = 'HLR_Code_CONUS_1k.tif'
    conus_code_path = raster_dir + conus_code_file
       
    #------------------------------ 
    # Create DEM clipped to CONUS
    #------------------------------
    regrid.regrid_geotiff(in_file=dem_path, out_file=conus_dem_path,
                   out_bounds=conus_bounds, GEO=False,
                   out_xres_m=None, out_yres_m=None,
                   RESAMPLE_ALGO='bilinear', REPORT=REPORT)

    #-----------------------------------------------------   
    # Also save new TIF file in RTG format with RTI file
    #-----------------------------------------------------
    if (MAKE_RTG):
        conus_dem_rtg_path = conus_dem_path.replace('.tif', '.rtg')               
        regrid.read_geotiff(in_file=conus_dem_path, REPORT=True,
                     MAKE_RTG=True, rtg_file=conus_dem_rtg_path, GEO=False,
                     zres=1, z_units='METERS', box_units='METERS',
                     data_source='HYDRO1K North America' )
                                    
    #----------------------------------------   
    # Create HLR Code grid clipped to CONUS
    #---------------------------------------- 
    regrid.regrid_geotiff(in_file=code_path, out_file=conus_code_path,
                   out_bounds=conus_bounds, GEO=False,
                   out_xres_m=None, out_yres_m=None,
                   RESAMPLE_ALGO='bilinear', REPORT=REPORT)
                   
    #-----------------------------------------------------   
    # Also save new TIF file in RTG format with RTI file
    #-----------------------------------------------------
    if (MAKE_RTG):
        conus_code_rtg_path = conus_code_path.replace('.tif', '.rtg')               
        regrid.read_geotiff(in_file=conus_code_path, REPORT=True,
                     MAKE_RTG=True, rtg_file=conus_code_rtg_path, GEO=False,
                     zres=1, z_units='NONE', box_units='METERS',
                     data_source='Wolock et al. (2004)' )
                 
    print('Finished.')
    print()
                      
#   create_hlr_conus_grids()
#---------------------------------------------------------------------
def get_hlr_outlet_info( SAVE_TO_FILE=True ):

    # Note: HLR outlet info can differ based on whether
    #       NWIS_ALL is True or not.
    ## usgs_dir  = get_usgs_dir( NWIS_ALL=True )
    ## info_path = usgs_dir + 'USGS_HLR_outlet_info.npy'
 
    new_dir   = dtu.get_new_data_dir( 'USGS_HLR' )
    info_path = new_dir + 'USGS_HLR_outlet_info.npy'
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
    hlr_path = dtu.get_new_tsv_filepath( 'USGS_HLR' )
    hlr_unit = open( hlr_path, 'r' )
    hlr_line = hlr_unit.readline()   # skip one-line header
    n_sites = dtu.get_n_records( 'USGS_HLR_outlets')
    ### n_sites = 9539
    
    #--------------------------------------    
    # Construct the hlr_outlet_info array
    #-----------------------------------------------
    # Numpy Unicode string array w/ up to 16 chars
    # Every element is initialized as null string
    # basin_id, lon, lat, hlr_code
    #-----------------------------------------------
    k = 0
    delim = '\t'  # tab character

    hlr_outlet_info = np.zeros((n_sites,4), dtype='<U16')
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
    
    if (SAVE_TO_FILE):
        np.save( info_path, hlr_outlet_info )
        print('Saved HLR outlet info to file:')
        print('  ' + info_path)
        print()
         
    return hlr_outlet_info
    
#   get_hlr_outlet_info()
#---------------------------------------------------------------------
def get_closest_hlr_outlet(usgs_id, usgs_site_info,
                           hlr_outlet_info, REPORT=False):

    #------------------------------------------------------------
    # This function uses usgs_site_info and hlr_outlet_info.
    # This function's caller should obtain these via functions:
    # get_usgs_site_info_dict() & get_hlr_outlet_info().
    #------------------------------------------------------------
    site_info = usgs.usgs_site_info[ usgs_id ]
    usgs_lon = site_info[ 'lon' ]
    usgs_lat = site_info[ 'lat' ]
    
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
def plot_color_legend( fig_xsize=8, fig_ysize=9,
               marker_size=65, font_size=24, SAVE_PNG=False ):

    #-----------------------------------------------------------
    # marker size for matplotlib.pyplot.plot is in points.
    # dot size (s) for matplotlib.pyplot.scatter is points^2.
    #-----------------------------------------------------------
    # If you increase fig xsize and ysize, you should also
    # increase marker_size and font_size.
    #-----------------------------------------------------------        
    start_time = time.time()   
    print('Running...')
   
    #-------------------------
    # Create a pyplot figure
    #-------------------------
    fig_config = {'figsize': (fig_xsize, fig_ysize)}
    fig = plt.figure(**fig_config)
    ax = plt.gca()
    ax.set_aspect('equal')  ########

    xmin = 0
    xmax = 1000
    ymin = 0
    ymax = 1000
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
      
    #--------------------------------------
    # Assign colors to the 20 HLR classes
    #--------------------------------------
    # Orange first
    rgb_color_map = {
    '20': (242, 0, 60),   '1':  (248, 89, 0),  '2':  (242, 136, 0),
    '3':  (242, 171, 0),  '4':  (239, 204, 0), '5':  (240, 234, 0),
    '6':  (177, 215, 0),  '7':  (0, 202, 36),  '8':  (0, 168, 119),
    '9':  (0, 167, 138),  '10': (0, 165, 156), '11': (0, 163, 172),
    '12': (0, 147, 175),  '13': (0, 130, 178), '14': (0, 110, 191),
    '15': (125, 0, 248),  '16': (159, 0, 197), '17': (185, 0, 166),
    '18': (208, 0, 129),  '19': (226, 0, 100) }
    #----------------------------------------------------------------
    # Red first
#     rgb_color_map = {
#     '1':  (242, 0, 60),   '2':  (248, 89, 0),  '3':  (242, 136, 0),
#     '4':  (242, 171, 0),  '5':  (239, 204, 0), '6':  (240, 234, 0),
#     '7':  (177, 215, 0),  '8':  (0, 202, 36),  '9':  (0, 168, 119),
#     '10': (0, 167, 138),  '11': (0, 165, 156), '12': (0, 163, 172),
#     '13': (0, 147, 175),  '14': (0, 130, 178), '15': (0, 110, 191),
#     '16': (125, 0, 248),  '17': (159, 0, 197), '18': (185, 0, 166),
#     '19': (208, 0, 129),  '20': (226, 0, 100) }

    #----------------------------------------------
    # For matplotlib, rgb colors must be in [0,1]
    #----------------------------------------------
    color_map = dict()
    for key in rgb_color_map.keys():
        color_map[ key ] = tuple(np.array( rgb_color_map[key] ) / 255)

    #----------------------------------------
    # If swb_g2 is an ndarray, you can use:
    #----------------------------------------------------------
    # For an even faster method see:
    # https://stackoverflow.com/questions/16992713/
    # translate-every-element-in-numpy-array-according-to-key
    #----------------------------------------------------------
    ## hlr_colors = np.vectorize(color_map.get)(hlr_codes)
    ## print('hlr_colors[0:19]', hlr_colors[0:19])

    #---------------------------------------
    # Add a legend. Can change x_L and y_T
    # and everything else will adjust.
    #---------------------------------------
    fs = font_size   # set by keyword
    ms = marker_size
    ## mkr = 'o'  # circle
    mkr = 's'  # square
    mec = 'black'  # mec = markeredgecolor for plot command
    mew = 0.5      # mew = markeredgewidth for plot command
    dx_text1 = -18
    dx_text2 = -33  # for 2-digit numbers
    dx = 180
    dy = 180
    #---------------------
    x1   = 140.0
    y_T  = 850.0
    #----------------------
    xt1  = x1 + dx_text1
    x2   = x1 + dx
    xt2  = x2 + dx_text1
    x3   = x2 + dx
    xt3a = x3 + dx_text1
    xt3  = x3 + dx_text2
    x4   = x3 + dx
    xt4  = x4 + dx_text2
    x5   = x4 + dx
    xt5  = x5 + dx_text2   
    #--------------------------------------------------------------------------
    plt.plot(x1,  y_T,    marker=mkr, markersize=ms, color=color_map['1'],
             mec=mec, mew=mew)
    plt.text(xt1, y_T,    '1', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x1,  y_T-dy, marker=mkr, markersize=ms, color=color_map['2'],
             mec=mec, mew=mew)           
    plt.text(xt1, y_T-dy, '2', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x1,  y_T-2*dy, marker=mkr, markersize=ms, color=color_map['3'],
             mec=mec, mew=mew)           
    plt.text(xt1, y_T-2*dy, '3', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x1,  y_T-3*dy, marker=mkr, markersize=ms, color=color_map['4'],
             mec=mec, mew=mew)           
    plt.text(xt1, y_T-3*dy, '4', fontsize=fs, ha='left', va='center')
    #----------------------
    # 2nd set of 4 colors
    #----------------------
    plt.plot(x2,  y_T,      marker=mkr, markersize=ms, color=color_map['5'],
             mec=mec, mew=mew)           
    plt.text(xt2, y_T,      '5', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x2,  y_T-dy,    marker=mkr, markersize=ms, color=color_map['6'],
             mec=mec, mew=mew)
    plt.text(xt2, y_T-dy,    '6', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x2,  y_T-2*dy, marker=mkr, markersize=ms, color=color_map['7'],
             mec=mec, mew=mew)           
    plt.text(xt2, y_T-2*dy, '7', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x2,  y_T-3*dy, marker=mkr, markersize=ms, color=color_map['8'],
             mec=mec, mew=mew)           
    plt.text(xt2, y_T-3*dy, '8', fontsize=fs, ha='left', va='center')
    #----------------------
    # 3rd set of 4 colors
    #----------------------   
    plt.plot(x3,  y_T,      marker=mkr, markersize=ms, color=color_map['9'],
             mec=mec, mew=mew)           
    plt.text(xt3a, y_T,     '9',  fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x3,  y_T-dy,   marker=mkr, markersize=ms, color=color_map['10'],
             mec=mec, mew=mew)           
    plt.text(xt3, y_T-dy,   '10', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x3,  y_T-2*dy, marker=mkr, markersize=ms, color=color_map['11'],
             mec=mec, mew=mew)
    plt.text(xt3, y_T-2*dy, '11', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x3,  y_T-3*dy, marker=mkr, markersize=ms, color=color_map['12'],
             mec=mec, mew=mew)           
    plt.text(xt3, y_T-3*dy, '12', fontsize=fs, ha='left', va='center')
    #----------------------
    # 4th set of 4 colors
    #----------------------
    plt.plot(x4,  y_T,      marker=mkr, markersize=ms, color=color_map['13'],
             mec=mec, mew=mew)           
    plt.text(xt4, y_T,      '13', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x4,  y_T-dy,   marker=mkr, markersize=ms, color=color_map['14'],
             mec=mec, mew=mew)           
    plt.text(xt4, y_T-dy,   '14', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x4,  y_T-2*dy, marker=mkr, markersize=ms, color=color_map['15'],
             mec=mec, mew=mew)           
    plt.text(xt4, y_T-2*dy, '15', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x4,  y_T-3*dy, marker=mkr, markersize=ms, color=color_map['16'],
             mec=mec, mew=mew)
    plt.text(xt4, y_T-3*dy, '16', fontsize=fs, ha='left', va='center')
    #----------------------
    # 5th set of 4 colors
    #----------------------
    plt.plot(x5,  y_T,    marker=mkr, markersize=ms, color=color_map['17'],
             mec=mec, mew=mew)           
    plt.text(xt5, y_T,    '17', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x5,  y_T-dy,   marker=mkr, markersize=ms, color=color_map['18'],
             mec=mec, mew=mew)           
    plt.text(xt5, y_T-dy,   '18', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x5,  y_T-2*dy, marker=mkr, markersize=ms, color=color_map['19'],
             mec=mec, mew=mew)           
    plt.text(xt5, y_T-2*dy, '19', fontsize=fs, ha='left', va='center')
    #--------------------------------------------------------------------------
    plt.plot(x5,  y_T-3*dy, marker=mkr, markersize=ms, color=color_map['20'],
             mec=mec, mew=mew)           
    plt.text(xt5, y_T-3*dy, '20', fontsize=fs, ha='left', va='center')
        
    #------------------------------------- 
    # Draw a rectangle around the legend
    #-------------------------------------
    width  = 940
    height = 760
    x0     = 30.0
    y0     = 200.0 
    rect = patches.Rectangle((x0,y0), width, height,
                   linewidth=1, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
       
    #--------------------------
    # Show the completed plot
    #--------------------------
    plt.show()

    #----------------------------    
    # Save figure to image file
    #----------------------------
    if (SAVE_PNG):
        fig.savefig(data_dir + 'SWB_GAGES2_CONUS.png')

    print('Finished.')
    
#   plot_color_legend()
#---------------------------------------------------------------------




