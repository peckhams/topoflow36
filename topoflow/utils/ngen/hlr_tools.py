#
# Copyright (c) 2023, Scott D. Peckham
#
# May 2023. Wrote: create_tsv,
#           get_polygon_points (all versions),
#           get_cols_and_rows, get_polygon_zvals,
#           check_lon_lats, get_zmin_lon_and_lat,
#           convert_coords)
# Jun 2023. Wrote convert_vat_adf_to_csv, get_bounding_box.
#           Moved many functions into shape_utils.py.
# Jul 2023. Modified to use newer shape_utils.py, which
#           contains many of the original functions from here.
#           Added get_hlr_data_dir, get_hlr_dem_dir.
#           Added get_repeated_value_list, to exclude them.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils import hlr_tools as hlr
#  >>> hlr.create_tsv(nf_max=100, tsv_file='new_hlr_na_100.tsv')
#  >>> hlr.create_tsv(nf_max=50000, tsv_file='new_hlr_na_all.tsv')
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
#    * 47479 features/basins (3648 more than 43391)
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
#  get_basin_repo_dir()    # Modify these as needed.
#  get_hlr_data_dir()
#  get_hlr_dem_dir()
#  get_repeated_values_list()
#  get_missing_basin_name()   # NOT READY YET
#
#  create_tsv()      # main function
#
#  get_cols_and_rows()
#  get_polygon_zvals()
#  check_lons_lats()
#  get_zmin_lon_and_lat()
#  get_zmin_lon_and_lat1()  # uses np.argmin()
#
#  convert_vat_adf_to_csv()
#
#---------------------------------------------------------------------
import numpy as np
from osgeo import ogr, osr
import json, time

from osgeo import gdal
from topoflow.utils import shape_utils as su

# For convert_vat_adf_to_csv()
import csv, sys, gdal

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
def get_hlr_data_dir():

    #-----------------------------------
    # Modify this directory as needed.
    #-----------------------------------
    repo_dir = get_basin_repo_dir()
    data_dir = repo_dir + 'HLR/'
    return data_dir
       
#   get_hlr_data_dir()
#---------------------------------------------------------------------
def get_hlr_dem_dir():

    #-----------------------------------
    # Modify this directory as needed.
    #-----------------------------------
    # repo_dir  = get_basin_repo_dir()
    dem_dir = '/Users/peckhams/DEMs/HYDRO1K/'
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
    #       Try to match lon/lat to a USGS station to get name.
    #------------------------------------------------------------
    pass
    
#   get_missing_basin_name()
#---------------------------------------------------------------------
def create_tsv( shape_dir=None, dem_dir=None,
                shape_file='hlrus4.shp', # shape_file='hlrus4.shp',
                ## prj_file='hlrus.prj',                     
                prj_file='na_dem.prj', dem_file='na_dem.bil',
                old_shapefile = False,  # old = original shapefile
                tsv_file='new_hlr_na.tsv',  # only North America
                #### ADD_BOUNDING_BOX=True, 
                nf_max=50, SWAP_XY=True, REPORT=True ):

    start_time = time.time()
    REPORT2 = True  # for basins outside North America
    # REPORT2 = REPORT
    print('Running...')

    if (shape_dir is None):
        shape_dir = get_hlr_data_dir() + 'hlrshape/'
    if (dem_dir is None):
        dem_dir  = get_hlr_dem_dir()

    ## prj_file = shape_file.replace('.shp', '.prj')
    
    shape_path = shape_dir + shape_file
    prj_path   = shape_dir + prj_file
    dem_path   = dem_dir   + dem_file
    tsv_path   = shape_dir + tsv_file
    shape_unit = su.open_shapefile( shape_path )
    layer      = shape_unit.GetLayer()

    #-------------------------------------------
    # HYDRO-1K DEM Info (Lambert Azimuthal EA)
    #-----------------------------------------------
    # The HYDRO-1K DEM used for the original HLR
    # work appears to have: ncols=9106, nrows=6855
    # based on the HLR metadata page.
    #-----------------------------------------------
    ncols = 9102
    nrows = 8384
    dx    = 1000.0  # meters
    dy    = 1000.0  # meters
    #-----------------------------
    # Use half grid cell offsets
    #-----------------------------
    xmin  = -4462500.0
    xmax  =  4639500.0
    ymin  = -3999500.0
    ymax  =  4384500.0
    #----------------------
    # From blw world file
    #----------------------
#     xmin  = -4462000.0
#     xmax  = xmin + (ncols * dx)
#     ymax  = 4384000.0
#     ymin  = ymax - (nrows * dy)
    #-----------------------------                 
    zmin  = -207
    zmax  = 6106    # Denali, AK
    nodata = -9999  # used for ocean grid cells

    dem_unit = open( dem_path, 'rb')
    tsv_unit = open( tsv_path, 'w') 

    n_features        = np.int32(0)
    zmin_matches      = np.int32(0)
    one_zmin_matches  = np.int32(0)
    all_good_matches  = np.int32(0)
    #--------------------------------
    two_try_count     = np.int32(0)
    empty_zvals_count = np.int32(0)
    code_zero_count   = np.int32(0)
    not_in_dem_count  = np.int32(0)
    one_cell_count    = np.int32(0)
    count_area_count  = np.int32(0)  # for old_shapefile only
   
    repeated_val_list = get_repeated_values_list()
    
    max_val = 48000 
    hlr_histogram   = np.zeros(21, dtype='int32')
    val_histogram   = np.zeros(max_val, dtype='int32')
    not_in_dem_vals = list()
    
    insert_key    = 'VALUE'
    new_headings  = ['OUTLON', 'OUTLAT', 'OUTCOL', 'OUTROW']
    new_headings += ['MINLON', 'MAXLON', 'MINLAT', 'MAXLAT']
    #--------------------------------------------------------------
#     new_headings = ['OUTLON', 'OUTLAT', 'OUTCOL', 'OUTROW']
#     if (ADD_BOUNDING_BOX):
#         new_headings += ['MINLON', 'MAXLON', 'MINLAT', 'MAXLAT']

    #-----------------------------------------    
    # Iterate over the features in the layer
    # GeometryType = 3 is POLYGON.
    #---------------------------------------------------------
    # Names in the shapefile attribute table are:
    # AREA, PERIMETER, HLRUS_COV_, HLRUS_COV1, VALUE, COUNT,
    # AQPERMNEW, SLOPE, TAVE, PPT, PET, SAND, PMPE, MINELE,
    # RELIEF, PFLATTOT, PFLATLOW, PFLATUP, HLR
    #---------------------------------------------------------
    for feature in layer:
        geometry   = feature.GetGeometryRef()
        attributes = feature.items()  # This is a Python dictionary
        # n_rings    = geometry.GetGeometryCount()

        #--------------------------------
        # Write header for new TSV file
        #--------------------------------
        if (n_features == 0):
            su.write_tsv_header( tsv_unit, attributes, insert_key=insert_key,
                  new_headings=new_headings )
              
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
        if (old_shapefile):
            hlr_area   = attributes['AREA']
            hlr_perim  = attributes['PERIMETER']
            if (hlr_area == 1e6 * hlr_count):
                count_area_matches += 1
            if (hlr_area == 1e6):
                one_cell_matches += 1
        hlr_value  = attributes['VALUE']  # watershed ID  
        hlr_count  = attributes['COUNT'] 
        hlr_zmin   = attributes['MINELE']
        hlr_relief = attributes['RELIEF']
        hlr_slope  = attributes['SLOPE']
        hlr_code   = attributes['HLR']
        if (hlr_code == 0):
            code_zero_matches += 1
        if (hlr_count == 1):
            # In new shapefile, (min,max) = (31, 23476)
            one_cell_matches += 1

        #----------------------------------------
        # Print attributes in attribute table ?
        #----------------------------------------
        if (REPORT):
            su.print_attributes( attributes, n_features )
 
        #----------------------------------------------
        # Skip 103 watersheds with "repeated VALUE"
        # that seem to have fragmented basin polygons
        #----------------------------------------------
        # Could return bounding box of the fragments
        # but then which outlet would we choose?
        #----------------------------------------------
        if (int(hlr_value) in repeated_val_list):
            print('Skipping repeated value =', hlr_value)
            continue
                          
        #---------------------------------- 
        # Get all points in this geometry
        #----------------------------------
        x1, y1 = su.get_polygon_points( feature )

        #----------------------------------
        # Round to nearest 1000 or 500 ??
        #----------------------------------------
        # Most x1 values end in 16.5 or 16.625,
        # and most y1 values end in 190.5.
        #----------------------------------------
#         x1 = np.around(x1 - 500, -3)
#         y1 = np.around(y1 - 500, -3)
#         x1 = np.around(x1, -3) + 500
#         y1 = np.around(y1, -3) - 500
#         x1 -= 16.5      # (old shapefile)  
#         y1 -= 190.5
#         x1 -= 16.6      # (new shapefile)  
#         y1 -= 190.4     
#         print('x1 =', x1)
#         print('y1 =', y1)
        # break

        #--------------------------------------------------        
        # Convert Lambert coords to Geo. lon/lat (WGS-84)
        #--------------------------------------------------
        x2,y2 = su.convert_coords(x1,y1, inPRJfile=prj_path,
                                  SWAP_XY=SWAP_XY, PRINT=False)
#         print('x2 =', x2)
#         print('y2 =', y2)
#         break

        #------------------------------ 
        # Get geographic bounding box
        #------------------------------
        minlon, maxlon, minlat, maxlat = su.get_bounding_box(x2, y2)
  
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
        cols, rows = get_cols_and_rows( x1, y1, xmin, ymax, dx, dy,
                         round_method='floor', print_c2_r2=False )
        ## break
         
        #----------------------------------------------
        # Need this because Hawaii is not in this DEM
        #----------------------------------------------
        lons, lats, cols, rows = check_lon_lats( cols, rows, ncols, nrows, x2, y2)       
        if (lons.size == 0):
            if (REPORT2):
                print('##### WARNING #####')
                print('Basin is not in the North America DEM.')
                print('Watershed ID (VALUE) =', hlr_value)
                print('It is likely in Hawaii.')
                print('Skipping to next basin.')
                print()
            not_in_dem_count += 1
            not_in_dem_vals.append( hlr_value )
            continue
               
        #---------------------------------------------------       
        # Get elevation (z) values for every polygon point
        # then get zmin and corresponding lon and lat.
        #---------------------------------------------------
        zvals = get_polygon_zvals( cols, rows, ncols, dem_unit,
                                   REPORT=REPORT )
        # This cannot happen now.
#         if (zvals.size == 0):
#             empty_zvals_count += 1
#             print('#### WARNING: zvals has size 0.')
#             print('#### VALUE =', hlr_value)
#             print()
#             continue
      
        zmin, zmin_arg, lon, lat, one_zmin = get_zmin_lon_and_lat( zvals,
              lons, lats, nodata, REPORT=REPORT)

        #--------------------------------------------------------        
        # If zmin doesn't match hlr_zmin, then try again with
        # round_method = 'ceil'.  Still don't get all matches.
        # For old_shapefile,
        # Got 168 matches out of first 200 polygons (16% fails)
        # Also got 20/20, 28/30, 42/50.
        #-------------------------------------------------------- 
        if (zmin != hlr_zmin):
            two_try_count += 1
            cols, rows = get_cols_and_rows( x1, y1, xmin, ymax, dx, dy,
                             round_method='ceil', print_c2_r2=False )
            lons, lats, cols, rows = check_lon_lats( cols, rows, ncols, nrows, x2, y2)       
            if (lons.size == 0):
                if (REPORT2):
                    print('##### WARNING #####')
                    print('Basin is not in the North America DEM.')
                    print('Watershed ID (VALUE) =', hlr_value)
                    print('It is likely in Hawaii.')
                    print('Skipping to next basin.')
                    print()
                not_in_dem_count += 1  # can't double count due to continue
                not_in_dem_vals.append( hlr_value )
                continue
            zvals = get_polygon_zvals( cols, rows, ncols, dem_unit,
                                       REPORT=REPORT )
            if (zvals.size == 0):
                empty_zvals_count += 1
                print('#### WARNING: zvals has size 0.')
                print('#### VALUE =', hlr_value)
                print()
                continue      
            zmin, zmin_arg, lon, lat, one_zmin = get_zmin_lon_and_lat( zvals,
                  lons, lats, nodata, REPORT=REPORT )
   
        #-----------------------------------        
        # Count when zmin matches hlr_zmin
        #-----------------------------------
        if (zmin == hlr_zmin):
            zmin_matches += 1
            if (one_zmin):
                one_zmin_matches += 1
        #------------------------------
        col = cols[ zmin_arg ]
        row = rows[ zmin_arg ]
           
        if (REPORT):
            print()
            print('####### Basin number:', n_features, '########')
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

        #--------------------------------
        # How many "all good matches" ?
        #--------------------------------
        if (old_shapefile):
            ALL_GOOD = ((hlr_area == 1e6 * hlr_count) and
                (hlr_code > 0) and (hlr_area > 1e6) and
                (zmin == hlr_zmin) and (one_zmin))
        else:
            ALL_GOOD = ((hlr_code > 0) and (hlr_count > 1) and
                (zmin == hlr_zmin) and (one_zmin))
 
        #---------------------------------------------               
        # Save info for good matches to new TSV file
        #---------------------------------------------
        if (ALL_GOOD):
            all_good_matches += 1

            su.write_tsv_line( tsv_unit, attributes, insert_key=insert_key,
                    new_values=[lon,lat,col,row,minlon,maxlon,minlat,maxlat])
            
            # How many basins with each HLR code?
            hlr_histogram[ hlr_code ]  += 1
            val_histogram[ hlr_value ] += 1

        n_features += 1
        if (n_features == nf_max):
            break
            
    #--------------------------------
    # Print info, close files, etc.
    #--------------------------------
    print('# features          =', n_features)
    print('# zmin matches      =', zmin_matches, 'of', nf_max)
    print('# one zmin matches  =', one_zmin_matches)
    print('# all good matches  =', all_good_matches, '(written to TSV)')
    #-------------------------------------------------------
    print('# two try count     =', two_try_count)
    print('# empty zvals count =', empty_zvals_count)
    print('# code zero count   =', code_zero_count)
    print('# not in dem count  =', not_in_dem_count)
    print('# one cell count    =', one_cell_count)
    if (old_shapefile):
        print('# count-area count  =', count_area_count)
        
    print()
    print('hlr_histogram =')
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
def get_cols_and_rows( x1, y1, xmin, ymax, dx, dy,
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

    #------------------------------------------
    # Note: First and last zvals are the same
    #------------------------------------------
    if (REPORT):
        print('zvals =', zvals)  #######
        print('n_pts =', zvals.size - 1)
        print()
    return zvals
 
#   get_polygon_zvals()
#---------------------------------------------------------------------
def check_lon_lats( cols, rows, ncols, nrows, x2, y2):

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

    #--------------------------------------------------
    # This version does not use np.argmin(), since it
    # only returns index of the first min value.
    #--------------------------------------------------
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
    
    #----------------------------------- 
    # Round to 4 places after decimal;
    # roughly 11 meters.  Can only do
    # in-place w/out for ndarrays.
    #-----------------------------------
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
    hlr_dir = get_hlr_data_dir()
    vat_dir = hlr_dir + 'arctar00000/hlrus/'
    
    if (vat_adf_file is None):
        vat_adf_file = vat_dir + 'vat.adf'
    if (csv_file is None):
        csv_file = vat_dir + 'vat_adf.csv'

    #------------------------------------------
    # Note: ogr is a vector library;
    #       gdal is a raster library
    #------------------------------------------
    ds  = gdal.Open( vat_adf_file )
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



