
# Copyright (c) 2023, Scott D. Peckham
#
# May 2023. Wrote: check_shapefile,
#           get_polygon_points (all versions),
#           get_cols_and_rows, get_polygon_zvals,
#           check_lon_lats, get_zmin_lon_and_lat,
#           convert_coords)
# Jun 2023. Wrote convert_vat_adf_to_csv, get_bounding_box.
#
# NOTE: Many code comments are for original (old) shapefile,
#       that has numerous issues. The "new" shapefile
#       comes from QGIS work with the raster version
#       of HLR dataset in arctar00000 folder.
#       For some reason, the same watershed ID "VALUE"
#       occurs in multiple rows, but with different outlet.
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils import hlr_tools as hlr
#  >>> hlr.check_shapefile()
#
#---------------------------------------------------------------------
#  check_shapefile()      # main function
#  get_polygon_points()
#  get_polygon_points1()  # also works
#  get_polygon_points2()  # doesn't work
#  get_polygon_points3()  # doesn't work
#  get_cols_and_rows()
#  get_polygon_zvals()
#  check_lons_lats()
#  get_zmin_lon_and_lat()
#  get_zmin_lon_and_lat1()  # uses np.argmin()
#  convert_coords()
#  get_bounding_box()
#
#  convert_vat_adf_to_csv()
#
#---------------------------------------------------------------------
import numpy as np
from osgeo import ogr, osr
import json, time

# For convert_vat_adf_to_csv()
import csv, sys
from osgeo import gdal

#---------------------------------------------------------------------
def check_shapefile( data_dir=None, dem_dir=None,
                     shape_file='hlrus4.shp',
                     ## shape_file='hlrus.shp',
                     ## prj_file='hlrus.prj',                     
                     prj_file='na_dem.prj',
                     dem_file='na_dem.bil',
                     old_shapefile = False,  # old = original shapefile
                     csv_file='new_hlr_na.csv',  # only North America
                     nf_max=50, REPORT=True ):

    start_time = time.time()
    REPORT2 = True  # for basins outside North America
    # REPORT2 = REPORT
    print('Running...')
    #--------------------------------------------
    # You will need to modify these directories
    # and get all the necessary files.
    #--------------------------------------------
    if (data_dir is None):
        data_dir  = '/Users/peckhams/Dropbox/NOAA_NextGen/'
        data_dir += '__NextGen_Example_Basin_Repo/'
        data_dir += 'HLR_Data_Sets/hlrshape/'
    if (dem_dir is None):
        dem_dir  = '/Users/peckhams/DEMs/HYDRO1K/'
    shape_path = data_dir + shape_file
    prj_path   = data_dir + prj_file
    dem_path   = dem_dir  + dem_file
    csv_path   = data_dir + csv_file 
    shape_unit = ogr.Open( shape_path )
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
    csv_unit = open( csv_path, 'w') 

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
   
    max_val = 48000 
    hlr_histogram   = np.zeros(21, dtype='int32')
    val_histogram   = np.zeros(max_val, dtype='int32')
    not_in_dem_vals = list()
    
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
        # Write header for new CSV file
        #--------------------------------
        if (n_features == 0):
            csv_header = ''
            for key in attributes.keys():
                csv_header += (key + ',')
                if (key == 'VALUE'):
                    csv_header += 'OUTLON,OUTLAT,'
                    csv_header += 'OUTCOL,OUTROW,'
                    csv_header += 'MINLON,MAXLON,MINLAT,MAXLAT,'
            csv_header = csv_header[:-1] + '\n' # remove comma, add newline
            csv_unit.write( csv_header )
              
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
        #   Claim is 43,931 basins, but there are
        #   47,579 basins in shapefile (3648 more).
        #--------------------------------------------
        #   659 basins have invalid HLR code of 0.
        #   These basins have COUNT=0, VALUE=-9999
        #--------------------------------------------
        #   Claim is pruning threshold of 100 km2,
        #   but a large number are much smaller. 
        #   3391 basins have area < 1e7 m2 = 10 km2
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
                
        #---------------------------------- 
        # Get all points in this geometry
        #----------------------------------
        x1, y1 = get_polygon_points( feature )
        # x1, y1 = get_polygon_points1( geometry ) 
        # x1, y1 = get_polygon_points2( feature )
        # x1, y1 = get_polygon_points3( feature )
        # p = geometry.GetPoints()  # doesn't work
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
        x2,y2 = convert_coords(x1,y1, inPRJfile=prj_path,
                               PRINT=False)
#         print('x2 =', x2)
#         print('y2 =', y2)
#         break

        #------------------------------ 
        # Get geographic bounding box
        #------------------------------
        minlon, maxlon, minlat, maxlat = get_bounding_box(x2, y2)
  
        #---------------------------------------------------------
        # The source DEM is the GTOPO30 HYDRO1K DEM and uses
        # a Lambert Azimuthal EqA projection as described above.
        # Can index this 1-km DEM with the Lambert xy coords
        # to get an elevation, then compare to MINELE attribute
        # as a way to get the outlet xy.  Or just find the xy
        # pair for each polygon with the lowest elevation.
        # But the DEM will be pretty big.  The HYDRO1K DEM for
        # North America includes Alaska, CONUS, PR, etc. but not
        # Hawaii. It has 9102 cols and 8384 rows.  Need to
        # convert Lambert xy coords to DEM column and row.
        #---------------------------------------------------------
        # See: http://www.ncgia.ucsb.edu/SANTA_FE_CD-ROM/
        #       sf_papers/verdin_kristine/santafe2.html
        # HYDRO1K DEM:  https://d9-wret.s3.us-west-2.amazonaws.com/
        #    assets/palladium/production/s3fs-public/atoms/files/
        #    HYDRO1kDocumentationReadMe.pdf
        #---------------------------------------------------------
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
            print('hlr_value  =', hlr_value)
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
        # Save info for good matches to new CSV file
        #---------------------------------------------
        if (ALL_GOOD):
            all_good_matches += 1
            n_vals = 0
            line = ''
            for val in attributes.values():
                line += str(val) + ','   ###### CHECK FOR LOSS
                n_vals += 1
                if (n_vals == 1):
                    line += str(lon) + ','
                    line += str(lat) + ','
                    line += str(col) + ','
                    line += str(row) + ','
                    #------------------------
                    line += str(minlon) + ','
                    line += str(maxlon) + ','
                    line += str(minlat) + ','
                    line += str(maxlat) + ','
                    n_vals += 8
            line = line[:-1] + '\n' # remove comma, add newline
            csv_unit.write( line )
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
    print('# all good matches  =', all_good_matches)
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
    print('# Repeated watershed IDs (VALUE) =', w.sum())
    print('vals =')
    print( i[w] )
    print('hist =')
    print( val_histogram[w] )
    print('not_in_dem_vals =')
    print( not_in_dem_vals )
    print()

    dem_unit.close()
    csv_unit.close()
    run_time = (time.time() - start_time)
    print('Run time =', run_time, '[sec]')
    print('Finished.')
    
#   check_shapefile()
#---------------------------------------------------------------------
def get_polygon_points( feature ):

    #---------------------------------------------------- 
    # Get all points in this geometry by exporting
    # feature to GeoJSON and converting json to dict().
    #----------------------------------------------------
    # Printing jstr shows that triple square brackets
    # are used for coordinates, so need subscript [0].
    # This is likely because polygons can have more
    # than one "ring" for holes (in general).
    #----------------------------------------------------
    jstr = feature.ExportToJson()
    # print(jstr)
    jdict = json.loads( jstr ) # python dictionary
    points = np.array( jdict['geometry']['coordinates'][0] )
    x1 = points[:,0]
    y1 = points[:,1]
    return x1, y1
    
#   get_polygon_points()
#---------------------------------------------------------------------
def get_polygon_points1( geometry ):

    #--------------------------------------------------
    # Note: About same speed as get_polygon_points().
    #--------------------------------------------------
    ring = geometry.GetGeometryRef(0)  # outer ring
    p  = ring.GetPoints()    # This is a list of xy tuples
    p2 = np.array(p)
    x1 = p2[:,0]
    y1 = p2[:,1]
    return x1, y1

#   get_polygon_points1()
#---------------------------------------------------------------------
# def get_polygon_points2( feature ):
# 
#     #----------------------------------
#     # Get all points in this geometry
#     # This version doesn't work.
#     #----------------------------------
#     n_points = geometry.GetPointCount()
#     if (n_points == 0):
#         print('###### WARNING ######')
#         print('No points in feature.')
#         continue
#     x1 = np.zeros( n_points, dtype='float64')
#     y1 = np.zeros( n_points, dtype='float64')
#     for i in range(n_points):
#         p = geometry.GetPoint(i)
#         x1[i] = p[0]
#         y1[i] = p[1]
#     return x1, y1
    
#   get_polygon_points2()
#---------------------------------------------------------------------
# def get_polygon_points3( feature ):
# 
#     #----------------------------------
#     # Get all points in this geometry 
#     #--------------------------------------------
#     # GetBoundary sort of works, but deprecated
#     #--------------------------------------------      
#     b = geometry.GetBoundary()  # deprecated
#     p = b.GetPoints()    # This is a list of xy tuples
#     #-------------------------------
#     # print( type(p) )   # list
#     # print( type(p[0])) # tuple
#     # print(p)
#     # print(p[0])
#     #------------------------------------------------
#     p2 = np.array(p)
#     if (p2.ndim < 2):
#         print('#### WARNING FOR FEATURE ####')
#         print('p  =', p)
#         print('p2 =', p2)
#         continue
#     x1 = p2[:,0]
#     y1 = p2[:,1]
# 
#     return x1, y1
    
#   get_polygon_points3()
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
def convert_coords(x1, y1, inEPSG=None, inPRJfile=None,
                   outEPSG=4326, PRINT=False):

    #------------------------------------------------------------
    # Note: ESPG=4326 is for WGS-84 (Geographic lon/lat)
    #------------------------------------------------------------
    # Read WKT (Well Known Text) from HLR shapefile's PRJ file.
    # Projection metadata is also provided in:
    # https://water.usgs.gov/GIS/metadata/usgswrd/XML/hlrus.xml
    #------------------------------------------------------------
    # Geographic bounding box for HLR DEM is given by:
    #    (minlon, maxlon) = (-180, -60)
    #    (minlat, maxlat) = (25, 75)
    # This includes Alaska and Hawaii; not just CONUS.
    # Can check in Google Earth via View > Grid option.
    #------------------------------------------------------------
    if (inPRJfile is not None):
        prj_unit = open(inPRJfile, 'r')
        hlr_wkt = prj_unit.read()
        prj_unit.close()

    #-----------------------------------
    # Create coordinate transformation
    #-----------------------------------
    srs_in = osr.SpatialReference()
    if (inEPSG is None):
        srs_in.ImportFromWkt( hlr_wkt )
    else:
        srs_in.ImportFromEPSG( inESPG )
    
    srs_out = osr.SpatialReference()
    srs_out.ImportFromEPSG( outEPSG )

    coordTransform = osr.CoordinateTransformation(srs_in, srs_out)

    #------------------------------------------
    # Is (x1,y1) a point or a basin polygon ?
    #------------------------------------------
    if (x1.size == 1):
        #----------------------
        # Transform the point
        #----------------------
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x1, y1)
        point.Transform( coordTransform )
        x2 = point.GetY()
        y2 = point.GetX()  #### Need to swap.
    else:
        #----------------------------------
        # Transform each point in polygon
        #----------------------------------
        x2 = np.zeros(x1.size, dtype='float64')
        y2 = np.zeros(x1.size, dtype='float64')
        for k in range(x1.size): 
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(x1[k], y1[k])
            point.Transform( coordTransform )
            x2[k] = point.GetY()
            y2[k] = point.GetX()  #### Need to swap.
            
        #----------------------------------
        # Transform each point in polygon
        #----------------------------------
#         ring = ogr.Geometry(ogr.wkbLinearRing)
#         for k in range(x1.size): 
#             point = ogr.Geometry(ogr.wkbPoint)
#             point.AddPoint(x1[k], y1[k])
#             point.Transform( coordTransform )
#             x2k = point.GetY()
#             y2k = point.GetX()  #### Need to swap.
#             x2[k] = x2k
#             y2[k] = y2k    
#             ring.AddPoint(x2k, y2k)
#         p  = ring.GetPoints()  # list of tuples
#         p2 = np.array(p)
#         x2 = p2[:,0]
#         y2 = p2[:,1]
        # x2 = ring.GetX()  # only gets x2[0]
        #------------------------------------------   
#         basin = ogr.Geometry(ogr.wkbPolygon)
#         basin.AddGeometry(ring)
#         basin.CloseRings()
#         p  = basin.GetPoints() # DOESN'T WORK
#       # This doesn't work
#       basin.Transform( coordTransform )
#       x2 = basin.GetX()
#       y2 = basin.GetY()
              
    #----------------------------------
    # Print input and output points ?
    #----------------------------------
    if (PRINT):
        print('Input:  x1, y1 =', x1, ',', y1 )
        print('Output: x2, y2 =', x2, ',', y2 )  # lon,lat
        print()

    return x2, y2

#   convert_coords()
#---------------------------------------------------------------------
def get_bounding_box(x2, y2, buffer=0.02):

    #---------------------------------------------------------
    # Note: Buffer is in decimal degrees, and 0.02 is
    #       roughly 2.2 km.  See Wikipedia: Decimal degrees.
    #       Recall that HYDRO1K DEM has 1 km grid cells.
    #---------------------------------------------------------
    minlon = x2.min() - buffer
    maxlon = x2.max() + buffer
    minlat = y2.min() - buffer
    maxlat = y2.max() + buffer

    #----------------------------------- 
    # Round to 4 places after decimal;
    # roughly 11 meters.  Can only do
    # in-place w/out for ndarrays.
    #-----------------------------------
    minlon = np.around(minlon, decimals=4)
    maxlon = np.around(maxlon, decimals=4)
    minlat = np.around(minlat, decimals=4)
    maxlat = np.around(maxlat, decimals=4)
 
    return minlon, maxlon, minlat, maxlat
    
#   get_bounding_box()
#---------------------------------------------------------------------
def convert_vat_adf_to_csv( vat_adf_file=None, csv_file=None):

    #-------------------------------------------------------
    # Note: This function was adapted from:
    # https://gis.stackexchange.com/questions/84700/
    # showing-only-raster-attribute-table-using-gdal/84729
    # See "__README_RAT_TO_CSV.txt." about error message.
    #-------------------------------------------------------
    data_dir  = '/Users/peckhams/Dropbox/NOAA_NextGen/'
    data_dir += '__NextGen_Example_Basin_Repo/HLR_Data_Sets/'
    data_dir += 'arctar00000/hlrus/'
    
    if (vat_adf_file is None):
        vat_adf_file = data_dir + 'vat.adf'
    if (csv_file is None):
        csv_file = data_dir + 'vat_adf.csv'

    ds  = gdal.Open(vat_adf_file)
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
    ds = None  # (Close vat.adf file)

#   convert_vat_adf_to_csv() 
#---------------------------------------------------------------------   



