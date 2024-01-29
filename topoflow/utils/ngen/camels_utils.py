
#  Copyright (c) 2023-2024, Scott D. Peckham
#
#  Jan 2024. Renamed all "ngen/utils" files to end in "_utils.py"
#            instead of "_tools.py".
#            Two versions of create_camels_bbox_file to compare.
#            Wrote get_camels_shapefile() and cleaned up.
#            Modified to use new data_utils.py.
#  Jul 2023. Original version. 
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import camels_utils as cu
#  >>> cu.create_camels_bbox_file2()
#  >>> cu.create_tsv()
#
#---------------------------------------------------------------------
#
#  get_camels_shapefile()
#  create_camels_bbox_file1()  # uses full-res shapefile
#  create_camels_bbox_file2()  # uses shapefile with HRUs
#  sort_camels_bbox_file()
#  find_extra_ids()
#  create_tsv()
#
#---------------------------------------------------------------------

import numpy as np
from osgeo import ogr, osr
# import json, sys

import time
from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import shape_utils as su

#---------------------------------------------------------------------
def get_camels_shapefile( ftype='full_res' ):

    camels_dir = dtu.get_data_dir( 'CAMELS' )
    
    if (ftype == 'full_res'):
        #------------------------------------------------
        # This contains a high-res version of the basin
        # boundary for each of the 671 CAMELS basins.
        #------------------------------------------------
        camels_dir += 'basin_set_full_res/'
        camels_file = camels_dir + 'HCDN_nhru_final_671.shp'
    elif (ftype == 'with_hrus'):
        #------------------------------------------------------
        # A given CAMELS basin may contain multiple HRUs.
        # This is a merged set of simplified basin boundaries
        # that include boundaries of 7844 individual HRUs.
        #------------------------------------------------------
        camels_dir += 'basin_timeseries_v1p2_metForcing_obsFlow/'
        camels_dir += 'basin_dataset_public_v1p2/'
        camels_dir += 'shapefiles/merge/'
        camels_file = camels_dir + 'basinset_gf_nhru.shp'
       
    return camels_file

#   get_camels_shapefile()
#---------------------------------------------------------------------
def create_camels_bbox_file1( tmp_bbox_file='camels_bbox.tmp',
                              bbox_file='camels_bbox.txt',
                              nb_max=700):
   
    #-------------------------------------------------------------
    # This function creates a semicolon-delimited text file,
    # similar to other camels files (e.g. camels_topo.txt),
    # that contains gauge_id, minlon, maxlon, minlat, & maxlat.
    # It reads from a full-res shapefile that has the basin
    # boundaries for all 671 CAMELS basins.
    # This is much faster than the similar function:
    #     create_camels_bbox_file2().
    #-------------------------------------------------------------
    # Since all CAMELS basins are included in GAGES2_SB3 and
    # GAGES2_ref, we could compare the bounding boxes.
    #-------------------------------------------------------------
    # Set nb_max to smaller value (e.g. 20) for testing.
    #-------------------------------------------------------------    
    shape_path  = get_camels_shapefile( 'full_res' )
    prj_path    = shape_path.replace('.shp', '.prj') 
    #---------------------------------------------------
    out_dir   = dtu.get_new_data_dir( 'CAMELS' )     
    bbox_unit = open( out_dir + tmp_bbox_file, 'w' )
 
    start_time = time.time()   
    shape_unit = ogr.Open( shape_path )
    layer      = shape_unit.GetLayer()
    n_basins   = np.int32(0)
    SWAP_XY    = True     #######
    delim      = ';'  # semi-colon delimiter
    minlon = 180.0
    maxlon = -180.0
    minlat = 90.0
    maxlat = -90.0

    #---------------------------------
    # Write header for new bbox_file
    #---------------------------------
    ## new_headings=['gauge_id','MINLON','MAXLON','MINLAT','MAXLAT']
    bbox_unit.write('gauge_id' + delim)
    bbox_unit.write('MINLON'   + delim)
    bbox_unit.write('MAXLON'   + delim)
    bbox_unit.write('MINLAT'   + delim)
    bbox_unit.write('MAXLAT'   + '\n')
  
    #-----------------------------------------    
    # Iterate over the features in the layer
    # GeometryType = 3 is POLYGON.
    #-----------------------------------------
    print('Working...')
    for feature in layer:
        geometry   = feature.GetGeometryRef()
        attributes = feature.items()  # This is a Python dictionary
        # n_rings    = geometry.GetGeometryCount()

        #-----------------------------        
        # Get the USGS site/gauge ID
        #-----------------------------
        gauge_id = str( attributes['hru_id'] ).strip()    #########
        if (len(gauge_id) == 7):
            gauge_id = '0' + gauge_id       
        
        #---------------------------------- 
        # Get all points in this geometry
        #----------------------------------
        x1, y1 = su.get_polygon_points( feature )
#         print('x1 =', x1)
#         print('y1 =', y1)
        # break

        #------------------------------------------       
        # Convert coords to Geo. lon/lat (WGS-84)
        #------------------------------------------
        x2,y2 = su.convert_coords(x1,y1, inPRJfile=prj_path,
                                  SWAP_XY=SWAP_XY, PRINT=False)
#         print('x2 =', x2)
#         print('y2 =', y2)
#         break        

        #---------------------------------------------------- 
        # Get geographic bounding box for this CAMELS basin
        #----------------------------------------------------
        minlon, maxlon, minlat, maxlat = su.get_bounding_box(x2, y2)
#         print()
#         print('minlon, maxlon =', minlon, maxlon)
#         print('minlat, maxlat =', minlat, maxlat)
        
        #-------------------------------------              
        # Write line of info to new CSV file
        #---------------------------------------------------------
        # su.get_bounding_box() rounds to 5 digits after decimal
        #---------------------------------------------------------   
        bbox_unit.write( gauge_id    + delim )
        bbox_unit.write( str(minlon) + delim )
        bbox_unit.write( str(maxlon) + delim )
        bbox_unit.write( str(minlat) + delim )
        bbox_unit.write( str(maxlat) + '\n' )
        n_basins += 1
        print('n_basins =', n_basins, 'of 671.')
        print('   gauge_id =', gauge_id)
        if (n_basins >= nb_max):
            break     

    #----------------------------------------
    # Close bbox file, then reopen and sort
    #----------------------------------------
    bbox_unit.close()
    sort_camels_bbox_file( in_bbox_file='camels_bbox.tmp',
                           out_bbox_file='camels_bbox.txt')
             
    #-------------------
    # Print final info
    #-------------------
    print('Processed', n_basins, 'CAMELS basins.') 
    run_time = (time.time() - start_time)
    print('Run time =', run_time, '[sec]')
    print('Finished.')
                           
#   create_camels_bbox_file1()
#---------------------------------------------------------------------
def create_camels_bbox_file2( tmp_bbox_file='camels_bbox.tmp',
                              bbox_file='camels_bbox.txt',
                              nb_max=700):

    #-------------------------------------------------------------
    # This function creates a semicolon-delimited text file,
    # similar to other camels files (e.g. camels_topo.txt),
    # that contains gauge_id, minlon, maxlon, minlat, & maxlat.
    # It reads from a merged shapefile that has the basin
    # boundaries for all HRUs within any given CAMELS basin.
    # However, it finds 677 vs. 671 basins
    #-------------------------------------------------------------
    # Set nb_max to smaller value (e.g. 20) for testing.
    #-------------------------------------------------------------    
    shape_path  = get_camels_shapefile( 'with_hrus' )
    prj_path    = shape_path.replace('.shp', '.prj') 
    #---------------------------------------------------
    out_dir   = dtu.get_new_data_dir( 'CAMELS')     
    bbox_unit = open( out_dir + tmp_bbox_file, 'w' )
 
    start_time = time.time()   
    shape_unit = ogr.Open( shape_path )
    layer      = shape_unit.GetLayer()
    n_hrus     = np.int32(0)
    n_basins   = np.int32(0)
    n_hrus_in_basin = np.int32(0)
    SWAP_XY    = True     #######
    delim      = ';'  # semi-colon delimiter
    last_gauge_id = '00000000'
    minlon = 180.0
    maxlon = -180.0
    minlat = 90.0
    maxlat = -90.0

    #---------------------------------
    # Write header for new bbox_file
    #---------------------------------
    ## new_headings=['gauge_id','MINLON','MAXLON','MINLAT','MAXLAT']
    bbox_unit.write('gauge_id' + delim)
    bbox_unit.write('MINLON'   + delim)
    bbox_unit.write('MAXLON'   + delim)
    bbox_unit.write('MINLAT'   + delim)
    bbox_unit.write('MAXLAT'   + '\n')
  
    #-----------------------------------------    
    # Iterate over the features in the layer
    # GeometryType = 3 is POLYGON.
    #-----------------------------------------
    print('Working...')
    for feature in layer:
        geometry   = feature.GetGeometryRef()
        attributes = feature.items()  # This is a Python dictionary
        # n_rings    = geometry.GetGeometryCount()

        gauge_id = attributes['GAGEID']           
        FIRST_ONE = (last_gauge_id == '00000000')
        # FIRST_ONE = (n_hrus == 0)
        
        #---------------------------------- 
        # Get all points in this geometry
        #----------------------------------
        x1, y1 = su.get_polygon_points( feature )
#         print('x1 =', x1)
#         print('y1 =', y1)
        # break

        #------------------------------------------       
        # Convert coords to Geo. lon/lat (WGS-84)
        #------------------------------------------
        x2,y2 = su.convert_coords(x1,y1, inPRJfile=prj_path,
                                  SWAP_XY=SWAP_XY, PRINT=False)
#         print('x2 =', x2)
#         print('y2 =', y2)
#         break        

        #-------------------------------------- 
        # Get geographic bounding box for one
        # of the HRUs in this CAMELS basin
        #--------------------------------------
        hru_minlon, hru_maxlon, hru_minlat, hru_maxlat = su.get_bounding_box(x2, y2)
#         print()
#         print('minlon, maxlon =', minlon, maxlon)
#         print('minlat, maxlat =', minlat, maxlat)
        n_hrus += 1
       
        if ((gauge_id == last_gauge_id) or FIRST_ONE):
            minlon = min(hru_minlon, minlon)
            maxlon = max(hru_maxlon, maxlon)
            minlat = min(hru_minlat, minlat)
            maxlat = max(hru_maxlat, maxlat)
        else:         
            #-------------------------------------              
            # Write line of info to new CSV file
            #---------------------------------------------------------
            # su.get_bounding_box() rounds to 5 digits after decimal
            #---------------------------------------------------------   
            bbox_unit.write( last_gauge_id + delim )
            bbox_unit.write( str(minlon)   + delim )
            bbox_unit.write( str(maxlon)   + delim )
            bbox_unit.write( str(minlat)   + delim )
            bbox_unit.write( str(maxlat)   + '\n' )
            n_basins += 1
            print('n_basins =', n_basins, 'of 671.')
            print('   gauge_id =', gauge_id)
            print('   n_hrus_in_basin =', n_hrus_in_basin)
            n_hrus_in_basin = np.int32(0)
            if (n_basins >= nb_max):
                break
            #-------------------------------------------
            # Use first HRU to initialize bounding box
            #-------------------------------------------
            minlon = hru_minlon
            maxlon = hru_maxlon
            minlat = hru_minlat
            maxlat = hru_maxlat

        n_hrus_in_basin += 1
        last_gauge_id = gauge_id      

    #----------------------------------------
    # Close bbox file, then reopen and sort
    #----------------------------------------
    bbox_unit.close()
    sort_camels_bbox_file( in_bbox_file='camels_bbox.tmp',
                           out_bbox_file='camels_bbox.txt')
             
    #-------------------
    # Print final info
    #-------------------
    print('Processed', n_hrus, 'HRUs.')
    print('Processed', n_basins, 'CAMELS basins.') 
    run_time = (time.time() - start_time)
    print('Run time =', run_time, '[sec]')
    print('Finished.')
   
#   create_camels_bbox_file2()
#---------------------------------------------------------------------
def sort_camels_bbox_file( in_bbox_file='camels_bbox.tmp',
                           out_bbox_file='camels_bbox.txt'):

    print('Sorting lines in new bbox file...')
    print('Reading lines from:', in_bbox_file)
    out_dir       = dtu.get_new_data_dir( 'CAMELS')
    in_bbox_unit  = open( out_dir + in_bbox_file,  'r' )      
    out_bbox_unit = open( out_dir + out_bbox_file, 'w' )
    
    lines    = in_bbox_unit.readlines()
    headings = lines[0]
    lines    = lines[1:]
    in_bbox_unit.close()
    print('Total number of lines =', len(lines)+1, '(incl. header line)' )

    #--------------------------------------     
    # Write sorted lines to new bbox file
    #--------------------------------------
    print('Writing sorted lines to:', out_bbox_file) 
    lines.sort()
    out_bbox_unit.write(headings)
    out_bbox_unit.writelines( lines )
    out_bbox_unit.close()
    print('Finished sorting by ID.')
    print()
 
#   sort_camels_bbox_file()
#---------------------------------------------------------------------
def find_extra_ids( CHECK_DUPES=False ):

    #--------------------------------------------------------------
    # Note: The bbox file created with create_camels_bbox_file2
    #       has 678 lines, while camels_name.txt and the others
    #       only have 672 lines (incl. 1 header line).  This
    #       function helps to find the extra gauge IDs.
    #--------------------------------------------------------------
    # Note: In the shapefile attribute table, the gauge ID:
    #       07067000 does occur first in a set of 77 HRUs that
    #       have hru_id_loc in {1,...,77}, and then in another
    #       set of 10 HRUs that have hru_id_loc in {1,...,10}.
    #       Some have gtype="polygon" and some have gtype=
    #       "multipolygon".
    #--------------------------------------------------------------
    # Results of this function:
    #
    # gauge IDs in bbox file not in name file:
    # ['06775500', '01150900', '06846500',
    # '09535100', '02081113', '03448942']
    #
    # gauge IDs in name file not in bbox file:
    # ['01333000']
    #
    # 01150900;-73.00165;-72.33164;43.64477;44.12633
    # 0208111310;-89.60984;-89.24013;30.44601;31.10538
    # 0344894205;-87.88285;-87.21651;35.18618;35.58351
    # 06775500;-105.93089;-105.83394;40.49394;40.56312
    # 06846500;-96.48799;-96.06662;39.86638;40.34738
    # 09535100;-112.1861;-111.8549;33.91465;34.11396
    #
    # 01333000;02;GREEN RIVER AT WILLIAMSTOWN, MA
    #--------------------------------------------------------------
    # Duplicate ID in bbox file (which one is right?):
    # One has 77 hrus, but they don't get highlighted in
    # QGIS with View > Zoom to Selection.  Other has 10 hrus.
    #
    # 07067000;-92.02911;-91.3148;36.94448;37.25004 (wider)
    # 07067000;-94.73446;-94.31482;32.74102;33.23056 (taller)
    #
    # 07067000;11;Current River at Van Buren, MO
    # See USGS website for this river:
    # https://waterdata.usgs.gov/monitoring-location/
    #    07067000/#parameterCode=00065&period=P7D
    #--------------------------------------------------------------
    # In QGIS, install the BoundingBox plugin (how to use?).
    # And/or right click on a point & select Copy Coordinate.
    #--------------------------------------------------------------
    dir1      = dtu.get_new_data_dir( 'CAMELS' )  
    dir2      = dtu.get_data_dir( 'CAMELS' )
    bbox_file = dir1 + 'camels_bbox.txt'
    name_file = dir2 + 'camels_name.txt'

    bbox_unit = open( bbox_file, 'r' )    
    lines1 = bbox_unit.readlines()
    bbox_unit.close()
    print('Number of IDs in bbox_file =', len(lines1)-1)
    print('Number of unique IDs in bbox_file =', len(set(lines1))-1)
    print()
    
    name_unit = open( name_file, 'r' )    
    lines2 = name_unit.readlines()
    name_unit.close()
    print('Number of IDs in name_file =', len(lines2)-1)
    print('Number of unique IDs in name_file =', len(set(lines2))-1)

    #--------------------------------    
    # Get USGS IDs in the bbox_file
    #--------------------------------        
    k = 0
    id_list1 = list()
    for line in lines1:
        id = lines1[k][:8]
        id_list1.append( id )
        k += 1
        # print('id =', id)

    #--------------------------------    
    # Get USGS IDs in the name_file
    #--------------------------------  
    k = 0
    id_list2 = list()
    for line in lines2:
        id = lines2[k][:8] 
        id_list2.append( id  )
        k += 1

    #----------------------------------
    # Find duplicate IDs in bbox_file
    #----------------------------------
    if (CHECK_DUPES):
        # This requires that id_list1 is sorted.
        print('Checking for duplicates...')
        last_id = '00000000'
        for id in id_list1:
            if (id == last_id):
                print('duplicate id =', id)
            last_id = id

#     k = 0
#     for id1 in id_list1[:-1]:
#         id2 = id_list2[k]
#         if (id2 != id1):
#             print('id1, id2 =', id1, id2)
#         k += 1

    #---------------------------
    # Compare the two ID lists
    #---------------------------
    print()   
    diff1 = list( set(id_list1) - set(id_list2))
    print('gauge IDs in bbox file not in name file:')
    print(diff1)
    print()

    diff2 = list( set(id_list2) - set(id_list1))
    print('gauge IDs in name file not in bbox file:')
    print(diff2)
    print()
   
    print('Finished.')
      
#   find_extra_ids()
#---------------------------------------------------------------------
def create_tsv( new_tsv_file = 'new_camels.tsv', delim='\t' ):

    #-------------------------------------------------------
    # Note: The file "camels_bbox.txt" was not part of the
    #       original CAMELS dataset, so it is in the same
    #       directory as the new TSV file to be created.
    #-------------------------------------------------------
    out_dir = dtu.get_new_data_dir( 'CAMELS' )
    in_dir  = dtu.get_data_dir( 'CAMELS' )

    #---------------------------------------------------
    # These CAMELS data files use ";" as the delimeter
    # All files start with "gauge_id" and are sorted
    #   by gauge_id.  All gauge IDs are 8 chars.
    #---------------------------------------------------
    name_file   = 'camels_name.txt'  # all in camels2_dir
    bbox_file   = 'camels_bbox.txt'  # must have 672 lines
    topo_file   = 'camels_topo.txt'  
    clim_file   = 'camels_clim.txt'
    soil_file   = 'camels_soil.txt'
    geol_file   = 'camels_geol.txt'

    #--------------------------------------------------    
    # Open text files to read & new csv file to write
    #--------------------------------------------------  
    name_unit = open( in_dir + name_file, 'r' )
    topo_unit = open( in_dir + topo_file, 'r' )
    clim_unit = open( in_dir + clim_file, 'r' )
    # soil_unit = open( in_dir + soil_file, 'r' )
    # geol_unit = open( in_dir + geol_file, 'r' )
    #---------------------------------------------------
    bbox_unit = open( out_dir + bbox_file, 'r' )   # SPECIAL
    tsv_unit  = open( out_dir + new_tsv_file, 'w' )

    #---------------------------    
    # Read the column headings
    #---------------------------
    name_line = name_unit.readline()
    name_head = name_line.replace(';', delim)
    name_head = name_head.replace('\n', '')
    #-------------------------------------------------    
    bbox_line = bbox_unit.readline()
    bbox_head = bbox_line.replace(';', delim)
    bbox_head = bbox_head[9:]   # exclude "gauge_id"
    bbox_head = bbox_head.replace('\n', '')
    #-------------------------------------------------    
    topo_line = topo_unit.readline()
    topo_head = topo_line.replace(';', delim)
    topo_head = topo_head[9:]   # exclude "gauge_id"
    topo_head = topo_head.replace('\n', '')
    #-------------------------------------------------     
    clim_line = clim_unit.readline()
    clim_head = clim_line.replace(';', delim)
    clim_head = clim_head[9:]   # exclude "gauge_id"
    clim_head = clim_head.replace('\n', '')
        
    #---------------------------------------- 
    # Write column headings to new tsv file
    #----------------------------------------    
    tsv_unit.write( name_head + delim )
    tsv_unit.write( bbox_head + delim )
    tsv_unit.write( topo_head + delim )
    tsv_unit.write( clim_head + '\n')

    #----------------------------- 
    # Write data to new tsv file
    #--------------------------------------------
    # Note: All gauge IDs are 8 chars.
    # Note: "gauge_id" heading is also 8 chars.
    # Note: All files are sorted by gauge ID.
    # Note: "readline" result ends w/ newline.
    #--------------------------------------------
    print('Working...')
    while (True):
        name_line = name_unit.readline()
        if (name_line == ''):
            break  # (reached end of file)
        name_line = name_line.replace(';', delim)
        tsv_unit.write( name_line[:-1] + delim )
        #-------------------------------------------        
        bbox_line = bbox_unit.readline()
        bbox_line = bbox_line.replace(';', delim)
        tsv_unit.write( bbox_line[9:-1] + delim )
        #-------------------------------------------        
        topo_line = topo_unit.readline()
        topo_line = topo_line.replace(';', delim)
        tsv_unit.write( topo_line[9:-1] + delim )
        #-------------------------------------------
        clim_line = clim_unit.readline()
        clim_line = clim_line.replace(';', delim)
        tsv_unit.write( clim_line[9:-1] + '\n')    # notice newline
        #-------------------------------------------
#         soil_line = soil_unit.readline()
#         soil_line = soil_line.replace(';', delim)
#         tsv_unit.write( soil_line[9:-1] + delim)        
#         #-------------------------------------------
#         geol_line = geol_unit.readline()
#         geol_line = geol_line.replace(';', delim)
#         tsv_unit.write( geol_line[9:-1] + '\n')
                       
    #------------------    
    # Close all files
    #------------------
    print('Finished.')
    name_unit.close()
    bbox_unit.close()
    topo_unit.close()
    clim_unit.close()
    # soil_unit.close()
    # geol_unit.close()
    #--------------------
    tsv_unit.close()
             
#   create_tsv()
#---------------------------------------------------------------------



