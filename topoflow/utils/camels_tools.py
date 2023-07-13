
# Copyright (c) 2023, Scott D. Peckham
#
# July 2023. 
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils import camels_tools as ct
#  >>> ct.create_camels_bbox_file()
#  >>> ct.create_tsv()
#
#---------------------------------------------------------------------
#
#  get_camels_data_dir()
#  create_camels_bbox_file()
#  sort_camels_bbox_file()
#  find_extra_ids()
#  create_tsv()
#
#---------------------------------------------------------------------

import numpy as np
from osgeo import ogr, osr
# import json, sys

import time
from topoflow.utils import shape_utils as su

#---------------------------------------------------------------------
def get_camels_data_dir(version='1.2'):

    if (version == '1.2'):
       camels_dir = '/Users/peckhams/Data/CAMELS/2016_v1.2/'
    elif (version == '2.0'):
       camels_dir = '/Users/peckhams/Data/CAMELS/2023_v2.0/'
    else:
       camels_dir  = '/Users/peckhams/Dropbox/NOAA_NextGen/'
       camesl_dir += '__NextGen_Example_Basin_Repo/CAMELS/'
       
    return camels_dir

#   get_camels_data_dir()
#---------------------------------------------------------------------
def create_camels_bbox_file( bbox_file='camels_bbox.txt',
                             nb_max=700):

    #-------------------------------------------------------------
    # This function creates a semicolon-delimited text file,
    # similar to other camels files (e.g. camels_topo.txt),
    # that contains gauge_id, minlon, maxlon, minlat, & maxlat.
    # It reads from a merged shapefile that has the basin
    # boundaries for all HRUs within any given CAMELS basin.
    #-------------------------------------------------------------
    # Change these directories as needed.
    # Set nb_max to smaller value (e.g. 20) for testing.
    #-------------------------------------------------------------    
    camels1_dir = get_camels_data_dir(version='1.2')
    camels2_dir = get_camels_data_dir(version='2.0')
    #------------------------------------------------------------- 
    shape_dir   = camels1_dir + 'shapefiles/merge/'
    shape_path  = shape_dir + 'basinset_gf_nhru.shp'
    prj_path    = shape_dir + 'basinset_gf_nhru.prj'
    #-------------------------------------------------------------     
    bbox_unit = open( camels2_dir + bbox_file, 'w' )
 
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
            bbox_unit.write( gauge_id    + delim )
            bbox_unit.write( str(minlon) + delim )
            bbox_unit.write( str(maxlon) + delim )
            bbox_unit.write( str(minlat) + delim )
            bbox_unit.write( str(maxlat) + '\n' )
            n_basins += 1
            print('n_basins =', n_basins, 'of 671.')
            print('   gauge_id =', gauge_id)
            print('   n_hrus_in_basin =', n_hrus_in_basin)
            n_hrus_in_basin = np.int32(0)
            if (n_basins >= nb_max):
                break
            #----------------------------
            # Reset bounding box coords
            #----------------------------
            minlon = hru_minlon
            maxlon = hru_maxlon
            minlat = hru_minlat
            maxlat = hru_maxlat

        n_hrus_in_basin += 1
        last_gauge_id = gauge_id      

    #--------------------------------
    # Print info, close files, etc.
    #--------------------------------
    print('Processed', n_hrus, 'HRUs.')
    print('Processed', n_basins, 'CAMELS basins.') 
    bbox_unit.close()
    run_time = (time.time() - start_time)
    print('Run time =', run_time, '[sec]')
    print('Finished.')
   
#   create_camels_bbox_file()
#---------------------------------------------------------------------
def sort_camels_bbox_file( bbox_file='camels_bbox0.txt',
         bbox_file2='camels_bbox.txt'):

    print('Reading lines from:', bbox_file)
    camels2_dir = get_camels_data_dir(version='2.0')  
    bbox_unit = open( camels2_dir + bbox_file, 'r' )
    lines    = bbox_unit.readlines()
    headings = lines[0]
    lines    = lines[1:]
    bbox_unit.close()
    print('len(lines) =', len(lines))

    #--------------------------------------     
    # Write sorted lines to new bbox file
    #--------------------------------------
    print('Writing sorted lines to:', bbox_file2) 
    lines.sort()
    bbox_unit2 = open( camels2_dir + bbox_file2, 'w' )
    bbox_unit2.write(headings)
    bbox_unit2.writelines( lines )
    bbox_unit2.close()
    print('Finished.')
    print()
 
#   sort_camels_bbox_file()
#---------------------------------------------------------------------
def find_extra_ids( CHECK_DUPES=False ):

    #--------------------------------------------------------------
    # Note: camels_bbox.txt has 678 lines, while camels_name.txt
    #       and the others only have 672 lines.  This function
    #       helps to find the extra gauge IDs.
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
    camels2_dir = get_camels_data_dir(version='2.0')
    # bbox_file = camels2_dir + 'camels_bbox_677.txt'
    bbox_file = camels2_dir + 'camels_bbox_672.txt'
    name_file = camels2_dir + 'camels_name.txt'

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
            
    k = 0
    id_list1 = list()
    for line in lines1:
        id = lines1[k][:8]
        id_list1.append( id )
        k += 1
        # print('id =', id)
 
    k = 0
    id_list2 = list()
    for line in lines2:
        id = lines2[k][:8] 
        id_list2.append( id  )
        k += 1

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

    #--------------------------------------------------------------
    # Note: CAMELS version 1.2 was released in 2016 and has much
    #       more data including basin shapefiles. (91 GB)
    #       CAMELS version 2.0 was released in July 2023 and
    #       does not include basin shapefiles. (7.5 GB)
    #--------------------------------------------------------------
    #       Change these directories as needed.
    #--------------------------------------------------------------
    camels1_dir = get_camels_data_dir(version='1.2')
    camels2_dir = get_camels_data_dir(version='2.0')
    tsv_dir     = get_camels_data_dir(version='repo')

    #---------------------------------------------------
    # These CAMELS data files use ";" as the delimeter
    # All files start with "gauge_id" and are sorted
    #   by gauge_id.  All gauge IDs are 8 chars.
    #---------------------------------------------------
    name_file   = 'camels_name.txt'  # all in camels2_dir
    bbox_file   = 'camels_bbox.txt'
    topo_file   = 'camels_topo.txt'  
    clim_file   = 'camels_clim.txt'
    soil_file   = 'camels_soil.txt'
    geol_file   = 'camels_geol.txt'

    #--------------------------------------------------    
    # Open text files to read & new csv file to write
    #--------------------------------------------------
    name_unit = open( camels2_dir + name_file, 'r' )
    bbox_unit = open( camels2_dir + bbox_file, 'r' )
    topo_unit = open( camels2_dir + topo_file, 'r' )
    clim_unit = open( camels2_dir + clim_file, 'r' )
    # soil_unit = open( camels2_dir + soil_file, 'r' )
    # geol_unit = open( camels2_dir + geol_file, 'r' )
    #---------------------------------------------------
    tsv_unit = open( tsv_dir + new_tsv_file, 'w' )

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

    #-----------------------------------------   
    # Don't forget about minlon,maxlon, etc.      ##############
    # Add afterwards w/ separate function?
    #----------------------------------------- 

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



