#---------------------------------------------------------------------
# Copyright (c) 2024, Scott D. Peckham

# May 2024. Wrote first version.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import subset as sub
#  >>> sub.read_collated_tsv_file()
#  >>>
#
#---------------------------------------------------------------------
#
#  read_collated_tsv_file()
#  get_basin_subset()
#
#---------------------------------------------------------------------

from topoflow.utils.ngen import data_utils as dtu
import numpy as np

#---------------------------------------------------------------------
def read_collated_tsv_file():

    repo_dir = dtu.get_repo_dir( OLD=False )
    tsv_file = repo_dir + '__Collated/collated_basins_all.tsv'
    tsv_unit = open( tsv_file, 'r' )
    delim    = '\t'  # tab char
    
    #------------------------------------
    # Get the column headings as a list
    #------------------------------------
    line = tsv_unit.readline()
    line = line[:-1]  # remove newline char
    headings = line.split( delim )
    # print('headings =')
    # print(headings)
    # print()

    #----------------------------------    
    # Get indices for certain columns
    #----------------------------------
    is_sb3_col = headings.index('Is_GAGES2_SB3')
    swb_col    = headings.index('SWB_Class')
    hlr_col    = headings.index('HLR_Code_Outlet')
    area_col   = headings.index('Area')

    #--------------------------------------------    
    # Prepare to count basins in each SWB class
    #--------------------------------------------
    swb_classes   = ['A1','A2','A3','B1','B2','B3',
        'C1', 'C2', 'D1', 'D2', 'D3']
    swb_counts = dict()
    for item in swb_classes:
        swb_counts[ item ] = 0
    #--------------------------------------------
    swb_area_mins  = dict()
    swb_area_maxes = dict()
    for item in swb_classes:
        swb_area_mins[ item ]  = 100000000.0
        swb_area_maxes[ item ] = 0.0
        
    #--------------------------------------------    
    # Prepare to count basins in each HLR class
    #--------------------------------------------
    hlr_counts = dict()
    for j in range(21):
        hlr_counts[ j ] = 0
    #--------------------------------------------
    n_basins = 0
    n_missing_areas = 0
          
    while (True):
        line = tsv_unit.readline()
        if (line == ''):
            break   # (reached end of file)

        line = line[:-1]  # remove newline char
        values = line.split( delim )
        is_sb3_val = values[ is_sb3_col ]
        if (is_sb3_val == 'Y'):
            swb_class = values[ swb_col ]
            swb_counts[ swb_class ] += 1
            #---------------------------------
            area_val = values[ area_col ].strip()
            ## print('area =', area_val)
            if ((area_val != '') and (area_val != '-9999')):
                area = np.float32( values[ area_col ] )
                old_min = swb_area_mins[ swb_class ]
                old_max = swb_area_maxes[ swb_class ]
                swb_area_mins[ swb_class ]  = min(old_min, area)
                swb_area_maxes[ swb_class ] = max(old_max, area)
            else:
                n_missing_areas += 1
            #---------------------------------
            hlr_code  = values[ hlr_col ]
            hlr_counts[ int(hlr_code) ] += 1
            #---------------------------------
            n_basins += 1

    #-------------------------            
    # Print some information
    #-------------------------
    for item in swb_classes:
        amin = swb_area_mins[ item ]
        amax = swb_area_maxes[ item ]
        msg  = 'swb_count(' + item + ') =' + str(swb_counts[item])
        msg += ', amin = ' + str(amin) + ', amax = ' + str(amax)
        print(msg)
    print()
    for j in range(21):
        print('hlr_count(' + str(j) + ') =', hlr_counts[j])
    print()
    print('n_basins =', n_basins)
    print('n_missing_areas =', n_missing_areas)
    print()

    #------------------------------------  
    # Close the collated basin TSV file
    #------------------------------------
    tsv_unit.close()
            
#   read_collated_tsv_file()
#---------------------------------------------------------------------
def get_basin_subset():

    pass
    
#   get_basin_subset()
#---------------------------------------------------------------------