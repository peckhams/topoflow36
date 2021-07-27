
import glob, os
import numpy as np

#-----------------------------------------------------------------------
def check_npz_files():

    ## os.chdir(
    npz_file_list = sorted( glob.glob( '*.npz' ) )
    n_files = len( npz_file_list )

    #---------------------------------    
    # Check "contents" of first file
    #---------------------------------
    data = np.load( npz_file_list[0] )
    print('contents of file0 =', data.files )
    var_name = data.files[0]
    
    for npz_file in npz_file_list:
        data = np.load( npz_file )
        vals = data[ var_name ]
        print('vmin, vmax =', vals.min(), ', ', vals.max() )

#   check_npz_files()
#-----------------------------------------------------------------------