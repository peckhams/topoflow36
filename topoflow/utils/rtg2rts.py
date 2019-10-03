#!/usr/bin/env python

# Copyright (c) 2009, Scott D. Peckham

#-----------------------------------------------------
# Example of use at a Unix prompt:
#
#    % ./rtg2rts.py pc_sb_nger_ PCTest2.rts 160 241
#-----------------------------------------------------

import glob
import os.path
import sys

import numpy
from . import rti_files  # (see /data/progs/csdms/python)
    
#-----------------------------------------------------------------------
def rtg2rts( RTG_prefix, new_RTS_file,
             nx=None, ny=None, dtype='float32' ):
     
    #--------------------------------------------------------------
    # Notes:  This function "bundles" a set of binary grid files
    #         (generic RTG format) with the same file_name prefix
    #         and embedded time index into a single RTS file.
    #--------------------------------------------------------------
    # Note:   The machine byte order on beach is "little" or LSB.
    #         However, the "bin" files that Mark wants to merge
    #         into a single file are apparently "big" or MSB.
    #         In cases like this, change SWAP_BYTES to True.
    #--------------------------------------------------------------    
    SWAP_BYTES = True   #################
    
    #--------------------
    # Open new RTS file
    #--------------------
    RTS_EXISTS = os.path.exists( new_RTS_file )
    if (RTS_EXISTS):
        print('SORRY, An RTS file with the name')
        print(new_RTS_file)
        print('already exists.')
        return

    RTG_list = glob.glob( RTG_prefix + '*' )
    RTG_list = numpy.sort( RTG_list )
    print('Number of RTG files to merge =', len(RTG_list))
    
    #--------------------------------
    # Try to get info from RTI file
    #--------------------------------
    if (nx == None) and (ny == None):
        RTG_file1 = RTG_list[0]
        RTI_file  = rti_files.try_to_find_rti_file( RTG_file1 )
        if (RTI_file != 'none'):
            info = rti_files.read_info( RTI_file )
            nx   = info.ncols
            ny   = info.nrows
        else:
            print('ERROR: Could not find RTI file and nx and ny not provided.')
            return
    else:
        #--------------------------------------------
        # For case when used as a script under Unix
        #--------------------------------------------
        nx = eval(nx)
        ny = eval(ny)
        
    RTS_unit = open( new_RTS_file, 'wb' )

    #------------------------------
    # Write RTG grids to RTS file
    #------------------------------    
    for RTG_file in RTG_list:
        print('Reading values from:', RTG_file)
        RTG_unit = open( RTG_file, 'rb')
        grid     = numpy.fromfile( RTG_unit, count=nx*ny, dtype=dtype )
        if (SWAP_BYTES):  grid.byteswap()
        RTG_unit.close()
        grid.tofile( RTS_unit )

    #---------------------
    # Close new RTS file
    #---------------------
    RTS_unit.close()
    print('Finished writing RTG files to RTS format.')
    print(' ')
    
#   rtg2rts()
#-----------------------------------------------------------------------
if (__name__ == "__main__"):
    #-----------------------------------------------------
    # Note: First arg in sys.argv is the command itself.
    #-----------------------------------------------------
    n_args = len(sys.argv)
    if (n_args < 3):
        print('ERROR: This utility requires an RTG_prefix')
        print('and a new RTS filename as arguments.')
        print('sys.argv =', sys.argv)
        print(' ')
    elif (n_args == 3):
        rtg2rts( sys.argv[1], sys.argv[2] )
    elif (n_args == 5):
        rtg2rts( sys.argv[1], sys.argv[2], nx=sys.argv[3], ny=sys.argv[4] )
    elif (n_args == 6):
        rtg2rts( sys.argv[1], sys.argv[2], nx=sys.argv[3], ny=sys.argv[4],
                 dtype=sys.argv[5] )
    else:
        print('ERROR: Invalid number of arguments.')
        
