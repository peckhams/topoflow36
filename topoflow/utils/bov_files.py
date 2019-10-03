

## Copyright (c) 2009, Scott D. Peckham
## October 2009

#-------------------------------------------------------------------

import sys
from . import rti_files

#-------------------------------------------------------------------
#  Functions:
#
#  write_info_as_bov()
#  write_bov_file()
#
#-------------------------------------------------------------------
def write_info_as_bov( file_name, info, var_name='X', time=0.0 ):

    #-----------------------------------------------------
    # Notes: BOV byte order is LITTLE or BIG.
    #        BOV data types are FLOAT (and what else ??)
    #-----------------------------------------------------
    prefix   = rti_files.get_file_prefix( file_name )
    bov_file = (prefix + '.bov')
    
    write_bov_file( bov_file, file_name, var_name, time,
                    info.ncols, info.nrows, nz=1,
                    data_type=info.data_type,
                    byte_order=sys.byteorder.upper(),
                    x_west_edge=info.x_west_edge,
                    y_south_edge=info.y_south_edge )
            
#   write_info_as_bov()   
#-------------------------------------------------------------------
def write_bov_file(bov_file, dat_file, var_name, time,
                   nx, ny, nz=1, centering="zonal",
                   data_type='FLOAT', byte_order='LITTLE',
                   x_west_edge=0.0, y_south_edge=0.0,
                   SILENT=True):

    bov_unit   = open(bov_file, 'w')
    dim_str    = str(nx) + ' ' + str(ny) + ' ' + str(nz)
    origin_str = str(x_west_edge) + ' ' + str(y_south_edge)
    origin_str = origin_str + ' 0.0'
    time_str   = str(time)
    
    bov_unit.write('TIME:          ' + time_str   + "\n")
    bov_unit.write('DATA_FILE:     ' + dat_file   + "\n")
    bov_unit.write('DATA_SIZE:     ' + dim_str    + "\n")
    bov_unit.write('DATA_FORMAT:   ' + data_type  + "\n")
    bov_unit.write('VARIABLE:      ' + var_name   + "\n")
    bov_unit.write('DATA_ENDIAN:   ' + byte_order + "\n")
    bov_unit.write('CENTERING:     ' + centering  + "\n")
    bov_unit.write('BRICK_ORIGIN:  ' + origin_str + "\n")
    bov_unit.write('BRICK_SIZE:    ' + dim_str    + "\n")

    bov_unit.close()

    if not(SILENT):
        print('BOV file written to: ')
        print(bov_file)
        print(' ')
        
#   write_bov_file()
#-------------------------------------------------------------------

