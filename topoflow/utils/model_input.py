
## Copyright (c) 2009-2021, Scott D. Peckham
## January 9, 2009
## April 29, 2009
## May 2010, Changed var_types from {0,1,2,3} to
#            {'Scalar', 'Time_Series', 'Grid'}, etc.
#-------------------------------------------------------------------

#  open_file()
#  read_next2()
#  read_next()
#  read_scalar()
#  read_grid()
#  close_file()

#-------------------------------------------------------------------
import numpy as np
import os.path
import netCDF4 as nc

#-------------------------------------------------------------------
def open_file(var_type, input_file, NGEN_CSV=False):

    #-----------------------------------------------------
    # Note:  This method's name cannot be "open" because
    #        it calls Python's built-in "open()" method.
    #-----------------------------------------------------
    # print 'var_type   =', var_type
    # print 'input_file =', input_file
    
    #--------------------------------------------
    # A scalar input value was provided already
    #--------------------------------------------
    file_unit = None
    if (var_type.lower() == 'scalar'):
        return file_unit
    if (input_file == ''):
        print('ERROR in model_input.open_file():')
        print('    Input file is null string.')
        # print '    variable type =' + var_type
        return file_unit

    #----------------------------------
    # Does input file exist locally ?
    #----------------------------------
    if not(os.path.exists(input_file)):
        print('ERROR in model_input.open_file():')
        print('    Could not find input file =')
        print('    ' + input_file)
        # print '   ' + input_file
        return file_unit
            
    if (var_type.lower() == 'time_series'):
        #-----------------------------------------
        # Input file contains a time series and
        # is ASCII text with one value per line.
        #-----------------------------------------
        file_unit = open(input_file, 'r')
        if (NGEN_CSV):
            #----------------------
            # Skip one header row
            #----------------------
            line = file_unit.readline()
    else:
        file_ext = os.path.splitext(input_file)[1]
        if file_ext in ('.nc','.nc4'):
            #--------------------------------------------
            # Input file is netcdf format
            #--------------------------------------------
            file_unit = nc.Dataset(input_file, mode = 'r')
        else:
            #--------------------------------------------
            # Input file contains a grid or grid stack
            # as row-major, binary file with no header.
            #--------------------------------------------
            file_unit = open(input_file, 'rb')

    #print('file_unit:')
    #print(file_unit)
    return file_unit

#   open_file()
#-------------------------------------------------------------------
def read_next2(self, var_name, rti, dtype='float32', factor=1.0):

#     exec( 'file_unit = self.' + var_name + '_unit' )
#     exec( 'var_type  = self.' + var_name + '_type' )

    ####### (2019-10-03, For Python 3)
    file_unit = eval('self.' + var_name + '_unit' )
    var_type  = eval('self.' + var_name + '_type' )
    
    if (var_type.lower() == 'scalar'): 
        #-------------------------------------------
        # Scalar value was entered by user already
        #-------------------------------------------
        data = None
    elif (var_type.lower() == 'time_series'):  
        #----------------------------------------------
        # Time series: Read scalar value from file.
        # File is ASCII text with one value per line.
        #----------------------------------------------
        data = read_scalar(file_unit, dtype)
    elif (var_type.lower() in ['grid', 'grid_sequence']):   
        #--------------------------------------
        # Single grid: Read grid from file
        # read_grid() accounts for byte order
        #----------------------------------------------
        # NB!  grid_type argument allows DEM to be
        # read for GW vars, which might not be FLOAT'
        #----------------------------------------------
        data = read_grid(file_unit, rti, dtype)
    else:
        raise RuntimeError('No match found for ' + var_type + '.')
        return None

    #---------------------------------------------
    # Multiply by a conversion or scale factor ?
    #---------------------------------------------
    if (factor != 1) and (data is not None):
        data = (data * factor)

    #------------------------------------------------------
    # Perform an "in-place" update of the variable and
    # preserve reference.  This should upcast, if needed,
    # from float32 to float64.  (11/15/16)
    #------------------------------------------------------
    self.update_var( var_name, data )

    #-----------------------------------------------------
    # Values must usually be read from file as FLOAT32
    # but then need to be returned as FLOAT64. (5/17/12)
    # But np.float64( None ) = NaN. (5/18/12)
    #-----------------------------------------------------
#     if (data is None):
#         return data
#     else:
#         return np.float64( data )

#   read_next2()
#-------------------------------------------------------------------
def read_next(file_unit, var_type, rti,
              dtype='float32', units_factor=1.0,
              NGEN_CSV=False,time_index=None):

    if (var_type.lower() == 'scalar'): 
        #-------------------------------------------
        # Scalar value was entered by user already
        #-------------------------------------------
        data = None
    elif (var_type.lower() == 'time_series'):  
        #----------------------------------------------
        # Time series: Read scalar value from file.
        # File is ASCII text with one value per line.
        #----------------------------------------------
        if (NGEN_CSV):
            col=1;  delim=','
        else:
            col=0;  delim=' '      
        data = read_scalar(file_unit, dtype, col=col,
                           delim=delim )
    elif (var_type.lower() in ['grid', 'grid_sequence']):   
        #--------------------------------------
        # Single grid: Read grid from file
        # read_grid() accounts for byte order
        #----------------------------------------------
        # NB!  grid_type argument allows DEM to be
        # read for GW vars, which might not be FLOAT'
        #----------------------------------------------
        data = read_grid(file_unit, rti, dtype, time_index)

    else:
        raise RuntimeError('No match found for ' + var_type + '.')
        return None

    #---------------------------------------------
    # Multiply by a conversion or scale factor ?
    #---------------------------------------------
    if (units_factor != 1) and (data is not None):
        data *= units_factor

    #-----------------------------------------------------
    # Values must usually be read from file as FLOAT32
    # but then need to be returned as FLOAT64. (5/17/12)
    # But np.float64( None ) = NaN. (5/18/12)
    #-----------------------------------------------------
    if (data is None):
        return data
    else:
        return np.float64( data )
        ## return data

#   read_next()
#-------------------------------------------------------------------
def read_scalar(file_unit, dtype='float32',
                col=0, delim=' '):

    #-------------------------------------------------
    # Note:  Scalar values are read from text files.
    #-------------------------------------------------
    # Note:  Added col and delim args for NextGen
    #        on 2022-11-18.
    #-------------------------------------------------
        
    #-------------------------------
    # Recheck if any of these work
    #-------------------------------
    # scalar = np.fromfile(file_unit, count=1, dtype=dtype)
    # scalar = np.fromfile(file_unit, count=1, dtype=dtype, sep="\n")
    # scalar = np.fromfile(file_unit, count=1, dtype=dtype, sep=" ")    
    # scalar = np.loadtxt(file_unit, dtype=dtype)
       
    line = file_unit.readline()
    line = line.strip()
##    print 'line =', line
    if (line != ''):
        words  = line.split( delim )
##        print 'type(words)   =', type(words)
##        print 'size(words)   =', size(words)
##        print 'words[0]      =', words[0]
##        print 'type(words[0] =', type(words[0])

        #---------------------------------       
        # Why was eval() used here?
        # For scientific notation maybe?
        # Or an old numpy version?
        #---------------------------------
        scalar = np.float32(words[ col ])
        ## scalar = np.float32(eval(words[ col ]))
    else:
        scalar = None

    return scalar

#   read_scalar()  
#-------------------------------------------------------------------
def read_grid(file_unit, rti, dtype='float32', time_index=None):

    #----------------------------------------------------
    # Note:  Read 2D grid from row-major, binary file.
    #        If there is no more data left to read from
    #        file_unit (end of file), then fromfile
    #        returns an empty array (size 0) with same
    #        dtype.  In this case, return None.
    #----------------------------------------------------

    if isinstance(file_unit,nc.Dataset):

        #determine name of netcdf variable (from ncgs.get_var_names())
        names    = list( file_unit.variables.keys() )
        dims     = list( file_unit.dimensions.keys() )
        names = [s for s in names if (s not in dims) and s != 'datetime']
        if len(names) != 1:
            print('ERROR more than one non dim variables in ',file_unit.filepath())

        #get var
        grid = file_unit.variables[names[0]][time_index,:,:]

    else:

        grid = np.fromfile(file_unit, count=rti.n_pixels, dtype=dtype)
        if (grid.size == 0):
            return None

        grid = np.reshape(grid, (rti.nrows, rti.ncols))
        if (rti.SWAP_ENDIAN):
            grid.byteswap(True)

    #-------------------------        
    return grid

#   read_grid()
#-------------------------------------------------------------------
def close_file(file_unit):

    #----------------------------------------------------------
    # Notes:  Each process module closes its own input files
    #         with a method called "close_input_files()."
    #         This function is included for completeness and
    #         may be used instead.
    #----------------------------------------------------------
    file_unit.close()

#   close_file()
#-------------------------------------------------------------------
#-------------------------------------------------------------------
##def open_file(var_type, input_file):
##
##    #-----------------------------------------------------
##    # Note:  This method's name cannot be "open" because
##    #        it calls Python's built-in "open()" method.
##    #-----------------------------------------------------
##    # print 'var_type   =', var_type
##    # print 'input_file =', input_file
##    
##    #--------------------------------------------
##    # A scalar input value was provided already
##    #--------------------------------------------
##    file_unit = None
##    if (var_type == 0):  return file_unit
##    if (input_file == ''):
##        print 'ERROR in model_input.open_file():'
##        print '    Input file is null string.'
##        # print '    variable type =' + var_type
##        return file_unit
##
##    #----------------------------------
##    # Does input file exist locally ?
##    #----------------------------------
##    if not(os.path.exists(input_file)):
##        print 'ERROR in model_input.open_file():'
##        print '    Could not find input file ='
##        print '    ' + input_file
##        # print '   ' + input_file
##        return file_unit
##            
##    if (var_type == 1):
##        #-----------------------------------------
##        # Input file contains a time series and
##        # is ASCII text with one value per line.
##        #-----------------------------------------
##        file_unit = open(input_file, 'r')
##    elif (var_type > 1):
##        #--------------------------------------------
##        # Input file contains a grid or grid stack
##        # as row-major, binary file with no header.
##        #--------------------------------------------
##        file_unit = open(input_file, 'rb')
##            
##    return file_unit
##
###   open_file()
###-------------------------------------------------------------------
##def read_next(file_unit, var_type, rti, \
##              dtype='float32', factor=1.0):
##
##    #-------------------------------------------------------
##    # (5/7/09) Allow "dtype" to be given using RTI types.
##    #-------------------------------------------------------
##    rti_types = ['BYTE','INTEGER','LONG','FLOAT','DOUBLE']
##    if (dtype.upper() in rti_types):
##        dtype_map = {'BYTE':'uint8', 'INTEGER':'int16',
##                     'LONG':'int32',
##                     'FLOAT':'float32', 'DOUBLE':'float64'}
##        dtype = dtype_map[dtype]
##
##
##    if (var_type == 0): 
##        #-------------------------------------------
##        # Scalar value was entered by user already
##        #-------------------------------------------
##        data = None
##    elif (var_type == 1):  
##        #----------------------------------------------
##        # Time series: Read scalar value from file.
##        # File is ASCII text with one value per line.
##        #----------------------------------------------
##        data = read_scalar(file_unit, dtype)
##    elif (var_type in [2,3]):   
##        #--------------------------------------
##        # Single grid: Read grid from file
##        # read_grid() accounts for byte order
##        #----------------------------------------------
##        # NB!  grid_type argument allows DEM to be
##        # read for GW vars, which might not be FLOAT'
##        #----------------------------------------------
##        data = read_grid(file_unit, rti, dtype)
##    else:
##        raise RuntimeError('No match found for "var_type".')
##        return None
##    
##    if (factor != 1) and (data is not None):
##        data = (data * factor)
##    return data
##
###   read_next()
#-------------------------------------------------------------------




