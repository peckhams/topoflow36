"""
This module reads and writes a series of data cubes to a NetCDF file.


=== Set up some parameters ===

Define the shape of the data cube.
>>> shape = (4, 5, 3)
>>> res = (1., 100., 200.)

Set the values of the data cube.
>>> var = np.random.random (shape)

Write cubes at these times.
>>> times = range (5)

This is the name of the variable that we will read/write.
>>> var_name = 'grain_size'


=== Write a cube ===

>>> import nccs_files
>>> f = nccs_files.nccs_file ()

Open the file for writing.
>>> f.open_new_file ('test_file.nc', shape=shape, res=res,
...                  dtype='float64',
...                  var_name=var_name,
...                  long_name="Sediment grain size",
...                  units_name="phe",
...                  comment="Cube of sediment grain size")
True


The name of the file may not be the same as what was specified in the call
to open_new_file.  If the file already exists, a number will be appended to
the file name (before the extension).

>>> file_name = f.file_name
>>> print file_name # doctest: +ELLIPSIS
test_file....nc

Write the variable
>>> f.add_cube (var, var_name)

Close the file.
>>> f.close_file ()

=== Read a cube ===

>>> import nccs_files
>>> f = nccs_files.nccs_file ()
>>> f.open_file (file_name)
True

Read variable from file and compare it to what we wrote.
>>> var_from_file = f.get_cube (var_name, 0)
>>> (var_from_file == var).all ()
True
>>> f.close_file ()

=== Write a series of cubes ===

Open the file for writing.
>>> f.open_new_file ('test_file.nc', shape=shape, res=res,
...                  dtype='float64',
...                  var_name=var_name,
...                  long_name="Sediment grain size",
...                  units_name="phe",
...                  comment="Cube of sediment grain size")
True

>>> file_name = f.file_name
>>> print file_name # doctest: +ELLIPSIS
test_file....nc

Write the variable
>>> for time in times:
...   f.add_cube (var+time, var_name)

Close the file.
>>> f.close_file ()

=== Read a series of cubes ===

>>> import nccs_files
>>> f = nccs_files.nccs_file ()
>>> f.open_file (file_name)
True

Read variable from file and compare it to what we wrote.
>>> values_match = True
>>> for time in times:
...   var_from_file = f.get_cube (var_name, time)
...   values_match &= (var_from_file == var+time).all ()
>>> values_match
True
>>> f.close_file ()

"""

import os
import sys
import time

import numpy as np
from . import bov_files
from . import file_utils
from . import rti_files

import netCDF4 as nc

#---------------------------------------------------------------------
# This class is for I/O of time-indexed 3D arrays to netCDF files.
#---------------------------------------------------------------------
#
#   unit_test()
#   save_nccs_cube()   # (not ready yet)
#
#   class nccs_file():
#
#       import_netCDF4()
#       open_file()
#       get_dtype_map()
#       open_new_file()
#       update_time_index()
#----------------------------
#       add_cube()
#       get_cube()
#----------------------------
#       close_file()
#       close()
#
#-------------------------------------------------------------------
def save_nccs_cube(nccs_file_name=None, rtg_file_name=None):

    nccs = nccs_file ()
    OK = nccs.open_file(nccs_file_name)
    if not OK:
        return

    var_name  = 'H'
    time_index = 200
    cube = nccs.get_cube(var_name, time_index)
    nccs.close()
    
    cube = np.array(cube)
    print('min(cube), max(cube) =', cube.min(), cube.max())

    rtg_unit = open(rtg_file_name, 'wb')
    cube.tofile(unit)
    rtg_unit.cube()

#   save_nccs_cube()
#-------------------------------------------------------------------
class nccs_file():

    #----------------------------------------------------------
    # Note:  ncgs = NetCDF Grid Stack (used by CSDMS)
    #----------------------------------------------------------
    def import_netCDF4(self):

        try:
            import netCDF4
            # print 'Imported netCDF4 version: ' + netCDF4.__version__
            return netCDF4
        except:
##            print ' '
##            print 'SORRY, Cannot write netCDF files because'
##            print 'the "netCDF4" package cannot be imported.'
##            print ' '
##            python_version = sys.version[:3]
##            if (python_version != '2.6'):
##                print 'Note that "PyNIO" is only installed for'
##                print 'Python version 2.6 on "beach".'
##                print 'The current Python version is:', python_version
##                print ' '
            return False
        
    #   import_netCDF4()
    #----------------------------------------------------------
    def open_file(self, file_name):
        
        #-------------------------
        # Open file to read only
        #-------------------------
        try:
            nccs_unit = nc.Dataset(file_name, mode='r')
            self.nccs_unit = nccs_unit
            return True
        except:
            return False
    
    #   open_file()
    #----------------------------------------------------------
    def get_dtype_map(self):

        #----------------------------------------
        # Possible settings for "dtype_code"
        #----------------------------------------------------
        # These two-char codes are used for netCDF4 package
        #----------------------------------------------------
        # See:  http://unidata.github.io/netcdf4-python/
        #----------------------------------------------------
        dtype_map = {'float64':'f8', 'float32':'f4',
                     'int64':'i8', 'int32':'i4',
                     'int16':'i2', 'int8':'i1',
                     'S|100':'S1'}  # ( ????? )       
        
        #-------------------------------------------------
        # These one-char codes are used for Nio in PyNIO
        #-------------------------------------------------
        # dtype_code = "d"  # (double, Float64)
        # dtype_code = "f"  # (float,  Float32)
        # dtype_code = "l"  # (long,   Int64)
        # dtype_code = "i"  # (int,    Int32)
        # dtype_code = "h"  # (short,  Int16)
        # dtype_code = "b"  # (byte,   Int8)
        # dtype_code = "S1" # (char)
        #-------------------------------------------
#         dtype_map = {'float64':'d', 'float32':'f',
#                         'int64':'l', 'int32':'i',
#                         'int16':'s', 'int8':'b',
#                         'S|100':'S1'}  # (check last entry)                      

        return dtype_map
    
    #   get_dtype_map()
    #----------------------------------------------------------
    def open_new_file(self, file_name, info=None,
                      time_info=None,
                      var_name='X',
                      long_name=None,
                      units_name='None',
                      dtype='float64',
                      time_units='minutes',
                      comment='',
                      shape=(1,1,1),
                      res=(1.,1.,1.),
                      MAKE_RTI=True, MAKE_BOV=False):

        #----------------------------
        # Does file already exist ?
        #----------------------------
        file_name = file_utils.check_overwrite( file_name )
        self.file_name = file_name
        
        #---------------------------------------
        # Check and store the grid information
        #---------------------------------------
        self.format     = 'nccs'
        self.file_name  = file_name
        self.time_index = 0
        self.var_name   = var_name
        self.shape      = shape
        self.res        = res
        
        if (long_name is None):
            long_name = var_name
        self.long_name  = long_name
        self.units_name = units_name
        self.dtype      = dtype

        #-----------------------------------
        # Save the two-char data type code
        #-----------------------------------
        dtype_map  = self.get_dtype_map()        
        dtype_code = dtype_map[ dtype.lower() ]        
        self.dtype_code = dtype_code
        
        #-------------------------------------
        # Open a new netCDF file for writing
        #-------------------------------------   
        try:
            ## format = 'NETCDF4'
            format = 'NETCDF4_CLASSIC'
            nccs_unit = nc.Dataset(file_name, mode='w', format=format)
            OK = True
        except:
            OK = False
            return OK

        #------------------------------------------------------------
        # Option to pre-fill with fill values
        # Set fill_value for a var with "var._Fill_Value = number"
        # For Nio was:  opt.PreFill = False # (for efficiency)
        #------------------------------------------------------------
        ncgs_unit.set_fill_off()
        # ncgs_unit.set_fill_on()

        #-------------------------------------
        # Prepare and save a history string
        #-------------------------------------
        # Sample output from time.asctime():
        #     "Thu Oct  8 17:10:18 2009"
        #-------------------------------------
        history = "Created using netCDF4 " + nc.__version__ + " on "
        history = history + time.asctime() + ". " 
        history = history + comment
        nccs_unit.history = history
        # print 'MADE IT PAST history BLOCK'
        
        #----------------------------------------------
        # Create grid dimensions nx and ny, plus time
        #----------------------------------------------
        # Without using "int()" here, we get this:
        #     TypeError: size must be None or integer
        #----------------------------------------------
        nccs_unit.createDimension('nz', self.shape[0])
        nccs_unit.createDimension('ny', self.shape[1])
        nccs_unit.createDimension('nx', self.shape[2])
        nccs_unit.createDimension('time', None)   # (unlimited dimension)
        # print 'MADE IT PAST create_dimension CALLS.'
        
        #-------------------------
        # Create a time variable
        #------------------------------------------
        #('d' = float64; must match in add_cube()
        #------------------------------------------
        tvar = nccs_unit.createVariable ('time', 'f8', ('time',))
        nccs_unit.variables['time'].units = time_units
        
        #--------------------------------
        # Create a variable in the file
        #----------------------------------
        # Returns "var" as a PyNIO object
        #----------------------------------
        var = nccs_unit.createVariable (var_name, dtype_code,
                                         ('time', 'nz', 'ny', 'nx'))

        #----------------------------------
        # Specify a "nodata" fill value ?
        #----------------------------------
        # var._Fill_Value = -9999.0    ## Used for pre-fill above ?
        
        #------------------------------------
        # Create attributes of the variable
        #------------------------------------
        nccs_unit.variables[var_name].long_name = long_name
        nccs_unit.variables[var_name].units = units_name
        #----------------------------------------------------
        nccs_unit.variables[var_name].dz = self.res[0]
        nccs_unit.variables[var_name].dy = self.res[1]
        nccs_unit.variables[var_name].dx = self.res[2]
        nccs_unit.variables[var_name].y_south_edge = 0.
        nccs_unit.variables[var_name].y_north_edge = self.res[1]*self.shape[1]
        nccs_unit.variables[var_name].x_west_edge = 0.
        nccs_unit.variables[var_name].x_east_edge = self.res[2]*self.shape[2]
        nccs_unit.variables[var_name].z_bottom_edge = 0.
        nccs_unit.variables[var_name].z_top_edge = self.res[0]*self.shape[0]
        #--------------------------------------------
        # This is how it is done in ncgs_files.py.
        #-------------------------------------------
#         ncgs_unit.variables[var_name].dx           = self.info.xres
#         ncgs_unit.variables[var_name].dy           = self.info.yres  ### (12/2/09)
#         ncgs_unit.variables[var_name].y_south_edge = self.info.y_south_edge
#         ncgs_unit.variables[var_name].y_north_edge = self.info.y_north_edge
#         ncgs_unit.variables[var_name].x_west_edge  = self.info.x_west_edge
#         ncgs_unit.variables[var_name].x_east_edge  = self.info.x_east_edge  
                
        self.nccs_unit = nccs_unit
        return OK
    
    #   open_new_file()
    #----------------------------------------------------------
    def add_cube(self, grid, var_name, time=None,
                 time_index=-1):

        #---------------------------------
        # Assign a value to the variable
        #-------------------------------------------
        # This syntax works for scalars and grids
        #-------------------------------------------
        # nc_unit.variables[var_name].assign_value( grid )
 
        #-------------------------------------
        # Can use time_index to overwrite an
        # existing grid vs. simple append.
        #-------------------------------------
        if time_index == -1:
            time_index = self.time_index
        if time is None:
            time = np.float64(time_index)
            
        #---------------------------------------
        # Write a time to existing netCDF file
        #---------------------------------------
        times = self.nccs_unit.variables['time']
        times[time_index] = time
        
        #---------------------------------------
        # Write a grid to existing netCDF file
        #---------------------------------------
        var = self.nccs_unit.variables[var_name]
        if np.ndim(grid) == 0:
            #-----------------------------------------------
            # "grid" is actually a scalar (dynamic typing)
            # so convert it to a grid before saving
            #-----------------------------------------------
            grid2 = grid + np.zeros([self.nz, self.ny, self.nx],
                                        dtype=self.dtype)
            var[time_index] = grid2.astype(self.dtype)
        else:
            var[time_index] = grid.astype(self.dtype)

        #---------------------------
        # Increment the time index
        #---------------------------            
        self.time_index += 1

    #   add_cube()
    #----------------------------------------------------------
    def get_cube(self, var_name, time_index):

        var = self.nccs_unit.variables[var_name]
        return var[time_index]
        
    #   get_cube()
    #-------------------------------------------------------------------
    def close_file(self):

        # self.ncgs_unit.sync()  ## (netCDF4 has no "flush")
        self.nccs_unit.close ()

    #   close_file()
    #-------------------------------------------------------------------
    def close(self):

        # self.ncgs_unit.sync()  ## (netCDF4 has no "flush")
        self.nccs_unit.close ()

    #   close()
    #-------------------------------------------------------------------
   
if __name__ == "__main__":
    import doctest
    doctest.testmod()
 
