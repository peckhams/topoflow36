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

#-------------------------------------------------------------------

import os
import sys
import time

import numpy as np
from . import bov_files
from . import file_utils
from . import rti_files
from . import time_utils

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
    # Note:  nccs = NetCDF Cube Stack (used by CSDMS)
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
        # dtype_code = "d"  # (double, float64)
        # dtype_code = "f"  # (float,  float32)
        # dtype_code = "l"  # (long,   int64)
        # dtype_code = "i"  # (int,    int32)
        # dtype_code = "h"  # (short,  int16)
        # dtype_code = "b"  # (byte,   int8)
        # dtype_code = "S1" # (char)
        #-------------------------------------------
#         dtype_map = {'float64':'d', 'float32':'f',
#                         'int64':'l', 'int32':'i',
#                         'int16':'s', 'int8':'b',
#                         'S|100':'S1'}  # (check last entry)                      

        return dtype_map
    
    #   get_dtype_map()
    #----------------------------------------------------------
    def open_new_file(self, file_name,
                      grid_info=None,
                      time_info=None,
                      z_values=[None],  # Added 2021-07-15
                      z_units='m',      # Added 2021-07-15
                      var_name='X',
                      long_name=None,
                      units_name=None,
                      dtype='float64',
                      time_units='minutes', time_res='60.0',
                      comment='',
                      MAKE_RTI=True, MAKE_BOV=False):

        #----------------------------
        # Does file already exist ?
        #----------------------------
        file_name = file_utils.check_overwrite( file_name )
        
        #-------------------------------------------
        # Check and store the grid information
        # ncgs_files.py defines & uses this function
        #---------------------------------------------
        # self.check_and_store_info( file_name, grid_info, var_name,
        #                           dtype, MAKE_RTI, MAKE_BOV )
                
        #---------------------------------------
        # Check and store the grid information
        #---------------------------------------
        self.format     = 'nccs'
        self.file_name  = file_name
        self.time_index = 0
        self.var_name   = var_name
        #-----------------------------------------------
        if (long_name is None): long_name = var_name
        self.long_name  = long_name
        self.units_name = units_name
        self.dtype      = dtype
        #-----------------------------------------------            
        self.z_values  = z_values
        self.z_units   = z_units
        nz             = np.size(z_values)

        #---------------------------------------------------
        # Create time metadata strings (Added: 2021-07-15)
        #---------------------------------------------------
        start_date     = time_info.start_date
        start_time     = time_info.start_time
        end_date       = time_info.end_date
        end_time       = time_info.end_time
        start_datetime = time_info.start_datetime
        end_datetime   = time_info.end_datetime
        dur_units      = time_units        
        duration = time_utils.get_duration( start_date, start_time,
                                            end_date, end_time,
                                            dur_units)
        self.start_datetime = start_datetime  # for add_cube()
       
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
            format = 'NETCDF4'  # (better string support)
            ## format = 'NETCDF4_CLASSIC'
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
        nccs_unit.set_fill_off()
        # nccs_unit.set_fill_on()
  
        #----------------------------------------------
        # Create grid dimensions X, Y, Z, and time
        #----------------------------------------------
        # Without using "int()" here, we get this:
        #     TypeError: size must be None or integer
        #----------------------------------------------
#         ncols = int(self.info.ncols)   # also works
#         nrows = int(self.info.nrows)   # also works
        ncols = int(grid_info.ncols)
        nrows = int(grid_info.nrows)
        nccs_unit.createDimension('Z', int(nz))
        nccs_unit.createDimension('Y', nrows)
        nccs_unit.createDimension('X', ncols)
        nccs_unit.createDimension('time', None)   # (unlimited dimension)
        nccs_unit.createDimension('bnds', 1)

        #-------------------------------
        # Create X_vector and Y_vector
        #-------------------------------
#         minlon = self.info.x_west_edge
#         maxlon = self.info.x_east_edge
#         minlat = self.info.y_south_edge
#         maxlat = self.info.y_north_edge
        minlon = grid_info.x_west_edge
        maxlon = grid_info.x_east_edge
        minlat = grid_info.y_south_edge
        maxlat = grid_info.y_north_edge
        
        #-----------------------------------------------------
        # Method 1 (maybe safer, if any bounds are wrong)
        # xres_deg = (self.info.xres / 3600.0)  # [deg]
        # yres_deg = (self.info.yres / 3600.0)  # [deg]
        # X_vector = np.arange( minlon, maxlon, xres_deg, dtype='float64')
        # Y_vector = np.arange( minlat, maxlat, yres_deg, dtype='float64' )
        #------------------------------------------------------
        # Method 2
        X_vector = np.linspace( minlon, maxlon, ncols, dtype='float64' )
        Y_vector = np.linspace( minlat, maxlat, nrows, dtype='float64' )
        Z_vector = np.float64( z_values )  # maybe not uniformly spaced
        
        #------------------------------------------------------ 
        # According to netCDF4 docs:
        # "The dimensions themselves are usually also defined
        # as variables, called *coordinate variables*."
        #------------------------------------------------------
        
        #-----------------------------------
        # Create coordinate variable, time
        #---------------------------------------------------
        #('f8' = float32; must match in add_grid()
        #------------------------------------------------------------
        # If using the NETCDF4 format (vs. NETCDF4_CLASSIC),
        # then for a fixed-length string (e.g. 2021-07-01 00:00:00)
        # the 2nd argument here can be 'S19', and for a variable-
        # length string it can just be 'str'.
        # Otherwise, you must use: netCDF4.stringtochar() to
        # convert the string array to a character array.
        #------------------------------------------------------------               
        tvar  = nccs_unit.createVariable('time', 'f8', ('time',))   
        dtvar = nccs_unit.createVariable('datetime', 'S19', ('time',))        
        
        #-----------------------------------------------------------
        # Create coordinate variables, time, X, Y and Z
        #-----------------------------------------------------------           
        Xvar = nccs_unit.createVariable('X', 'f8', ('X',))
        Yvar = nccs_unit.createVariable('Y', 'f8', ('Y',))
        Zvar = nccs_unit.createVariable('Z', 'f8', ('Z',))

        #-----------------------------------------
        # Create a computed variable in the file
        #-----------------------------------------
        # Vars must be in the order shown
        #----------------------------------
        var = nccs_unit.createVariable(var_name, dtype_code,
                                      ('time', 'Z', 'Y', 'X'))

        #----------------------------------
        # Specify a "nodata" fill value ?
        #----------------------------------
        # var._Fill_Value = -9999.0    ## Used for pre-fill above ?
        
        #-------------------------------------
        # Prepare and save a history string
        #-------------------------------------
        # Sample output from time.asctime():
        #     "Thu Oct  8 17:10:18 2009"
        #-------------------------------------
        history = "Created using netCDF4 " + nc.__version__ + " on "
        history = history + time.asctime() + ". "
        
        #---------------------------------------------------       
        # Create title, summary and other metadata strings
        #---------------------------------------------------
        title = 'Cube stack for variable: ' + long_name
        tf_version = str(tf_utils.TF_Version_Number())
        summary  = 'This file contains a stack of 3D spatial grids, indexed '
        summary += 'by time for the single variable: ' + long_name + '.'
        email = 'Scott.Peckham@colorado.edu'
        date_created = str( datetime.date.today() )
        naming_authority = 'edu.isi.workflow'
        if (comment == ''):
            comment = 'Created by TopoFlow version ' + tf_version + '.'
        else:
            history += comment

        #-----------------------------------------
        # Save some global attributes (metadata)
        #-----------------------------------------
        nccs_unit.title             = title
        nccs_unit.summary           = summary
        nccs_unit.comment           = comment
        nccs_unit.history           = history
        nccs_unit.creator_email     = email
        nccs_unit.date_created      = date_created 
        nccs_unit.naming_authority  = naming_authority   
        nccs_unit.geospatial_bounds_crs = '+init=epsg:4979'
        bounds = [minlon, minlat, maxlon, maxlat]   #### MINT order
        nccs_unit.geospatial_bounds = bounds 

        #---------------------------------------
        # Save attributes of geospatial var X
        #---------------------------------------
        nccs_unit.variables['X'][:] = X_vector
        nccs_unit.variables['X'].long_name = 'longitude'
        nccs_unit.variables['X'].units = 'degrees_east'       
        nccs_unit.variables['X'].geospatial_lon_min = minlon
        nccs_unit.variables['X'].geospatial_lon_max = maxlon
        nccs_unit.variables['X'].xres_arcsec = grid_info.xres
        ## nccs_unit.variables['X'].xres_arcsec = self.info.xres
        # Will need something like this for UTM coords later
        # nccs_unit.variables['X'].x_west_edge = minlon
        # nccs_unit.variables['X'].x_east_edge = maxlon  
        
        #---------------------------------------
        # Save attributes of geospatial var Y
        #---------------------------------------
        nccs_unit.variables['Y'][:] = Y_vector
        nccs_unit.variables['Y'].long_name = 'latitude'
        nccs_unit.variables['Y'].units = 'degrees_north'     
        nccs_unit.variables['Y'].geospatial_lat_min = minlat
        nccs_unit.variables['Y'].geospatial_lat_max = maxlat
        nccs_unit.variables['Y'].yres_arcsec = grid_info.yres
        ## nccs_unit.variables['Y'].yres_arcsec = self.info.yres
        # Will need something like this for UTM coords later
        # nccs_unit.variables['Y'].y_south_edge  = minlat
        # nccs_unit.variables['Y'].y_north_edge  = maxlat

        #---------------------------------------
        # Save attributes of geospatial var Z
        # Assume z-axis points downward from land surface
        #--------------------------------------------------
        nccs_unit.variables['Z'][:] = Z_vector
        nccs_unit.variables['Z'].long_name = 'depth'
        nccs_unit.variables['Z'].units = z_units     
        ## nccs_unit.variables['Z'].zres = zres  # can be variable
        nccs_unit.variables[var_name].z_top    = 0.0
        nccs_unit.variables[var_name].z_bottom = z_bottom
                
        #------------------------------------------
        # Save attributes of coordinate var, time
        #----------------------------------------------------
        # Note!  add_values_at_IDs() builds a time vector
        #        but starts at 0 and doesn't include dates.
        #        Recall time is an unlimited dimension.
        #----------------------------------------------------
        nccs_unit.variables['time'].long_name = 'time'
        nccs_unit.variables['time'].units = time_units
        nccs_unit.variables['time'].time_units = time_units
        nccs_unit.variables['time'].time_coverage_resolution = time_res  
        nccs_unit.variables['time'].time_coverage_start = start_datetime 
        nccs_unit.variables['time'].time_coverage_end = end_datetime 
        nccs_unit.variables['time'].time_coverage_duration = duration

        #-----------------------------------------
        # Save attributes of extra var, datetime
        #-------------------------------------------------------------
        # Note: The duration is determined from start_datetime,
        #       end_datetime and dur_units.  However, if the model
        #       is run for a shorter time period, then the T_vector
        #       computed here may be longer than the time var vector.
        #       This results in nodata values "--" in time vector.
        #-------------------------------------------------------------
        # Note!  add_profile() now builds a datetime vector.
        #        Recall time is an unlimited dimension.
        #-------------------------------------------------------------
        time_dtype = time_utils.get_time_dtype( time_units )
        nccs_unit.variables['datetime'].long_name = 'datetime' 
        nccs_unit.variables['datetime'].units = time_dtype
                            
        #---------------------------------------
        # Save attributes of the main variable
        #---------------------------------------
        svo_name = svo_names.get_svo_name( var_name )
        nccs_unit.variables[var_name].svo_name  = svo_name 
        nccs_unit.variables[var_name].long_name = long_name
        nccs_unit.variables[var_name].units     = units_name
        nccs_unit.variables[var_name].n_cubes   = 0 
                                       
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
        # if time is None:
        #    time = np.float64(time_index)
            
        #----------------------------------------------
        # Write current time to existing netCDF file
        # Recall that time has an unlimited dimension
        #----------------------------------------------
        times = self.nccs_unit.variables[ 'time' ]
        times[ time_index ] = time

        #-------------------------------------------------
        # Write current datetime to existing netCDF file
        # Recall that time has an unlimited dimension
        # Datetime strings have netCDF4 type 'S19'
        #-------------------------------------------------
        datetime = time_utils.get_current_datetime(
                              self.start_datetime,
                              time, time_units='minutes')
        datetimes = self.nccs_unit.variables[ 'datetime' ]
        datetimes[ time_index ] = str(datetime)
        
        #---------------------------------------
        # Write a grid to existing netCDF file
        #---------------------------------------
        var = self.nccs_unit.variables[var_name]
        var.n_cubes += 1
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

        # self.nccs_unit.sync()  ## (netCDF4 has no "flush")
        self.nccs_unit.close ()

    #   close_file()
    #-------------------------------------------------------------------
    def close(self):

        # self.nccs_unit.sync()  ## (netCDF4 has no "flush")
        self.nccs_unit.close ()

    #   close()
    #-------------------------------------------------------------------
   
if __name__ == "__main__":
    import doctest
    doctest.testmod()
 
