
# S.D. Peckham
# Sept 2014 (new version to use netCDF4)
# June 2010 (streamlined a bit more)
# December 2, 2009 (updated open_new_file to use "info")
# October 13, 2009
# Nov. 2019 (changed how netCDF is written in open_new_file())

import os
import sys
import time
import datetime

import numpy as np
from . import bov_files
from . import file_utils
from . import rti_files
from . import svo_names
from . import tf_utils

import netCDF4 as nc

#-------------------------------------------------------------------

#   unit_test()
#   save_ncgs_frame()    ## (12/7/09)
#
#   class ncgs_file():
#
#       import_netCDF4()
#       open_file()
#       check_and_store_info()   # (12/2/09)
#       get_dtype_map()
#       open_new_file()
#       add_grid()
#       get_var_names()       # 2019-11-21
#       get_var_long_name()   # 2019-11-21
#       get_var_units()       # 2019-11-21
#       get_grid()
#       close_file()
#       close()
#
#-------------------------------------------------------------------
def unit_test(nx=4, ny=6, n_grids=6, VERBOSE=False,
              file_name="NCGS_Grid_Test.nc"):

    print('Running unit_test()...')

    #-------------------------------------
    # Make instance of ncgs_file() class
    #-------------------------------------
    ncgs = ncgs_file()
    xres_sec = 900
    yres_sec = 900
    var_name = "depth"

    info = rti_files.make_info( file_name, ncols=nx, nrows=ny,
                                xres=xres_sec, yres=yres_sec,
                                x_west_edge=0.0, x_east_edge=1.0,
                                y_south_edge=0.0, y_north_edge=1.5 )
                                
    OK = ncgs.open_new_file( file_name, info,
                             dtype='float32',
                             var_name=var_name,
                             long_name="depth of water",
                             units_name="meters",
                             comment="Created by TopoFlow 3.6.")
    if not(OK):
        print('ERROR during open_new_file().')
        return
    
    grid = np.arange(nx * ny, dtype='Float32')
    grid = grid.reshape( (ny, nx) )

    #-----------------------------------------------
    # (6/10/10) Can use this to test new ability
    # of add_grid() to convert from scalar to grid
    #-----------------------------------------------
    # grid = np.float32(0)
    
    #----------------------------------
    # Add some test grids to the file
    #----------------------------------
    print('Writing grids to NCGS file...')
    for time_index in range(n_grids):
        ncgs.add_grid( grid, var_name )  
        ## ncgs.add_grid( grid, var_name, time_index )                         
        grid = (grid + 1)
        
    if (VERBOSE):
        print(self.ncgs_unit)  # (print a summary)

    ncgs.close_file()
    print('Finished writing NCGS file: ' + file_name)
    print(' ')

    #---------------------------------------------
    # Re-open the file and read grids one-by-one 
    #------------------------------------------------
    # NB!  If errors, check if there are multiple
    # copies of the test file with a number suffix.
    #------------------------------------------------
    # ncgs = ncgs_file()  # not necessary
    OK = ncgs.open_file( file_name )
    if not(OK): return
    print('Reading grids from NCGS file: ')
    
    for time_index in range(n_grids):
        grid = ncgs.get_grid(var_name, time_index)
        print('grid[' + str(time_index) + '] = ')
        print(grid)
        print('-----------------------------------------------')
    ncgs.close_file()    
    print('Finished reading NCGS file: ' + file_name)
    print(' ')
    
#   unit_test()
#-------------------------------------------------------------------
def save_ncgs_frame(ncgs_file_name=None, rtg_file_name=None):

    ncgs = ncgs_file()
    OK = ncgs.open_file( ncgs_file_name )
    if not(OK): return

    grid_name  = 'H'
    time_index = 200
    grid = ncgs.get_grid( grid_name, time_index )
    ncgs.close()
    
    grid = np.array( grid )
    print('min(grid), max(grid) =', grid.min(), grid.max())

    rtg_unit = open( rtg_file_name, 'wb' )
    grid.tofile( unit )
    rtg_unit.close()

#   save_ncgs_frame()
#-------------------------------------------------------------------
class ncgs_file():

    #------------------------------------------------------
    # Note:  ncgs = NetCDF Grid Stack (used by CSDMS)
    #
    #        (10/9/10) Added check_netcdf() function in
    #        model_output.py that each component can call
    #        in its open_output_files() method.
    #------------------------------------------------------
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
            ncgs_unit = nc.Dataset(file_name, mode='r')
            self.ncgs_unit = ncgs_unit
            ### return ncgs_unit
            return True
        except:
            print('ERROR: Could not open file:')
            print( '   ' + file_name )
            print( 'Current working directory =')
            print( '   ' + os.getcwd() )
            return False
    
    #   open_file()
    #----------------------------------------------------------
    def check_and_store_info(self, file_name, info=None,
                             grid_name='UNKNOWN',
                             dtype='float32',
                             MAKE_RTI=True, MAKE_BOV=False):

        #-----------------------------------------------------
        # Note: This object (self) may be new or it may have
        #       been used previously.  In the latter case,
        #       "info" should still be available in "self".
        #       We only need info if MAKE_RTI or MAKE_BOV.
        #-----------------------------------------------------
        self.format     = 'NCGS'
        self.file_name  = file_name
        self.time_index = 0
        self.grid_name  = grid_name
        
        #-----------------------------------------------------
        # This was used by rts_files.check_and_store_info()
        # but is not appropriate here because we need to
        # know nx, ny, dx and dy for the netCDF file.
        #-----------------------------------------------------        
        ### if not(MAKE_RTI or MAKE_BOV): return
        
        #---------------------------------
        # Was "info" argument provided ?
        #---------------------------------
        NEW_INFO = True
        if (info is None):
            try:
                info    = self.info
                self.nx = info.ncols  ###
                self.ny = info.nrows
                NEW_INFO = False
                ## print 'Found info in state.'
            except:
                #------------------------------------------
                # Try to find RTI file to copy info from.
                # Don't create a new RTI file.
                #------------------------------------------
                RTI_file = rti_files.try_to_find_rti_file( file_name )
                if (RTI_file != 'none'):
                    print('Reading info from: ' + RTI_file)
                    info = rti_files.read_info( RTI_file )
                else:
                    print('ERROR during open_new_file():')
                    print('   Could not find RTI file and "info"')
                    print('   argument was not provided.')
                    print(' ')
                    return

        #-----------------------------
        # Update "info" as necessary
        #-----------------------------
        info.grid_file   = file_name
        info.data_type   = rti_files.get_rti_data_type( dtype )
        info.data_source = 'TopoFlow 3.0'
        info.gmin        = -9999.0
        info.gmax        = -9999.0
        
        #---------------------------------------
        # If new "info" was provided, store it
        #---------------------------------------
        if (NEW_INFO):
            self.info = info
            self.nx   = info.ncols
            self.ny   = info.nrows           
            ## print 'Stored new info in state.'

        #-------------------
        # Write RTI file ?
        #-------------------
        if (MAKE_RTI):
            prefix   = rti_files.get_file_prefix( file_name )
            RTI_file = (prefix + '.rti')
            rti_files.write_info( RTI_file, info )
            # print 'Wrote grid info to: ' + RTI_file   ######
            
        #-------------------
        # Write BOV file ?
        #-------------------
        if (MAKE_BOV):
            bov_files.write_info_as_bov( file_name, info, grid_name)
                                         ###  time )
        
    #   check_and_store_info()
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
#     def open_new_file0(self, file_name, info=None,
#                       var_name='X',
#                       long_name=None,
#                       units_name='None',
#                       dtype='float32',
#                       ### dtype='float64'
#                       time_units='minutes',
#                       comment='',
#                       MAKE_RTI=True, MAKE_BOV=False):
# 
#         #----------------------------
#         # Does file already exist ?
#         #----------------------------
#         file_name = file_utils.check_overwrite( file_name )
#         self.file_name = file_name
#         
#         #---------------------------------------
#         # Check and store the grid information
#         #---------------------------------------
#         self.check_and_store_info( file_name, info, var_name,
#                                    dtype, MAKE_RTI, MAKE_BOV )
#         if (long_name is None): long_name = var_name
#         self.long_name  = long_name
#         self.units_name = units_name
#         self.dtype      = dtype
# 
#         #-----------------------------------
#         # Save the two-char data type code
#         #-----------------------------------
#         dtype_map  = self.get_dtype_map()
#         dtype_code = dtype_map[ dtype.lower() ]
#         self.dtype_code = dtype_code
# 
#         #-------------------------------------
#         # Open a new netCDF file for writing
#         #-------------------------------------        
#         try:
#             ## format = 'NETCDF4'
#             format = 'NETCDF4_CLASSIC'
#             ncgs_unit = nc.Dataset(file_name, mode='w', format=format)
#             OK = True
#         except:
#             OK = False
#             return OK
# 
#         #------------------------------------------------------------
#         # Option to pre-fill with fill values
#         # Set fill_value for a var with "var._Fill_Value = number"
#         # For Nio was:  opt.PreFill = False # (for efficiency)
#         #------------------------------------------------------------
#         ncgs_unit.set_fill_off()
#         # ncgs_unit.set_fill_on()
#         
#         #-------------------------------------
#         # Prepare and save a history string
#         #-------------------------------------
#         # Sample output from time.asctime():
#         #     "Thu Oct  8 17:10:18 2009"
#         #-------------------------------------
#         history = "Created using netCDF4 " + nc.__version__ + " on "
#         history = history + time.asctime() + ". " 
#         history = history + comment
#         ncgs_unit.history = history
#         # print 'MADE IT PAST history BLOCK'
#         
# ##        print 'nx =', self.info.ncols
# ##        print 'ny =', self.info.nrows
# ##        print 'dx =', self.info.xres
# ##        print 'dy =', self.info.yres
# ##        print ' '
#         
#         #----------------------------------------------
#         # Create grid dimensions nx and ny, plus time
#         #----------------------------------------------
#         # Without using "int()" here, we get this:
#         #     TypeError: size must be None or integer
#         #----------------------------------------------
#         ncgs_unit.createDimension('nx', int(self.info.ncols))
#         ncgs_unit.createDimension('ny', int(self.info.nrows))
#         ncgs_unit.createDimension('time', None)   # (unlimited dimension)
#         # print 'MADE IT PAST create_dimension CALLS.'
#         
#         #-------------------------
#         # Create a time variable
#         #------------------------------------------
#         #('d' = float64; must match in add_grid()
#         # In netCDF4 package, use 'f8' vs. 'd'.
#         #------------------------------------------
#         tvar = ncgs_unit.createVariable('time', 'f8', ('time',))
#         tvar.units = time_units
#         ### ncgs_unit.variables['time'].units = time_units
#         
#         #--------------------------------
#         # Create a variable in the file
#         #--------------------------------
#         var = ncgs_unit.createVariable(var_name, dtype_code,
#                                         ('time', 'ny', 'nx'))
# 
#         #----------------------------------
#         # Specify a "nodata" fill value ?
#         #----------------------------------
#         # var._Fill_Value = -9999.0    ## Used for pre-fill above ?
#         
#         #-------------------------------------------
#         # Create a separate, scalar "time stamp" ?
#         #-------------------------------------------
#         # t = nc_unit.createVariable('time', dtype_code, ('time'))
#         
#         #----------------------------------
#         # Specify a "nodata" fill value ?
#         #----------------------------------
# #         var._FillValue = -9999.0    ## Was used for Nio.
#         
#         #------------------------------------
#         # Create attributes of the variable
#         #------------------------------------
#         ncgs_unit.variables[var_name].long_name    = long_name
#         ncgs_unit.variables[var_name].units        = units_name
#         ncgs_unit.variables[var_name].dx           = self.info.xres
#         ncgs_unit.variables[var_name].dy           = self.info.yres  ### (12/2/09)
#         ncgs_unit.variables[var_name].y_south_edge = self.info.y_south_edge
#         ncgs_unit.variables[var_name].y_north_edge = self.info.y_north_edge
#         ncgs_unit.variables[var_name].x_west_edge  = self.info.x_west_edge
#         ncgs_unit.variables[var_name].x_east_edge  = self.info.x_east_edge        
#         
#         self.ncgs_unit = ncgs_unit
#         return OK
#     
#     #   open_new_file0()
    #----------------------------------------------------------
    def open_new_file(self, file_name,
                      grid_info=None,
                      time_info=None,
                      var_name='Z',
                      long_name=None,
                      units_name=None,
                      dtype='float32',
                      ### dtype='float64'
                      time_units='minutes', time_res='60.0',
                      comment='',
                      MAKE_RTI=True, MAKE_BOV=False):

        #----------------------------
        # Does file already exist ?
        #----------------------------
        file_name = file_utils.check_overwrite( file_name )
        self.file_name = file_name
        
        #---------------------------------------
        # Check and store the grid information
        #---------------------------------------
        self.check_and_store_info( file_name, grid_info, var_name,
                                   dtype, MAKE_RTI, MAKE_BOV )
        if (long_name is None): long_name = var_name
        self.long_name  = long_name
        self.units_name = units_name
        self.dtype      = dtype

        #---------------------------------------------
        # Create time metadata strings  (2020-01-14)
        #---------------------------------------------
        # time.asctime() = 'Fri Jan 10 10:48:48 2020'
        # str(datetime.date.today().ctime()) =
        #      'Fri Jan 10 00:00:00 2020' 
        # str(datetime.date.today()) = '2010-01-10'
        # str(datetime.datetime.now()) =
        #   '2020-01-14 12:35:32.087911'
        #---------------------------------------------------
        # str(datetime.datetime.now().ctime()) =
        #   'Tue Jan 14 12:36:48 2020'
        #---------------------------------------------------
        # x = datetime.datetime(2018, 9, 15, 12, 45, 35)
        # str(x) = '2018-09-15 12:45:35'
        # x.strftime("%b %d %Y %H:%M:%S") =
        #    Sep 15 2018 12:45:35
        #---------------------------------------------------
#         conversion_factor_map = { 'years': 31536000, 
#         'days': 86400, 'hours': 3600, 'minutes': 60, 'seconds': 1 }
#         factor = conversion_factor_map[ time_units ]
#         time_res_sec = factor * int(time_res)
#         time_res_sec_str = str(time_res_sec)
        #----------------------------------------------------------
        start_date = time_info.start_date
        start_time = time_info.start_time
        end_date   = time_info.end_date
        end_time   = time_info.end_time
        #-----------------------------------------------
        start_datetime = start_date + ' ' + start_time
        end_datetime   = end_date   + ' ' + end_time
        dur_units      = time_units        
        duration = tf_utils.get_duration( start_date, start_time,
                                          end_date, end_time,
                                          dur_units)
                                          ## time_res_sec_str) 
        
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
            ncgs_unit = nc.Dataset(file_name, mode='w', format=format)
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

        #----------------------------------------------
        # Create grid dimensions X, Y, and time
        #----------------------------------------------
        # Without using "int()" here, we get this:
        #     TypeError: size must be None or integer
        #----------------------------------------------
#         ncols = int(self.info.ncols)   # also works
#         nrows = int(self.info.nrows)   # also works
        ncols = int(grid_info.ncols)
        nrows = int(grid_info.nrows)
        ncgs_unit.createDimension('X', ncols)
        ncgs_unit.createDimension('Y', nrows)
        ncgs_unit.createDimension('time', None)   # (unlimited dimension)
        ncgs_unit.createDimension('bnds', 1)

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

        #------------------------------------------------------ 
        # According to netCDF4 docs:
        # "The dimensions themselves are usually also defined
        # as variables, called *coordinate variables*."
        #------------------------------------------------------
        #('d' = float64; must match in add_grid()
        # In netCDF4 package, use 'f8' vs. 'd'.
        #-----------------------------------------------------------
        # Create coordinate variables, time, X, and Y
        #-----------------------------------------------------------        
        tvar = ncgs_unit.createVariable('time', 'f8', ('time',))     
        Xvar = ncgs_unit.createVariable('X', 'f8', ('X',))
        Yvar = ncgs_unit.createVariable('Y', 'f8', ('Y',))

        #-----------------------------------------
        # Create a computed variable in the file
        #-----------------------------------------
        # Note:  Y must come before X here !
        #------------------------------------------
        var = ncgs_unit.createVariable(var_name, dtype_code,
                                        ('time', 'Y', 'X'))

        #----------------------------------
        # Specify a "nodata" fill value ?
        #----------------------------------
        # var._Fill_Value = -9999.0   ## Used for pre-fill above ?
        # var._FillValue = -9999.0    ## Was used for Nio.

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
        title = 'Grid stack for variable: ' + long_name
        tf_version = str(tf_utils.TF_Version_Number())
        summary  = 'This file contains a stack of spatial grids, indexed '
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
        ncgs_unit.title             = title
        ncgs_unit.summary           = summary
        ncgs_unit.comment           = comment
        ncgs_unit.history           = history
        ncgs_unit.creator_email     = email
        ncgs_unit.date_created      = date_created 
        ncgs_unit.naming_authority  = naming_authority   
        ncgs_unit.geospatial_bounds_crs = '+init=epsg:4979'
        bounds = [minlon, minlat, maxlon, maxlat]   #### MINT order
        ncgs_unit.geospatial_bounds = bounds 
        
        #---------------------------------------
        # Save attributes of geospatial var X
        #---------------------------------------
        ncgs_unit.variables['X'][:] = X_vector
        ncgs_unit.variables['X'].long_name = 'longitude'
        ncgs_unit.variables['X'].units = 'degrees_east'       
        ncgs_unit.variables['X'].geospatial_lon_min = minlon
        ncgs_unit.variables['X'].geospatial_lon_max = maxlon
        ncgs_unit.variables['X'].xres_sec = grid_info.xres
        ## ncgs_unit.variables['X'].xres_sec = self.info.xres
        # Will need something like this for UTM coords later
        # ncgs_unit.variables['X'].x_west_edge = minlon
        # ncgs_unit.variables['X'].x_east_edge = maxlon  
        
        #---------------------------------------
        # Save attributes of geospatial var Y
        #---------------------------------------
        ncgs_unit.variables['Y'][:] = Y_vector
        ncgs_unit.variables['Y'].long_name = 'latitude'
        ncgs_unit.variables['Y'].units = 'degrees_north'     
        ncgs_unit.variables['Y'].geospatial_lat_min = minlat
        ncgs_unit.variables['Y'].geospatial_lat_max = maxlat
        ncgs_unit.variables['Y'].yres_sec = grid_info.yres
        ## ncgs_unit.variables['Y'].yres_sec = self.info.yres
        # Will need something like this for UTM coords later
        # ncgs_unit.variables['Y'].y_south_edge  = minlat
        # ncgs_unit.variables['Y'].y_north_edge  = maxlat

        #------------------------------------------
        # Save attributes of coordinate var, time
        #------------------------------------------
        ncgs_unit.variables['time'].long_name = 'time'
        ncgs_unit.variables['time'].units = time_units
        ncgs_unit.variables['time'].time_coverage_resolution = time_res    
        ncgs_unit.variables['time'].time_coverage_start = start_datetime 
        ncgs_unit.variables['time'].time_coverage_end = end_datetime 
        ncgs_unit.variables['time'].time_coverage_duration = duration
                    
        #---------------------------------------
        # Save attributes of the main variable
        #---------------------------------------
        svo_name = svo_names.get_svo_name( var_name )
        ncgs_unit.variables[var_name].svo_name= svo_name 
        ncgs_unit.variables[var_name].long_name = long_name
        ncgs_unit.variables[var_name].units     = units_name
        ncgs_unit.variables[var_name].n_grids   = 0   ##########
        #----------------------------------------------------------------
#         ncgs_unit.variables[var_name].valid_min     = valid_min
#         ncgs_unit.variables[var_name].valid_max     = valid_max
#         ncgs_unit.variables[var_name].valid_range   = valid_range
#         ncgs_unit.variables[var_name].missing_value = missing_value
#         ncgs_unit.variables[var_name].fill_value    = fill_value                         

        self.ncgs_unit = ncgs_unit
        return OK
    
    #   open_new_file()
    #----------------------------------------------------------
    def add_grid(self, grid, grid_name, time=None,
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
        if (time_index == -1):
            time_index = self.time_index
        if (time is None):
            time = np.float64( time_index )
            
        #---------------------------------------
        # Write a time to existing netCDF file
        #---------------------------------------
        # times = self.ncgs_unit.variables[ 'time' ]
        # times[ time_index ] = time  ############################ CHECK
        
        #---------------------------------------
        # Write a grid to existing netCDF file
        #---------------------------------------
        var = self.ncgs_unit.variables[ grid_name ]
        var.n_grids += 1  ##########
        if (np.ndim(grid) == 0):
            #-----------------------------------------------
            # "grid" is actually a scalar (dynamic typing)
            # so convert it to a grid before saving
            #-----------------------------------------------
            grid2 = grid + np.zeros([self.ny, self.nx],
                                       dtype=self.dtype)
            var[ time_index ] = grid2.astype(self.dtype)
        else:
            var[ time_index ] = grid.astype(self.dtype)

        #---------------------------
        # Increment the time index
        #---------------------------            
        self.time_index += 1

        #-------------------------------------------------
        # 12/2/09:  netCDF is supposed to take care of
        # byteorder transparently.  However, we need to
        # make sure we don't byteswap in the function
        # "model_output.save_as_grid_to_file()" when the
        # output format is netCDF.
        #-------------------------------------------------        
##        if (sys.byteorder == 'big'):
##            var[time_index] = grid
##        else:
##            grid2 = grid.copy()
##            var[time_index] = grid2.byteswap() 
##        self.time_index += 1
        
    #   add_grid()
    #----------------------------------------------------------
    def get_var_names(self):
    
        var_dict = self.ncgs_unit.variables
        return list( var_dict.keys() )

    #   get_var_names()
    #----------------------------------------------------------
    def get_var_long_name(self, var_name ):

        var = self.ncgs_unit.variables[ var_name ]
        return var.long_name 
            
    #   get_var_long_name()
    #----------------------------------------------------------
    def get_var_units(self, var_name ):

        var = self.ncgs_unit.variables[ var_name ]
        return var.units
            
    #   get_var_units()
    #----------------------------------------------------------
    def get_grid(self, var_name, time_index):

        var = self.ncgs_unit.variables[ var_name ]
        return var[ time_index ]
        
    #   get_grid()
    #-------------------------------------------------------------------
    def close_file(self):

        # self.ncgs_unit.sync()  ## (netCDF4 has no "flush")
        self.ncgs_unit.close()

    #   close_file()
    #-------------------------------------------------------------------
    def close(self):

        # self.ncgs_unit.sync()  ## (netCDF4 has no "flush")
        self.ncgs_unit.close()

    #   close()
    #-------------------------------------------------------------------
    
