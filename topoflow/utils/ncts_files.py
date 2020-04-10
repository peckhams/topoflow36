
# This new verion (6/11/10) hasn't been tested yet.
# Can't run unit tests on my MacPro w/o Nio.
#---------------------------------------------------

# S.D. Peckham
# Sept 2014 (new version to use netCDF4)
# May, June 2010
# Nov 2019  (MINT netCDF compliance)

import os
import sys
import time
import datetime

import numpy as np
from . import file_utils
from . import rti_files
from . import svo_names
from . import tf_utils

import netCDF4 as nc

#-------------------------------------------------------------------
# This class is for I/O of time series data to netCDF files.
#-------------------------------------------------------------------
#
#   unit_test1()
#   unit_test2()
#
#   save_as_text()   # (not ready yet)
#   get_dtype_map()
#   get_dtype_codes()      # 2020-01-26 (separate function)
#
#   class ncts_file():
#
#       import_netCDF4()
#       open_file()
#       open_new_file()
#       update_time_index()
#-----------------------------
#       add_value()
#       get_value()
#-----------------------------
#       values_at_IDs()
#       add_values_at_IDs()
#-----------------------------
#       add_series()
#       get_var_names()        # 2019-11-21
#       get_var_long_name()    # 2019-11-24
#       get_var_units()        # 2019-11-24
#       get_var_lons()         # 2019-11-24
#       get_var_lats()         # 2019-11-24
#       get_series()
#-----------------------------
#       close_file()
#       close()
#
#-------------------------------------------------------------------
def unit_test1(n_values=10, VERBOSE=False,
               file_name="NCTS_Series_Test.nc"):

    #--------------------------------------------------------
    # Notes: This test uses add_value() and get_value() to
    #        add and retrieve a time series to/from a file,
    #        one value at a time.
    #--------------------------------------------------------
    print(' ')
    print('Running unit_test1()...')

    #-------------------------------------
    # Make instance of ncts_file() class
    #-------------------------------------
    ncts = ncts_file()
    var_names = ['Q_3_6', 'Q_4_5']

    info = rti_files.make_info( file_name, ncols=10, nrows=12,
                                xres=900.0, yres=900.0,
                                x_west_edge=0.0, x_east_edge=1.0,
                                y_south_edge=0.0, y_north_edge=1.5 )
    ncts.info = info
                                    
    OK = ncts.open_new_file( file_name,
                             var_names=var_names,
                             long_names=['volumetric discharge'],
                             units_names=['m3 s-1'],
                             dtypes=['float32'],
                             comment="Created by TopoFlow 3.6.")
                             ## time_long_name='time',
                             ## time_units_name="minutes")

    ###########################################################
    # WHAT ABOUT UNITS AND LONG_NAME for the TIME VALUES ??
    ###########################################################
    
    if not(OK):
        print('ERROR during open_new_file().')
        return

    series = np.sqrt(np.arange( n_values, dtype='Float32'))
    times  = np.arange( n_values, dtype='Float32') * 60.0
    
    #--------------------------
    # Add time series to file
    #--------------------------
    print('Writing values to NCTS file...')
    for time_index in range(n_values):
        time  = times[ time_index ]
        value = series[ time_index ]
        ncts.add_value( value, var_names[0], time )
        #----------------------------------------
        ncts.update_time_index()
    if (VERBOSE):
        print(self.ncts_unit)  # (print a summary)

    ncts.close_file()
    print('Finished writing ncts file: ' + file_name)
    print(' ')

    #--------------------------------------------
    # Re-open the file and read the time series 
    #--------------------------------------------
    OK = ncts.open_file( file_name )
    if not(OK): return
    print('Reading values from ncts file: ')
    
    for time_index in range(n_values):
        value, time = ncts.get_value(var_names[0], time_index)
        ti_str = str(time_index)
        t_str  = 'time[' + ti_str + '], '
        v_str  = 'value[' + ti_str + '] = '
        print((t_str + v_str), time, value)
        ## print '-----------------------------------------------'

    #-----------------
    # Close the file
    #-----------------
    ncts.close_file()    
    print('Finished reading ncts file: ' + file_name)
    print(' ')
    
#   unit_test1()
#-------------------------------------------------------------------
def unit_test2(n_values=10, VERBOSE=False,
               file_name="NCTS_Series_Test.nc"):

    #--------------------------------------------------------
    # Notes: This test uses add_series() and get_series() to
    #        add and retrieve a time series to/from a file,
    #        all values at once.
    #--------------------------------------------------------
    print(' ')
    print('Running unit_test2()...')

    #-------------------------------------
    # Make instance of ncts_file() class
    #-------------------------------------
    ncts = ncts_file()
    var_name = "depth"

    OK = ncts.open_new_file( file_name,
                             var_names=[var_name],
                             long_names=["depth of water"],
                             units_names=["meters"],
                             dtypes=['float32'],
                             time_units='minutes',
                             comment="Created by TopoFlow 3.0.")

    ###############################################
    # WHAT ABOUT LONG_NAME for the TIME VALUES ??
    ###############################################
    
    if not(OK):
        print('ERROR during open_new_file().')
        return

    series = np.sqrt(np.arange( n_values, dtype='Float32'))
    times  = np.arange( n_values, dtype='Float32') * 60.0
    
    #--------------------------
    # Add time series to file
    #--------------------------
    print('Writing values to NCTS file...')
    ncts.add_series( series, var_names[0], times )
    #--------------------------------------------
    ncts.update_time_index( step=n_values )
        
    if (VERBOSE):
        print(self.ncts_unit)  # (print a summary)

    ncts.close_file()
    print('Finished writing ncts file: ' + file_name)
    print(' ')

    #--------------------------------------------
    # Re-open the file and read the time series 
    #--------------------------------------------
    OK = ncts.open_file( file_name )
    if not(OK): return
    print('Reading values from ncts file: ')
    series, times = ncts.get_series( var_names[0] )
    
    for n in range(n_values):
        time   = times[n]
        value  = series[n]
        ti_str = str(n)
        t_str  = 'time[' + ti_str + '], '
        v_str  = 'value[' + ti_str + '] = '
        print((t_str + v_str), time, value)
        ## print '-----------------------------------------------'

    #-----------------
    # Close the file
    #-----------------
    ncts.close_file()    
    print('Finished reading ncts file: ' + file_name)
    print(' ')
    
#   unit_test2()
#-------------------------------------------------------------------
def save_as_text(ncts_file_name=None, text_file_name=None):

    ncts = ncts_file()
    OK = ncts.open_file( ncts_file_name )
    if not(OK): return

    var_name  = 'H'
    data = ncts.get_series( var_name )
    ncts.close()
    
    data = np.array( data )
    print('min(data), max(data) =', data.min(), data.max())

    text_unit = open( text_file_name, 'w' )
    data.tofile( unit )  ###### CHECK THIS #######
    text_unit.close()

#   save_as_text()
#----------------------------------------------------------
def get_dtype_map():

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
def get_dtype_codes( dtypes, var_names ):

    #---------------------------------------------
    # Create array of dtype codes from dtypes
    # for multiple time series (i.e. columns).
    #---------------------------------------------
    n_vars = len(var_names)
    dtype_map   = get_dtype_map()
    dtype_codes = []
    if (len(dtypes) == n_vars):
        for dtype in dtypes:
           dtype_code = dtype_map[ dtype.lower() ]
           dtype_codes.append( dtype_code )
    else:
        dtype = dtypes[0]
        dtype_code = dtype_map[ dtype.lower() ]
        for k in range(n_vars):
            dtype_codes.append( dtype_code ) 
                    
    return dtype_codes               

#   get_dtype_codes()
#-------------------------------------------------------------------
class ncts_file():

    #----------------------------------------------------------
    # Note:  ncts = NetCDF Time Series (used by CSDMS)
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
            ncts_unit = nc.Dataset(file_name, mode='r')
            self.ncts_unit = ncts_unit
            ### return ncts_unit
            return True
        except:
            print('ERROR: Could not open file:')
            print( '   ' + file_name )
            print( 'Current working directory =')
            print( '   ' + os.getcwd() )
            return False
    
    #   open_file()
    #----------------------------------------------------------
    def open_new_file(self, file_name,
                      grid_info=None,
                      time_info=None,
                      var_names=['Z_2_3'],
                      long_names=['None'],
                      units_names=['None'],
                      dtypes=['float32'],
                      ### dtypes=['float64'],
                      time_units='minutes',
                      time_res='60.0',
                      comment=''):
              
        #----------------------------
        # Does file already exist ?
        #----------------------------
        file_name = file_utils.check_overwrite( file_name )
                
        #---------------------------------------
        # Check and store the time series info
        #---------------------------------------
        self.format     = 'ncts'
        self.file_name  = file_name
        self.time_index = 0
        if (long_names[0] is None):
            long_names = var_names

        #######################################################             
        # Assume for now that var_names only differ by the
        # appended row and column info, so only need one each
        # of svo_name, long_name and units_name.
        # First, strip trailing row and column numbers.
        #######################################################
        s  = var_names[0]
        p1 = s.rfind('_')
        s  = s[:p1]
        p2 = s.rfind('_')
        short_name = s[:p2]
        svo_name   = svo_names.get_svo_name( short_name )   
        long_name  = long_names[0]
        units_name = units_names[0]
        
        #-------------------------------------------
        # Need this to compute grid cell lat & lon
        #-------------------------------------------
        xres_deg = (grid_info.xres / 3600.0)
        yres_deg = (grid_info.yres / 3600.0)
        minlon   = grid_info.x_west_edge
        maxlon   = grid_info.x_east_edge
        minlat   = grid_info.y_south_edge
        maxlat   = grid_info.y_north_edge

        #-------------------------------------------
        # We may not need to save these in self.
        # I don't think they're used anywhere yet.
        #-------------------------------------------
        self.var_names   = var_names           
        self.long_names  = long_names
        self.units_names = units_names
        self.time_units  = time_units
        self.dtypes      = dtypes

        #---------------------------------------------
        # Create time metadata strings  (2020-01-14)
        #---------------------------------------------
        # str(datetime.datetime.now()) =
        #   '2020-01-14 12:35:32.087911'
        #-------------------------------------------------
        # x = datetime.datetime(2018, 9, 15, 12, 45, 35)
        # str(x) = '2018-09-15 12:45:35'
        #-------------------------------------------------
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
        #------------------------------------------------
        start_datetime = start_date + ' ' + start_time
        end_datetime   = end_date   + ' ' + end_time
        dur_units      = time_units        
        duration = tf_utils.get_duration( start_date, start_time,
                                          end_date, end_time,
                                          dur_units)
                                          
        #---------------------------------------------
        # Create array of dtype codes from dtypes
        # for multiple time series (i.e. columns).
        #---------------------------------------------
        dtype_codes = get_dtype_codes( dtypes, var_names )
        self.dtype_codes = dtype_codes
            
        #-------------------------------------
        # Open a new netCDF file for writing
        #-------------------------------------
        try:
            ## format = 'NETCDF4'
            format = 'NETCDF4_CLASSIC'
            ncts_unit = nc.Dataset(file_name, mode='w', format=format)
            OK = True
        except:
            OK = False
            return OK

        #------------------------------------------------------------
        # Option to pre-fill with fill values
        # Set fill_value for a var with "var._Fill_Value = number"
        # For Nio was:  opt.PreFill = False # (for efficiency)
        #------------------------------------------------------------
        ncts_unit.set_fill_off()
        # ncts_unit.set_fill_on()
         
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
        title = 'Time series data for variable: ' + long_name
        tf_version = str(tf_utils.TF_Version_Number())
        summary  = 'This file contains one or more time series for '
        summary += 'the single variable: ' + long_name + ', at '
        summary += 'model grid cells specified in an outlets file. '
        summary += 'Short var names have form:  symbol_row_col.'
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
        ncts_unit.title             = title
        ncts_unit.summary           = summary
        ncts_unit.comment           = comment
        ncts_unit.history           = history
        ncts_unit.creator_email     = email
        ncts_unit.date_created      = date_created 
        ncts_unit.naming_authority  = naming_authority   
        ncts_unit.geospatial_bounds_crs = '+init=epsg:4979'
        ## bounds = [minlon, minlat, maxlon, maxlat]   #### MINT order
        ## ncts_unit.geospatial_bounds = bounds 
                       
        #------------------------------------------------
        # Create an unlimited time dimension (via None)
        #------------------------------------------------
        # Without using "int()" for length, we get this:
        #     TypeError: size must be None or integer
        #------------------------------------------------
        ncts_unit.createDimension("time", None)

        #------------------------------------------
        # Save attributes of coordinate var, time
        #---------------------------------------------------
        #('f' = float32; must match in add_values_at_IDs()
        #---------------------------------------------------
        # NB! Can't use "time" vs. "tvar" here unless we
        #     add "import time" inside this function.
        #---------------------------------------------------
        tvar = ncts_unit.createVariable('time', 'f8', ("time",))
        ncts_unit.variables['time'].units = time_units
        ncts_unit.variables['time'].time_coverage_resolution = time_res    
        ncts_unit.variables['time'].time_coverage_start = start_datetime 
        ncts_unit.variables['time'].time_coverage_end = end_datetime 
        ncts_unit.variables['time'].time_coverage_duration = duration
                
        #-----------------------------------
        # Create variables using var_names
        #---------------------------------------------------
        # NB! The 3rd argument here (dimension), must be a
        #     tuple.  If there is only one dimension, then
        #     we need to add a comma, as shown.
        #---------------------------------------------------
        for k in range(len(var_names)):
            var_name = var_names[k]
            var = ncts_unit.createVariable(var_name, dtype_codes[k], ("time",))
        
            #-----------------------------------------
            # Create attributes of the main variable
            #-----------------------------------------
            # ncts_unit.variables[var_name].standard_name = standard_names[k] 
            # ncts_unit.variables[var_name].long_name = long_names[k]
            # ncts_unit.variables[var_name].units     = units_names[k] 
            #-------------------------------------------------------------
            ncts_unit.variables[var_name].svo_name  = svo_name             
            ncts_unit.variables[var_name].long_name = long_name
            ncts_unit.variables[var_name].units     = units_name       
            ncts_unit.variables[var_name].n_values  = 0   ##########
            #-------------------------------------------------------------
            # Compute & save geospatial info
            #----------------------------------------------------
            # NOTE:  var_name can have "_", so index from right
            #----------------------------------------------------
            ## print('var_name =', var_name)
            p   = var_name.split('_')
            row = np.int16( p[-2] ) 
            col = np.int16( p[-1] )
            lon = minlon + (col * xres_deg)
            lat = minlat + (row * yres_deg)
            ncts_unit.variables[var_name].geospatial_lon = lon
            ncts_unit.variables[var_name].geospatial_lat = lat            
            #----------------------------------------------------------------           
#         ncts_unit.variables[var_name].valid_min     = valid_min
#         ncts_unit.variables[var_name].valid_max     = valid_max
#         ncts_unit.variables[var_name].valid_range   = valid_range
#         ncts_unit.variables[var_name].missing_value = missing_value
#         ncts_unit.variables[var_name].fill_value    = fill_value  

            #----------------------------------
            # Specify a "nodata" fill value ?
            #----------------------------------
            # var._Fill_Value = -9999.0    ## Used for pre-fill above ?
            
        self.ncts_unit = ncts_unit
        return OK
    
    #   open_new_file()
    #----------------------------------------------------------
    def update_time_index(self, step=1): 

        #---------------------------------------------------
        # We shouldn't update clock in every add_value()
        # call because different values (for same time)
        # may be written with different add_value() calls.
        #---------------------------------------------------
        
        #------------------------------------
        # Increment the internal time index
        #------------------------------------
        self.time_index += step

    #   update_time_index()
    #----------------------------------------------------------
    def add_value(self, value, var_name, time=None,
                  time_index=-1):

        #---------------------------------------------------
        # Note: "time_index" allows insertion/overwrite
        #       of a value at a particular location.
        #---------------------------------------------------
        # This syntax works for scalars and grids
        # nc_unit.variables[var_name].assign_value( value )
        #---------------------------------------------------

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
        times = self.ncts_unit.variables[ 'time' ]
        times[ time_index ] = time
        
        #---------------------------------------------
        # Write a data value to existing netCDF file
        #---------------------------------------------
        values = self.ncts_unit.variables[ var_name ]
        values[ time_index ] = value
        self.ncts_unit.variables[ var_name ].n_values += 1
        
        ####################################################
        # We shouldn't update clock in every add_value()
        # call because different values (for same time)
        # may be written with different add_value() calls.
        ####################################################
        #------------------------------------
        # Increment the internal time index
        #------------------------------------
        # self.time_index += 1
        
        #-------------------------------------------------
        # 12/2/09:  netCDF is supposed to take care of
        # byteorder transparently.  However, we need to
        # make sure we don't byteswap in the function
        # "model_output.save_value_to_file()" when the
        # output format is netCDF.
        #-------------------------------------------------        
##        if (sys.byteorder == 'big'):
##            var[time_index] = value
##        else:
##            value2 = value.copy()
##            var[time_index] = value2.byteswap() 
##        self.time_index += 1
        
    #   add_value()
    #----------------------------------------------------------
    def get_value(self, var_name, time_index):

        values = self.ncts_unit.variables[ var_name ]
        times  = self.ncts_unit.variables[ 'time' ]
        return (values[ time_index ], times[ time_index ])
        
    #   get_value()
    #-------------------------------------------------------------------
    def values_at_IDs(self, var, IDs):

        #----------------------------------------------------------
        # Notes:  If "var" is a grid, subscript with self.IDs to
        #         get a 1D array of values.  If "var" is scalar,
        #         return a vector with the scalar value repeated
        #         once for each ID in self.IDs.
        #----------------------------------------------------------
        
        #---------------------------------
        # Is variable a grid or scalar ?
        #---------------------------------
        if (np.ndim(var) > 0):
            return np.float32( var[ IDs ] )
        else:
            #-----------------------------------------------------
            # (3/16/07) Bug fix.  This gets used in case of q0,
            # which is a scalar when INFIL_ALL_SCALARS is true.
            # Without this, don't get a value for every ID.
            #-----------------------------------------------------
            n_IDs  = np.size(IDs[0])
            vector = np.zeros( n_IDs, dtype='Float32')
            return (vector + np.float32(var)) 
        
    #   values_at_IDs()
    #-------------------------------------------------------------------
    def add_values_at_IDs(self, time, var, var_name, IDs,
                          time_index=-1):

        #---------------------------------------------------
        # Note: Here "var" is typically a grid and IDs are
        #       (row,col) subscripts into the grid.  A set
        #       of variable names are constructed from the
        #       actual "var_name" (e.g. "Q") and the
        #       row and column.  Note that we must have
        #       called open_new_file() with these same
        #       var_names.
        #---------------------------------------------------
        # Note: "time_index" allows insertion/overwrite
        #       of a value at a particular location.
        #---------------------------------------------------
        # This syntax works for scalars and grids
        # nc_unit.variables[var_name].assign_value( value )
        #---------------------------------------------------

        #-------------------------------------
        # Can use time_index to overwrite an
        # existing grid vs. simple append.
        #-------------------------------------
        if (time_index == -1):
            time_index = self.time_index

        #---------------------------------------------
        # Write current time to existing netCDF file
        #---------------------------------------------
        times = self.ncts_unit.variables[ 'time' ]
        times[ time_index ] = time
        
        #--------------------------------------------
        # Write data values to existing netCDF file
        #--------------------------------------------
        vals  = self.values_at_IDs( var, IDs )
        rows  = IDs[0]
        cols  = IDs[1]
        n_IDs = np.size(rows)
        for k in range(n_IDs):
            #----------------------------------------
            # Construct var_name of form:  Q[24,32]
            # or, if necessary, Q_24_32
            #----------------------------------------
            row_str  = '_' + str(rows[k])
            col_str  = '_' + str(cols[k])
            #--------------------------------------------------
            # Must match with model_output.open_new_ts_file()
            #--------------------------------------------------
            ## row_str = '[' + str(rows[k]) + ','
            ## col_str = str(cols[k]) + ']'
            
            vname  = var_name + row_str + col_str
            values = self.ncts_unit.variables[ vname ]
            values[ time_index ] = vals[k]
            values.n_values += 1
        
        #---------------------------
        # Increment the time index
        #---------------------------
        self.time_index += 1
         
    #   add_values_at_IDs() 
    #-------------------------------------------------------------------
    def add_series(self, values, var_name, times,
                   time_index=-1):

        #-----------------------------------------------------
        # Note: "time_index" allows insertion/overwrite
        #       of a time series at a particular location.
        #-----------------------------------------------------
        # This syntax works for scalars and grids
        # nc_unit.variables[var_name].assign_value( values )
        #-----------------------------------------------------


        #-------------------------------------
        # Can use time_index to overwrite an
        # existing grid vs. simple append.
        #-------------------------------------
        if (time_index == -1):
            time_index = self.time_index

        #---------------------------------------------
        # Write a data value to existing netCDF file
        #---------------------------------------------
        series = self.ncts_unit.variables[ var_name ]
        series[:] = values

        ######################################################
        # WE SHOULDN'T update clock in every add_value()
        # call because different vars (e.g. the time)
        # must be written with different add_value() calls.
        ######################################################
        #------------------------------------
        # Increment the internal time index
        #------------------------------------
        # self.time_index += np.size(values)
        
    #   add_series()
    #----------------------------------------------------------
    def get_var_names(self):
    
        var_dict = self.ncts_unit.variables
        return list( var_dict.keys() )

    #   get_var_names()
    #----------------------------------------------------------
    def get_var_long_name(self, var_name ):

        var = self.ncts_unit.variables[ var_name ]
        return var.long_name 
            
    #   get_var_long_name()
    #----------------------------------------------------------
    def get_var_units(self, var_name ):

        var = self.ncts_unit.variables[ var_name ]
        return var.units

    #   get_var_units()
    #----------------------------------------------------------
    def get_var_lons(self):
        
        var_names = self.get_var_names()
        var_names = var_names[1:]    # exclude 'time'
        lons = []
        for name in var_names:
            var = self.ncts_unit.variables[ name ]
            lons.append( var.geospatial_lon )
        return lons
 
    #   get_var_lons()
    #----------------------------------------------------------
    def get_var_lats(self):
        
        var_names = self.get_var_names()
        var_names = var_names[1:]    # exclude 'time'
        lats = []
        for name in var_names:
            var = self.ncts_unit.variables[ name ]
            lats.append( var.geospatial_lat )
        return lats
 
    #   get_var_lats()
    #----------------------------------------------------------
    def get_series(self, var_name):

        #----------------------------------------------------------
        # Note:  Attributes such as long_name and units can
        #        be obtained as series.units or series.long_name.
        #        The actual values can be obtained from:
        #        values = series[:], or np.array(series)
        #----------------------------------------------------------
        series = self.ncts_unit.variables[ var_name ]
        times  = self.ncts_unit.variables[ 'time' ]
        return (series, times)
               
    #   get_series()
    #-------------------------------------------------------------------
    def close_file(self):

        # self.ncts_unit.sync()  ## (netCDF4 has no "flush")
        self.ncts_unit.close()

    #   close_file()
    #-------------------------------------------------------------------
    def close(self):

        # self.ncts_unit.sync()  ## (netCDF4 has no "flush")
        self.ncts_unit.close()

    #   close()
    #-------------------------------------------------------------------
    
