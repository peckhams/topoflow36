
# This new verion (6/11/10) hasn't been tested yet.
# Can't run unit tests on my MacPro w/o Nio.
#---------------------------------------------------

# S.D. Peckham
# Sept 2014 (new version to use netCDF4)
# May, June 2010

import os
import sys
import time

import numpy as np
from . import file_utils

import netCDF4 as nc

#-------------------------------------------------------------------
# This class is for I/O of time series data to netCDF files.
#-------------------------------------------------------------------
#
#   unit_test1()
#   unit_test2()
#   save_as_text()   # (not ready yet)
#
#   class ncts_file():
#
#       import_netCDF4()
#       open_file()
#       get_dtype_map()
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
    var_names = ['depth']

    OK = ncts.open_new_file( file_name,
                             var_names=var_names,
                             long_names=["depth of water"],
                             units_names=["meters"],
                             dtypes=['float32'],
                             comment="Created by TopoFlow 3.0.")
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
    def open_new_file(self, file_name,
                      var_names=['X'],
                      long_names=[None],
                      units_names=['None'],
                      dtypes=['float32'],
                      ### dtypes=['float64'],
                      time_units='minutes',
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
        # Create array of dtype codes from dtypes
        # for multiple time series (i.e. columns).
        #---------------------------------------------
        dtype_map   = self.get_dtype_map()
        dtype_codes = []
        if (len(dtypes) == len(var_names)):
            for dtype in dtypes:
               dtype_code = dtype_map[ dtype.lower() ]
               dtype_codes.append( dtype_code )
        else:
            dtype = dtypes[0]
            dtype_code = dtype_map[ dtype.lower() ]
            for k in range(len(var_names)):
                dtype_codes.append( dtype_code )                
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
        history = history + comment
        ncts_unit.history = history
                
        #------------------------------------------------
        # Create an unlimited time dimension (via None)
        #------------------------------------------------
        # Without using "int()" for length, we get this:
        #     TypeError: size must be None or integer
        #------------------------------------------------
        ncts_unit.createDimension("time", None)

        #-------------------------
        # Create a time variable
        #---------------------------------------------------
        #('f' = float32; must match in add_values_at_IDs()
        #---------------------------------------------------
        # NB! Can't use "time" vs. "tvar" here unless we
        #     add "import time" inside this function.
        #---------------------------------------------------
        tvar = ncts_unit.createVariable('time', 'f8', ("time",))
        ncts_unit.variables['time'].units = time_units
        
        #-----------------------------------
        # Create variables using var_names
        #-----------------------------------
        # Returns "var" as a PyNIO object
        #---------------------------------------------------
        # NB! The 3rd argument here (dimension), must be a
        #     tuple.  If there is only one dimension, then
        #     we need to add a comma, as shown.
        #---------------------------------------------------
        for k in range(len(var_names)):
            var_name = var_names[k]
            var = ncts_unit.createVariable(var_name, dtype_codes[k], ("time",))
        
            #------------------------------------
            # Create attributes of the variable
            #------------------------------------
            ncts_unit.variables[var_name].long_name = long_names[k]
            ncts_unit.variables[var_name].units     = units_names[k]        

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
    
