
# S.D. Peckham
# May 2010

import os
import sys
import time

import numpy as np
from . import file_utils
# import rti_files   # (not used)

import netCDF4 as nc

#---------------------------------------------------------------------
# This class is for I/O of time-indexed 1D profiles to netCDF files.
#---------------------------------------------------------------------
#
#   unit_test()
#   unit_test2()
#   save_as_text()   # (not ready yet)
#
#   class ncps_file():
#
#       import_netCDF4()
#       open_file()
#       get_dtype_map()
#       open_new_file()
#       update_time_index()
#-------------------------------
#       add_profile()
#       get_profile()
#-------------------------------
#       profiles_at_IDs()
#       add_profiles_at_IDs()
#-------------------------------
#       close_file()
#       close()
#
#-------------------------------------------------------------------
def unit_test(n_times=5, nz=10, VERBOSE=False,
              file_name="NCPS_Profile_Test.nc"):

    #--------------------------------------------------------
    # Notes: This test uses add_profile() and get_profile()
    #        to add and retrieve a time-indexed set of 1D
    #        profiles to/from a file.  An example would be
    #        a set of soil-moisture profiles that vary with
    #        both depth, z, and time.
    #--------------------------------------------------------
    print(' ')
    print('Running unit_test()...')

    #-------------------------------------
    # Make instance of ncps_file() class
    #-------------------------------------
    ncps = ncps_file()
    var_names = ['theta']
    z_values = np.arange( nz, dtype='Float64' )
    z_units  = 'm'
    
    OK = ncps.open_new_file( file_name,
                             z_values=z_values,
                             z_units=z_units,
                             var_names=var_names,
                             long_names=['soil_water_content'],
                             units_names=['none'],
                             dtypes=['float64'],
                             time_units='minutes',
                             comment="Created by TopoFlow 3.0.")
                          
    ###############################################
    # WHAT ABOUT LONG_NAME for the TIME VALUES ??
    ###############################################
    
    if not(OK):
        print('ERROR during open_new_file().')
        return

    profile  = np.exp(-0.1 * z_values)
    times    = np.arange( n_times, dtype='Float64') * 0.1
    
    #-----------------------------------
    # Add a series of profiles to file
    #-----------------------------------
    print('Writing profiles to ncps file...')
    for time_index in range(n_times):
        time  = times[ time_index ]
        ncps.add_profile( profile, var_names[0], time )
        #----------------------------------------------
        ncps.update_time_index()
        profile += 1    ## (make profile change in time)
    if (VERBOSE):
        print(self.ncps_unit)  # (print a summary)

    ncps.close_file()
    print('Finished writing ncps file: ' + file_name)
    print(' ')

    #-----------------------------------------
    # Re-open the file and read the profiles
    #-----------------------------------------
    OK = ncps.open_file( ncps.file_name )
    if not(OK): return
    print('Reading values from ncps file: ')
    
    for time_index in range(n_times):
        profile, time = ncps.get_profile(var_names[0], time_index)

        ti_str = str(time_index)
        print('time[' + ti_str + '] =', time)
        print('profile[' + ti_str + '] =', profile)
        print('-----------------------------------------------')

    #-----------------
    # Close the file
    #-----------------
    ncps.close_file()    
    print('Finished reading ncps file: ' + file_name)
    print(' ')
    
#   unit_test()
#-------------------------------------------------------------------
def unit_test2(n_times=5, nz=10, VERBOSE=False,
               file_name="NCPS_Profile_Test2.nc"):

    #--------------------------------------------------------
    # Notes: This test uses add_profile() and get_profile()
    #        to add and retrieve a time-indexed set of 1D
    #        profiles to/from a file.  An example would be
    #        a set of soil-moisture profiles that vary with
    #        both depth, z, and time.
    #--------------------------------------------------------
    print(' ')
    print('Running unit_test2()...')

    #-------------------------------------
    # Make instance of ncps_file() class
    #-------------------------------------
    ncps = ncps_file()

    var_name   = 'theta'
    z_values   = np.arange( nz, dtype='Float64' )
    z_units    = 'm'
    IDs        = ([1,2,3], [1,2,3])
    var_names  = ['theta_1_1', 'theta_2_2', 'theta_3_3']
    long_names = ['soil_water_content_profile_at_1_1',
                  'soil_water_content_profile_at_2_2',
                  'soil_water_content_profile_at_3_3']
    units_names = ['none', 'none', 'none']
    dtypes      = ['float64']
    # dtypes      = ['float64', 'float64', 'float64']
    
    OK = ncps.open_new_file( file_name,
                             z_values=z_values,
                             z_units=z_units,
                             var_names=var_names,
                             long_names=long_names,
                             units_names=units_names,
                             dtypes=dtypes,
                             time_units='minutes',
                             comment="Created by TopoFlow 3.0.")
                          
    ###############################################
    # WHAT ABOUT LONG_NAME for the TIME VALUES ??
    ###############################################
    
    if not(OK):
        print('ERROR during open_new_file().')
        return

    profile = np.exp(-0.1 * z_values)
    times   = np.arange( n_times, dtype='Float64') * 0.1
    print('z_values =', z_values)
    print('profile  =', profile)
    print('times    =', times)
    print(' ')
    
    ny = 5
    nx = 5
    var = np.zeros([nz,ny,nx], dtype='float64')
    for k in range(nz):
        var[k,:,:] = profile[k]
    
    #-----------------------------------
    # Add a series of profiles to file
    #-----------------------------------
    print('Writing profiles to ncps file...')
    for time_index in range(n_times):
        time  = times[ time_index ]
        ncps.add_profiles_at_IDs(var, var_name, IDs, time )
        #-------------------------------------------------
        # Don't need to update_time_index, done already.
        #-------------------------------------------------
        #### ncps.update_time_index()
        var += 1    ## (make profiles change in time)
    if (VERBOSE):
        print(self.ncps_unit)  # (print a summary)

    ncps.close_file()
    print('Finished writing ncps file: ' + file_name)
    print(' ')
    
    #-----------------------------------------
    # Re-open the file and read the profiles
    #-----------------------------------------
    OK = ncps.open_file( ncps.file_name )
    if not(OK): return
    print('Reading values from ncps file: ')
    
    for time_index in range(n_times):
        profile, time = ncps.get_profile(var_names[0], time_index)

        ti_str = str(time_index)
        print('time[' + ti_str + '] =', time)
        print('profile[' + ti_str + '] =', profile)
        print('-----------------------------------------------')

    #-----------------
    # Close the file
    #-----------------
    ncps.close_file()    
    print('Finished reading ncps file: ' + file_name)
    print(' ')
    
#   unit_test2()
#-------------------------------------------------------------------
def save_as_text(ncps_file_name=None, text_file_name=None):

    ncps = ncps_file()
    OK = ncps.open_file( ncps_file_name )
    if not(OK): return

    var_name  = 'theta'
    data = ncps.get_profile( var_name )
    ncps.close()
    
    data = np.array( data )
    print('min(data), max(data) =', data.min(), data.max())

    text_unit = open( text_file_name, 'w' )
    data.tofile( unit )  ###### CHECK THIS #######
    text_unit.close()

#   save_as_text()
#-------------------------------------------------------------------
class ncps_file():

    #----------------------------------------------------------
    # Note:  ncps = NetCDF Time Series (used by CSDMS)
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
            ncps_unit = nc.Dataset(file_name, mode='r')
            self.ncps_unit = ncps_unit
            ### return ncps_unit
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
                      z_values=np.arange(10),
                      z_units='m',
                      var_names=['X'],
                      long_names=[None],
                      units_names=['None'],
                      dtypes=['float64'],
                      time_units='minutes',
                      comment=''):

        #----------------------------------------------------
        # Notes: It might be okay to have "nz" be an
        #        unlimited dimension, like "time".  This
        #        would mean replacing "int(profile_length)"
        #        with "None".
        #----------------------------------------------------

        #----------------------------
        # Does file already exist ?
        #----------------------------
        file_name = file_utils.check_overwrite( file_name )
        self.file_name = file_name
        
        #---------------------------------------
        # Check and store the time series info
        #---------------------------------------
        self.format     = 'ncps'
        self.file_name  = file_name
        self.time_index = 0
        if (long_names[0] is None):
            long_names = var_names
        #-------------------------------------------            
        self.z_values  = z_values
        self.z_units   = z_units
        nz             = np.size(z_values)
        #-------------------------------------------
        # We may not need to save these in self.
        # I don't think they're used anywhere yet.
        #-------------------------------------------
        self.var_names   = var_names 
        self.long_names  = long_names
        self.units_names = units_names
        self.dtypes      = dtypes

        #---------------------------------------------
        # Create array of Nio type codes from dtypes
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
            ncps_unit = nc.Dataset(file_name, mode='w', format=format)
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
        ncps_unit.history = history
        # print 'MADE IT PAST history BLOCK'
       
        #------------------------------------------------
        # Create an unlimited time dimension (via None)
        #------------------------------------------------
        # Without using "int()" here, we get this:
        #     TypeError: size must be None or integer
        #------------------------------------------------
        ncps_unit.createDimension('nz', int(nz))
        ncps_unit.createDimension('time', None)

        #-------------------------
        # Create a time variable
        #---------------------------------------------------
        #('f' = float32; must match in add_values_at_IDs()
        #---------------------------------------------------
        # NB! Can't use "time" vs. "tvar" here unless we
        #     add "import time" inside this function.
        #---------------------------------------------------
        tvar = ncps_unit.createVariable('time', 'f8', ('time',))
        ncps_unit.variables['time'].units = time_units

        #--------------------------------------
        # Create a distance/depth variable, z
        #--------------------------------------
        zvar = ncps_unit.createVariable('z', 'd', ('z',))
        zvar[ : ] = z_values  # (store the z-values)
        ncps_unit.variables['z'].units = z_units
        
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
            var = ncps_unit.create_variable(var_name, dtype_codes[k],
                                            ("time", "nz"))
        
            #------------------------------------
            # Create attributes of the variable
            #------------------------------------
            ncps_unit.variables[var_name].long_name = long_names[k]
            ncps_unit.variables[var_name].units     = units_names[k]        

            #----------------------------------
            # Specify a "nodata" fill value ?
            #----------------------------------
            # var._Fill_Value = -9999.0    ## Used for pre-fill above ?
            
        self.ncps_unit = ncps_unit
        return OK
    
    #   open_new_file()
    #----------------------------------------------------------
    def update_time_index(self, step=1): 

        #-----------------------------------------------------
        # We shouldn't update clock in every add_profile()
        # call because different profiles (for same time)
        # may be written with different add_profile() calls.
        #-----------------------------------------------------
        
        #------------------------------------
        # Increment the internal time index
        #------------------------------------
        self.time_index += step

    #   update_time_index()
    #----------------------------------------------------------
    def add_profile(self, profile, var_name, time=None,
                    time_index=-1):
    
        #-----------------------------------------------------
        # Note: "time_index" allows insertion/overwrite
        #       of a profile at a particular time index.
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
        if (time is None):
            time = np.float64( time_index )

        #---------------------------------------------
        # Write a data value to existing netCDF file
        #---------------------------------------------
        profiles = self.ncps_unit.variables[ var_name ]
        profiles[ time_index ] = profile
        #------------------------------------------------
        times = self.ncps_unit.variables[ 'time' ]
        times[ time_index ] = time

        ######################################################
        # We shouldn't update clock in every add_profile()
        # call because different profiles (for same time)
        # may be written with different add_profile() calls.
        ######################################################
        #------------------------------------
        # Increment the internal time index
        #------------------------------------
        # self.time_index += np.size(values)
        
    #   add_profile()
    #----------------------------------------------------------
    def get_profile(self, var_name, time_index):

        profiles = self.ncps_unit.variables[ var_name ]
        times    = self.ncps_unit.variables[ 'time' ]
        return (profiles[ time_index ], times[ time_index ])
        
    #   get_profile()
    #-------------------------------------------------------------------
    def profiles_at_IDs(self, var, IDs):

        #---------------------------
        # Get the dimensions, etc.
        #---------------------------
        ndims = np.ndim(var)
        dtype = self.dtypes[0]
        nz    = var.shape[0]
        # nz  = np.size(var, 0)   # (also works)
        n_IDs = np.size(IDs[0])
        profiles = np.zeros([n_IDs, nz], dtype=dtype)
 
        if (ndims == 1):    
            #------------------------------
            # Variable is a 1D profile,
            # and is the same for all IDs
            #------------------------------
            for k in range(n_IDs):
                profiles[k, :] = var.astype(dtype)
        else:    
            #---------------------------------
            # Variable is a 3D array; return
            # a profile for each ID
            #---------------------------------         
            for k in range(nz):
                layer = var[k,:,:].astype(dtype)
                profiles[:, k] = layer[ IDs ]
                
        return profiles

    #   profiles_at_IDs()
    #-------------------------------------------------------------------
    def add_profiles_at_IDs(self, var, var_name, IDs, time=None,
                            time_index=-1):

        #-------------------------------------
        # Can use time_index to overwrite an
        # existing grid vs. simple append.
        #-------------------------------------
        if (time_index == -1):
            time_index = self.time_index
        if (time is None):
            time = np.float64( time_index )
            
        #---------------------------------------------
        # Write current time to existing netCDF file
        #---------------------------------------------
        times = self.ncps_unit.variables[ 'time' ]
        times[ time_index ] = time
        
        #--------------------------------------------
        # Write data values to existing netCDF file
        #--------------------------------------------
        profiles = self.profiles_at_IDs( var, IDs )
        rows     = IDs[0]
        cols     = IDs[1]
        n_IDs    = np.size(rows)
        for k in range(n_IDs):
            #----------------------------------------
            # Construct var_name of form:  Q[24,32]
            # or, if necessary, Q_24_32
            #----------------------------------------
            row_str  = '_' + str(rows[k])
            col_str  = '_' + str(cols[k])
            #--------------------------------------------------
            # Must match with model_output.open_new_ps_file()
            #--------------------------------------------------
            ## row_str = '[' + str(rows[k]) + ','
            ## col_str = str(cols[k]) + ']'
            
            vname  = var_name + row_str + col_str
            profile_series = self.ncps_unit.variables[ vname ]
            profile_series[ time_index ] = profiles[k,:]

            ## print 'added profile =', profiles[k,:]  ###########
            
        #---------------------------
        # Increment the time index
        #---------------------------
        self.time_index += 1

    #   add_profiles_at_IDs()
    #-------------------------------------------------------------------
    def close_file(self):

        # self.ncts_unit.sync()  ## (netCDF4 has no "flush")
        self.ncps_unit.close()

    #   close_file()
    #-------------------------------------------------------------------
    def close(self):

        # self.ncts_unit.sync()  ## (netCDF4 has no "flush")
        self.ncps_unit.close()

    #   close()
    #-------------------------------------------------------------------
    
