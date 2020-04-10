
# S.D. Peckham
# May 2010
# Jan 2020.   Added new "MINT netCDF" metadata.

import os
import sys
import time
import datetime

import numpy as np
from . import file_utils
# from . import rti_files   # (not used for unit_test() yet.
from . import svo_names
from . import tf_utils

import netCDF4 as nc

#---------------------------------------------------------------------
# This class is for I/O of time-indexed 1D profiles to netCDF files.
#---------------------------------------------------------------------
#
#   unit_test()
#   unit_test2()
#
#   save_as_text()   # (not ready yet)
#   get_dtype_map()
#   get_dtype_coes()      # 2020-01-26  (separate function)
#
#   class ncps_file():
#
#       import_netCDF4()
#       open_file()

#       open_new_file()
#       update_time_index()
#-------------------------------
#       add_profile()
#       get_profile()
#       profiles_at_IDs()
#       add_profiles_at_IDs()
#--------------------------------------------
#       get_var_names()        # 2020-01-26
#       get_var_long_name()    # 2020-01-26
#       get_var_units()        # 2020-01-26
#       get_var_lons()         # 2020-01-26
#       get_var_lats()         # 2020-01-26
#--------------------------------------------
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
                             comment="Created by TopoFlow 3.6.")
                          
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
                             comment="Created by TopoFlow 3.6.")
                          
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
            return False
        
    #   import_netCDF4()
    #----------------------------------------------------------
    def open_file(self, file_name):
      
        #-------------------------
        # Open file to read only
        #-------------------------
#         try:
        ncps_unit = nc.Dataset(file_name, mode='r')
        self.ncps_unit = ncps_unit
        ### return ncps_unit
        return True
#         except:
#             print('ERROR: Could not open file:')
#             print( '   ' + file_name )
#             print( 'Current working directory =')
#             print( '   ' + os.getcwd() )
#             return False
    
    #   open_file()
    #----------------------------------------------------------
    def open_new_file(self, file_name,
                      grid_info=None,
                      time_info=None,
                      z_values=np.arange(10),
                      z_units='m',
                      var_names=['q_2_3'],
                      long_names=[None],
                      units_names=['None'],
                      dtypes=['float32'],
                      ## dtypes=['float64'],
                      time_units='minutes',
                      time_res='60.0',
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
        ## print('######## dtype_codes =', dtype_codes)
           
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
        ncps_unit.set_fill_off()
        # ncps_unit.set_fill_on()
        
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
        title = 'Profile series data for variable: ' + long_name
        tf_version = str(tf_utils.TF_Version_Number())
        summary  = 'This file contains one or more profile series for '
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
        ncps_unit.title             = title
        ncps_unit.summary           = summary
        ncps_unit.comment           = comment
        ncps_unit.history           = history
        ncps_unit.creator_email     = email
        ncps_unit.date_created      = date_created 
        ncps_unit.naming_authority  = naming_authority   
        ncps_unit.geospatial_bounds_crs = '+init=epsg:4979'
        ## bounds = [minlon, minlat, maxlon, maxlat]   #### MINT order
        ## ncps_unit.geospatial_bounds = bounds 
               
        #------------------------------------------------
        # Create an unlimited time dimension (via None)
        #------------------------------------------------
        # Without using "int()" here, we get this:
        #     TypeError: size must be None or integer
        #------------------------------------------------
        ncps_unit.createDimension('time', None)

        #------------------------------------------
        # Save attributes of coordinate var, time
        #---------------------------------------------------
        #('f' = float32; must match in add_values_at_IDs()
        #---------------------------------------------------
        # NB! Can't use "time" vs. "tvar" here unless we
        #     add "import time" inside this function.
        #---------------------------------------------------
        tvar = ncps_unit.createVariable('time', 'f8', ("time",))
        ncps_unit.variables['time'].units = time_units
        ncps_unit.variables['time'].time_coverage_resolution = time_res    
        ncps_unit.variables['time'].time_coverage_start = start_datetime 
        ncps_unit.variables['time'].time_coverage_end = end_datetime 
        ncps_unit.variables['time'].time_coverage_duration = duration

        #--------------------------------------
        # Create an "z" dimension.
        # Create a distance/depth variable, z
        #--------------------------------------
        ncps_unit.createDimension('z', int(nz))
        zvar = ncps_unit.createVariable('z', 'f4', ('z',))
        zvar[ : ] = z_values  # (store the z-values)
        ncps_unit.variables['z'].units = z_units
        
        #-----------------------------------
        # Create variables using var_names
        #---------------------------------------------------
        # NB! The 3rd argument here (dimension), must be a
        #     tuple.  If there is only one dimension, then
        #     we need to add a comma, as shown.
        #---------------------------------------------------
        for k in range(len(var_names)):
            var_name = var_names[k]
            var = ncps_unit.createVariable(var_name, dtype_codes[k],
                                            ("time", "z"))

            #-----------------------------------------
            # Create attributes of the main variable
            #-----------------------------------------
            # ncps_unit.variables[var_name].standard_name = standard_names[k] 
            # ncps_unit.variables[var_name].long_name = long_names[k]
            # ncps_unit.variables[var_name].units     = units_names[k] 
            #-------------------------------------------------------------
            ncps_unit.variables[var_name].svo_name  = svo_name             
            ncps_unit.variables[var_name].long_name = long_name
            ncps_unit.variables[var_name].units     = units_name       
            ncps_unit.variables[var_name].n_profiles = 0
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
            ncps_unit.variables[var_name].geospatial_lon = lon
            ncps_unit.variables[var_name].geospatial_lat = lat            
            #----------------------------------------------------------------           
#         ncps_unit.variables[var_name].valid_min     = valid_min
#         ncps_unit.variables[var_name].valid_max     = valid_max
#         ncps_unit.variables[var_name].valid_range   = valid_range
#         ncps_unit.variables[var_name].missing_value = missing_value
#         ncps_unit.variables[var_name].fill_value    = fill_value
            
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
        z        = self.ncps_unit.variables[ 'z' ]
        
        return (profiles[ time_index ], z, times[ time_index ])
        
    #   get_profile()
    #----------------------------------------------------------
    def get_profiles(self, var_name):

        profiles = self.ncps_unit.variables[ var_name ]
        times    = self.ncps_unit.variables[ 'time' ]
        z        = self.ncps_unit.variables[ 'z' ]
        
        return (profiles, z, times)
        
    #   get_profiles()
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
            profile_series.n_profiles += 1
            
            ## print 'added profile =', profiles[k,:]  ###########
            
        #---------------------------
        # Increment the time index
        #---------------------------
        self.time_index += 1

    #   add_profiles_at_IDs()
    #----------------------------------------------------------
    def get_var_names(self):
    
        var_dict = self.ncps_unit.variables
        return list( var_dict.keys() )

    #   get_var_names()
    #----------------------------------------------------------
    def get_var_long_name(self, var_name ):

        var = self.ncps_unit.variables[ var_name ]
        return var.long_name 
            
    #   get_var_long_name()
    #----------------------------------------------------------
    def get_var_units(self, var_name ):

        var = self.ncps_unit.variables[ var_name ]
        return var.units

    #   get_var_units()
    #----------------------------------------------------------
    def get_var_lons(self):
        
        var_names = self.get_var_names()
        var_names = var_names[2:]    # exclude 'time' & 'z'
        lons = []
        for name in var_names:
            var = self.ncps_unit.variables[ name ]
            lons.append( var.geospatial_lon )
        return lons
 
    #   get_var_lons()
    #----------------------------------------------------------
    def get_var_lats(self):
        
        var_names = self.get_var_names()
        var_names = var_names[2:]    # exclude 'time' & 'z'
        lats = []
        for name in var_names:
            var = self.ncps_unit.variables[ name ]
            lats.append( var.geospatial_lat )
        return lats
 
    #   get_var_lats()
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
    
