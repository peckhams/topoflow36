#
#  We should use "initialize_scalar()" for all scalar assignments.
#  See the Notes for that method.  Search for "np.float64(0".
#
#-----------------------------------------------------------------------      
#  Copyright (c) 2009-2023, Scott D. Peckham
#
#  Sep 2023. Improved how read_config_file() handles boolean vars.
#            Added save_edge_ids(); from d8_base.py/get_edge_ids();
#            called from read_grid_info().
#  Nov 2022. Minor changes to support NOAA NextGen.
#  Oct 2022. Full compliance with BMI v2.1.  Cleaned up.
#  May 2020. Revamped version of set_directories()
#  Sep 2019. Updates to allow TF to run in Python 3, including
#            print() and exec() functions.  For exec(), see:
#            https://docs.python.org/3.1/library/functions.html#exec 
#            https://docs.python.org/3/reference/executionmodel.html
#            Optional 2nd argument is a dictionary of globals, such
#               as the one returned by the builtin globals().
#            Optional 3rd argument (seems to require 2nd) is a
#               dictionary of locals, such as returned by locals().
#            I pass an empty dictionary "{}" for 2nd arg in most
#               cases to emphasize when globals aren't needed.
#            NB!  In Py3.x, you can't use exec() to set locals in a
#            function, unless the name was already assigned to in that function.
#            However, in most of these cases we can use "eval()" instead.
#            And it is okay to set a variable into "self" directly.
#  Nov 2016. Added new BMI version 2 functions.
#            In get_grid_shape(),   ordering is now [nz, ny, nx].
#            In get_grid_spacing(), ordering is now [dz, dy, dx].
#            In get_grid_origin(),  ordering is now [z0, y0, x0].
#  Sep 2014. New initialize_basin_vars(), using outlets.py.
#            Removed obsolete functions.
#  Jan 2013. Added "initialize_scalar()" method.
#  Feb 2012. Complete BMI, starting from CSDMS_base.py
#            This now takes the place of CSDMS_base.py.
#  Nov 2011. Cleanup and conversion to BMI function names.
#  May 2010. initialize_config_vars(), cleanup, etc.
#  Aug 2009. Created, as CSDMS_base.py.
#
#-----------------------------------------------------------------------
#
#  Notes:  This file defines a "base class" with a BMI (Basic Model
#          Interface) for CSDMS "process" components.  These methods
#          allow a component to "fit into" a CMI harness (IMPL file)
#          that allows it to be used in the CSDMS/CCA framework.
#
#          The BMI interface is designed to be completely framework
#          independent, so this class contains no methods related to
#          CCA framework concepts such as ports.
#
#          Some "private" utility methods are defined at the end.
#
#-----------------------------------------------------------------------
#
#  unit_test()
#
#  class BMI_component
#
#      __init__()
#
#      -------------------------------------------------
#      BMI Functions (Many will need to be overridden)
#      -------------------------------------------------
#      get_bmi_version()
#      initialize()
#      update()          # override
#      update_frac()     # not BMI function
#      update_until()
#      finalize()
#      --------------------------
#      get_component_name()
#      get_input_item_count()
#      get_output_item_count()
#      get_input_var_names()      # override
#      get_output_var_names()     # override
#      --------------------------
#      get_var_grid()             # override
#      get_var_type()
#      get_var_units()            # override
#      get_var_itemsize()
#      get_var_nbytes()
#      get_var_location()
#      --------------------------
#      get_current_time()
#      get_start_time()
#      get_end_time()
#      get_time_units()
#      get_time_step()
#      --------------------------
#      get_value()
#      get_value_ptr()
#      get_value_at_indices()
#      set_value()
#      set_value_at_indices()
#      --------------------------
#      get_grid_type()
#      get_grid_rank()
#      get_grid_size()
#      get_grid_shape()
#      get_grid_spacing()
#      get_grid_origin()
#      get_grid_x()
#      get_grid_y()
#      get_grid_z()

#      ---------------------------
#      Non-BMI Utility Functions
#      ---------------------------
#      get_var_name()            # (override;  Not part of BMI)
#      get_status()              # Not part of BMI
#      get_attribute()           # Deprecated in BMI version 2
#      get_grid_attribute()      # NEW ADDITION TO BMI ??
#      ------------------------
#      run_model()               # (not required for BMI)
#      initialize_time_vars()
#      update_time()
#      check_finished()
#      -------------------------
#      print_time_and_value()
#      get_run_time_string()
#      print_run_time()
#      print_final_report()          # (6/30/10)
#      print_traceback()             # (10/10/10)
#      -------------------------
#      read_grid_info()              # read grid info from RTI file
#      save_edge_ids()               # (9/11/23)
#      read_path_info()              # (2/12/17)
#      read_time_info()              # (1/14/20)
#      read_config_file()            # (5/17/10, 5/9/11)
#      initialize_config_vars()      # (5/6/10)
#      set_computed_input_vars       # (5/6/10) over-ridden by each comp.
#      initialize_basin_vars()       # (9/19/14) New version that uses outlets.py.
#      -------------------------
#      prepend_directory()           # (may not work yet)
#      set_directories()             # revamped, renamed (5/3/20)
#      -------------------------
#      initialize_scalar()           # (2/5/13, for ref passing)
#      initialize_grid()
#      initialize_var()
#      update_var()                  # (11/15/16)
#      is_scalar()
#      is_vector()
#      is_grid()
#
#-----------------------------------------------------------------------

import numpy as np
import os, os.path
import sys
import time
import traceback        # (10/10/10)

#--------------------------------------------
# (5/14/10. Can't be here because basins.py
# has import BMI_base at top (circular).
# See initialize_basin_vars() below.
#--------------------------------------------
# import basins

from . import outlets          ## (9/19/14)
from . import pixels
from . import rti_files
from . import tf_utils
from . import time_utils

#-----------------------------------------------------------------------
def unit_test():

    c = BMI_component()
    print('Instantiated BMI component.')

#   unit_test()
#-----------------------------------------------------------------------
# def test1():
# 
#     c = BMI_component()
#     var_name = 'T_air_file'
#     cmd = 'c.' + var_name + " = ''"
#     exec( cmd, {}, locals() )
#     print( "c.T_air_file = " + c.T_air_file )
#     print( "locals() = ")
#     print( locals() )
#     ## print( "globals() = ")
#     ## print( globals() )    
#     
# #   test1()
#-----------------------------------------------------------------------
class BMI_component:

    def __init__(self):

        self.SILENT      = True    # (new default: 11/16/11)        
        ## self.DEBUG    = True
        self.DEBUG       = False   # a "VERBOSE" setting
        if (self.DEBUG):
            self.SILENT  = False
        self.REPORT      = False
        self.DONE        = False
        self.SKIP_ERRORS = False
        self.status      = 'created'   # (OpenMI 2.0 conventions)

        self.in_directory     = None
        self.out_directory    = None
        self.site_prefix      = None
        self.case_prefix      = None
        self.cfg_prefix       = None       ###### (9/18/14)
        self.comp_status      = 'Enabled'
        self.OVERWRITE_OK     = False   # Now in path_info.cfg; 2022-02-16
              
        # NB! This probably won't work here, because a driver
        #     may be instantiated later that then changes the
        #     current working directory.
##        self.cfg_directory  = os.getcwd()
##        self.cfg_prefix     = 'Case5'
        
    #   __init__()
    #-------------------------------------------------------------------
    def get_bmi_version(self):
    
        return '2.1'
        
    #   get_bmi_version()
    #-------------------------------------------------------------------
    # BMI methods for fine-grained control of model
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver"):

        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file
        # Will use cfg_file from above.
        #-----------------------------------------------
        ## self.set_constants()
        self.initialize_config_vars()  # calls set_directories().
        self.read_grid_info()
        self.initialize_basin_vars()
        self.initialize_time_vars()

        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()

        #--------------------------------------------
        # Set any input variables that are computed
        #---------------------------------------------------
        # NOTE:  This must be called in initialize(), AFTER
        #        calling read_input_files().
        #---------------------------------------------------
        # self.set_computed_input_vars()
        
##        self.open_output_files()

        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self):

        #------------------------------------------------------
        # Note: This method must be overridden by each model.
        #------------------------------------------------------
        # update model vars    
        self.update_time()
        
    #   update()
    #-------------------------------------------------------------------
    def update_frac(self, time_frac):
        """Update model by a fraction of a time step."""
        
        time_step = self.get_time_step()
        self.time_step = time_frac * time_step
        self.update()
        self.time_step = time_step

    #   update_frac()
    #-------------------------------------------------------------------
    def update_until(self, later_time):
        """Update model until a particular time."""

        n_steps = (later_time - self.get_current_time()) / self.get_time_step()
        for k in range(int(n_steps)):
            self.update()
        self.update_frac(n_steps - int(n_steps))
        
    #   update_until()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalizing'  # (OpenMI 2.0 convention)
        self.close_input_files()    ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'  # (OpenMI 2.0 convention)
        
    #   finalize()
    #-------------------------------------------------------------------
    def get_component_name(self):

        #----------------------------------------------------------
        # Note: This method must be overridden by each component.
        #----------------------------------------------------------    
        return 'BMI_base'

    #   get_component_name()
    #-------------------------------------------------------------------
    def get_input_item_count(self):
    
        input_names = self.get_input_var_names()
        return len(input_names)
        
    #   get_input_item_count()
    #-------------------------------------------------------------------
    def get_output_item_count(self):
    
        output_names = self.get_output_var_names()
        return len(output_names)
        
    #   get_output_item_count()
    #-------------------------------------------------------------------
    def get_input_var_names(self):

        #-------------------------------------------------------------
        # Note: This method must be overridden by each model.
        #-------------------------------------------------------------
        # Note: There may be a way to retrieve this list auto-
        #       matically using code similar to self.has_variable().
        #-------------------------------------------------------------
        # Note:  NumPy dtype of returned array will be "|Sn",
        #        where n is the length of the longest string.
        #-------------------------------------------------------------
        items = ['None']
        return np.array( items )   # (string array vs. list)
    
    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):

        #------------------------------------------------------
        # Note: This method must be overridden by each model.
        #------------------------------------------------------
        items = ['None']
        return np.array( items )   # (string array vs. list)
    
    #   get_output_var_names()    
    #-------------------------------------------------------------------
    def get_var_grid(self, long_var_name):    
    
        pass   ##################################

    #   get_var_grid()
    #-------------------------------------------------------------------
    def get_var_type(self, long_var_name):

        #--------------------------------------------------------
        # Notes: I think we should use the same type names that
        #        NumPy "ndarrays" store as "dtype" because they
        #        are unambiguous.  We assume here that all vars
        #        are NumPy types with dtype defined. e.g.
        #        uint8, int16, int32, float32, float64, etc.
        #--------------------------------------------------------

        # CSDMS BMI Python Example method        
        # return str(self.get_value_ptr(long_var_name).dtype)
        
        var_name = self.get_var_name( long_var_name )  # (2/20/12)

        dtype = 'unknown'  # (For exec in Python 3.x.)
        try:
            dtype = eval( "self." + var_name + ".dtype" )
        except:
            dtype = 'unknown'
        return str(dtype)       # (need str() here)

        #----------------------------------------------------------------
##        # This should also work.
##        exec( "HAS_DTYPE = hasattr( self." + var_name + ", 'dtype')", {}, locals() )
##        if (HAS_DTYPE):
##            exec( "dtype = self." + var_name + ".dtype", {}, locals() )
##        else:
##            dtype = 'unknown'   ###############
##        return dtype
    
    #   get_var_type()    
    #-------------------------------------------------------------------
    def get_var_units(self, long_var_name):

        #------------------------------------------------------
        # Note: This method must be overridden by each model.
        #------------------------------------------------------
        # Should define the "map" just once in __init__().
        #------------------------------------------------------
        # EXAMPLE:
        self._var_units_map = {
            'Q':     'm3/s',
            'Qs':    'm3/s', ########
            'A':     'km2',
            'S':     'none',  ## (or 'm/m')
            'z':     'm',
            'dz_dt': 'm/yr',
            'DEM':   'm',
            #--------------------------
            'T_air':     'degrees_C',
            'T_surf':    'degrees_C',
            'rho_H2O':   'kg/m3',
            'rho_snow':  'kg/m3',
            'Q_sum':     'W/m2',
            'h_swe':     'm',
            'h_snow':    'm' }
        var_name = self.get_var_name( long_var_name )
        return self._var_units_map[ var_name ]

        #---------------------------------------------------
        # Most components map long name directly to units.
        #---------------------------------------------------
        ## return self._var_units_map[ long_var_name ]

    #   get_var_units()
    #-------------------------------------------------------------------    
#     def get_var_rank(self, long_var_name):
# 
#         #----------------------------------------------------
#         # In BMI 2.0, this was replaced by get_grid_rank().
#         #---------------------------------------------------- 
#         var_name = self.get_var_name( long_var_name )
#         
#         try:
#             rank = eval( "np.ndim(self." + var_name + ")" )
#         except:
#             rank = -1
#         return rank
#     
#         ## return np.int32( rank )
# 
#     #   get_var_rank()
    #-------------------------------------------------------------------
    def get_var_itemsize(self, long_var_name):

        # CSDMS BMI Python Example method
        ## return np.dtype(self.get_var_type(name)).itemsize
        
        var_name = self.get_var_name( long_var_name )

        try:
            itemsize = eval( "self." + var_name + ".itemsize" )
        except:
            itemsize = -1
        return itemsize

    #   get_var_itemsize()
    #-------------------------------------------------------------------
    def get_var_nbytes(self, long_var_name):

        # CSDMS BMI Python Example method
        # return self.get_value_ptr(long_var_name).nbytes
        
        var_name = self.get_var_name( long_var_name )

        try:
            nbytes = eval( "self." + var_name + ".nbytes" )
        except:
            nbytes = -1
        return nbytes

    #   get_var_nbytes()    
    #-------------------------------------------------------------------
    def get_var_location(self):    

        #-----------------------------------------------------
        # Note: This method should be overridden, as needed.
        #-----------------------------------------------------   
        return 'node'
        
    #   get_var_location()
    #-------------------------------------------------------------------
    # BMI methods to get time-related info
    #-------------------------------------------------------------------
    def get_current_time(self):

        return self.time
    
    #   get_current_time()
    #-------------------------------------------------------------------
    def get_start_time(self):

        #--------------------------------------------------------
        # Notes: Most TF and Erode components do not specify
        #        a particular start time.  An exception is the
        #        TF Meteorology component, which we could check
        #        for here and process differently.
        #--------------------------------------------------------
        return np.float64( 0  )
    
    #   get_start_time()
    #-------------------------------------------------------------------
    def get_end_time(self):

        #--------------------------------------------------
        # Does model have a "stop_time" attribute ?
        # Even if it does, it may not be honored exactly.
        #--------------------------------------------------
        if (hasattr( self, 'stop_time' )):
            return np.float64( self.stop_time )

        #---------------------------------
        # Can we compute the stop_time ?
        #---------------------------------
        dt_type    = self.get_attribute( 'time_step_type' )
        COMPUTABLE = (dt_type == 'fixed') and (hasattr(self, 'n_steps'))
        
        if (hasattr( self, 'stop_code' )):
            if (self.stop_code != 0):
                COMPUTABLE = False  # (overrides COMPUTABLE above)

        if (COMPUTABLE):
            return np.float64( self.n_steps * self.dt )
        else:  
            print('##############################################')
            print(' ERROR: Unable to compute model stop_time.')
            print('##############################################')
            print(' ')
            return np.finfo("d").max
            ### return np.float64( -1 )           
            
    #   get_end_time()     
    #-------------------------------------------------------------------        
    def get_time_units(self):

        #-----------------------------------------------------------
        # Notes: The get_attribute() method should always provide
        #        a model's time units.  However, some components
        #        have an initialize_time_vars() method that saves
        #        a units string in self.time_units.  These strings
        #        are lower-case, e.g. 'seconds', 'minutes',
        #        or 'years'.  They should conform to UDUNITS.
        #-----------------------------------------------------------
        unit_str = self.get_attribute( 'time_units' )
        return unit_str
        ## return self.time_units  
    
    #   get_time_units()
    #------------------------------------------------------------------- 
    def get_time_step(self):

        return np.float64( self.dt )
    
    #   get_time_step()   
    #-------------------------------------------------------------------
    # BMI variable getter and setter methods
    #-------------------------------------------------------------------
    # Note:  Before v2.0, these functions contained the substring
    #        "values" vs. "value" and there was no "dest" argument.
    #-------------------------------------------------------------------    
    # def get_values( self ):
    # def get_values_ptr( self ):
    # def get_values_at_indices( self ):
    # def set_values( self ):
    # def set_values_at_indices( self ):
    #-------------------------------------------------------------------                     
    def get_value(self, long_var_name, dest):

        #------------------------------------------------------- 
        # Note: The value returned by getattr() will have rank
        #       and data type that goes with long_var_name.
        #------------------------------------------------------- 
        var_name = self.get_var_name( long_var_name )

        # CSDMS BMI Python Example method
        ## dest[:] = self.get_value_ptr(long_var_name).flatten()
            
        try:
            dest[:] = getattr(self, var_name)
        except:
            print('ERROR in BMI_base.get_values()')
            print('    for var_name =', var_name)
            print('    Returning 0.')

            #-----------------------------------------           
            # Return array of zeros with right shape
            #-----------------------------------------
            try:
                shape = eval( "np.shape(self." + var_name + ")" )
            except:
                shape = ()    # numpy scalar
            dtype = self.get_var_type( long_var_name )
            value = np.zeros( shape, dtype=dtype)
            dest[:] = value
        return dest
        
    #   get_value()
    #-------------------------------------------------------------------                     
    def get_value_ptr(self, long_var_name):

        var_name = self.get_var_name( long_var_name )    
        return getattr(self, var_name)

    #   get_value_ptr()
    #-------------------------------------------------------------------
    def get_value_at_indices(self, long_var_name, dest, IDs):

        #---------------------------------------------------------
        # Notes: This function was tested in the new Diversions
        #        component on (2/18/10).
        #---------------------------------------------------------
        # CSDMS BMI Python Example method
        dest[:] = self.get_value_ptr(long_var_name).take(IDs)
        return dest

        #--------------------------------------------------------        
#         var_name = self.get_var_name( long_var_name )
#         var_IDs_name = var_name + '.flat[IDs]'
#         result = getattr(self, var_IDs_name)
#         dest[:] = np.array(result, dtype='float64')
#         return dest
        #--------------------------------------------------------
#         var_name = self.get_var_name( long_var_name )                    
#         try:
#             var_IDs_name = var_name + '.flat[IDs]'
#             result = getattr(self, var_IDs_name)  ## (2/19/13)
#             dest[:] = np.array(result, dtype='float64')
#         except:
#             print('ERROR in BMI_base.get_value_at_indices().')
#             print('    Returning zeros.')
#             dtype = self.get_var_type( long_var_name )
#             dest[:] = np.zeros(len(IDs), dtype=dtype)
#         return dest
        
    #   get_values_at_indices()
    #-------------------------------------------------------------------
    def set_value(self, long_var_name, value):

        #---------------------------------------------------------------
        # Notes: The "var_name" string cannot contain a ".". (5/17/12)
        #---------------------------------------------------------------
        # (2/7/13) We are now using 0D numpy arrays as a way to
        # produce "mutable scalars" that allow a component with a
        # reference to the scalar to see changes to its value.
        # But we can't apply np.float64() to the value as we did
        # before or it destroys the reference.
        # See BMI_base.initialize_scalar() for more information.
        #---------------------------------------------------------------
        
        # CSDMS BMI Python Example method       
        # val = self.get_value_ptr(long_var_name)
        # val[:] = value.reshape(val.shape)
        #----------------------------------------------- 
        var_name = self.get_var_name( long_var_name )
        setattr( self, var_name, value )
 
    #   set_value()
    #-------------------------------------------------------------------
    def set_value_at_indices(self, long_var_name, IDs, value):

        #--------------------------------------------------------
        # Notes: This function was tested in the new Diversions
        #        component on (2/18/10).
        #--------------------------------------------------------

        # CSDMS BMI Python Example method        
        # val = self.get_value_ptr(long_var_name)
        # val.flat[IDs] = value
        
        var_name = self.get_var_name( long_var_name )
        var_IDs_name = var_name + '.flat[IDs]'   ##############
        setattr( self, var_IDs_name, value )
        
    #   set_values_at_indices()   
    #-------------------------------------------------------------------
    # BMI grid information methods
    #-------------------------------------------------------------------
    def get_grid_rank( self, grid_id ):
    
        if (grid_id == 0):
            return 2    # 2D array, shape of DEM
        else:
            return 0    # 0D array or scalar
        
    #   get_grid_rank()
     #-------------------------------------------------------------------
    def get_grid_size( self, grid_id ):

        if (grid_id == 0):
            if not(hasattr( self, 'grid_info' )):
                self.read_grid_info()
            info = self.grid_info
            return (info.ny * info.nx)        
        else:
            return 1
    
    #   get_grid_size()
    #-------------------------------------------------------------------
    def get_grid_type( self, grid_id ):

        if (grid_id == 0):
            if not(hasattr( self, 'grid_info' )):
                self.read_grid_info()
            info = self.grid_info
            if (info.pixel_geom == 1):        
                return 'uniform_rectilinear'
            else:
                ########################################################
                # Need to then implement more grid info functions !!!
                ########################################################
                return 'structured_quadrilateral'
        else:
            return 'scalar'

    #   get_grid_type()
    #-------------------------------------------------------------------
    def get_grid_shape(self, grid_id, shape):

        if (grid_id == 0):
            if not(hasattr( self, 'grid_info' )):
                self.read_grid_info()
            info = self.grid_info
            #---------------------------------------
            # Note: Ordering must be [nz, ny, nx].
            #---------------------------------------
            shape[:] = np.array( [info.ny, info.nx] )
        else:   # scalar
            pass  ############################
            ## shape[:] = ()   # must be array
        
        return shape
    
    #   get_grid_shape()
    #-------------------------------------------------------------------
    def get_grid_spacing(self, grid_id, spacing):

        #-------------------------------------------------------
        # Note: Ordering must be [dz, dy, dx].
        #-------------------------------------------------------
        # Note: xres and yres could have angular units like
        #       arcseconds, if (info.pixel_geom == 0).  In
        #       that case, xres (or dx) will vary with latitude
        #       and this method isn't really appropriate.
        #--------------------------------------------------------
        if (grid_id == 0):
            if not(hasattr( self, 'grid_info' )):
                self.read_grid_info()
            info = self.grid_info
            spacing[:] = np.array( [info.yres, info.xres] )
        else:  # scalar
            pass  ############################
            ## spacing[:] = 0
        return spacing

    #   get_grid_spacing()
    #-------------------------------------------------------------------
    def get_grid_origin(self, grid_id, origin):

        if (grid_id == 0): 
            if not(hasattr( self, 'grid_info' )):
                self.read_grid_info()
            info = self.grid_info
            #---------------------------------------
            # Note: Ordering must be [z0, y0, x0].
            #---------------------------------------
            origin[:] = np.array( [info.y_south_edge, info.x_west_edge] )
        else:  # scalar
            pass  ############################
            ## origin[:] = 0
             
        return origin

    #   get_grid_origin()    
    #-------------------------------------------------------------------
    def get_grid_x(self):
    
        raise NotImplementedError("get_grid_x") 

    #   get_grid_x()
    #------------------------------------------------------------------- 
    def get_grid_y(self):
    
        raise NotImplementedError("get_grid_y") 

    #   get_grid_y()
    #------------------------------------------------------------------- 
    def get_grid_z(self):
    
        raise NotImplementedError("get_grid_z")
        
    #   get_grid_z()
    #-------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------  
    # Non-BMI (private or convenience) Functions  
    #-------------------------------------------------------------------
    def get_var_name(self, long_var_name):

        """Map long_var_name to short, internal var name."""
        #------------------------------------------------------
        # Note: This method must be overridden by each model.
        #------------------------------------------------------
        # Should define the "map" just once in __init__().
        #------------------------------------------------------
        # EXAMPLE:
        self._var_name_map = {
            'model__time_step':'dt' }

        return self._var_name_map[ long_var_name ]
        
    #   get_var_name()
    #-------------------------------------------------------------------        
    def get_status(self):

        #-----------------------------------------------------
        # Notes: Return component status as a string.  The
        #        possible return values are from OpenMI 2.0:
        #
        #        created, initializing, initialized,
        #        updating, updated, finalizing, finalized,
        #        failed (could add "stopped").
        #-----------------------------------------------------
        return self.status

    #   get_status()
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        #-------------------------------------------------------------
        # Note: This dictionary must be provided for each component.
        #-------------------------------------------------------------
#         self._att_map = {
#             'model_name':         'Model_Name',
#             'version':            '1.0',
#             'author_name':        'Scott D. Peckham',
#             'grid_type':          'uniform',
#             'time_step_type':     'fixed',
#             'step_method':        'explicit',
#             #------------------------------------------------------
#             'comp_name':          'Model_Name',
#             'model_family':       'TopoFlow',
#             'cfg_template_file':  'Model_Name.cfg.in',
#             'cfg_extension':      '_Model_Name.cfg',
#             'cmt_var_prefix':     '/ModelName/Input/Var/',
#             'gui_xml_file':       '/home/csdms/cca/model_name/1.0/src/share/cmt/gui/Model_Name.xml',
#             'dialog_title':       'ProcessName: Component Parameters',
#             'time_units':         'seconds' }

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print('###################################################')
            print(' ERROR: Could not find attribute: ' + att_name)
            print('###################################################')
            print(' ')

    #   get_attribute()
    #-------------------------------------------------------------------
    def get_grid_attribute(self, long_var_name, att_name):

        amap = {'x_units':'meters', 'y_units':'meters',
                'z_units':'meters', 'ellipsoid':'None',
                'datum':'None', 'projection':'None'}
                ## 'utm_zone':'None'}

        try:
            return amap[ att_name.lower() ]
        except:
            print('######################################################')
            print(' ERROR: Could not find grid attribute: ' + att_name)
            print('######################################################')
            print(' ')            

    #   get_grid_attribute()
    #-------------------------------------------------------------------
    def run_model(self, cfg_directory=None, cfg_prefix=None,
                  n_steps=5):

        #--------------------------------------------------
        # NOTE: This method is not called from the go()
        #       method in the IMPL file, but it is often
        #       called by a unit_test().
        #--------------------------------------------------

        #--------------------------------------------------
        # All components, including this one (the driver)
        # will look in the CWD for their CFG file.
        #--------------------------------------------------
        if (cfg_directory is not None):
            os.chdir( cfg_directory )
        self.cfg_prefix = cfg_prefix

        #-----------------------------------------
        # Initialize the model run (driver mode)
        #---------------------------------------------
        # This will set in_directory, out_directory,
        # site_prefix and case_prefix
        #---------------------------------------------
        self.initialize( cfg_prefix=cfg_prefix, mode='driver' )
            
        #----------------------------------------------------------- 
        # Note:  If the number of timesteps is specified in a
        #        component's CFG file and is then saved by
        #        "read_cfg_file()" as "n_steps" then we should
        #        honor that setting.  Otherwise we use the n_steps
        #        argument.
        #-----------------------------------------------------------
        if not(hasattr( self, 'stop_code' )):
           self.stop_code = 0
        #--------------------------------------
        if (hasattr(self, 'n_steps')):
            ## print 'NUMBER OF STEPS =', self.n_steps  ####
            n_steps = self.n_steps
        # else:
            # (Use the n_steps argument.)

        #-------------------------------------------        
        # Note: __init__() sets self.DONE to False
        #-------------------------------------------
        while not(self.DONE):
            if not(self.SKIP_ERRORS):
                #-------------------------------------------
                # Exceptions will not be caught and
                # Python error messages will be displayed.
                #-------------------------------------------
                if (self.DEBUG):
                    print('time_index =', self.time_index)
                self.update()
                # self.update( -1 )  # (BMI later: use own dt)
            else:   
                try:
                    self.update()
                except:
                    print('ERROR in run_model() method at:')
                    print('   time_index =', self.time_index)
                    self.status = 'failed'
                    self.DONE = True

            #------------------------------------------------
            # If the model has set self.DONE = True, then
            # stop, even if we haven't reached n_steps yet.
            #------------------------------------------------
            self.check_finished()   # (see below: 2/13/12)

        #-------------------------
        # Finalize the model run
        #-------------------------
        self.finalize()

    #   run_model()
    #-------------------------------------------------------------------
    def initialize_time_vars(self, units='seconds'):

        #--------------------------------------------------------
        # NOTE:  Since time is not an input or output variable
        #        passed between components, it may not be a
        #        a good idea to use "initialize_scalar()" here.
        #        The data type is then "ndarray" (0D).
        #        See notes for:  initialize_scalar()
        #        This may be the cause of a strange problem
        #        when "time" is written to a netCDF file.
        #        Now, np.float64() is applied before writing.
        #        CHECK if this affects time interpolation.
        #--------------------------------------------------------
        
        #------------------
        # Start the clock
        #------------------
        self.start_time = time.time()
        
        #--------------------------------
        # Initialize the time variables
        #--------------------------------
        self.time_units = units.lower()
        self.time_index = np.int32(0)
        self.time       = self.initialize_scalar(0, dtype='float64') 
        self.DONE       = False

        #--------------------------------------------------------
        # 2022-04-08.  Result is: <class> 'ndarray' (0D array)
        # See notes for: BMI.initialize_scalar().
        #--------------------------------------------------------
        ## print('In BMI.initialize_time_vars:')
        ## print('    type(self.time) =', type(self.time))

        #-------------------------------------------
        # (2/5/12) Set default stopping method:
        # (0=n_steps, 1=stop_time, 2=steady-state)
        #-------------------------------------------
        # Don't do this in set_constants().
        #-------------------------------------------
        if not(hasattr(self, 'stop_code')):
            self.stop_code = 0

        #---------------------------------------
        # Make sure n_steps is defined for new
        # use in update_time().  Was not set
        # for diversions_base. (11/14/11)
        #---------------------------------------
        if not(hasattr(self, 'n_steps')):
            self.n_steps = 1
 
        #--------------------------
        # Time conversion factors
        #--------------------------
        self.sec_per_year = np.float64(365) * 24 * 3600
        self.min_per_year = np.float64(365) * 24 * 60
        
        #-------------------------------------------
        # For backward compatibility with TopoFlow
        #------------------------------------------------
        # Use 1D array (mutable) vs. scalar (immutable)
        # to permit passing as a reference.  (2/4/13)
        #------------------------------------------------
        self.time_sec = self.initialize_scalar(0, dtype='float64')
        self.time_min = self.initialize_scalar(0, dtype='float64')
            
        #--------------------------------------------
        # For print_time_and_value() function below
        #--------------------------------------------
        # Substract 100 seconds so we'll always
        # print values at time zero. (6/29/10)
        #--------------------------------------------
        self.last_print_time = time.time() - 100.0
        
##        self.last_check_time  = time.time()  # (for user interrupt)
##        self.last_plot_time   = np.float32(0.0)   ### CHECK ###
        
    #   initialize_time_vars()
    #-------------------------------------------------------------------
    def update_time(self, dt=-1):

        #--------------------------------------------------
        # Note: The BMI update() method has a dt argument
        #       that is passed to this routine.
        #--------------------------------------------------
        
        #---------------------
        # Increment the time
        #---------------------
        self.time_index += 1
        if (dt == -1):
            self.time += self.dt  # (use same units as dt)
        else:
            self.time += dt

        #------------------------------------------------
        # Compute and store the current time in seconds
        # (time_sec) and in minutes (time_min).
        #------------------------------------------------
        # Added 'weeks' and 'months': 2022-02-17
        #------------------------------------------------
        if (self.time_units == 'seconds'):
            self.time_sec = self.time         # [seconds]
        elif (self.time_units == 'minutes'):
            self.time_sec = self.time * np.float64(60)
            ## self.time_min = self.time
        elif (self.time_units == 'hours'):
            self.time_sec = self.time * np.float64(3600)
        elif (self.time_units == 'days'):
            self.time_sec = self.time * np.float64(3600) * 24
        elif (self.time_units == 'weeks'):
            self.time_sec = self.time * np.float64(3600) * 24 * 7
        elif (self.time_units == 'months'):
            #---------------------------------------
            # NOTE!  This is just an approximation
            #---------------------------------------
            self.time_sec = self.time * self.sec_per_year / 12
        elif (self.time_units == 'years'):
            #-----------------------------------
            # Used by GC2D and Erode (12/4/09)
            #-----------------------------------
            self.time_sec = self.time * self.sec_per_year  ####
        else:
            print('### WARNING: time_units =', self.time_units, 'not recognized.')
            self.time_sec = self.time

        #------------------------------------------
        # Compute & store current time in minutes
        #------------------------------------------
        self.time_min = self.time_sec / np.float64(60)   # [minutes]

        #---------------------------------------------------
        # Are we DONE yet?  BMI_base.py contains a basic
        # version of check_finished(), but components that
        # inherit from BMI_base.py may override it.
        #---------------------------------------------------
        # Moved here for new framework approach. (5/18/12)
        #---------------------------------------------------
        self.check_finished()
                
    #   update_time()
    #-------------------------------------------------------------------
    def check_finished(self):

        #---------------------------------------------------------
        # Note: If self.DONE has already been set to True by
        #       another function or component, this function
        #       preserves that setting (see below).
        #---------------------------------------------------------
        #       TINY_DZ can occur either because dt required for
        #       stability is really small or because we have
        #       converged to a steady-state landscape.
        #---------------------------------------------------------
        #       Moved here from erode_base.py on 2/13/12.
        #---------------------------------------------------------

        #---------------------------------------------------
        # If component isn't the driver, return. (2/20/12)
        #---------------------------------------------------
        ### if (self.mode == 'nondriver'):
        if (self.mode != 'driver'):
            return
        
        if (self.stop_code == 0):
            #---------------------
            # Stop after n_steps
            #---------------------
            TIMES_UP = (self.time_index >= self.n_steps)
        elif (self.stop_code == 1):
            #-----------------------
            # Stop after stop_time 
            #-----------------------
            TIMES_UP = (self.time >= self.stop_time)
        elif (self.stop_code == 2):
            #-----------------------------------------
            # Stop if "steady-state", but only works
            # as written here for global timesteps. 
            #-----------------------------------------
            TIMES_UP = (self.dz_max < self.dz_tolerance)

        self.DONE = (self.DONE or TIMES_UP)

    #   check_finished()
    #-------------------------------------------------------------------
    def print_time_and_value(self, value, var_name='Q_out',
                             units_name='[m^3/s]',
                             interval=2.0,
                             PRINT_INDEX=False):

        #-----------------------------------------
        # (8/2/10) Print message about interval.
        #-----------------------------------------
        if (self.time_index == 0):
            istr = str(interval)
            print('Will print values every ' + istr + ' seconds.')
            
        #---------------------------------------------------
        # Note: Print the model time, in minutes, and the
        #       current value of "var", at the specified
        #       real-time "interval" (in seconds).
        #---------------------------------------------------
        # Note: Plotting hydrograph at same interval is
        #       generally too infrequent.
        #---------------------------------------------------
        elapsed_time = (time.time() - self.last_print_time)
        if (elapsed_time > interval):
            if (self.time_units == 'seconds'):
                cur_time = self.time_min
                time_units_str = ' [min]'
            else:
                cur_time = self.time
                time_units_str = ' [' + self.time_units + ']' 
            #------------------------
            # Build the time string
            #------------------------
            time_str = 'Time = ' + ("%10.2f" % cur_time)
            time_str = time_str + time_units_str

            #------------------------
            # Build the value string
            #------------------------
            if (value < 1e6):
                val_str   = ("%10.5f" % value)
                units_str = units_name
            else:
                val_str = ("%10.5f" % (value / 1e6))
                units_str = 'x 10^6 ' + units_name          
            #--------------------------------------
            var_str  = var_name + ' = ' + val_str
            var_str  = var_str  + ' ' + units_str                 
            # var_str  = var_name + ' = ' + ("%10.5f" % var)
            # var_str  = var_str  + ' ' + units_name          
            #----------------------------------------     
            print((time_str + ',  ' + var_str))
            #-----------------------------------------------------
            if (PRINT_INDEX):
                index = (self.time_index + 1)  # (starts at 0)
                print('n =', index, 'of', self.n_steps)
            #-----------------------------------------------------                
            self.last_print_time = time.time()

    #   print_time_and_value()
    #-------------------------------------------------------------------
    def get_run_time_string(self, proc_name='component',
                            sec_digits=4, seconds=None,
                            SILENT=None):
                            ### SUB_PROCESS=False) 

        #------------------------------------------------------
        # If "seconds" argument is only provided for testing.
        # You can provide this value to make sure that the
        # minuts, hours, days, etc. are computed correctly.
        #------------------------------------------------------
        if (seconds == None):    
            finish  = time.time()
            seconds = (finish - self.start_time)

        #-------------------------------------------------
        # Could later add this to help gauge performance
        # for all models.  It is currently coded into
        # erode_d8_local.finalize(). But not all models
        # have self.time in "years". (10/4/11)
        #-------------------------------------------------
##        finish         = time.time()
##        run_time_secs  = (finish - self.start_time)
##        run_time_hours = (run_time_secs / 3600.0)
##        sim_time_years = self.time
##        years_per_hour = (sim_time_years / run_time_hours)
##        print ' '
##        print 'Years per hour =', years_per_hour
        
        #----------------------------------
        # Compute minutes, hours and days
        #----------------------------------
        dec_part  = (seconds % np.float32(1.0))     #(Save decimal part)
        days      = np.int32(seconds) / np.int32(86400)
        secs_left = np.int32(seconds) % np.int32(86400)
        hours     = (secs_left / np.int32(3600))
        secs_left = (secs_left % np.int32(3600))
        minutes   = (secs_left / np.int32(60))
        seconds   = (secs_left % np.int32(60))
        #-----------------------------------------
        #hours     = long(seconds)  /  3600L
        #secs_left = long(seconds) mod 3600L
        #minutes   = (secs_left  /  60L)
        #seconds   = (secs_left mod 60L)
        
        #----------------------------
        # Construct the time string
        #----------------------------
        time_string = ''
        #--------------------------------------------------------
        if (days > 0):    
            if (days > 1):    
                e0 = ' days, '
            else:    
                e0 = ' day, '
            time_string += str(days) + e0
        #--------------------------------------------------------
        if (hours > 0):    
            if (hours > 1):    
                e1 = ' hours, '
            else:    
                e1 = ' hour, '
            time_string += str(hours) + e1
        #--------------------------------------------------------
        if (minutes > 0):    
            if (minutes > 1):    
                e2 = ' minutes, '
            else:    
                e2 = ' minute, '
            time_string += str(minutes) + e2
        
        #-----------------------------------------
        # Default is 4 digits after the decimal.
        #-----------------------------------------
        dec_pastr = ('.' + str(dec_part)[2:2+sec_digits])
        time_string += str(seconds) + dec_pastr + ' seconds.'

        return time_string
            
##            if (SUB_PROCESS):    
##                PART1 = '>> '
##            else:    
##                PART1 = ''
##            print (PART1 + 'Run time for ' + procname + ' = ')
##            print (PART1 + time_string)
##            print ' '
     
    #   get_run_time_string()
    #-------------------------------------------------------------------
    def print_run_time(self, proc_name='component',
                       sec_digits=4, seconds=None,
                       SILENT=None):
                       ### SUB_PROCESS=False) 


        print(('Run time for ' + proc_name + ' = '))
        print(self.get_run_time_string())
        print(' ')
     
    #   print_run_time()
    #-------------------------------------------------------------------    
    def print_final_report(self, comp_name='BMI component',
                           mode='nondriver'):

        if (mode == 'nondriver'):
            print(comp_name + ': Finished.')
            return

        if not(hasattr( self, 'in_directory' )):
            return   # (2/20/12)
        
        #-------------------
        # Print the report
        #-------------------
        hline = ''.ljust(60, '-')
        print(hline)
        print(comp_name)
        print(time.asctime())
        print(' ')
        print('Input directory:      ' + self.in_directory)
        print('Output directory:     ' + self.out_directory)
        print('Site prefix:          ' + self.site_prefix)
        print('Case prefix:          ' + self.case_prefix)
        print(' ')

        #-----------------------------------
        # Construct sinulation time string
        #-----------------------------------
        sim_units    = ' [' + self.time_units + ']'
        sim_time_str = str(self.time) + sim_units
        
        #----------------------------
        # Construct run time string
        #----------------------------
        run_time_str = self.get_run_time_string() 

        print('Simulated time:      ' + sim_time_str)
        print('Program run time:    ' + run_time_str)
        print(' ')
        print('Number of timesteps: ' + str(self.time_index))
        print('Process timestep:    ' + str(self.dt) + sim_units)
        print('Number of columns:   ' + str(self.nx))
        print('Number of rows:      ' + str(self.ny))
        print(' ')
        print('Finished. (' + self.case_prefix + ')')
        print(' ')

        ## finish_str = ': Finished. (' + self.case_prefix + ')'
        ## print finish_str
        ## print comp_name + finish_str
        
    #   print_final_report()
    #-------------------------------------------------------------------    
    def print_traceback(self, caller_name='TopoFlow'):

        print('################################################')
        print(' ERROR encountered in ' + caller_name)
        print('       Please check your input parameters.')
        print('################################################')
        
        traceback.print_exc()
        
    #   print_traceback()
    #-------------------------------------------------------------------
    def read_grid_info(self):

        #-------------------------------------------------------
        # Note: This isn't a BMI method, but is currently used
        #       by all of the TopoFlow and Erode components.
        #       Notice that it uses "rti_files" and "pixels".
        #-------------------------------------------------------
        # Read grid info from an RTI file that is
        # in the current working directory.
        #------------------------------------------
        if (self.DEBUG):
            print('Process component: Reading grid info...')
        
        rti_file = (self.site_prefix + '.rti')
        self.grid_info_file = self.topo_directory + rti_file
        info = rti_files.read_info( self.grid_info_file, SILENT=True )
        if (info == None):
            #-------------------------------------------------------
            # The RTI file is usually in "in_directory" and
            # uses "site_prefix", but when it needs to be
            # created, as for Erode (see erode_base.py), then
            # it uses "out_directory" and "case_prefix". (2/17/13)
            #-------------------------------------------------------
            print('ERROR: RTI file: ' + rti_file + ' not found in input directory:')
            ## print('In the input_directory:')
            print( self.topo_directory )
            #-------------------------------------------------------
            print('RTI file names usually begin with site_prefix')
            print('   and are stored in the input (or topo) directory.')
            print('In rare cases, they begin with the case_prefix')
            print('   and are stored in the output directory.')
            print('Looking in the output directory...')
            print()
            rti_file = (self.case_prefix + '.rti') 
            self.grid_info_file = self.out_directory + rti_file
            info = rti_files.read_info( self.grid_info_file, SILENT=True )
            if (info == None):
                print('ERROR: RTI file: ' + rti_file + ' not found in output directory.')
                ## print('In the output_directory:')
                print( self.out_directory )
                sys.exit()   ####################
                
        # print '##### In BMI_base.read_grid_info():'
        # print '##### in_directory   =', self.in_directory
        # print '##### grid_info_file =', self.grid_info_file
        
        #----------------------
        # Convenient synonyms
        #-----------------------------------------------------
        # Note that "info" has additional derived attributes
        # such as: n_pixels, bpe, grid_size and SWAP_ENDIAN.
        #-----------------------------------------------------
        self.rti = info
        self.nx  = info.ncols
        self.ny  = info.nrows

        #------------------------------------------------
        # (2/17/12) Start to phase out "rti" in favor
        #           of "grid_info".  See methods below.
        #------------------------------------------------
        self.grid_info    = self.rti
        self.grid_info.nx = self.rti.ncols
        self.grid_info.ny = self.rti.nrows
        
        #------------------------------------------------
        # Get grid cell areas, "da", which is either a
        # scalar (if same for all grid cells) or a grid
        # with default units of "m^2".
        #------------------------------------------------
        self.da = pixels.get_da( info )

        #-----------------------------------------        
        # Save the edge IDs to self (2023-09-11)
        #-----------------------------------------
        self.save_edge_ids()
        
    #   read_grid_info()
    #-------------------------------------------------------------------
    def save_edge_ids(self):

        #-----------------------------------------------
        # Note:  Copied from d8_base.py to BMI_base.py
        #        on 2023-09-11, for mass balance use,
        #        available to every BMI component.
        #-----------------------------------------------
        
        #------------------------------------------
        # Get IDs of edge pixels, making sure not
        # to double-count the corner pixels
        #------------------------------------------
        nx    = self.nx
        ny    = self.ny
        T_IDs = np.arange(nx, dtype='int32')
        B_IDs = T_IDs + (ny - 1) * nx
        L_IDs = (1 + np.arange(ny - 2, dtype='int32')) * nx
        R_IDs = L_IDs + (nx - 1)
        edge_IDs = np.concatenate([T_IDs, B_IDs, L_IDs, R_IDs])

        #-------------------------------------------
        # Save IDs as a tuple of row indices and
        # calendar indices, "np.where" style
        #-------------------------------------------        
        self.edge_IDs = divmod(edge_IDs, nx)   ##  NB! (row, col)

        #------------------------------------------
        # Save IDs as a 1D array of long-integer,
        # calendar-style indices
        #------------------------------------------
        # self.edge_IDs = edge_IDs
        
    #   save_edge_ids()
    #-------------------------------------------------------------------
    def read_path_info(self):

        #------------------------------------------------------------
        # Note:  In earlier versions of TopoFlow, every CFG file
        #        had: in_directory, out_directory, site_prefix and
        #        case_prefix at the top.  However, these are
        #        typically the same for all CFG files, so to
        #        prevent editing this block in every file, this
        #        method reads those 4 variables from a single file
        #        called: [case_prefix]_path_info.cfg.
        #        Since read_config_info() uses these 4 variables,
        #        this method must be called before that from the
        #        initialize_config_vars() method.
        #------------------------------------------------------------
        if (self.DEBUG):
            print('Reading path info config file.')
        # print( 'In read_path_info(), cfg_file =')
        # print( self.cfg_file )         

        #----------------------------------------------
        # Construct name of path info CFG file.
        # All we have to work with is cfg_file, which
        # was passed to initialize.  cfg_prefix and
        # cfg_directory have not been set yet.
        # cfg_directory is set in set_directories().
        # cfg_extension is the extension for the
        #   current component's CFG file.
        #-----------------------------------------------
        cfg_extension = self.get_attribute( 'cfg_extension' )
        cfg_prefix = self.cfg_file.replace( cfg_extension, '' )
        path_info_cfg = (cfg_prefix + '_path_info.cfg')
            
        #------------------------
        # Does CFG file exist ?
        #------------------------
        if (self.DEBUG):
            print('path_info_cfg =', path_info_cfg)
        if not(os.path.exists(path_info_cfg)):
            print('WARNING: path info CFG file not found:')
            print('         ' + path_info_cfg)
            return

        #---------------------------------------
        # Open file_info CFG file to read data
        #---------------------------------------
        cfg_unit = open( path_info_cfg, 'r' )

        #--------------------------------------------------
        # Recall that a "blank line", with just a (hidden)
        # newline character will not be null and will
        # have len(line) = 1.
        #--------------------------------------------------
        while (True):
            line  = cfg_unit.readline()
            if (line == ''):
                break                  # (reached end of file)
            
            COMMENT = (line[0] == '#')
            #--------------------------------------------
            # Using "|" as a delimiter means we can use
            # " ", "," or "=" in the filenames.
            #--------------------------------------------
            # Added ".strip()" on 10/25/11.
            #--------------------------------------------
            words   = line.split('|')  # (split on pipe)
            if (len(words) == 4) and not(COMMENT):
                var_name = words[0].strip()
                value    = words[1].strip()
                var_type = words[2].strip()  # (e.g. 'double', 'int', 'string', etc.)

                #---------------------------------------------
                # Note: Need to include locals() for "value"
                #---------------------------------------------
                exec( "self." + var_name + " = value", {}, locals() )

        #-------------------------------------
        # Check the directories and prefixes
        #-------------------------------------
        cfg_unit.close()
        # print '#### CALLING set_directories()...'
        self.set_directories()  # (Moved here: 2020-01-21)
        
    #   read_path_info()
    #-------------------------------------------------------------------
    def read_time_info(self):

        #--------------------------------------------------------------
        # Note:  Since we don't know which components a user will
        #        select, start_date, start_time, end_date & end_time
        #        needs to be available to every component.  This is
        #        needed to write time metadata to netCDF output files.
        #--------------------------------------------------------------
        #        This is virtually identical to read_path_info().
        #        Also see read_grid_info() in this file.
        #--------------------------------------------------------------        
        if (self.DEBUG):
            print('Reading time info config file.')

#         print ('In read_file_info, cfg_directory = ' + self.cfg_directory)
#         print ('In read_file_info, cfg_file   = ' + self.cfg_file)
#         print ('In read_file_info, cfg_prefix = ' + self.cfg_prefix)

        #---------------------------------------------------
        # Create an empty time_info object with dummy info
        # Will add things from "*_time_info.cfg" file.
        #---------------------------------------------------
#         class bunch:
#             def __init__(self, **kwds):
#                 self.__dict__.update(kwds)
#         time_info = bunch( dum = 0 )
        #-----------------------------------------
        class time_info_class:
            pass
        time_info = time_info_class()
#         time_info.dum = 0

        #---------------------------------------
        # Construct name of path info CFG file
        #--------------------------------------
        cfg_extension = self.get_attribute( 'cfg_extension' )
        cfg_prefix    = self.cfg_file.replace( cfg_extension, '' )
        cfg_file = (cfg_prefix + '_time_info.cfg')

        #------------------------
        # Does CFG file exist ?
        #------------------------
        if (self.DEBUG):
            print('cfg_file =', cfg_file)
        if not(os.path.exists(cfg_file)):
            print('WARNING: time info CFG file not found:')
            print('         Using arbitrary defaults.')
            print('         ' + cfg_file)
            #------------------
            # Set to defaults?
            #------------------
            time_info.start_date = '2020-01-01'
            time_info.start_time = '00:00:00'
            time_info.end_date   = '2020-01-02'
            time_info.end_time   = '00:00:00'
            time_info.duration   = 1440.0
            time_info.duration_units = 'minutes'
            self.time_info = time_info
            return

        #---------------------------------------
        # Open file_info CFG file to read data
        #---------------------------------------
        cfg_unit = open( cfg_file, 'r' )

        #--------------------------------------------------
        # Recall that a "blank line", with just a (hidden)
        # newline character will not be null and will
        # have len(line) = 1.
        #--------------------------------------------------
        while (True):
            line  = cfg_unit.readline()
            if (line == ''):
                break                  # (reached end of file)
            
            COMMENT = (line[0] == '#')
            #----------------------------------------------
            # Using "|" as a delimiter means we can use
            # " ", "," or "=" in the filenames.
            #----------------------------------------------
            # Added ".strip()" on 10/25/11.
            #----------------------------------------------
            # var_names will be:  start_date, start_time,
            #     end_date, end_time
            #----------------------------------------------            
            words   = line.split('|')  # (split on pipe)
            if (len(words) == 4) and not(COMMENT):
                var_name = words[0].strip()
                value    = words[1].strip()
                var_type = words[2].strip()  # (e.g. 'double', 'int', 'string', etc.)

                #---------------------------------------------
                # Note: Need to include locals() for "value"
                #---------------------------------------------
                ## exec( "self." + var_name + " = value", {}, locals() )
                exec( "time_info." + var_name + " = value", {}, locals() )
    
        #----------------------------------------------------    
        # Compute duration, save time_info object into self
        #----------------------------------------------------
        # Note that "time_res" will vary by component and
        # is not included here.
        #----------------------------------------------------        
        dur_units = 'minutes'  #######
        duration  = time_utils.get_duration( time_info.start_date, time_info.start_time,
                                             time_info.end_date, time_info.end_time,
                                             dur_units )
        time_info.duration = duration  # (2020-04-29) 
        time_info.duration_units = dur_units

        time_info.start_datetime = time_info.start_date + ' ' + time_info.start_time
        time_info.end_datetime   = time_info.end_date   + ' ' + time_info.end_time
                   
        self.time_info = time_info
        cfg_unit.close()  #######

    #   read_time_info()
    #-------------------------------------------------------------------
    def read_config_file(self):

        #----------------------------------------------------------
        # Notes: This version reads configuration settings from
        #        a new type of CFG file that only has var_name,
        #        value and type; more like key-value. (5/9/11)
        #        The GUI is no longer generated using a CFG file.
        #----------------------------------------------------------
        #        All of the TopoFlow and Erode components now use
        #        the same type of CFG file and all of them use
        #        this method to read settings from the file.
        #----------------------------------------------------------
        if not(self.SILENT):
            print('Reading config file into component state.')

        # print '######## In read_config_file(): cfg_file =', self.cfg_file
        # print '######## In read_config_file(): cfg_prefix =', self.cfg_prefix

        #############################################################
        # No longer needed; cfg_file is passed to BMI.initialize()
        # and saved as self.cfg_file. (9/17/14)
        #############################################################
        #---------------------------
        # Get name of the cfg_file
        #---------------------------
#         cfg_extension = self.get_attribute( 'cfg_extension' )  # (10/26/11)
#         cfg_directory = (os.getcwd() + os.sep)
#         file_name     = (self.cfg_prefix + cfg_extension)
#         self.cfg_file = (cfg_directory + file_name) 

        #------------------------
        # Does CFG file exist ?
        #------------------------
        if (self.DEBUG):
            print('cfg_file =', self.cfg_file)
        if not(os.path.exists(self.cfg_file)):
            print('WARNING: cfg_file not found:')
            print('         ' + self.cfg_file)
            return

        #-----------------------------
        # Open CFG file to read data
        #-----------------------------
        cfg_unit = open( self.cfg_file, 'r' )
        last_var_name = ''

        #-----------------------------------------
        # Save user input into component's state
        #--------------------------------------------------
        # Recall that a "blank line", with just a (hidden)
        # newline character will not be null and will
        # have len(line) = 1.
        #--------------------------------------------------
        while (True):
            line  = cfg_unit.readline()
            if (line == ''):
                break                  # (reached end of file)
            
            COMMENT = (line[0] == '#')
            #--------------------------------------------
            # Using "|" as a delimiter means we can use
            # " ", "," or "=" in the filenames.
            #--------------------------------------------
            # Added ".strip()" on 10/25/11.
            #--------------------------------------------
            words   = line.split('|')  # (split on pipe)
            if (len(words) == 4) and not(COMMENT):
                var_name = words[0].strip()
                value    = words[1].strip()
                var_type = words[2].strip()  # (e.g. 'double', 'int', 'string', etc.)
                ## help_str = words[3]
                READ_SCALAR   = False
                READ_FILENAME = False

                # For debugging
#                 print('var_name = ' + var_name)
#                 print('value    = ' + value)
#                 print('var_type = ' + var_type)
#                 print('----------------------------')
                
                #------------------------------------------------
                # For layered variables (e.g. soil layers) call
                # initialize_layered_vars().  (11/15/16)
                #------------------------------------------------
                if (var_name == 'n_layers'):
                    self.n_layers = np.int32( value )  # (will be repeated below)
                    self.initialize_layer_vars()

                #----------------------------------------------
                # Does var_name end with an array subscript ?
                # This is the case for soil layer variables.
                #----------------------------------------------
                p1 = var_name.rfind('[')
                p2 = var_name.rfind(']')
                if (p1 > 0) and (p2 > p1):
                    var_base  = var_name[:p1]
                    subscript = var_name[p1:p2+1]
                    var_name_file_str = var_base + '_file' + subscript
                else:
                    var_base = var_name
                    var_name_file_str = var_name + '_file'

                #----------------------------------------------------------
                # If previous variable ended with "_type", then assume
                # its purpose was to set the data type for this variable.
                #----------------------------------------------------------
                if (last_var_name.startswith(var_base + '_type')):
                    #--------------------------------------------------
                    # For Python 2.7, this worked with 2nd arg = "{}".
                    # For 3.7, didn't work with "{}" or "globals()".
                    # NB!  You can't use exec() to set locals in a
                    # function, unless the name was *already* assigned
                    # to in that function.  Can use eval(), though.
                    #--------------------------------------------------
                    ## type_choice = ''   # (For exec in Python 3.x.)
                    ## exec( "type_choice = self." + last_var_name, {}, locals() )
                    #---------------------------------------------------------------
                    ### type_choice = eval( "self." + last_var_name )  ## OK TOO.
                    type_choice = last_value  ### (2019-10-03, for Python 3()
                    if (type_choice.lower() == 'scalar'):
                        #--------------------
                        # For Scalar option
                        #--------------------
                        ### print('#### var_name_file_str = ' + var_name_file_str )
                        exec( "self." + var_name_file_str + " = ''", {}, locals() )
                        READ_SCALAR = True
                    elif (type_choice.lower() == 'time_series'):
                        #------------------------------------------
                        # For Time Series option (read from file)
                        #------------------------------------------
                        # dtype will default to 'float64'
                        ## exec( 'self.' + var_name + '= self.initialize_scalar(0)', {}, locals() )

                        exec( "self." + var_name + " = 0.0", {}, locals() )  # (just a placeholder)
                        READ_FILENAME = True
                    else:
                        #--------------------------------------------------------
                        # For Grid & Grid Sequence options (read from file)
                        # Can't call initialize_grid(), don't know nx & ny yet.
                        #--------------------------------------------------------
                        ### exec( 'self.' + var_name + '= self.initialize_grid(0)', {}, locals() )

                        exec( "self." + var_name + " = 0.0", {}, locals() )  # (just a placeholder)
                        READ_FILENAME = True

                #-----------------------------------           
                # Read a value of type "var_type"
                #-----------------------------------
                # Convert scalars to numpy scalars
                #-----------------------------------
                if (var_type in ['double', 'float']):
                    #-----------------------------------------
                    # Save value of a double or float scalar
                    #-----------------------------------------
                    value = np.float64( value )   # (convert string to float64)
                    ## cmd = 'self.' + var_name + ' = value'
                    cmd = 'self.' + var_name + '= self.initialize_scalar( value )'
                    exec( cmd, {}, locals() )

                    #------------------------
                    # For testing (5/18/12)
                    #------------------------
#                     print 'var_name =', var_name
#                     print 'var_type =', var_type
#                     print 'value    =', value
#                     exec( 'dtype = self.' + var_name + '.dtype', {}, locals() )
#                     print 'dtype    =', dtype
#                     print '---------------------------------'

                elif (var_type in ['long', 'int']):
                    #----------------------------------------
                    # Save value of a short or long integer
                    #----------------------------------------
                    value = np.int32( value )  # (convert string to int32)
                    dtype = 'int32'
                    ## cmd = 'self.' + var_name + ' = value'
                    cmd = 'self.' + var_name + '= self.initialize_scalar( value, dtype=dtype )'
                    exec( cmd, {}, locals() )

                elif (var_type == 'string'):
                    #-----------------------------------------
                    # Get the value string for this var_name
                    #----------------------------------------------------
                    # If string has a placeholder filename prefix, then
                    # expand it here.  Need to use original "var_name"
                    # without appending "_file" until assignment.
                    #----------------------------------------------------
                    # case_str = '<case_prefix>'
                    # site_str = '<site_prefix>'
                    case_str = '[case_prefix]'
                    site_str = '[site_prefix]'
                    #---------------------------------
                    s = value
                    if (s[:13] == case_str):
                        value_str = (self.case_prefix + s[13:])
                    elif (s[:13] == site_str):
                        value_str = (self.site_prefix  + s[13:])
                    else:
                        value_str = s

                    #-----------------------------------------------
                    # If value_str indicates True or False, then
                    # convert to Python boolean.  This is a more
                    # general method than before for handling
                    # yes/no booleans like PRECIP_ONLY that don't
                    # start with 'SAVE_'. Note that both 'Yes' and
                    # 'No' evaluate to True, regardless of case.
                    # 2023-09-01.
                    #-----------------------------------------------
                    # Should we include "1" and "0" if string??
                    # As integers, they evaluate to True & False,
                    # as occurs with FLOOD_OPTION & ATTENUATE
                    # in the channel component CFG file.
                    #-----------------------------------------------
                    true_list  = ['yes', 'true', 'on', '1']
                    false_list = ['no', 'false', 'off', '0']
                    bool_list  = true_list + false_list  # join lists
                    if (s.lower() in bool_list):
                        VALUE_SET = True
                        if (s.lower() in true_list):
                            exec( "self." + var_name + " = True", {}, locals() )
                        else:
                            exec( "self." + var_name + " = False", {}, locals() )
                    else:
                        VALUE_SET = False

                    #-----------------------------------------------
                    # If var_name starts with "SAVE_" and value is
                    # Yes or No, then convert to Python boolean.
                    # See more general method above. 2023-09-01
                    #-----------------------------------------------                   
#                     if (var_name[:5] == 'SAVE_'):
#                         VALUE_SET = True
#                         if (s.lower() == 'yes'):
#                             exec( "self." + var_name + " = True", {}, locals() )
#                         elif (s.lower() == 'no'):
#                             exec( "self." + var_name + " = False", {}, locals() )
#                         else:
#                             VALUE_SET = False
#                     else:
#                         VALUE_SET = False
                    #----------------------------------------------------------
                    if not(VALUE_SET):
                        if (READ_FILENAME):
                            cmd = 'self.' + var_name_file_str + ' = value_str'
                            # print('var_name_file_str = ' + var_name_file_str)
                            # print('cmd = ' + cmd)
                            # print('==========================================')
                            exec( cmd, {}, locals() )  ## For Python 3
                            ### exec( cmd, globals(), locals() )  ## For Python 3
                        elif (READ_SCALAR):
                            ##############################################################
                            # NOTE!  var_type = 'string', so not sure why this is here.
                            ##############################################################
                            ## exec( "self." + var_name + " = np.float64(value_str)")
                            exec( "self." + var_name + " = value_str", {}, locals() )
                        else:
                            exec( "self." + var_name + " = value_str", {}, locals() )
                else:
                    print('ERROR in BMI_base.read_config_file().')
                    print('   Unsupported data type = ' + var_type + '.')
                    print(' ')

                last_var_name = var_name
                last_value    = value     #### (2019-10-03, for Python 3)

        cfg_unit.close()

    #   read_config_file()
    #-------------------------------------------------------------------
    def initialize_config_vars(self, READ_GRID_INFO=True):
   
        #--------------------------------------------------------------
        # Note: EMELI calls bmi.initialize() with full cfg_file name.
        #       bmi.initialize() saves cfg_file into self.
        #       Until read_path_info() below, we don't know:
        #          in_dir, out_dir, site_prefix or case_prefix.
        #-------------------------------------------------------------
                      
        #--------------------------------------
        # Read input variables from CFG files
        #--------------------------------------
        self.read_path_info()   ###### (2016-02-12)
        self.read_time_info()   ###### (2020-01-14)
        #---------------------------------------------
        # May not be needed by all components.
        # Later put grid_info in a CFG vs. RTI file?
        # Was called by each "*_base.py" before.
        #---------------------------------------------
        if (READ_GRID_INFO):
            self.read_grid_info()   # (stores rti in self, adds "da")        
        
        # print '#### CALLING read_config_file()...'
        self.read_config_file()
        # print '#### AFTER read_config_file():'
        # print '#### in_directory  =', self.in_directory
        # print '#### out_directory =', self.out_directory
        # print ' '

        #-------------------------------------
        # Check the directories and prefixes
        #-------------------------------------
        # Defaults are set in __init__() and
        # some may have just been read in.
        #-------------------------------------
        # print '#### CALLING set_directories()...'
        self.set_directories()  # (Moved up: 10/25/11)

        #-------------------------------------------------------
        # A few things to note:
        #
        # (1) Each component's "read_cfg_file()" method
        #     needs to know where to find its CFG file.
        #     We _could_ store the paths to all of a project's
        #     CFG files in the sync_file (created by driver).
        #     However, assuming that all of the CFG files are
        #     in the CWD is simpler and keeps them together.
        #
        # (2) If a user is running an example w/o the GUI,
        #     then it is simplest if the CFG files needed for
        #     the example are stored together with the input
        #     data for the example in a central place on
        #     beach.  It is likely that the user won't have
        #     write permission there, which is why we need a
        #     separate output directory for model output.
        #     We also want to avoid the need for users to
        #     create copies of example CFG files in their own
        #     directories. (But they will have read permission
        #     and can make copies if they want.)
        #
        # (3) Component CFG files represent model input and
        #     may contain input parameters and/or the names
        #     of input files, such as initial-value grids.
        #     If a user is not running an example, then they
        #     will need to create an appropriate set of input
        #     files using whatever tools we give them.
        #
        # (4) Except when running examples (or using someone
        #     else's input files) the directories for input
        #     and output will typically be the same.
        #
        #-------------------------------------------------------        #
        # self.mode is set to "driver" or "nondriver" in the
        # initialize() method.
        #-------------------------------------------------------
        
        #--------------------------------------------
        # Set any input variables that are computed
        #----------------------------------------------
        # NOTE:  This must be called in initialize(), 
        #        AFTER calling read_input_files().
        #----------------------------------------------
        # self.set_computed_input_vars()
     
    #   initialize_config_vars()
    #-------------------------------------------------------------------
    def set_computed_input_vars( self):

        pass
    
    #   set_computed_input_vars()
    #-------------------------------------------------------------------
    def initialize_basin_vars( self ):

        #------------------------------------------------------------
        # This saves outlet_file, outlet_IDs, outlet_ID, n_outlets,
        # basin_area and basin_relief into self.  (9/19/14)
        #------------------------------------------------------------
        outlets.read_outlet_file( self )

    #   initialize_basin_vars()    
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def prepend_directory(self, file_list, INPUT=True):

        #-----------------------------------------------------------
        # (2020-05-03) This was previously used by just a few
        # components, such as channels_base.  But is not ideal
        # because it uses both eval and exec.  So phase out?
        #-----------------------------------------------------------        
        # (11/14/11) Call from a component's open_input_file() 
        # method something like this:
        #   self.prepend_directory( ['slope_file', 'width_file'] ) 
        #-----------------------------------------------------------
        if (INPUT):
            dir_part = " = self.in_directory + "
        else:
            dir_part = " = self.out_directory + "

        for file_str in file_list:
            self_part = "self." + file_str
#             filename = ''   # (For exec in Python 3.x.)
#             exec( 'filename = ' + self_part, {}, locals() )
            filename = eval( self_part )   ### For Python 3.
            if (filename != ''):
                exec( self_part + dir_part + self_part, {}, locals() ) 

    #   prepend_directory
    #-------------------------------------------------------------------
    def set_directories(self, REPORT=False):

        #--------------------------------------------
        # Note:  This was cleaned up on 2020-05-03.
        #--------------------------------------------       
        
        #--------------------------------------------------
        # These should have been set by read_path_info().
        #--------------------------------------------------        
#         if (self.in_directory == None):
#             self.in_directory = os.getcwd() + os.sep
#         #-----------------------------------------------
#         if (self.out_directory == None):
#             self.out_directory = os.getcwd() + os.sep
#         #-----------------------------------------------
#         if (self.site_prefix == None):
#             self.site_prefix = self.cfg_prefix
#         #-----------------------------------------------
#         if (self.case_prefix == None):
#             self.case_prefix = self.cfg_prefix
        
        ## REPORT = True  # (to debug)
        if (REPORT):
            print('Before calling set_directories()...')
            print('cfg_file      =', self.cfg_file)
            print('cur_directory =', os.getcwd() + os.sep )
            print('in_directory  =', self.in_directory)
            print('out_directory =', self.out_directory)
#             print('site_prefix   =', self.site_prefix)
#             print('case_prefix   =', self.case_prefix)
            print()

        #--------------------------------------    
        # Set the cfg_directory from cfg_file
        # It is usually not set before now.
        #--------------------------------------
        if not(hasattr(self, 'cfg_directory')):
            if (self.cfg_file is not None):
                cfg_dir = os.path.dirname(self.cfg_file)
                #---------------------------------------
                # Need this when cfg_file is local, as
                # in a Jupyter notebook (2020-05-28)
                #---------------------------------------
                if (cfg_dir == ''):
                    cfg_dir = os.getcwd()
                cfg_dir = (cfg_dir + os.sep)
                ## print 'cfg_dir =', cfg_dir
        else:
            cfg_dir = self.cfg_directory

        #----------------------------------------------------
        # First, expand user path abbreviations like "~"
        #----------------------------------------------------
        # Note that os.path.expanduser() does not expand
        # relative paths like ".", "..", "./". (2/13/12)
        # And it does NOT remove trailing path separator.
        #--------------------------------------------------
        in_dir  = os.path.expanduser( self.in_directory  )
        out_dir = os.path.expanduser( self.out_directory )

        #------------------------------------------------------------
        # In CFG files, the input directory is often set to ".",
        # which indicates the directory that contains the CFG file.
        # However, once "os.chdir()" is called, "." will expand to
        # something else.  Also, we want to avoid calling
        # "os.chdir()", because it creates problems for finding
        # package paths.  So to address these issues, we expand
        # relative paths like "." & ".." here. (Update: 2020-05-03)
        #------------------------------------------------------------
        # Can use os.path.abspath() or os.path.realpath() to
        # expand relative paths like ".", "..", "../..", etc.
        # os.path.realpath() also expands symbolic links.
        # These DO NOT expand "~".
        # Note that trailing path separator may be removed.
        # cfg_dir already ends in os.sep.
        #------------------------------------------------------------
        if (in_dir[0] == '.'):
            in_dir  = os.path.abspath( cfg_dir + in_dir  )
        if (out_dir[0] == '.'):
            out_dir = os.path.abspath( cfg_dir + out_dir )

        #--------------------------------------------------
        # Add trailing separator to directory, if missing
        # Note: Must come AFTER os.path.realpath calls.
        #--------------------------------------------------
        if (in_dir[-1] != os.sep):
            in_dir += os.sep
        #----------------------------
        if (out_dir[-1] != os.sep):
            out_dir += os.sep

        #-----------------------------------------------------
        # New idea of separate "static var" directories
        # in input_directory, called "soil", "topo", "met".
        # If present, these can be used to keep the input
        # files better organized.
        #-----------------------------------------------------
        # Most components have "open_input_files()" methods
        # where in input directory is used. A soil directory
        # would mostly be used by infil components.
        # A topo directory would mostly be used by the
        # channels and meteorology components.
        #-----------------------------------------------------
        # Could add a "snow_dir", and maybe a directory for
        # every process: chan, infil, evap, snow, etc.
        #-----------------------------------------------------        
        topo_dir = in_dir + '__topo' + os.sep
        soil_dir = in_dir + '__soil' + os.sep
        met_dir  = in_dir + '__met'  + os.sep
        misc_dir = in_dir + '__misc' + os.sep
        #----------------------------------------
        if not(os.path.exists( topo_dir )):
            topo_dir = in_dir
        if not(os.path.exists( soil_dir )):
            soil_dir = in_dir
        if not(os.path.exists( met_dir )):
            met_dir  = in_dir  # (see note below)
        if not(os.path.exists( misc_dir )):
            misc_dir = in_dir
        #------------------------------------------
        # Not static, so use cfg_dir not in_dir ?
        #------------------------------------------       
        # if not(os.path.exists( met_dir )):
        #     met_dir = cfg_dir

        #---------------------------------------------------- 
        # These could also be set in the path_info CFG file    
        #----------------------------------------------------        
#         if not(hasattr(self, 'soil_directory')):
#             soil_dir = in_dir
#         if not(hasattr(self, 'topo_directory')):
#             topo_dir = in_dir
#         if not(hasattr(self, 'met_directory')):
#             met_dir = cfg_dir
        #-----------------------------------------------------

        #------------------------------------------------
        # Save fully expanded directory names into self
        #------------------------------------------------
        self.cfg_directory  = cfg_dir
        self.in_directory   = in_dir
        self.out_directory  = out_dir
        self.topo_directory = topo_dir
        self.soil_directory = soil_dir
        self.met_directory  = met_dir
        self.misc_directory = misc_dir                                         
        if (REPORT):
            print('After calling set_directories()...')
            print('cfg_directory  =', cfg_dir)
            print('in_directory   =', in_dir)
            print('out_directory  =', out_dir)
            print('topo_directory =', topo_dir)
            print('soil_directory =', soil_dir)
            print('met_directory  =', met_dir)
            print('misc_directory =', misc_dir)
#             print('site_prefix   =', self.site_prefix)
#             print('case_prefix   =', self.case_prefix)
            print()
               
    #   set_directories()
    #-------------------------------------------------------------------
    def initialize_scalar(self, value=0.0, dtype='float64'):

        #------------------------------------------------------------
        # Note:  See utils/scalar_tests.py for more info.
        #------------------------------------------------------------
        # If we create scalars as 0D numpy arrays, then:
        #    - Return var is MUTABLE and other components with
        #      a reference to it WILL see it change.
        #    - Values will print without square brackets.
        #    - Trying to subscript it will generate the error:
        #      "IndexError: 0-d arrays can't be indexed"
        #    - In-place updates will preserve the reference, e.g.:
        #         x.fill(3), x += 1, x *= 2
        #         Use optional "out" argument to numpy ufuncs:
        #            e.g. np.sin(x, out=x)
        #    - Non in-place updates will BREAK the reference, e.g.:
        #         x = x+1, x = x*2, x = np.sin(x)
        #    - Data type changes will BREAK the reference, and may
        #      require more or less memory, e.g.:
        #         x = np.float32(x)
        #    - We CAN compare this type of scalar to others with
        #      "<", ">", "==", etc.
        #    - We can get the data type with ".dtype".
        #    - The data type will be "numpy.ndarray".
        #    - numpy.rank(x) will return 0.
        #------------------------------------------------------------        
        return np.array(value, dtype=dtype)

        #------------------------------------------------------------      
        # If we create scalars as 1D numpy arrays, then:
        #    - They are really 1D arrays with one element.
        #    - Return var is MUTABLE and other components with
        #      a reference to it WILL see it change.
        #    - Values will print WITH square brackets unless we
        #      ask for element 0 "[0]".
        #    - In-place updates will preserve the reference, e.g.:
        #         x[0]=3, x[:]=3, x[0] += 1, x[0] *= 2
        #         Use optional "out" argument to numpy ufuncs:
        #            e.g. np.sin(x, out=x)
        #    - Non in-place updates will BREAK the reference, e.g.:
        #         x = x+1, x = x*2, x = np.sin(x)
        #    - Data type changes will BREAK the reference, and may
        #      require more or less memory, e.g.:
        #         x = np.float32(x)
        #    - We CAN compare this type of scalar to others with
        #      "<", ">", "==", etc.
        #    - We can get the data type with ".dtype".
        #    - The data type will be "numpy.ndarray".
        #    - numpy.rank(x) will return 1.
        #------------------------------------------------------------ 
        ## return np.array([value], dtype=dtype)

        #------------------------------------------------------------
        # If we create scalars as typed numpy scalars, then:
        #    -  Return var is now IMMUTABLE and other components
        #       with a reference to it WILL NOT see it change.
        #    - We CAN compare this type of scalar to others with
        #      "<", ">", "==", etc.
        #    - We can get the data type with ".dtype".
        #    - The data type will be "numpy.float64".
        #    - numpy.rank(x) will return 0.
        #------------------------------------------------------------
        # return np.float64( value )
        
    #   initialize_scalar()
    #-------------------------------------------------------------------
    def initialize_grid(self, value=0.0, dtype='float64'):

        #----------------------------------------------------
        # Written on 11/13/16 to match initialize_scalar().
        #----------------------------------------------------
        grid = np.zeros([self.ny, self.nx], dtype=dtype)
        if (value != 0.0):
            grid += value

        return grid

    #   initialize_grid()
    #-------------------------------------------------------------------
    def initialize_var(self, var_type, value=0.0, dtype='float64'):

        #--------------------------------------------------------
        # Note:  if (var_type.lower() == 'scalar'), then do not
        #        call this function;  var already initialized
        #        by read_config_file().
        #--------------------------------------------------------
        if (var_type.lower() == 'time_series'):
            var = self.initialize_scalar( value=value, dtype=dtype)
        else:
            var = self.initialize_grid( value=value, dtype=dtype)

        return var

    #   initialize_var()
    #-------------------------------------------------------------------
    def update_var(self, var_name, value, factor=1.0):

        #---------------------------------------------------------------
        # Note:  Added factor keyword on 2022-02-15.
        #---------------------------------------------------------------
        # Note:  Vars must be initialized first, e.g. using either
        #        initialize_scalar() or initialize_grid() in BMI_base.
        #---------------------------------------------------------------
        if (value is None):
            return

        #--------------------------------------------------
        # Support for new "factor" keyword (2022-02-15)
        # This function is called from read_input_files()
        # and it is okay to modify the passed "value".
        #--------------------------------------------------
        if (factor != 1.0):
            value *= factor
        
        if (self.is_scalar( var_name )):
            #----------------------------------------------------
            # Update the value of a scalar without breaking the
            # reference so other components will see it change.
            #----------------------------------------------------
            # See Notes for initialize_scalar() above.
            # This is needed for type "Time Series".
            #----------------------------------------------------
            exec( 'self.' + var_name + '.fill( value )', {}, locals() )
        else:
            #----------------------------------------------
            # Update the value of a grid (2D, 3D) without
            # breaking the reference (or making a copy).
            #----------------------------------------------
            exec( 'self.' + var_name + '[:] = value', {}, locals() )

    #   update_var()
    #-------------------------------------------------------------------
    # These are for convenience;  not part of BMI.
    #-------------------------------------------------------------------
    def is_scalar(self, var_name):

        #--------------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        # globals() is needed for "np", & locals() for var_name
        #--------------------------------------------------------
#         n = -1   # (For exec in Python 3.x.) 
#         exec("n = np.ndim(self." + var_name + ")", {}, locals() )      
#         ## exec("n = np.ndim(self." + var_name + ")", globals(), locals() )       
#         return (n == 0)
    
        n = eval( "np.ndim(self." + var_name + ")" )
        return (n == 0)
        
    #   is_scalar()
    #-------------------------------------------------------------------
    def is_vector(self, var_name):

        #--------------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        # globals() is needed for "np", & locals() for var_name
        #--------------------------------------------------------
#         n = -1   # (For exec in Python 3.x.) 
#         exec("n = np.ndim(self." + var_name + ")", {}, locals() )       
#         ## exec("n = np.ndim(self." + var_name + ")", globals(), locals() )       
#         return (n == 1)

        n = eval( "np.ndim(self." + var_name + ")" )
        return (n == 1)
            
    #   is_vector()
    #-------------------------------------------------------------------
    def is_grid(self, var_name):

        #------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        #------------------------------------------------ 

        #-------------------------------------------------
        # (9/29/09) This might be causing a problem with
        # the c++ bindings for this CCA component.
        #-------------------------------------------------         
##        exec("type_str = str(type(self." + var_name + "))", {}, locals() )
##        p1 = type_str.find("ndarray")
##        p2 = type_str.find("float")
##        if (p1 == -1) and (p2 == -1):
##            print 'ERROR: type(' + var_name + ') =' + type_str
##            return False
        #-------------------------------------------------
        # (9/29/09) This might be causing a problem with
        # the c++ bindings for this CCA component.
        #-------------------------------------------------        
##        if ("ndarray" not in type_str) and \
##           ("float" not in type_str):
##            print 'ERROR: type(' + var_name + ') =' + type_str
##            return False
        #------------------------------------------------------- 
        # globals() is needed for "np", & locals() for var_name
        #--------------------------------------------------------
#         n = -1   # (For exec in Python 3.x.)
#         exec("n = np.ndim(self." + var_name + ")", {}, locals() )       
#         ## exec("n = np.ndim(self." + var_name + ")", globals(), locals() )
#         return (n == 2)

        n = eval( "np.ndim(self." + var_name + ")" )
        return (n == 2)
        
    #   is_grid()
    #---------------------------------------------------------------------

    
    
