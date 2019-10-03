
## Copyright (c) 2001-2016, Scott D. Peckham
## January 2009  (converted from IDL)
## May, July, August, October 2009
## March 2010 (modified to support local timestepping)
## September 2011  (Updated unit tests.)
## October 2011 (prepared for use as a CMT component)
## January 2012 (Added set_default_config_vars().)
## January 2012 (Removed call to update_noflow_IDs() and
##               removed GLOBAL implementations of many methods.
##
#####################################################
# (3/4/10)  Replaced ".flow_grid" with ".d8_grid"
#####################################################
## NB!  TF expects d8.codes vs. d8.d8_grid
## NB!  Need to update a few more functions to use g = "code_list".

#############################################################
## Note:  channel_base.py calls many of these but then
##        saves the results among its state variables.
#############################################################

import numpy as np

import os        # (for os.chdir(), in unit_test())
import os.path
import time

from topoflow.utils import BMI_base
from topoflow.utils import fill_pits
from topoflow.utils import model_output
from topoflow.utils import pixels
from topoflow.utils import rtg_files

#-------------------------------------------
# For use outside of the TopoFlow package.
#-------------------------------------------
# import BMI_base
# import fill_pits
# import model_output     # (added: 11/8/11)
# import pixels
# import rtg_files
## import tf_utils  # (not used now)

#---------------------------------------------------------------------
#
#   class d8_base
#
#       get_attribute()           # (NOT READY YET)  #############
#       get_input_var_names()
#       get_output_var_names()
#       get_var_name()
#       get_var_units()
#-----------------------------
#       set_constants()
#       set_default_config_vars()   # (1/17/12)
#----------------------------------
#       initialize()
#       update()
#       set_computed_input_vars()
#       initialize_computed_vars()
#----------------------------------
#       get_pixel_dimensions()
#       get_flow_code_list()
#       get_flow_code_list_opps()
#       get_valid_code_map()
#       get_ID_grid()
#       get_parent_inc_map()
#       get_edge_IDs()
#       get_not_edge_grid()
#       resolve_array_cycle()
#       get_resolve_array()
#----------------------------------
#       update_parent_ID_grid()     # (can handle periodic BCs) 
#       update_parent_IDs()         # (used by erosion_base.update_slope_grid())
#       update_flow_from_IDs()
#       update_flow_to_IDs()
#       update_noflow_IDs()
#----------------------------------
#       read_flow_grid()
#       update_flow_grid()
#          start_new_d8_codes()
#          break_flow_grid_ties()
#          link_flats()
#----------------------------------
#       update_flow_width_grid()
#       update_flow_length_grid()
#       update_area_grid()          # (added on 10/28/09)
#----------------------------------
#       change_extension_to_rtg()   # (9/20/11)    ( NB!  RTS_FILES = True )
#       open_input_files()          # (11/8/11) 
#       read_input_files()          # (11/8/11)
#       close_input_files()         # (11/8/11)
#---------------------------------
#       update_outfile_names()      # (9/20/11)
#       open_output_files()         # (11/8/11)
#       write_output_files()        # (11/8/11)
#       close_output_files()        # (11/8/11)
#       save_grids()                # (9/20/11)
#       save_pixel_values()         # (11/8/11)
#
#----------------------------------
#       update_area_grid_OLD()      # (can't handle periodic BCs)
#       get_flow_width_grid()       # OBSOLETE ?
#       get_flow_length_grid()      # OBSOLETE ?
#
#-----------------------------------------------------------------------
class d8_component( BMI_base.BMI_component ):

    _input_var_names = [ 'land_surface__elevation' ]  # (DEM)
    
    #------------------------------------------------------------
    # Note: no_flow_IDs, parent_IDs, w1, p1, etc. are currently
    #       not 1D arrays but tuples with array of rows and an
    #       array of cols.  Therefore they can't be accessed by
    #       another component using current interface functions
    #       for getting and setting 1D arrays or "vectors".
    #------------------------------------------------------------        
    # items = ['ID_grid', 'd8_grid', 'A', 'ds', 'dw',
    # 'parent_IDs', 'edge_IDs',
    # 'w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w7', 'w8',
    # 'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8',
    # 'noflow_IDs', nx, ny, dx, dy, dd, da ]
       
    _output_var_names = []
    
    _var_name_map = {
        'land_surface__elevation': 'DEM' }

    _var_units_map = {
        'land_surface__elevation': 'm' }

    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_input_var_names(self):

        #--------------------------------------------------------
        # Note: These are currently variables needed from other
        #       components vs. those read from files or GUI.
        #--------------------------------------------------------   
        return self._input_var_names
    
    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):
 
        return self._output_var_names
    
    #   get_output_var_names()
    #-------------------------------------------------------------------
    def get_var_name(self, long_var_name):
            
        return self._var_name_map[ long_var_name ]

    #   get_var_name()
    #-------------------------------------------------------------------
    def get_var_units(self, long_var_name):

        return self._var_units_map[ long_var_name ]
   
    #   get_var_units()
    #-------------------------------------------------------------------
##    def get_var_type(self, long_var_name):
##
##        #---------------------------------------
##        # So far, all vars have type "double",
##        # but use the one in BMI_base instead.
##        #---------------------------------------
##        return 'float64'
##    
##    #   get_var_type()                    
    #-------------------------------------------------------------------
    def set_constants(self):

        #------------------------
        # Define some constants
        #------------------------
        self.dt     = 1.0   # (needs to be defined)
        self.nodata = np.float32(-9999)

        #----------------------------------------
        # (9/19/11) Put these here because they
        # don't seem to be anywhere else.
        #----------------------------------------
        self.BREAK_TIES = True
        # self.LINK_FLATS = True       # (this is read from CFG file)
        # self.FILL_PITS_IN_Z0 = True  # (this is read from CFG file)
        
    #   set_constants()  
    #-------------------------------------------------------------------
    def set_default_config_vars(self):

        #-------------------------------------------------------------
        # Notes:  D8 components are used as utility classes in 
        #         components like Erode_D8_Global and DEM_Smoother,
        #         but are not then connected through a port in the
        #         CMT.  In these cases, a CFG file is not created
        #         for the D8 component, and this function is needed
        #         in order to set the config vars.  This is done
        #         before trying to find and read the missing CFG
        #         file (which gives a warning but continues).
        #         There is one entry here for each config var that
        #         is set in the CFG file.
        #-------------------------------------------------------------
        self.comp_status      = 'Enabled'
        if not(hasattr(self, 'in_directory')):
            self.in_directory = '~/CMT_Output/'  # (Is this used ??)
        self.out_directory    = '~/CMT_Output/'
        if not(hasattr(self, 'site_prefix')):
            self.site_prefix = 'NOT_SET'
        self.case_prefix      = 'NOT_SET'
        self.n_steps          = 1
        self.dt               = 1.0
        if not(hasattr(self, 'DEM_file')):
            self.DEM_file = 'NOT_SET'
        if not(hasattr(self, 'A_units')):
            # (1/17/11) Need "km^2" for DEM_Smoother ??
            self.A_units          = 'km^2'
        self.LINK_FLATS       = 1
        if not(hasattr(self, 'FILL_PITS_IN_Z0')):
            # (1/17/12) DEM_file is overwritten if this is set to 1.
            # Currently fails for KY_Sub, maybe due to INT, FLOAT data type.
            self.FILL_PITS_IN_Z0 = 0  
        self.LR_PERIODIC      = 0
        self.TB_PERIODIC      = 0
        self.save_grid_dt     = 1.0
        self.SAVE_AREA_GRIDS  = False
        self.area_gs_file     = ''
        self.SAVE_CODE_GRIDS  = False
        self.code_gs_file     = ''
        self.SAVE_DS_GRIDS    = False
        self.ds_gs_file       = ''
        self.SAVE_DW_GRIDS    = False
        self.dw_gs_file       = ''
        self.save_pixels_dt   = 1.0
        self.pixel_file       = 'NOT_SET'
        self.SAVE_AREA_PIXELS = False
        self.area_ts_file     = ''
        self.SAVE_CODE_PIXELS = False
        self.code_ts_file     = ''
        self.SAVE_DS_PIXELS   = False
        self.ds_ts_file       = ''
        self.SAVE_DW_PIXELS   = False
        self.dw_ts_file       = ''

    #   set_default_config_vars()
    #-----------------------------------------------------------------
    def initialize(self, cfg_file=None, mode='nondriver',
                   SILENT=False, REPORT=True):

        #--------------------------------------------------------
        # Note:  This function calls functions that compute a
        #        variety of D8 variables and saves them in its
        #        state.  This "d8 state" can be embedded within
        #        another class such as the "channels_base" or
        #        "erosion_base" class.
        #--------------------------------------------------------
        if not(SILENT):
            print('D8 component: Initializing...')

        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode

        if (cfg_file == None):
            #-----------------------------------------------------------------
            # Changed from using site_prefix to case_prefix on 2/11/2017
            # since this is the only CFG file that uses site_prefix.
            # See matching change in initialize_d8_vars in channels_base.py.
            #-----------------------------------------------------------------
            cfg_extension = self.get_attribute( 'cfg_extension' )  #########
            # print 'cfg_extension =', cfg_extension
            # print 'case_prefix = ', self.case_prefix

            
            filename      = self.case_prefix + cfg_extension
            #filename      = self.site_prefix + cfg_extension
            cfg_file      = self.in_directory + filename
            ## self.cfg_file   = os.path.join( os.getcwd(), filename )
            ## self.cfg_prefix = self.site_prefix
        self.cfg_file = cfg_file
        
        ## print 'In d8_base.initialize(), cfg_file =', cfg_file  #########

        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        # print '### Calling set_constants()...'
        self.set_constants()
        #----------------------------------------------------------------
        # (1/17/12) EXPERIMENT to solve problem with missing CFG files.
        #----------------------------------------------------------------
        # print '### Calling set_default_config_vars()...'
        self.set_default_config_vars()
        #----------------------------------------------------------------
        # print '### Calling initialize_config_vars()...'
        self.initialize_config_vars()
        # print '### Calling read_grid_info()...'
        self.read_grid_info()   # (also gets & stores self.da)
        self.initialize_basin_vars()    # (uncommented on 11/8/11.)

        #-------------------------------------------
        # This must come before "Disabled" test ??
        #-------------------------------------------
        # print '### Calling initialize_time_vars()...' 
        self.initialize_time_vars()
        ## self.initialize_time_vars( units='years' )

        #--------------------------------------------
        # Convert units for da from m^2 to km^2 ??
        #--------------------------------------------
        # Better to do this in "update_area_grid()"
        # See "unit_test()".
        #--------------------------------------------        
##        if ('km' in self.A_units.lower()):
##            self.da = self.da / 1e6
  
        #-----------------------------------------------
        # These return values that don't depend on the
        # flow grid and don't change, so they should
        # simply be stored for subsequent use.
        #-----------------------------------------------
        # Values that do depend on the flow grid are
        # computed by calls in the update() method.
        #-----------------------------------------------        
        self.get_pixel_dimensions(SILENT=SILENT, REPORT=REPORT)
        self.get_flow_code_list()
        self.get_flow_code_list_opps()
        self.get_ID_grid(SILENT=SILENT)
        self.get_parent_inc_map()
        self.get_edge_IDs(SILENT=SILENT)
        self.get_not_edge_grid(SILENT=SILENT)
        #-------------------------------------------------
        self.get_resolve_array(SILENT=SILENT)   ######
        self.get_valid_code_map(SILENT=SILENT)
        
        #-----------------------------
        # Initialize dw, ds, A, etc.
        #-----------------------------
        ## self.initialize_time_vars()  # (done above now)
        self.initialize_computed_vars(SILENT=SILENT, REPORT=REPORT)
        
        self.open_output_files()    # (added on 11/8/11) ###################
        self.status = 'initialized'
        
##        finish = time.time()
##        print 'Run time for initialize =', (finish - start), ' [secs]'
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, time=None, DEM=None,
               SILENT=True, REPORT=False):

        self.status = 'updating'  # (OpenMI 2.0 convention)
        if not(SILENT):
            print('D8 component: Updating...')
        # if (self.mode == 'main'):
        #     self.print_time_and_value(self.z_outlet, 'z_out', '[m]')
        
        #-------------------------
        # Update computed values
        #------------------------------------------------------
        # Note that ds grid is used to compute updated slopes
        #------------------------------------------------------
        ## fill_pits.fill_pits()   ## pass DEM to here? ## 

        self.update_flow_grid(DEM=DEM,
                              SILENT=SILENT, REPORT=REPORT)
        self.update_parent_ID_grid()
        self.update_parent_IDs()     # (needed for gradients)
        self.update_flow_from_IDs()
        self.update_flow_to_IDs()
        #-----------------------------------------------------------
        # Next line was removed because it was hurting performance
        # of erode_d8_global.py and erode_d8_local.py even though
        # "noflow_IDs" were not being used. (1/25/12)
        #-----------------------------------------------------------
        ### self.update_noflow_IDs()
        self.update_flow_width_grid(SILENT=SILENT, REPORT=REPORT)   # (dw)
        self.update_flow_length_grid(SILENT=SILENT, REPORT=REPORT)  # (ds)
        self.update_area_grid(SILENT=SILENT, REPORT=REPORT)

        #-------------------------------------------
        # Read from files as needed to update vars 
        #-------------------------------------------
        # if (self.time_index > 0):
        #     self.read_input_files()

        #----------------------------------------------
        # Use this for saving D8 flow and area grids
        #------------------------------------------------
        # Write user-specified data to output files ?
        #------------------------------------------------
        self.write_output_files( time )    # (uncommented on 11/8/11)

        #------------------------
        # Check computed values
        #------------------------
##        if (OK):
##            self.status = 'updated'  # (OpenMI 2.0 convention)
##        else:
##            self.status = 'failed'
##            self.DONE   = True

        #------------------------
        # Update internal clock
        #------------------------
        self.update_time()
        self.status = 'updated'  # (OpenMI 2.0 convention)
            
    #   update()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        self.LR_PERIODIC  = (self.LR_PERIODIC != 0)
        self.TB_PERIODIC  = (self.TB_PERIODIC != 0)
        self.ALL_PERIODIC = (self.LR_PERIODIC and self.TB_PERIODIC)

        #-------------------------------------
        # These are set elsewhere now.
        # See d8_local.initialize_d8_vars()
        #-------------------------------------
        ## self.LINK_FLATS = (self.LINK_FLATS != 0)
        ## self.BREAK_TIES = (self.BREAK_TIES != 0)

        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        # self.save_grid_dt   = maximum(self.save_grid_dt,   self.dt)
        # self.save_pixels_dt = maximum(self.save_pixels_dt, self.dt)
        
    #   set_computed_input_vars()
##    #-------------------------------------------------------------------
    def initialize_computed_vars(self, DOUBLE=False,
                                 SILENT=True, REPORT=False):

        self.RT3_TEST = False  # (default setting)
        nx = self.nx  # (Local synonyms)
        ny = self.ny

        if (DOUBLE):
            self.dw = np.zeros([ny, nx], dtype='Float64')
            self.ds = np.zeros([ny, nx], dtype='Float64')
            self.A  = np.zeros([ny, nx], dtype='Float64')
        else:
            self.dw = np.zeros([ny, nx], dtype='Float32')
            self.ds = np.zeros([ny, nx], dtype='Float32')
            self.A  = np.zeros([ny, nx], dtype='Float32')

        #-----------------------------------
        # Read DEM from DEM_file (11/7/11)
        #----------------------------------------------
        # This was removed from start_new_d8_codes().
        #
        # NB! start_new_d8_codes() has a DEM argument
        # that is used by erode_d8_global.py, etc.
        #----------------------------------------------
        self.DEM_file = (self.in_directory + self.DEM_file)  ######## (11/11/16, QUICK FIX)
        if not( os.path.exists(self.DEM_file) ):
            print('ERROR: Could not find DEM file:')
            print(self.DEM_file)
            return
        DEM = rtg_files.read_grid( self.DEM_file, self.rti, SILENT=False )
        if not(SILENT):
            print('   min(DEM), max(DEM) =', DEM.min(), DEM.max())
        self.DEM = DEM

        #---------------------------------------------------
        # Option to fill pits in the initial DEM (11/7/11)
        #---------------------------------------------------
        if (self.FILL_PITS_IN_Z0):
            print(' ')
            print('Filling pits in initial DEM...')
            fill_pits.fill_pits( self.DEM, 'FLOAT', self.nx, self.ny,
                                 SILENT=SILENT )
            #--------------------------
            # Save new DEM to a file. 
            #-------------------------- 
            rtg_files.write_grid( self.DEM, self.DEM_file, self.rti)

            ###################################################
            #  THE ABOVE CALL OVERWRITES THE ORIGINAL DEM !!
            ###################################################

    #   initialize_computed_vars() 
    #---------------------------------------------------------------------
    def get_pixel_dimensions(self, DOUBLE=False,
                             SILENT=True, REPORT=False):

        if not(SILENT):
            print('Computing pixel dimensions...')

        dx, dy, dd = pixels.get_sizes_by_row(self.rti, METERS=True)
        self.dx = dx
        self.dy = dy
        self.dd = dd

        if not(DOUBLE):
            self.dx = np.float32(self.dx)
            self.dy = np.float32(self.dy)
            self.dd = np.float32(self.dd)

        #------------------------------------------------
        # Get grid cell areas, "da", which is either a
        # scalar (if same for all grid cells) or a grid
        #------------------------------------------------
        # self.da = pixels.get_da( self.rti )
        #---------------------------------------
        # Note that "da" is saved into self by
        # self.read_grid_info().
        #---------------------------------------
        if (REPORT):
            w = 8
            dx_str  = str(self.dx.min()).ljust(w) + ', '
            dx_str += str(self.dx.max()).ljust(w) + ' [m]'
            print('    min(dx), max(dx) = ' + dx_str)
            #--------------------------------------------------           
            dy_str  = str(self.dy.min()).ljust(w) + ', '
            dy_str += str(self.dy.max()).ljust(w) + ' [m]'
            print('    min(dy), max(dy) = ' + dy_str)
            #--------------------------------------------------
            dd_str  = str(self.dd.min()).ljust(w) + ', '
            dd_str += str(self.dd.max()).ljust(w) + ' [m]'
            print('    min(dd), max(dd) = ' + dd_str)
            #--------------------------------------------------
            da_str  = str(self.da.min()).ljust(w) + ', '
            da_str += str(self.da.max()).ljust(w) + ' [m^2]'
            print('    min(da), max(da) = ' + da_str)
      
##            print '    min(dx) =', self.dx.min(),  ' [m]'
##            print '    max(dx) =', self.dx.max(),  ' [m]'
##            #------------------------------------------------
##            print '    min(dy) =', self.dy.min(),  ' [m]'
##            print '    max(dy) =', self.dy.max(),  ' [m]'
##            #------------------------------------------------
##            print '    min(dd) =', self.dd.min(),  ' [m]'
##            print '    max(dd) =', self.dd.max(),  ' [m]'
##            #------------------------------------------------
##            print '    min(da) =', self.da.min(),  ' [m^2]'
##            print '    max(da) =', self.da.max(),  ' [m^2]'

    #   get_pixel_dimensions()
    #---------------------------------------------------------------------
    def get_flow_code_list(self, ARC=False):

        #--------------------------------------------------
        # Notes: RT/Jenson84 flow codes  = | 64 128 1 |
        #                                  | 32  x  2 |
        #                                  | 16  8  4 |

        #        ARC/INFO codes          = | 32 64 128 |
        #                                  | 16  x   1 |
        #                                  |  8  4   2 |
        #--------------------------------------------------------
        # (9/23/11) Looks like we could use uint8 vs. int16
        #           int16 here.  Then see get_valid_code_map().
        #--------------------------------------------------------
        if not(ARC):   
            self.code_list = np.int16([1, 2, 4, 8, 16, 32, 64, 128])
        else:    
            self.code_list = np.int16([128, 1, 2, 4, 8, 16, 32, 64])
        
    #   get_flow_code_list()
    #---------------------------------------------------------------------
    def get_flow_code_list_opps(self, ARC=False):

        #--------------------------------------------------------
        # (9/23/11) Looks like we could use uint8 vs. int16
        #           int16 here.  Then see get_valid_code_map().
        #--------------------------------------------------------
        if not(ARC):    
            self.code_opps     = np.int16([16, 32, 64, 128, 1, 2, 4, 8])
            self.code_opps_3x3 = np.int16([[4,8,16],[2,0,32],[1,128,64]])
        else:    
            self.code_opps     = np.int16([8, 16, 32, 64, 128, 1, 2, 4])
            self.code_opps_3x3 = np.int16([[2,4,8],[1,0,16],[128,64,32]])
        
    #   get_flow_code_list_opps()
    #---------------------------------------------------------------------
    def get_valid_code_map(self, SILENT=True):

        #-----------------------------------------------------------------
        # Notes: This map is used near the end of update_flow_grid()
        #        to set any invalid flow code to 0, which signifies
        #        that the flow direction for that grid cell is undefined.
        #-----------------------------------------------------------------
        self.valid_code_map = np.zeros([256], dtype='UInt8')
        self.valid_code_map[ self.code_list ] = np.uint8( self.code_list )

    #   get_valid_code_map()
    #---------------------------------------------------------------------
    def get_ID_grid(self, SILENT=True):

        #-----------------------------------------------------
        # Get a grid which for each grid cell contains the
        # calendar-style index (or ID) of that grid cell.
        #-----------------------------------------------------
        if not(SILENT):
            print('Computing pixel IDs...')
        
        nx = self.nx
        ny = self.ny
        self.ID_grid = np.reshape(np.arange(nx*ny, dtype='Int32'), [ny, nx])
        
    #   get_ID_grid() 
    #---------------------------------------------------------------------
    def get_parent_inc_map(self):

        #--------------------------------------------
        # Note: parent_ID = ID + incs[flow_code].
        #--------------------------------------------
        # But this doesn't work with periodic BC's.
        #--------------------------------------------
        nx   = self.nx
        incs = np.array([-nx + 1, 1, nx + 1, nx, nx - 1,
                            -1, -nx - 1, -nx], dtype='Int32')
        
        MAP = np.zeros([129], dtype='Int32')
        MAP[self.code_list] = incs
        
        self.inc_map = MAP
        
    #   get_parent_inc_map()       
    #-------------------------------------------------------------------
    def get_edge_IDs(self, SILENT=True):

        if not(SILENT):
            print('Computing edge pixel IDs...')

        #------------------------------------------
        # Get IDs of edge pixels, making sure not
        # to double-count the corner pixels
        #------------------------------------------
        nx    = self.nx
        ny    = self.ny
        T_IDs = np.arange(nx, dtype='Int32')
        B_IDs = T_IDs + (ny - 1) * nx
        L_IDs = (1 + np.arange(ny - 2, dtype='Int32')) * nx
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
        
    #   get_edge_IDs()
    #---------------------------------------------------------------------
    def get_not_edge_grid(self, SILENT=True):

        if not(SILENT):
            print('Computing "not_edge" grid...')

        self.not_edge_grid = np.ones([self.ny, self.nx],
                                        dtype='UInt8')
        self.not_edge_grid[ self.edge_IDs ] = 0
        
##        self.not_edge_grid[:, 0]      = 0
##        self.not_edge_grid[:, nx - 1] = 0
##        self.not_edge_grid[0, :]      = 0
##        self.not_edge_grid[ny - 1, :] = 0

    #   get_not_edge_grid()
    #---------------------------------------------------------------------
    def resolve_array_cycle(self, w, t):

        #--------------------------------------------------------
        #NOTES:  This function takes a set of "byte-coded ties"
        #        and returns all of the "byte-coded ties" that are
        #        equivalent up to a rotation.  This ensures that
        #        consistent tie-breaker directions are assigned
        #        in producing the array called "resolve".
        #        w is the array to be "Cycled."
        #        t is the number of cycles.
        #        v is the "Cycled" ARRAY.
        #--------------------------------------------------------
        v = w
        for j in range(1, t+1):
            v = (v * 2)
            v += (v > 255)
            v = np.bitwise_and(v, 255)
            ##### v = np.logical_and(v, 255)   
        return v
        
    #   resolve_array_cycle()
    #---------------------------------------------------------------------
    def get_resolve_array(self, SILENT=True, REPORT=False):

        #------------------------------------------------------
        # NOTES:  RT/SJ D8 flow direction codes are given by:

        #       |64  128  1|
        #       |32   x   2|
        #       |16   8   4|
        #
        # (9/20/11) Confirmed that this agrees with RT3.
        #------------------------------------------------------
        if not(SILENT):
            print('Computing "resolve array"...')
            
        resolve = np.zeros([256], dtype='UInt8')
        
        #-----------------------------------
        # High-symmetry groups are:
        # {17,34,68,136},{51,102,204,153},
        # {119,238,193,187},{85,170},{255}
        #-----------------------------------
        resolve[[17, 85]]                  = 1
        resolve[[34, 102, 119, 238, 170]]  = 2
        resolve[[136, 153, 187, 221, 255]] = 8
        resolve[68]  = 4
        resolve[51]  = 32
        resolve[204] = 128
        
        #-----------------------
        # Resolve corner cases
        #-----------------------
        resolve[[1, 5, 69]]   = 1
        resolve[[4, 20, 21]]  = 4
        resolve[[16, 80, 84]] = 16
        resolve[[64, 65, 81]] = 64
        
        w = np.uint8([2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23])              #(12)
        w = np.concatenate((w,[27,31,35,38,39,43,47,54,55,66,67,70]))           #(12)
        w = np.concatenate((w,[71,74,75,78,79,82,83,86,87,91,95,103]))          #(12)
        w = np.concatenate((w,[107,111,134,135,138,139,142,143,150,151,155]))   #(11)
        w = np.concatenate((w,[159,166,167,171,175,183,191,206,207,223]))       #(10)
        
        #----------------------------------------
        # Resolve rotationally equivalent cases
        #----------------------------------------
        resolve[w] = 2
        resolve[self.resolve_array_cycle(w, 2)] = 8
        resolve[self.resolve_array_cycle(w, 4)] = 32
        resolve[self.resolve_array_cycle(w, 6)] = 128

        self.resolve = resolve

        #--------------
        # For testing
        #--------------
        if (REPORT):
            print('resolve =')
            print(np.reshape(resolve, (16,16)))
            print('sum(resolve) =', np.sum(resolve))
            print(' ')
        
    #   get_resolve_array()
    #---------------------------------------------------------------------
    def update_parent_ID_grid(self, SILENT=True):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: update_parent_ID_grid() not implemented.')

    #   update_parent_ID_grid()
    #---------------------------------------------------------------------
    def update_parent_IDs(self):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: update_parent_IDs() not implemented.')
        
    #   update_parent_IDs()
    #---------------------------------------------------------------------
    def update_flow_from_IDs(self):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: update_flow_from_IDs() not implemented.')
            
    #   update_flow_from_IDs()
    #---------------------------------------------------------------------
    def update_flow_to_IDs(self):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: update_flow_to_IDs() not implemented.')

    #   update_flow_to_IDs()
    #-------------------------------------------------------------------
    def update_noflow_IDs(self):

        #-------------------------------------------------------
        # (1/25/12) d8_global.py has an implementation of this,
        # but erode_d8_global.py doesn't use it and it affects
        # performance.  The update() method no longer calls
        # this method.
        #-------------------------------------------------------
        print('ERROR: update_nowflow_IDs() not implemented.')

    #   update_noflow_IDs()
    #-------------------------------------------------------------------
    def read_flow_grid(self, SILENT=True):

        #----------------------------------------------------
        # Read a grid of D8 flow codes, same size as DEM.
        #----------------------------------------------------
        if not(SILENT):
            print('Reading D8 flow grid...')
        code_file = (self.directory +
                     self.data_prefix + '_flow.rtg')
        self.d8_grid = rtg_files.read_grid(code_file, self.rti,
                                           RTG_type='BYTE')

    #   read_flow_grid()
    #-------------------------------------------------------------------
    def update_flow_grid(self, DEM=None, SILENT=True, REPORT=False):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: update_flow_grid() not implemented.')
 
    #   update_flow_grid()
    #-------------------------------------------------------------------    
    def start_new_d8_codes(self, DEM=None,
                           z=None, IDs=None,
                           SILENT=True, REPORT=False):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: start_new_d8_codes() not implemented.')

    #   start_new_d8_codes()
    #-------------------------------------------------------------------    
    def break_flow_grid_ties(self, SILENT=True):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: break_flow_grid_ties() not implemented.')         

    #   break_flow_grid_ties()
    #-------------------------------------------------------------------    
    def link_flats(self, SILENT=True):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: link_flats() not implemented.')
        
    #   link_flats()    
    #-------------------------------------------------------------------    
    def update_flow_width_grid(self, DOUBLE=False, METHOD2=False,
                               SILENT=True, REPORT=False):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: update_flow_width_grid() not implemented.')

    #   update_flow_width_grid()
    #-------------------------------------------------------------------
    def update_flow_length_grid(self, DOUBLE=False,
                                SILENT=True, REPORT=False):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: update_flow_length_grid() not implemented.')

    #   update_flow_length_grid()
    #-------------------------------------------------------------------
    def update_area_grid(self, SILENT=True, REPORT=False):

        #----------------------------------------------
        # (1/25/12) d8_global.py and d8_local.py each
        # have their own implementation.
        #----------------------------------------------
        print('ERROR: update_area_grid() not implemented.')
                
    #   update_area_grid()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------  
    def change_extension_to_rtg(self, filename):

        p = filename.find('.')
        if (p == -1):
            return filename
        
        extension = filename[p:]
        
        if (extension != '.rtg'):
            filename = filename[0:p] + '.rtg'

        return filename
    
    #   change_extension_to_rtg()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        pass

    #   open_input_files()
    #------------------------------------------------------------------- 
    def read_input_files(self):

        pass

    #   read_input_files()
    #-------------------------------------------------------------------  
    def close_input_files(self):

        # This is only needed because called by finalize().
        pass

    #   close_input_files()
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        RTS_FILES = True
        # RTS_FILES = False   ###########################
        if (RTS_FILES):
            self.code_gs_file = self.change_extension_to_rtg( self.code_gs_file )
            self.area_gs_file = self.change_extension_to_rtg( self.area_gs_file )            
            self.dw_gs_file   = self.change_extension_to_rtg( self.dw_gs_file )
            self.ds_gs_file   = self.change_extension_to_rtg( self.ds_gs_file )
            
        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.code_gs_file = (self.out_directory + self.code_gs_file)
        self.area_gs_file = (self.out_directory + self.area_gs_file)
        self.dw_gs_file   = (self.out_directory + self.dw_gs_file) 
        self.ds_gs_file   = (self.out_directory + self.ds_gs_file) 
        #-------------------------------------------------------------
        self.code_ts_file = (self.out_directory + self.code_ts_file)
        self.area_ts_file = (self.out_directory + self.area_ts_file) 
        self.dw_ts_file   = (self.out_directory + self.dw_ts_file) 
        self.ds_ts_file   = (self.out_directory + self.ds_ts_file) 

    #   update_outfile_names()
    #-------------------------------------------------------------------  
    def open_output_files(self):

        #-----------------------------------------------------------
        # (1/26/12) This component uses var_names "area" and 
        # "code" where Erode components use "A" and perhaps "D8".
        # Make sure that var_name passed to open_new_gs_file(),
        # open_new_ts_file(), add_grid(), etc.
        # matches the substring in: "self.*_gs_file" or else
        # open_new_gs_file() method in model_output.py won't work.
        #-----------------------------------------------------------
        # open_new_gs_file() has a "dtype" keyword that defaults
        # to "float32".  Flow codes have dtype = "uint8".
        #-----------------------------------------------------------
        model_output.check_netcdf()    # (test import and info message)
        self.update_outfile_names()

        #--------------------------------------
        # Open new files to write grid stacks
        #--------------------------------------
        if (self.SAVE_CODE_GRIDS):
            model_output.open_new_gs_file( self, self.code_gs_file, self.rti,
                                           dtype='uint8',
                                           var_name='code',
                                           long_name='D8 flow direction codes',
                                           units_name='none')

        if (self.SAVE_AREA_GRIDS):
            model_output.open_new_gs_file( self, self.area_gs_file, self.rti,
                                           var_name='area',
                                           long_name='D8 contributing area',
                                           units_name='km^2')

        if (self.SAVE_DS_GRIDS):
            model_output.open_new_gs_file( self, self.ds_gs_file, self.rti,
                                           var_name='ds',
                                           long_name='flow length in D8 direction',
                                           units_name='m')

        if (self.SAVE_DW_GRIDS):
            model_output.open_new_gs_file( self, self.dw_gs_file, self.rti,
                                           var_name='dw',
                                           long_name='flow width orthogonal to D8 direction',
                                           units_name='m')


        #--------------------------------------
        # Open new files to write time series
        #--------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_CODE_PIXELS):
            model_output.open_new_ts_file( self, self.code_ts_file, IDs,
                                           dtype='uint8',
                                           var_name='code',
                                           long_name='D8 flow direction codes',
                                           units_name='none')

        if (self.SAVE_AREA_PIXELS):
            model_output.open_new_ts_file( self, self.area_ts_file, IDs,
                                           var_name='area',
                                           long_name='D8 contributing area',
                                           units_name='km^2')

        if (self.SAVE_DS_PIXELS):
            model_output.open_new_ts_file( self, self.ds_ts_file, IDs,
                                           var_name='ds',
                                           long_name='flow length in D8 direction',
                                           units_name='m')

        if (self.SAVE_DW_PIXELS):
            model_output.open_new_ts_file( self, self.dw_ts_file, IDs,
                                           var_name='dw',
                                           long_name='flow width orthogonal to D8 direction',
                                           units_name='m')

        #-------------------------------------
        # Save FLOAT version of original DEM
        # as the rawDEM for the new DEM
        #-------------------------------------
##        if (self.rti.SWAP_ENDIAN):
##            np.array(np.float32(self.z0), copy=0).byteswap(True)
##        new_rawDEM_unit = open( self.new_rawDEM_file, 'wb' )
##        np.float32(self.z0).tofile( new_rawDEM_unit )
##        new_rawDEM_unit.close()

    #   open_output_files()
    #-------------------------------------------------------------------  
    def write_output_files(self, time_seconds=None):

        #---------------------------------------------------------
        # Notes:  This function was write_output_filestten to use only model
        #         time (maybe from a caller) in seconds, and
        #         the save_grid_dt and save_pixels_dt parameters
        #         read from CFG file
        #
        #         read_cfg_file() makes sure that all of
        #         the "save_dts" are larger than or equal to the
        #         process dt.
        #---------------------------------------------------------

        #-----------------------------------------
        # Allows time to be passed from a caller
        #-----------------------------------------
        if (time_seconds is None):
            time_seconds = self.time_sec
        model_time = int(time_seconds)

        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
        if (model_time % int(self.save_grid_dt) == 0):
            self.save_grids()
        if (model_time % int(self.save_pixels_dt) == 0):
            self.save_pixel_values()

    #   write_output_files()
    #-------------------------------------------------------------------  
    def close_output_files(self):

        if (self.SAVE_CODE_GRIDS): model_output.close_gs_file( self, 'code')
        if (self.SAVE_AREA_GRIDS): model_output.close_gs_file( self, 'area')
        if (self.SAVE_DS_GRIDS):   model_output.close_gs_file( self, 'ds')
        if (self.SAVE_DW_GRIDS):   model_output.close_gs_file( self, 'dw')
        #---------------------------------------------------------------------
        if (self.SAVE_CODE_PIXELS): model_output.close_ts_file( self, 'code')
        if (self.SAVE_AREA_PIXELS): model_output.close_ts_file( self, 'area')
        if (self.SAVE_DS_PIXELS):   model_output.close_ts_file( self, 'ds')
        if (self.SAVE_DW_PIXELS):   model_output.close_ts_file( self, 'dw')

    #   close_output_files()   
    #-------------------------------------------------------------------
    def save_grids(self, SILENT=True, REPORT=False):

        if not(SILENT):
            print('Saving D8 grids specified in CFG file...')

        #--------------------------------------------------------------
        # Note: filenames are set in the CFG file.  The out_directory
        #       is prepended when self.open_output_files() calls
        #       self.update_outfile_names().  The latter method
        #       should never be called from this method.
        #--------------------------------------------------------------
        
        #------------------------------
        # Save D8 flow direction grid
        #------------------------------
        if (self.SAVE_CODE_GRIDS):
            model_output.add_grid( self, self.d8_grid, 'code', self.time )
            #--------------------
            # For one-time use.
            #--------------------
            # rtg_files.write_grid( self.d8_grid, self.code_gs_file, self.rti,
            #                       RTG_type='BYTE')

        #---------------------------------
        # Save D8 contributing area grid
        #---------------------------------
        if (self.SAVE_AREA_GRIDS):
            model_output.add_grid( self, self.A, 'area', self.time )
            #--------------------
            # For one-time use.
            #--------------------
            # rtg_files.write_grid( self.A, self.area_gs_file, self.rti)
 
        #------------------------------
        # Save D8 flow width grid, dw
        #------------------------------
        if (self.SAVE_DW_GRIDS):
            model_output.add_grid( self, self.dw, 'dw', self.time )
            #--------------------
            # For one-time use.
            #--------------------
            # rtg_files.write_grid( self.dw, self.dw_gs_file, self.rti)

        #-------------------------------
        # Save D8 flow length grid, ds
        #-------------------------------
        if (self.SAVE_DS_GRIDS):
            model_output.add_grid( self, self.ds, 'ds', self.time )
            #--------------------
            # For one-time use.
            #--------------------
            # rtg_files.write_grid( self.ds, self.ds_gs_file, self.rti)

    #   save_grids()
    #-------------------------------------------------------------------
    def save_pixel_values(self):
         # save_time_series_data(self): 

        IDs  = self.outlet_IDs
        time = self.time         #####   years ??

        #-------------
        # New method
        #-------------
        if (self.SAVE_CODE_PIXELS):
            model_output.add_values_at_IDs( self, time, self.d8_grid, 'code', IDs )

        if (self.SAVE_AREA_PIXELS):
            model_output.add_values_at_IDs( self, time, self.A, 'area', IDs )

        if (self.SAVE_DS_PIXELS):
            model_output.add_values_at_IDs( self, time, self.ds, 'ds', IDs )

        if (self.SAVE_DW_PIXELS):
            model_output.add_values_at_IDs( self, time, self.dw, 'dw', IDs )

    #   save_pixel_values()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
##    def update_area_grid_OLD(self, SILENT=True, REPORT=False):
##
##        #------------------------------------------------------
##        # Notes: This version cannot handle periodic boundary
##        #        conditions.
##        #------------------------------------------------------
##        # Notes: Idea is to find the pixels whose area has
##        #        not yet been assigned, and to recursively
##        #        assign areas to those for which all children
##        #        have been assigned areas.
##
##        #        g = [1, 2, 4, 8, 16, 32, 64, 128]
##        #        h = [16, 32, 64, 128, 1, 2, 4, 8]
##
##        #        Notice that the first and last pixel of each
##        #        line are assigned an area of zero.
##
##        #        DIRECTION CODES:      SUBSCRIPTS:
##        #        ----------------      -----------
##        #           |64 128 1|          |6  7  0|    (i1)
##        #           |32  x  2|          |5  x  1|    (i2)
##        #           |16  8  4|          |4  3  2|    (i3)
##        #------------------------------------------------------
##        if not(SILENT):    
##            print 'Updating upstream area grid...'
##        
##        #------------------
##        # Initialize vars
##        #------------------------------------------------
##        # initialize_computed_vars() initializes self.A
##        #------------------------------------------------
##        self.A = minimum(self.A, 0)   # (reset to all zeros)
##        n_reps = int32(0)
##
##        #-----------------        
##        # Local synonyms
##        #-----------------
##        nx = self.nx
##        ny = self.ny
##        h  = self.code_opps
##        # g  = self.code_list  # (not used here)\
##        
##        #-----------------------------------------
##        # da was stored by self.read_grid_info()
##        #-----------------------------------------
##        pixel_area = (self.da  / 1e6) # [m^2 -> km^2]
##        
##        #-------------------------------------------------
##        # Use "flattened versions" of grids for indexing
##        #-------------------------------------------------
##        A    = self.A.flat
##        dirs = self.d8_grid.flat  # (these are "iterators")
##        
##        while (True):
##        
##            STILL_ACTIVE = False
##            n_reps += 1
##            
##            #-------------------------------------------------------
##            # Added test for (dirs ne 0), so that edge
##            # and nodata pixels won't be considered "READY".
##            #-------------------------------------------------------
##            # This led to an infinite loop bug because the edge
##            # pixels were "ready", then LINECHANGE=1b, then their
##            # areas were changed back to zero in caller which made
##            # them "ready" again.
##            #-------------------------------------------------------
##            unknown   = np.where(np.logical_and((A == 0), (dirs != 0)))           
##            IDs       = unknown[0]
##            n_unknown = np.size( IDs )
##            
##            if (n_unknown != 0):    
##                #-----------------------------------------
##                # Initialize all unknown pixels as READY
##                #-----------------------------------------
##                READY = np.zeros([n_unknown], dtype='UInt8') + 1
##                
##                #--------------------------------------
##                # Upstream areas of 8 neighbor pixels
##                #--------------------------------------
##                a0 = A[IDs - nx + 1]   #(upper-right)
##                a1 = A[IDs + 1]        #(right)
##                a2 = A[IDs + nx + 1]   #(lower-right)
##                a3 = A[IDs + nx]       #(bottom)
##                a4 = A[IDs + nx - 1]   #(lower-left)
##                a5 = A[IDs - 1]        #(left)
##                a6 = A[IDs - nx - 1]   #(upper-left)
##                a7 = A[IDs - nx]       #(top)
##                
##                #----------------------------------
##                # Flow codes of 8 neighbor pixels
##                #----------------------------------
##                d0 = dirs[IDs - nx + 1]   #(upper-right)
##                d1 = dirs[IDs + 1]        #(right)
##                d2 = dirs[IDs + nx + 1]   #(lower-right)
##                d3 = dirs[IDs + nx]       #(bottom)
##                d4 = dirs[IDs + nx - 1]   #(lower-left)
##                d5 = dirs[IDs - 1]        #(left)
##                d6 = dirs[IDs - nx - 1]   #(upper-left)
##                d7 = dirs[IDs - nx]       #(top)
##                
##                #----------------------------------
##                # A pixel is not READY if any of
##                # it's children have unknown area
##                #----------------------------------
##                # Pixels w/ no children are READY
##                # unless they are "edge" pixels
##                #----------------------------------
##                w0 = np.where(np.logical_and((d0 == h[0]), (a0 <= 0)))
##                nw0 = np.size(w0[0])
##                if (nw0 != 0):    
##                    READY[w0] = np.uint8(0)
##                w1 = np.where(np.logical_and((d1 == h[1]), (a1 <= 0)))
##                nw1 = np.size(w1[0])
##                if (nw1 != 0):    
##                    READY[w1] = np.uint8(0)
##                w2 = np.where(np.logical_and((d2 == h[2]), (a2 <= 0)))
##                nw2 = np.size(w2[0])
##                if (nw2 != 0):    
##                    READY[w2] = np.uint8(0)
##                w3 = np.where(np.logical_and((d3 == h[3]), (a3 <= 0)))
##                nw3 = np.size(w3[0])
##                if (nw3 != 0):    
##                    READY[w3] = np.uint8(0)
##                w4 = np.where(np.logical_and((d4 == h[4]), (a4 <= 0)))
##                nw4 = np.size(w4[0])
##                if (nw4 != 0):    
##                    READY[w4] = np.uint8(0)
##                w5 = np.where(np.logical_and((d5 == h[5]), (a5 <= 0)))
##                nw5 = np.size(w5[0])
##                if (nw5 != 0):    
##                    READY[w5] = np.uint8(0)
##                w6 = np.where(np.logical_and((d6 == h[6]), (a6 <= 0)))
##                nw6 = np.size(w6[0])
##                if (nw6 != 0):    
##                    READY[w6] = np.uint8(0)
##                w7 = np.where(np.logical_and((d7 == h[7]), (a7 <= 0)))
##                nw7 = np.size(w7[0])
##                if (nw7 != 0):    
##                    READY[w7] = np.uint8(0)
##                
##                #----------------------------------------
##                # If a pixel is my child and I'm READY,
##                # then add it's value to mine.
##                #----------------------------------------
##                WR      = np.where(READY)
##                n_ready = np.size(WR[0])
##                WR      = WR[0]   ##### (Need this) #####
##                
##                if (n_ready != 0):    
##                    #---------------------------
##                    # This also assigns areas
##                    # to pixels w/ no children
##                    #---------------------------
##                    ### PASS_CHANGE  = True
##                    STILL_ACTIVE = True
##                    UWR = IDs[WR]         #####
##                    A[UWR] = pixel_area
##                    
##                    w0  = np.where(d0[WR] == h[0])
##                    nw0 = np.size(w0[0])
##                    if (nw0 != 0):    
##                        A[UWR[w0]] += a0[WR[w0]]
##                    w1  = np.where(d1[WR] == h[1])
##                    nw1 = np.size(w1[0])
##                    if (nw1 != 0):    
##                        A[UWR[w1]] += a1[WR[w1]]
##                    w2  = np.where(d2[WR] == h[2])
##                    nw2 = np.size(w2[0])
##                    if (nw2 != 0):    
##                        A[UWR[w2]] += a2[WR[w2]]
##                    w3  = np.where(d3[WR] == h[3])
##                    nw3 = np.size(w3[0])
##                    if (nw3 != 0):    
##                        A[UWR[w3]] += a3[WR[w3]]
##                    w4  = np.where(d4[WR] == h[4])
##                    nw4 = np.size(w4[0])
##                    if (nw4 != 0):    
##                        A[UWR[w4]] += a4[WR[w4]]
##                    w5  = np.where(d5[WR] == h[5])
##                    nw5 = np.size(w5[0])
##                    if (nw5 != 0):    
##                        A[UWR[w5]] += a5[WR[w5]]
##                    w6  = np.where(d6[WR] == h[6])
##                    nw6 = np.size(w6[0])
##                    if (nw6 != 0):    
##                        A[UWR[w6]] += a6[WR[w6]]
##                    w7  = np.where(d7[WR] == h[7])
##                    nw7 = np.size(w7[0])
##                    if (nw7 != 0):    
##                        A[UWR[w7]] += a7[WR[w7]]
##                
##                #------------------------
##                # Are we finished now ?
##                #------------------------
##                if (n_unknown != 0) and (n_ready != n_unknown):   
##                    UNFINISHED = True
##                else:    
##                    UNFINISHED = False
##                
##            else:    
##                UNFINISHED = False
##
##            if not(UNFINISHED) or not(STILL_ACTIVE):
##                break
##        
##        if (UNFINISHED):    
##            print 'Upstream area not defined for all pixels.'
##
##        #-----------------------------
##        # Save area grid for testing
##        #-----------------------------
##        SAVE_TEST = False
##        if (SAVE_TEST):
##            file_unit = open('00_AREA_GRID_TEST.rtg', 'wb')
##            if (self.rti.SWAP_ENDIAN):
##                grid = self.A.copy()
##                grid.byteswap(True)
##                grid.tofile(file_unit)
##            else:
##                self.A.tofile( file_unit )
##            file_unit.close()
##        
##        #------------------
##        # Optional report
##        #------------------
##        if (REPORT):
##            A_str = str(self.A.min()) + ', ' + str(self.A.max())
##            print '    min(A), max(A) = ' + A_str + ' [km^2]'
##            print '    Number of iterations = ' + str(n_reps)
##    
##    #   update_area_grid_OLD()
##    #-------------------------------------------------------------------
##    #-------------------------------------------------------------------    
##    def get_flow_width_grid(self, DOUBLE=False, METHOD2=False,
##                            SILENT=True, REPORT=False):
##
##        #-------------------------------------------------------------
##        # NOTES: This routine returns the flow widths for each
##        #        pixel in the DEM.  The extra width of flowing to a
##        #        diagonal pixel is taken into account, as well as the
##        #        lat/lon-dependence for DEMs with fixed-angle pixels
##        #        (Geographic lat/lon).
##
##        #        METHOD2 version ensures that the sum of all flow
##        #        widths around a pixel is equal to 2*(dx + dy), but
##        #        is incorrect for case of a plane and others.
##
##        #        Flow widths are zero where (flow grid eq 0).
##        
##        #        Flow widths are for entire pixel and are appropriate
##        #        for overland or subsurface flow.
##        #------------------------------------------------------------- 
##        #        Is this only used by Seepage function now ???
##        #-------------------------------------------------------------
##        # NB!    np.where returns a 2-tuple, but with "empty"
##        #        second part when applied to a 1D array.  This means
##        #        that nw = size(w[0]) still works.
##        #-------------------------------------------------------------        
##        if not(SILENT):
##            print 'Computing flow width grid...'
##
##        #-----------------------
##        # Get pixel dimensions
##        #-----------------------
##        dx, dy, dd = pixels.get_sizes_by_row(self.rti, METERS=True)
##
##        #-------------------------
##        # Double or Float type ?
##        #-------------------------
##        if (DOUBLE):    
##            dw = np.zeros([self.ny, self.nx], dtype='Float64')
##        else:    
##            dw = np.zeros([self.ny, self.nx], dtype='Float32')
##            dx = np.float32(dx)
##            dy = np.float32(dy)
##            dd = np.float32(dd)
##
##        #----------------------------------------------
##        # Initialize to default value that is used
##        # for pixels with flow code of zero.  This is
##        # done to avoid "divide by zero" errors.
##        #----------------------------------------------
##        dw += dx[0]
##        ## dw = dw + dx.min()
##        
##        for row in xrange(self.ny):
##            g = self.d8_grid[row,:]
##            
##            #----------------
##            # Diagonal flow
##            #----------------
##            wd  = np.where(np.logical_or(np.logical_or(np.logical_or((g == 1), (g == 4)), \
##                                              (g == 16)), (g == 64)))
##            nwd = np.size(wd[0])
##            if (nwd != 0):    
##                if not(METHOD2):    
##                    dw[row, wd] = dd[row]
##                else:    
##                    dw[row, wd] = (dx[row] + dy[row]) / 4
##            
##            #---------------------
##            # East and west flow
##            #---------------------
##            wh  = np.where(np.logical_or((g == 2), (g == 32)))
##            nwh = np.size(wh[0])
##            if (nwh != 0):    
##                dw[row, wh] = dy[row]
##                if (METHOD2):    
##                    dw[row, wh] = dw[row,wh] / 2
##            
##            #-----------------------
##            # North and south flow
##            #-----------------------
##            wv  = np.where(np.logical_or((g == 8), (g == 128)))
##            nwv = np.size(wv[0])
##            if (nwv != 0):    
##                dw[row, wv] = dx[row]
##                if (METHOD2):    
##                    dw[row, wv] = dw[row, wv] / 2
##
##        #------------------
##        # Optional report
##        #------------------
##        if (REPORT):
##            print '    min(dw) = ' + str(dw.min()) + '  [m]'
##            print '    max(dw) = ' + str(dw.max()) + '  [m]'
##
##        self.dw = dw
##
##    #   get_flow_width_grid()
##    #-------------------------------------------------------------------
##    def get_flow_length_grid(self, DOUBLE=False,
##                             SILENT=True, REPORT=False):
##
##        #-------------------------------------------------------------
##        # NOTES: This routine returns the flow lengths for each
##        #        pixel in the DEM.  The extra length of flowing to a
##        #        diagonal pixel is taken into account, as well as the
##        #        latitude-dependence for DEMs with fixed-angle pixels
##        #        (Geographic lat/lon).
##
##        #        The only difference between this and the Flow_Widths
##        #        function is that roles of dx and dy are switched.
##
##        #        Flow lengths are set to dx[0] for the pixels where
##        #       (flow grid eq 0), such as on the edges of the DEM.
##        #-------------------------------------------------------------
##        if not(SILENT):
##            print 'Computing flow length grid...'
##
##        #-----------------------
##        # Get pixel dimensions
##        #-----------------------
##        dx, dy, dd = pixels.get_sizes_by_row(self.rti, METERS=True)
##            
##        #-------------------------
##        # Double or Float type ?
##        #-------------------------
##        if (DOUBLE):    
##            ds = np.zeros([self.ny, self.nx], dtype='Float64')
##        else:    
##            ds = np.zeros([self.ny, self.nx], dtype='Float32')
##            dx = np.float32(dx)
##            dy = np.float32(dy)
##            dd = np.float32(dd)
##        
##        #----------------------------------------------
##        # Initialize to default value that is used
##        # for pixels with flow code of zero.  This is
##        # done to avoid "divide by zero" errors.
##        #----------------------------------------------
##        ds += dx[0]
##        ## ds += dx.min()
##        
##        for row in xrange(self.ny):
##            g = self.d8_grid[row,:]
##            
##            #----------------
##            # Diagonal flow
##            #----------------
##            wd  = np.where(np.logical_or(np.logical_or(np.logical_or((g == 1), (g == 4)), \
##                                              (g == 16)), (g == 64)))
##            nwd = np.size(wd[0])
##            if (nwd != 0):    
##                ds[row, wd] = dd[row]
##            
##            #---------------------
##            # East and west flow
##            #---------------------
##            wh  = np.where(np.logical_or((g == 2), (g == 32)))
##            nwh = np.size(wh[0])
##            if (nwh != 0):    
##                ds[row, wh] = dx[row]
##            
##            #-----------------------
##            # North and south flow
##            #-----------------------
##            wv  = np.where(np.logical_or((g == 8), (g == 128)))
##            nwv = np.size(wv[0])
##            if (nwv != 0):    
##                ds[row, wv] = dy[row]
##
##        #------------------
##        # Optional report
##        #------------------
##        if (REPORT):
##            print '    min(ds) = ' + str(ds.min()) + '  [m]'
##            print '    max(ds) = ' + str(ds.max()) + '  [m]'
##
##        self.ds = ds
##
##    #   get_flow_length_grid()
##    #-------------------------------------------------------------------

