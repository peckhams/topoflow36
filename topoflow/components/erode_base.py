
## Copyright (c) 2001-2013, Scott D. Peckham
##
## Jan 2013.  Fixed bug in create_initial_DEM() (PLANE option)
## Jan 2009  (converted from IDL)
## September, October, November 2009
## February 2010 (local time-stepping)
## August 2010 (updated to use new model_output, conventions, etc.)
## Sep to Nov 2011 (miscellaneous)
## Jan 16, 2012 (Fixed bug: self.time_step_type vs. self.FIXED_STEPS.)

## NOTE: D8 components no longer compute "noflow_IDs" by default
##       even though d8_global.py has a method to do so.  It is
##       only used here for the dt method that uses update_min_dz_up().

## See: FILL_PITS_IN_Z0   (in initialize_DEM())  ##################

#-----------------------------------------------------------------------
#
#  class erosion_component
#
#      set_constants()
#--------------------------
#      initialize()
#      update()           # (in erode_d8_local.py & erode_d8_global.py)
#      finalize()
#      set_computed_input_vars()
#----------------------------------
#      initialize_DEM()
#      create_initial_DEM()
#      initialize_boundary_conditions()
#      initialize_computed_vars()
#----------------------------------
#      update_R()
#      update_R_integral()
#      update_U()
#      update_U_integral()
#      update_base_level()
#      update_DEM_edge_values()     ###
#      update_d8_vars()
#-----------------------------------------------
# See erode_d8_local.py and erode_d8_global.py
#-----------------------------------------------
#####      update()
#####      update_slope_grid()
#####      update_Q_grid()
#####      update_Qs_grid()
#####      update_dz_dt()
#####      update_DEM()
#--------------------------------
#      fill_pits_in_DEM()           ###
#      update_DEM_min_and_max()
#--------------------------------
#      update_dt_grid()
#      update_dt_grid_method1()     ## (3/9/10, based on original)
#      update_min_dz_up_grid()      ## (2/19/10)
#      update_dt_grid_method2()     ## (2/9/10, based on min_dz_up)
#      update_dt_grid_local1()      ## (2/19/10)
#      update_n_grid()              ## (2/19/10)
#      update_n_grid2()             ## (2/19/10)  COMMENTED OUT
#--------------------------------
#      print_mins_and_maxes()       ## (3/12/10)
#      check_stability()
#      check_finished()
#      check_steady_state()         ### (not written yet)
#--------------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#--------------------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#      print_time_and_value()       ## (override one in BMI_base.py)
#
#---------------------------
#  THESE ARE NOT USED YET
#---------------------------
#  Get_Timestep()
#  Stable_Timestep()
#
#-----------------------------------------------------------------------

import numpy as np
import os.path
import time
#----------------------
# Could also do this.
#----------------------
# from numpy import where, logical_and, logical_or

from topoflow.components import d8_global

from topoflow.utils import BMI_base
from topoflow.utils import fill_pits
from topoflow.utils import midpoints
from topoflow.utils import model_input
from topoflow.utils import model_output
from topoflow.utils import pixels
from topoflow.utils import rtg_files
from topoflow.utils import rti_files

#-------------------------------------------
# For use outside of the TopoFlow package.
#-------------------------------------------
##import d8_global
###------------------
##import BMI_base
##import fill_pits  #####
##import midpoints  #####
##import model_input
##import model_output
##import pixels
##import rtg_files
##import rti_files

## from model_output import *   ## ELIMINATE THIS

#-----------------------------------------------------------------------
class erosion_component( BMI_base.BMI_component ):

    #-------------------------------------------------------------------
    def set_constants(self):

        #---------------------------------------------
        # (2/5/12) Set default stopping method:
        # (0=n_steps, 1=stop_time, 2=steady-state)
        # This is now done in initialize_time_vars()
        # in CSDMS_base.py. (2/14/12)
        #---------------------------------------------
        ## self.stop_code = 0
        ## self.stop_time = 1000.0  # [years]

        #------------------------
        # Define some constants
        #------------------------
        self.secs_per_year = np.float64(31536000)
        self.mm_per_m      = np.float64(1000)
        self.dz_tolerance  = np.float64(1e-6)
        ## self.dz_tolerance  = np.float64(1e-5)
        ## self.dz_tolerance  = np.float64(1e-4)
        ## self.dz_tolerance  = np.float64(1e-3)
  
        #--------------------------------
        # Set constants related to time
        #--------------------------------
        self.dt_too_small = np.float64(1e-2)
        # dt_too_big = 5d  # (should there be a max?)
        self.dt_limit = 1e20  # [years]
        self.dt = 1.0  # [year]

        ##########################################
        # CFL_factor is now set in the CFG file.
        ##########################################
        ## self.CFL_factor = np.float64(0.2)

        #--------------------------------------
        # Removed this option from CFG files.
        # See erode_d8_global.update().
        #--------------------------------------
        self.FILL_PITS = False

        #-------------------------------------------
        # (11/20/11) Need to avoid frequent string
        # comparisons for best performance.
        #-------------------------------------------
        self.time_step_type = self.get_attribute( 'time_step_type' )
        self.FIXED_STEPS = (self.time_step_type.lower() == 'fixed')
        self.DRIVER      = (self.mode == 'driver')
        
    #   set_constants()
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        #--------------------------------------------------------
        # Note: CSDMS_base.run_model() changes to cfg_directory.
        #       And now get_neighbors_test() does this also.
        #--------------------------------------------------------       
        if not(SILENT):
            print('Erosion component: Initializing...')

        self.status = 'initializing'  # (OpenMI 2.0 convention)
        self.mode   = mode
        self.cfg_file = cfg_file

        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()
        self.initialize_config_vars()

##        print '### In erode_base.initialize():'
##        print '### out_directory =', self.out_directory
##        print ' '
        
        #-------------------------------------------------------
        # Makes more sense to generate RTI structure,
        # instead of reading a separate RTI file. (11/7/11)
        # Note: basins.py needs RTI file for read_grid_info().
        # Note: read_grid_info() is in BMI_base.py.
        #-------------------------------------------------------
        # self.read_grid_info() 
        rti = rti_files.make_info( grid_file=self.z0_file,
                                   ncols=self.nx, nrows=self.ny,
                                   xres=self.dx, yres=self.dy )
        self.rti = rti
        self.da = pixels.get_da( rti )
        RTI_file = self.out_directory + self.case_prefix + '.rti'
        rti_files.write_info( RTI_file, rti )
        #--------------------------------------------------
        # Note: BMI_base.initialize_basin_vars() calls
        # basins.initialize() which then calls
        # BMI_base.read_grid_info().  It looks for an
        # RTI file in input directory first, then in the
        # output directory. (11/5/13)
        #--------------------------------------------------
        self.initialize_basin_vars()
        
        #-------------------------------------------
        # This must come before "Disabled" test ??
        #-------------------------------------------
        self.initialize_time_vars( units='years' )
          
        #---------------------------------------------
        # Create an RTI file with grid info which
        # can be used by this and other components ?
        #---------------------------------------------
        if (self.make_z0_method != 'READ_FILE'):
            out_prefix   = self.out_directory + self.case_prefix
            DEM_file     = (out_prefix + '_2D-z0.bin')
            RTI_file     = (out_prefix + '.rti')
            self.z0_file = DEM_file
            self.rti = rti_files.make_info( DEM_file, self.nx,
                                            self.ny, self.dx, self.dy )
            rti_files.write_info(RTI_file, self.rti)

        #----------------------------------
        # Has component been turned off ?
        #----------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print('Erosion component: Disabled in CFG file.')
            self.SAVE_Z_GRIDS  = False    # (It is True by default.)
            self.SAVE_Z_PIXELS = False    # (It is True by default.)
            self.DONE = True
            self.status = 'initialized'  # (OpenMI 2.0 convention) 
            return
                                                  
        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        # Can't move read_input_files() to start of
        # update(), since initial values needed here.
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()

        #-----------------------
        # Initialize variables
        #---------------------------------------------------------
        # Do we need initialize_DEM() after initialize_d8_vars()
        # for erode_d8_global.py ??  #### CHECK THIS ####
        #---------------------------------------------------------
##        self.initialize_d8_vars()  # (depend on D8 flow grid)
##        self.initialize_DEM()
##        self.initialize_boundary_conditions()
##        self.initialize_computed_vars()
        
        #-----------------------
        # Initialize variables
        #----------------------------------------------------
        # Need initialize_DEM() before initialize_d8_vars()
        # for erode_d8_local.py.
        #----------------------------------------------------
        self.initialize_DEM()      # (must come before d8_vars)
        self.initialize_d8_vars()  # (depend on D8 flow grid)
        self.initialize_boundary_conditions()
        self.initialize_computed_vars()

        #--------------------------------------------
        # New safety precaution. (2/14/12)
        # Make sure rtg_files.write_grid() function
        # does not change data type fo 'float32'.
        #--------------------------------------------
        print('Data type of initial DEM is:', str(self.DEM.dtype))
        print(' ')

        #--------------------------------------------------
        # Make sure self.Q_ts_file is not NULL (12/22/05) 
        # This is only output file that is set by default
        # and is still NULL if user hasn't opened the
        # output var dialog for the channel process.
        #--------------------------------------------------
        if (self.SAVE_Z_PIXELS and (self.z_ts_file == '')):    
            self.z_ts_file = (self.case_prefix + '_0D-z.txt')
        
        self.open_output_files()
        self.status = 'initialized'  # (OpenMI 2.0 convention) 
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, time=None, SILENT=None, REPORT=None):

        #------------------------------------------------------
        # Note: erode_d8_local.py and erode_d8_global.py each
        #       have their own version of this function.
        #------------------------------------------------------
        pass
    
    #   update()
    #-------------------------------------------------------------------
##    def update(self, time=None, SILENT=None, REPORT=None):        
##        if (SILENT == None): SILENT=self.SILENT
##        if (REPORT == None): REPORT=self.REPORT
##
##        #### if not(SILENT) and (self.time_index == 0):
##        if (self.time_index == 0):
##            print 'Erosion component: Processing...'
##        self.status = 'updating'  # (OpenMI 2.0 convention)
##
##        #---------------------------------
##        # Print dz_max to track progress
##        #---------------------------------
##        if (self.mode == 'main'):
##            self.print_time_and_value(self.dz_max, 'dz_max', '[m]',
##                                      interval=5.0, PRINT_INDEX=True)
##            
##        #--------------------------------------
##        # Update values from other components
##        #--------------------------------------
##        self.update_R()
##        self.update_R_integral()
##        self.update_U()
##        self.update_U_integral()
##        
##        #-------------------------
##        # Update computed values
##        #-------------------------
##        self.update_base_level()
##        self.update_DEM_edge_values()   ######
##        if (self.FILL_PITS):
##            self.fill_pits_in_DEM(SILENT=SILENT)    ############
##
##        #--------------------------------------------
##        # Update the D8 flow grid and all vars that
##        # depend on it, including D8 area grid.
##        #--------------------------------------------
##        self.update_d8_vars(SILENT=SILENT, REPORT=REPORT)  #########
##        self.update_slope_grid(SILENT=SILENT, REPORT=REPORT)
##        self.update_Q_grid(SILENT=SILENT, REPORT=REPORT)
##        self.update_Qs_grid(SILENT=SILENT, REPORT=REPORT)
##        self.update_dz_dt_grid(SILENT=SILENT, REPORT=REPORT)
##        
##        self.update_DEM(SILENT=SILENT, REPORT=REPORT)
##        self.update_DEM_min_and_max()
##
##        ########################################
##        # CAN THE UPDATED DEM HAVE PITS ??
##        # IF NOT, DON'T CALL FILL_PITS.
##        ########################################
##        
##        #------------------------
##        # Check computed values
##        #------------------------
##        OK = self.check_stability()
##
##        #-------------------------------------------
##        # Read from files as needed to update vars 
##        #-----------------------------------------------------
##        # NB! This is currently not needed for the "erosion
##        # process" because values don't change over time and
##        # read_input_files() is called by initialize().
##        #-----------------------------------------------------
##        # if (self.time_index > 0):
##        #     self.read_input_files()
##
##        #----------------------------------------------
##        # Write user-specified data to output files ?
##        #----------------------------------------------
##        self.write_output_files( time )
##        if (OK):
##            self.status = 'updated'  # (OpenMI 2.0 convention)
##        else:
##            self.status = 'failed'
##            self.DONE   = True
##
##        #------------------------
##        # Update internal clock
##        #------------------------
##        self.update_time()
##        ## print 'time_index =', self.time_index
##        
##        #-------------------------------------------
##        # Check for steady-state condition instead
##        #-------------------------------------------
##        self.check_finished()   ######################
##        
##    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalizing'  # (OpenMI)   
        self.close_input_files()   ##  TopoFlow input "data streams"
        self.close_output_files()

        #-----------------------------
        # Save final DEM to DEM_file        # (write a separate function ###########)
        #-----------------------------
##        DEM_unit = open(self.final_DEM_file, 'wb')
##        if (self.rti.SWAP_ENDIAN):
##            final_DEM = self.DEM.byteswap(True)
##            final_DEM.tofile(DEM_unit)
##        else:
##            self.DEM.tofile( DEM_unit )
##        DEM_unit.close()
    
        #----------------------
        # Print final message
        #----------------------
        print(' ')
        print('min(dt), max(dt) =', self.dt_grid.min(), self.dt_grid.max())
        print(' ')
        print('min(z), max(z)   =', self.DEM.min(), self.DEM.max())
        print(' ')
        print('dz_max_vec[0]    =', self.dz_max_vec[0])
        print('dz_max_vec[nt-1] =', self.dz_max_vec[-1])
        print('dz_max_vec.min() =', self.dz_max_vec.min())
        print('dz_max_vec.max() =', self.dz_max_vec.max())
        print(' ')
        self.print_final_report( comp_name='Erode 3.1 (11/15/11)' )


        self.status = 'finalized'  # (OpenMI)
        
    #   finalize()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        self.U_mpyr       = self.U   / self.mm_per_m
        self.BLR_mpyr     = self.BLR / self.mm_per_m 
        #--------------------------------------------------------------------
        self.FLAT         = (self.make_z0_method.upper() == 'FLAT')
        self.PLANE        = (self.make_z0_method.upper() == 'PLANE')
        self.CORNER_PLANE = (self.make_z0_method.upper() == 'CORNER_PLANE')
        self.READ_FILE    = (self.make_z0_method.upper() == 'READ_FILE')
        #--------------------------------------------------------------------
        self.GAUSSIAN     = (self.noise_method.upper() == 'GAUSSIAN')
        self.MIDPOINTS    = (self.noise_method.upper() == 'MIDPOINTS')
        self.NO_NOISE     = (self.noise_method.upper() == 'NO_NOISE')
        #--------------------------------------------------------------------        
        self.BOTTOM       = (self.BC_method.upper() == 'BOTTOM')
        self.RIGHT        = (self.BC_method.upper() == 'RIGHT')
        self.CORNER       = (self.BC_method.upper() == 'CORNER')
        self.FOUR_SIDES   = (self.BC_method.upper() == 'FOUR_SIDES')

        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt   = np.maximum(self.save_grid_dt,   self.dt)
        self.save_pixels_dt = np.maximum(self.save_pixels_dt, self.dt)
        
    #   set_computed_input_vars()
    #-------------------------------------------------------------------
##    def initialize_d8_vars(self):
##
##        #------------------------------------------------------
##        # Note: erode_d8_local.py and erode_d8_global.py each
##        #       have their own version of this function.
##        #------------------------------------------------------
##        pass
##    
##    #   initialize_d8_vars()        
    #-------------------------------------------------------------
    def initialize_DEM(self, SILENT=False):

        #-------------------------------------------------------
        # Notes: This function initializes the DEM either from
        #        self.z0_file or using create_initial_DEM().
        #-------------------------------------------------------
        ### if (self.z0_file != ''):
        if (self.make_z0_method == 'READ_FILE'):
            #---------------------------------------
            # Read inital elevation grid from file
            #---------------------------------------
            DEM_unit  = open(self.z0_file, 'rb')
            file_size = os.path.getsize(DEM_unit.name)
            dbl_size  = self.nx * self.ny * np.int32(8)
            flt_size  = self.nx * self.ny * np.int32(4)
            int_size  = self.nx * self.ny * np.int32(2)
            if (file_size == dbl_size):
                RTG_type = 'DOUBLE'
            elif (file_size == flt_size):
                RTG_type = 'FLOAT'
            elif (file_size == int_size):
                RTG_type = 'INTEGER'
            else:
                print('ERROR in initialize_DEM().')
                print('   Cannot determine DEM data type.')
                return
            
            self.DEM = rtg_files.read_grid( self.z0_file, self.rti,
                                            RTG_type=RTG_type ) 
            if (RTG_type != 'FLOAT'):
                self.DEM = np.float32(self.DEM)
            if not(SILENT):
                print('Read DEM from:', self.z0_file)
        else:    
            #--------------------------------
            # Create initial elevation grid
            #--------------------------------
            self.create_initial_DEM(SILENT=SILENT)

            #-------------------------------------
            # See self.update_DEM_edge_values().
            #----------------------------------------------------------
            # 01/13/06. Otherwise we end up with very large values on
            # the indicated boundary that prevent good use of color
            #----------------------------------------------------------
            # self.update_DEM_edge_values()

            #------------------------------------
            # Save the initial DEM to a file ??
            #------------------------------------
            #### Check for overwrite here !!!  ##############################
##            DEM_unit = open(self.z0_file, 'wb')
##            if (self.rti.SWAP_ENDIAN):   ## (rti undefined so far)
##                self.DEM.byteswap(True)
##            self.DEM.tofile(DEM_unit)
##            DEM_unit.close()

        #---------------------------------------------------
        # Option to fill pits in the initial DEM (2/23/10)
        #---------------------------------------------------
        FILL_PITS_IN_Z0 = False
        ## FILL_PITS_IN_Z0 = True
        if (FILL_PITS_IN_Z0):
            print('Filling pits in initial DEM...')
            self.fill_pits_in_DEM(SILENT=True)
            #-------------------------------------------
            # Save new DEM to a file.  This overwrites
            # the one saved by create_initial_DEM()
            #------------------------------------------- 
            rtg_files.write_grid( self.DEM, self.z0_file, self.rti)

        self.DEM_min = np.nanmin( self.DEM )
        self.DEM_max = np.nanmax( self.DEM )
        if not(SILENT):
            print('Initial (z_min, z_max) =', self.DEM_min, self.DEM_max)
            print(' ')
            
    #   initialize_DEM()
    #-------------------------------------------------------------
    def create_initial_DEM(self, SILENT=False):

        #------------------------------------------------------
        # Notes: This routine allows a "z0_method" and a
        #        "noise_method" to be combined to generate
        #        an initial surface (z0) with optional noise.
        #
        #        z0_methods:    FLAT, PLANE, CORNER_PLANE
        #        noise_methods: GAUSSIAN, MIDPOINTS

        #        Noise grid is generated first, which is
        #        more efficient in the FLAT case.

        #        "seed" is read from CFG file and should be
        #        a 4- or 5-digit integer (e.g. 36421).
        #        If the same seed is used, the same sequence
        #        of random numbers is generated, which allows
        #        for reproducible results and comparisons.
        #------------------------------------------------------
        nx = self.nx  # (local synonyms)
        ny = self.ny

        #---------------------------------------
        # Uncorrelated Gaussian random numbers
        # which essentially models white noise
        #---------------------------------------
        # When sigma or factor = 1, then range
        # is pretty much -3 to 3.
        #---------------------------------------
        if (self.GAUSSIAN):    
            np.random.seed( self.seed )
            self.DEM = np.random.normal(loc=0.0, scale=1.0, size=(ny, nx))
               #(mean = 0.0, stddev = 1.0)
            #-----------------------------------
            # factor = (1 / float64(3))               #(-3,3) -> (-1,1) 
            # factor = factor * float64(300)
            factor = self.noise_scale
            
            #-----------------------------------
            # Slope should dominate over noise
            #-----------------------------------
            #*** factor = factor * (slope / 2d)
            self.DEM = factor * self.DEM
            #------------------------------
            # if (PLANE OR CORNER_PLANE) then mean=0.0 else mean=2.0
            # stddev = 4.0  ;(1.0)
            # DEM = (stddev * DEM) + mean
            if not(SILENT):
                print('Created initial DEM by GAUSSIAN method.')
        elif (self.MIDPOINTS):
            #-----------------------------------------------------------
            # Use "midpoint displacement" to create a fractal surface
            # with correlation structure as in MARSSIM (A. Howard).
            # See midpoints.py in code directory.
            #-----------------------------------------------------------
            nn = max(nx, ny)
            n_levels = np.ceil(np.log(nn-1) / np.log(2))
            n_levels = np.int16( n_levels )
            surf = midpoints.make_fractal_surface( n_levels, H=1.5,
                                                   scale=self.noise_scale,
                                                   seed=self.seed,
                                                   SILENT=True)
            self.DEM = surf[0:ny, 0:nx]
            if not(SILENT):
                print('Created initial DEM by MIDPOINTS method.')
        else:    
            self.DEM = float32(0)

        #----------------------------
        # Construct x and y grids ?
        #----------------------------
        if (self.PLANE or self.CORNER_PLANE):
            #---------------------------
            # Option 1 for x & y grids 
            #---------------------------
            x_vals = self.dx * np.arange( nx )
            y_vals = self.dy * (ny - 1 - np.arange( ny ))
            x,y = np.meshgrid( x_vals, y_vals )

            #----------------------------------------------------
            # Note: We haven't called initialize_d8_vars() yet.
            #----------------------------------------------------
            # IDs  = self.d8.ID_grid

            #---------------------------
            # Option 2 for x & y grids 
            #---------------------------
            ## ramp = np.arange(nx*ny, dtype='Int32')
            ## IDs  = np.reshape( ramp, [ny, nx] )
            ## cols = (IDs % nx)
            ## rows = (IDS / nx)
            ## x    = (self.dx * cols)
            ## y    = (self.dy * (ny - 1 - rows))

        #--------------------------------------
        # Inclined plane tilted toward bottom
        #--------------------------------------
        if (self.PLANE):
            z    = (self.z0_plane_dz_dx * x) + \
                   (self.z0_plane_dz_dy * y)
            self.DEM += z

        #-------------------------------------------------
        # Inclined plane tilted toward lower left corner
        #-------------------------------------------------
        if (self.CORNER_PLANE):
            a = (self.z0_plane_S / sqrt(np.float32(2)))
            b = a
            z = (a * x) + (b * y)
            #--------------------------------------------
##            z    = (self.z0_plane_dz_dx * x) + \
##                   (self.z0_plane_dz_dy * y)
            self.DEM += z

        #----------------------------
        # Make sure type is FLOAT ?
        #-----------------------------------
        # Changed to "float64". (2/14/12)
        #----------------------------------
        # self.DEM = np.float32( self.DEM )
        self.DEM = np.float64( self.DEM )
        
        #-------------------------
        # Save new DEM to a file
        #--------------------------------------------------
        # Uses "write_grid() *function*, not class method
        # in rtg_files.py.  That function has an RTG_type
        # keyword and converts type before writing.
        #--------------------------------------------------
        rtg_files.write_grid( self.DEM, self.z0_file, self.rti,
                              RTG_type='DOUBLE')  ## (2/15/12)

    #   create_initial_DEM()
    #-------------------------------------------------------------
    def initialize_boundary_conditions(self):

        nx = self.nx   # (local synonyms)
        ny = self.ny
        ID_type = 'Int32'
        
        #------------------------------------
        # Use bottom row/edge as base level
        #------------------------------------
        if (self.BOTTOM):    
            self.base_IDs = np.arange(nx, dtype=ID_type) + (nx * (ny - 1))
            #*** above_row = (base_IDs - nx)
            #*** base_IDs  = [above_row, base_IDs]
            
            #---------------------------------
            # Change values in the top row ?
            #----------------------------------
            # Copy values so slope & fluxes
            # are zero at top edge for PLANE1
            #----------------------------------
            top_IDs  = np.arange(nx, dtype=ID_type)
            row2_IDs = np.arange(nx, dtype=ID_type) + nx
            self.DEM[top_IDs] = self.DEM[row2_IDs]
            #*** self.DEM[top_IDs] = 0.0
        
        #-------------------------------
        # Use right edge as base level
        #-------------------------------
        if (self.RIGHT):    
            self.base_IDs = nx * (np.arange(ny, dtype=ID_type) + 1)
            self.base_IDs = self.base_IDs - 1
            #*** prev_col = (base_IDs - 1L)
            #*** base_IDs = [prev_col, base_IDs]
            
            #---------------------------------
            # Change values on the left edge
            #---------------------------------
            left_IDs = np.arange(ny, dtype=ID_type) * nx
            #** emax = np.maximum(DEM, /NAN)
            #** DEM[left_IDs] = (emax * 0.2)
            self.DEM[left_IDs] = np.float32(0)
        
        #-----------------------------------
        # Use all four sides as base level
        #-----------------------------------
        if (self.FOUR_SIDES):    
            T_IDs = np.arange(nx, dtype=ID_type)
            B_IDs = T_IDs + nx * (ny - 1)
            L_IDs = nx * (np.arange(ny - 2, dtype=ID_type) + 1)
            R_IDs = L_IDs + (nx - 1)
            self.base_IDs = np.concatenate((T_IDs, B_IDs, L_IDs, R_IDs))
        
        #---------------------------------------
        # Use bottom left corner as base level
        #---------------------------------------
        if (self.CORNER):    
            ID1 = nx * (ny - 1)    # (lower left pixel)
            ID2 = ID1 - nx         # (just above ID1)
            ID3 = ID1 - (2 * nx)   # (just above ID2)
            self.base_IDs = np.concatenate(( ID1, ID1 + 1, ID1 + 2,
                                          ID2, ID2 + 1, ID2 + 2,
                                          ID3, ID3 + 1, ID3 + 2 ))
        
        #--------------------------------------
        # Set the initial base-level height ?
        #--------------------------------------
        if (self.BOTTOM or self.RIGHT or self.CORNER or self.FOUR_SIDES):    
            #-------------------------------------
            # Subtracting 1 here is not good for
            # continuing on from a previous run
            #-------------------------------------
            self.base_level = np.nanmin(self.DEM)
        else:    
            self.base_level = np.float32(0)
        
    #   initialize_boundary_conditions()  
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        self.dt = self.initialize_scalar( 0.1, dtype='float32') # [years]
        ## self.dt = float32(0.1)  # [years]

        #-------------------------
        # For mass balance check
        #-------------------------
        self.vol_R = self.initialize_scalar( 0, dtype='float64')
        self.vol_U = self.initialize_scalar( 0, dtype='float64')

        self.dz_max_vec = np.zeros([self.n_steps], dtype='Float32')
        self.dz_max     = self.initialize_scalar( -9999, dtype='float32')
        self.dz_min     = self.initialize_scalar(  9999, dtype='float32')

        #---------------------------------------------------
        # (11/15/11) Added these for write_output_files().
        # Initial values cause first grid to be saved.
        #---------------------------------------------------
        self.last_grid_time  = -(self.save_grid_dt   + 1)
        self.last_pixel_time = -(self.save_pixels_dt + 1)
        # self.last_model_time = -1.0  # (Don't use 0.)
        
    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def update_R(self):

        #------------------------------------------------------
        # Note: self.R is currently set by read_config_file()
        #------------------------------------------------------
        return

        ##################################################
        ##################################################
        ##  CONVERT UNITS FROM [m/s] to [m/yr] BELOW !!
        ##################################################
        ##################################################
    
        #----------------------------------------
        # Compute the "excess rainrate", R.
        # Each term must have same units: [m/s]
        # Sum = net gain/loss rate over pixel.
        #----------------------------------------------------
        # R can be positive or negative.  If negative, then
        # water is removed from the surface at rate R until
        # surface water is consumed.
        #--------------------------------------------------------------
        # P  = precip_rate   [m/s]  (converted by read_input_data()).
        # SM = snowmelt rate [m/s]
        # GW = seep rate     [m/s]  (water_table intersects surface)
        # ET = evap rate     [m/s]
        # IN = infil rate    [m/s]
        # MR = icemelt rate  [m/s]
        #--------------------------------------------------------------        
        P_rain = self.P_rain
        
        #--------------
        # For testing
        #--------------        
##        print '(Pmin,  Pmax)  =', P.min(),  P.max()
##        print '(SMmin, SMmax) =', SM.min(), SM.max()
##        print '(GWmin, GWmax) =', GW.min(), GW.max()
##        print '(ETmin, ETmax) =', ET.min(), ET.max()
##        print '(INmin, INmax) =', IN.min(), IN.max()
##        print '(MRmin, MRmax) =', MR.min(), MR.max()
##        # print '(Hmin,  Hmax)  =', H.min(), H.max()
##        print ' '

        self.R = P_rain
        ## self.R = (P_rain + SM + GW + MR) - (ET + IN)
            
    #   update_R()
    #-------------------------------------------------------------------
    def update_R_integral(self):

        #-----------------------------------------------
        # Update mass total for R, sum over all pixels
        #-----------------------------------------------   
        ## volume = np.double(self.R * self.da * self.dt)  # [m^3]
        volume = np.float64(self.R * self.da * self.dt)  # [m^3]
        if (volume.size == 1):
            self.vol_R += (volume * self.rti.n_pixels)
        else:
            self.vol_R += np.sum(volume)

    #   update_R_integral()
    #-------------------------------------------------------------------
    def update_U(self):

        #------------------------------------------------------
        # Note: self.U is currently set by read_config_file()
        #------------------------------------------------------
        return
        
    #   update_U()     
    #-------------------------------------------------------------------
    def update_U_integral(self):

        #-----------------------------------------------
        # Update mass total for U, sum over all pixels
        #-----------------------------------------------
        ## volume = np.double(self.U * self.da * self.dt)  # [m^3]
        volume = np.float64(self.U * self.da * self.dt)  # [m^3]
        if (volume.size == 1):
            self.vol_U += (volume * self.rti.n_pixels)
        else:
            self.vol_U += np.sum(volume)

    #   update_U_integral() 
    #-------------------------------------------------------------------
    def update_base_level(self):

        ###################################################
        # NB! As written, this makes all base_IDs have
        # the same elevation.  Probably not what we want.
        # No longer called in erode_d8_global.py.
        ###################################################
        
        #--------------------------------------
        # Lower base level for bottom row or
        # rightmost column or LL corner, etc.
        #--------------------------------------
        if (self.BOTTOM or self.RIGHT or self.CORNER or self.FOUR_SIDES):    
            #---------------------------------------
            # NB!  Inside loop since dt is dynamic
            #---------------------------------------
            # BLR_mpyr has units of meters/yr.
            #---------------------------------------
            drop = (self.BLR_mpyr * self.dt)
            self.base_level -= drop
            self.DEM.flat[ self.base_IDs ] = self.base_level
            
        #-------------------------------
        # Maintain boundary condition   (Used for Test5)
        # e.g. zero out all four edges
        #-------------------------------
        # nx = self.nx
        # ny = self.ny
        # self.DEM[0,:]     = 0.0    # (left edge)
        # self.DEM[nx-1, :] = 0.0    # (right edge)
        # self.DEM[:, 0]    = 0.0    # (top edge)
        # self.DEM[:, ny-1] = 0.0    # (bottom edge)
    
    #   update_base_level()
    #-------------------------------------------------------------------
    def update_DEM_edge_values(self):
        
        #-------------------------------------------
        # 01/16/06.  Adjust DEM edge values since
        # they can't erode and will otherwise stay
        # big and skew the color stretches.
        #-------------------------------------------
        if (self.BOTTOM):
            self.DEM[0,:] = self.DEM[1,:] + np.float32(0.01)
        if (self.RIGHT):
            self.DEM[:,0] = self.DEM[:,1] + np.float32(0.01)

##        if (self.BOTTOM):
##            self.DEM[:, 1] = np.nanmax( self.DEM )
##            self.DEM[:, 0 ]= self.DEM[:, 1] + 0.01
##        if (self.RIGHT):
##            self.DEM[1, :] = np.nanmax( self.DEM )
##            self.DEM[0, :] = self.DEM[1, :] + 0.01

    #   update_DEM_edge_values()
    #-------------------------------------------------------------------
    def update_d8_vars(self, SILENT=True, REPORT=False,
                       SAVE_RTG=False):

        #--------------------------------------------------------
        # Note: This is written so that it works for both
        #       erode_d8_local.py and erode_d8_global.py,
        #       because each has its own, embedded d8.update().
        #--------------------------------------------------------
        # Update the D8 flow grid and all vars that
        # depend on it, including D8 area grid.
        #---------------------------------------------
        # Area grid units are either 'm^2' or 'km^2'
        # based on a setting in "*_d8.cfg" file.
        # All length units are given in meters.
        #---------------------------------------------
        # d8.update() needs a depression-filled DEM
        # and can later get it from a CCA port.
        #---------------------------------------------        
        self.d8.update( self.time, DEM=self.DEM,
                        SILENT=SILENT, REPORT=REPORT )

        #----------------------------------------
        # Erode model needs A_units to be "m^2"
        #----------------------------------------
        if (self.d8.A_units == 'km^2'):
            self.d8.A = self.d8.A * 1e6   # [km^2 -> m^2]
            
        #-----------------------------
        # Save grid as an RTG file ?
        #-----------------------------
        if (SAVE_RTG):
            d8_file = (self.case_prefix + '_flow.rtg')
            rtg_files.write_grid( self.d8.d8_grid, d8_file, self.rti,
                                  RTG_type='BYTE')
            area_file = (self.case_prefix + '_area.rtg')
            rtg_files.write_grid( self.d8.A, area_file, self.rti)

    #   update_d8_vars()
##    #-------------------------------------------------------------------
      # Note:  These next few functions are implemented differently
      #        by erode_d8_local.py and erode_d8_global.py.
##    #-------------------------------------------------------------------
##    def update_slope_grid(self, SILENT=True, REPORT=False):
##
##    #   update_slope_grid()
##    #-------------------------------------------------------------------
##    def update_Q_grid(self, SILENT=True, REPORT=False):
##
##    #   update_Q_grid()
##    #-------------------------------------------------------------------
##    def update_Qs_grid(self, SILENT=True, REPORT=False):
##
##    #   update_Qs_grid()
##    #-------------------------------------------------------------------
##    def update_dz_dt_grid(self, SILENT=True, REPORT=False):
##                
##    #   update_dz_dt_grid()
##    #-------------------------------------------------------------------
##    def update_DEM(self, SILENT=True, REPORT=False):
##
##    #   update_DEM()   
    #-------------------------------------------------------------------
    def fill_pits_in_DEM(self, SILENT=True):

        if not(SILENT):    
            print('Filling depressions in DEM...')

##        DEM_before = self.DEM.copy()  # (For testing)
        
        fill_pits.fill_pits(self.DEM, 'FLOAT', self.nx, self.ny,
                            SILENT=SILENT)

##        #--------------
##        # For testing
##        #--------------
##        w  = where(DEM_before != self.DEM)
##        nw = w[0].size
##        print 'Number of pixels changed by fill_pits =', nw
        
    #   fill_pits_in_DEM()
    #-------------------------------------------------------------------
    def update_DEM_min_and_max(self, REPORT=False):
        
        self.DEM_min = np.nanmin( self.DEM )
        self.DEM_max = np.nanmax( self.DEM )

        #------------------
        # Optional report
        #------------------
        if (REPORT):
            z_str = str(self.DEM_min) + ', ' + str(self.DEM_max)
            print('    min(z), max(z) = ' + z_str + ' [m]')
            
    #   update_DEM_min_and_max()
    #-------------------------------------------------------------------
    def update_dt_grid(self, SILENT=True, REPORT=False,
                       SAVE_RTG=False):

        self.update_dt_grid_method1( SILENT=SILENT, REPORT=REPORT,
                                     SAVE_RTG=SAVE_RTG )

##        self.update_dt_grid_method2( SILENT=SILENT, REPORT=REPORT,
##                                     SAVE_RTG=SAVE_RTG )
        
    #   update_dt_grid()
    #-------------------------------------------------------------------
    def update_dt_grid_method1(self, SILENT=True, REPORT=False,
                               SAVE_RTG=False):

        #------------------------------------------------------------
        # Notes: Compute dt value that each pixel needs to
        #        satisfy an experimental stability condition.
        #
        #        Idea is that no D8 "kid pixel" can contribute
        #        enough sediment to its D8 parent to make its
        #        parent have a higher elevation than it does.
        #
        #        dt < fac * da * (del_z / del_Qs)
        #
        #        fac    = factory of safety
        #        DEM    = current elevation grid  [m]
        #        Qs     = K * Q^m * S^n = sed. discharge [m^3/yr]
        #        da     = pixel area grid  [m^2]
        #        del_z  = elevation drops to D8 parent pixels [m]
        #        del_Qs = sed. discharge "drops" to D8 parents
        #------------------------------------------------------------
        # Notes: fac is a factory of safety.  Decreased fac from
        #        0.8 to 0.4 on 8/10/04.  Tried factor of 0.6 on
        #        8/12/04, but for m = n = 1, one or more pits would
        #        form near the outlet so decreased back to 0.4.
        #------------------------------------------------------------     
        if not(SILENT):    
            print('Updating dt_grid...')

        #-----------------------
        # Initialize some vars
        #-----------------------
        dt_too_small = np.float32(1E-2)
        # dt_too_big = 5d  # (should there be a max?)

        #----------------------------------------------------------- 
        # (8/17/10) RT Time Profiles at random pixels shows that
        # state vars are still noisy (oscillating) for some pixels
        # when we use 0.4, but much less so than when using 1.0.
        # Reduced default to 0.3.
        #-----------------------------------------------------------         
        # fac = np.float64(0.4)   # (less than 1)
        fac = np.float64(0.3)   # (less than 1)
        
        #------------------------------------------------------------- 
        # (8/17/10) With default parameters and n_steps=5000,
        # using fac = 1.0 also works.  Despite very similar
        # drainage patterns the base level drops to -56.9028
        # vs. -24.5874.  Max elevations are 2.02733 and 2.9543,
        # respectively, and distribution of elevations is different.
        # Total simulated time is 541213 vs. 218060.
        # These comments are for an erode_d8_global simulation.
        #-------------------------------------------------------------
        # fac = np.float64(1)
        
        #---------------------------
        # Compute downstream drops
        #----------------------------------------------------
        # Get a nonzero value if a pixel's elevation is
        # greater than than that of its parent pixel.
        #----------------------------------------------------
        # #####  THIS SHOULD BE TRUE EXCEPT FOR "FLATS" #######
        #----------------------------------------------------
        # Computing pIDs as: "self.d8.parent_IDs" works for
        # erode_d8_global but doesn't for erode_d8_local.
        #----------------------------------------------------        
        pIDs  = divmod(self.d8.parent_ID_grid, self.nx)
        del_z = np.maximum((self.DEM - self.DEM[ pIDs ]), 0)
        
        #------------------------
        # Experiment 2: 8/10/04
        #----------------------------------------------------
        # Get a nonzero value if the amount leaving a pixel
        # is more than the amount leaving its parent pixel
        #----------------------------------------------------
        del_Qs = np.maximum((self.Qs - self.Qs[pIDs]), 0)
        
        #--------------------------------------------------
        # Initialize dt_grid to a very long dt value.
        # This becomes the value at non-deposition sites.
        #--------------------------------------------------
        ###########################################################
        ###########################################################
        # NB! Original version had default of 1 year vs. 1e9 !!!
        ###########################################################
        ###########################################################
        self.dt_grid = np.zeros([self.ny, self.nx], dtype='float64')  # [years]
        self.dt_grid += 1e9  # [years]
        ## self.dt_grid += 10000.0 # [years]  # (this is a mid-range value)

        ## self.dt_grid = ones([self.ny, self.nx], dtype='float64')  # [years]
        
        #----------------------------
        # Find the deposition sites
        #----------------------------
        wd = np.where(np.logical_and((del_Qs > 0), (del_z > 0)))
        nd = wd[0].size 
        if (nd != 0):
            #------------------------------
            # Compute stable dt, in years
            #------------------------------
            term = fac * (del_z[wd] / del_Qs[wd])
            if (self.da.size == 1):
                self.dt_grid[wd] = term * self.da 
            else:
                self.dt_grid[wd] = term * self.da[wd]
        else:
            if not(SILENT):
                print('-------------------------------------------')
                print(' WARNING: There are no deposition sites.')
                print('-------------------------------------------')
                print(' ')
            
        #--------------------------
        # Save the min and max dt
        #--------------------------
        self.dt_min = self.dt_grid.min()
        self.dt_max = self.dt_grid.max()
        ## if not(SILENT):
        if (self.DEBUG):
            print('#########################################')
            print(' dt_min =', np.around(self.dt_min, 2))
            print(' dt_max =', np.around(self.dt_max, 2))
            print(' in update_dt_grid() (del_z method)')
            print('#########################################')
        
        #-------------------------------
        # Don't let dt get too small ?
        #-----------------------------------------------------
        # dt_min must be less than 1e-4 for the case m=1.5,
        # n=1.0, K=1.0. Making K smaller allows dt_too_small
        # to be bigger. Now K is smaller.
        #-----------------------------------------------------
        # Fixed units issue in Qs_Grid so should be able to
        # reduce dt_too_small.
        #-----------------------------------------------------
        dt_grid_min = np.nanmin(self.dt_grid)
        if (dt_grid_min < dt_too_small):   
            print('******************************************')
            print(' Aborting: Stable dt is too small.')
            print(' Computed dt = ' + str(dt_grid_min))
            print('******************************************')
            print(' ')
            sys.exit()

        #------------------
        # Optional report
        #------------------
        if (REPORT):
            # dt_grid_min = np.nanmin(self.dt_grid)
            dt_grid_max = np.nanmax(self.dt_grid)
            dt_str     = str(dt_grid_min)  + ', ' + str(dt_grid_max)
            del_z_str  = str(del_z.min())  + ', ' + str(del_z.max())
            del_Qs_str = str(del_Qs.min()) + ', ' + str(del_Qs.max())
            print('    min(dt),     max(dt)     = ' + dt_str + ' [yrs]')
            print('    min(del_z),  max(del_z)  = ' + del_z_str)
            print('    min(del_Qs), max(del_Qs) = ' + del_Qs_str)
            print(' ')

        #-----------------------------
        # Save grid as an RTG file ?
        #-----------------------------
        if (SAVE_RTG):
            RTG_file = (self.case_prefix + '_dt.rtg')
            rtg_files.write_grid( self.dt_grid, RTG_file, self.rti )
            
    #   update_dt_grid_method1()
    #-------------------------------------------------------------------
    def update_min_dz_up_grid(self, SILENT=True, REPORT=False,
                              SAVE_RTG=False):

        #-------------------------------------------------------------
        # Notes: Compute elevation drop from lowest uphill neighbor.
        #        Any pixel without an uphill neighbor, such as peaks
        #        and ridges will be assigned a value of 9999.0.
        #
        #        Pixels with a flow code of zero may either be on
        #        the edge or at bottom of a pit.  We need to know
        #        min_dz_up for the pit pixels.
        #
        #        Flats should be okay (min_dz_up == 0).  We expect
        #        multi-pixel depressions to fill from the edges
        #        inward, because Qs will be zero wherever S = 0.
        #-------------------------------------------------------------
        if not(SILENT):    
            print('Updating min_dz_up_grid...')

        #--------------
        # For testing
        #--------------
##        print 'n1 =', self.d8.n1
##        print 'n2 =', self.d8.n2
##        print 'n3 =', self.d8.n3
##        print 'n4 =', self.d8.n4
##        print 'n5 =', self.d8.n5
##        print 'n6 =', self.d8.n6
##        print 'n7 =', self.d8.n7
##        print 'n8 =', self.d8.n8
##        n0 = self.d8.noflow_IDs[0].size
##        print 'n0 =', n0 
##        print 'n_total  =', (self.d8.n1 + self.d8.n2 + self.d8.n3 +
##                             self.d8.n4 + self.d8.n5 + self.d8.n6 +
##                             self.d8.n7 + self.d8.n8 + n0)
##        print 'n_pixels =', self.rti.n_pixels
##        print ' '
        
        #-------------------------------------------------------
        # Compute smallest elevation drop from any "kid pixel"
        #-------------------------------------------------------
        # We may need to re-initialize with 9999's here each
        # time this function is called.  Not sure yet.
        #-------------------------------------------------------
        ## min_dz_up = self.min_dz_up_grid
        min_dz_up = np.zeros((self.ny, self.nx), dtype='float64')
        min_dz_up += 9999.0
        #-------------------------------------------------------------
        if (self.d8.n1 != 0):
            dz1 = self.DEM[self.d8.w1] - self.DEM[self.d8.p1]            
            min_dz_up[self.d8.p1] = np.minimum( min_dz_up[self.d8.p1], dz1)
        #-----------------------------------------------------------------------
        if (self.d8.n2 != 0):
            dz2 = self.DEM[self.d8.w2] - self.DEM[self.d8.p2]
            min_dz_up[self.d8.p2] = np.minimum( min_dz_up[self.d8.p2], dz2)
        #-----------------------------------------------------------------------
        if (self.d8.n3 != 0):
            dz3 = self.DEM[self.d8.w3] - self.DEM[self.d8.p3]
            min_dz_up[self.d8.p3] = np.minimum( min_dz_up[self.d8.p3], dz3)
        #-----------------------------------------------------------------------
        if (self.d8.n4 != 0):
            dz4 = self.DEM[self.d8.w4] - self.DEM[self.d8.p4]
            min_dz_up[self.d8.p4] = np.minimum( min_dz_up[self.d8.p4], dz4)
        #-----------------------------------------------------------------------
        if (self.d8.n5 != 0):
            dz5 = self.DEM[self.d8.w5] - self.DEM[self.d8.p5]
            min_dz_up[self.d8.p5] = np.minimum( min_dz_up[self.d8.p5], dz5)
        #-----------------------------------------------------------------------
        if (self.d8.n6 != 0):
            dz6 = self.DEM[self.d8.w6] - self.DEM[self.d8.p6]
            min_dz_up[self.d8.p6] = np.minimum( min_dz_up[self.d8.p6], dz6)
        #-----------------------------------------------------------------------
        if (self.d8.n7 != 0):
            dz7 = self.DEM[self.d8.w7] - self.DEM[self.d8.p7]
            min_dz_up[self.d8.p7] = np.minimum( min_dz_up[self.d8.p7], dz7)
        #-----------------------------------------------------------------------
        if (self.d8.n8 != 0):
            dz8 = self.DEM[self.d8.w8] - self.DEM[self.d8.p8]
            min_dz_up[self.d8.p8] = np.minimum( min_dz_up[self.d8.p8], dz8)

        #-----------------------------------
        # Is min_dz_up negative anywhere ?
        # This should be impossible.
        #-----------------------------------
##        w_neg  = np.where(min_dz_up < 0)
##        n_neg  = w_neg[0].size
##        if (n_neg != 0):
##            print 'WARNING: min_dz_up < 0 at', n_neg, 'locations.'

        #-------------------------------
        # Is min_dz_up zero anywhere ?
        #--------------------------------------------------------
        # This should also be impossible if we're not doing
        # any "D8 flat resolution" such as "iterative linking".
        # If it did happen, then we'd probably want to look at
        # the other kids and find the "min positive dz up."
        #--------------------------------------------------------
        # This will occur if a neighbor pixel has the same
        # elevation (but flows toward this one) even if other
        # "neighbor kids" have positive drops.
        #--------------------------------------------------------
        w_zero = np.where(min_dz_up == 0)
        n_zero = w_zero[0].size
        if (n_zero != 0):
            print('WARNING: min_dz_up = 0 at', n_zero, 'locations.')

        #-----------------------------------
        # Monitor the number of pit pixels
        #-----------------------------------
        ### min_dz_up[ self.d8.noflow_IDs ] = 9999.0
        if not(SILENT):
            n0 = self.d8.noflow_IDs[0].size
            n_edge = 2*(self.nx + self.ny) - 4
            print('Number of no-flow pixels =', n0)
            print('Number in interior       =', (n0 - n_edge))
        #--------------------------------------------
        # Checked these on 2/23/10 and all are zero
        #--------------------------------------------
##        S0  = self.S[ self.d8.noflow_IDs ]
##        A0  = self.d8.A[ self.d8.noflow_IDs ]
##        Q0  = self.Q[ self.d8.noflow_IDs ]
##        Qs0 = self.Qs[ self.d8.noflow_IDs ]
##        print 'At noflow pixels:'
##        print '   Smin,  Smax  =', S0.min(), S0.max()
##        print '   Amin,  Amax  =', A0.min(), A0.max()
##        print '   Qmin,  Qmax  =', Q0.min(), Q0.max()
##        print '   Qsmin, Qsmax =', Qs0.min(), Qs0.max()
        #-----------------------------------------------------
        # Check dz_dt for edge_IDs and noflow_IDs (interior)
        #-----------------------------------------------------
        if not(SILENT):
            dz_dt0 = self.dz_dt[ self.d8.noflow_IDs ]
            w0     = np.where(dz_dt0 == 0)
            nw0    = w0[0].size
            print('At noflow pixels:')
            print('   dz_dt_min, dz_dt_max    =', dz_dt0.min(), dz_dt0.max())   
            print('   number of pixels with 0 =', nw0)
            dz_dte = self.dz_dt[ self.d8.edge_IDs ]
            print('At edge pixels:')
            print('   dz_dt_min, dz_dt_max    =', dz_dte.min(), dz_dte.max())
        
        #--------------
        # For testing
        #--------------
        if not(SILENT):
            w = np.where(min_dz_up != 9999.0)
            nw = w[0].size
            if (nw == 0):
                min_dz_up_max = 9999.0
            else:
                min_dz_up_max = min_dz_up[w].max()
            print('min_dz_up: min, max:', min_dz_up.min(), min_dz_up_max)
            print('DEM: min, max:', self.DEM.min(), self.DEM.max())
        
        #---------------------------
        # Store the result in self
        #---------------------------
        self.min_dz_up_grid = min_dz_up

        #-----------------------------
        # Save grid as an RTG file ?
        #-----------------------------
        if (SAVE_RTG):
            RTG_file = (self.case_prefix + '_min_dz_up.rtg')
            rtg_files.write_grid( min_dz_up, RTG_file, self.rti )
            # if not(SILENT):
            print('Saved min_dz_up_grid to:')
            print('   ' + RTG_file)
            print(' ')
                
    #   update_min_dz_up_grid()
    #-------------------------------------------------------------------
##    def update_min_dz_up_grid2(self, SILENT=True, REPORT=False):
##   
##        #-----------------------------------------------------------
##        # NB!  A neighbor pixel with (dz > 0) may not flow towards
##        #      the "center" pixel, so this approach won't work.
##        #-------------------------------------------------------------
##        # Notes: Compute elevation drop from lowest uphill neighbor.
##        #-------------------------------------------------------------
##        z = self.DEM
##        
##        #----------------------------------
##        # Elevations of 8 neighbor pixels
##        #----------------------------------
##        z1 = np.roll(np.roll(z, 1, axis=0), -1, axis=1)   # (upper-right)
##        z2 = np.roll(z, -1, axis=1)                          # (right)
##        z3 = np.roll(np.roll(z, -1, axis=0), -1, axis=1)  # (lower-right)
##        z4 = np.roll(z, -1, axis=0)                          # (bottom)
##        z5 = np.roll(np.roll(z, -1, axis=0), 1, axis=1)   # (lower-left)
##        z6 = np.roll(z, 1, axis=1)                           # (left)
##        z7 = np.roll(np.roll(z, 1, axis=0), 1, axis=1)    # (upper-left)
##        z8 = np.roll(z, 1, axis=0)                           # (top)
##
##        #---------------------------------------------------
##        # Grid of elevation drops *from* 8 neighbor pixels
##        #---------------------------------------------------
##        # NB!  A neighbor pixel with (dz > 0) may not flow
##        #      towards the "center" pixel!
##        #---------------------------------------------------
##        dz1 = np.maximum(z1-z, 0)
##        dz2 = np.maximum(z2-z, 0)
##        dz3 = np.maximum(z3-z, 0)
##        dz4 = np.maximum(z4-z, 0)
##        dz5 = np.maximum(z5-z, 0)
##        dz6 = np.maximum(z6-z, 0)
##        dz7 = np.maximum(z7-z, 0)
##        dz8 = np.maximum(z8-z, 0)
##
##        min_dz_up = np.minimum(dz1, dz2)
##        min_dz_up = np.minimum(min_dz_up, dz3)
##        min_dz_up = np.minimum(min_dz_up, dz4)
##        min_dz_up = np.minimum(min_dz_up, dz5)
##        min_dz_up = np.minimum(min_dz_up, dz6)
##        min_dz_up = np.minimum(min_dz_up, dz7)
##        min_dz_up = np.minimum(min_dz_up, dz8)
##        
##        #-------------------------------
##        # Is min_dz_up zero anywhere ?
##        #-------------------------------
##        w_neg  = np.where(min_dz_up < 0)
##        n_neg  = w_neg[0].size
##        if (n_neg != 0):
##            print 'WARNING: min_dz_up < 0 at', n_neg, 'locations.'
##        #-------------------------------------------------------------
##        w_zero = np.where(min_dz_up == 0)
##        n_zero = w_zero[0].size
##        if (n_zero != 0):
##            print 'WARNING: min_dz_up = 0 at', n_zero, 'locations.'
##
##        #--------------
##        # For testing
##        #--------------
##        w = where(min_dz_up != 9999.0)
##        print 'min_dz_up: min, max:', min_dz_up.min(), min_dz_up[w].max()
##        print 'DEM: min, max:', self.DEM.min(), self.DEM.max()
##        
##        #---------------------------
##        # Store the result in self
##        #---------------------------
##        self.min_dz_up_grid = min_dz_up
##        
##    #   update_min_dz_up_grid2()
    #-------------------------------------------------------------------
    def update_dt_grid_method2(self, SILENT=True, REPORT=False,
                               SAVE_RTG=False):

        #---------------------------------------------------
        # Note: Compute dt value that each pixel needs to
        #       satisfy CFL condition.
        #---------------------------------------------------
        # Where min_dz_up_grid = 0, we have S=0 and Qs=0.
        #---------------------------------------------------
        # Note: fac is a factory of safety that should be
        #       applied to the entire dt_grid.
        #---------------------------------------------------
        if not(SILENT):    
            print('Updating dt_grid...')
            
        #----------------------------------------
        # This function requires min_dz_up_grid,
        # so compute/update it now.
        #----------------------------------------
        self.update_min_dz_up_grid(SILENT=SILENT, REPORT=REPORT,
                                   SAVE_RTG=SAVE_RTG)
        
        ### fac          = np.float64(0.0002)  # (less than 1)
        ### fac          = np.float64(0.1)
        fac          = np.float64(1.0)  # (OK for global ??)
        self.dt_grid = fac * self.min_dz_up_grid / self.dz_dt
        
        #-------------------------------------------
        # Find places where dt is not well-defined.
        #-------------------------------------------
        w_bad = np.where(np.logical_or(self.min_dz_up_grid == 9999.0,
                                 self.dz_dt <= 0.0))
        n_bad = w_bad[0].size
        if (n_bad != 0):
            self.dt_grid[w_bad] = 0.0

        #--------------------------------------------------
        # Set those places to the largest stable timestep
        #--------------------------------------------------
        dt_max = self.dt_grid.max()
        if (dt_max == 0):
            dt_max = 1.0  #######
        #-------------------------------------------
        # Include where min_dz_up_grid = 0.
        #-------------------------------------------
        w_bad = np.where(self.dt_grid == 0)
        n_bad = w_bad[0].size
        if (n_bad != 0):
            self.dt_grid[w_bad] = dt_max

        self.dt_min = self.dt_grid.min()  #######
        self.dt_max = dt_max

        if not(SILENT):
            print('##########################################')
            print(' dt_min =', np.around(self.dt_min, 4))
            print(' in update_dt_grid() (min_dz_up method)')
            print('##########################################')

        #-----------------------------
        # Save grid as an RTG file ?
        #-----------------------------
        if (SAVE_RTG):
            RTG_file = (self.case_prefix + '_dt.rtg')
            rtg_files.write_grid( self.dt_grid, RTG_file, self.rti )
            
    #   update_dt_grid_method2()
    #-------------------------------------------------------------------
    def update_dt_grid_local1(self, SILENT=True, REPORT=False,
                              SAVE_RTG=False):

        #------------------------------------------------------------
        # Notes: Compute dt value that each pixel needs to
        #        satisfy an experimental stability condition.
        #
        #        Idea is that no D8 "kid pixel" can contribute
        #        enough sediment to its D8 parent to make its
        #        parent have a higher elevation than it does.
        #
        #        dt < fac * da * (del_z / del_Qs)
        #
        #        fac    = factory of safety
        #        DEM    = current elevation grid  [m]
        #        Qs     = K * Q^m * S^n = sed. discharge [m^3/yr]
        #        da     = pixel area grid  [m^2]
        #        del_z  = elevation drops to D8 parent pixels [m]
        #        del_Qs = sed. discharge "drops" to D8 parents
        #------------------------------------------------------------
        # Notes: fac is a factory of safety.  Decreased fac from
        #        0.8 to 0.4 on 8/10/04.  Tried factor of 0.6 on
        #        8/12/04, but for m = n = 1, one or more pits would
        #        form near the outlet so decreased back to 0.4.
        #------------------------------------------------------------     
        if not(SILENT):    
            print('Updating dt_grid...')

        #-----------------------
        # Initialize some vars
        #-----------------------
        fac          = np.float64(0.4)   # (less than 1)
        dt_too_small = np.float64(1E-2)
        # dt_too_big = 5d  # (should there be a max?)
        
        #---------------------------
        # Compute downstream drops
        #----------------------------------------------------
        # Get a nonzero value if a pixel's elevation is
        # greater than than that of its parent pixel.
        #----------------------------------------------------
        # #####  THIS SHOULD TRUE EXCEPT FOR "FLATS" #######
        #----------------------------------------------------
        IDs    = self.d8.IDs
        pIDs   = self.d8.parent_ID_grid.flat[ IDs ]
        z_IDs  = self.DEM.flat[ IDs ]
        z_pIDs = self.DEM.flat[ pIDs ]
        del_z  = np.maximum((z_IDs - z_pIDs), 0)
        
        #------------------------
        # Experiment 2: 8/10/04
        #----------------------------------------------------
        # Get a nonzero value if the amount leaving a pixel
        # is more than the amount leaving its parent pixel
        #----------------------------------------------------
        Qs_IDs  = self.Qs.flat[ IDs ]
        Qs_pIDs = self.Qs.flat[ pIDs ]
        del_Qs  = np.maximum((Qs_IDs - Qs_pIDs), 0)

        ###################################################
        # THIS PART VIOLATES THE IDEA OF A "LOCAL" CALC.
        ###################################################
        ###################################################
        #-----------------------------------------------
        # Old version had default of 1 year, but there
        # is no particular reason
        #-----------------------------------------------
        self.dt_grid = np.zeros([self.ny, self.nx], dtype='float64')  # [years]
        self.dt_grid += 10000.0 # [years]
        
        #----------------------------
        # Find the deposition sites
        #----------------------------
        wd = np.where(np.logical_and((del_Qs > 0), (del_z > 0)))
        nd = wd[0].size 
        if (nd != 0):
            #------------------------------
            # Compute stable dt, in years
            #------------------------------
            term = fac * (del_z[wd] / del_Qs[wd])
            if (self.da.size == 1):
                self.dt_grid.flat[IDs[wd]] = term * self.da 
            else:
                self.dt_grid.flat[IDs[wd]] = term * self.da.flat[IDs[wd]]
        else:
            print('-------------------------------------------')
            print(' WARNING: There are no deposition sites.')
            print('-------------------------------------------')
            print(' ')
            
        #--------------------------
        # Save the min and max dt
        #--------------------------
        self.dt_min = self.dt_grid.min()
        self.dt_max = self.dt_grid.max()
        if not(SILENT):
            print('#########################################')
            print(' dt_min =', np.around(self.dt_min, 2))
            print(' in update_dt_grid() (del_z method)')
            print('#########################################')
        
        #-------------------------------
        # Don't let dt get too small ?
        #-----------------------------------------------------
        # dt_min must be less than 1e-4 for the case m=1.5,
        # n=1.0, K=1.0. Making K smaller allows dt_too_small
        # to be bigger. Now K is smaller.
        #-----------------------------------------------------
        # Fixed units issue in Qs_Grid so should be able to
        # reduce dt_too_small.
        #-----------------------------------------------------
        dt_grid_min = np.nanmin(self.dt_grid)
        if (dt_grid_min < dt_too_small):   
            print('******************************************')
            print(' Aborting: Stable dt is too small.')
            print(' Computed dt = ' + str(dt_grid_min))
            print('******************************************')
            print(' ')
            sys.exit()

        #------------------
        # Optional report
        #------------------
        if (REPORT):
            # dt_grid_min = np.nanmin(self.dt_grid)
            dt_grid_max = np.nanmax(self.dt_grid)
            dt_str     = str(dt_grid_min)  + ', ' + str(dt_grid_max)
            del_z_str  = str(del_z.min())  + ', ' + str(del_z.max())
            del_Qs_str = str(del_Qs.min()) + ', ' + str(del_Qs.max())
            print('    min(dt),     max(dt)     = ' + dt_str + ' [yrs]')
            print('    min(del_z),  max(del_z)  = ' + del_z_str)
            print('    min(del_Qs), max(del_Qs) = ' + del_Qs_str)
            print(' ')

        #-----------------------------
        # Save grid as an RTG file ?
        #-----------------------------
        if (SAVE_RTG):
            RTG_file = (self.case_prefix + '_dt.rtg')
            rtg_files.write_grid( self.dt_grid, RTG_file, self.rti )
            
    #   update_dt_grid_local1()
    #-------------------------------------------------------------------
    def update_n_grid(self, SILENT=True, REPORT=False, STORE=True,
                      SAVE_RTG=False):

        #--------------------------------------------------
        # Notes: Bin the pixels according to the timestep
        #        they require to satisfy a CFL condition.
        #        This version is based on dt_min.
        #--------------------------------------------------
        # (1) Find dt_min = min(dt_grid) = the smallest
        #                   dt required for stability.
        #
        # (2) Let x = (dt_grid / dt_min) >= 1.
        #--------------------------------------------------        
        self.dt_min = self.dt_grid.min()
        x = (self.dt_grid / self.dt_min)

        #----------------------------------------------------
        # Assign value of n to every pixel that satisfies:
        #
        #       n * dt_min <= dt_grid < (n+1) * dt_min
        #  <=>  n <= x < (n+1)
        #
        # Then pixels in n_grid with a given value of n
        # will satisfy CFL for dt = (n+1) * dt_min.
        #------------------------------------------------------
        # n=1  =>  1 * dt_min <= dt_grid < 2 * dt_min
        # n=2  =>  2 * dt_min <= dt_grid < 3 * dt_min
        # n=3  =>  3 * dt_min <= dt_grid < 4 * dt_min
        #------------------------------------------------------
        ## self.n_grid = np.floor(x)

        #-------------------------------------------------------
        # Assign value of n to every pixel that satisfies:
        #
        #       (4^(n-1)) * dt_min <= dt_grid < (4^n) * dt_min
        #  <=>  (n-1) <= log(x)/log(4) < n
        #
        # Then pixels in n_grid with a given value of n
        # will satisfy CFL for dt = (4^(n-1)) * dt_min.
        #-------------------------------------------------------
        # n=1  =>  1  * dt_min <= dt_grid < 4  * dt_min
        # n=2  =>  4  * dt_min <= dt_grid < 16 * dt_min
        # n=3  =>  16 * dt_min <= dt_grid < 64 * dt_min
        #-------------------------------------------------------
        # Same idea works if we replace 4 by some other
        # positive integer.  Try changing this value to get
        # about the same number in each bin ?
        #-------------------------------------------------------
        # np.floor( np.log(x) / np.log(4), a)
        a = np.floor( np.log(x) / np.log(4))
        n_grid = (1 + a).astype('int32') 
        ## n_grid = 1 + np.floor(np.log(x)/np.log(4))

        #------------------------------------------------
        # Store the pixel IDs for all of the "n groups"
        #------------------------------------------------
        start       = 0
        n_max       = n_grid.max()
        n_pixels    = self.rti.n_pixels
        group_IDs   = np.empty(n_pixels,  dtype='int32')        
        group_start = np.zeros(n_max + 1, dtype='int32')
        group_count = np.zeros(n_max + 1, dtype='int32')
        
        for n in range(1, n_max+1):
            #-----------------------------------------------------
            # Should we use n_grid.flatten() or ravel(n_grid) ??
            #-----------------------------------------------------
            ## group_n = np.where(self.n_grid.flatten() == n)
            group_n = np.where( np.ravel(n_grid) == n )
            count   = group_n[0].size
            #-----------------------------------------------
            # Notes: When (count == 0), IDs[i:i] = [].
            #        group_n is a tuple.
            #        group_n[0] is an nd_array.
            #-----------------------------------------------
            group_IDs[start: start + count] = group_n[0]
            group_start[n] = start
            group_count[n] = count
            start         += count

        if (start != n_pixels):
            print('---------------------------------------------')
            print('ERROR: Sum over all n groups is not')
            print('       equal to number of pixels in grid.')
            print('---------------------------------------------')

        #---------------------------------
        # Store "n group" info in self ?
        #---------------------------------
        if (STORE):
            self.n_grid      = n_grid
            self.n_max       = n_max
            self.group_IDs   = group_IDs
            self.group_start = group_start
            self.group_count = group_count
        else:
            #-----------------------------------------
            # This is used to see how n_grid changes
            # inside of the update_DEM() method.
            #-----------------------------------------
            self.new_n_grid  = n_grid

        #------------------
        # Optional report
        #------------------
        ## REPORT = True   ###########
        ## REPORT = STORE  ##########
        if (REPORT):
            n_min = n_grid.min()
            print('n_min, n_max =', n_min, n_max)
            #------------------------------------------
            print('group_count =')
            print(group_count[1:])  # (skip 0 for n=0)
            #------------------------------------------
            n_pixels = self.nx * self.ny
            T_global = n_pixels * 4.0**(n_max - 1)
            n_vec    = np.arange(1,n_max+1, dtype='int64')
            M_vec    = group_count[1:]
            T_local  = np.sum( M_vec * 4**(n_max - n_vec) )
            ratio    = T_global / T_local
            print('T_global =', T_global)
            print('T_local  =', T_local)
            print('T_global / T_local =', ratio)
            print('------------------------------------------------------------')
            
        #-----------------------------
        # Save grid as an RTG file ?
        #-----------------------------
        if (SAVE_RTG):
            RTG_file = (self.case_prefix + '_n.rtg')
            rtg_files.write_grid( n_grid, RTG_file, self.rti )
            
        #------------------
        # Optional report
        #------------------
##        REPORT = True   ###########
##        if (REPORT):
##            n_min = self.n_grid.min()
##            n_max = self.n_grid.max()
##            print 'n_min, n_max =', n_min, n_max
##            #----------------------------------------
##            w1 = np.where(self.n_grid == 1)
##            print 'n=1 for', size(w1[0]), 'pixels'
##            w2 = np.where(self.n_grid == 2)
##            print 'n=2 for', size(w2[0]), 'pixels'
##            w3 = np.where(self.n_grid == 3)
##            print 'n=3 for', size(w3[0]), 'pixels'
##            w4 = np.where(self.n_grid == 4)
##            print 'n=4 for', size(w4[0]), 'pixels'
##            w5 = np.where(self.n_grid == 5)
##            print 'n=5 for', size(w5[0]), 'pixels'
##            w6 = np.where(self.n_grid == 6)
##            print 'n=6 for', size(w6[0]), 'pixels'
##            w7 = np.where(self.n_grid == 7)
##            print 'n=7 for', size(w7[0]), 'pixels'
##            w8 = np.where(self.n_grid == 8)
##            print 'n=8 for', size(w8[0]), 'pixels'
            
    #   update_n_grid()
    #-------------------------------------------------------------------
##    def update_n_grid2(self, SILENT=True, REPORT=False):
##
##        #---------------------------------------------------
##        # Notes: Bin the pixels according to the timestep
##        #        they require to satisfy a CFL condition.
##        #        This version is based on dt_max.
##        #---------------------------------------------------
##        # In the general case, n_values need not increase
##        # in the downstream direction.  However, if the
##        # slope exponent (n) equals 1, then perhaps this
##        # will be the case.
##        #---------------------------------------------------
##        # (1) Find dt_max = max(dt_grid) = the biggest
##        #                   dt required for stability.
##        #
##        # (2) Let x = (dt_grid / dt_max) <= 1.
##        #---------------------------------------------------
##        self.dt_max = self.dt_grid.max()
##        x = (self.dt_grid / self.dt_max)
##
##        #----------------------------------------------------
##        # Assign value of n to every pixel that satisfies:
##        #
##        #       dt_max / (n+1) < dt_grid <= (dt_max / n)
##        #  <=>  1/(n+1) < x <= (1/n)
##        #
##        # Then pixels in n_grid with a given value of n
##        # will satisfy CFL for dt = (dt_max / n).
##        #------------------------------------------------------
##        # n=1  =>  (1/2) * dt_max < dt_grid <= dt_max
##        # n=2  =>  (1/3) * dt_max < dt_grid <= (1/2) * dt_max
##        # n=3  =>  (1/4) * dt_max < dt_grid <= (1/3) * dt_max
##        #------------------------------------------------------
##        ## self.n_grid = np.floor(1.0 / x)
##
##        #----------------------------------------------------
##        # Assign value of n to every pixel that satisfies:
##        #
##        #       dt_max / (4^n) < dt_grid <= (dt_max / (4^(n-1))
##        #  <=>  (n-1) <= -log(x)/log(4) < n
##        #
##        # Then pixels in n_grid with a given value of n
##        # will satisfy CFL for dt = (dt_max / 4^(n-1)).
##        #---------------------------------------------------------
##        # n=1  =>  (1/4) * dt_max < dt_grid <= dt_max
##        # n=2  =>  (1/16) * dt_max < dt_grid <= (1/4)  * dt_max
##        # n=3  =>  (1/64) * dt_max < dt_grid <= (1/16) * dt_max
##        #--------------------------------------------------------
##        # Same idea works if we replace 4 by some other
##        # positive integer.  Try changing this value to get
##        # about the same number in each bin ?
##        #-------------------------------------------------------
##        self.n_grid = np.floor(1.0 - np.log(x)/np.log(4.0))
##        
##        REPORT = True   ###########
##        if (REPORT):
##            n_min = self.n_grid.min()
##            n_max = self.n_grid.max()
##            print 'n_min, n_max =', n_min, n_max
##            #----------------------------------------
##            w1 = np.where(self.n_grid == 1)
##            print 'n=1 for', size(w1[0]), 'pixels'
##            w2 = np.where(self.n_grid == 2)
##            print 'n=2 for', size(w2[0]), 'pixels'
##            w3 = np.where(self.n_grid == 3)
##            print 'n=3 for', size(w3[0]), 'pixels'
##            w4 = np.where(self.n_grid == 4)
##            print 'n=4 for', size(w4[0]), 'pixels'
##            w5 = np.where(self.n_grid == 5)
##            print 'n=5 for', size(w5[0]), 'pixels'
##            w6 = np.where(self.n_grid == 6)
##            print 'n=6 for', size(w6[0]), 'pixels'
##            w7 = np.where(self.n_grid == 7)
##            print 'n=7 for', size(w7[0]), 'pixels'
##            w8 = np.where(self.n_grid == 8)
##            print 'n=8 for', size(w8[0]), 'pixels'
##            
##    #   update_n_grid2()
    #-------------------------------------------------------------------
    def print_mins_and_maxes(self, step, DEM=False, A=False,
                             S=False, Q=False, QS=False,
                             DZ_DT=False, DZ=False, DT=False):

##        print 'step =', step, '; DEM ='
##        print np.around(self.DEM, decimals=0)
##        print ' '
        
        if (DEM):
            vmin = self.DEM.min()
            vmax = self.DEM.max()
            vstr = '; zmin, zmax ='

        if (A):      
            vmin = self.d8.A.min()
            vmax = self.d8.A.max()
            vstr = '; Amin, Amax ='

        if (S):      
            vmin = self.S.min()
            vmax = self.S.max()
            vstr = '; Smin, Smax ='

        if (Q):      
            vmin = self.Q.min()
            vmax = self.Q.max()
            vstr = '; Qmin, Qmax ='
            
        if (QS):      
            vmin = self.Qs.min()
            vmax = self.Qs.max()
            vstr = '; Qsmin, Qsmax ='

        if (DZ):      
            vmin = self.dz.min()
            vmax = self.dz.max()
            vstr = '; dz_min, dz_max ='
            
        if (DT):      
            vmin = self.dt_grid.min()
            vmax = self.dt_grid.max()
            vstr = '; dt_min, dt_max ='

        if (DZ_DT):      
            vmin = self.dz_dt.min()
            vmax = self.dz_dt.max()
            vstr = '; dz_dt_min, dz_dt_max ='
            
        print('step =', step, vstr, vmin, vmax)

    #   print_mins_and_maxes()
    #-------------------------------------------------------------------
    def check_stability(self):
        
        #----------------------------------
        # Check for one type of stability
        # (This may be obsolete now.)
        #----------------------------------
        ## if (self.dz_max > float32(200)):    # (before 4/15/10)
        if (self.dz_max > np.float32(1000)):    
            print('************************************************')
            print('Program aborted because  dz > 1000.')
            print('Time step or K is probably too large')
            print('for this set of input parameters.')
            print('   dx = ' + str(self.dx) + ' [meters]')
            print('   dy = ' + str(self.dy) + ' [meters]')
            print('   dt = ' + str(self.dt) + ' [years]')
            print('   K  = ' + str(self.K))
            print('   m  = ' + str(self.m))
            print('   n  = ' + str(self.n))
            print('************************************************')
            print(' ')
            return False
            ## sys.exit()
        else:
            return True
        
    #   check_stability()
    #-------------------------------------------------------------------
    def check_finished(self):

        #---------------------------------------------------------
        # Note: If self.DONE has already been set to True by
        #       another function or component, this function
        #       preserves that setting (see below).
        #---------------------------------------------------------
        #       CSDMS_base.run_model() also uses self.DONE as
        #       well as self.n_steps. 
        #---------------------------------------------------------
        #       TINY_DZ can occur either because dt required for
        #       stability is really small or because we have
        #       converged to a steady-state landscape.
        #---------------------------------------------------------
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

        #----------------------
        # Used before 2/5/12.
        #----------------------
        # TIMES_UP  = (self.time_index >= self.n_steps)
        ## TINY_DZ   = (self.dz_max < self.dz_tolerance)
        # self.DONE = (self.DONE or TIMES_UP)

##        self.DONE = (self.DONE or TIMES_UP or TINY_DZ)

##        if (TINY_DZ):
##            tol_str = str(self.dz_tolerance)
##            print '### WARNING: dz_max < ' + tol_str + '.'
##            ## print 'Aborting since dz_max < ' + tol_str + '.'
##            ## self.DONE = True
##            print '      dz_max =', self.dz_max, '.'

    #   check_finished()
    #-------------------------------------------------------------------
##    def check_steady_state(self):
##
##    #   check_steady_state()
    #-------------------------------------------------------------------
    def open_input_files(self):

        if (self.make_z0_method == 'READ_FILE'):
            self.z0_unit = model_input.open_file(self.z0_type, self.z0_file)

    #   open_input_files()        
    #-------------------------------------------------------------------  
    def read_input_files(self):

        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        if (self.make_z0_method == 'READ_FILE'):
            self.z0 = model_input.read_next(self.z0_unit, self.z0_type, self.rti)
            self.z0_unit.close()   #########
        
##        slopes = model_input.read_next(self.slope_unit, self.slope_type, rti)
##        if (slopes is not None): self.slopes = slopes

    #   read_input_files()     
    #-------------------------------------------------------------------  
    def close_input_files(self):

        if (self.make_z0_method == 'READ_FILE'):
            self.z0_unit.close()
            ### if (self.z0_file != ''): self.z0_unit.close()

    #   close_input_files()
     #-------------------------------------------------------------------
    def update_outfile_names(self):

        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.z_gs_file     = (self.out_directory + self.z_gs_file)
        self.S_gs_file     = (self.out_directory + self.S_gs_file)
        self.A_gs_file     = (self.out_directory + self.A_gs_file) 
        self.Q_gs_file     = (self.out_directory + self.Q_gs_file)
        self.Qs_gs_file    = (self.out_directory + self.Qs_gs_file)
        self.dz_gs_file    = (self.out_directory + self.dz_gs_file)
        self.dt_gs_file    = (self.out_directory + self.dt_gs_file) 
        self.dz_dt_gs_file = (self.out_directory + self.dz_dt_gs_file)
##        self.n_gs_file     = (self.out_directory + self.n_gs_file)
##        self.T_gs_file     = (self.out_directory + self.T_gs_file)
        #----------------------------------------------------------------
        self.z_ts_file     = (self.out_directory + self.z_ts_file)
        self.S_ts_file     = (self.out_directory + self.S_ts_file)
        self.A_ts_file     = (self.out_directory + self.A_ts_file)
        self.Q_ts_file     = (self.out_directory + self.Q_ts_file)
        self.Qs_ts_file    = (self.out_directory + self.Qs_ts_file)
        self.dz_ts_file    = (self.out_directory + self.dz_ts_file)
        self.dt_ts_file    = (self.out_directory + self.dt_ts_file)
        self.dz_dt_ts_file = (self.out_directory + self.dz_dt_ts_file)
##        self.n_ts_file     = (self.out_directory + self.n_ts_file)
##        self.T_ts_file     = (self.out_directory + self.T_ts_file)

     
##        self.z_gs_file     = (self.case_prefix + '_2D-z.rts')
##        self.S_gs_file     = (self.case_prefix + '_2D-S.rts')
##        self.A_gs_file     = (self.case_prefix + '_2D-A.rts')
##        self.Q_gs_file     = (self.case_prefix + '_2D-Q.rts')
##        self.Qs_gs_file    = (self.case_prefix + '_2D-Qs.rts')
##        self.dz_gs_file    = (self.case_prefix + '_2D-dz.rts')
##        self.dt_gs_file    = (self.case_prefix + '_2D-dt.rts')
##        self.dz_dt_gs_file = (self.case_prefix + '_2D-dz_dt.rts')
##        self.n_gs_file     = (self.case_prefix + '_2D-n.rts')
##        self.T_dt_gs_file  = (self.case_prefix + '_2D-Tnext.rts')
##        #-----------------------------------------------------------
##        self.z_ts_file     = (self.case_prefix + '_0D-z.txt')
##        self.S_ts_file     = (self.case_prefix + '_0D-S.txt')
##        self.A_ts_file     = (self.case_prefix + '_0D-A.txt')
##        self.Q_ts_file     = (self.case_prefix + '_0D-Q.txt')
##        self.Qs_ts_file    = (self.case_prefix + '_0D-Qs.txt')
##        self.dz_ts_file    = (self.case_prefix + '_0D-dz.txt')
##        self.dt_ts_file    = (self.case_prefix + '_0D-dt.txt')
##        self.dz_dt_ts_file = (self.case_prefix + '_0D-dz_dt.txt')
##        self.n_ts_file     = (self.case_prefix + '_0D-n.txt')
##        self.T_ts_file     = (self.case_prefix + '_0D-Tnext.txt')
        
    #   update_outfile_names()     
    #-------------------------------------------------------------------  
    def open_output_files(self):

        model_output.check_netcdf()    # (test import and info message)
        self.update_outfile_names()
        
        #--------------------------------------
        # Open new files to write grid stacks
        #--------------------------------------
        if (self.SAVE_Z_GRIDS):   
            model_output.open_new_gs_file( self, self.z_gs_file, self.rti,
                                           var_name='z',
                                           long_name='elevation grid',
                                           units_name='m',
                                           time_units='years')
            
        if (self.SAVE_S_GRIDS):    
            model_output.open_new_gs_file( self, self.S_gs_file, self.rti,
                                           var_name='S',
                                           long_name='slope grid',
                                           units_name='m/m',
                                           time_units='years')
        
        if (self.SAVE_A_GRIDS):    
            model_output.open_new_gs_file( self, self.A_gs_file, self.rti,
                                           var_name='A',
                                           long_name='contributing area grid',
                                           units_name='km^2',
                                           time_units='years')

        if (self.SAVE_Q_GRIDS):    
            model_output.open_new_gs_file( self, self.Q_gs_file, self.rti,
                                           var_name='Q',
                                           long_name='water discharge grid',
                                           units_name='m^3/s',
                                           time_units='years')

        if (self.SAVE_QS_GRIDS):    
            model_output.open_new_gs_file( self, self.Qs_gs_file, self.rti,
                                           var_name='Qs',
                                           long_name='sediment discharge grid',
                                           units_name='m^3/s',
                                           time_units='years')        

        if (self.SAVE_DZ_GRIDS):    
            model_output.open_new_gs_file( self, self.dz_gs_file, self.rti,
                                           var_name='dz',
                                           long_name='elevation increment grid',
                                           units_name='m',
                                           time_units='years')
            
        if (self.SAVE_DT_GRIDS):    
            model_output.open_new_gs_file( self, self.dt_gs_file, self.rti,
                                           var_name='dt',
                                           long_name='local timestep grid',
                                           units_name='yr',
                                           time_units='years')

        if (self.SAVE_DZ_DT_GRIDS):    
            model_output.open_new_gs_file( self, self.dz_dt_gs_file, self.rti,
                                           var_name='dz_dt',
                                           long_name='z rate of change grid',
                                           units_name='m/yr',
                                           time_units='years')

        # (2/15/12)
##        if (self.SAVE_UPCNT_GRIDS):    
##            model_output.open_new_gs_file( self, self.upcnt_gs_file, self.rti,
##                                           var_name='up_count',
##                                           long_name='cell update count grid',
##                                           units_name='none',
##                                           time_units='years')
            
##        if (self.SAVE_N_GRIDS):    
##            model_output.open_new_gs_file( self, self.n_gs_file, self.rti,
##                                           var_name='n',
##                                           long_name='timestep bin number',
##                                           units_name='none',
##                                           time_units='years')

##        if (self.SAVE_TN_GRIDS):    
##            model_output.open_new_gs_file( self, self.T_next_gs_file, self.rti,
##                                           var_name='T_next',
##                                           long_name='processing time grid',
##                                           units_name='yr',
##                                           time_units='years')
            
        #---------------------------------------
        # Open text files to write time series
        #---------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_Z_PIXELS):
            model_output.open_new_ts_file( self, self.z_ts_file, IDs,
                                           var_name='z',
                                           long_name='elevation grid',
                                           units_name='m',
                                           time_units='years')
            
        if (self.SAVE_S_PIXELS):
            model_output.open_new_ts_file( self, self.S_ts_file, IDs,
                                           var_name='S',
                                           long_name='slope grid',
                                           units_name='m/m',
                                           time_units='years')
            
        if (self.SAVE_A_PIXELS):
            model_output.open_new_ts_file( self, self.A_ts_file, IDs,
                                           var_name='A',
                                           long_name='contributing area grid',
                                           units_name='km^2',
                                           time_units='years')  
            
        if (self.SAVE_Q_PIXELS):
            model_output.open_new_ts_file( self, self.Q_ts_file, IDs,
                                           var_name='Q',
                                           long_name='water discharge grid',
                                           units_name='m^3/s',
                                           time_units='years')
        if (self.SAVE_QS_PIXELS):
            model_output.open_new_ts_file( self, self.Qs_ts_file, IDs,
                                           var_name='Qs',
                                           long_name='sediment discharge grid',
                                           units_name='m^3/s',
                                           time_units='years')
            
        if (self.SAVE_DZ_PIXELS):
            model_output.open_new_ts_file( self, self.dz_ts_file, IDs,
                                           var_name='dz',
                                           long_name='elevation increment grid',
                                           units_name='m/m',
                                           time_units='years')
            
        if (self.SAVE_DT_PIXELS):
            model_output.open_new_ts_file( self, self.dt_ts_file, IDs,
                                           var_name='dt',
                                           long_name='local timestep grid',
                                           units_name='yr',
                                           time_units='years')  
            
        if (self.SAVE_DZ_DT_PIXELS):
            model_output.open_new_ts_file( self, self.dz_dt_ts_file, IDs,
                                           var_name='dz_dt',
                                           long_name='z rate of change grid',
                                           units_name='m/yr',
                                           time_units='years')

        ## (2/15/12)
##        if (self.SAVE_UPCNT_PIXELS):
##            model_output.open_new_ts_file( self, self.upcnt_ts_file, IDs,
##                                           var_name='up_count',
##                                           long_name='cell update count grid',
##                                           units_name='none',
##                                           time_units='years')
            
##        if (self.SAVE_N_PIXELS):
##            model_output.open_new_ts_file( self, self.n_ts_file, IDs,
##                                           var_name='n',
##                                           long_name='timestep bin number grid',
##                                           units_name='none',
##                                           time_units='years')

##        if (self.SAVE_TN_PIXELS):
##            model_output.open_new_ts_file( self, self.T_next_ts_file, IDs,
##                                           var_name='T_next',
##                                           long_name='processing time grid',
##                                           units_name='yr',
##                                           time_units='years')
 
    #   open_output_files()
    #-------------------------------------------------------------------  
    def write_output_files(self, time=None):

        #---------------------------------------------------------
        # Notes:  This function was written to use only model
        #         time (maybe from a caller) in seconds, and
        #         the save_grid_dt and save_pixels_dt parameters
        #         read by read_config_file().
        #
        #         set_computed_input_vars() makes sure that all
        #         of the "save_dts" are larger than or equal to
        #         the process dt.
        #---------------------------------------------------------
##        print '===> dt      =', self.dt
##        print '===> T_clock =', self.T_clock
##        print '===> Time    =', self.time
        
        #-----------------------------------------
        # Allows time to be passed from a caller
        #-----------------------------------------
        if (time is None):
            time = self.time
        else:
            time = np.float64(time)
        
        #------------------------------------- 
        # Avoid string comparisons for speed
        #------------------------------------- 
        if not(self.FIXED_STEPS): 
        ## if (self.time_step_type != 'fixed'):
            #------------------------------------------------
            # (11/15/11). Save computed values based on the
            # elapsed time since last saved.  Is this the
            # best we can do to honor "save_grid_dt" ??
            #------------------------------------------------
            elapsed_grid_time = (time - self.last_grid_time)
            if (elapsed_grid_time > self.save_grid_dt):
                ## print '#### Writing frame at time =', self.time
                self.save_grids()
                self.last_grid_time = time.copy()   ## (2/7/13)
            #-----------------------------------------------------
            elapsed_pixel_time = (time - self.last_pixel_time)
            if (elapsed_pixel_time > self.save_pixels_dt):
                self.save_pixel_values()
                self.last_pixel_time = time.copy()  ## (2/7/13)
        else:
            #---------------------------------------------
            # (11/15/11)  This does not work as intended
            # for the case of adaptive timesteps.           
            #---------------------------------------------
            # Save computed values at sampled times
            #----------------------------------------
            model_time = round(time)   # (11/15/11)
            ## model_time = int(time)
            if (model_time % int(self.save_grid_dt) == 0):
                # print '#### Writing frame at time =', self.time
                self.save_grids()
            if (model_time % int(self.save_pixels_dt) == 0):
                self.save_pixel_values()

        #--------------------------        
        # SHOULD BE OBSOLETE NOW.
        #---------------------------------------------------
        # For Erode-D8-local, we need this, too. (11/15/1)
        #---------------------------------------------------
##        if (model_time == self.last_model_time):
##            return
##        self.last_model_time = model_time

    #   write_output_files()
    #-------------------------------------------------------------------  
##    def write_output_files(self, time=None):
##
##        #---------------------------------------------------------
##        # Notes:  This function was written to use only model
##        #         time (maybe from a caller) in seconds, and
##        #         the save_grid_dt and save_pixels_dt parameters
##        #         read by read_config_file().
##        #
##        #         read_config_file() makes sure that all of
##        #         the "save_dts" are larger than or equal to the
##        #         process dt.
##        #---------------------------------------------------------
##        
##        #-----------------------------------------
##        # Allows time to be passed from a caller
##        #-----------------------------------------
##        if (time is None):
##            time = self.time
##        model_time = int(time)
##        
##        #----------------------------------------
##        # Save computed values at sampled times
##        #----------------------------------------
##        if (model_time % int(self.save_grid_dt) == 0):
##            # print '#### Writing frame at time =', self.time
##            self.save_grids()
##        if (model_time % int(self.save_pixels_dt) == 0):
##            self.save_pixel_values()
##        
##    #   write_output_files()
    #-------------------------------------------------------------------  
    def close_output_files(self):

        if (self.SAVE_Z_GRIDS):      model_output.close_gs_file( self, 'z')   
        if (self.SAVE_S_GRIDS):      model_output.close_gs_file( self, 'S')   
        if (self.SAVE_A_GRIDS):      model_output.close_gs_file( self, 'A')   
        if (self.SAVE_Q_GRIDS):      model_output.close_gs_file( self, 'Q')
        if (self.SAVE_QS_GRIDS):     model_output.close_gs_file( self, 'Qs')
        if (self.SAVE_DZ_GRIDS):     model_output.close_gs_file( self, 'dz')
        if (self.SAVE_DT_GRIDS):     model_output.close_gs_file( self, 'dt')
        if (self.SAVE_DZ_DT_GRIDS):  model_output.close_gs_file( self, 'dz_dt')
##        if (self.SAVE_UPCNT_GRIDS):  model_output.close_gs_file( self, 'up_count')
                
##        if (self.SAVE_N_GRIDS):      model_output.close_gs_file( self, 'n')
##        if (self.SAVE_TN_GRIDS):     model_output.close_gs_file( self, 'T_next')
        #-------------------------------------------------------------------------
        if (self.SAVE_Z_PIXELS):     model_output.close_ts_file( self, 'z')   
        if (self.SAVE_S_PIXELS):     model_output.close_ts_file( self, 'S')    
        if (self.SAVE_A_PIXELS):     model_output.close_ts_file( self, 'A')    
        if (self.SAVE_Q_PIXELS):     model_output.close_ts_file( self, 'Q')
        if (self.SAVE_QS_PIXELS):    model_output.close_ts_file( self, 'Qs')
        if (self.SAVE_DZ_PIXELS):    model_output.close_ts_file( self, 'dz')
        if (self.SAVE_DT_PIXELS):    model_output.close_ts_file( self, 'dt')
        if (self.SAVE_DZ_DT_PIXELS): model_output.close_ts_file( self, 'dz_dt')
##        if (self.SAVE_UPCNT_PIXELS): model_output.close_ts_file( self, 'up_count')
                
##        if (self.SAVE_N_PIXELS):     model_output.close_ts_file( self, 'n')
##        if (self.SAVE_TN_PIXELS):    model_output.close_ts_file( self, 'T_next')
        
    #   close_output_files()              
    #-------------------------------------------------------------------  
    def save_grids(self):
        
        #---------------------------------------------------
        # Notes:  Each variable is saved as a grid whether
        #         it is a scalar or already a (2D) grid.
        #---------------------------------------------------
        if (self.SAVE_Z_GRIDS):
            model_output.add_grid( self, self.DEM, 'z' )
            
        if (self.SAVE_S_GRIDS):
            model_output.add_grid( self, self.S, 'S' )
            
        if (self.SAVE_A_GRIDS):
            model_output.add_grid( self, self.d8.A, 'A' )

        if (self.SAVE_Q_GRIDS):
            model_output.add_grid( self, self.Q, 'Q' )   

        if (self.SAVE_QS_GRIDS):
            model_output.add_grid( self, self.Qs, 'Qs' )

        if (self.SAVE_DZ_GRIDS):
            model_output.add_grid( self, self.dz, 'dz' )
            
        if (self.SAVE_DT_GRIDS):
            model_output.add_grid( self, self.dt_grid, 'dt' )

        if (self.SAVE_DZ_DT_GRIDS):
            model_output.add_grid( self, self.dz_dt, 'dz_dt' )   

##        if (self.SAVE_UPCNT_GRIDS):
##            model_output.add_grid( self, self.update_count, 'up_count' )
            
##        if (self.SAVE_N_GRIDS):
##            model_output.add_grid( self, self.n_grid, 'n' )
##            
##        if (self.SAVE_T_GRIDS):
##            model_output.add_grid( self, self.T_next, 'T_next' )            
            
    #   save_grids()
    #-------------------------------------------------------------------  
    def save_pixel_values(self):   ##### save_time_series_data(self)  #######

        IDs  = self.outlet_IDs
        time = self.time
        
        if (self.SAVE_Z_PIXELS):
            model_output.add_values_at_IDs( self, time, self.DEM, 'z', IDs )
                    
        if (self.SAVE_S_PIXELS):
            model_output.add_values_at_IDs( self, time, self.S, 'S', IDs )
            
        if (self.SAVE_A_PIXELS):
            model_output.add_values_at_IDs( self, time, self.d8.A, 'A', IDs )
            
        if (self.SAVE_Q_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Q, 'Q', IDs )
            
        if (self.SAVE_QS_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Qs, 'Qs', IDs )

        if (self.SAVE_DZ_PIXELS):
            model_output.add_values_at_IDs( self, time, self.dz, 'dz', IDs )
                    
        if (self.SAVE_DT_PIXELS):
            model_output.add_values_at_IDs( self, time, self.dt_grid, 'dt', IDs )
            
        if (self.SAVE_DZ_DT_PIXELS):
            model_output.add_values_at_IDs( self, time, self.dz_dt, 'dz_dt', IDs )

##        if (self.SAVE_UPCNT_PIXELS):
##            model_output.add_values_at_IDs( self, time, self.update_count, 'up_count', IDs )
            
##        if (self.SAVE_N_PIXELS):
##            model_output.add_values_at_IDs( self, time, self.n_grid, 'n', IDs )
##            
##        if (self.SAVE_TN_PIXELS):
##            model_output.add_values_at_IDs( self, time, self.T_next, 'T_next', IDs )
            
    #   save_pixel_values()
    #-------------------------------------------------------------------
    def print_time_and_value(self, var, var_name='dz_max',
                             units_name='[m]', interval=5.0):
                             ####  PRINT_INDEX=False):

        #------------------------------------------------------------
        # Note: Print the model time and the current value of
        #       "var" and perhaps n_steps.
        #------------------------------------------------------------
        # (2/9/12) This overrides CSDMS_base.print_time_and_value()
        # and uses self.stop_code, used by erode_d8_global.py and
        # erode_d8_local.py.
        #------------------------------------------------------------

        #--------------------------------
        # Print message about interval.
        #--------------------------------
        if (self.time_index == 0):
            istr = str(interval)
            print('Will print values every '+ istr + ' seconds.')

        elapsed_time = (time.time() - self.last_print_time)
        if (elapsed_time <= interval):
            return

        #--------------------------------------------
        # Get time info and build time units string
        #--------------------------------------------
        index = (self.time_index + 1)        # (starts at 0)
        if (self.time_units == 'seconds'):
            cur_time = self.time_min
            time_units_str = ' [min]'
        else:
            cur_time = self.time
            time_units_str = ' [' + self.time_units + ']'

        #----------------------------
        # Build the variable string
        #----------------------------
        var_str = var_name + ' = ' + ("%10.5f" % var)
        var_str += ' ' + units_name          

        #------------------------
        # Build the time string
        #------------------------
        msg1 = 'Time = ' + ("%10.2f" % cur_time)
        msg2 = '  n = ' + str(index)
        #-----------------------------------------
        if (self.stop_code == 0): 
            #-------------------------------------
            # Will stop when (n >= self.n_steps)
            #-------------------------------------
            msg1 += time_units_str
            msg1 += ', ' + var_str
            #-----------------------------------
            msg2 += ' of ' + str(self.n_steps)
            print(msg1)
            print(msg2)
        elif (self.stop_code == 1): 
            #------------------------------------------
            # Will stop when (time >= self.stop_time)
            #------------------------------------------
            msg1 += ' of ' + ("%10.2f" % self.stop_time) 
            msg1 += time_units_str
            print(msg1)
            print(msg2)
            print('  ' + var_str) 
        else:
            print(msg1)
            print(msg2)
            print('  ' + var_str) 

        #-------------------------
        # Update last_print_time
        #-------------------------
        self.last_print_time = time.time()
    
    #   print_time_and_value()
    #-------------------------------------------------------------------


#---------------------------------------------------------------------------------
# (2/15/12) These are not used by erode_d8_local.py but could possibly
# still be used by erode_d8_global.py.  They are functions not class methods.
#---------------------------------------------------------------------------------
##def Get_Timestep(cmin, nmin, dx, Rmax, Amax, Smax):
##
##    #------------------------------------------------------------
##    # Notes: The Courant condition: v < dx/dt is used together
##    #        with the following equations to compute a stable
##    #        timestep:
##
##    #        (1)  Q = R * A = v * w * d
##    #        (2)  v = d^(2/3) * S^(1/2) / n
##    #        (2b) d = (n v )^(3/2) * S^(-3/4)
##
##    #        Combining (1) and (2) we get:
##    #
##    #        (3) v  = (R * A / w)^(2/5) * S^(3/10) * n^(-3/5)
##    #        (4) w  = c * dx     (c <= 1)
##    #        (5) v  < dx / dt
##
##    #        Combining these and solving for dt we get:
##
##    #        (6) dt < [c^(2/5) * n^(3/5) * dx^(7/5)] /
##    #                 [(R * A)^(2/5) * S^(3/10)]
##
##    #        Use cmin, nmin, dx_min, Rmax, Amax, Smax.
##    #------------------------------------------------------------
##    numer = (cmin ** np.float64(0.4)) * (nmin ** np.float64(0.6)) * dx ** np.float64(1.4)
##    denom = (Rmax * Amax) ** np.float64(0.4) * Smax ** np.float64(0.3)
##    dt = (numer / denom)
##    
##    return dt
##    
###   Get_Timestep()
###-----------------------------------------------------------------------
##def Stable_Timestep(A, S, dx=None, dy=None, dt=None, R=None,
##                    theta=None, m=None, n=None, k=None):
##
##    #--------------------------------------------------------
##    # Notes: This routine is based on a similarity idea for
##    #        finding a stable timestep using the fact that
##    #        the model was stable for a previous set of
##    #        parameters.  Recall that:  Qs = K Q^m S^n,
##    #        Q = R A^theta, qs = Q/dw
##
##    #        K R^m A^(theta * m) S^n dt / (ds * dw) = const
##
##    #        We also assume that dx=dy and:
##    #            ds0/ds = dw0/dw = dx0/dx = dy0/dy
##    #--------------------------------------------------------
##
##    #-----------------------
##    # Stable parameter set
##    #-----------------------
##    n_params = 2
##    dx0 = np.float32(40.0)  #(meters)
##    dy0 = np.float32(40.0)  #(meters)
##    dt0 = np.float32(10.0)  #(years)
##    R0 = np.float32(1.0)   #(m/year)
##    m0 = np.float32(1.0)
##    n0 = np.float32(1.0)
##    k0 = np.float32(0.01)
##    theta0 = np.float32(1.0)
##    P0 = np.float32(100.0)  #(Not known as well, but max(A*S).)
##    
##    #-------------------
##    # Keyword defaults
##    #-------------------
##    if (dx in [0,None]):    
##        dx = dx0
##    if (dy in [0,None]):    
##        dy = dy0
##    if (dt in [0,None]):    
##        dt = dt0
##    if (R in [0,None]):    
##        R = R0
##    if (m in [0,None]):    
##        m = m0
##    if (n in [0,None]):    
##        n = n0
##    if (k in [0,None]):    
##        k = k0
##    if (theta in [0,None]):    
##        theta = theta0
##    
##    #------------------------------------
##    # Get max value of the A-S function
##    #------------------------------------
##    grid = (A ** (theta * m)) * (S ** n)
##    P = np.nanmax(grid)
##    
##    #---------------
##    # New timestep
##    #---------------
##    dt = (dx / dx0) * (dx / dx0) * (k0 / k) * ((R0 ** m0) / (R ** m))
##    dt = dt * (P0 / P) * dt0
##    
##    return np.int16(dt)
##    
###   Stable_Timestep()
#-----------------------------------------------------------------------
