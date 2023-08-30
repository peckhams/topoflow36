#
# NB! Several GC2D parameters that are not used yet have been
#     "disabled" in set_gc2d_parameters() and in the GUI info 
#     file, such as "lapse_rate".
#
#-----------------------------------------------------------------------
# Copyright (c) 2009-2023, Scott D. Peckham
#
# Aug 2023.  Renamed H, Zi, Zb, MR throughout.
#            H = h_ice, Zi = z_ice, Zb = z_bed, MR = mr_ice.
#            MR conflicts with MR in snow_base for SVO maps
#            so added imr_gs_file, imr_ts_file synonyms.
# May 2020.  Added disable_all_output().
# Jan 2013.  Revised handling of input/output names.
# Oct 2012.  CSDMS Standard Names and BMI.
# May 2010   Changes to initialize() and read_cfg_file().
# Sep 2009.  Updates
# Aug 2009.  Initial conversion from Kessler Matlab version.
#
#-----------------------------------------------------------------------
#
#  class ice_component    (inherits from CSDMS_base)
#
#      get_component_name()
#      get_attribute()            # (10/26/11)
#      get_input_var_names()      # (5/15/12)
#      get_output_var_names()     # (5/15/12)
#      get_var_name()             # (5/15/12)
#      get_var_units()            # (5/15/12)
#      --------------------------
#      set_constants()
#      initialize_vars_matlab()
#      initialize_vars_rtg()
#      initialize_computed_vars()
#      save_matlab_dem_as_rtg()
#      --------------------------
#      initialize()
#      update()
#      finalize()
#      --------------------------
#      set_computed_input_vars()
#      set_gc2d_parameters()
#      -----------------------------
#      update_*()                ## (not used yet)
#      ------------------------
#      open_input_files()        ## (not used yet)
#      read_input_files()        ## (not used yet)
#      close_input_files()       ## (not used yet)
#      ------------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()     #####
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
#-----------------------------------------------------------------------

import numpy as np
import os      # (will be needed later for os.chdir)

from topoflow.components import gc2d

from topoflow.utils import BMI_base
from topoflow.utils import model_input
from topoflow.utils import model_output
from topoflow.utils import rtg_files

#-----------------------------------------------------------------------
class ice_component( BMI_base.BMI_component ):

    #--------------------------------------
    # Should we move these into gc2d.py ?
    #-------------------------------------------------------------------
    # NB! GC2D has a VARIABLE_DT_TOGGLE that determines whether
    #     time_step_type is "fixed" or "adaptive".  The default
    #     is VARIABLE_DT_TOGGLE = 0.
    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Ice_GC2D_Valley_Glacier',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',    ###### SEE NOTE ABOVE.
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'IceGC2D',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Ice_GC2D.cfg.in',
        'cfg_extension':      '_ice_valley_glacier.cfg',     # (matches old get_cfg_extension())
        # 'cfg_extension':      '_ice_gc2d_valley_glacier.cfg',
        'cmt_var_prefix':     '/IceGC2D/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Ice_GC2D.xml',
        'dialog_title':       'Ice: GC2D Valley Glacier Parameters',
        'time_units':         'years' }

    #--------------------------------------------------
    # Could include: 'bedrock_surface__elevation' OR
    # 'land_surface__elevation'.  Now read from file.
    #--------------------------------------------------
    _input_var_names = []  # (4/18/13)
    ### _input_var_names = ['']

    #------------------------------------------------
    # Should we use "glacier__thickness" or maybe
    # "glacier__depth", or something else ?
    #------------------------------------------------
    # Note that GC2D uses "conserveIce" to check
    # mass balance.  There is a line in the code:
    #    iceVolumeLast = conserveIce * dx * dy
    # that shows it is "grid_sum_of_thickness", or
    # the sum of the thickness in every grid cell.
    #------------------------------------------------
    # GC2D computes many other vars that could be
    # included here; see gc2d.py.
    #---------------------------------------------------
    # NB! "domain_time_integral_of_melt_volume_flux"
    #     vs. "basin_cumulative_ice_meltwater_volume".
    #     But new name doesn't connote "over basin".
    #---------------------------------------------------   
    _output_var_names = [
        'glacier_ice__domain_time_integral_of_melt_volume_flux', # vol_MR
        'glacier_ice__melt_volume_flux',     # mr_ice
        'glacier_top_surface__elevation',    # z_ice
        'glacier_ice__thickness',            # h_ice
        ## 'glacier_ice__grid_sum_of_thickness', # conserveIce
        'model_grid_cell__x_length',  # dx
        'model_grid_cell__y_length',  # dy
        'model__time_step' ]          # dt
    
    _var_name_map = {
        # 'bedrock_surface__elevation': 'Zb',  # (provide under 2 names?)
        'land_surface__elevation': 'Zb',
        #---------------------------------------------------------------
        'glacier_ice__domain_time_integral_of_melt_volume_flux': 'vol_MR',
        'glacier_ice__melt_volume_flux': 'mr_ice',
        'glacier_top_surface__elevation': 'z_ice',
        ## 'glacier_ice__grid_sum_of_thickness': 'conserveIce',
        'glacier_ice__thickness': 'h_ice',
        'model_grid_cell__x_length': 'dx',
        'model_grid_cell__y_length': 'dy',
        'model__time_step': 'dt' }
    
    _var_units_map = {
        # 'bedrock_surface__elevation': 'm',  # (provide under 2 names?)
        'land_surface__elevation': 'm',
        #---------------------------------------------------------------
        'glacier_ice__domain_time_integral_of_melt_volume_flux': 'm3',
        'glacier_ice__melt_volume_flux': 'm s-1',
        'glacier_top_surface__elevation': 'm',
        'glacier_ice__thickness': 'm',
        ## 'glacier_ice__grid_sum_of_thickness': 'm',
        'model_grid_cell__x_length': 'm',
        'model_grid_cell__y_length': 'm',
        'model__time_step': 'yr' }   ### (changed from 's' on 9/11/14.)

    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )
 
    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Ice_GC2D_Valley_Glacier'

    #   get_component_name()        
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print('###################################################')
            print(' ERROR: Could not find attribute: ' + att_name)
            print('###################################################')
            print(' ')

    #   get_attribute()
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

        #----------------------------------
        # Define some constants (12/3/09)
        #----------------------------------
        self.sec_per_year = np.float64(3600) * 24 * 365  # [secs]
        self.mps_to_mmph = np.float64(3600000)
        self.mmph_to_mps = (np.float64(1) / np.float64(3600000))
        #--------------------------------------------------------
        self.rho_H2O      = np.float64(1000)     # [kg/m**3]
        self.rho_ice      = np.float64(917)      # [kg/m**3]
        self.Cp_ice       = np.float64(2060)     # [J/(kg * K)]
        self.Qg           = np.float64(1.575e6)  # [(J/yr)/m**2]
        self.grad_Tz      = np.float64(-0.0255)  # [none ??]
        self.g            = np.float64(9.81)     # [m/s**2]
        self.DONE         = False  # (see initialize_time_vars())
        
    #   set_constants()
    #-------------------------------------------------------------------
    def initialize_vars_matlab(self):

        print('Reading GC2D input from MatLab file...')
        DEM_file = (self.in_directory + 'Animas_200.mat')
        ( H, Zb, Zi, dx, dy ) = gc2d.load_state( DEM_file,
                                                 RESTART_TOGGLE=0,
                                                 INIT_COND_TOGGLE=self.Toggles.INIT_COND_TOGGLE,
                                                 GENERIC_ICE_TOGGLE=self.Toggles.GENERIC_ICE_TOGGLE)
        ny, nx = Zb.shape
        #------------------
        self.h_ice = H
        self.z_bed = Zb
        self.z_ice = Zi
        self.dx = dx
        self.dy = dy
        self.nx = nx
        self.ny = ny

    #   initialize_vars_matlab()
    #-------------------------------------------------------------------
    def initialize_vars_rtg(self):

        print('Reading GC2D input from grid file...')

        #---------------------
        # Read the DEM as Zb
        #---------------------
        DEM_file = (self.in_directory + self.DEM_file)
        self.z_bed = rtg_files.read_grid(DEM_file, self.rti,
                                    RTG_type=self.rti.data_type)

        #-----------------------------------
        # Read or define ice depth grid, H
        #-----------------------------------------------------
        # Could search for an "H0_file" in current directory
        # to allow restarting from a previous run, etc.
        #-----------------------------------------------------
        if (self.H0_file.upper() == "NONE"):
            self.h_ice  = np.zeros((self.ny, self.nx), dtype='float32')
        else:
            H0_file = (self.in_directory + self.H0_file)
            self.h_ice  = rtg_files.read_grid(H0_file, self.rti,
                                         RTG_type='FLOAT')

        
        self.z_ice = (self.z_bed + self.h_ice)
        
    #   initialize_vars_rtg()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        #------------------------------------
        # Get initial values of H, Zb, etc.
        #-------------------------------------------------------
        # Can use ice_base.save_matlab_dem_as_rtg() to convert
        #-------------------------------------------------------
        ## self.initialize_vars_matlab()
        self.initialize_vars_rtg()
        
        #-----------------------        
        # Initialize more vars
        #-----------------------
        self.mr_ice      = np.zeros((self.ny, self.nx), dtype='float64')
        self.MR          = self.mr_ice   ## synonym
        self.meltrate    = self.mr_ice   # (useful synonym ref)
        self.vol_MR      = self.initialize_scalar(0, dtype='float64')
        self.dt_min      = np.float64(99999.0)

        ################################################
        # Need to check what "conserveIce" is exactly.
        # Is it the total ice mass ?  (2/6/13)
        ################################################        
        self.conserveIce = np.float64(0)
        
        #--------------------------------------------
        # This is set by read_cfg_file(), but is
        # not the value that GC2D uses internally.
        #--------------------------------------------
        ## self.dt = 1.0 * self.sec_per_year   # [seconds]

    #   initialize_computed_vars()  
    #-------------------------------------------------------------------
    def save_matlab_dem_as_rtg(self, prefix='Animas_200'):

        import scipy.io
        import rti_files
        
        print('Saving MatLab file DEM to RTG file...')
        mat_file = (prefix + '.mat')
        rtg_file = (prefix + '_DEM.rtg')
        
        vars = scipy.io.loadmat( mat_file )
        cellsize = np.float64(vars['cellsize'])
        easting  = np.float64(vars['easting'])
        northing = np.float64(vars['northing'])
        topo     = np.float64(vars['topo'])
        ny, nx   = topo.shape
        dx       = cellsize
        dy       = cellsize
        print('   (nx, ny) =', nx, ny)
        print('   (dx, dy) =', dx, dy)
        file_unit = open(rtg_file, 'wb')
        topo = np.float32(topo)
        data_type = 'FLOAT'
##        topo = np.float64(topo)
##        data_type = 'DOUBLE'
        topo.tofile(file_unit)
        file_unit.close()

##        class info:
##            grid_file    = rtg_file
##            data_source  = 'Converted from Animas_200.mat, M. Kessler'
##            ncols        = nx
##            nrows        = ny
##            data_type    = 'FLOAT'
##            byte_order   = 'LSB'
##            pixel_geom   = 1       # (ASSUMES fixed-length, e.g. UTM)
##            xres         = dx      # [meters]
##            yres         = dy      # [meters]
##            zres         = 1.0     #######################  CHECK THIS
##            z_units      = 'METERS'
##            y_south_edge = northing.min()
##            y_north_edge = y_south_edge + (ny * dy)
##            x_west_edge  = easting.min()
##            x_east_edge  = x_west_edge  + (nx * dx)
##            box_units    = 'METERS'
##            gmin         = topo.min()
##            gmax         = topo.max()
##            UTM_zone     = 'UNKNOWN'

        #------------------------------------
        # Create an RTI file with grid info
        #------------------------------------
        rti_files.make_info( rtg_file, ncols=nx, nrows=ny,
                             xres=dx, yres=dy,
                             data_source= 'Converted from Animas_200.mat, M. Kessler',
                             y_south_edge=northing.min(),
                             x_west_edge=easting.min(),
                             gmin=topo.min(), gmax=topo.max())
                             ### box_units='METERS', z_units='METERS')
        RTI_filename = (prefix + '.rti')
        rti_files.write_info( RTI_filename, info )
        print('Finished.')
        print(' ')
        
    #   save_matlab_dem_as_rtg()    
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        self.SILENT = SILENT
        #---------------------------------------------------
        # When a user clicks on a component's "run" button,
        # that component's "run_model()" method is called
        # which sets "mode" to "main" and passes the mode
        # to initialize().  If this initialize() method is
        # called by another component (a driver), then the
        # mode has the default setting of "nondriver".  In
        # this case, the directory, etc. is set by the
        # caller, so after "load_user_input()" below, we
        # call "set_directory()" again to override and new
        # directory setting from the user.
        #---------------------------------------------------
        if not(self.SILENT):
            print(' ')
            print('Ice component: Initializing...')

        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file

        #--------------------------------
        # Valley glacier or ice sheet ?
        #--------------------------------
        cfg_extension  = self.get_attribute( 'cfg_extension' )
        self.ICE_SHEET = ('ice_sheet' in cfg_extension.lower())
        
        #-------------------------------------------
        # Save enumerations defined in gc2d module
        # These are used by read_cfg_file().
        #-------------------------------------------
        self.MassBalance  = gc2d.MassBalance   # (enumeration)
        self.BoundaryCond = gc2d.BoundaryCond  # (enumeration)
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()        # (12/03/09)
        self.initialize_config_vars() 
        # self.read_grid_info()    # NOW IN initialize_config_vars()
        self.initialize_basin_vars()  # (5/14/10)
        #-----------------------------------------
        # This must come before "Disabled" test.
        #-----------------------------------------
        self.initialize_time_vars( units='years' )
        
        self.dx = self.rti.xres     # (Now assumes UTM coords)
        self.dy = self.rti.yres
        
        self.set_gc2d_parameters()
        
        if (self.comp_status.lower() == 'disabled'):
            if not(self.SILENT):
                print('Ice component: Disabled in CFG file.')
            self.disable_all_output()
            self.mr_ice   = self.initialize_scalar(0, dtype='float64')
            self.vol_MR   = self.initialize_scalar(0, dtype='float64')
            self.meltrate = self.mr_ice
            self.dt_min   = np.float64(99999.0)
            self.DONE     = True
            self.status   = 'initialized'
            return
  
##        #---------------------------------------------
##        # Open input files needed to initialize vars 
##        #---------------------------------------------
##        self.open_input_files()
##        self.read_input_files()
##
##        #-----------------------
##        # Initialize variables
##        #-----------------------
####        self.check_input_types()
        self.initialize_computed_vars()
        
        self.open_output_files()
        self.status = 'initialized'
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, dt=-1.0):

        #--------------------------------
        # Has component been disabled ?
        #-----------------------------------------------------
        # Note: self.mr_ice already set to 0 by initialize()
        #-----------------------------------------------------
        if (self.comp_status.lower() == 'disabled'):
            # Note: self.status should be 'initialized'.
            return

        #----------------------------------------
        # Read next met vars from input files ?
        #-----------------------------------------------------       
        # Note: read_input_files() is called by initialize()
        # and those values must be used for the "update"
        # calls before reading new ones.
        #-----------------------------------------------------
#         if (self.time_index > 0):
#             self.read_input_files()
        
        #-------------------------
        # Update computed values
        #-------------------------
        self.status = 'updating'
        # print '### CALLING GC2D.update()...'
        (dt, t, H, Zi, MR, conserveIce) = gc2d.update( self.time, self.h_ice,
                                               self.z_bed, self.dx,
                                               self.dy, self.mr_ice,
                                               ### self.dy, self.meltrate,
                                               self.conserveIce,
                                               ICEFLOW_TOGGLE=self.ICEFLOW_TOGGLE,
                                               ICESLIDE_TOGGLE=self.ICESLIDE_TOGGLE,
                                               VARIABLE_DT_TOGGLE=self.VARIABLE_DT_TOGGLE,
                                               dtDefault=self.dt_max,  ###
                                               dtMax=self.dt_max )
##                                           ICEFLOW_TOGGLE=self.Toggles.ICEFLOW_TOGGLE,
##                                           ICESLIDE_TOGGLE=self.Toggles.ICESLIDE_TOGGLE,
##                                           VARIABLE_DT_TOGGLE=self.Toggles.VARIABLE_DT_TOGGLE,
##                                           dtDefault=self.Parameters.dtDefault,
##                                           dtMax=self.Parameters.dtMax )

        #----------------------------------
        # Update computed values (FUTURE)
        #----------------------------------
        # self.status = 'updating'
        # self.update_rectangle()
        # self.update_boundary_values()
        # self.update_basal_shear_stress()
        # self.update_deformation_velocity()
        # self.update_sliding_velocity()
        # self.update_total_velocity()
        # self.update_dH_dt()
        # self.update_specific_discharge()
        # self.update_mass_balance()   # or update_bz()
        #### self.update_avalanches()
        #### self.update_calving()
        # self.update_thickness()
        # self.update_melt_rate()
        # self.update_time()
        
        if (self.mode == 'driver'):
            self.print_time_and_value(H.max(), 'H_max', '[m]')
              
        #--------------
        # For testing
        #--------------
##        print 'type(H)  =', type(H)
##        print 'type(Zi) =', type(Zi)
##        print 'type(MR) =', type(MR)
##        print ' '
        
        #------------------------------------------
        # Save computed vars in component's state
        #------------------------------------------
        self.h_ice       = H
        self.z_ice       = Zi
        self.mr_ice      = MR
        self.conserveIce = conserveIce
        
##        self.update_meltrate()
##        self.update_meltrate_integral()

        #------------------------
        # Update internal clock
        #---------------------------------------------------
        # GC2D's update() function updates the time and
        # optionally uses an adaptive timestep (in years).
        #---------------------------------------------------
        self.dt = dt
        self.dt_min = np.minimum(dt, self.dt_min)   ###

        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # Components use own self.time_sec by default.
        #-----------------------------------------------
        self.write_output_files()
        ## self.write_output_files( time_seconds )

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        self.update_time( dt )
        self.status = 'updated'  # (OpenMI)

    #   update()
    #-------------------------------------------------------------------
    def update_rectangle(self):
    
        """Only simulate on sub-rectangle that contains ice.
           Inputs: H, Zb, COMPRESS_TOGGLE, RESTART_TOGGLE, THERMAL_TOGGLE
           Outputs: Zi, compression_ration, COMPRESSED_FLAG"""

        #---------------------------------------------------------------------
        # Notes: It seems like this might have a bug, because rectangle
        #        with ice is dynamic and Zi, Zb & H get redefined each time.
        #        Better to handle some other way?
        #---------------------------------------------------------------------
        if (self.COMPRESS_TOGGLE) and (self.h_ice.max() > 1) and \
           (self.RESTART_TOGGLE != 2):
            H_FullSpace  = self.h_ice.copy()
            Zb_FullSpace = self.z_bed.copy()
        
            if (self.THERMAL_TOGGLE):
                Ts_FullSpace = self.Ts.copy()  # Ts = surface temp. ?
                Tb_FullSpace = self.Tb.copy()  # Tb = bed temp. ?
                Tm_FullSpace = self.Tm.copy()  # Tm = ??? temp.

            #[indrw,indcl] = find(H ~= 0);    # Matlab?
            indrw, indcl = np.where( self.h_ice != 0 )

            mxrw, mxcl = self.z_bed.shape

            mnrw = max( 0    , min(indrw) - 2 )
            mxrw = min( mxrw , max(indrw) + 2 )
            mncl = max( 0    , min(indcl) - 2 )
            mxcl = min( mxcl , max(indcl) + 2 )

            self.h_ice = H [ mnrw:mxrw , mncl:mxcl ]  ##### H not defined yet
            self.z_bed = Zb[ mnrw:mxrw , mncl:mxcl ]  ##### Zb not defined yet
            ## Zi = Zb + max( H, 0 )
            ## Zi = Zb + np.choose( H<0 , (H,0) )
            self.z_ice = self.z_bed + np.maximum(H, 0)

            if (self.THERMAL_TOGGLE):
                Ts = Ts[ mnrw:mxrw , mncl:mxcl ]
                Tb = Tb[ mnrw:mxrw , mncl:mxcl ]
                Tm = Tm[ mnrw:mxrw , mncl:mxcl ]

            ny, nx          = self.h_ice.shape
            mx_ny, mx_nx    = Zb_FullSpace.shape
            ny, nx          = Zb.shape
            self.compression_ratio = (mx_nx * mx_ny) / (nx * ny)
            self.COMPRESSED_FLAG   = 1
        else:
            ## Zi = Zb + max( H, 0 ) # included for restarts
            ## Zi = Zb + np.choose( H<0 , (H,0) )
            self.z_ice = self.z_bed + np.maximum(self.h_ice, 0)
            self.compression_ratio = 1.
            self.COMPRESSED_FLAG   = 0

        ### return ( Zi , compression_ratio , COMPRESSED_FLAG )
    
    #   update_rectangle()
    #-------------------------------------------------------------------
    def update_boundary_values(self):
    
        """Modify boundary cells to enforce boundary conditions.
           Inputs: H, Zb, Zi
           Outputs: H_ext, Zb_ext, Zi_ext"""
        # ( H_ext , Zb_ext , Zi_ext ) = set_bc( H , Zb , Zi )
        # ( dZidxX , dZidyY ) = difference_grid( Zi_ext , dx , dy )
        pass
    
    #   update_boundary_values()
    #-------------------------------------------------------------------
#     def update_Zi_gradient(self):
#             
#     #   update_Zi_gradient()    
    #-------------------------------------------------------------------
    def update_basal_shear_stress(self):

        """Update the basal shear stress.
           Inputs: H_ext, Zi_ext
           Outputs: H_ext, Zb_ext, Zi_ext"""
        #( ( xcmpnt , ycmpnt ) , ( taubX , taubY ) , ( HX , HY ) ) = 
        #    basal_shear_stress( H_ext , Zi_ext , dx=dx , dy=dy )
        pass
        
    #   update_basal_shear_stress()
    #-------------------------------------------------------------------
    def update_deformation_velocity(self):

        """Update ice velocity due to deformation.
           Inputs: taubX, taubY, HX, HY, xcmpnt, ycmpnt
           Outputs: UdxX, UdyY  """
           
#         if (ICEFLOW_TOGGLE):
#             ( UdxX , UdyY ) = iceflow( taubX , taubY , HX , HY , xcmpnt , ycmpnt )
#         else:
#             UdxX = np.zeros( xcmpnt.shape , dtype='float64' )    #### INEFFICIENT TO HAVE INSIDE LOOP
#             UdyY = np.zeros( ycmpnt.shape , dtype='float64' )
        pass
          
    #   update_deformation_velocity()
    #-------------------------------------------------------------------
    def update_sliding_velocity(self):

        """Update ice sliding velocity.
           Inputs: taubX, taubY, xcmpnt, ycmpnt
           Outputs: UsxX, UsyY  """
            
#         if (ICESLIDE_TOGGLE):
#             ( UsxX , UsyY ) = ice_sliding( taubX , taubY , xcmpnt , ycmpnt )
#         else:
#             #### INEFFICIENT TO HAVE INSIDE LOOP
#             UsxX = np.zeros( xcmpnt.shape , dtype='float64' )
#             UsyY = np.zeros( ycmpnt.shape , dtype='float64' )
        pass
           
    #   update_sliding_velocity()
    #-------------------------------------------------------------------
    def update_total_velocity(self):
 
        """Sum all contributions to ice motion.
           Inputs: UdxX, UdyY, UsxX, UsyY
           Outputs: UxX, UyY"""

        # ( UxX , UyY ) = sum_ice_motion( UdxX , UdyY , UsxX , UsyY )
        pass
        
    #   update_total_velocity()
    #-------------------------------------------------------------------
    def update_dH_dt(self):

        # ( dHdt , ( qxX , qyY ) ) =
        # mass_conservation( H_ext , UxX , UyY , HX , HY ,
        #                    dZidxX , dZidyY , dx=dx , dy=dy )
        pass
                                                      
    #   update_dH_dt()
    #-------------------------------------------------------------------
    def update_specific_discharge(self):
    
        pass  # (code from mass_conservation() .)
        
    #   update_specific_discharge()
    #-------------------------------------------------------------------
    def update_bz(self):

        # ( Bxy , ELA ) = mass_balance( Zi , t )
        pass
        
    #   update_bz()
    #-------------------------------------------------------------------
    def update_timestep(self):

#         if (VARIABLE_DT_TOGGLE):
#             dt = get_timestep( H , Zi_ext , Zi , dHdt , Bxy )
#         else:
#             #### INEFFICIENT TO HAVE INSIDE LOOP
#             dt = dtDefault
        pass
            
    #   update_timestep()
    #-------------------------------------------------------------------
    def update_thickness(self):

#         ( t, H, Zi, meltrate, conserveIce ) =
#              update_vars( H, Zb, Zi, Bxy, qxX, qyY, dHdt, t, dt,
#                           conserveIce, dx=dx, dy=dy )
        pass
                                                   
    #   update_thickness()
    #-------------------------------------------------------------------
    def update_melt_rate(self):
    
        pass

    #   update_melt_rate()
    #-------------------------------------------------------------------
    def update_time(self):
    
        pass

    #   update_time()
    #-------------------------------------------------------------------
    def finalize(self):

        #--------------------------------
        # Has component been disabled ?
        #--------------------------------
        if (self.comp_status.lower() == 'disabled'):
            # Note: self.status should be 'initialized'.
            return

        self.status = 'finalizing' 
        ## self.close_input_files()   ##  TopoFlow input "data streams"
        self.close_output_files()

        if not(self.SILENT):
            self.print_final_report(comp_name='Ice component')

        if (self.DEBUG):
            print('Ice: dt_min =', self.dt_min, '### Smallest dt ###')
        self.status = 'finalized'  # (OpenMI)
                
    #   finalize()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

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
    def set_gc2d_parameters(self):

        #------------------------------------------------------
        # Note: Assume input parameters have been read from a
        #       CFG file by read_cfg_file().  In order for
        #       the modules in gc2d_funcs.py to use these new
        #       values, they must be copied into its internal
        #       data structures.
        #------------------------------------------------------
        # gc2d.Parameters.dt =  # (convert from secs to years)  ######

        #----------------------------------------------------------------------
        gc2d.Parameters.tMax              = self.t_max
        gc2d.Parameters.dtMax             = self.dt_max
        gc2d.Parameters.dtDefault         = self.dt_max  ### (need this!)
        gc2d.Parameters.MinGlacThick      = self.min_glacier_thick
        gc2d.Parameters.glensA            = self.glens_A
        gc2d.Parameters.B                 = self.B
        gc2d.Parameters.UsChar            = self.char_sliding_vel
        gc2d.Parameters.taubChar          = self.char_tau_bed
        gc2d.Parameters.DepthToWaterTable = self.depth_to_water_table
        gc2d.Parameters.MaxFloatFraction  = self.max_float_fraction
        gc2d.Parameters.Hpeff             = self.Hp_eff
        gc2d.Parameters.initELA           = self.init_ELA
        gc2d.Parameters.ELAStepSize       = self.ELA_step_size
        gc2d.Parameters.ELAStepInterval   = self.ELA_step_interval
        gc2d.Parameters.gradBz            = self.grad_Bz
        gc2d.Parameters.maxBz             = self.max_Bz
        ## gc2d.Parameters.lapseRate         = self.lapse_rate    #############
        gc2d.Parameters.tmin              = self.spinup_time
        ## gc2d.Parameters.angleOfRepose     = self.angle_of_repose
        ## gc2d.Parameters.avalancheFreq     = self.avalanche_freq
        gc2d.Parameters.seaLevel          = self.sea_level
        ## gc2d.Parameters.calvingCoef       = self.calving_coeff
        gc2d.Parameters.c                 = self.Cp_ice
        gc2d.Parameters.Qg                = self.geothermal_heat_flux
        gc2d.Parameters.gradTz            = self.geothermal_gradient

        #--------------------------
        # Now set all the toggles
        #--------------------------
##        gc2d.Toggles.GUISTART_TOGGLE     = self.GUISTART_TOGGLE
##        gc2d.Toggles.SAVE_TOGGLE         = self.SAVE_TOGGLE
##        gc2d.Toggles.PLOT_TOGGLE         = self.PLOT_TOGGLE 
##        gc2d.Toggles.REPORT_TOGGLE       = self.REPORT_TOGGLE
        
        ## gc2d.Toggles.COMPRESS_TOGGLE     = self.COMPRESS_TOGGLE
        gc2d.Toggles.VARIABLE_DT_TOGGLE  = self.VARIABLE_DT_TOGGLE
        gc2d.Toggles.INIT_COND_TOGGLE    = self.INIT_COND_TOGGLE
        gc2d.Toggles.GENERIC_ICE_TOGGLE  = self.GENERIC_ICE_TOGGLE      
        gc2d.Toggles.ICEFLOW_TOGGLE      = self.ICEFLOW_TOGGLE
        gc2d.Toggles.ICESLIDE_TOGGLE     = self.ICESLIDE_TOGGLE       
        ## gc2d.Toggles.THERMAL_TOGGLE      = self.THERMAL_TOGGLE
        gc2d.Toggles.FREEZEON_TOGGLE     = self.FREEZE_ON_TOGGLE  ######
        ## gc2d.Toggles.AVALANCHE_TOGGLE    = self.AVALANCHE_TOGGLE
        ## gc2d.Toggles.CALVING_TOGGLE      = self.CALVING_TOGGLE
        ## gc2d.Toggles.ERODE_TOGGLE        = self.ERODE_TOGGLE

        #------------------------------------------------------------
        # These toggles are different because they use enumerations
        #------------------------------------------------------------
        # NB! While we can change Toggles in our own "instance" of
        #     gc2d, then are not used by the functions within
        #     gc2d_funcs.py.  So this next part doesn't work.
        #------------------------------------------------------------
        MBT = eval("self.MassBalance." + self.MASS_BALANCE_TOGGLE)
        ## print '### MASS_BALANCE_TOGGLE = MBT =', MBT
        gc2d.Toggles.MASS_BALANCE_TOGGLE = MBT
        ## print '### gc2d.Toggles.MASS_BALANCE_TOGGLE =', gc2d.Toggles.MASS_BALANCE_TOGGLE
        
        gc2d.Toggles.WEST_BC_TOGGLE      = eval("self.BoundaryCond." + self.WEST_BC_TOGGLE)
        gc2d.Toggles.EAST_BC_TOGGLE      = eval("self.BoundaryCond." + self.EAST_BC_TOGGLE)
        gc2d.Toggles.SOUTH_BC_TOGGLE     = eval("self.BoundaryCond." + self.SOUTH_BC_TOGGLE)
        gc2d.Toggles.NORTH_BC_TOGGLE     = eval("self.BoundaryCond." + self.NORTH_BC_TOGGLE)
        #------------------------------------------------------------
##        gc2d.Toggles.MASS_BALANCE_TOGGLE = self.MASS_BALANCE_TOGGLE
##        
##        gc2d.Toggles.WEST_BC_TOGGLE      = self.WEST_BC_TOGGLE
##        gc2d.Toggles.EAST_BC_TOGGLE      = self.EAST_BC_TOGGLE
##        gc2d.Toggles.SOUTH_BC_TOGGLE     = self.SOUTH_BC_TOGGLE
##        gc2d.Toggles.NORTH_BC_TOGGLE     = self.NORTH_BC_TOGGLE
        
        #--------------------------------------------
        # Save Parameters, Toggles etc. in state ??
        #--------------------------------------------
##        self.Parameters   = gc2d.Parameters    # (structure)
##        self.InputParams  = gc2d.InputParams   # (structure)
##        self.OutputParams = gc2d.OutputParams  # (structure)
##        self.Toggles      = gc2d.Toggles       # (structure)

    #   set_gc2d_parameters()
    #-------------------------------------------------------------------
    def update_outfile_names(self):

        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.hi_gs_file = (self.out_directory + self.hi_gs_file)
        self.zi_gs_file = (self.out_directory + self.zi_gs_file)
        self.mr_gs_file = (self.out_directory + self.mr_gs_file)
        #---------------------------------------------------------
        self.hi_ts_file = (self.out_directory + self.hi_ts_file)
        self.zi_ts_file = (self.out_directory + self.zi_ts_file)
        self.mr_ts_file = (self.out_directory + self.mr_ts_file)

        #-------------------------------------------------
        # 2023-08-28.  Create synonym for mr_gs_file to
        # avoid SVO name mapping conflict (svo_names.py)
        #-------------------------------------------------
        self.imr_gs_file = self.mr_gs_file
        self.imr_ts_file = self.mr_ts_file

    #   update_outfile_names()
    #-------------------------------------------------------------------  
    def disable_all_output(self):

        self.SAVE_HI_GRIDS  = False
        self.SAVE_ZI_GRIDS  = False
        self.SAVE_MR_GRIDS  = False
        #-----------------------------
        self.SAVE_HI_PIXELS = False
        self.SAVE_ZI_PIXELS = False
        self.SAVE_MR_PIXELS = False
                  
    #   disable_all_output()    
    #-------------------------------------------------------------------  
    def open_output_files(self):

        model_output.check_netcdf( SILENT=self.SILENT )
        self.update_outfile_names()
        
        #----------------------------------
        # Open files to write grid stacks
        #----------------------------------
        if (self.SAVE_HI_GRIDS):
            model_output.open_new_gs_file( self, self.hi_gs_file, self.rti,
                                           var_name='hi',
                                           long_name='ice_depth',
                                           units_name='m')
            
        if (self.SAVE_ZI_GRIDS):
            model_output.open_new_gs_file( self, self.zi_gs_file, self.rti,
                                           var_name='zi',
                                           long_name='ice_free_surface',
                                           units_name='m')
            
        if (self.SAVE_MR_GRIDS):
            model_output.open_new_gs_file( self, self.imr_gs_file, self.rti,
                                           var_name='imr',
                                           long_name='ice_meltrate',
                                           units_name='mm/hr')
                                               #####  units_name='m/s')
            
        #---------------------------------------
        # Open text files to write time series
        #---------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_HI_PIXELS):
            model_output.open_new_ts_file( self, self.hi_ts_file, IDs,
                                           var_name='hi',
                                           long_name='ice_depth',
                                           units_name='m',
                                           time_units='years')
            
        if (self.SAVE_ZI_PIXELS):
            model_output.open_new_ts_file( self, self.zi_ts_file, IDs,
                                           var_name='zi',
                                           long_name='ice_free_surface',
                                           units_name='m',
                                           time_units='years')
            
        if (self.SAVE_MR_PIXELS):
            model_output.open_new_ts_file( self, self.imr_ts_file, IDs,
                                           var_name='imr',
                                           long_name='ice_meltrate',
                                           units_name='mm/hr',
                                           time_units='years')
            
    #   open_output_files()
    #-------------------------------------------------------------------
    def write_output_files(self, time_seconds=None):

        #---------------------------------------------------------
        # Notes:  This function was written to use only model
        #         time (maybe from a caller) in seconds, and
        #         the save_grid_dt and save_pixels_dt parameters
        #         read by read_cfg_file().
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

        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
##        if ((self.time_index % self.grid_save_step) == 0):
##             self.save_grids()
##        if ((self.time_index % self.pixel_save_step) == 0):
##             self.save_pixel_values()

    #   write_output_files()
    #-------------------------------------------------------------------
    def close_output_files(self):
    
        if (self.SAVE_HI_GRIDS):  model_output.close_gs_file( self, 'hi') 
        if (self.SAVE_ZI_GRIDS):  model_output.close_gs_file( self, 'zi') 
        if (self.SAVE_MR_GRIDS):  model_output.close_gs_file( self, 'imr')  
        #-----------------------------------------------------------------
        if (self.SAVE_HI_PIXELS): model_output.close_ts_file( self, 'hi')
        if (self.SAVE_ZI_PIXELS): model_output.close_ts_file( self, 'zi')
        if (self.SAVE_MR_PIXELS): model_output.close_ts_file( self, 'imr') 

    #   close_output_files()   
    #-------------------------------------------------------------------  
    def save_grids(self):
        
        #------------------------------------------------------
        # Notes:  The RTG function will work whether argument
        #         is a scalar or already a 2D grid.
        #------------------------------------------------------
        if (self.SAVE_HI_GRIDS):
            model_output.add_grid( self, self.h_ice, 'hi', self.time_min )

        if (self.SAVE_ZI_GRIDS):
            model_output.add_grid( self, self.z_ice, 'zi', self.time_min )

        if (self.SAVE_MR_GRIDS):
            MR_mmph = self.mr_ice * self.mps_to_mmph
            model_output.add_grid( self, MR_mmph, 'imr', self.time_min )

    #   save_grids()            
    #-------------------------------------------------------------------  
    def save_pixel_values(self):

        IDs  = self.outlet_IDs
        time = self.time  # (Note that "self.time_units" = 'years')
         
        if (self.SAVE_HI_PIXELS):
            model_output.add_values_at_IDs( self, time, self.h_ice,
                                            'hi', IDs )

        if (self.SAVE_ZI_PIXELS):
            model_output.add_values_at_IDs( self, time, self.z_ice,
                                            'zi', IDs )

        if (self.SAVE_MR_PIXELS):
            MR_mmph = self.mr_ice * self.mps_to_mmph
            model_output.add_values_at_IDs( self, time, MR_mmph,
                                            'imr', IDs )
            
    #   save_pixel_values()
    #-------------------------------------------------------------------

    
        
