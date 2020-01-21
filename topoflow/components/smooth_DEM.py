
#-------------------------------------------------------------------
# Note: We can now compute a new D8 flow grid and area grid for
#       the new DEM and repeat this process until the flow grid
#       no longer changes.  Need to use d8_global.py (which
#       computes flow and area grids; used by Erode) instead of
#       tf_d8_base.py as done here.  The result will be a DEM that
#       satisfies Flint's law exactly.  For steady-state rainfall,
#       it will also satisfy a slope-discharge power law exactly.
#       We can then use this tool to explore how the network
#       topology and geometry change as we vary c and p.
#
#       Note that if S = c * A^p, Q = R * A, w = w1 * Q^b, and
#       q = q1 * S^gamma, then we have:
#           p = (1 - b)/gamma
#           c = [R^(1-b) / (w1 * q1)]^(1/gamma)
#       If b = 1/2 and gamma = -1, then p = -1/2.
#
#-------------------------------------------------------------------
import numpy as np
import os
import scipy.optimize

from topoflow.components import d8_global

from topoflow.utils import cfg_files
from topoflow.utils import BMI_base
from topoflow.utils import file_utils   # (for count_lines())
from topoflow.utils import model_output
from topoflow.utils import rtg_files
from topoflow.utils import rti_files

# import matplotlib.pyplot

#-------------------------------------------------------------------------
#   smooth_DEM.py
#
#   Copyright (c) 2005-2012, Scott D. Peckham
#   Created:   May 2004
#   Modified:  Jul-Aug 2005
#   Converted from IDL to Python:  July 2010
#   Worked on: read_profile_data() and find_best_fit_c_and_p(). (1/18/12)
#   Updated:  Aug. 2018
#
#-------------------------------------------------------------------------
#
#   unit_test()
#   curve_fit_test()
#
#   class DEM_smoother
#
#       get_component_name()
#       get_attribute()
#       set_constants()
#       initialize()
#       update()
#       finalize()
#       ----------------------
#       get_gui_info()
#       get_cfg_extension()
#       build_filenames()           ##### OSBOLETE SOON
#       initialize_d8_vars()
#       ----------------------
#       update_d8_vars()
#       update_slopes()           (Step 3)
#       update_DEM()              (Step 4)
#       ------------------------
#       read_profile_data()       (Step 1)
#       find_best_fit_c_and_p()   (Step 2)
#       ------------------------
#       open_input_files()
#       read_input_files()
#       close_input_files()
#       ------------------------
#       update_outfile_names()
#       open_output_files()
#       write_output_files
#       close_output_files()
#       save_grids()
#       save_pixel_values()
#
#-------------------------------------------------------------------------
def unit_test():

    c = DEM_smoother()
    # c.CCA   = False
    c.DEBUG = True

##    # cfg_directory = '/data/sims/erode/small_ky/'
##    cfg_directory = 'Applications/Erode/Data/Small_KY/'
##    cfg_prefix    = 'Case1'
##    c.site_prefix = 'Small'

    cfg_directory = '/home/csdms/models/erode/0.5/share/data/KY_Sub/'
    # cfg_directory = 'Applications/Erode/Data/KY_Sub/'
    cfg_prefix    = 'Smooth1'    ### 'Case1'
    c.site_prefix = 'KY_Sub'
   
    #-------------------------------------------- 
    # Note: n_steps must be read from CFG file;
    #       setting it here gets over-ridden.
    #-------------------------------------------- 
    c.run_model(cfg_directory=cfg_directory,
                cfg_prefix=cfg_prefix)
    
    ## c.initialize()
    ## c.update()

#   unit_test()
#-------------------------------------------------------------------------
def curve_fit_test():

    #------------------------------------------------------------
    # Notes:  This test function shows that the function:
    #         find_best_fit_c_and_p() works, but typically
    #         does not give the p-value to high accuracy.
    #------------------------------------------------------------

    #------------------------
    # Starting on a divide
    # and moving downstream
    #------------------------
    #** x0   = 0.001  # (IDL: doesn't converge)
    #** x0   = 0.01   # (IDL: doesn't converge)
    #** x0 = 0.1      # (IDL: converges; large stderr)
    x0   = float64(1)
    x    = arange(100, dtype='Float64') + x0     # (distance [km];  NB! x[0]=x0)
    xmin = x.min()
    xmax = x.max()
    Amax = float64(625)                # [km^2]
    ca   = Amax / xmax ** float64(2)   # [unitless]
    A    = ca * x ** 2                 # (area [km^2])
    
    #--------------------------------
    # If eps is small, then expect:
    # p = (b - 1)/2  or  b = (1 + 2p)
    #--------------------------------
    #b = -1.0d   ;(p = -1.00)
    #b = -0.9d   ;(p = -0.95)
    #b = -0.7d   ;(p = -0.85)
    #b = -0.5d    ;(p = -0.75)
    b = -float64(0.3)   #(p = -0.65)     ;(closest to actual for KY_Sub?)
    #b = -0.1d   ;(p = -0.55)
    
    #------------------------------------------
    # Make sure that z[x0] = z0.  Note that
    # if x0=0, then alog(0) will occur in fit
    # and fitting procedure will fail.
    #------------------------------------------
    z0 = np.float64(600)
    z  = z0 * (x - x0 + float64(1)) ** b     # (elevation [meters])
    
    #** eps = 1e-6
    #** z   = z0 * (x + eps)^b     ;(elevation [meters])
    #** z   = z / (1d + eps)^b     ;(so that z[0] = z0)
    #** z   = -1d * 0.01d * alog(x + 1d)   ;(elevation [meters])
    
    #---------------------------------
    # Doesn't perform well for these
    #---------------------------------
    #** z = 600d - (5.9d * x^2d)
    #** z = 600d - (5.9d * x^0.5)
    #** z = 600d - (5.9d * x)
    
    #------------------------------------
    # Reverse the vectors so that we
    # start at outlet and move upstream
    #-----------------------------------------------------------
    # Must use FLIPUD(x) vs. ROT90(x,-2) to reverse 1D arrays.
    #-----------------------------------------------------------
    x2 = np.flipud( x )
    A2 = np.flipud( A )
    z2 = np.flipud( z )
    
    #--------------------------
    # Find the best-fit curve
    #--------------------------
    c, p = best_slope_area_curve_fit( A2, z2 )
    print('best-fit c =', c)
    print('best-fit p =', p)
    #-----------------------------------
    zfit = slope_area_function( A, c, p )    # (z0 and ds via "data")
    print('zfit[0] = ', zfit[0])
    
    #----------------------------------
    # Print expected curve-fit values
    #----------------------------------
    pe = (b - float64(1)) / float64(2)
    ce = absolute((z0 * b) / (ca ** pe))   #(abs since S>0 by convention)
    print('ce =', ce)
    print('pe =', pe)
    print(' ')
    
    #---------------------------
    # Create a plot to compare
    # fitted curve to original
    #---------------------------
##    matplotlib.pyplot.figure(1, figsize=(800/80.0, 600/80.0), dpi=80)
##    matplotlib.pyplot.show()
##    loadct(39, silent=True)  #####
##    black = int16(0)
##    white = int16(255)
##    red = int16(220)
##    matplotlib.pyplot.plot(x2, z2, color='k')
##    matplotlib.pyplot.axes(axisbg='w')
##    matplotlib.pyplot.show()
##    oplot(x2, yfit, psym=-1, color=red)  ####
    
#   curve_fit_test()
#-------------------------------------------------------------------------
class DEM_smoother( BMI_base.BMI_component ):

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_DEM_Smoother'

    #   get_component_name()  
    #-----------------------------------------------------------------
    # Note: Do not define an __init__() method here.  It will
    #       override things needed from CSDMS_base.__init__()
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        map = {'comp_name':          'DEMSmoother',
               'version':            '0.5',
               'model_name':         'DEM_Smoother',
               'model_family':       'Erode',
               'cfg_template_file':  'DEM_Smoother.cfg.in',
               'cfg_extension':      '_dem_smoother.cfg', 
               'cmt_var_prefix':     '/DEMSmoother/Input/Var/',
               'gui_xml_file':       '/home/csdms/cca/erode/0.5/src/share/cmt/gui/DEM_Smoother.xml',
               'dialog_title':       'DEM Profile Smoother Parameters',
               'time_step_type':     'fixed',
               'time_units':         'years',
               'mesh_type':          'uniform',
               'author_name':        'Scott Peckham'}

        try:
            return map[ att_name.lower() ]
        except:
            print('###################################################')
            print(' ERROR: Could not find attribute: ' + att_name)
            print('###################################################')
            print(' ')

    #   get_attribute()
    #-------------------------------------------------------------------
    def set_constants(self):
       
        self.nodata = np.float32(-9999)

        #----------------------------------------------
        # Maybe set constants "c" and "p" this way ??
        # Or maybe read them as input parameters ??
        #----------------------------------------------
##        self.read_profile_data()
##        self.find_best_fit_c_and_p()

    #   set_constants()   
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        self.comp_name  = 'DEM Smoother component'
        if not(SILENT):
            print(' ')
            print(self.comp_name + ': Initializing...')
        
        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()
        self.initialize_config_vars()
        # self.build_filenames()   ##########
        # self.read_grid_info()    # NOW IN initialize_config_vars()
        self.initialize_basin_vars()
        #-----------------------------------------
        # This must come before "Disabled" test.
        #-----------------------------------------
        self.initialize_time_vars()

        #------------------------------------------------
        # Append in_directory to input files. (1/17/12)
        #------------------------------------------------
        self.DEM_file     = self.in_directory + self.DEM_file
        self.profile_file = self.in_directory + self.profile_file

        #----------------------------------------------
        # Maybe set constants "c" and "p" this way ??
        # Or maybe read them as input parameters ??
        #----------------------------------------------
        if (self.FIT_C_AND_P):
            print('Finding best-fit c and p from:')
            print('    ' + self.profile_file)
            self.read_profile_data()
            self.find_best_fit_c_and_p()
        
        #----------------------------------
        # Has component been turned off ?
        #----------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print(self.comp_name + ': Disabled in CFG file.')
            self.DONE = True
            self.status = 'initialized'  # (OpenMI 2.0 convention) 
            return
        else:
            self.DONE = False
            
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
        #-----------------------
        self.initialize_d8_vars()  # (depend on D8 flow grid)
        # self.initialize_computed_vars()
        
        self.open_output_files()
        self.status = 'initialized'  # (OpenMI 2.0 convention)
                   
    #   initialize()  
    #-------------------------------------------------------------------
    def update(self, time_seconds=None):

        #---------------------------------------------
        # Note that u and d from previous time step
        # must be used on RHS of the equations here.
        #---------------------------------------------
        self.status = 'updating'  # (OpenMI 2.0 convention)
        
##        if (self.mode == 'driver'):
##            self.print_time_and_value(self.Q_outlet, 'Q_out', '[m^3/s]')
        # self.print_time_and_value(self.Q_outlet, 'Q_out', '[m^3/s]')
            
        #-------------------------
        # Update computed values
        #-------------------------
        self.update_d8_vars()    # (compute new D8 flow and area grids)
        self.update_slopes()     # (compute new slopes from D8 areas)
        self.update_DEM()

        #-------------------------------------------
        # Read from files as needed to update vars 
        #--------------------------------------------------------
        # NB! This is currently not needed because values don't
        # change over time and read_input_files() is called by
        # initialize().
        #--------------------------------------------------------
        # if (self.time_index > 0):
        #     self.read_input_files()

        #------------------------------------------------------
        # Update internal clock *before* write_output_files()
        # because we're going to save original DEM, too, with
        # a time of zero.
        #------------------------------------------------------
        self.update_time()
        
        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        self.write_output_files( time_seconds )

        ### self.update_time()   # (after write_output_files()
        OK = True   ##### (can be used for some test)
        if (OK):
            self.status = 'updated'  # (OpenMI 2.0 convention)
        else:
            self.status = 'failed'
            self.DONE   = True
            
    #   update()   
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalizing'  # (OpenMI)
        self.close_input_files()   ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'  # (OpenMI)

        print('(c, p) = ' + str(self.c) + ', ' + str(self.p))
        print(' ')
        
        self.print_final_report(comp_name=self.comp_name)
        
    #   finalize()
    #-------------------------------------------------------------------
##    def build_filenames(self):
##
##        #--------------------------------------------------------
##        # Note: These could all be read from an input file, or
##        #       we could just prompt for prefix and new_prefix.
##        #--------------------------------------------------------
##        if (hasattr(self, 'site_prefix')):
##            prefix = self.site_prefix
##            self.DEM_file = prefix + '_DEM.rtg'
##        else:
##            prefix, extension = file_utils.get_prefix_and_extension( self.DEM_file )
##        
##        #--------------------------------------------
##        # Build input filenames from site_prefix ??
##        #--------------------------------------------
####        if (self.DEM_file == None):    
####            self.DEM_file = prefix + '_DEM.rtg'
##
##        ####################################
##        self.profile_file    = None
##        self.area_file       = None
##        self.flow_file       = None
##        self.new_DEM_file    = None
##        self.new_RTI_file    = None
##        self.new_slope_file  = None
##        self.new_rawDEM_file = None
##        self.new_flow_file   = None
##        ####################################
##        
##        if (self.profile_file == None):    
##            self.profile_file = prefix + '_prof1.txt'
##            
##        if (self.area_file == None):    
##            self.area_file = prefix + '_area.rtg'
##            #----------------------------------------------------
##            # D-infinity areas may not be monotonic increasing,
##            # and then the slopes won't decrease downstream.
##            #----------------------------------------------------
##            ### area_file = prefix + '_dinf-area.rtg'
##            
##        if (self.flow_file == None):    
##            self.flow_file = prefix + '_flow.rtg'
##
##        #----------------------------------------------
##        # Build output filenames from site_prefix ??
##        #----------------------------------------------
##        new_prefix = (prefix + '2')  #####
##        if (self.new_DEM_file == None):
##            self.new_DEM_file = new_prefix + '_DEM.rtg'
##
##        if (self.new_RTI_file == None):
##            self.new_RTI_file = new_prefix + '.rti'
##            
##        if (self.new_slope_file == None):
##            self.new_slope_file = new_prefix + '_slope.rtg'
##
##        if (self.new_rawDEM_file == None):
##            self.new_rawDEM_file = new_prefix + '_rawDEM.rtg'
##
##        if (self.new_flow_file == None):
##            self.new_flow_file = new_prefix + '_flow.rtg'
##        
##    #   build_filenames()    
    #-------------------------------------------------------------------    
    def initialize_d8_vars(self):

        #---------------------------------------------
        # Compute and store a variety of (static) D8
        # flow grid variables.  Embed structure into
        # the current component.
        #---------------------------------------------
        self.d8 = d8_global.d8_component()
        
        self.d8.DEBUG = False  # (make sure self tests are OFF)
        
        ################################################
        # (5/13/10)  Do next lines here for now, until
        # the d8 cfg_file includes site prefix.
        # Same is done in GW_base.py.
        ################################################
        # (1/17/12) Note that d8_base.py now has a new
        # method called: set_default_config_vars()
        # that is used to intialize vars in cases
        # where there is no "*_d8_global.cfg" file.
        # It is called in d8.initialize().
        ################################################
        self.d8.DEM_file        = self.DEM_file       # (1/17/12) in_directory already prepended?
        self.d8.FILL_PITS_IN_Z0 = 0                   # (1/17/12)
        self.d8.A_units         = 'km^2'              # (1/17/12) May be needed.

        #--------------------------------------------------         
        # D8 component builds its cfg filename from these  
        #--------------------------------------------------      
        self.d8.site_prefix  = self.site_prefix
        self.d8.case_prefix  = self.case_prefix   ## (07/12/18)
        self.d8.in_directory = self.in_directory
        self.d8.initialize( cfg_file=None,  
                            SILENT=self.SILENT,
                            REPORT=self.REPORT )

        #---------------------------------------------------------
        # Need "A_units" to be km^2 to compare to RT area grid
        # so override setting in the CFG file, needed for Erode.
        #---------------------------------------------------------
##        if (self.DEBUG):
##            self.d8.A_units = 'km^2'   #####
    
    #   initialize_d8_vars()
    #-------------------------------------------------------------------
    def update_d8_vars(self, SILENT=True, REPORT=False,
                       SAVE_RTG=False):

        #---------------------------------------------
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
    #-------------------------------------------------------------------------  
    def update_slopes(self):

        Amin = np.nanmin( self.d8.A )
        Amax = np.nanmax( self.d8.A ) 
        print('Min(A), Max(A) = ', Amin, ', ', Amax)

        #-------------------------------------------
        # Compute new slope grid from D8 area grid
        # using best-fit Flint's law parameters.
        #-------------------------------------------
        # S[0]=0  and  S[Inf]=0
        # Smax = (1-exp(-1)) * Ah^p * c
        #-----------------------------------

        #----------------------------------------------------
        # Idea to produce convex hilltops at Ah = hillslope
        # scale.  Otherwise get singularity at A=0.
        #----------------------------------------------------
        # self.S = self.c * (A**self.p) * (1.0 - exp(-A/Ah))
   
        #-------------------------------------------
        # Raising zero to a negative power produces
        # a "divide by zero" error message.
        # Also can't use "float64" for S.
        #-------------------------------------------
        ## self.S = self.c * (self.d8.A ** self.p)

        self.S = np.zeros( self.d8.A.shape, dtype='float32')
        wpos = where( self.d8.A > 0 )
        if (wpos[0].size > 0):
            self.S[wpos] = self.c * (self.d8.A[wpos] ** self.p)
        wneg = where (self.d8.A <= 0 ) 
        if (wneg[0].size > 0):
            self.S[wneg] = 0
        
        #----------------------------------------------------------
        # Note: These should match the slopes computed from the
        #       D8 area grid using Flint's law, but remember that
        #       we're using the original D8 flow directions which
        #       may not be strictly valid for the new DEM.
        #----------------------------------------------------------
        ## dz = (self.new_DEM - self.new_DEM.flat[ self.d8.parent_ID_grid ]) 
        ## self.S = (dz / self.d8.ds)

    #   update_slopes()
    #-------------------------------------------------------------------------
    def update_DEM(self):

        #------------------------------------------------------------------
        # NOTES: This routine uses a slope-area relationship, an area
        #        grid and a D8 flow grid to create a new, smoother DEM
        #        from an old one.  The reason for wanting to do some-
        #        thing like this is that slopes in channels are very
        #        poorly resolved by local methods, even though slopes
        #        on hillslopes can be computed reasonably well.

        #        It operates on a principle of "raster recursion" and
        #        should be adaptable to the computation of other
        #        recursively-defined quantities.  That is, it:

        #        (1) initializes the raster file, and
        #        (2) makes several passes through the file (line by
        #            line), looking for pixels whose _parent_ has a
        #            known value, and
        #        (3) assigns these pixels a value which is determined
        #            from the value of the parent.
        #
        #        Note that self.d8.ds has units of meters.
        #------------------------------------------------------------------
        #
        #   This routine is based on the one used to compute flow
        #   distances to a set of masked pixels (e.g. pixels with
        #   undefined flow codes.)  They use a type of raster recursion,
        #   but unlike most others, they work upstream, from parents to
        #   their kids, instead of downstream from kids to parent pixels.
        #
        #------------------------------------------------------------------ 
        info = self.rti
        ## info = rti_files.read_info( self.DEM_file )
        nx = info.ncols
        ny = info.nrows
        ## byte_order = info.byte_order

        #--------------------------------------------
        # Read the DEM, area and D8 flow code grids
        #--------------------------------------------
        # self.DEM   = rtg_files.read_grid( DEM_file,  info, RTG_type=info.data_type )
        # self.areas = rtg_files.read_grid( area_file, info, RTG_type='FLOAT' )
        # self.codes = rtg_files.read_grid( flow_file, info, RTG_type='BYTE' )
        
        #----------------------------------------------
        # Get a where-style tuple of parent pixel IDs
        #-------------------------------------------------------
        # Find the pixels that flow to a nodata or edge pixel;
        # the parents of these pixels have a flow code of 0.
        #-------------------------------------------------------
        pIDs         = self.d8.parent_IDs  # (where-style tuple)
        parent_codes = self.d8.d8_grid[ pIDs ]
        w  = np.where(parent_codes == 0)
        nw = w[0].size  # (much faster)
        ## nw = np.size( w[0] )
        #---------------------------------------------
        # OLD METHOD that can't handle nodata pixels
        #---------------------------------------------
        ## w = where(codes == 0)
        ## w = where(logical_and(codes == 0, DEM > self.nodata) )
        ## nw = w[0].size  
        #--------------------------------------
        # Are there any pixels to work with ?
        #--------------------------------------
        if (nw == 0):    
            print('ERROR: ')
            print('No pixels to initialize iteration.')
            print(' ')
            return

        #-------------------------------------------------
        # Initialize values in new DEM to be same as in
        # old DEM for pixels whose parent flow code = 0
        # and nodata value otherwise
        #-------------------------------------------------
        self.DEM = zeros([ny, nx], dtype='Float32') + self.nodata
        self.DEM[ w ] = self.z0[ w ]
        #----------------------------------------------------------------
        self.flow_dist = zeros([ny, nx], dtype='Float32') + self.nodata
        self.flow_dist[ w ] = 0
        #----------------------------------------------------------------        
        n_reps = np.int32(0)
        DONE   = False
        
        #------------------------------------------
        # Iteratively assign new elevation values
        #------------------------------------------
        while True:
            
            STILL_ACTIVE = False
            IDs = where( self.DEM == self.nodata )  # (tuple)
            n_IDs = IDs[0].size
            ## n_IDs = np.size(IDs[0])  # (much slower)
            n_reps += 1
            
            if (n_IDs != 0):    
                #-------------------------------------
                # Get elevations of D8 parent pixels
                #-------------------------------------
                ## dvals  = self.d8.d8_grid[ IDs ]  # (not needed)
                pIDs   = self.d8.parent_ID_grid[ IDs ]
                p_elev = self.DEM.flat[ pIDs ]
                p_dist = self.flow_dist.flat[ pIDs ]   ####
                
                #-------------------------------------
                # If D8 parent elevation is known,
                # then assign elevations to D8 kids.
                #-------------------------------------
                wp = where( p_elev != self.nodata )
                n_assigned = wp[0].size  # (much faster)
                ## n_assigned = size(wp[0])  # (much slower)
                if (n_assigned != 0):
                    #----------------------------------------------
                    # Get "calendar-style" indices of "ready IDs"
                    #----------------------------------------------
                    ID_rows   = IDs[0]
                    ID_cols   = IDs[1]
                    ID_vals   = (ID_rows * nx) + ID_cols
                    ready_IDs = ID_vals[ wp ]
                    
                    #--------------------------------
                    # Compute new slopes from areas
                    #--------------------------------
                    S_vals = self.S.flat[ ready_IDs ]
                    
                    #--------------------------------------
                    # Get upstream areas of parent's kids
                    # and compute new slopes from areas
                    #--------------------------------------
                    #### Avals = self.d8.A.flat[ ready_IDs ]  # (later on)
                    ## A_vals = self.areas.flat[ ready_IDs ]
                    ## S_vals = self.c * (A_vals ** self.p)
                    
                    #-----------------------------------
                    # S(0)=0  and  S(Inf)=0
                    # Smax = (1-exp(-1)) * Ah^p * c
                    #-----------------------------------
                    #** S_vals = c * (A_vals^p) * (1.0 - exp(-A_vals/Ah))
                    
                    #-----------------------------------
                    # Try to capture convex hillslopes
                    # with a second power-law curve.
                    #-------------------------------------------------------------
                    # Can force continuity, but can't match derivatives or
                    # get p2=p1. This can be seen with a figure.  We expect
                    # 0 < p2 < 1, so we'll just fix p2 and compute c2 from cont.
                    #-------------------------------------------------------------
                    # ww = where(A_vals < Ah)
                    # nww = ww[0].size
                    # if (nww != 0):
                    #     Smax = c * (Ah**p)
                    #** p2 = 0.1
                    #** p2 = 0.5
                    #   p2 = 0.8
                    #** p2 = 1
                    #** p2 = 2
                    #** p2 = 4
                    #   c2   = Smax / Ah**p2
                    #   S_vals[ww] = c2 * (A_vals[ww]**p2)
                    
                    #------------------------------------------
                    # Update the new, smooth elevation values
                    #---------------------------------------------------------
                    # Note: Since D8 areas always increase downstream, the
                    # slopes computed from Flint's law must always decrease.
                    #---------------------------------------------------------
                    ds_vals = self.d8.ds.flat[ ready_IDs ]  # [meters]
                    dz_vals = S_vals * ds_vals  # [meters]
                    self.DEM.flat[ ready_IDs ] = (p_elev[wp] + dz_vals)
                    STILL_ACTIVE = True

                    #-------------------------------------
                    # Compute the flow distances to edge
                    #-------------------------------------
                    self.flow_dist.flat[ ready_IDs ] = (p_dist[wp] + ds_vals)
                    
                #------------------------
                # Are we finished yet ?
                #------------------------
                DONE = (n_assigned == n_IDs)
            
            if (DONE or not(STILL_ACTIVE)):  break

        #--------------------------
        # Compute DEM min and max
        #--------------------------
        self.zmin = np.nanmin(self.DEM)
        self.zmax = np.nanmax(self.DEM)

        #--------------------------------------------------
        # Adjust the values by a distance-weighted amount
        # so that max of new DEM will be same as old
        #-------------------------------------------------
        # wmax  = where( self.DEM == self.zmax )
        # dmax  = self.flow_dist[ (wmax[0][0], wmax[1][0]) ]
        # del_z = (self.flow_dist / dmax)*(self.zmax - self.z0max)
        # self.DEM = self.DEM - del_z
        
        #-------------------------------------------------
        # Scale the values by a distance-weighted factor
        # so that max of new DEM will be same as old
        #-------------------------------------------------
        # factor = (1 - (self.flow_dist / dmax)) * 
        # self.DEM = self.DEM * factor
        
        #----------------------
        # Print final message
        #----------------------
        ## if (self.REPORT):
##        print 'Finished with new DEM. '
##        print ' '
        print('Number of iterations = ' + str(n_reps))
        print('Min/Max of orig. DEM = ' + \
                     str(self.z0min) + ', ' + str(self.z0max))
        print('Min/Max of new  DEM  = ' + \
                     str(self.zmin) + ', ' + str(self.zmax))
        print(' ')
        
    #   update_DEM()            
    #-------------------------------------------------------------------------  
    def read_profile_data(self, n_header=None):

        #--------------------------------------------------------
        # Notes: This routine gets pixel IDs for a main channel
        #        streamline from profile_file and uses them to
        #        get elevations, areas and pixel-to-pixel flow
        #        lengths along the main channel for use by the
        #        best_slope_area_curve_fit routine.
        #--------------------------------------------------------  
        if (n_header == None):
            n_header = np.int16(6)

        #------------------------------
        # Count lines in profile file
        #------------------------------
        n_lines = file_utils.count_lines( self.profile_file, SILENT=True )
        n_lines = (n_lines - n_header)
        #-------------------------------
        dist = np.zeros([n_lines], dtype='Float64')  ## 1/16/12
        elev = np.zeros([n_lines], dtype='Float64')  ## 1/16/12
        cols = np.zeros([n_lines], dtype='Int32')
        rows = np.zeros([n_lines], dtype='Int32')

        #-----------------------------
        # Open file to read IDs and
        # skip over the header lines
        #-----------------------------
        file_unit = open(self.profile_file, 'r')
        cfg_files.skip_header( file_unit, n_lines=n_header )

        #----------------------------------
        # Read the column and row vectors
        #-----------------------------------------------------
        # Profile file has: distance, elevation, column, row
        #-----------------------------------------------------
        dtype_list = ['float64','float64','int32', 'int32']
        for k in range(n_lines):
            var_list = cfg_files.read_list( file_unit, dtype_list=dtype_list )
            dist[k] = var_list[0]  ## 1/16/12
            elev[k] = var_list[1]  ## 1/16/12
            cols[k] = var_list[2]
            rows[k] = var_list[3]

        #---------------------
        # Close profile_file
        #---------------------
        file_unit.close()

        #--------------------------------------------
        # Read the DEM, area and D8 flow code grids
        #-------------------------------------------------
        # 1/16/12. Should we add area_file and flow_file
        # to the CFG file?  It already has DEM_file.
        #-------------------------------------------------
        dp = (self.in_directory + self.site_prefix)
        DEM_file  = dp + '_DEM.rtg'   ## 1/16/12
        area_file = dp + '_area.rtg'  ## 1/16/12
        #--------------------------------------------
        info  = self.rti
        DEM   = rtg_files.read_grid( self.DEM_file,  info, RTG_type=info.data_type )
        areas = rtg_files.read_grid( area_file, info, RTG_type='FLOAT' )
        ## ds = rtg_files.read_grid( ds_file, info, RTG_type='FLOAT' )

        ######### Done by read_input_files() ??
        
        #---------------------------------------
        # Only used for Flow_Lengths function.
        #---------------------------------------
        # flow_file = self.site_prefix + '_flow.rtg'  ## 1/16/12        
        # codes = rtg_files.read_grid( flow_file, info, RTG_type='BYTE' )
        #-----------------------------------------------------
        # Compute the along-channel flow lengths (ds vector)
        #-----------------------------------------------------
        # ds = Flow_Lengths(codes, RTI_file, METERS=True, DOUBLE=True)  ########

        #------------------------------------------------------
        # Construct ds vector from distances in profile_file.
        # First distance is always zero.
        # Also need to convert from km to meters.
        # Size of "diffs" is one less than size of "dist".
        #------------------------------------------------------
        diffs          = np.diff( dist )
        # print 'size(dist)  =', dist.size
        # print 'size(diffs) =', diffs.size
        ds_profile      = np.zeros( dist.size, dtype='Float64' )
        ds_profile[:-1] = diffs
        ds_profile[-1]  = diffs[-2]  ######################  NOT STRICTLY CORRECT
        ds_profile      = ds_profile * 1000.0   # [meters]
    
        #------------------------------------------
        # Construct calendar-style streamline IDs
        #------------------------------------------
        ncols = np.size(DEM, 1)
        IDs   = (ncols * rows) + cols

        #-------------------------------------
        # Get the profile elevations & areas
        #-------------------------------------
        ### z_profile  = elev     # (Use this instead ?? 1/16/12)
        z_profile  = DEM.flat[ IDs ]     # [meters]
        A_profile  = areas.flat[ IDs ]   # [km^2]
        # ds_profile = ds.flat[ IDs ]    # [meters]

        #-------------------------------------
        # Reverse the vectors so that values
        # start at outlet and work upstream
        #-----------------------------------------------------------
        # Must use FLIPUD(x) vs. ROT90(x,-2) to reverse 1D arrays.
        #-----------------------------------------------------------
        self.A_profile  = np.flipud( A_profile )
        self.z_profile  = np.flipud( z_profile )
        self.ds_profile = np.flipud( ds_profile )

    #   read_profile_data()
    #-------------------------------------------------------------------------
    def find_best_fit_c_and_p(self, weights=None, REPORT=True):
                              ## itmax=None, tol=None ):

        #------------------------------------------------------------
        # Notes: These notes are for the original IDL version.
        #
        #        This function uses IDL's CURVEFIT function and the
        #        procedure slope_area_curve (above) to find the
        #        best-fit parameters for fitting the data vectors
        #        A and z.

        #        x and y can have as few as 3 unique points, but
        #        must contain 4 elements each to avoid an error
        #        from IDL.  The 3rd value can simply be repeated.
        #        Initial guesses are required for all of the power
        #        curve parameters (a,c,p) and the choice of these
        #        has a big impact on whether CURVEFIT converges to
        #        a solution.  Some experimenting is necessary but
        #        the initial guess for p must be a large negative
        #        number like -10 if p is expected to be negative
        #        and a small positive number like 0.001 if p is
        #        expected to be positive ??

        #        The array of flags, fitvars, determines which
        #        parameters are fixed and which ones to find, we
        #        don't need to find z0, but we need to pass it.
        #------------------------------------------------------------
        A  = self.A_profile
        z  = self.z_profile
        ds = self.ds_profile

        #---------------------------------------------
        # Set weights for the curve fitting function
        #---------------------------------------------
        if (weights == None):
            #-----------------------------------------------
            # Use equal weights; gives smaller stderr
            # but further from computed p value. A
            # leading constant didn't matter for CURVEFIT.
            #-----------------------------------------------
            # weights = np.ones( A.size )

            #----------------------------------------------
            # Poisson statistical weighting, gives a
            # smaller stderr, but further from computed p
            #----------------------------------------------
            # weights = 1 / z

            #------------------------------------------------
            # Weight by contributing areas: improved fit.
            # Since c and p values are used for entire DEM,
            # and since the number of streamlines that pass
            # through a given pixel is proportional to the
            # contributing area, A, this makes some sense.
            #------------------------------------------------
            weights = A
            # weights = (A ** 1.1)
            # weights = np.sqrt(A)  ;(good compromise ?)
            # weights = (A ** 0.75)

            #---------------------------------------------
            # Combination of previous two methods, gives
            # worst stderr but closest to computed p.
            # Note that xdata=A, ydata=z in curve fit.
            #---------------------------------------------
            # weights = (A / z)
            
        w0  = where(weights == 0)
        nw0 = w0[0].size
        if (nw0 != 0):
            weights[w0] = np.float64(1)

        #------------------------------------------
        # Points used to generate initial guesses
        #------------------------------------------
        z0 = z[0]
        # z1 = z[1]
        z2 = z[-1]
        
        #---------------------------------------------
        # Need initial guesses with good convergence
        # properties; not necessarily close to value
        # (These worked well for IDL's CURVEFIT.)
        # Can't use p0 since keyword to curve_fit().
        #---------------------------------------------
        pg = np.float64( -0.5 )
        cg = (z2 - z0) / np.sum(np.float64(ds * (A ** pg)))

        #-------------------------------------------------------------
        # Define fitting function needed by scipy.optimize.curve_fit.
        # First argument is only allowed independent variable, while
        # remaining arguments are the fitting parameters.  Note that
        # ds (a vector) and z0 are treated as known values; the
        # curve_fit() function does not allow them as arguments.
        # Recall that S = c * A^p, and z = z0 + cumsum(ds * S).
        # We also want the first element of the estimate to be z0,
        # so we prepend 0 to dz.
        #
        # It looks like the first argument needs to be a scalar
        # in order for this to find best fit for both c and p.
        #-------------------------------------------------------------
        def fit_function(AA, cc, pp):
            dz = cc * np.float64( ds * ( AA ** pp ) )
            dz = np.concatenate(( [0], dz[:-1] ))
            return z0 + np.cumsum( dz )
        
        #--------------------------------------------------
        # Define "L2_error" function, also called E(c,p).
        #--------------------------------------------------
        def L2_error( params ):  ###, *weights ):
            cc = params[0]
            pp = params[1]
            nz = z.size
            dz = cc * np.float64( ds * ( A ** pp ) )
            dz = np.concatenate(( [0], dz[:-1] ))
            zf = z0 + np.cumsum(dz)
            #------------------------------------------------
            # Experiment:  Weighting by contributing area.
            # This gives a lower p-value, but seems to give
            # better results when applied to entire DEM.
            #------------------------------------------------
            weights = A
            return np.sqrt( np.sum( weights*(z - zf)**2 ) / nz)
            # if (weights == None):
            #     return np.sqrt( np.sum( (z - zf)**2 ) / nz)
            # else: 
            #     return np.sqrt( np.sum( weights*(z - zf)**2 ) / nz)
  
        #----------------------------------------------------
        # Define "Fk(p)" function used by c1(p) and c2(p). 
        #----------------------------------------------------
        def Fk_function( k, p ):
            if (k == 0): return 0.0
            A_vals  = A[1: k+1]
            ds_vals = ds[1: k+1]
            return np.sum( (A_vals**p) * ds_vals )

        #----------------------------------------------------
        # Define "Fk(p)" function used by c1(p) and c2(p). 
        #----------------------------------------------------
        def Fkd_function( k, p ):
            if (k == 0): return 0.0
            A_vals  = A[1: k+1]
            ds_vals = ds[1: k+1]
            return np.sum( (A_vals**p) * np.log(A_vals) * ds_vals )

        #----------------------------------------------------
        # Define "c1(p)" function from d/dp[ E(c,p) ] = 0.
        #----------------------------------------------------
        def c1_function( p ):
            nz = z.size
            Fk_vals  = np.zeros( nz, dtype='float64' )
            Fkd_vals = np.zeros( nz, dtype='float64' )
            for k in range( nz ):
                Fk_vals[ k ]  = Fk_function( k, p )
                Fkd_vals[ k ] = Fkd_function( k, p ) 
            top = np.sum( (z - z0) * Fkd_vals )
            bot = np.sum( Fk_vals * Fkd_vals )
            return (top / bot)
 
        #----------------------------------------------------
        # Define "c2(p)" function from d/dc[ E(c,p) ] = 0.
        #----------------------------------------------------
        def c2_function( p ):
            nz = z.size
            Fk_vals  = np.zeros( nz, dtype='float64' )
            Fkd_vals = np.zeros( nz, dtype='float64' )
            for k in range( nz ):
                Fk_vals[ k ]  = Fk_function( k, p )
                Fkd_vals[ k ] = Fkd_function( k, p ) 
            top = np.sum( (z - z0) * Fk_vals )
            bot = np.sum( Fk_vals ** 2)
            return (top / bot)

        #-------------------------------------------------
        # Define "c_diff(p)" function  (for root finder)
        # Best c and p should be where c1(p) = c2(p).
        #-------------------------------------------------
        def c_diff( p ):
            return ( c1_function(p) - c2_function(p) )

        #-------------------------------
        # Define "c_diff2(p)" function 
        #-------------------------------
        def c_diff2( p ):
            return ( c1_function(p) - c2_function(p) )**2

        #---------------------------------------------------------------
        # Use scipy.optimize.fmin() to find best-fit parameters
        # by finding parameters that minimize the L2 error.
        # This uses the Nelder-Mead downhill simplex algorithm.
        #---------------------------------------------------------------
        # See: http://docs.scipy.org/doc/scipy/reference/optimize.html
        #---------------------------------------------------------------
        # If (disp=True), convergence messages are printed.
        # If (retall=True), best_params contains a list of solutions. 
        #-------------------------------------------------------------
        xtol    = 1e-12   # (tolerance in best parameters)
        maxiter = 300     # (max number of iterations)
        best_guesses = np.array((cg, pg))    # (an nd_array)

        #-----------------------------------------------------------
        # Each of these methods works, with very similar results,
        # including c, p, maxerr, E(c,p), c_1(p) and c_2(p).
        #-----------------------------------------------------------
        # Note that L2_error() now uses "weights".  It was
        # found previously with IDL's CURVEFIT that best results
        # were obtained with (weights = A), which causes downstream
        # points/pixels to have greater influence. This makes some
        # sense since the number of distinct streamlines/profiles
        # that pass through a given pixel is proportional to its
        # contributing area.  It also causes the max in the new
        # DEMs to have a much more reasonable value, even though
        # the fit to the main channel profile used to find c and
        # p has a greater error. 
        #-----------------------------------------------------------
        results = scipy.optimize.fmin( L2_error, best_guesses,
                    xtol=xtol, maxiter=maxiter,
                    disp=True, retall=True )
        #------------------------------------------------------------------
        # results = scipy.optimize.fmin_powell( L2_error, best_guesses,
        #             xtol=xtol, maxiter=maxiter,
        #             disp=True, retall=True )

        #------------------------------------------------------------
        # This experimental method also worked, but resulted in
        # larger maxerr and stderr, even though c1(p) and c2(p)
        # were closer to equal.  Note that func(a) and func(b) must
        # have opposite signs and they did for KY_Sub when a=-1.0,
        # b=1.0, as shown.  Also took longer to run.
        #------------------------------------------------------------
        # best_p = scipy.optimize.brentq( c_diff, -1.0, 1.0,  
        #             xtol=xtol, maxiter=maxiter, disp=True )
        # best_c = c1_function( best_p )
        # best_pair = np.array( [best_c, best_p] )
        # results = ( best_pair, best_pair )

        #-----------------------------------------------------------
        # Experimental method.  Didn't work with c_diff2 above.
        #-----------------------------------------------------------
        # p_guess = np.array( pg )
        # results = scipy.optimize.fmin( c_diff2, p_guess, 
        #             xtol=xtol, maxiter=maxiter, disp=True, retall=True )
        # best_p = results[0]
        # best_c = c1_function( best_p )
        # best_pair = np.array( best_c, best_p )
        # results[0] = best_pair

        #-----------------------------------------------------------
        # Didn't work with the default settings, as shown here.
        # DISP keyword not suppported in SciPy 0.9.
        #-----------------------------------------------------------
        # best_params = scipy.optimize.anneal( L2_error, best_guesses,
        #                 feps=xtol, maxiter=maxiter )
        # results = [ best_params, best_params ]

        #--------------------------------------------------------------------
        # This method requires a function for the derivative, "fprime"
        #--------------------------------------------------------------------
        # results = scipy.optimize.fmin_ncg( L2_error, best_guesses,
        #             fprime= ????????,
        #             avextol=xtol, maxiter=maxiter, disp=True, retall=True )

        #--------------------------------------------------------------------
        # These methods didn't give similar results; p did not change from
        # its initial value.  Also, they don't allow the XTOL keyword,
        # but tried the GTOL keyword.
        #--------------------------------------------------------------------
        # results = scipy.optimize.fmin_cg( L2_error, best_guesses,
        #             gtol=xtol, maxiter=maxiter, disp=True, retall=True )
        #--------------------------------------------------------------------
        # results = scipy.optimize.fmin_bfgs( L2_error, best_guesses,
        #             gtol=xtol, maxiter=maxiter, disp=True, retall=True )
        #--------------------------------------------------------------------
        print(' ')      # (after converence message)
        best_params = results[0]
        pair_list   = results[1]
        self.c      = best_params[0]
        self.p      = best_params[1]

        if (REPORT):
            print('List of (c,p) pairs:')
            for pair in pair_list:
                print('    (c,p) =', pair)
            print(' ')
 
        # Note: minimize() is not available in SciPy 0.9.
        # best_params, info = scipy.optimize.minimize( L2_error, best_guesses,
        #                                              method='Nelder-Mead')
   
        #-------------------------------------------------------------
        # Use scipy.optimize.curve_fit() to find best-fit parameters.
        # It uses nonlinear least squares to fit a function to data.
        #-------------------------------------------------------------
        # http://docs.scipy.org/doc/scipy/reference/generated/
        #        scipy.optimize.curve_fit.html
        # Uses the Levenburg-Marquardt algorithm implemented as:
        #      scipy.optimize.leastsq()
        # Additional keyword arguments are passed directly to that
        # algorithm. See help(scipy.optimize.leastsq) for more info
        # on keywords such as:
        #     maxfev: max number of iterations
        #     ftol:   Relative error desired in the sum of squares.
        #     xtol:   Relative error desired in the approximate solution.
        #     ier:    An integer information flag. (returned)
        #     mesg:   An error message string.     (returned)
        #
        # popt, pcov = scipy.optimize.curve_fit(f, xdata, ydata,
        #                             p0=None, sigma=None, **kw)
        #
        # Keywords not expected by curve_fit() are passed directly
        # to the underlying leastsq() function.
        #-------------------------------------------------------------
        maxfev = 300                    # (used for IDL's CURVEFIT)
        xtol   = np.float64( 1e-10 )
        # xtol   = np.float64( 1e-20 ) # (used for IDL's CURVEFIT)
        # kwargs = { "maxfev":maxfev, "xtol":xtol }  # (Works, but not needed.)
        # I don't know how to get values returned in keywords.
        # This doesn't work: kwargs = { "ier":None, "mesg":None }
        # This doesn't work: kwargs = { "ier":0, "mesg":'NONE' }
        # best_guesses = [cg, pg]   # (a list)
        # best_guesses = (cg, pg)   # (a tuple)

        # best_guesses = np.array((cg, pg))  # (an nd_array)
        # xdata = A
        # ydata = z
        # best_params, best_cov = scipy.optimize.curve_fit( fit_function,
        #                                        xdata, ydata,
        #                                        p0=best_guesses, ## p0=None,
        #                                        ## sigma=weights, ## sigma=None,
        #                                        maxfev=maxfev,
        #                                        xtol=xtol )
        #                                        ## **kwargs )
        # self.c = best_params[0]
        # self.p = best_params[1]

        # ier  = kwargs['ier']
        # mesg = kwargs['mesg']
        # print 'ier =', ier
        # print 'mesg =', mesg
        # ier  = 1 
        # mesg = 'NOT_SET'

        #--------------------------
        # Compute error estimates
        #--------------------------
        nz     = z.size
        zfit   = fit_function( A, self.c, self.p )
        maxerr = np.max( np.absolute( z - zfit ))
        stderr = np.sqrt( np.sum( (z - zfit)**2  )/ nz )
        # stderr = np.sqrt( np.sum( (z - zfit)**2  )/(nz - 1))  # (Bessel's correction?)

        #--------------------------------
        # Print comparison of zfit to z 
        #--------------------------------
        if (REPORT):
            for k in range( len(z) ):
                print('(z, zfit, diff) =', z[k], ',', zfit[k], ',', (z[k]-zfit[k]))
            print(' ')
            # print 'A =', A
            # print ' '
            # print 'ds =', ds
            # print ' '

        #---------------------------
        # Print an optional report
        #---------------------------
        if (REPORT):
            print('--------------------------------------')
            print(' Least squares curve fit to profile')
            print(' weighted by contributing area') 
            print('--------------------------------------')
            print('z(A)   = z0 + np.cumsum(dz(A))')
            print('dz(A)  = [0, ds * S(A)]')
            print('S(A)   = c * (A^p)')
            print('--------------------------------------')
            print('z0     = ' + str(z0))
            print('zmin,   zmax =', np.nanmin(z),  ',', np.nanmax(z))
            print('Amin,   Amax =', np.nanmin(A),  ',', np.nanmax(A))
            print('dsmin, dsmax =', np.nanmin(ds), ',', np.nanmax(ds))
            print('--------------------------------------')
            print('c0     = ' + str(cg))
            print('p0     = ' + str(pg))
            print('--------------------------------------')
            print('c      = ' + str(self.c))
            print('p      = ' + str(self.p))
            print('maxerr = ' + str(maxerr)) 
            print('stderr = ' + str(stderr))    # (same as E(c,p)
            print('--------------------------------------')
            print('E(c,p) = ' + str( L2_error( best_params ) ))
            print('c_1(p) = ' + str( c1_function( self.p ) ))
            print('c_2(p) = ' + str( c2_function( self.p ) ))
            print('--------------------------------------')
            print(' ')

            #--------------------------------
            # Print status of the curve fit
            #-----------------------------------------------------
            # IDL's CURVEFIT provided information about whether
            # the algorithm converged or not and the number of
            # iterations.  scipy.optimize.leastsq() provides
            # similar information in "ier" and "mesg".
            #-----------------------------------------------------
            # good_codes = [1,2,3,4]
            # if (ier not in good_codes):
            #     print 'Error: ' + mesg
            # else:
            #     print 'Message: ' + mesg
            # print ' '

        #---------------------------------------------------
        # Use IDL's CURVEFIT() to find best-fit parameters
        #-------------------------------------------------------------------
        # Result = CURVEFIT( X, Y, Weights, A [, Sigma] [, CHISQ=variable]
        #                    [, /DOUBLE] [, FITA=vector]
        #                    [, FUNCTION_NAME=string] [, /NODERIVATIVE]
        #                    [, ITER=variable] [, ITMAX=value]
        #                    [, STATUS={0 | 1 | 2}] [, TOL=value]
        #                    [, YERROR=variable] )
        #-------------------------------------------------------------------
        # This is how CURVEFIT would be used:
        #
        # params  = [c0, p0, z0]
        # fitvars = [1, 1, 0]
        # zfit = curvefit(A, z, weights, params, sigma, DOUBLE=True,
        #                 FUNCTION_NAME='IDL_fit_function', TOL=tol,
        #                 ITMAX=itmax, YERROR=stderr, FITA=fitvars,
        #                 STATUS=status, ITER=iter)
        # c    = params[0]  ; (these are passed back out)
        # p    = params[1]
        # zfit = IDL_fit_function( A, c, p )  # (z0 and ds via "data")
        # zfit = z0 + (c * cumsum(ds * A**p))
        
        # if (status == 0):    
        #     print 'Curve fit was successful!'
        #     print 'Number of iterations = ' + str(iter)
        # elif (status == 1):    
        #     print 'Curve fit failed. Chi square was '
        #     print 'increasing without bounds.'
        # elif (status == 2):    
        #     print 'Curve fit failed to converge after'
        #     print str(itmax) + ' iterations.'
        # else:
        #     raise RuntimeError('no match found for expression')
        # print ' '
        #---------------------------------------------------------------------
        
    #   find_best_fit_c_and_p()       
    #-------------------------------------------------------------------------
##    def IDL_fit_function(A, c, p):
##                         ### params, z, partials):
##
##        #-------------------------------------------------------------
##        # Notes:  For use with IDL's CURVEFIT() function
##        #
##        #         CUMULATIVE keyword to TOTAL gives partial sums and
##        #         returns a vector vs. a scalar.
##        #-------------------------------------------------------------
##        # NB!     z0 is the elevation of the parent pixel of the
##        #         outlet pixel.  It is not the elevation of the
##        #         outlet pixel and A is max (not zero) at the outlet
##        #         pixel.
##        #-------------------------------------------------------------
##        # NB!     Procedures called by IDL's CURVEFIT must conform
##        #         to strict rules.  This means that we have no way
##        #         to pass an additional vector argument like ds.
##        #         However, we can use a COMMON block, as done here.
##        #-------------------------------------------------------------   
##        ds = common_block.ds_profile
##        z0 = common_block.z0
##
##        z = z0 + (c * np.cumsum( float64(ds * A ** p) ))
##
##        return z
##
##        #------------------------------
##        # Compute partial derivatives
##        #---------------------------------
##        # n_params() refers to number of
##        # arguments to this procedure.
##        #---------------------------------
##    ##    if (n_params >= 4):    
##    ##        dz_dc = np.cumsum(double(ds * A ** p))
##    ##        dz_dp = c * np.cumsum(double(ds * log(A) * A ** p))
##    ##        nA = np.size(A)
##    ##        dz_dz0 = zeros([nA], dtype='Float64') + 1.0
##    ##        partials = array([array([dz_dc]), array([dz_dp]), array([dz_dz0])])
##    ##    
##    ##    return (A, params, z, partials)
##
##    #   IDL_fit_function() 
    #-------------------------------------------------------------------------
    def open_input_files(self):

        pass

    #   open_input_files()
    #------------------------------------------------------------------------- 
    def read_input_files(self):

        #----------------------------------------
        # Get name of the info file and read it
        #----------------------------------------
        info = self.rti
        # info = rti_files.read_info( self.DEM_file )

        #-----------------------
        # Read the initial DEM
        #-----------------------
        self.z0  = rtg_files.read_grid( self.DEM_file,  info,
                                        RTG_type=info.data_type )
        self.DEM = self.z0.copy()

        #---------------------------------
        # Store original DEM min and max
        #---------------------------------
        self.z0min = np.nanmin( self.z0 )
        self.z0max = np.nanmax( self.z0 )
        
        #------------------------------------------------
        # Could read these, but now we use d8_global.py
        # to compute them to allow evolution.
        #------------------------------------------------
        # self.areas = rtg_files.read_grid( self.area_file, info, RTG_type='FLOAT' )
        # self.codes = rtg_files.read_grid( self.flow_file, info, RTG_type='BYTE' )
        
    #   read_input_files()
    #-------------------------------------------------------------------------  
    def close_input_files(self):

        pass

    #   close_input_files()
    #-------------------------------------------------------------------------
    def update_outfile_names(self):

        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.z_gs_file  = (self.out_directory + self.z_gs_file)
        self.D8_gs_file = (self.out_directory + self.D8_gs_file)
        self.S_gs_file  = (self.out_directory + self.S_gs_file) 
        self.A_gs_file  = (self.out_directory + self.A_gs_file) 
        #---------------------------------------------------------
        self.z_ts_file  = (self.out_directory + self.z_ts_file)
        self.D8_ts_file = (self.out_directory + self.D8_ts_file) 
        self.S_ts_file  = (self.out_directory + self.S_ts_file) 
        self.A_ts_file  = (self.out_directory + self.A_ts_file)
        
##        self.new_DEM_file    = (self.out_directory + self.new_DEM_file)
##        self.new_rawDEM_file = (self.out_directory + self.new_rawDEM_file) 
##        self.new_slope_file  = (self.out_directory + self.new_slope_file) 
##        self.new_flow_file   = (self.out_directory + self.new_flow_file) 

    #   update_outfile_names()     
    #-------------------------------------------------------------------  
    def open_output_files(self):

        model_output.check_netcdf()    # (test import and info message)
        self.update_outfile_names()

        #--------------------------------------
        # Open new files to write grid stacks
        #--------------------------------------
        # open_new_gs_file() has a "dtype" keyword that defaults
        # to "float32".  Flow codes have dtype = "uint8".
        #-----------------------------------------------------------
        if (self.SAVE_Z_GRIDS):   
            model_output.open_new_gs_file( self, self.z_gs_file, self.rti,
                                           var_name='z',
                                           long_name='elevation',
                                           units_name='m')
            
        if (self.SAVE_D8_GRIDS):    
            model_output.open_new_gs_file( self, self.D8_gs_file, self.rti,
                                           dtype='uint8',
                                           var_name='D8',
                                           long_name='D8 flow direction codes',
                                           units_name='none')
            
        if (self.SAVE_S_GRIDS):    
            model_output.open_new_gs_file( self, self.S_gs_file, self.rti,
                                           var_name='S',
                                           long_name='surface slope',
                                           units_name='m/m')
        
        if (self.SAVE_A_GRIDS):    
            model_output.open_new_gs_file( self, self.A_gs_file, self.rti,
                                           var_name='A',
                                           long_name='D8 contributing area',
                                           units_name='km^2')

            
        #--------------------------------------
        # Open new files to write time series
        #--------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_Z_PIXELS):  
            model_output.open_new_ts_file( self, self.z_ts_file, IDs,
                                           var_name='z',
                                           long_name='elevation',
                                           units_name='m')
                                          
        if (self.SAVE_D8_PIXELS):
            model_output.open_new_ts_file( self, self.D8_ts_file, IDs,
                                           dtype='uint8',
                                           var_name='D8',
                                           long_name='D8 flow direction codes',
                                           units_name='none')
                                          
        if (self.SAVE_S_PIXELS):    
            model_output.open_new_ts_file( self, self.S_ts_file, IDs,
                                           var_name='S',
                                           long_name='surface slope',
                                           units_name='m/m')
            
        if (self.SAVE_A_PIXELS):    
            model_output.open_new_ts_file( self, self.A_ts_file, IDs,
                                           var_name='A',
                                           long_name='D8 contributing area',
                                           units_name='km^2')

        #-------------------------------------
        # Save FLOAT version of original DEM
        # as the rawDEM for the new DEM
        #-------------------------------------
##        if (self.rti.SWAP_ENDIAN):
##            array(float32(self.z0), copy=0).byteswap(True)
##        new_rawDEM_unit = open( self.new_rawDEM_file, 'wb' )
##        float32(self.z0).tofile( new_rawDEM_unit )
##        new_rawDEM_unit.close()

##        self.new_DEM_unit    = open(self.new_DEM_file,    'wb')
##        self.new_slope_unit  = open(self.new_slope_file,  'wb')
##        self.new_rawDEM_unit = open(self.new_rawDEM_file, 'wb')
##        self.new_flow_unit   = open(self.new_flow_file,   'wb')

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

##        SWAP_ENDIAN = self.rti.SWAP_ENDIAN
##        
##        #-----------------------
##        # Save new DEM in file
##        #-----------------------
##        if (SWAP_ENDIAN):
##            array(float32(self.DEM), copy=0).byteswap(True)
##        float32(self.DEM).tofile( self.new_DEM_unit )
##        #-----------------------------
##        # Create RTI file for new DEM
##        #------------------------------
##        info = self.rti
##        info.data_type = 'FLOAT'
##        #info.DEM_file  = str(self.new_DEM_unit.name)
##        rti_files.write_info( self.new_RTI_file, info )
##        
##        #--------------------------------------
##        # Create and save new slope grid file
##        #-----------------------------------------
##        # Subpixel sinuosity, if any, is applied
##        # later in Route_Flow.  Both ds and the
##        # pID_grid were computed above.
##        #-----------------------------------------
##        ## slopes = (self.new_DEM - self.new_DEM[pID_grid]) / ds
##        if (SWAP_ENDIAN):
##            array(float32(self.S), copy=0).byteswap(True)
##        float32(self.S).tofile( self.new_slope_unit )
## 
##        #------------------------------------
##        # Save D8 flow grid of original DEM
##        # as the flow grid of the new DEM
##        #-----------------------------------------
##        # Check that flow grid hasn't changed ??      ;**********************
##        #-----------------------------------------
##        if (SWAP_ENDIAN):
##            array(self.d8.d8_grid, copy=0).byteswap(True)
##        self.d8.d8_grid.tofile( self.new_flow_unit )
##        # self.d8.d8_codes.tofile( self.new_flow_unit )
        
    #   write_output_files()
    #-------------------------------------------------------------------  
    def close_output_files(self):

##        self.new_DEM_unit.close()
##        self.new_slope_unit.close()
##        self.new_rawDEM_unit.close()
##        self.new_flow_unit.close()

        if (self.SAVE_Z_GRIDS):  model_output.close_gs_file( self, 'z')   
        if (self.SAVE_D8_GRIDS): model_output.close_gs_file( self, 'D8')  
        if (self.SAVE_S_GRIDS):  model_output.close_gs_file( self, 'S')   
        if (self.SAVE_A_GRIDS):  model_output.close_gs_file( self, 'A')
        #---------------------------------------------------------------------
        if (self.SAVE_Z_PIXELS):  model_output.close_ts_file( self, 'z')   
        if (self.SAVE_D8_PIXELS): model_output.close_ts_file( self, 'D8')    
        if (self.SAVE_S_PIXELS):  model_output.close_ts_file( self, 'S')    
        if (self.SAVE_A_PIXELS):  model_output.close_ts_file( self, 'A')
        
    #   close_output_files()
    #-------------------------------------------------------------------  
    def save_grids(self):
        
        #-----------------------------------
        # Save grid stack to a netCDF file
        #---------------------------------------------
        # Note that add_grid() methods will convert
        # var from scalar to grid now, if necessary.
        #---------------------------------------------        
        if (self.SAVE_Z_GRIDS):
            if (self.time_index == 0):
                #--------------------------------------
                # Save original DEM as the first grid
                #--------------------------------------
                model_output.add_grid( self, self.z0, 'z', 0.0 )
            model_output.add_grid( self, self.DEM, 'z', self.time_min )
            
        if (self.SAVE_D8_GRIDS):
            model_output.add_grid( self, self.d8.d8_grid, 'D8', self.time_min )
            
        if (self.SAVE_S_GRIDS):
            model_output.add_grid( self, self.S, 'S', self.time_min )

        if (self.SAVE_A_GRIDS):
            model_output.add_grid( self, self.d8.A, 'A', self.time_min )     

    #   save_grids()
    #-------------------------------------------------------------------  
    def save_pixel_values(self):   ##### save_time_series_data(self)  #######
        
        IDs  = self.outlet_IDs
        time = self.time_min       #####

        #-------------
        # New method
        #-------------
        if (self.SAVE_Z_PIXELS):
            model_output.add_values_at_IDs( self, time, self.DEM, 'z', IDs )
                    
        if (self.SAVE_D8_PIXELS):
            model_output.add_values_at_IDs( self, time, self.d8.d8_grid, 'D8', IDs )
            
        if (self.SAVE_S_PIXELS):
            model_output.add_values_at_IDs( self, time, self.S, 'S', IDs )
            
        if (self.SAVE_A_PIXELS):
            model_output.add_values_at_IDs( self, time, self.d8.A, 'A', IDs )
        
    #   save_pixel_values()
    #------------------------------------------------------------------- 

