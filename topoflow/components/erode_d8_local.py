
############################################################
## 2/28/12.  Changed "test1" in update_neighbors().
##           Now using update_dt_grid6().
##           Wrote check_zone_info() method to run various
##           tests on values at zone_IDs.  For example, it
##           checks that A=0, Qs=0, Q=0 when d8_code=0.
##
## 2/29/12.  No longer putting cells and the heap that have
##           dt >= dt_limit in update_T_next().  Also, in
##           update_active_ID(), some cells have (d8_code=0)
##           when drawn from the heap. And update_area_grid()
##           is still reporting a backflow problem.
##
## 1/13/12.  Check paths in all unit tests.
## 
## In CSDMS_base.update_time(), and elsewhere, string
##    comparisons may be reducing performance.  Removed
##    some others on 11/20/11.
##
## In d8_local.update_area_grid(), needed to use ".copy()"
## in several places; might be faster now. (11/20/11)
##
## To profile the code, do this at the Python prompt:
##    >>> import cProfile   # (a built-in profiler module)
##    >>> cProfile.run('unit_test()')
##
## Note: There is a slim chance of cases where max_reps = 
##       2 * max(nx, ny) used in update_area_grid() in
##       d8_local.py is not large enough.
##
############################################################

## Copyright (c) 2001-2013, Scott D. Peckham
##
#-------------------------------------------------------------
##
## January 2009  (converted Erode-D8-Global from IDL)
## September, October, November 2009
## February, March 2010 (local time-stepping)
## August 2010 (local time-stepping based on DES algorithm.)
## Sept. to November 2010 (ongoing R&D, new CFL condition)
## August-September 2011 (more testing, get_neighbors() bug fix)
## October 2011 (update_fluxes2; changes to update_neighbors.)
## November 2011 (another bug in update_fluxes2(), etc.)
##
#-------------------------------------------------------------
##
## Feb. 2012  Changes to work on all"zone_IDs" at once,
##            which clarifies the logic and greatly reduces
##            slow calls to: setdiff1d(), in1d() and unique().
##
#-------------------------------------------------------------
## Things to test:
##
## + Check why AREA_GRID_WARNING is encountered within
##   d8_local.update_area_grid(), with "backflow".  This
##   loop may be costly (before it breaks out).
##   ====>  Needed for update_flow_grid() to update nbr_IDs
##          as well as IDs. Messages gone now. (3/1/12)
##
## + We should add a total mass balance test that compares
##   total initial mass, minus whatever mass has flowed off
##   of the grid, to the total current or final mass.  This
##   may require integrating all pixels up to current time.
## + Check performance of update_active_ID() over time as
##   the size of the heap grows. While loop runs longer to
##   find cells that have not been retracted.
##
## + Check how cells are put on and "removed" from the heap,
##   especially those with a flow code of zero, like pits.
##   ====> Looked into this on 3/1/12.

## + What should we do in update_neighbors() if dz_cap
##   and dz_target have *opposite* signs ?
##   ====> Added code to update them also. (10/6/11)
##
## + Is update_fluxes() updating fluxes for some cells that
##   it shouldn't be?  ====> See new update_fluxes2().
##
## + Should update_fluxes() be called before the call to
##   update_dz_dt_grid()in update() ? (Not clear from paper.)
##   ====>  Doesn't seem to make any difference. (9/15/11)
##   ====>  Does make a difference now. (11/11/11)
##
## + Should we use A_IDs or new_A_IDs in update_fluxes() ?
##   ====>  Using new_A_IDs causes many "stalling" errors.
##
## + We may be able to use a CFL_factor > 0.2 now, after
##   fixing a bug in get_neighbors().
##
## + Is it okay to call update_dt_grid() before calling
##   update_dz_target() due to how implemented ?
##
## + Do we need a new version of update_dz_target() ?
##
## + Do we still need to set dz_target to zero for pits
##   so that they are always considered "affected" ?
##
## + Right now, cells with (pID == 0) or (d8_code == 0)
##   are "deactivated", meaning that they can never be
##   put on the heap and scheduled for future update.
##   They have their dt_grid value set to dt_limit.
##   These cells are always updated when one of their
##   neighbors is updated, and this is the only way for
##   them to get updated.
##
## + Should we use a pID of -1 instead of 0, or just flow
##   code values of 0 to identify pits, edges, etc. ?
##   Iteration to find downstream cells (as now done in
##   update_area_grid()?) won't work if we use -1.
##
## + Should run the 2 unit tests in d8_local.py again.
##   =====> Updated and fixed them, also in d8_global.py.
##   =====> They agree with RT3 except with how flats are
##          linked.  (TREYNOR, KY_SUB, BEAVER)
##   =====> Make sure LINK_FLATS=1 in new CFG file.
##
#-------------------------------------------------------------
## NB!  LINK_FLATS is set to False in initialize_d8_vars().
##      We don't need to fill pits or link flats because
##      this model fills all depressions dynamically.
##      However, there is a pit-filling option flag called
##      FILL_PITS_IN_Z0, in erode_base.initialize_DEM().
#-------------------------------------------------------------
## NB!  The following is done every time update_dz_dt_grid()
##      is called.  This may be costly but may also be
##      needed if base_IDs are allowed on the heap.
## 
##      self.dz_dt.flat[ self.base_IDs ] = -self.BLR_mpyr
#-------------------------------------------------------------
## NB!  In general, base_IDs will be a subset of edge_IDs,
##      but for now they're assumed to be the same.
#-------------------------------------------------------------

## See: self.SAVE_DZ_GRIDS = False   # (only use dz_target now)
        
#-----------------------------------------------------------------------
#
#  class erosion_component   # (inherits from erode_base.py)
#
#      get_component_name()
#      get_attribute()              ## (10/27/11)
#      initialize_d8_vars()
#      get_interior_pit_IDs()       ## (uses not_edge_grid)
#      initialize_computed_vars()
#----------------------------------
#      initialize_not_edge_grid()   ## (2/28/12)
#      initialize_slope_grid()
#      initialize_Q_grid()
#      initialize_Qs_grid()
#      ### initialize_D_grid_old()
#      initialize_D_grid()          ## (11/15/11)
#      initialize_dz_dt_grid()
#      ### initialize_dt_grid4()    ## (10/11/10)
#      initialize_dt_grid5()        ## (11/1/10)
#      initialize_dz_target_grid2()
#      initialize_T_last_grid()
#      initialize_T_next_grid()
#      save_initial_grids()         ## (3/7/12, like save_final_grids())
#-------------------------------
#      initialize_zone_IDs()        ## (2/10/12)
#      check_zone_info()            ## (2/28/12)
#      reset_zone_IDs()             ## (2/28/12)
#      append_zone_IDs()            ## (2/28/12)
#      append_all_now_IDs()         ## (2/28/12)
#-------------------------------
#      update()
#      write_output_files()         ## (11/15/11)
#      finalize()
#-------------------------------
#      update_active_ID()           ## (4/16/10)
#      update_T_clock()
#      update_DEM()                 ## (2/22/10)
#      update_entire_DEM()          ## (10/21/10)
#-------------------------------
#      get_neighbors()              ## (Moved to d8_local.py)
#-------------------------------
#      update_neighbors()           ## (8/23/10)
#      ### update_fluxes()              ## (9/13/10, separate)
#      ### update_fluxes2()             ## (9/16/11)
#      ### update_fluxes3()             ## (11/11)
#      update_fluxes4()             ## (2/10/12)  ### For new method.
#      update_flux_cap_grid()       ## (4/23/10)
#      update_T_last_grid()         ## (9/3/10, separate)
#-------------------------------
#      update_d8_vars()             ## (in erode_base.py )
#      update_slope_grid()
#      update_Q_grid()
#      update_Qs_grid()
#      update_D_grid_old()
#      update_D_grid()              ## (11/14/11)
#      update_dz_dt_grid()          ## (4/15/10)
#      update_mins_and_maxes()      ## (1/25/12)
#      update_grid_mins_and_maxes() ## (3/4/12)
#-------------------------------
#      ### update_dt_grid4()            ## (10/11/10,  OBSOLETE ??)
#      ### update_dt_grid5()            ## (11/1/10)
#      update_dt_grid6()            ## (2/28/12) #######################
#      update_dz_target_grid2()     ## (9/7/10)
#      update_T_next_grid()         ## (4/14/10)
#-------------------------------
#      save_update_count_grid()     ## (11/18/10)
#      save_d8_grids()              ## (11/18/10)
#      save_final_grids()           ## (9/2/11)

#-----------------------------------------------------------------------

import heapq
import numpy as np
import os          # for os.chdir in get_neighbors_test().
import os.path
import sys
import time
#----------------------
# Could also do this.
#----------------------
# from numpy import where, logical_and, logical_or

from topoflow.components import d8_local
from topoflow.components import erode_base
from topoflow.utils      import rtg_files

#-------------------------------------------
# For use outside of the TopoFlow package.
#-------------------------------------------
##import d8_local
##import erode_base
##import rtg_files

#-----------------------------------------------------------------------
class erosion_component( erode_base.erosion_component ):

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Erode_D8_Local'

    #   get_component_name() 
    #-------------------------------------------------------------------
    # Update this to match erode_d8_global.py.
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        map = {'comp_name':          'ErodeLocal',
               'version':            '0.5',
               'model_name':         'Erode_D8_Local',
               'model_family':       'Erode',
               'cfg_template_file':  'Erode_Local.cfg.in',
               'cfg_extension':      '_erode_local.cfg',
               'cmt_var_prefix':     '/ErodeLocal/Input/Var/',
               'gui_xml_file':       '/home/csdms/cca/erode/0.5/src/share/cmt/gui/Erode_Local.xml',
               'dialog_title':       'LEM: Erode DES Method Parameters',
               'time_step_type':     'local',
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
    def initialize_d8_vars(self):
        
        #---------------------------------------------
        # Compute and store a variety of (static) D8
        # flow grid variables.  Embed structure into
        # the "erosion_base" component.
        #---------------------------------------------
        self.d8 = d8_local.d8_component()   ### NOTE: LOCAL HERE ###

        #################################################
        # (5/13/10)  Do next line here for now, until
        # the d8 cfg_file includes site prefix.
        # Same is done in GW_base.py.
        #################################################
        # (1/23/12) Note that d8_base.py now has a new
        # method called: set_default_config_vars()
        # that is used to intialize vars in cases
        # where there is no "*_d8_global.cfg" file.
        # It is called in d8.initialize().
        ################################################
        self.d8.in_directory    = self.in_directory   # (1/23/12)
        self.d8.site_prefix     = self.site_prefix
        self.d8.FILL_PITS_IN_Z0 = 0                   # (1/23/12)
        self.d8.A_units         = 'm^2'               # (1/23/12) May be needed.        

        self.d8.initialize( cfg_prefix=self.cfg_prefix,
                            SILENT=self.SILENT,
                            REPORT=self.REPORT )

        #---------------------------------------------------
        # (11/14/11) Moved this here from update_d8_vars()
        # but not tested yet because we always use "m^2".
        # Note that self.da always has units of meters and
        # pixel_area may be converted to "km^2" in
        # d8_local.update_area_grid().
        #---------------------------------------------------
        # Erode model needs A_units to be "m^2"
        #----------------------------------------
        if ('km' in self.d8.A_units.lower()):
            self.d8.A = self.d8.A * 1e6   # [km^2 -> m^2]
            self.d8.A_units == 'm^2'  # (should this be here?)

        #---------------------------------------------------
        # This overrides settings from D8_Local CFG file.
        # We don't need to link flats in the Erode models
        # because pits are filled dynamically.  Note that
        # S=0 for flats, so outfluxes Q and Qs will also
        # be zero and deposition will therefore occur.
        #---------------------------------------------------
        self.d8.LINK_FLATS = False
        self.d8.BREAK_TIES = True
        print('In erode_d8_local.initialize_d8_vars():')
        print('   LINK_FLATS =', self.d8.LINK_FLATS)
        print('   BREAK_TIES =', self.d8.BREAK_TIES)
        print(' ')
        
        #------------------------------------------
        # Need this here, too, for now.
        # Note: IDs = None means use entire grid.
        #------------------------------------------
        self.update_d8_vars(IDs=None, SILENT=self.SILENT,
                            REPORT=self.REPORT)
        ### self.save_d8_grids()  # (save the initial grids)
         
    #   initialize_d8_vars()
    #-------------------------------------------------------------------
    def get_interior_pit_IDs(self, FLATS=False,
                             SILENT=True, REPORT=False):

        pit_IDs = self.d8.pit_IDs

        #-------------------------------------------
        # If FLATS keyword set, include those also
        #-------------------------------------------
        if (FLATS and (self.d8.flat_IDs.size > 0)):
            pit_IDs = np.concatenate(( pit_IDs, self.d8.flat_IDs ))

        #-----------------------------------------------
        # Exclude pits on the edges: new way (2/29/12)
        #-----------------------------------------------
        w = np.where( self.not_edge_grid.flat[ pit_IDs ] )
        if (w[0].size > 0):
            return pit_IDs[ w ]
        else:
            return np.array([])

        #-------------------------------------
        # Exclude pits on the edges: old way
        #-------------------------------------
##        nx = self.nx
##        ny = self.ny
##        pit_rows, pit_cols = divmod( pit_IDs, nx )
##
####        print '##### BEFORE excluding edge pixels:'
####        print 'pit_rows =', pit_rows
####        print 'pit_cols =', pit_cols
####        print ' '
##        
##        w1 = np.where( np.logical_and( pit_rows > 0, pit_rows < (ny-1) ) )
##        if (w1[0].size > 0):
##            pit_rows = pit_rows[ w1 ]
##            pit_cols = pit_cols[ w1 ]
##        else:
##            return np.array([])
##        
##        w2 = np.where( np.logical_and( pit_cols > 0, pit_cols < (nx-1) ) )
##        if (w2[0].size > 0):
##            pit_rows = pit_rows[ w2 ]
##            pit_cols = pit_cols[ w2 ]
##        else:
##            return np.array([])
##            
####        print '##### AFTER excluding edge pixels:'
####        print 'pit_rows =', pit_rows
####        print 'pit_cols =', pit_cols
####        sys.exit()
##        
##        return (pit_rows * nx) + pit_cols
        
    #   get_interior_pit_IDs()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self, DOUBLE=True):

        #------------------------------------------------------
        # Notes: Since self.d8.IDs is initialized to all IDs,
        #        original grid values here get overwritten.
        #
        #        initialize_DEM() is in erode_base.py and is
        #        called by initialize() there.
        #------------------------------------------------------
        self.active_ID_count = 0
        self.max_active_ID_count = 0
        self.last_active_ID = np.int32(0)
        self.last_last_active_ID = np.int32(0)
        
        #-------------------------------------
        # Override setting in erode_base.py.
        # (11/19/10) Not used right now.
        #-------------------------------------
        self.dz_tolerance = 1e-15    

        #-----------------------------------------------------
        # (11/15/11) Added these for write_output_files().
        # Initial values cause first grid to be saved.
        # "time_step_type" is in erode_base.set_constants().
        #-----------------------------------------------------
        self.last_grid_time  = -(self.save_grid_dt   + 1)
        self.last_pixel_time = -(self.save_pixels_dt + 1)
        # self.last_model_time = -1.0  # (Don't use 0.)

        #------------------------------------------------
        # This determines whether all computed vars are
        # 4-byte floats or 8-byte floats (doubles).
        # Notice use of ".astype(dtype)" below.
        #------------------------------------------------
        # dtype of DEM is set by create_initial_DEM()
        # in erode_base.py to "float64". (2/14/12)
        #------------------------------------------------
        if (DOUBLE):
            dtype = 'Float64'
        else:
            dtype = 'Float32'

        #-----------------------
        # Initialize some vars
        #-----------------------
        self.vol_R = np.float32(0)  # (for mass balance)
        self.vol_U = np.float32(0)  # (for mass balance)
        self.vol_R = self.vol_R.astype(dtype)
        self.vol_U = self.vol_U.astype(dtype)
        
        self.dz_max_vec = np.zeros([self.n_steps], dtype=dtype)
        self.dz_max     = np.float32(-9999).astype(dtype)
        self.dz_min     = np.float32(9999).astype(dtype)

        #----------------------------------------------------
        # This is done by initialize_DEM() in erode_base.py
        #----------------------------------------------------
        ## self.DEM_min = np.nanmin( self.DEM )
        ## self.DEM_max = np.nanmax( self.DEM )
        
        #--------------------------------
        # Need these here, too, for now
        #--------------------------------
        self.update_R()
        self.update_R_integral()
        self.update_U()
        self.update_U_integral()
        
        #-------------------------------
        # Initialize the "state" grids
        #-------------------------------
        self.initialize_not_edge_grid()       # (2/28/12)
        self.initialize_slope_grid()
        self.initialize_Q_grid()
        self.initialize_Qs_grid()
        self.initialize_D_grid()              # (9/7/10)
        self.initialize_dz_dt_grid()
        ## self.initialize_dt_grid4()         # (10/18/10)
        self.initialize_dt_grid5()            # (11/5/10)
        self.initialize_dz_target_grid2()     # (9/7/10)
        self.initialize_T_last_grid()
        self.initialize_T_next_grid()
        self.initialize_zone_IDs()   # (2/10/12  ###################)
        
        self.flux_cap_grid = np.zeros([self.ny, self.nx], dtype=dtype)

        self.update_count = np.zeros([self.ny, self.nx], dtype='Int32')

        #------------------------------------
        # (3/7/12) Save grids for debugging
        #------------------------------------
        self.save_initial_grids()
        
        #---------------------------------------
        # Initialize IDs to all pixels in grid
        #---------------------------------------
        self.T_clock = np.float32(0).astype(dtype)

        #----------------------------------
        # Override settings from CFG file
        #----------------------------------
        self.SAVE_N_GRIDS  = False
        self.SAVE_N_PIXELS = False
        self.SAVE_DZ_GRIDS = False   # (only use dz_target now)
     
    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def initialize_not_edge_grid(self, SILENT=True, REPORT=False):

        #----------------------------------------------------------
        # Note: Added this on 2/28/12.  Edge (or boundary)
        #       cells should be in "wait_IDs" while all other
        #       cells with (d8_code == 0) should be in "now_IDs".
        #----------------------------------------------------------
        #       This is used by initialize_dz_target_grid2().
        #----------------------------------------------------------        
        nx = self.nx
        ny = self.ny
        self.not_edge_grid = np.ones([ny, nx], dtype='bool') # (True)
        self.not_edge_grid[0, :]    = False
        self.not_edge_grid[ny-1, :] = False
        self.not_edge_grid[:, 0]    = False
        self.not_edge_grid[:, nx-1] = False

        ## print 'not_edge_grid ='
        ## print self.not_edge_grid
        
    #   initialize_not_edge_grid()
    #-------------------------------------------------------------------
    def initialize_slope_grid(self, SILENT=True, REPORT=False):

        #----------------------------------------------------
        # Notes: Make sure that d8 component is initialized
        #        with the same directory, site_prefix, etc.
        #        as this component.  Otherwise, the "shape"
        #        of DEM and ds, A, etc. won't match.
        #----------------------------------------------------
        if not(SILENT):    
            print('Initializing slope grid...')
   
        #--------------------------------------------------
        # Compute slope (rise/run) toward D8 parent pixel
        #--------------------------------------------------       
        pID_grid = self.d8.parent_ID_grid
        pIDs     = divmod( pID_grid, self.nx )
        self.S   = ((self.DEM - self.DEM[pIDs]) / self.d8.ds)

##        w = np.where(self.DEM < self.DEM[pIDs])
##        n_bad = w[0].size
##        if (n_bad != 0):
##            print '   Number of uphill parents =', n_bad
            
        #---------------------------------------------------
        # Set slope to zero wherever the D8 flow direction
        # is undefined.  This makes Qs(Q,S) = 0.
        #---------------------------------------------------
        ## w  = np.where(self.d8.parent_ID_grid == 0)
        w  = np.where(self.d8.d8_grid == 0)
        nw = w[0].size
        if (nw != 0):    
            self.S[w] = 0

        #-------------------------
        # Check for bad S values
        #-------------------------
        w2  = np.where(self.S < 0)
        nw2 = w2[0].size
        if (nw2 != 0):    
            # self.S[w2] = 1e-8
            print('   Negative values found in slope grid.')
            print('   Found ' + str(nw2) + ' negative values.')
            
            #---------------------------
            # Return the cols and rows
            #---------------------------
            rows = w2[0]
            cols = w2[1]
            print('   cols =')
            print(cols)
            print('   rows =')
            print(rows)
            sys.exit()   #####################

        #------------------------------------
        # Save the min and max slope values
        #------------------------------------
        self.S_min = np.nanmin( self.S )
        self.S_max = np.nanmax( self.S )
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            S_str = str(self.S_min) + ', ' + str(self.S_max)
            print('    min(S), max(S) = ' + S_str + ' [m/m]')

    #   initialize_slope_grid()
    #-------------------------------------------------------------------
    def initialize_Q_grid(self, SILENT=True, REPORT=False):

        #--------------------------------------------------------------
        # NOTES:  Q = annual discharge [m^3 / yr] (*leaving* a cell)
        #         R = geomorphically effective rainrate [meters / yr]
        #             (unless p ne 0, then R = coefficient)
        #         A = contributing area [meters^2]
        #--------------------------------------------------------------
        # Note: Q is zero wherever A is zero, and A will be zero
        #       wherever the D8 flow code is zero.
        #--------------------------------------------------------------        
        if not(SILENT):    
            print('Initializing discharge grid...')
        
        if (self.p != 1):    
            self.Q = self.R * (self.d8.A ** self.p)
        else:    
            #----------------------
            # Don't convert units
            #----------------------
            #** R2 = R / self.secs_per_year   ;[m/yr -> m/s]
            self.Q = self.R * self.d8.A

        #--------------------------------------------------------
        # Test whether there are any pixels where flow code is
        # zero but area grid isn't.  This would be bad because
        # flow code is zero in pits, and A > 0 would lead to a
        # nonzero discharge (water and sed.) *out* of the cell.
        #--------------------------------------------------------
        # Tested for 100 timesteps on 2/22/10.
        #--------------------------------------------------------
##        w  = np.where(np.logical_and(self.d8.d8_grid == 0, self.d8.A != 0))
##        nw = w[0].size
##        if (nw != 0):
##            print 'WARNING: There are places where flow grid is'
##            print '         zero and area grid is nonzero.'
##            print ' '

        #----------------------------------------
        # Save the min and max discharge values
        #----------------------------------------
        self.Q_min = np.nanmin( self.Q )
        self.Q_max = np.nanmax( self.Q )
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            Q_str = str(self.Q_min) + ', ' + str(self.Q_max)
            print('    min(Q), max(Q) = ' + Q_str + ' [m^3/yr]')

    #   initialize_Q_grid()
    #-------------------------------------------------------------------
    def initialize_Qs_grid(self, SILENT=True, REPORT=False):

        #--------------------------------------------------------
        #NOTES:  Qs = annual sed. discharge [m^3/ yr]
        #        Q  = annual discharge [m^3 / yr]
        #        S  = slope [unitless]
        #        k  = coefficient [(m^3 / yr)^(1 - m)]
        #        m  = discharge exponent (usually in [1,2])
        #        n  = slope exponent (usually in [1,2])

        #        The standard formula Qs = k * Q^m * S^n is for
        #        the case where Q and Qs have units of m^3/sec.
        #        The extra factor below makes the formula valid
        #        for the case where units of both are m^3/yr.
        #        That is, Qs' = fac * kf * (Q')^m * S^n.

        #        The Slope_Grid function checks for negative
        #        slopes and either adjusts them or aborts.
        #--------------------------------------------------------
        # NB!    S and A are zero on edges, so Q and Qs are.
        #--------------------------------------------------------        
        if not(SILENT):
            print('Initializing sed. discharge grid...')
        
        fac1    = self.secs_per_year ** (np.float64(1) - self.m)
        fac2    = self.K * (self.Q ** self.m)
        fac3    = (self.S ** self.n)
        self.Qs = fac1 * fac2 * fac3     #[m^3 / year]

        #---------------------------------
        # Save the min and max Qs values
        #---------------------------------
        self.Qs_min = np.nanmin( self.Qs )
        self.Qs_max = np.nanmax( self.Qs )
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            Qs_str = str(self.Qs_min) + ', ' + str(self.Qs_max)
            print('    min(Qs), max(Qs) = ' + Qs_str + ' [m^3/yr]')

    #   initialize_Qs_grid()
    #-------------------------------------------------------------------
##    def initialize_D_grid_old(self, SILENT=True, REPORT=False):
##
##        #------------------------------------------------------------
##        # Note: The LEM has the form:
##        #          div( D * grad(f) ) = f_t - U
##        #       where D is the diffusion coefficient, and given by:
##        #          D = (qs / S),   where qs = Qs / w.
##        #       Here, Qs(S,A) is computed as in last function, and
##        #       we can take w = dw (computed by D8 component).
##        #
##        #       For linear diffusion and constant D, the CFL
##        #       condition is:    dt < (dx^2) / (2 * D).
##        #
##        #       However, Omelchenko and Karimabadi (2006, p. 186)
##        #       discuss an example 1D problem where instead of 2*D
##        #       in the denominator, they use D_(i-1/2) + D_(i+1/2).
##        #
##        #       For our case, it may be better to use something
##        #       similar, e.g. (D[IDs] + D[pIDs]).
##        #------------------------------------------------------------
##        # NOTE: We have A=0 and S=0 where D8 code is zero, as in
##        #       pits and flats.  As a result, we have Q=0 and Qs=0,
##        #       which means that D will currently be undefined
##        #       anywhere the D8 code is zero.
##        #------------------------------------------------------------
##        if not(SILENT):    
##            print 'Initializing D grid...'
##        
##        self.D = (self.Qs / (self.d8.dw * self.S))
##        self.D.flat[ self.base_IDs ] = 0
##        
##        #--------------------------------
##        # Save the min and max D values
##        #--------------------------------
##        self.D_min = np.nanmin( self.D )
##        self.D_max = np.nanmax( self.D )
##
##
##        w = np.where( self.D > 0 ) ####################  FOR TESTING
##        D_min_pos = np.nanmin( self.D[w] )
##
##        print 'D = (Qs / (dw * S))'
##        print 'D = 0 where D8 code = 0.'     
##        print 'min(D), max(D) =', self.D_min, self.D_max
##        print 'min positive D-value =', D_min_pos
##        print ' '
##        
##    #   initialize_D_grid_old()
    #-------------------------------------------------------------------
    def initialize_D_grid(self, SILENT=True, REPORT=False):

        #------------------------------------------------------
        # Note: For linear diffusion and constant D, the CFL
        #       condition is:    dt < (dx^2) / (2 * D).
        #       For our new CFL stability condition, we have:
        #       D = (Qs / delta_z).
        #------------------------------------------------------
        
        #---------------------------------------------
        # Compute downstream elevation drop, between
        # each cell in DEM and its D8 patent
        #---------------------------------------------
        # del_z will always be positive as long as
        # LINK_FLATS = False. (flats have D8 code 0)
        #---------------------------------------------
        pID_grid = self.d8.parent_ID_grid
        pIDs     = divmod( pID_grid, self.nx )
        del_z    = (self.DEM - self.DEM[ pIDs ])

        #---------------------------------------
        # Note that (del_z[0, 0] == 0), since
        # (pID[0,0] == 0), so (D[0,0] == NaN).
        #---------------------------------------
        self.D = (self.Qs / del_z )
        self.D.flat[ self.base_IDs ] = 0
        
        #--------------------------------
        # Save the min and max D values
        #--------------------------------
        self.D_min = np.nanmin( self.D )
        self.D_max = np.nanmax( self.D )

        w = np.where( self.D > 0 )     ################  FOR TESTING
        D_min_pos = np.nanmin( self.D[w] )

        print('D = (Qs / delta_z).')
        print('D = 0 where D8 code = 0.')
        print('min(D), max(D) =', self.D_min, self.D_max)
        print('min positive D-value =', D_min_pos)
        print(' ')
        
    #   initialize_D_grid()
    #-------------------------------------------------------------------
    def initialize_dz_dt_grid(self, SILENT=True, REPORT=False):
        
        #-------------------------------------------------------------
        # Notes: dV_dt is the net rate at which the volume in a
        #        grid cell is changing.  It is computed as the
        #        difference between all of the fluxes into the grid
        #        cell and flux out of the grid cell. (D8-based here)
        #        dz_dt = (dV_dt / da)  [m/yr].
        #
        #        (dz_dt > 0) => deposition at that cell
        #        (dz_dt < 0) => erosion at that cell
        #
        #        Qs    = kf * Q^mf * S^nf = sed. discharge [m^3/yr]
        #        da    = pixel area grid  [m^2]
        #-------------------------------------------------------------        
        if not(SILENT):    
            print('Initializing dz/dt grid...')

        #----------------------------------------------
        # Initialize dz_dt with the sediment flux out
        # of each pixel divided by the pixel area.
        # Pixel area, self.da, may be scalar or grid.
        # Should be 0 for pixels with flow code of 0.
        #----------------------------------------------
        # Note that in timestep, dt,  dz = dz_dt * dt
        #----------------------------------------------
        flux_out   = self.Qs / self.da
        self.dz_dt = -flux_out   # [m^3 / year]

        #-----------------------------------------
        # Add the contribution to to uplift, U.
        # self.U has units of [mm/yr], but
        # self.U_mpyr has units of [m/yr].
        #-----------------------------------------
        print('Uplift rate [m/yr] =', self.U_mpyr)
        self.dz_dt += self.U_mpyr
        
        #-----------------------------------------
        # Add contributions from neighbor pixels
        #-----------------------------------------
        if (self.d8.n1 != 0):    
            self.dz_dt[ self.d8.p1 ] += flux_out[ self.d8.w1 ]
        if (self.d8.n2 != 0):    
            self.dz_dt[ self.d8.p2 ] += flux_out[ self.d8.w2 ]
        if (self.d8.n3 != 0):    
            self.dz_dt[ self.d8.p3 ] += flux_out[ self.d8.w3 ]
        if (self.d8.n4 != 0):    
            self.dz_dt[ self.d8.p4 ] += flux_out[ self.d8.w4 ]
        if (self.d8.n5 != 0):    
            self.dz_dt[ self.d8.p5 ] += flux_out[ self.d8.w5 ]
        if (self.d8.n6 != 0):    
            self.dz_dt[ self.d8.p6 ] += flux_out[ self.d8.w6 ]
        if (self.d8.n7 != 0):    
            self.dz_dt[ self.d8.p7 ] += flux_out[ self.d8.w7 ]
        if (self.d8.n8 != 0):    
            self.dz_dt[ self.d8.p8 ] += flux_out[ self.d8.w8 ]

        #----------------------------------------------
        # This is only necessary for "completeness"
        # since elevations at base_IDs are set by the
        # update_base_level() function.
        #----------------------------------------------        
        # Use Base-level Lowering Rate (BLR) to set
        # dz_dt on edges.
        # self.BLR_mpyr has units of [m/yr].
        #----------------------------------------------
        # base_IDs are defined in the function:
        # erode_base.initialize_boundary_conditions()
        #----------------------------------------------
        print('Base-level lowering rate [m/yr] =', self.BLR_mpyr)
        print(' ')
        self.dz_dt.flat[ self.base_IDs ] = -self.BLR_mpyr

        #----------------------------------------------------
        # If we have U>0 and BLR=0, we need this (9/30/10)
        #----------------------------------------------------
        ## self.dz_dt.flat[ self.base_IDs ] = 0
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            dz_dt_min = self.dz_dt.min()
            dz_dt_max = self.dz_dt.max()
            print('(dz_dt_min, dz_dt_max) =', dz_dt_min, dz_dt_max)
            print(' ')
            
        #--------------------------------------
        # For testing.  Check where dz_dt = 0
        #--------------------------------------
##        w  = np.where(self.dz_dt == 0)
##        nw = w[0].size
##        if (nw != 0):
##            print 'dz_dt grid is zero at', nw, 'pixel(s).'
##            if (nw <= 20):
##                for k in xrange(nw):
##                    print '(row, col) =', w[0][k], w[1][k]
##            self.dz_dt[w] = self.BLR_mpyr  # (make positive)
        
    #   initialize_dz_dt_grid()
    #-------------------------------------------------------------------
##    def initialize_dt_grid4(self, SILENT=True, REPORT=False):
## 
##        #---------------------------------------------------------------
##        # Notes: This version is based on the following rule:
##        #
##        #        "No pixel can become lower than its D8 parent in
##        #        one timestep.  A one-pixel pit cannot be raised
##        #        higher than any of its neighbors (or perhaps just
##        #        its D8 kids) in one timestep."
##        #
##        #        Applying this rule everywhere has the consequence
##        #        that no pixel can become higher than any of its
##        #        D8 child pixels in one timestep.
##        #
##        #        z2[ID]  = z[ID]  + z_dot[ID]  * dt
##        #        z2[pID] = z[pID] + z_dot[pID] * dt
##        #
##        #        where pID is the D8 parent pixel for pixel at ID.
##        #        Note that same dt is used for both. If we now
##        #        subtract these 2 equations, we find that:
##        #
##        #        If (z_dot[pID] > z_dot[ID]), then we need:
##        #            dt < (z[ID] - z[pID]) / (z_dot[pID] - z_dot[ID])
##        #        In this case, z[ID] and z[pID] are converging,
##        #        so the slope, S[ID], is getting smaller.
##        
##        #        If (z_dot[pID] = z_dot[ID]), then we need:
##        #            dt < Infinity
##        #        In this case, z[ID] and z[pID] are maintaining
##        #        the same difference (and S[ID]) over time.
##        #
##        #        If (z_dot[pID] < z_dot[ID]), then we need:
##        #            dt > -(z[ID] - z[pID]) / (z_dot[ID] - z_dot[pID])
##        #            dt > (z[ID] - z[pID]) / (z_dot[pID] - z_dot[ID])
##        #            dt > 0  (since both diffs above > 0)
##        #        In this case, z[ID] and z[pID] are diverging,
##        #        so the slope, S[ID], is getting larger.
##        #---------------------------------------------------------------
##        # NB!    If we use CFL_factor = 1, then we'll have
##        #        z2[ID] = z2[pID] instead of z2[ID] > z2[pID].
##        #---------------------------------------------------------------
##        if not(SILENT):
##            print 'Initializing dt_grid...'
##
##        #-----------------------
##        # Initialize some vars
##        #--------------------------------------------------
##        # In erode_d8_global.py, a value of 1e-2 is used.
##        #--------------------------------------------------
##        self.dt_too_small = np.float64(1e-4)   ## (need 1e-3 or less, 9/8/100
##        # self.dt_too_small = np.float32(1e-6)
##        # self.dt_too_big = 5   # (should there be a max?)
##        self.dt_limit = 1e20  # [years]
##        
##        #------------------------------------------
##        # Set a "CFL safety factor" (< 1)
##        #-------------------------------------------------
##        # (9/13/10). Using 0.9 works pretty well, but
##        # using 0.88 & smaller causes it to "get stuck".
##        #-------------------------------------------------
##        ## self.CFL_factor     = np.float64(0.2)
##        ## self.CFL_factor     = np.float64(0.4)
##        ## self.CFL_factor     = np.float64(0.5)
##        self.CFL_factor     = np.float64(0.7)         ## somewhat optimal ?
##        ## self.CFL_factor     = np.float64(0.8)
##        ## self.CFL_factor     = np.float64(0.85)
##        ## self.CFL_factor     = np.float64(0.88)  # (threshold)
##        ## self.CFL_factor     = np.float64(0.89)
##        ## self.CFL_factor     = np.float64(0.9)
##        ## self.CFL_factor     = np.float64(0.95)
##        ## self.CFL_factor     = np.float64(1)
##        self.min_allowed_dz = np.float64(1e-6)      ######
##        ## self.min_allowed_dz = np.float64(1e-8)   ######
##        self.max_n_now = np.int32(0)
##        
##        #------------------------------------
##        # Compute downstream elevation drop
##        #------------------------------------------------
##        # del_z will always be positive as long as D8
##        # flow code is not zero and LINK_FLATS = False.
##        #------------------------------------------------
##        # NB! del_z will be wrong wherever pID=0 since
##        # it will be computed from DEM[0].  Places with
##        # pID=0 are fixed further down.
##        #------------------------------------------------
##        pIDs  = self.d8.parent_ID_grid
##        del_z = (self.DEM - self.DEM.flat[ pIDs ])
##        
##        #-------------------------------------------------
##        # See Notes above. This can be pos, neg or zero.
##        #---------------------------------------------------
##        # NB! del_dz_dt will be wrong wherever pID=0 since
##        # it will be computed from dz_dt[0].  Places with
##        # pID=0 are fixed further down.
##        #---------------------------------------------------
##        del_dz_dt = (self.dz_dt - self.dz_dt.flat[ pIDs ])
##       
##        #------------------------------
##        # Compute stable dt, in years
##        #-----------------------------------------------------
##        # If (del_dz_dt < 0), z[ID] & z[pID] are converging.
##        # If (del_dz_dt > 0), z[ID] & z[pID] are diverging.
##        # If (del_dz_dt = 0), z[ID] & z[pID] are neither.
##        #-----------------------------------------------------
##        self.dt_grid = np.zeros( [self.ny, self.nx], dtype='float64' )
## 
##        w1 = np.where( del_dz_dt != 0 )
##        if (w1[0].size > 0):
##            abs_del_dz_dt = np.absolute( del_dz_dt[w1] )
##            TL_term  = 0.0  # (only in the beginning)
##            dt_equal = (del_z[w1] + TL_term) / abs_del_dz_dt
##            self.dt_grid[ w1 ] = self.CFL_factor * dt_equal
##
####            TL_term  = 0.0    # (only in the beginning)
####            dt_equal = (del_z[w1] + TL_term) / del_dz_dt[w1]
####            dt_equal = np.absolute( dt_equal )   #########
####            self.dt_grid[ w1 ] = self.CFL_factor * dt_equal
##            #--------------------------------------------------
##            # (10/20/10) This is only happening at some pits
##            # and their values are set separately below.
##            #--------------------------------------------------
####            w0 = np.where( dt_equal < 0 )
####            if (w0[0].size > 0):
####                print '############  In initialize_dt_grid4():'
####                print '############  dt_equal < 0 at IDs:'
####                w0_rows = w1[0][w0]
####                w0_cols = w1[1][w0]
####                w0_IDs  = (w0_rows * self.nx) + w0_cols
####                print w0_IDs
####                print ' '
##
##        #-----------------------------------------------------
##        # (10/21/10) Pits might match some of this, so make
##        # sure it comes before pit handling section.
##        #-----------------------------------------------------
##        # (10/22/10) Initial dt_grid has no dt_limit values
##        # in the interior, so this may only happen on edges.
##        #-----------------------------------------------------
##        w2 = np.where( del_dz_dt == 0 )
##        if (w2[0].size > 0):
##            self.dt_grid[ w2 ] = self.dt_limit
##            
##        w3 = np.where( del_z <= 0)
##        if (w3[0].size > 0):
##            self.dt_grid[ w3 ] = self.dt_limit
##            
##        #----------------------------------------------------
##        # Deactivate all pixels that have a flow code of 0.
##        # Their z values can then only be changed when they
##        # are synchronized by one of their neighbors.
##        #----------------------------------------------------
##        w4 = np.where( pIDs == 0 )
##        if (w4[0].size > 0):
##            self.dt_grid[ w4 ] = self.dt_limit
##      
##        #------------------------------------------------------------
##        # For each pit, set its dt-value so that it will fill in
##        # one timestep. It will fill to the level of the pixel that
##        # is its lowest neighbor at that future time. (10/18/10)
##        #------------------------------------------------------------
##        # NB! Pits are lower than *all* of their neighbors, but
##        #     some of them may not flow into the pit.  We must have
##        #     (del_z > 0)  for the pit's D8 kids, but we could have
##        #     (del_z <= 0) for other neighbors.
##        #------------------------------------------------------------
##        # Two issues: Pit may not become *higher* than any of its
##        # D8 kids (in one step), BUT pit should only be filled to
##        # the point where it is as high as its lowest neighbor.
##        # Are both issues covered as written ?? (10/19/10)
##        #------------------------------------------------------------
##        # (10/19/10) Checked that pits are being identified and
##        # saved correctly in self.d8.pit_IDs.  However, the list
##        # includes "pits" on the edges.
##        #------------------------------------------------------------
##        print 'Number of pits =', self.d8.pit_IDs.size
##        print '#### Pit IDs ='
##        print self.d8.pit_IDs
##        print ' '
##
##        ## interior_pit_IDs = self.get_interior_pit_IDs( FLATS=True )        
##        interior_pit_IDs = self.get_interior_pit_IDs( FLATS=False )
##        print 'Number of interior pits =', interior_pit_IDs.size
##        print '#### Interior pit IDs ='
##        print interior_pit_IDs
##        print ' '
##        
##        ## for ID in self.d8.pit_IDs:
##        for ID in interior_pit_IDs:
##            #---------------------------------------------------
##            # (10/20/10) Results seem very similar now whether
##            #  we use D8 kids or all 8 neighbors.  COMPARE.
##            #---------------------------------------------------
##            ## kid_IDs = self.d8.get_kid_IDs( ID )
##            kid_IDs = self.d8.get_neighbors( ID )
##            if (kid_IDs.size > 0):
##                # Use this if we want to include flats in with pits. (10/20/10)
####                wcon = np.where( np.logical_and( del_dz_dt.flat[ kid_IDs ] < 0,
####                              self.DEM.flat[ kid_IDs ] > self.DEM.flat[ID] ) )
##                wcon = np.where( del_dz_dt.flat[ kid_IDs ] < 0 ) # (converging)
##                if (wcon[0].size > 0):
##                    kid_IDs = kid_IDs[ wcon ]
##                    abs_del_dz_dt = np.absolute(del_dz_dt.flat[ kid_IDs ])
##                    del_z_term = del_z.flat[ kid_IDs ]
##                    TL_term    = 0.0   # (only in the beginning) 
##                    dt_equal = (del_z_term + TL_term) / abs_del_dz_dt
##                    self.dt_grid.flat[ ID ] = dt_equal.min()
##                    if (ID in self.d8.flat_IDs):
##                        self.dt_grid.flat[ ID ] *= self.CFL_factor
##                    ## self.dt_grid.flat[ ID ] = dt_equal.min() * self.CFL_factor
##                    ## self.dt_grid.flat[ ID ] = dt_equal.min() * 1.1
##                    #---------------------------------------------
##                    # (10/22/10) EXPERIMENT.  If one of pit's kids
##                    # is active_ID, then change kid's dt-value.
##                    # Only need this in update_dt_grid4().
##                    #---------------------------------------------
####                    if (self.active_ID in kid_IDs):
####                        self.dt_grid.flat[ self.active_ID ] /= self.CFL_factor
##
##                    #---------------------------------------------
##                    # (10/22/10) EXPERIMENT. Make sure kid's are
##                    # updated AFTER their pit parent.
##                    #---------------------------------------------
####                    epsilon = 1e-08
####                    self.dt_grid.flat[ kid_IDs ] = np.maximum( dt_equal.min() + epsilon,
####                                                                  self.dt_grid.flat[ kid_IDs ] )
##                    #----------------------------------------------------
##                    # (10/22/10) EXPERIMENT. Deactivate pits neighbors.
##                    # Seems to be a bad idea.
##                    #----------------------------------------------------
##                    ## self.dt_grid.flat[ kid_IDs ] = self.dt_limit
##                    
##                    ## wm = np.where( dt_equal == dt_equal.min() )
##                    ## self.dt_grid.flat[ ID ] = dt_equal[ wm ]
##                    #----------------------------------------------------
##                    # (10/20/10) This hasn't been triggered yet.
##                    #----------------------------------------------------                    
##                    if (dt_equal.min() < 0):
##                        print '############  In initialize_dt_grid4():'
##                        print '############  dt_equal < 0 at PIT ID =', ID
##                        print ' '
##
####        #-------------------------------------------------------------
####        # (10/22/10) Where dt is very small, divide it by CFL_factor 
####        #-------------------------------------------------------------
####        w5 = np.where( self.dt_grid < 1e-1 )
####        if (w5[0].size > 0):
######            self.dt_grid[ w5 ] = np.maximum( self.dt_grid[w5], self.dt_grid[ pIDs[w5] ] )
####            self.dt_grid[ w5 ] = self.dt_grid[ w5 ] / self.CFL_factor
##        
##        #------------------------------------------
##        # Deactivate pixels where dt is too small
##        #------------------------------------------
##        # Commented out on: 10/21/10.
####        w5 = np.where( self.dt_grid < 1e-3 )
####        if (w5[0].size > 0):
####            self.dt_grid[ w5 ] = self.dt_limit
##          
##        #---------------------------------------------------------
##        # If we have U>0 and BLR=0, then we need this. (9/20/10)
##        #---------------------------------------------------------
##        #### if (self.U_mpyr > 0) and (self.BLR_mpyr == 0):
##        ## self.dt_grid.flat[ self.base_IDs ] = self.dt_limit
##            
##        #--------------------------
##        # Save the min and max dt
##        #--------------------------
##        self.dt_min = np.nanmin( self.dt_grid )
##        self.dt_max = np.nanmax( self.dt_grid )
####        self.dt_min = self.dt_grid.min()
####        self.dt_max = self.dt_grid.max()
##        if not(SILENT):
##            print '#########################################'
##            print ' dt_min =', np.around(self.dt_min, 2)
##            print ' in initialize_dt_grid()'
##            print '#########################################'
##
##        #-------------------------------
##        # Don't let dt get too small ?
##        #-----------------------------------------------------
##        # dt_min must be less than 1e-4 for the case m=1.5,
##        # n=1.0, K=1.0. Making K smaller allows dt_too_small
##        # to be bigger. Now K is smaller.
##        #-----------------------------------------------------
##        # Fixed units issue in Qs_Grid so should be able to
##        # reduce dt_too_small.
##        #-----------------------------------------------------
##        ## if (self.dt_min <= 0):
##        if (self.dt_min < self.dt_too_small):
##            print '### WARNING: dt < ' + str(self.dt_too_small)
##            print '        dt = ' + str(self.dt_min)
##            print ' '
##            if (self.dt_min <= 0):
##                sys.exit()
##            
##        print 'min(del_z), max(del_z) =', del_z.min(), del_z.max()
##        print 'min(del_z_dot), max(del_z_dot) =', del_dz_dt.min(), del_dz_dt.max()
##        print 'min(dt), max(dt) =', self.dt_grid.min(), self.dt_grid.max()
##        wm = np.where( self.dt_grid != self.dt_grid.max() )
##        print 'max(dt_used)     =', self.dt_grid[wm].max()
##        print ' '
##        
##    #   initialize_dt_grid4()
    #-------------------------------------------------------------------
    def initialize_dt_grid5(self, SILENT=True, REPORT=False):
 
        #---------------------------------------------------------------
        # Notes: This version is based on the following stability
        #        condition:
        #
        #        dt < dx * dy * [z(k) - z(p)] / [2 * Qs(k)]
        #
        #        This condition can be derived from the fact that
        #        the rate of sediment transport between cells k and
        #        p, Qs(k), must be less than the rate one gets by
        #        taking half of the sediment volume at cell k that
        #        is higher than cell p and moving it to cell p in
        #        one timestep.
        #
        #        This condition ensures a well-defined timestep
        #        everywhere except at pits and flats.  It may,
        #        however, allow the elevation at cell p to exceed
        #        the elevation at cell k during a (global) timestep
        #        if there is a sediment flux from other neighbors
        #        into cell k and sufficiently low flux out of cell p.
        #
        #        It could be that with global timesteps, z(k) will
        #        be raised by inflow from its neighbors instead of
        #        lowered by the amount used here.  But the result
        #        should still be valid.
        #---------------------------------------------------------------
        if not(SILENT):
            print('Initializing dt_grid...')

        #####################################
        # SOME OF THESE ARE ALSO SET IN:
        #    erode_base.set_constants()
        #####################################
        
        #-----------------------
        # Initialize some vars
        #--------------------------------------------------
        # In erode_d8_global.py, a value of 1e-2 is used.
        #--------------------------------------------------
        self.max_n_now    = np.int32(0)
        self.dt_too_small = np.float64(1e-4)   ## (need < 1e-3, 9/8/10)
        self.dt_limit     = np.float64(1e20)   # [years]
        self.max_dz_for_pits = np.float64( 0.01 )  # [meters], (3/7/15)

        
        # self.dt_too_small = np.float64(1e-6)
        # self.dt_too_big   = 5   # (should there be a max?)
        
        #####################################################
        # (11/20/11) CHECK SOON IF THIS IS STILL USED !!
        #####################################################
        self.min_allowed_dz = np.float64(1e-6) 
        ## self.min_allowed_dz = np.float64(1e-8)
        ## self.max_allowed_dz = np.float64(0.001)
        
        ###################################################
        # CFL_factor is now set in the CFG file.
        # erode_base.initialize() calls
        # erode_d8_local.initialize_config_vars() AND
        # erode_d8_local.initialize_computed_vars() calls
        # erode_d8_local.initialize_dt_grid5().
        ###################################################        
        
 
        #------------------------------------
        # Compute downstream elevation drop
        #------------------------------------------------
        # del_z will always be positive as long as D8
        # flow code is not zero and LINK_FLATS = False.
        #------------------------------------------------
        # NB! del_z will be wrong wherever pID=0 since
        # it will be computed from DEM[0].  Places with
        # pID=0 are fixed further down.
        #------------------------------------------------
        pIDs  = self.d8.parent_ID_grid
        del_z = (self.DEM - self.DEM.flat[ pIDs ])
        
        #-------------------------------------
        # This should be equivalent (3/1/12)
        #-------------------------------------
##        pID_grid = self.d8.parent_ID_grid
##        pIDs     = divmod( pID_grid, self.nx )
##        del_z    = (self.DEM - self.DEM[ pIDs ])
        
        #------------------------------
        # Compute stable dt, in years
        #------------------------------
        self.dt_grid = np.zeros( [self.ny, self.nx], dtype='float64' )
 
        w1 = np.where( self.Qs > 0 )
        if (w1[0].size > 0):
            if (self.da.size == 1):
                da = self.da
            else:
                da = self.da[ w1 ]
            ## da = (self.rti.xres * self.rti.yres)
            dt_equal = (da / 2) * del_z[w1] / self.Qs[w1]
            self.dt_grid[ w1 ] = self.CFL_factor * dt_equal

            
        #------------------------------------------------
        # Wherever pIDs = 0, z[pIDs]=z[0] and we could
        # easily have (del_z < 0). But pIDs are handled
        # in the next block, so we shouldn't need this.
        # It was commented out on 2/28/12.
        #------------------------------------------------
##        w2 = np.where( del_z < 0)
##        if (w2[0].size > 0):
##            self.dt_grid[ w2 ] = self.dt_limit
##            print '###########################################'
##            print ' ERROR in initialize_dt_grid5():'
##            print '       del_z < 0 at', w2[0].size, 'cells.'
##            print '###########################################'
##            print ' '
##            sys.exit()

        #----------------------------------------------------
        # For pits and flats, use dz_dt and max_dz_for_pits
        # to set dt.  They are now put on heap. (3/7/12)
        #----------------------------------------------------
        w2 = np.where( self.d8.d8_grid == 0 )
        if (w2[0].size > 0):
            abs_dz_dt = np.absolute( self.dz_dt[ w2 ] )
            self.dt_grid[ w2 ] = self.max_dz_for_pits / abs_dz_dt
     
        #-----------------------------------------------------
        # Deactivate all edge pixels (but not pits & flats).
        # Their z values can then only be changed when they
        # are synchronized by one of their neighbors.
        #-----------------------------------------------------
        self.dt_grid.flat[ self.base_IDs ] = self.dt_limit

        
        #----------------------------------------------------
        # Deactivate all pixels that have a flow code of 0.
        # Their z values can then only be changed when they
        # are synchronized by one of their neighbors.
        # (Used before 3/7/12 experiment.)
        #----------------------------------------------------
##        ## w3 = np.where( pIDs == 0 )
##        w3 = np.where( self.d8.d8_grid == 0 )  ### (2/28/12; same?)
##        if (w3[0].size > 0):
##            self.dt_grid[ w3 ] = self.dt_limit
##            ## self.dt_grid[ w3 ] = 999

            
        #------------------------------------------------------------
        # For each pit, set its dt-value so that it will fill in
        # one timestep. It will fill to the level of the pixel that
        # is its lowest neighbor at that future time. (10/18/10)
        #------------------------------------------------------------
        # NB! Pits are lower than *all* of their neighbors, but
        #     some of them may not flow into the pit.  We must have
        #     (del_z > 0)  for the pit's D8 kids, but we could have
        #     (del_z <= 0) for other neighbors.
        #------------------------------------------------------------
        # Two issues: Pit may not become *higher* than any of its
        # D8 kids (in one step), BUT pit should only be filled to
        # the point where it is as high as its lowest neighbor.
        # Are both issues covered as written ?? (10/19/10)
        #------------------------------------------------------------
        # (10/19/10) Checked that pits are being identified and
        # saved correctly in self.d8.pit_IDs.  However, the list
        # includes "pits" on the edges.
        #------------------------------------------------------------
        print('Number of pits =', self.d8.pit_IDs.size)
        print('Pit IDs =')
        print(self.d8.pit_IDs)
        print(' ')
        
        interior_pit_IDs = self.get_interior_pit_IDs( FLATS=False )
        print('Number of interior pits =', interior_pit_IDs.size)        
          
        #---------------------------------------------------------
        # If we have U>0 and BLR=0, then we need this. (9/20/10)
        #---------------------------------------------------------
        #### if (self.U_mpyr > 0) and (self.BLR_mpyr == 0):
        ## self.dt_grid.flat[ self.base_IDs ] = self.dt_limit
            
        #--------------------------
        # Save the min and max dt
        #----------------------------------------------
        # (1/25/12) nanmin() and nanmax() are slower,
        # and it seems we don't need them here.
        #----------------------------------------------
        self.dt_min = self.dt_grid.min()
        self.dt_max = self.dt_grid.max()
##        self.dt_min = np.nanmin( self.dt_grid )
##        self.dt_max = np.nanmax( self.dt_grid )

        if not(SILENT):
            print('#########################################')
            print(' dt_min =', np.around(self.dt_min, 2))
            print(' in initialize_dt_grid()')
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
        ## if (self.dt_min <= 0):
        if (self.dt_min < self.dt_too_small):
            print('### WARNING: dt < ' + str(self.dt_too_small))
            print('        dt = ' + str(self.dt_min))
            print(' ')
            if (self.dt_min <= 0):
                sys.exit()
            
        print('min(del_z), max(del_z) =', del_z.min(), del_z.max())
        ## print 'min(del_z_dot), max(del_z_dot) =', del_dz_dt.min(), del_dz_dt.max()
        print('min(dt), max(dt) =', self.dt_grid.min(), self.dt_grid.max())
        wm = np.where( self.dt_grid != self.dt_grid.max() )
        print('max(dt_used)     =', self.dt_grid[wm].max())
        print(' ')

        ##########################################
##        print 'CFL_factor.dtype =', self.CFL_factor.dtype
##        print 'da.shape       =', da.shape
##        print 'da.dtype       =', da.dtype
##        print 'da.min()       =', da.min()
##        print 'da.max()       =', da.max()
##        print 'dt_limit.dtype =', self.dt_limit.dtype
##        print 'del_z.shape    =', del_z.shape
##        print 'del_z.dtype    =', del_z.dtype
##        print 'Qs.shape       =', self.Qs.shape
##        print 'Qs.dtype       =', self.Qs.dtype
##        print 'dt_grid.shape  =', self.dt_grid.shape
##        print 'dt_grid.dtype  =', self.dt_grid.dtype   
                    
    #   initialize_dt_grid5()                  
    #-------------------------------------------------------------------
    def initialize_dz_target_grid2(self, SILENT=True, REPORT=False):
      
        #----------------------------------------
        # Use dt_grid to determine dz_target
        # Note: dt_grid already has CFL_factor.
        #----------------------------------------------
        # (9/12/10) Allow dz_target to keep its sign.
        #----------------------------------------------
        self.dz_target_grid = self.dt_grid * self.dz_dt

        #######################################################
        # (2/28/12) Set dz_target to zero for all pixels
        # that have (d8_code == 0) that are NOT on the edge.
        # Cells with (d8_code == 0) are never drawn from the
        # heap.  They get updated only by their neighbors.
        # The edge cells should be in "wait_IDs" and all of
        # the others should be in "now_IDs".
        # See update_dz_target_grid2() also.
        #######################################################
        test1 = (self.d8.d8_grid == 0)
        test2 = self.not_edge_grid
        w = np.where( np.logical_and( test1, test2) )
        if (w[0].size > 0):
            self.dz_target_grid[ w ] = 0.0

        #------------------------------------------------------
        # This block was used before 2/28/12. It excludes the
        # edges but may not get same IDs as above.
        #------------------------------------------------------
        # Set dz_target to 0 for pits so they will always be
        # considered "affected", to keep their dt_grids, etc.
        # up-to-date.  (11/4/10)
        #------------------------------------------------------
        ## interior_pit_IDs = self.get_interior_pit_IDs()
##        interior_pit_IDs = self.get_interior_pit_IDs(FLATS=True)
##        self.dz_target_grid.flat[ interior_pit_IDs ] = 0.0
            
    #   initialize_dz_target_grid()
    #-------------------------------------------------------------------
    def initialize_T_last_grid(self, SILENT=True, REPORT=False):

        if not(SILENT):
            print('Initializing T_last grid...')

        self.T_last = np.zeros([self.ny, self.nx], dtype='Float64')
        
    #   initialize_T_last_grid()
    #-------------------------------------------------------------------
    def initialize_T_next_grid(self, SILENT=True, REPORT=False):

        if not(SILENT):
            print('Initializing T_next grid...')
            
        self.T_next = self.dt_grid.copy()

        #--------------------------------
        # Prepare to use priority queue
        #--------------------------------
        T_heap = list(zip( self.T_next.flat, self.d8.ID_grid.flat ))
        heapq.heapify( T_heap )
        self.T_heap = T_heap
        self.n_heap = len( T_heap )
        print('Initial size of heap is:', self.n_heap)
            
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            s       = np.argsort(self.T_next, axis=None)
            T_sort  = self.T_next.flat[ s ]
            ID_sort = self.d8.ID_grid.flat[ s ]
            print('T_sort[0:10]  =', T_sort[0:10])
            print('ID_sort[0:10] =', ID_sort[0:10])

        ## Tn_min = self.T_next.min()
        ## Tn_max = self.T_next.max()
        ## print 'Tn_min, Tn_max =', Tn_min, Tn_max
            
    #   initialize_T_next_grid()
    #-------------------------------------------------------------------
    def save_initial_grids(self, SILENT=True, REPORT=False,
                           ABORT=False):
 
        #---------------------------------------
        # Save some initial grids as RTG files
        #--------------------------------------------------
        # This is called from initialize_computed_vars()
        # after all grids have been initialized.
        #--------------------------------------------------
        dt0_file = (self.out_directory + self.case_prefix + '_dt0.rtg')
        rtg_files.write_grid( self.dt_grid, dt0_file, self.rti,
                              RTG_type='DOUBLE')
        print('#######################################################')
        print('Saved initial dt_grid to: ')
        print('     ' + dt0_file)
        #-------------------------------------------------------------------
        pIDs      = self.d8.parent_ID_grid
        del_z     = (self.DEM - self.DEM.flat[ pIDs ])
        delz_file = (self.out_directory + self.case_prefix + '_delz0.rtg')
        rtg_files.write_grid( del_z, delz_file, self.rti,
                              RTG_type='DOUBLE')
        print('Saved initial del_z grid to:')
        print('     ' + delz_file)
        #-------------------------------------------------------------------
        Qs0_file = (self.out_directory + self.case_prefix + '_Qs0.rtg')
        rtg_files.write_grid( self.Qs, Qs0_file, self.rti,
                              RTG_type='DOUBLE')
        print('Saved initial Qs grid to: ')
        print('     ' + Qs0_file)
        #-------------------------------------------------------------------
        dzdt0_file = (self.out_directory + self.case_prefix + '_dzdt0.rtg')
        rtg_files.write_grid( self.dz_dt, dzdt0_file, self.rti,
                              RTG_type='DOUBLE')
        print('Saved initial dz/dt grid to: ')
        print('     ' + dzdt0_file)
        #-------------------------------------------------------------------
        dztr0_file = (self.out_directory + self.case_prefix + '_dztr0.rtg')
        rtg_files.write_grid( self.dz_target_grid, dztr0_file, self.rti,
                              RTG_type='DOUBLE')
        print('Saved initial dz_target grid to: ')
        print('     ' + dztr0_file)
        #-------------------------------------------------------------------
        print('#######################################################')
        print(' ')

        #------------------
        # Option to abort  
        #------------------
        if (ABORT):
            sys.exit() ########

    #   save_initial_grids()
    #-------------------------------------------------------------------
    def initialize_zone_IDs(self, SILENT=True, REPORT=False):

        #-----------------------------------------------------
        # (2/10/12) New method that uses zone_IDs, etc.
        # We're using fixed-size arrays here to avoid slow
        # calls to np.concatenate().  The largest number
        # of now_IDs in previous runs was around 14.
        # zone_IDs is union of all_now_IDs and all_wait_IDs.
        # This is called by initialize_computed_vars().
        #-----------------------------------------------------
        # Type (int32) must match pIDs in d8_base.py.
        #-----------------------------------------------------
        # If we take out "test4" in update_neighbors(), then
        # active_ID_count for ID = 905 goes up to 147 and
        # all_now_IDs gets up to 453 at one point. (3/1/12)
        #-----------------------------------------------------        
        # nmax = 10000
        # nmax = 1000  # (still failed)
        nmax = 100

        self.zone_IDs    = np.zeros( nmax, dtype='int32' )
        self.all_now_IDs = np.zeros( nmax, dtype='int32' )
##        self.all_wait_IDs = np.zeros( nmax, dtype='int32' )
        
    #   initialize_zone_IDs()
    #-------------------------------------------------------------------
    def check_zone_info(self, SILENT=True, REPORT=False):

        #--------------------------------------------------
        # Note: New test function for zone_IDs. (2/28/12)
        #--------------------------------------------------
        w1 = np.where( self.zone_IDs != -1 )
        print('==========================================')
        print(' Checking information for zone_IDs:')
        print(' ')
        print(' active_ID         =', self.active_ID)
        print(' size(zone_IDs)    =', w1[0].size)
        # print ' zone_IDs_index    =', self.zone_IDs_index
        if (w1[0].size == 0):
            return
        w1b = np.where( self.all_now_IDs != -1 )
        print(' size(all_now_IDs) =', w1b[0].size)
        # print ' all_now_IDs_index =', self.all_now_IDs_index

        #----------------------------
        # Print the actual zone_IDs
        #----------------------------
        print(' ')
        print(' zone_IDs =')
        print(self.zone_IDs[ w1 ])
        print(' ')
        
        #-------------------------
        # Get values at zone_IDs
        #-------------------------
        IDs      = self.zone_IDs[ w1 ]  ####
        d8_codes = self.d8.d8_grid.flat[ IDs ]
        pIDs     = self.d8.parent_ID_grid.flat[ IDs ]
        A        = self.d8.A.flat[ IDs ]
        S        = self.S.flat[ IDs ]
        Q        = self.Q.flat[ IDs ]
        Qs       = self.Qs.flat[ IDs ]
        dt       = self.dt_grid.flat[ IDs ]
        dz_targ  = self.dz_target_grid.flat[ IDs ]
        dz_cap   = self.flux_cap_grid.flat[ IDs ]
        PASS     = True   # (default)
        
        #-------------------------------------
        # Check values where (d8_codes == 0)
        #------------------------------------- 
        w2 = np.where( d8_codes == 0 )
        print(' d8_code = 0 at', w2[0].size, 'cells.')
        #----------------------------------------------- 
        w3 = np.where( pIDs[ w2 ] != 0 )
        if (w3[0].size > 0):
            print(' FAIL: pIDs!=0 where d8_codes=0')
            print('       at', w3[0].size, 'cells.')
            PASS = False
        #-----------------------------------------------    
        w4 = np.where( A[ w2 ] != 0 )
        if (w4[0].size > 0):
            print(' FAIL: A!=0 where d8_codes=0')
            print('       at', w4[0].size, 'cells.')
            PASS = False
        #-----------------------------------------------          
        w5 = np.where( S[ w2 ] != 0 )
        if (w5[0].size > 0):
            print(' FAIL: S!=0 where d8_codes=0')
            print('       at', w5[0].size, 'cells.')
            PASS = False
        #----------------------------------------------- 
        w6 = np.where( Q[ w2 ] != 0 )
        if (w6[0].size > 0):
            print(' FAIL: Q!=0 where d8_codes=0')
            print('       at', w6[0].size, 'cells.')
            PASS = False
        #-----------------------------------------------         
        w7 = np.where( Qs[ w2 ] != 0 )
        if (w7[0].size > 0):
            print(' FAIL: Qs!=0 where d8_codes=0')
            print('       at', w7[0].size, 'cells.')
            PASS = False
        #----------------------------------------------------
        # See detailed comments in update_dz_target_grid2()
        # that resulted from this test.  It seems okay that
        # that can happen, but I don't think it was being
        # handled correctly before 2/29/12.  ##############
        #----------------------------------------------------
##        w8 = np.where( dt[ w2 ] != self.dt_limit )
##        if (w8[0].size > 0):
##            print ' FAIL: dt!=dt_limit where d8_codes=0'
##            print '       at', w8[0].size, 'cells.'
##            PASS = False
        #-------------------------------------------------
        # dz_targ can now be nonzero on edges (2/28/12).
        # This is triggered for active_ID = 134.
        #-------------------------------------------------        
##        w9 = np.where( dz_targ[ w2 ] != 0 )
##        if (w9[0].size > 0):
##            print ' FAIL: dz_targ!=0 where d8_codes=0'
##            print '       at', w9[0].size, 'cells.'
##            PASS = False
        #-------------------------------------------------
        # dz_cap will often be nonzero for pits, flats
        # and edges.  While all pits and flats are now
        # included in "now_IDs", edge_IDs are not.
        # Pits and flats will therefore have been reset
        # to 0 before this test is called by update().
        #-------------------------------------------------        
##        wA = np.where( dz_cap[ w2 ] != 0 )
##        if (wA[0].size > 0):
##            print ' FAIL: dz_cap!=0 where d8_codes=0'
##            print '       at', wA[0].size, 'cells.'
##            PASS = False
            
       
        #------------------------------------
        # Compute downstream elevation drop
        #------------------------------------------------
        # del_z will always be positive as long as
        # LINK_FLATS = False. (flats have D8 code 0)
        #------------------------------------------------
        # NB! del_z will be wrong wherever pID=0 since
        # it will be computed from DEM[0].
        #------------------------------------------------
        z_IDs  = self.DEM.flat[ IDs ]
        z_pIDs = self.DEM.flat[ pIDs ]
        del_z  = (z_IDs - z_pIDs)
        wY     = np.where( d8_codes != 0 )
        wZ     = np.where( del_z[ wY ] <= 0 )
        if (wZ[0].size > 0):
            print(' FAIL: del_z <= 0 where d8_code > 0')
            print('       at', wZ[0].size, 'cells.')
            PASS = False

        #---------------------------
        # Abort if any test fails.
        #---------------------------
        if not(PASS):
            sys.exit()
        else:
            print(' All tests PASSED.')
            print(' ')
            
    #   check_zone_info()
    #-------------------------------------------------------------------
    def reset_zone_IDs(self, ID):

        #-------------------------------------------------
        # Note: zone_IDs and all_now_IDs are initialized
        #       in "initialize_zone_IDs()".
        #-------------------------------------------------
        # It may be faster to set to -1 just up to the
        # associated index.
        #-------------------------------------------------
        
        #---------------------------
        # Reset the zone_IDs array
        #---------------------------
        self.zone_IDs[:]    = -1
        self.zone_IDs[0]    = ID        
        self.zone_IDs_index = np.int32( 1 )

        #------------------------------
        # Reset the all_now_IDs array
        #------------------------------
        self.all_now_IDs[:]    = -1
        self.all_now_IDs[0]    = ID
        self.all_now_IDs_index = np.int32( 1 )

        #-------------------------------------------------
        # This is older; maybe same as all_now_IDs_index
        #-------------------------------------------------
        self.total_n_now = 1
        
        #-------------------------------
        # Reset the all_wait_IDs array
        # (no longer needed)
        #-------------------------------
##        self.all_wait_IDs[:]    = -1
##        self.all_wait_IDs_index = np.int32( 0 )

    #   reset_zone_IDs()
    #-------------------------------------------------------------------
    def append_zone_IDs(self, nbr_IDs):

        #-----------------------------------------------        
        # Note:  Append nbr_IDs to the zone_IDs array.
        #        zone_IDs is initialized with -1's.
        #-----------------------------------------------
        j = self.zone_IDs_index
        n_nbrs = nbr_IDs.size
        top_index = (j + n_nbrs)  # (this is correct)
        if (top_index > self.zone_IDs.size):
            print('############################################')
            print(' ERROR in append_zone_IDs():')
            print('       Splash zone exceeded max size of:')
            print('      ', self.zone_IDs.size, ',')
            print('       at time_index =', self.time_index)
            print(' ')
            sys.exit()

        self.zone_IDs[ j: top_index ] = nbr_IDs
        self.zone_IDs_index += n_nbrs

        #----------------------------------------------------
        # Avoiding "concatenate()", which tends to be slow.
        #----------------------------------------------------
        ## self.zone_IDs = np.concatenate( (self.zone_IDs, nbr_IDs) )

    #   append_zone_IDs()
    #-------------------------------------------------------------------
    def append_all_now_IDs(self, now_IDs):

        #-------------------------------------------------        
        # Note:  Append now_IDs to the all_now_IDs array.
        #        all_now_IDs is initialized with -1's.
        #-------------------------------------------------
        j = self.all_now_IDs_index
        n_NOW = now_IDs.size #######
        top_index = (j + n_NOW)      # (this is correct)
        if (top_index > self.all_now_IDs.size):
            print('############################################')
            print(' ERROR in append_all_now_IDs():')
            print('       all_now_IDs exceeded max size of:')
            print('      ', self.all_now_IDs.size, '.')
            print('       time_index =', self.time_index)
            print(' ')
            sys.exit()    
        self.all_now_IDs[ j: top_index ] = now_IDs
        self.all_now_IDs_index += n_NOW
        self.total_n_now += n_NOW  # (same as all_now_IDs_index)
        
        #---------------------------------------------
        # Keep track of the largest n_NOW ever seen.
        #---------------------------------------------
        self.max_n_now = np.maximum( self.total_n_now,
                                     self.max_n_now )

        #----------------------------------------------------
        # Avoiding "concatenate()", which tends to be slow.
        #----------------------------------------------------
        ## self.all_now_IDs = np.concatenate( (self.all_now_IDs, now_IDs)

    #   append_all_now_IDs()
    #-------------------------------------------------------------------
    def update(self, time=None, SILENT=True, REPORT=False):
     
        if not(SILENT):
            print('Erosion component: Updating...')
        self.status = 'updating'  # (OpenMI 2.0 convention)

        #---------------------------------
        # Print dz_max to track progress
        #---------------------------------
        ### if (self.mode == 'driver'):   # (maybe slower)
        if (self.DRIVER):
            self.print_time_and_value( self.dz_max, 'dz_max', '[m]',
                                       interval=5.0 )

        #-----------------------------------------------------------
        # (11/15/11) Moved this here when I realized that due to
        # recursive nature of update_neighbors(), output was not
        # being written when it was supposed to (skipped, etc.).
        #-----------------------------------------------------------
        # This will write DEM, etc. before anything happens since
        # we haven't called update_T_clock() or update_time() yet.
        #-----------------------------------------------------------
        # (11/15/11)  Note that dt can be 0, but active_ID will
        # be different.  This happens when multiple pixels get
        # scheduled for same update time, T_next.
        #-----------------------------------------------------------        
        self.write_output_files()
        
        #--------------------------------------------------
        # Update the active_ID (where T = T_min)
        # Call update_T_clock() AFTER update_active_ID().
        #--------------------------------------------------        
        self.update_active_ID( SILENT=SILENT, REPORT=REPORT )
        self.update_T_clock()
        ### self.total_n_now = 1  # (moved to reset_zone_IDs)
        
        #--------------------------------------
        # Update values from other components
        #----------------------------------------------
        # NB! This shouldn't be called this often !!!
        #     That is, it must be rewritten to look
        #     at self.T_clock.
        #----------------------------------------------
##        self.update_R()
##        self.update_R_integral()
##        self.update_U()
##        self.update_U_integral()

        #------------------------------
        # Update the DEM at active_ID
        #-----------------------------------------------
        # This uses dz_dt_grid and (T_clock - T_last).
        #-----------------------------------------------
        ## ID = self.active_ID
        ID = np.array([self.active_ID])  # (must be array)
        self.update_DEM(ID, SILENT=SILENT, REPORT=REPORT)
        self.flux_cap_grid.flat[ ID ] = 0
        self.update_T_last_grid(ID, SILENT=SILENT, REPORT=REPORT) # (9/3/10)

        #--------------------------------------------
        # Reset the zone_IDs and all_now_IDs arrays
        #--------------------------------------------
        self.reset_zone_IDs( ID )
        
        #----------------------------------------------------
        # Find all cells in the "splash zone" associated
        # with active_ID and store their IDs in zone_IDs.
        #----------------------------------------------------
        # Cells in this zone consist of two types: those
        # that need to be fully updated (and rescheduled)
        # NOW (now_IDs), and those that store their net dz
        # change in the flux capacitor and can WAIT (called
        # wait_IDs) until their scheduled update time.
        #----------------------------------------------------
        # update_neighbors() recursively builds these lists
        # of IDs, avoiding duplicate entries.  It does this
        # by looking at dz_target, dz_cap and d8_codes.
        #----------------------------------------------------        
        # The DEM, D8 vars (code, A, S), fluxes (Q, Qs) and
        # dz_dt are always updated for *all* zone_IDs.
        # D8 vars are actually updated for zone_IDs as well
        # as *all* of their neighbors.
        #
        ##### This is necessary to avoid backflow and loops
        ##### in the D8 code grid, but could have other
        ##### unknown side effects. (3/1/12)
        #----------------------------------------------------
        # Note that dz/dt is always computed using just
        # (fully) updated D8 flow direction codes.  So we
        # don't need to worry about a flow direction
        # changing before its corresponding flux rate has
        # been used.
        #####################################################       
        self.update_neighbors( ID )  # (update_zone_IDs ??)
        w1          = np.where( self.zone_IDs != -1 )
        zone_IDs    = self.zone_IDs[w1]
        w2          = np.where( self.all_now_IDs != -1 )
        all_now_IDs = self.all_now_IDs[w2]

        #----------------------------------------------       
        # Since np.setdiff1d() is now used within
        # update_neighbors, we shouldn't need to call
        # unique() here. (2/14/12)
        #----------------------------------------------
##        w1 = np.where( self.zone_IDs != -1 )
##        zone_IDs = np.unique( self.zone_IDs[w1] )
##        w2 = np.where( self.all_now_IDs != -1 )
##        all_now_IDs = np.unique( self.all_now_IDs[w2] )

        #------------------------------------------------------
        # We could call update_DEM() and update_T_last_grid()
        # just once here for all zone_IDs.  However, we still
        # need to call update_flux_cap_grid() inside of
        # update_neighbors() in order to build out zone_IDs.
        # This would remove similar calls just for ID above.
        # Might need update_T_last_grid() inside of update_
        # neighbors() also, since it is used to compute
        # dz_cap. (3/3/12)
        #------------------------------------------------------
##        self.update_DEM( zone_IDs, SILENT=SILENT,
##                         REPORT=REPORT)
##        self.flux_cap_grid.flat[ zone_IDs ] = 0
##        self.update_T_last_grid( zone_IDs, SILENT=SILENT,
##                                REPORT=REPORT)
        
        #--------------------------------------------------
        # Synchronize all neighbor cells in "splash zone"
        #--------------------------------------------------
        # This uses current DEM to update D8 vars (A, S
        # and flow_code) as well as the fluxes, Q and Qs.
        #--------------------------------------------------        
        self.update_fluxes4( zone_IDs, SILENT=SILENT,
                             REPORT=REPORT )

        #-----------------------------------------------
        # This uses current fluxes and flow directions.
        #-----------------------------------------------
        self.update_dz_dt_grid( zone_IDs, SILENT=SILENT,
                                REPORT=REPORT)

        ########################################################
        # (3/5/12) Now that DEM, D8_codes and fluxes have been
        # fully updated, we can compute new dt_grid values for
        # *all* cells in zone_IDs.  Any cells for which the
        # new dt-value is less than the current value should
        # perhaps be appended to all_now_IDs right here.
        # Otherwise they may be unstable at their scheduled time.
        #
        # Should we also "grow out" from these new IDs and
        # include those as well?  Note that a cell cannot
        # need a new dt unless its del_z or Qs has changed.
        #
        # Maybe we don't even need dz_target and dz_cap any
        # more.  We could just grow all_now_IDs by computing
        # updated dt-values and appending all cells that now
        # have a lower value.  But we'd need to move update_
        # fluxes() back into update_neighbors() in order to
        # compute updated dt-values for region growing.
        ########################################################
        
        #---------------------------------------------
        # Update dt values, but only for "now_IDs".
        # This uses DEM (for del_z), d8_code and Qs.
        #---------------------------------------------
        self.update_dt_grid6( all_now_IDs, SILENT=SILENT,
                              REPORT=REPORT)


        
        #--------------------------------------------
        # Update dz_target, but only for "now_IDs".
        # This uses dz_dt_grid and dt_grid.
        #--------------------------------------------
        self.update_dz_target_grid2( all_now_IDs, SILENT=SILENT,
                                     REPORT=REPORT)
        
        #-----------------------------------------------------
        # Schedule next update time, but only for "now_IDs".
        # Pairs with (T_next, ID) are put on a heap.
        #-----------------------------------------------------
        # When update_active_ID() draws an ID from the heap,
        # it discards "retracted" events where the T_next
        # value drawn from the heap is less than the T_next
        # value in the T_next grid.  So when we update the
        # T_next grid for a cell, we are retracting any
        # previously scheduled event for that cell.
        #-----------------------------------------------------
        self.update_T_next_grid( all_now_IDs, SILENT=SILENT, REPORT=REPORT)
        
        #----------------------------------------------
        # Update min and max values of various grids.
        # (1/25/12) This may have a significant cost,
        # especially in we use nanmin() and nanmax().
        #----------------------------------------------
        # self.update_mins_and_maxes()
        
        #------------------------
        # Check computed values
        #--------------------------------------------------
        # The method is inherited from erode_base.py, and
        # just checks if dz_max exceeds a limit (1000).
        #--------------------------------------------------
        OK = self.check_stability()

        #-------------------------------------------
        # Read from files as needed to update vars 
        #-----------------------------------------------------
        # NB! This is currently not needed for the "erosion
        # process" because values don't change over time and
        # read_input_files() is called by initialize().
        #-----------------------------------------------------
        # if (self.time_index > 0):
        #     self.read_input_files()

        #--------------------------------------------------------
        # NB! Since update_neighbors() is recursive,
        # putting this here causes some times to be
        # missed, etc.  So moved it to update_DEM(). (11/15/11)
        #--------------------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # Force writing grids at every step.
        # This will kill performance and should only
        # be used for testing and debugging.
        #----------------------------------------------
        ## self.write_output_files(0)       # (force every time)
        ## self.write_output_files()        # (BEFORE: 11/15/11)
        ## self.write_output_files( time )  # (AFTER:  11/15/11)

        if (OK):
            self.status = 'updated'  # (OpenMI 2.0 convention)
        else:
            self.status = 'failed'
            self.DONE   = True

        #------------------------
        # Update internal clock
        #------------------------
        self.update_time()
        ## print 'time_index =', self.time_index
        
        #-------------------------------------
        # Check if the model run is finished
        #-------------------------------------
        self.check_finished()    # (sets self.DONE)

        ####################################################
        # Run some diagnostic tests on zone_IDs (2/28/12)
        ####################################################
        ## self.check_zone_info()
                        
        #--------------
        # For testing
        #--------------
##        step = self.time_index
##        if ((step % 1) == 0):
##            print 'active_ID   =', self.active_ID
##            print 'size(IDs)   =', self.d8.IDs.size
##            print 'self.dt     =', self.dt, ' [years]'
##            print 'self.dz_max =', self.dz_max, ' [m]'
##            self.print_mins_and_maxes(step, DEM=True)
##            self.print_mins_and_maxes(step, dz_dt=True)
##            self.print_mins_and_maxes(step, A=True)
##            self.print_mins_and_maxes(step, S=True)
##            # self.print_mins_and_maxes(step, Q=True)
##            # self.print_mins_and_maxes(step, Qs=True)
##            print ' '
            
    #   update()
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

        #######################################################
        # (11/20/11) This string comparison may be costly.
        #######################################################
        ## if (self.time_step_type != 'fixed'):
        if not(self.FIXED_STEPS):
            #------------------------------------------------
            # (11/15/11). Save computed values based on the
            # elapsed time since last saved.  Is this the
            # best we can do to honor "save_grid_dt" ??
            #------------------------------------------------
            elapsed_grid_time = (time - self.last_grid_time)
            if (elapsed_grid_time > self.save_grid_dt):
                # print '#### Writing frame at time =', self.time
                self.save_grids()
                self.last_grid_time = time
            #-----------------------------------------------------
            elapsed_pixel_time = (time - self.last_pixel_time)
            if (elapsed_pixel_time > self.save_pixels_dt):
                self.save_pixel_values()
                self.last_pixel_time = time
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

    #   write_output_files()
    #-------------------------------------------------------------------
    def finalize(self):

        #--------------------------------------------------------------
        # Note: This overrides finalize() in erode_base.py. (10/4/11)
        #--------------------------------------------------------------
        self.status = 'finalizing'  # (OpenMI)   
        self.close_input_files()    # Input "data streams"
        self.close_output_files()

        ## if (self.mode == 'driver'):
        if (self.DRIVER):
            self.print_time_and_value(self.dz_max, 'dz_max', '[m]',
                                      interval=0.0 )
        self.update_entire_DEM()    ## (10/21/10)
        ###########################################
        # (11/8/10) We should also update fluxes
        # and the area grid, not just the DEM.
        ###########################################
        self.write_output_files(0)
        
        #-------------------------------
        # Save the "update_count" grid
        #-------------------------------
        self.save_update_count_grid()

        #--------------------------------------
        # Save all other final grids (9/2/11)
        #--------------------------------------
        self.save_final_grids()

        #----------------------------------------------
        # Added on 10/4/11 to help gauge performance.
        # See CSDMS_base.print_final_report.
        # This could eventually go there because it
        # would be useful to any model.  But not all
        # of them have self.time in "years".
        #----------------------------------------------
        finish         = time.time()
        run_time_secs  = (finish - self.start_time)
        sim_time_years = self.time
        years_per_sec  = (sim_time_years / run_time_secs)
        print(' ')
        print('Years per second =', years_per_sec)

        #------------------------------------------
        # Update mins and maxes of various grids.
        #------------------------------------------
        self.update_grid_mins_and_maxes()  ## (3/4/12)
        
        #------------------------------------------
        # Print report info for the Local version
        #------------------------------------------
        print(' ')
        print('CFL_factor =', self.CFL_factor)
        print(' ')
        print('Final heap size  =', self.n_heap)        
        print(' ')
        print('Max number of "now_IDs" =', self.max_n_now)
        print(' ')
        print('Max active_ID_count =', self.max_active_ID_count)
        print(' ')
##        print 'min(dt), max(dt) =', self.dt_min, self.dt_max, ' [yr]'
##        print ' '
        print('D = (Qs / delta_z).')
        print('D = 0 where D8 code = 0.')        
        print('min(D), max(D)   =', self.D_min, self.D_max, ' [m^2/yr]')
        print('min positive D-value =', self.D_min_pos)
        print(' ')
        print('S = 0 where D8 code = 0.')
        print('min(S),  max(S)  =', self.S_min,  self.S_max)
        print('min positive S-value =', self.S_min_pos)
        print(' ')
        print('Q = 0 where D8 code = 0.')
        print('min(Q),  max(Q)  =', self.Q_min,  self.Q_max,  ' [m^3/yr]')
        print('min positive Q-value =', self.Q_min_pos)
        print(' ')
        print('Qs = 0 where D8 code = 0.')
        print('min(Qs), max(Qs) =', self.Qs_min, self.Qs_max, ' [m^3/yr]')
        print('min positive Qs-value =', self.Qs_min_pos)
        print(' ')
        
        #----------------------------------------------
        # Print report info common with Global version
        #-----------------------------------------------
        print('min(dt), max(dt) =', self.dt_min, self.dt_max, ' [yr]')
        print(' ')
        print('min(z), max(z)   =', self.DEM_min, self.DEM_max, ' [m]')
        print(' ')
        print('dz_max_vec[0]    =', self.dz_max_vec[0])
        print('dz_max_vec[nt-1] =', self.dz_max_vec[-1])
        print('dz_max_vec.min() =', self.dz_max_vec.min())
        print('dz_max_vec.max() =', self.dz_max_vec.max())
        print(' ')
        self.print_final_report(comp_name='Erode-D8-Local 0.9 (3/4/12)',
                                mode='driver')  ## NEED THIS !!

        self.status = 'finalized'  # (OpenMI)

    #   finalize()
    #-------------------------------------------------------------------
    def update_active_ID(self, SILENT=True, REPORT=False):
        
        #------------------------------------------------------
        # Note:  This draws the next smallest update time and
        #        its corresponding ID from a priority queue.
        #        The new ID is called "active_ID".
        #------------------------------------------------------
        while (True):
            #----------------------------------------------
            # Use a heap and priority queue to find T_min
            #----------------------------------------------
            pair   = heapq.heappop( self.T_heap )
            min_ID = pair[1]
            T_min  = pair[0]
            self.n_heap -= 1
            #-----------------------------------------------
            # Skip over this ID and get the next smallest
            # T-value if this one was processed before its
            # scheduled time and retracted.
            #-----------------------------------------------
            # (3/2/12) Changed this test to use "<" which
            # may be needed due to roundoff errors.
            #-----------------------------------------------
            # RETRACTED = (T_min != self.T_next.flat[ min_ID ])
            RETRACTED = (T_min < self.T_next.flat[ min_ID ])
            #--------------------------------------
            # Comment this out to use next block.
            #--------------------------------------
            if not(RETRACTED):
                break

            #--------------------------------------------------
            # (3/3/12) This experiment did not appear to
            # succeed. For n=1, at pixels with IDs = 8029
            # and 2546, abrupt jumps occurred.  Both cells
            # started as higher points located within large
            # depressions, eventually "engulfed".  ID 8029
            # jumped to a peak when just before it seemed
            # to have similar elevation to neighboring cells.
            #
            # Is there any harm with drawing a pit or flat
            # from the heap?  Earlier updates generally seem
            # better than later updates.
            #--------------------------------------------------
            # (3/2/12) Experiment to disallow active_IDs
            # that have (d8_code == 0).  This should not
            # violate mass conservation because the cell
            # will be updated later on with same dz/dt.
            # This gets triggered a lot, so took out the
            # "Discarded pit" message.
            #--------------------------------------------------
            # NOTE: update_T_next() does not put cells on the
            # heap that have (d8_code == 0).  However, cells
            # drawn from the heap can have (d8_code == 0),
            # due to changes that occurred around them since
            # they were scheduled.  These must be pits and
            # flats because a cell cannot be changed to become
            # an edge cell.
            #--------------------------------------------------
##            d8_code = self.d8.d8_grid.flat[ min_ID ]
##            PIT_OR_FLAT = ( d8_code == 0 )
####            if (PIT_OR_FLAT):
####                print '### Discarded pit from heap:', min_ID
##            if not(RETRACTED) and not(PIT_OR_FLAT):
##                break

        #------------------------------------------------------
        # Cells with (d8_code == 0) were found to occur quite
        # often for 100x100 test on 2/29/12.  For example,
        # the cell with ID = 2375 and (col, row) = (75, 23)
        # is the higher of 2 cells in a 2-pixel pit where
        # (76, 23) is lower.  It starts with a nonzero D8
        # code, but then when the neighbor pit fills, it
        # becomes a pit with (d8_code == 0).
        #------------------------------------------------------        
        # If a pixel was not a pit when scheduled, then:
        # (1) It could have (dz/dt < 0).  This could cause
        #     update_DEM() to make it an even deeper pit.
        # (2) It could have (dz/dt > 0).  This could cause
        #     update_DEM() to make it a peak. ########
        #------------------------------------------------------
##        d8_code = self.d8.d8_grid.flat[ min_ID ]
##        if (d8_code == 0):
##            print '### WARNING: D8 code = 0 for active_ID =', min_ID
##            ## sys.exit()
            
        #----------------------------
        # Store active_ID and T_min
        #----------------------------
        self.active_ID = min_ID
        self.T_min     = T_min
        #--------------------------------------
##        row, col = divmod( min_ID, self.nx )
##        self.active_row = row
##        self.active_col = col

        #--------------------------------------
        # Update "active_ID count" (10/22/10)
        #--------------------------------------
        if (self.active_ID == self.last_active_ID):
            self.active_ID_count += 1
            self.max_active_ID_count = np.maximum( self.max_active_ID_count, \
                                                      self.active_ID_count )
        else:
            self.active_ID_count = 0
            self.last_active_ID = self.active_ID
            
        #--------------------------------------
        # Update "active_ID count" (10/26/10)
        # This version includes oscillations.
        #--------------------------------------
##        if (self.active_ID == self.last_active_ID) or \
##           (self.active_ID == self.last_last_active_ID):
##            self.active_ID_count += 1
##        else:
##            self.active_ID_count = 0
##            self.last_last_active_ID = self.last_active_ID
##            self.last_active_ID = self.active_ID
            
        #------------------
        # Optional report
        #------------------
        if (self.DEBUG):
        ## if (True):
            print(' time_index =', self.time_index)
            print(' Active ID  =', min_ID)
            print(' T_last[ID] =', self.T_last.flat[ min_ID ])
            print(' T_next[ID] =', self.T_next.flat[ min_ID ])
            print(' Heap size  =', self.n_heap)
            print('-------------------------------------')

        ### print ' Active ID =', min_ID  ################

    #   update_active_ID()
    #-------------------------------------------------------------------
    def update_T_clock(self, SILENT=True, REPORT=False):

        #--------------------
        # Update clock time
        #--------------------
        T_clock_last = self.T_clock.copy()
        self.T_clock = self.T_min

        #-----------------------------------------------------
        # Update the "clock dt", which is no longer the
        # same as the dt used to compute dz in update_DEM().
        # This dt is used for reporting elapsed time and by
        # update_base_level().
        #-----------------------------------------------------
        self.dt = (self.T_clock - T_clock_last)

        #----------------------------------------------------
        # (11/15/11) This happens a lot, but active_ID is
        # different each time.  So it just seems to mean
        # that many pixels are scheduled for the same time.
        # Order of processing them could affect result.
        #----------------------------------------------------
##        if (self.dt == 0.0):
##            print '===> dt = 0, active_ID =', self.active_ID

        #-----------------------------------------------------
        # (3/1/12) This should never happen.  It would mean
        # that there's a problem with Python's "heapq".
        # ##### Haven't tried this yet.
        #-----------------------------------------------------
        # According this this URL, there is some kind of
        # problem with Python's "heapq':
        #   http://brannerchinese.wordpress.com/2011/11/30/
        #          apparent-error-in-pythons-priority-queue/
        #-----------------------------------------------------        
        if (self.dt < 0):
            print('#####################################')
            print(' ERROR in update_T_clock():')
            print('       Time is going backwards!')
            print('       (Possible bug in heapq.)')
            print('       active_ID =', self.active_ID)
            print('       dt        =', self.dt)
            print('#####################################')
            print(' ')
            sys.exit()
            
        #------------------
        # Optional report
        #------------------
        REPORT = False
        if (REPORT):
            # print 'T_clock =', self.T_clock
            ID = self.active_ID
            row, col = divmod(ID, self.nx)
            S = self.S.flat[ ID ]
            A = self.d8.A.flat[ ID ]
            Qs = self.Qs.flat[ ID ]
            z = self.DEM.flat[ ID ]
            print('T_clock =', self.T_clock, ', (col,row, S, A, Qs, z) =', col, row, S, A, Qs, z)
            
##        print 'T_clock_last =', T_clock_last
##        print 'T_clock      =', self.T_clock
##        print 'dt           =', self.dt
##        print '----------------------------------------------'   
##        if (self.dt < 0):
##            print '##############  WARNING: dt < 0 !'
##            print 'dt =', self.dt
        
    #   update_T_clock()
    #-------------------------------------------------------------------
    def update_DEM(self, IDs, SILENT=True, REPORT=False):

        #------------------------------------------------------------
        # Notes: DEM   = current elevation grid  [m]
        #        Qs    = K * Q^m * S^n = sed. discharge [m^3/yr]
        #        da    = pixel area grid  [m^2]
        #        dt    = channel flow timestep  [s]  (returned)
        #        del_z = elevation drops to parent pixels [m]
        #        dz    = net elevation change due to sed. flux [m]
        #        U     = tectonic uplift rate [mm/year]
        #        w1    = IDs of D8 "flow-from" pixels
        #        p1    = IDs of D8 "flow-to" pixels

        #        NB!  Don't want elevations to change at pixels
        #        where base level is fixed, as with flow to sea.
        #------------------------------------------------------------
        # Where dz_dt <= 0 (eroding), no constraint on dt
        # is required for stability ??  (OLD COMMENT)
        #------------------------------------------------------------
        # Where dz_dt > 0 (depositing), dt must satisfy
        # CFL-type condition for stability.  (OLD COMMENT)
        #------------------------------------------------------------
        # Once z-value has been changed for the pixels in a set
        # of IDs, need to update the flow directions for the
        # pixels in the set as well as any of their D8 "kids".
        # Wherever the D8 flow code changes, the area grid of
        # that pixel's parent pixel must be updated.  For area
        # grid updates, use the recursive approach (starting
        # with pixels with no D8 kids and working downstream).
        # However, we only need to reinitialize the area grid
        # values for pixels where area could have changed; any
        # area grid values upstream of these pixels will still
        # be valid.  Slope grid values must also be updated for
        # any pixel whose D8 flow code has changed.  Once the
        # area and slope grids have been updated, the Q and Qs
        # grids then need to be updated.
        #------------------------------------------------------------
        if not(SILENT):    
            print('Updating DEM...')

        #--------------------------------------------------
        # This dt is NOT the elapsed time since ANY pixel
        # was updated but should instead be the elapsed
        # time since THIS pixel was updated.
        #--------------------------------------------------
        # Note that T_clock grows monotonically, but this
        # dt value could be very large if the cell has
        # been waiting on the heap for a while or was
        # "deactivated" (with dt = dt_limit) a while ago.
        #-------------------------------------------------
        ## dt = self.dt_grid.flat[ IDs ]  # (WRONG !)
        dt = (self.T_clock - self.T_last.flat[ IDs ])
        dz = self.dz_dt.flat[ IDs ] * dt
            
        ## print 'dz, dt =', dz[0], dt[0]
         
        #----------------
        # For debugging
        #----------------
        if (self.DEBUG):
        ## if (True):
        ## if (self.time_index > 660):
            zval = self.DEM.flat[ self.active_ID ]
            print(self.time_index, ':', self.active_ID, ':', self.T_clock, ': z=', zval)
            print('     dz =', dz, ',', 'dt =', dt)
            print(' ')

        ## print 'T_clock, T_last, ID =', self.T_clock, \
        ##        self.T_last.flat[IDs], IDs
        
        ## dz_dt = self.dz_dt.flat[ IDs ]
        ## print 'dt, dz, dz_dt =', dt, dz, dz_dt  #############
        
        #-----------------------------------------------
        # Impose a "minimum allowed dz" value
        # otherwise eventually dt goes to 0? (8/25/10)
        #-----------------------------------------------
        # (10/20/10) Tried again but not useful ?
        #-----------------------------------------------        
        ## dz = np.maximum( self.min_allowed_dz, dz )
        ## dt = np.maximum( self.min_allowed_dt, dt )
        
        #------------------------
        # Update the DEM at IDs
        #------------------------
        self.DEM.flat[ IDs ] += dz

        #------------------------------------------------------
        # (11/18/10) Count number of times each ID is updated
        #------------------------------------------------------
        self.update_count.flat[ IDs ] += 1
        
        #------------------------------------------------------
        # Raise pit parent to level of active ID ?? (11/4/10)
        #------------------------------------------------------
        # NB! The get_interior_pit_IDs() function only
        # returns pits that are included in current IDs !!
        # So it returns empty array for ID=137, even though
        # its parent (pID = 136) is a pit.
        #------------------------------------------------------
        # Filling the pit abruptly like this may alter the
        # overall evolution of the landscape.
        #------------------------------------------------------        
####        if (self.active_ID == self.last_active_ID):
##        ID = self.active_ID
##        if (ID in IDs):      
##            pID = self.d8.parent_ID_grid.flat[ ID ]    
##            ## if (self.is_a_pit(pID)):
##            if (self.is_an_interior_pit(pID)):
####                print '##########################################'
####                print ' Raising pit elevation to:', self.DEM.flat[ ID ]
####                print ' ID, pID    =', ID, ', ', pID
####                print ' time_index =', self.time_index
####                print ' Initial z[ID] =', self.DEM.flat[ ID ]
####                print ' z[pID] =', self.DEM.flat[ pID ]
####                print '##########################################'
##                ## self.DEM.flat[ pID ] = self.DEM.flat[ ID ]
##                self.DEM.flat[ pID ] = self.DEM.flat[ ID ] + 1e-8   #########

                    
##        if (self.time_index >= 280):
##            self.DEBUG  = True
##            self.SILENT = False
            
##        print 'IDs =', IDs
##        print 'New DEM values =', self.DEM.flat[ IDs ]
##        print ' '

##        if (32 in IDs):
##            print '############## Updated z[32] at:', self.T_clock
            
##        if (32 in IDs) or (51 in IDs):
##        if (self.time_index > 580):
##            pIDs = self.d8.parent_ID_grid.flat[ IDs ]
####            dz_dt_pID = self.dz_dt.flat[ pIDs ]
####            dz_dt_ID  = self.dz_dt.flat[ IDs]
##            print '---------------------------------------------'
##            print 'active_ID  =', self.active_ID
##            print 'time_index =', self.time_index
##            print 'T_clock    =', self.T_clock
##            print 'IDs   =', IDs
##            print 'pIDs  =', pIDs
##            print 'z[12] =', self.DEM.flat[ 12 ]
##            print 'z[32] =', self.DEM.flat[ 32 ]
##            print 'z[51] =', self.DEM.flat[ 51 ]
##            print 'dz_dt[12] =', self.dz_dt.flat[ 12 ]
##            print 'dz_dt[32] =', self.dz_dt.flat[ 32 ]
##            print 'dz_dt[51] =', self.dz_dt.flat[ 51 ]
##            print 'dz    =', dz
##            print 'dt    =', dt
##            print ' '

##        if (self.time_index > 706):
##            pIDs = self.d8.parent_ID_grid.flat[ IDs ]
####            dz_dt_pID = self.dz_dt.flat[ pIDs ]
####            dz_dt_ID  = self.dz_dt.flat[ IDs]
##            print '---------------------------------------------'
##            print 'active_ID  =', self.active_ID
##            print 'time_index =', self.time_index
##            print 'T_clock    =', self.T_clock
##            print 'IDs   =', IDs
##            print 'pIDs  =', pIDs
##            print 'z[12] =', self.DEM.flat[ 12 ]
##            print 'z[32] =', self.DEM.flat[ 32 ]
##            print 'z[51] =', self.DEM.flat[ 51 ]
##            print 'dz_dt[12] =', self.dz_dt.flat[ 12 ]
##            print 'dz_dt[32] =', self.dz_dt.flat[ 32 ]
##            print 'dz_dt[51] =', self.dz_dt.flat[ 51 ]
##            print 'dz    =', dz
##            print 'dt    =', dt
##            print ' '
        
        #-------------------------
        # Update DEM min and max
        #-------------------------
        zmin = self.DEM.flat[ IDs ].min()
        zmax = self.DEM.flat[ IDs ].max()
        #--------------------------------------------------
        # It seems that nanmin() and nanmax() are pretty
        # slow.  If everything is working, shouldn't need
        # to worry about NaNs. (1/25/12)
        #--------------------------------------------------
        # zmin = np.nanmin( self.DEM.flat[ IDs ] )
        # zmax = np.nanmax( self.DEM.flat[ IDs ] )
        
        #-------------------------------------------
        # Check whether zmax ever increases, which
        # should never occur (10/14/11)
        #-------------------------------------------
        if (zmax > self.DEM_max):
            w = np.where( self.DEM.flat[ IDs ] == zmax )  # (or argmax)
            max_ID  = IDs[ w[0] ]
            max_pID = self.d8.parent_ID_grid.flat[ max_ID ]
            del_z   = self.DEM.flat[ max_ID] - self.DEM.flat[ max_pID ]
            dt_here = self.T_clock - self.T_last.flat[ max_ID ]
            dz_here = self.dz_dt.flat[ max_ID ] * dt_here
            print('########################################')
            print(' ERROR: zmax just increased.')
            print('     old_zmax   =', self.DEM_max)
            print('     new_zmax   =', zmax)
            print('     time_index =', self.time_index)
            print('     ID         =', max_ID)
            print('     del_z[ID]  =', del_z)
            print('     Qs[ID]     =', self.Qs.flat[ max_ID ])
            print('     dz_dt[ID]  =', self.dz_dt.flat[ max_ID ])
            print('     dt[ID]     =', dt_here)
            print('     dz[ID]     =', dz_here)
            print('     dz_targ[ID]=', self.dz_target_grid.flat[ max_ID ])
            print('########################################')            
        self.DEM_min = np.minimum( self.DEM_min, zmin )
        self.DEM_max = np.maximum( self.DEM_max, zmax )
        

        #-----------------------------------------------------------
        # Compute dz and the largest dz obtained during this
        # timestep. This is used to check stability & convergence.
        #-----------------------------------------------------------
        # (4/14/10) Needed to add "abs" because otherwise
        #           simulation is aborted if dz_max < 0. 
        #--------------------------------------------------
        abs_dz      = np.absolute( dz )     ## (9/30/10)
        # self.dz_max = np.nanmax( abs_dz )   ## (9/30/10)
        self.dz_max = abs_dz.max()   # (1/25/12. See Notes about nanmax().)
        ## self.dz_max = np.nanmax(dz)
        ## self.dz_max = np.absolute( self.dz_max )
        if (self.stop_code == 0):
            self.dz_max_vec[ self.time_index - 1 ] = self.dz_max

        #------------------------------------
        # This should never occur (9/2/10).
        #------------------------------------
        if (self.DEBUG and (self.dz_max == 0)):
            print('ERROR: dz_max = 0 at time index:', self.time_index)
            ## print '       IDs =', IDs
            ## sys.exit()
        ## print 'dz_max =', self.dz_max   ####################
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            self.dz_min = np.nanmin( dz )
            dz_str = str(self.dz_min)   + ', ' + str(self.dz_max)
            print('    min(dz), max(dz) = ' + dz_str + ' [m]')


        #-----------------------------------------------------------
        # Experiment: 10/6/11. Update D8 codes and slopes whenever
        # DEM is updated, but only update fluxes Q and Qs when DES
        # algorithm says to do so.  Haven't used this yet, though,
        # because looking at how update_fluxes2() works, I don't
        # think it can be done here.
        #-----------------------------------------------------------
        # We are trying to address a bug where dz_dt remains > 0
        # for a peak pixel.  Since update_dz_dt_grid() uses D8
        # flow codes of all 8 neighbors, maybe they need to be
        # updated after every call to update_DEM() instead of only
        # in update_fluxes2().
        #-----------------------------------------------------------
##        self.update_d8_vars( IDs, SILENT=SILENT, REPORT=REPORT)
##        self.update_slope_grid( IDs, SILENT=SILENT, REPORT=REPORT)
        
    #   update_DEM()
    #-------------------------------------------------------------------
    def update_entire_DEM(self, SILENT=True, REPORT=False):

        #------------------------------------------------------------
        # Notes: DEM   = current elevation grid  [m]
        #        Qs    = K * Q^m * S^n = sed. discharge [m^3/yr]
        #        da    = pixel area grid  [m^2]
        #        dt    = channel flow timestep  [s]  (returned)
        #        del_z = elevation drops to parent pixels [m]
        #        dz    = net elevation change due to sed. flux [m]
        #        U     = tectonic uplift rate [mm/year]
        #        w1    = IDs of D8 "flow-from" pixels
        #        p1    = IDs of D8 "flow-to" pixels

        #        NB!  Don't want elevations to change at pixels
        #        where base level is fixed, as with flow to sea.
        #------------------------------------------------------------
        # Where dz_dt <= 0 (eroding), no constraint on dt
        # is required for stability ??
        #------------------------------------------------------------
        # Where dz_dt > 0 (depositing), dt must satisfy
        # CFL-type condition for stability.
        #------------------------------------------------------------
        # Once z-value has been changed for the pixels in a set
        # of IDs, need to update the flow directions for the
        # pixels in the set as well as any of their D8 "kids".
        # Wherever the D8 flow code changes, the area grid of
        # that pixel's parent pixel must be updated.  For area
        # grid updates, use the recursive approach (starting
        # with pixels with no D8 kids and working downstream).
        # However, we only need to reinitialize the area grid
        # values for pixels where area could have changed; any
        # area grid values upstream of these pixels will still
        # be valid.  Slope grid values must also be updated for
        # any pixel whose D8 flow code has changed.  Once the
        # area and slope grids have been updated, the Q and Qs
        # grids then need to be updated.
        #------------------------------------------------------------
        if not(SILENT):    
            print('Updating entire DEM...')

        #--------------------------------------------------
        # This dt is NOT the elapsed time since ANY pixel
        # was updated but should instead be the elapsed
        # time since THIS pixel was updated.
        #--------------------------------------------------
        # Note that T_clock grows monotonically, but this
        # dt value will be very large if the cell has
        # been waiting on the heap for a while.
        #-------------------------------------------------
        ## dt = self.dt_grid  # (WRONG !)  
        dt = (self.T_clock - self.T_last)
        dz = self.dz_dt * dt
        
        #------------------------
        # Update the DEM at IDs
        #------------------------
        self.DEM += dz
        
        #----------------------------------------------------
        # (9/2/11) Count number of times each ID is updated
        #----------------------------------------------------
        self.update_count += 1    # see update_DEM()
        
        #-------------------------
        # Update DEM min and max
        #-------------------------
        zmin = self.DEM.min()
        zmax = self.DEM.max()
##        zmin = np.nanmin( self.DEM )  # (nanmax is slow, 11/25/12)
##        zmax = np.nanmax( self.DEM )
        self.DEM_min = np.minimum( self.DEM_min, zmin )
        self.DEM_max = np.maximum( self.DEM_max, zmax )
             
        #-----------------------------------------------------------
        # Compute dz and the largest dz obtained during this
        # timestep. This is used to check stability & convergence.
        #-----------------------------------------------------------
        # (4/14/10) Needed to add "abs" because otherwise
        #           simulation is aborted if dz_max < 0. 
        #--------------------------------------------------
        abs_dz      = np.absolute( dz )     ## (9/30/10)
        ## self.dz_max = np.nanmax( abs_dz )   ## (9/30/10)
        self.dz_max = abs_dz.max()
        ## self.dz_max = np.nanmax(dz)
        ## self.dz_max = np.absolute( self.dz_max )
        if (self.stop_code == 0):
            self.dz_max_vec[ self.time_index - 1 ] = self.dz_max

        #------------------------------------
        # This should never occur (9/2/10).
        #------------------------------------
        if (self.DEBUG and (self.dz_max == 0)):
            print('ERROR: dz_max = 0 at time index:', self.time_index)
            ## print '       IDs =', IDs
            ## sys.exit()
        ## print 'dz_max =', self.dz_max   ####################
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            self.dz_min = np.nanmin( dz )
            dz_str = str(self.dz_min)   + ', ' + str(self.dz_max)
            print('    min(dz), max(dz) = ' + dz_str + ' [m]')
           
    #   update_entire_DEM()
    #-------------------------------------------------------------------    
    def update_neighbors(self, IDs, SILENT=True, REPORT=False):

        #------------------------------------------------------------
        # Note: Maybe rename this to "update_zone_IDs()" ?
        #------------------------------------------------------------
        # O&K(2006) paper says to call the event synchronization
        # procedure for all neighbors of IDs. This involves:
        #    (0) find the neighbors (we do several at once)
        #    (1) compute dz = dz/dt * (T_clock - T_last)
        #    (2) add dz to flux_capacitor (dz_cap)
        #    (3) add dz to z
        #    (4) set T_last = T_clock
        #    (5) determine which neighbor cells were affected
        #        enough that they need to be updated NOW.
        #        These have (abs_dz_cap >= dz_target), but we
        #          include pixels with D8 codes of 0 also.
        #    (6) for "now_IDs" that must be updated now:
        #        (a) retract the pending event
        #        (b) zero out the flux capacitor (dz_cap)
        #        (c) update all neighbors (of now_IDs)
        #        (d) schedule a new event
        #            (i)   update fluxes between IDs and now_IDs ??
        #            (ii)  update dz_dt (for now_IDs)
        #            (iii) update dt and dz_target
        #            (iv)  update T_next
        #    (7) for "wait_IDs" that don't need to be updated yet:
        #        (a) update z-values (implied by (1))
        #        (b) update fluxes between IDs and wait_IDs, using
        #            up-to-date z-values for IDs and wait_IDs
        #        (c) update dz/dt for wait_IDs
        #
        # Note: Boundary cells always have dz/dt = BLR (base-level
        # lowering rate.  Since they have D8 codes of 0, they can
        # only be updated if one of their neighbors is updated.
        # They also have A, S, Q and Qs all equal to zero.
        #------------------------------------------------------------
        
        #------------------------------------------------------
        # Find the set of unique neighbors to the set of IDs
        #------------------------------------------------------
        # See unit test called: d8_local.get_neighbors_test()
        #------------------------------------------------------     
        nbr_IDs = self.d8.get_neighbors( IDs )

        #----------------------------------------------------
        # Remove nbr_IDs already in self.zone_IDs (2/14/12)
        # Should no longer need to use unique() on zone_IDs
        # or all_now_IDs in update().
        #----------------------------------------------------
        # setdiff1d returns sorted, unique nbr_IDs that are
        # not in self.zone_IDs.
        #----------------------------------------------------
        nbr_IDs = np.setdiff1d( nbr_IDs, self.zone_IDs )
        if (nbr_IDs.size == 0):
            return

##        print 'nbr_IDs ='
##        print nbr_IDs
##        print ' '
        
##        print 'nbr_IDs.shape =', nbr_IDs.shape
##        print 'nbr_IDs.size  =', nbr_IDs.size

        
        ##########################################
        # Maybe we should run this test again.
        ##########################################
        if (self.DEBUG):
            print('#### In update_neighbors():')
            print('#### start IDs    =', IDs)
            print('#### neighbor IDs =', nbr_IDs)
            
        #---------------------------------
        # Append new nbr_IDs to zone_IDs
        #---------------------------------
        self.append_zone_IDs( nbr_IDs )

        #-------------------------------------------------
        # Do this for all neighbors, even if on boundary
        #-------------------------------------------------
        # Note that A, S, Q, Qs all equal 0 on boundary,
        # and dz/dt equals base-level lowering rate.
        #------------------------------------------------------
        # NB! update_DEM() uses T_last and T_clock to compute
        # dz.  So if we pick up previously processed IDs when
        # update_neighbors() calls itself, that should be OK
        # because the dt and computed dz will both be zero.
        # (This note no longer applies. 2/28/12)
        #--------------------------------------------------------      
        self.update_DEM( nbr_IDs, SILENT=SILENT, REPORT=REPORT )
  
        #-------------------------------------------------
        # Update the flux capacitor and last update time
        #------------------------------------------------------------
        # The flux capacitor grid stores the sum of all the dz's
        # that have occurred since the last "full" update when the
        # cell was included in now_IDs.  If it gets as big as
        # dz_target, then that cell is included in now_IDs and the
        # flux capacitor gets reset to 0.
        #------------------------------------------------------------        
        # NB! This should only be called when synchronizing a
        #     neighbor cell, so cannot be merged into update_DEM().
        #------------------------------------------------------------
        self.update_flux_cap_grid( nbr_IDs, SILENT=SILENT, REPORT=REPORT )
        self.update_T_last_grid( nbr_IDs, SILENT=SILENT, REPORT=REPORT )

        #-----------------------------------------------------
        # (9/14/11) Both dz_target and dz_cap carry signs.
        # If they have the *same* sign AND if abs(dz_cap)
        # >= abs(dz_target), then we need to update z-values
        # now vs. at the future, scheduled time.  If they
        # have *opposite* signs, then we should probably
        # also consider them "affected" and update them now.
        #
        # (10/6/11) This may be the cause of a bug that
        # allows a peak to have dz/dt > 0.  That is, it must
        # have had D8 kids when it was scheduled in order to
        # get a dz/dt > 0.  But then it became a peak, with
        # (dz_cap > 0) and (dz_target < 0).  If abs(dz_cap)
        # < abs(dz_target), then it would not be updated now
        # and when drawn from heap the peak would grow since
        # dz/dt > 0.  Can something similar happen for a
        # pixel that starts as a peak (with dz/dt < 0) and
        # then becomes a pixel with D8 kids and (dz/dt > 0)?
        # This situation seems less likely.
        #
        #########################################################
        #
        # (2/28/12) We still have a rare problem where a pit
        # suddenly becomes a peak with an adjacent pit.
        # Recall that pits are never drawn from the heap and
        # only get updated by neighbor updates.  As a pit,
        # the cell would have (dz_target = 0) (by convention),
        # (dz/dt > 0) and (d8_code = 0).
        #
        # Is it possible that the neighbor who becomes the pit
        # adjacent to the new peak is causing the problem?
        # Maybe it gives a large contribution to its (formerly)
        # pit neighbor.
        #
        # Pits in nbr_IDs, which have (d8_code == 0), will
        # usually have (dz_cap > 0) from the call above to
        # update_flux_cap_grid(), but also have their dz_targ
        # set to zero by initialize_dz_target_grid2() and by
        # update_dz_target_grid2(). Edges are similar, except
        # their dz_targ is NOT set to 0 by initialize_dz_target_grid2().
        #
        # Even if a pit has (dz_cap == 0) for some reason, we
        # may want to add it to now_IDs so that *its* neighbors
        # can be checked and possibly added to now_IDs.
        #
        #########################################################

        #-----------------------------------------------------
        # Determine which neighbors need to be updated NOW.
        #-----------------------------------------------------
        # Could there be other neighbors that need to be
        # updated now, such as those for which D8 code or
        # D8 area has changed ??  (2/28/12)
        ######################################################
        dz_cap      = self.flux_cap_grid.flat[ nbr_IDs ]
        dz_targ     = self.dz_target_grid.flat[ nbr_IDs ]
        abs_dz_cap  = np.absolute( dz_cap )
        abs_dz_targ = np.absolute( dz_targ )
        sgn_dz_cap  = np.sign( dz_cap )    # (10/6/11)
        sgn_dz_targ = np.sign( dz_targ ) 

        #------------------------------------------------------------
        # (2/29/12) check_zone_info() showed that there are cases
        # when (d8_code == 0) and (dt != dt_limit).  Because of
        # how update_dz_target_grid2() works, this means that they
        # did not get assigned (dz_target = 0) and did not get
        # included with now_IDs.
        #------------------------------------------------------------        
        # test1 = (abs_dz_targ == 0)    # (D8 code equals zero)
        
        #---------------------------------------------------------
        # (2/28/12) All cells with (d8_code == 0), except edges
        # should be included in "now_IDs". If edges are included
        # then first edge cell encountered causes all edge cells
        # to be added to zone_IDs.
        #
        #    => test1 is True for all pits/flats.
        #    => test1 is False for edges
        #---------------------------------------------------------
        # (3/7/12) Leaving this test in, even though we are now
        # assigning dt-values to pits/flats and putting on heap.
        #---------------------------------------------------------        
        d8_codes = self.d8.d8_grid.flat[ nbr_IDs ]
        not_edge = self.not_edge_grid.flat[ nbr_IDs ]
        test1    = np.logical_and(d8_codes == 0, not_edge)       

        #------------------------------------------------------------
        # All cells that have changed more than the amount they
        # were predicted to change at their next scheduled update
        # time should be included in "now_IDs".  dz_target is set
        # to zero for pits and flats and very large for edges.
        #
        #    => test2 is True  for pits/flats with (dz_cap != 0).
        #    => test2 is False for edges.  Since (dz_dt == -BLR)
        #          on edges, we have (dz_targ < dz_cap < 0).
        #------------------------------------------------------------        
        test2 = (abs_dz_cap >= abs_dz_targ) # (predictor size wrong)
        NOW   = np.logical_or( test1, test2 )
        
        #------------------------------------------------------------
        # All cells that have changed by an amount that has a
        # different sign than was predicted to occur at their next
        # update time should be included in "now_IDs".  This picks
        # up some cells with (abs_dz_cap < abs_dz_targ) that did
        # not pass test2.  But because of the sign difference,
        # the prediction is off by even more than dz_targ.
        # Keep in mind also that sign(0) = 0.
        #
        #    => test3 is True for pits/flats with (dz_cap != 0).
        #    => test3 is False for edges.  Since (dz_dt == -BLR)
        #          on edges, we have (dz_targ < dz_cap < 0).
        #------------------------------------------------------------        
        test3 = (sgn_dz_cap != sgn_dz_targ) # (predictor sign wrong)
        NOW   = np.logical_or( NOW, test3 )


##        T_last  = self.T_last.flat[ nbr_IDs ]
##        dt_vals = self.dt_grid.flat[ nbr_IDs ]
##        test3B  = (self.T_clock - T_last) < dt_vals
##        NOW     = np.logical_and( NOW, test3B )
        
        #----------------------------------------------------------
        # Note that test1, test2 and test3 are all True for pits
        # and flats and are all False for edges.  So even with
        # test4, edges will never be included in "now_IDs".
        #---------------------------------------------------------- 
        # Neighbors with (dz_cap == 0) should not be included in
        # "now_IDs".  This happens because it is fairly common for
        # different active_IDs to be scheduled for the same update
        # time. (See update_T_clock().)  This, in turn, happens
        # because when slope-exponent, n, equals 1, then dt_equal
        # is a function of A and ds.  Therefore, any cells with
        # the same A and same ds get the same dt-value.  This
        # issue seems to be related to the one-pixel instability,
        # because it doesn't happen if we take out test4, even
        # though zone_IDs then becomes very large.
        #
        # (3/2/12) I ran a test with n=1.1, which should result
        # in far fewer identical dt-values.  But it still had the
        # one-pixel instability (at least twice!).
        #
        # It seems that test4 is excluding some cells from
        # "now_IDs" that should be included.  Maybe their only
        # path of connectivity to the splash zone is through a
        # pixel that has (dz_cap == 0) ??
        #
        # Note that dz_cap was just updated above, so a cell that
        # has (dz_cap == 0) at this point must either have
        # (dz_dt == 0) OR (T_last == T_clock).  The latter can
        # only happen if it was in the splash zone of another
        # active_ID with the same T_next as the current one.
        # But it could have been a "wait_ID" vs. a "now_ID" and
        # it would still have (T_last == T_clock) now.  If it
        # was a "wait_ID", then it could have a nonzero dz_cap
        # from its last update, but this would not have been
        # augmented by update_flux_cap_grid() since it now has
        # T_last = T_clock.
        #
        # Will anything *new* happen if we include a pixel with
        # (dz_cap == 0) in now_IDs again?  Apparently so, because
        # we didn't encounter one-pixel instabilities.
        #
        # Any cell in the splash zone of a previous active_ID that
        # happens to have the same update time as the current
        # active_ID will still have (T_last == T_clock).  This
        # causes dz_cap to be set to 0 in update_flux_cap_grid().
        # Without this test, zone_IDs can grow to a large size
        # and exceed its max size.  See initialize_zone_IDs().
        #
        # We *could* process all active_IDs with the same update
        # time as a group.  But they are often spatially disjoint
        # and happen because dt = f(A,ds), as noted above.
        #----------------------------------------------------------
        # NB!!!  This test uses AND.  All others used OR.
        #-----------------------------------------------------------
        # (3/3/12) Tried test4 using T_last vs. dz_cap.  It seemed
        # to run faster (about 1500 steps per 5 secs vs. 1300) but
        # encountered "zmax" error early on at time_index = 43787.
        # Does this mean some cells have (dz/dt == 0) ??
        #-----------------------------------------------------------        
        ## test4 = (abs_dz_cap > 0)   # (nonzero change in z)
        test4 = (dz_cap != 0)
        ## T_last = self.T_last.flat[ nbr_IDs ]  # (see note above)
        ## test4  = (T_last != self.T_clock)   # (just updated)
        NOW   = np.logical_and( NOW, test4 )

        #-----------------------------------------------
        # (3/5/12) Could use this to ensure that edges
        # stay excluded from all_now_IDs.  Then test1
        # could just be (d8_codes == 0).
        #-----------------------------------------------
        ## NOW   = np.logical_and( NOW, not_edge )


        
##        if (self.time_index > 54170):
##            w = np.where( self.zone_IDs != -1 )
##            print '======================================================='
##            print 'zone_IDs =', self.zone_IDs[ w ]
##            print 'nbr_IDs  =', nbr_IDs
##            print 'd8_codes =', d8_codes
##            print 'not_edge =', not_edge
##            print 'dz_cap   =', dz_cap
##            print 'dz_targ  =', dz_targ
##            print 'dz_dt    =', self.dz_dt.flat[ nbr_IDs ]
##            print 'T_clock  =', self.T_clock
##            print 'T_last   =', self.T_last.flat[ nbr_IDs ]
##            print 'dt_vals  =', self.dt_grid.flat[ nbr_IDs ]
##            print ' '
            
        #---------------------------------------------------
        # We may be able to avoid this call to where. Note
        # that array[NOW] = array[w1] in Python.  However,
        # we need to know w1[0].size, etc.  (2/16/12)
        #---------------------------------------------------
        ### w1 = NOW  ## ????
        w1 = np.where( NOW )    # (update these now)

        #-------------------------------------
        # NO LONGER NEED wait_IDs. (2/14/12)
        #-------------------------------------
##        WAIT  = np.invert( NOW )
##        w2    = np.where( WAIT )   # (update these later)


        #------------------------------------------------
        # Should be equivalent to the 10/6/11 approach.
        #------------------------------------------------
        # NOW   = np.logical_or( test1, test2 )
        # NOW   = np.logical_and( NOW, abs_dz_cap > 0 )
        # WAIT  = np.invert( NOW )
        # w1 = np.where( NOW )
        # w2 = np.where( WAIT )
        
        #-------------------------------
        # Approach used until: 10/6/11
        #---------------------------------------------------
        # Note that if (abs_dz_cap == 0) then capacitor is
        # empty and no need to update now.  The case of
        # (abs_dz_targ == 0) is handled correctly, but is
        # a bit confusing as written.
        #---------------------------------------------------
        # w1 = np.where( logical_and(abs_dz_cap > 0,
        #                               abs_dz_cap >= abs_dz_targ ) )
        # w2 = np.where( logical_or(abs_dz_cap == 0,
        #                              abs_dz_cap  < abs_dz_targ ) )

        #--------------------
        # Original approach
        #--------------------
        # w1 = np.where( abs_dz_cap >= abs_dz_targ )
        # w2 = np.where( abs_dz_cap  < abs_dz_targ )

        
        #-----------------------------------------------------------
        # Recursive call to update neighbors of affected neighbors
        #-----------------------------------------------------------
        # O&K(2006) paper says that for neighbors that must be
        # updated now, do the following:
        #     (1) retract the pending event
        #     (2) zero out the flux capacitor
        #     (3) update all neighbors (of these neighbors)
        #     (4) schedule a new event
        #         (a) update fluxes between IDs and now_IDs ??
        #         (b) update dz_dt (for original IDs)
        #         (c) update dt and dz_target
        #         (d) update T_next
        # In our implementation, a pending event is retracted by
        # changing its scheduled time in T_next.  So (1) is part
        # of (4), but the order doesn't matter because we only
        # call update_active_ID() from update().
        #-----------------------------------------------------------
        n_NOW = w1[0].size
        if (n_NOW > 0):
            now_IDs = nbr_IDs[ w1 ]  # (must be updated now)

            #---------------------------------
            # Append now_IDs to all_now_IDs.
            #---------------------------------
            self.append_all_now_IDs( now_IDs )
            
            #--------------------------------------------
            # Flux cap grid values were only used to
            # find now_IDs, so can set it back to zero.
            # They should *not* be reset to zero for
            # the "wait_IDs" (i.e. nbrs not in now_IDs).
            #--------------------------------------------
            # NB! Look carefully at whether it is okay
            # to set dz_cap = 0 here, instead of after
            # the call to update_neighbors().
            #############################################
            self.flux_cap_grid.flat[ now_IDs ] = 0
            self.update_neighbors( now_IDs )

        #-------------------------------------------------------------
        # O&K(2006) paper says that for neighbors that do NOT need
        # to be updated now ("slightly affected"), do the following:
        #    (1) update z-values for IDs and wait_IDs
        #    (2) update fluxes between IDs and wait_IDs using
        #        just updated z-values at IDs and wait_IDs
        #    (3) update dz/dt for wait_IDs using
        #        just-updated fluxes between IDs and wait_IDs.
        # Note: These neighbors keep their old update time, T_next.
        # Note: In order to update the fluxes Q and Qs, we must
        #       first update D8 codes, A and S.
        #-------------------------------------------------------------
        # NO LONGER NEED THE wait_IDs. (2/14/12)
        #-------------------------------------------------------------
##        n_WAIT = w2[0].size
##        if (n_WAIT > 0):
##            wait_IDs = nbr_IDs[ w2 ]   # (IDs to update later)
##
##            #-----------------------------------
##            # Append wait_IDs to all_wait_IDs.
##            # It is initialized in update().
##            #-----------------------------------
##            j = self.all_wait_IDs_index
##            self.all_wait_IDs[ j: j + n_WAIT ] = wait_IDs  # (top_index is correct)
##            self.all_wait_IDs_index += n_WAIT

    #   update_neighbors()
    #-------------------------------------------------------------------
    def update_fluxes(self, IDs, nbr_IDs, SILENT=True, REPORT=False):

        #---------------------------------------------------
        # Note: This has been replaced by update_fluxes2()
        #       and may now be obsolete.
        #---------------------------------------------------
        
        #-----------------------------------------------
        # Combine IDs and their nbr_IDs into one group
        #-----------------------------------------------
        IDs2 = np.concatenate( (IDs, nbr_IDs) )
            
        #--------------------------------------------------------
        # Update the D8 flow grid and all vars that depend on
        # it, including the D8 area grid, A, and slope grid, S.
        #--------------------------------------------------------
        # The set of cells where the area grid changes is a
        # superset of "IDs";  it includes other, downstream
        # cells.  These are saved by the d8_update_area_grid()
        # function as d8.new_A_IDs.
        #---------------------------------------------------------
        # The D8 flow grid gives the directions of fluxes and
        # Q=Q(A) and Qs=Qs(A,S) give the magnitudes.
        #--------------------------------------------------------      
        self.update_d8_vars( IDs2, SILENT=SILENT, REPORT=REPORT)
        self.update_slope_grid( IDs2, SILENT=SILENT, REPORT=REPORT)

        #---------------------------------------------------------
        # (11/3/10) With A_IDs = IDs, stalls at ID = 311.
        #           With A_IDs = new_A_IDs, stalls at 131.
        # ID=131 is a pit.  ID=311 is in 2-pixel pit with 290.
        #---------------------------------------------------------
        # (11/18/10) If we use A_IDs instead of new_A_IDs, the
        # speed and final results are now very similar.
        # Differences when we use A_IDs include:
        #  (1) No more cases of stalling and deactivation.
        #  (2) Overall dt_min goes from 0.0009627 to 86.459 years.
        #  (3) Average number of updates per cell: 253 vs. 247
        #  (4) Final heap size = 1003 vs. 1007.
        #  (5) min(z), max(z) = -17.4571, 2.65324 vs.
        #                       -17.3981, 2.6528
        #---------------------------------------------------------
        # (9/12/11). If we use new_A_IDs, we seem to run the
        # risk of changing a Q or Qs-value before it has been
        # used, which may cause a mass balance problem.  But
        # maybe not, since dz_dt is used for z-value updates.
        #---------------------------------------------------------
        # (9/15/11) Tried new_A_IDs again and got a lot of
        # "Progress has stalled" error messages.
        #---------------------------------------------------------
        ### A_IDs = self.d8.new_A_IDs
        A_IDs = IDs2

        #--------------------------------------------
        # Update fluxes Q and Qs for cells in A_IDs
        #--------------------------------------------
        self.update_Q_grid(  A_IDs, SILENT=SILENT, REPORT=REPORT)
        self.update_Qs_grid( A_IDs, SILENT=SILENT, REPORT=REPORT)

        #-----------------------------------------------------
        # Update the diffusion coefficient, D = (Qs / del_z)
        #-----------------------------------------------------
        # Added back in on (11/14/11), with new version.
        #-----------------------------------------------------
        self.update_D_grid( flux_IDs, SILENT=SILENT, REPORT=REPORT)

    #   update_fluxes()
    #-------------------------------------------------------------------
    def update_fluxes2(self, IDs, nbr_IDs, SILENT=True, REPORT=False):

        #----------------------------------------------------------
        # Note: This version only updates fluxes Q and Qs for
        #       cells that have flow BETWEEN IDs and nbr_IDs.
        #       So if a cell in IDs flows to another in IDs,
        #       their fluxes will NOT be updated and similarly
        #       for a cell in nbr_IDs that flows to another
        #       in nbr_IDs.  It is not yet clear whether this
        #       will be different than update_fluxes(). (9/16/11)
        #----------------------------------------------------------
        #       Note that "in1d()" function is new in NumPy 1.4.
        #----------------------------------------------------------
        #       In 2D, there will typically be *no* flux between
        #       a cell and one of its 8 neighbors.  So a change
        #       in flux between them only occurs if:
        #       (1) D8 neighbor used to flow to ID, but not now.
        #       (2) D8 neighbor didn't used to, but does now.
        #       (3) ID used to flow to D8 neighbor, but not now.
        #       (4) ID didn't used to, but does now.
        #----------------------------------------------------------

        ###########################################
        #    EXPERIMENT (11/14/11)
        ###########################################
        ## self.update_fluxes3( IDs, nbr_IDs )
        ## return
        
        #--------------------------------------------
        # Find IDs that flow to nbr_IDs and nbr_IDs
        # that flow to IDs BEFORE update_d8_vars().
        #--------------------------------------------
        p_IDs     = self.d8.parent_ID_grid.flat[ IDs ]
        p_nbr_IDs = self.d8.parent_ID_grid.flat[ nbr_IDs ]
        mask1     = np.in1d( p_IDs, nbr_IDs )
        mask2     = np.in1d( p_nbr_IDs, IDs )
        flux_IDs1 = np.concatenate( (IDs[mask1], nbr_IDs[mask2]) )
        
        #--------------------------------------------------------
        # Update the D8 flow grid and all vars that depend on
        # it, including the D8 area grid, A, and slope grid, S.
        #--------------------------------------------------------
        # The set of cells where the area grid changes is a
        # superset of "IDs";  it includes other, downstream
        # cells.  These are saved by the d8.update_area_grid()
        # function as d8.new_A_IDs.
        #--------------------------------------------------------
        # The D8 flow grid gives the directions of fluxes and
        # Q=Q(A) and Qs=Qs(A,S) give the magnitudes.
        #--------------------------------------------------------
        all_IDs = np.concatenate( (IDs, nbr_IDs) )
        self.update_d8_vars( all_IDs, SILENT=SILENT, REPORT=REPORT )
        self.update_slope_grid( all_IDs, SILENT=SILENT, REPORT=REPORT )

        #--------------------------------------------
        # Find IDs that flow to nbr_IDs and nbr_IDs
        # that flow to IDs AFTER update_d8_vars().
        #--------------------------------------------
        ########################################################
        # Need to get new p_IDs, etc. after updating D8 vars.
        # Added next 2 lines to fix bug on 11/10/11.
        # But putting call to update_fluxes2() in update() did
        # fix some (or maybe all) problems.
        ########################################################
        p_IDs     = self.d8.parent_ID_grid.flat[ IDs ]
        p_nbr_IDs = self.d8.parent_ID_grid.flat[ nbr_IDs ]
        ########################################################
        mask3     = np.in1d( p_IDs, nbr_IDs )     
        mask4     = np.in1d( p_nbr_IDs, IDs )
        flux_IDs2 = np.concatenate( (IDs[mask3], nbr_IDs[mask4]) )

        #------------------------------------------
        # Get the set of all cells that have flow
        # between IDs and nbr_IDs.
        #------------------------------------------
        flux_IDs  = np.concatenate( (flux_IDs1, flux_IDs2) )
        flux_IDs  = np.unique( flux_IDs )
            
        #----------------------------------------------
        # Update fluxes Q and Qs for cells that have
        # flow between IDs and nbr_IDs.
        #----------------------------------------------
        self.update_Q_grid(  flux_IDs, SILENT=SILENT, REPORT=REPORT )
        self.update_Qs_grid( flux_IDs, SILENT=SILENT, REPORT=REPORT )

        #-----------------------------------------------------
        # Update the diffusion coefficient, D = (Qs / del_z)
        #-----------------------------------------------------
        # Added back in on (11/14/11), with new version.
        #-----------------------------------------------------
        self.update_D_grid( flux_IDs, SILENT=SILENT, REPORT=REPORT)

    #   update_fluxes2()
    #-------------------------------------------------------------------
    def update_fluxes3(self, IDs, nbr_IDs, SILENT=True, REPORT=False):

        #----------------------------------------------------------
        # Note: Based on comparison to Erode-D8-Global, it seems
        #       that update_fluxes2() may still not be right, even
        #       though there is no obvious bug now.  This is
        #       because networks in large depressions don't seem
        #       to develop as quickly (based on simulated time)
        #       as they do in Erode-D8-Global. (11/14/11)
        #
        #       This version is based on the idea that we always
        #       call update_fluxes() just before update_dz_dt(),
        #       and the latter only uses fluxes at IDs and at
        #       D8 neighbors that flow *towards* IDs. (11/14/11)
        
        #       This version only updates fluxes Q and Qs for
        #       cells that have flow BETWEEN IDs and nbr_IDs.
        #       So if a cell in IDs flows to another in IDs,
        #       their fluxes will NOT be updated and similarly
        #       for a cell in nbr_IDs that flows to another
        #       in nbr_IDs.  It is not yet clear whether this
        #       will be different than update_fluxes(). (9/16/11)
        #----------------------------------------------------------
        #       Note that "in1d()" function is new in NumPy 1.4.
        #       It tests whether each element of a 1D array is
        #       also present in a second array.
        #       ### Try setting assume_unique=True to speed up.
        #----------------------------------------------------------
  
        #--------------------------------------
        # Find nbr_IDs that flow *toward* IDs
        # BEFORE update_d8_vars().
        #--------------------------------------
        p_nbr_IDs = self.d8.parent_ID_grid.flat[ nbr_IDs ]
        mask1     = np.in1d( p_nbr_IDs, IDs )
        flux_IDs1 = nbr_IDs[ mask1 ]
        
        #--------------------------------------------------------
        # Update the D8 flow grid and all vars that depend on
        # it, including the D8 area grid, A, and slope grid, S.
        #--------------------------------------------------------
        # The set of cells where the area grid changes is a
        # superset of "IDs";  it includes other, downstream
        # cells.  These are saved by the d8.update_area_grid()
        # function as d8.new_A_IDs.
        #--------------------------------------------------------
        # The D8 flow grid gives the directions of fluxes and
        # Q=Q(A) and Qs=Qs(A,S) give the magnitudes.
        #--------------------------------------------------------
        all_IDs = np.concatenate( (IDs, nbr_IDs) )
        self.update_d8_vars( all_IDs, SILENT=SILENT, REPORT=REPORT )
        self.update_slope_grid( all_IDs, SILENT=SILENT, REPORT=REPORT )

        #--------------------------------------
        # Find nbr_IDs that flow *toward* IDs
        # AFTER update_d8_vars().
        #--------------------------------------
        p_nbr_IDs = self.d8.parent_ID_grid.flat[ nbr_IDs ]
        mask2     = np.in1d( p_nbr_IDs, IDs )
        flux_IDs2 = nbr_IDs[ mask2 ]

        #-------------------------------------
        # Get the set of all cells that have
        # (or had) flow from nbr_IDs to IDs.
        #-------------------------------------
        flux_IDs  = np.concatenate( (flux_IDs1, flux_IDs2) )
        flux_IDs  = np.concatenate( (flux_IDs, IDs) )  ## THIS TOO ???
        flux_IDs  = np.unique( flux_IDs )
            
        #----------------------------------------------
        # Update fluxes Q and Qs for cells that have
        # flow between IDs and nbr_IDs.
        #----------------------------------------------
        self.update_Q_grid(  flux_IDs, SILENT=SILENT, REPORT=REPORT )
        self.update_Qs_grid( flux_IDs, SILENT=SILENT, REPORT=REPORT )

        #-----------------------------------------------------
        # Update the diffusion coefficient, D = (Qs / del_z)
        #-----------------------------------------------------
        # Added back in on (11/14/11), with new version.
        #-----------------------------------------------------
        self.update_D_grid( flux_IDs, SILENT=SILENT, REPORT=REPORT)

    #   update_fluxes3()
    #-------------------------------------------------------------------
    def update_fluxes4(self, IDs, SILENT=True, REPORT=False):

        #---------------------------------------------------------
        # Update the D8 flow grid and all vars that depend on
        # it, including the D8 area grid, A, and slope grid, S.
        #---------------------------------------------------------
        # The D8 area grid is not just changed at IDs.
        # The set of cells where the area grid changes is a
        # superset of "IDs" that includes other, downstream
        # cells.  These are saved by the d8.update_area_grid()
        # function as d8.new_A_IDs.
        #
        # Note that update_d8_vars() is always called just
        # before update_Qs_grid() and update_dz_dt_grid(), so
        # updating the area grid at downstream IDs should never
        # be a problem.  However, saved grids of A, Q and Qs may
        # appear to be "out of sync" at these downstream IDs,
        # since they are snapshots.
        #---------------------------------------------------------
        # The D8 flow grid gives the directions of fluxes and
        # Q=Q(A) and Qs=Qs(A,S) give the magnitudes.
        #---------------------------------------------------------
        self.update_d8_vars( IDs, SILENT=SILENT, REPORT=REPORT )
        self.update_slope_grid( IDs, SILENT=SILENT, REPORT=REPORT )

        #----------------------------------------------
        # Update fluxes Q=Q(A) and Qs=Qs(A,S) for IDs
        #----------------------------------------------
        self.update_Q_grid(  IDs, SILENT=SILENT, REPORT=REPORT )
        self.update_Qs_grid( IDs, SILENT=SILENT, REPORT=REPORT )

        #-----------------------------------------------------
        # Update the diffusion coefficient, D = (Qs / del_z)
        #-----------------------------------------------------
        self.update_D_grid( IDs, SILENT=SILENT, REPORT=REPORT)        

    #   update_fluxes4()        
    #-------------------------------------------------------------------
    def update_flux_cap_grid(self, IDs, SILENT=True, REPORT=False):

        #------------------------------------------------------------
        # NB! This should only be called when synchronizing a
        #     neighbor cell, so cannot be merged into update_DEM().
        #------------------------------------------------------------
        # NB! CFL_factor is used to define dt_grid values, which
        #     are used to set T_next grid.  Don't use it here.
        #------------------------------------------------------------
        dt = (self.T_clock - self.T_last.flat[ IDs ])
        dz = (self.dz_dt.flat[ IDs ] * dt)

        #-----------------------------------------------
        # Impose a "minimum allowed dz" value
        # otherwise eventually dt goes to 0. (8/25/10)
        #-----------------------------------------------
        # (10/20/10) Tried again but not useful ?
        #-----------------------------------------------
        ## dz = np.maximum( self.min_allowed_dz, dz )
        ## dt = np.maximum( self.min_allowed_dt, dt )
        
        #----------------------------
        # Update the flux capacitor
        #----------------------------
        self.flux_cap_grid.flat[ IDs ] += dz
        
    #   update_flux_cap_grid()
    #-------------------------------------------------------------------
    def update_T_last_grid(self, IDs, SILENT=True, REPORT=False):

        self.T_last.flat[ IDs ] = self.T_clock
        
    #   update_T_last_grid() 
    #-------------------------------------------------------------------
    # def update_d8_vars(self):
    #
    #     (Inherited from erode_base.py)
    #     (But doesn't support IDs, so replace with one below.)
    
    #     update_d8_vars()
    #-------------------------------------------------------------------
    def update_d8_vars(self, IDs=None, SILENT=True, REPORT=False):

##        if not(SILENT):    
##            print 'Updating D8 variables...'
            
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
        #--------------------------------------------------
        # Note: d8.update() will use all IDs if IDs=None,
        # and this happens when initialize_d8_vars()
        # calls it.
        #--------------------------------------------------        
##        if (IDs is not None):
##            d8_before = self.d8.d8_grid.flat[ IDs ]
##            A_before  = self.d8.A.flat[ IDs ]


        #--------------------------------------------------
        # Perform updates only at IDs, but this results
        # in "backflow" situations where some of the pIDs
        # have a D8 parent (pID) in IDs. (2/29/12)
        #--------------------------------------------------
##        self.d8.IDs = IDs
##        self.d8.update( self.time, DEM=self.DEM,
##                        SILENT=SILENT, REPORT=REPORT )

        #-----------------------------------------------
        # This still doesn't catch all backflow cases.
        # It may be no different than next block.
        #-----------------------------------------------            
##        if (self.time_index == 0):
##            self.d8.IDs = IDs
##            self.d8.update( self.time, DEM=self.DEM,
##                            SILENT=SILENT, REPORT=REPORT )
##        else:
##            pIDs = self.d8.parent_ID_grid.flat[ IDs ]
##            affected_IDs = np.concatenate( (IDs, pIDs) )
##            affected_IDs = np.unique( affected_IDs )
##            self.d8.IDs = affected_IDs
##            self.d8.update( self.time, DEM=self.DEM,
##                            SILENT=SILENT, REPORT=REPORT )
            
        #---------------------------------------------------
        # Find the "grandparent IDs" (pIDs of pIDs) for
        # IDs and check if any gpIDs are also in IDs.
        # If so, the D8 flow codes for the corresponding
        # pIDs also need to be updated. (2/29/12)
        #---------------------------------------------------
        # This doesn't catch all "backflow" cases, though.
        #---------------------------------------------------
##        self.d8.IDs = IDs
##        self.d8.update( self.time, DEM=self.DEM,
##                        SILENT=SILENT, REPORT=REPORT )
##        if (self.time_index > 0):
##            pIDs  = self.d8.parent_ID_grid.flat[ IDs ]
##            gpIDs = self.d8.parent_ID_grid.flat[ pIDs ]
##            mask  = np.in1d( gpIDs, IDs )
##            bpIDs = pIDs[mask]  # ("bad" parent IDs)
##            if (bpIDs.size != 0):
##                ### print '###### FOUND bad parent IDs: updating...'
##                self.d8.IDs = bpIDs
##                self.d8.update( self.time, DEM=self.DEM,
##                                SILENT=SILENT, REPORT=REPORT )
        
        #---------------------------------------------------
        # Find the "grandparent IDs" (pIDs of pIDs) for
        # IDs and check if any gpIDs are also in IDs.
        # If so, the D8 flow codes for the corresponding
        # pIDs also need to be updated. (2/29/12)
        #---------------------------------------------------
        # This doesn't catch all "backflow" cases, though.
        #---------------------------------------------------
##        self.d8.IDs = IDs
##        self.d8.update( self.time, DEM=self.DEM,
##                        SILENT=SILENT, REPORT=REPORT )
##        if (self.time_index > 0):
##            pIDs  = self.d8.parent_ID_grid.flat[ IDs ]
##            gpIDs = self.d8.parent_ID_grid.flat[ pIDs ]
##            mask  = np.in1d( gpIDs, IDs )
##            #----------------------------------------------
##            # Next way to get "bad" parent IDs may be
##            # a bit different than the following 2 lines.
##            #----------------------------------------------
##            # bpIDs = pIDs[mask]  # ("bad" parent IDs)
##            bIDs  = gpIDs[ mask ]  # (these are ones in IDs)
##            bpIDs = self.d8.parent_ID_grid.flat[ bIDs ] 
##            if (bpIDs.size != 0):
##                ### print '###### FOUND bad parent IDs: updating...'
##                self.d8.IDs = bpIDs
##                self.d8.update( self.time, DEM=self.DEM,
##                                SILENT=SILENT, REPORT=REPORT )

        #-------------------------------------------
        # Update D8 info for all neighbors of IDs.
        #-------------------------------------------
##        if (self.time_index == 0):
##            self.d8.IDs = IDs
##        else:
##            #----------------------------------------
##            # Notice use of new INCLUDE_IDS option.
##            #----------------------------------------
##            nbr_hood_IDs = self.get_neighbors( IDs, INCLUDE_IDS=True )
##            self.d8.IDs  = nbr_hood_IDs
##            #--------------------------------------------------
##            # This shouldn't work because nbr_IDs does not
##            # include IDs, but it stops all error messages
##            # from update_area_grid(). Time goes faster, too?
##            #--------------------------------------------------
##            # nbr_IDs = self.get_neighbors( IDs )
##            # self.d8_IDs = nbr_IDs
##        self.d8.update( self.time, DEM=self.DEM,
##                        SILENT=SILENT, REPORT=REPORT )

        #----------------------------------------------------
        # (3/1/12) Moved this new "neighborhood" stuff into
        # d8_local.py where it belongs, along with the
        # get_neighbors() function.
        #----------------------------------------------------
        self.d8.IDs = IDs
        self.d8.update( self.time, DEM=self.DEM,
                        SILENT=SILENT, REPORT=REPORT )
        
 
##        if (IDs is not None):
##            d8_after = self.d8.d8_grid.flat[ IDs ]
##            A_after  = self.d8.A.flat[ IDs ]
##            if (d8_before[0] != d8_after[0]):
##                ## print 'D8_CODES_CHANGED       =', self.d8.D8_CODES_CHANGED
##                print 'd8 code: before, after =', d8_before, d8_after
##                print 'A value: before, after =', A_before, A_after
##                print '------------------------------------------------'
         
        #-----------------------------------
        # Check Amin and Amax, for testing
        #-----------------------------------
        if (self.DEBUG):
            print('#### Amin, Amax =', self.d8.A.min(), self.d8.A.max())

    #   update_d8_vars()
    #-------------------------------------------------------------------
    def update_slope_grid(self, IDs, SILENT=True, REPORT=False):

        #-----------------------------------------------------
        # Notes: Added IDs argument on 8/23/10.
        #
        #        *********  CHECK THIS  **********
        #        After DEM[active_ID] has been changed, the
        #        D8 flow codes of all 8 neighbor pixels may
        #        change, not just those of the D8 kids.  So
        #        their slopes may change as well.
        #-----------------------------------------------------
        # Notes: Make sure that d8 component is initialized
        #        with the same directory, site_prefix, etc.
        #        as this component.  Otherwise, the "shape"
        #        of DEM and ds, A, etc. won't match.
        #-----------------------------------------------------
        if not(SILENT):    
            print('Updating slope grid...')

        #---------------------------------------
        # Update slope (rise/run) at active_ID
        # and for all of its 8 neighbors
        #---------------------------------------
        pIDs   = self.d8.parent_ID_grid.flat[ IDs ]
        z_IDs  = self.DEM.flat[ IDs ]
        z_pIDs = self.DEM.flat[ pIDs ] 
        ds     = self.d8.ds.flat[ IDs ]
        self.S.flat[ IDs ] = (z_IDs - z_pIDs) / ds
        
        #--------------------------------------------
        # Set S=0 for IDs that have no D8 parent
        # (pits, flats, and edges), so Qs(Q,S) = 0.
        #--------------------------------------------
        w = np.where( pIDs == 0 )
        if (w[0].size > 0):
            self.S.flat[ IDs[w] ] = 0

        #------------------------------------------
        # Update the min and max values in S grid
        # See update_mins_and_maxes().
        #------------------------------------------

    #   update_slope_grid()
    #-------------------------------------------------------------------
    def update_Q_grid(self, IDs, SILENT=True, REPORT=False):

        #--------------------------------------------------------------
        # Notes:  If the flow code of any D8 kid has changed, then
        #         A[active_ID] and A[pID] will have been changed, but
        #         not the A-values at the D8 kid_IDs.  Q = Q(A,R)
        #--------------------------------------------------------------
        # Notes:  Q = annual discharge [m^3 / yr] (*leaving* a cell)
        #         R = geomorphically effective rainrate [meters / yr]
        #             (unless p ne 1, then R = coefficient)
        #         A = contributing area [meters^2]
        #--------------------------------------------------------------
        # Note: Q is zero wherever A is zero, and A will be zero
        #       wherever the D8 flow code is zero.
        #--------------------------------------------------------------        
        if not(SILENT):    
            print('Updating discharge grid...')
 
        #------------------------------------
        # Update discharge, Q, at active_ID
        # and for all of its 8 neighbors
        #------------------------------------
        A = self.d8.A.flat[ IDs ]

        #-------------------------------------
        # Get rainrate as scalar or 1D array
        #-------------------------------------
        if (self.R.size == 1):
            R = self.R
        else:
            R = self.R.flat[ IDs ]
        
        if (self.p != 1):
            Q = R * (A ** self.p)
        else:    
            Q = R * A

        self.Q.flat[ IDs ] = Q
      
        #--------------------------------------------------------
        # Test whether there are any pixels where flow code is
        # zero but area grid isn't.  This would be bad because
        # flow code is zero in pits, and A > 0 would lead to a
        # nonzero discharge (water and sed.) *out* of the cell.
        #--------------------------------------------------------
        # Tested for 100 timesteps on 2/22/10.
        #--------------------------------------------------------
##        w  = np.where(logical_and(self.d8.d8_grid == 0, self.d8.A != 0))
##        nw = w[0].size
##        if (nw != 0):
##            print 'WARNING: There are places where flow grid is'
##            print '         zero and area grid is nonzero.'
##            print ' '

        #------------------------------------------
        # Update the min and max values in S grid
        # See update_mins_and_maxes().
        #------------------------------------------
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            print('    discharge-area exponent =', self.p)
            print('    geomorphic rainrate     =', self.R)
            Q_str = str(np.nanmin(Q)) + ', ' + str(np.nanmax(Q))
            print('    min(Q), max(Q) = ' + Q_str + ' [m^3/yr]')
            A_str = str(np.nanmin(A)) + ', ' + str(np.nanmax(A))
            print('    min(A), max(A) = ' + A_str + ' [m^2]')

    #   update_Q_grid()
    #-------------------------------------------------------------------
    def update_Qs_grid(self, IDs, SILENT=True, REPORT=False):

        #---------------------------------------------------------
        # Notes: Since Qs = Qs(Q(A), S), we need to update Qs
        #        for any pixel that has had A or S updated.
        #        This includes active_ID (A and S updated), its
        #        D8 kids (S updated) and its parent (A updated).
        #---------------------------------------------------------
        # Notes: Qs = annual sed. discharge [m^3/ yr]
        #        Q  = annual discharge [m^3 / yr]
        #        S  = slope [unitless]
        #        K  = coefficient [(m^3 / yr)^(1 - m)]
        #        m  = discharge exponent (usually in [1,2])
        #        n  = slope exponent (usually in [1,2])

        #        The standard formula Qs = k * Q^m * S^n is for
        #        the case where Q and Qs have units of m^3/sec.
        #        The extra factor below makes the formula valid
        #        for the case where units of both are m^3/yr.
        #        That is, Qs' = fac * kf * (Q')^m * S^n.

        #        update_slope_grid() checks for negative
        #        slopes and either adjusts them or aborts.
        #---------------------------------------------------------
        if not(SILENT):
            print('Updating sed. discharge grid...')
        
        #------------------------------------
        # Update discharge, Q, at active_ID
        # and for all of its 8 neighbors
        #------------------------------------
        Q = self.Q.flat[ IDs ]
        S = self.S.flat[ IDs ]

        #-------------------------------------------------------
        # We could compute fac1 in initialize() and save it
        # but it probably wouldn't affect speed much. (9/2/11)
        #-------------------------------------------------------
        fac1 = self.secs_per_year ** (np.float64(1) - self.m)
        fac2 = self.K * (Q ** self.m)
        fac3 = (S ** self.n)
        Qs   = fac1 * fac2 * fac3     # [m^3 / year]
        #--------------------------
        self.Qs.flat[ IDs ] = Qs

        #------------------------------------------
        # Update the min and max values in S grid
        # See update_mins_and_maxes().
        #------------------------------------------
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            Qs_str = str(np.nanmin(Qs)) + ', ' + str(np.nanmax(Qs))
            print('    min(Qs), max(Qs) = ' + Qs_str + ' [m^3/yr]')

    #   update_Qs_grid()
    #-------------------------------------------------------------------
##    def update_D_grid_old(self, IDs, SILENT=True, REPORT=False):
##
##        #------------------------------
##        # THIS IS CURRENTLY NOT USED.
##        #------------------------------------------------------------
##        # Note: The LEM has the form:
##        #          div( D * grad(f) ) = f_t - U
##        #       where D is the diffusion coefficient, and given by:
##        #          D = (qs / S),   where qs = Qs / w.
##        #       Here, Qs(S,A) is computed as in last function, and
##        #       we can take w = dw (computed by D8 component).
##        #
##        #       For linear diffusion and constant D, the CFL
##        #       condition is:    dt < (dx^2) / (2 * D).
##        #
##        #       However, Omelchenko and Karimabadi (2006, p. 186)
##        #       discuss an example 1D problem where instead of 2*D
##        #       in the denominator, they use D_(i-1/2) + D_(i+1/2).
##        #
##        #       For our case, it may be better to use something
##        #       similar, e.g. (D[IDs] + D[pIDs]).
##        #------------------------------------------------------------
##        # NOTE: We have A=0 and S=0 where D8 code is zero, as in
##        #       pits and flats.  As a result, we have Q=0 and Qs=0,
##        #       which means that D will currently be undefined
##        #       anywhere the D8 code is zero.
##        #------------------------------------------------------------
##        if not(SILENT):    
##            print 'Updating D grid...'
##        
##        Qs = self.Qs.flat[ IDs ]
##        S  = self.S.flat[ IDs ]
##        w  = self.d8.dw.flat[ IDs ]
##        
##        self.D.flat[ IDs ] = (Qs / (w * S))
##
##        #-------------------------------------------
##        # Update the min and max D value (in grid)
##        #-------------------------------------------
##        D_min = np.nanmin( self.D.flat[ IDs ] )
##        D_max = np.nanmax( self.D.flat[ IDs ] )
##        if (np.isfinite(D_min)):
##            self.D_min = np.minimum( self.D_min, D_min )
##        if (np.isfinite(D_max)):
##            self.D_max = np.maximum( self.D_max, D_max )
##
##    #   update_D_grid_old()
    #-------------------------------------------------------------------
    def update_D_grid(self, IDs, SILENT=True, REPORT=False):

        #------------------------------------------------------------
        # Note: For linear diffusion and constant D, the CFL
        #       condition is:    dt < (dx^2) / (2 * D).
        #       Based on the new stability condition used in
        #       update_dt_grid5(), we have D = (Qs / del_z).
        #       Qs units are [m^3/yr], so D units are [m^2/yr].
        #------------------------------------------------------------
        #       This is not used for computations, but the spread
        #       of D values after a long time is of interest with
        #       regard to the benefits of the algorithm.
        #------------------------------------------------------------        
        if not(SILENT):    
            print('Updating D grid...')
 
        #------------------------------------
        # Compute downstream elevation drop
        #---------------------------------------------
        # del_z will always be positive as long as
        # LINK_FLATS = False. (flats have D8 code 0)
        #---------------------------------------------
        pIDs   = self.d8.parent_ID_grid.flat[ IDs ]
        z_IDs  = self.DEM.flat[ IDs ]
        z_pIDs = self.DEM.flat[ pIDs ]
        #--------------------------------
        del_z  = (z_IDs - z_pIDs)
        Qs = self.Qs.flat[ IDs ]
        #--------------------------------
        self.D.flat[ IDs ] = (Qs / del_z)  # [m^2/year]
        
        #------------------------------------------
        # Update the min and max values in D grid
        # See update_mins_and_maxes().
        #------------------------------------------
            
    #   update_D_grid()  
    #-------------------------------------------------------------------
    def update_dz_dt_grid(self, IDs, SILENT=True, REPORT=False):
        
        #-------------------------------------------------------------
        # Notes: This updates dz/dt at the specified IDs.
        #        dV_dt is the net rate at which the volume in a
        #        grid cell is changing.  It is computed as the
        #        difference between all of the fluxes into the grid
        #        cell and flux out of the grid cell. (D8-based here)
        #        dz_dt = (dV_dt / da)  [m/yr].
        #
        #        (dz_dt > 0) => deposition at that cell
        #        (dz_dt < 0) => erosion at that cell
        #
        #        Qs = K * Q^m * S^n = sed. discharge [m^3/yr]
        #        da = pixel area grid  [m^2]
        #-------------------------------------------------------------        
        if not(SILENT):    
            print('Updating dz/dt grid...')

        #-----------------        
        # Local synonyms
        #-----------------
        nx = self.nx
        ny = self.ny
        h  = self.d8.code_opps
        
        rows, cols = divmod( IDs, nx )
        
        #----------------------------------------------
        # Initialize dz_dt with the sediment flux out
        # of each pixel divided by the pixel area.
        # Should be 0 for pixels with flow code of 0.
        #----------------------------------------------
        # Note that in timestep, dt,  dz = dz_dt * dt
        #----------------------------------------------
        Qs = self.Qs.flat[ IDs ]
        if (self.da.size == 1):
            da = self.da
        else:
            da = self.da.flat[ IDs ]
        flux_out = Qs / da
        self.dz_dt.flat[ IDs ] = -flux_out   # [m / year]

        #-----------------------------------------
        # Add the contribution due to uplift, U.
        # self.U has units of [mm/yr], but
        # self.U_mpyr has units of [m/yr].
        #-----------------------------------------
        self.dz_dt.flat[ IDs ] += self.U_mpyr  # (must be scalar)
        
        #---------------------------------------------
        # Flow codes of 8 neighbor pixels (periodic)
        #--------------------------------------------------------------
        # In d8_local.start_new_codes(), d8_grid is set to zero on
        # left and right edges if not(LR_PERIODIC).  Similarly, it
        # is set to zero on top and bottom edges if not(TB_PERIODIC).
        # This ensures that flow from an opposite side does not
        # enter a cell unless periodic BCs are being used.
        #--------------------------------------------------------------
##        d0 = self.d8.d8_grid[(rows-1) % ny, (cols+1) % nx]   # (upper-right)
##        d1 = self.d8.d8_grid[rows,          (cols+1) % nx]   # (right)
##        d2 = self.d8.d8_grid[(rows+1) % ny, (cols+1) % nx]   # (lower-right)
##        d3 = self.d8.d8_grid[(rows+1) % ny,  cols]           # (bottom)
##        d4 = self.d8.d8_grid[(rows+1) % ny, (cols-1) % nx]   # (lower-left)
##        d5 = self.d8.d8_grid[rows,          (cols-1) % nx]   # (left)
##        d6 = self.d8.d8_grid[(rows-1) % ny, (cols-1) % nx]   # (upper-left)
##        d7 = self.d8.d8_grid[(rows-1) % ny,  cols]           # (top)

        #---------------------------------------------
        # Flow codes of 8 neighbor pixels (periodic)
        #------------------------------------------------
        # This is noticeably faster because things like
        # (cols+1)%nx are computed once and reused 3x.
        # Results from cProfile show this method's time
        # dropping from 14.6 to 12.1 seconds.
        #------------------------------------------------
        cols_R = (cols + 1) % nx
        cols_L = (cols - 1) % nx
        rows_U = (rows - 1) % ny
        rows_D = (rows + 1) % ny
        #--------------------------
        d0 = self.d8.d8_grid[rows_U,  cols_R]   # (upper-right)
        d1 = self.d8.d8_grid[rows,    cols_R]   # (right)
        d2 = self.d8.d8_grid[rows_D,  cols_R]   # (lower-right)
        d3 = self.d8.d8_grid[rows_D,  cols]     # (bottom)
        d4 = self.d8.d8_grid[rows_D,  cols_L]   # (lower-left)
        d5 = self.d8.d8_grid[rows,    cols_L]   # (left)
        d6 = self.d8.d8_grid[rows_U,  cols_L]   # (upper-left)
        d7 = self.d8.d8_grid[rows_U,  cols]     # (top)
        
        #-----------------------------------------
        # Add contributions from neighbor pixels
        #----------------------------------------------------
        # This part moves sediment downstream, and always
        # results in deposition, or a positive contribution
        # to dz/dt.  See Note above on periodic BCs.
        #----------------------------------------------------
        w0 = np.where( d0 == h[0] )    # (from upper-right)
        if (w0[0].size != 0):
            ## r0 = (rows[w0]-1) % ny
            ## c0 = (cols[w0]+1) % nx
            r0 = rows_U[ w0 ]
            c0 = cols_R[ w0 ]
            self.dz_dt.flat[ IDs[w0] ] += self.Qs[ r0,c0 ]/da
        #-------------------------------------------------------
        w1 = np.where( d1 == h[1] )    # (from right)
        if (w1[0].size != 0):
            ## r1 =  rows[w1]
            ## c1 = (cols[w1]+1) % nx
            r1 = rows[ w1 ]
            c1 = cols_R[ w1 ]
            self.dz_dt.flat[ IDs[w1] ] += self.Qs[ r1,c1 ]/da
        #-------------------------------------------------------
        w2 = np.where( d2 == h[2] )    # (from lower-right)
        if (w2[0].size != 0):
            ## r2 = (rows[w2]+1) % ny
            ## c2 = (cols[w2]+1) % nx
            r2 = rows_D[ w2 ]
            c2 = cols_R[ w2 ]
            self.dz_dt.flat[ IDs[w2] ] += self.Qs[ r2,c2 ]/da
        #-------------------------------------------------------
        w3 = np.where( d3 == h[3] )    # (from below)
        if (w3[0].size != 0):
            ## r3 = (rows[w3]+1) % ny
            ## c3 =  cols[w3]
            r3 = rows_D[ w3 ]
            c3 = cols[ w3 ]
            self.dz_dt.flat[ IDs[w3] ] += self.Qs[ r3,c3 ]/da
        #-------------------------------------------------------
        w4 = np.where( d4 == h[4] )    # (from lower-left)
        if (w4[0].size != 0):
            ## r4 = (rows[w4]+1) % ny
            ## c4 = (cols[w4]-1) % nx
            r4 = rows_D[ w4 ]
            c4 = cols_L[ w4 ]
            self.dz_dt.flat[ IDs[w4] ] += self.Qs[ r4,c4 ]/da
        #-------------------------------------------------------
        w5 = np.where( d5 == h[5] )    # (from left)
        if (w5[0].size != 0):
            ## r5 =  rows[w5]
            ## c5 = (cols[w5]-1) % nx
            r5 = rows[ w5 ]
            c5 = cols_L[ w5 ]
            self.dz_dt.flat[ IDs[w5] ] += self.Qs[ r5,c5 ]/da
        #-------------------------------------------------------
        w6 = np.where( d6 == h[6] )    # (from upper left)
        if (w6[0].size != 0):
            ## r6 = (rows[w6]-1) % ny
            ## c6 = (cols[w6]-1) % nx
            r6 = rows_U[ w6 ]
            c6 = cols_L[ w6 ]
            self.dz_dt.flat[ IDs[w6] ] += self.Qs[ r6,c6 ]/da
        #-------------------------------------------------------
        w7 = np.where( d7 == h[7] )    # (from above)
        if (w7[0].size != 0):
            ## r7 = (rows[w7]-1) % ny
            ## c7 =  cols[w7]
            r7 = rows_U[ w7 ]
            c7 = cols[ w7 ]
            self.dz_dt.flat[ IDs[w7] ] += self.Qs[ r7,c7 ]/da

        #--------------------------------------------------------
        # This is only done for "completeness" (as when viewing
        # values in a grid stack).  The values aren't actually
        # used, since elevations at base_IDs are set by the
        # update_base_level() function.
        #--------------------------------------------------------        
        # Use Base-level Lowering Rate (BLR) to set
        # dz_dt on edges.
        # self.BLR_mpyr has units of [m/yr].
        #----------------------------------------------
        # base_IDs are defined in the function:
        # erode_base.initialize_boundary_conditions()
        #----------------------------------------------        
        self.dz_dt.flat[ self.base_IDs ] = -self.BLR_mpyr
        
        ### self.dz_dt.flat[ self.base_IDs ] -= self.BLR_mpyr  ############ (9/30/10)

        #------------------------------------------------
        # If U>0 and BLR=0, then we need this (9/30/10)
        #------------------------------------------------
        ## self.dz_dt.flat[ self.base_IDs ] = 0
        
        #----------------------
        # This may be faster.
        #----------------------
##        edge_IDs = np.intersect1d_nu( IDs, self.base_IDs )
##        self.dz_dt.flat[ edge_IDs ] = -self.BLR_mpyr

        #----------------------------------------------
        # Update the min and max values in dz_dt grid
        # See update_mins_and_maxes().
        #----------------------------------------------
            
        #--------------------------------------
        # For testing.  Check where dz_dt = 0
        #--------------------------------------
##        w  = np.where(self.dz_dt == 0)
##        nw = w[0].size
##        if (nw != 0):
##            print 'dz_dt grid is zero at', nw, 'pixel(s).'
##            if (nw <= 20):
##                for k in xrange(nw):
##                    print '(row, col) =', w[0][k], w[1][k]
##            self.dz_dt[w] = self.BLR_mpyr  # (make positive)
        
    #   update_dz_dt_grid()
    #-------------------------------------------------------------------
    def update_mins_and_maxes(self, IDs, REPORT=False):

        #-----------------------------------------------------
        # (1/25/12) Introduced this when cProfile showed the
        # cost of nanmin(), nanmax(), etc.
        #
        # Note. ".min()" and ".max()" may be faster than
        # np.nanmin() or even np.min().  Check this.
        #-----------------------------------------------------

        ######## THIS IS STILL DONE IN update_DEM() #######
        
        #-----------------------------------------
        # Update the min and max value in z grid
        #-----------------------------------------
        zmin = self.DEM.flat[ IDs ].min()
        zmax = self.DEM.flat[ IDs ].max()
##        zmin = np.nanmin( self.DEM.flat[ IDs ] )
##        zmax = np.nanmax( self.DEM.flat[ IDs ] )
        #-------------------------------------------
        # Check whether zmax ever increases, which
        # should never occur (10/14/11)
        #-------------------------------------------
        if (zmax > self.DEM_max):
            print('######################################')
            print(' ERROR: zmax just increased at')
            print(' time_index =', self.time_index)
            print(' from', self.DEM_max, 'to', zmax, '.')
            if (IDs.size == 1):
                print(' at ID =', IDs[0])
            print('######################################')            
        self.DEM_min = np.minimum( self.DEM_min, zmin )
        self.DEM_max = np.maximum( self.DEM_max, zmax )
        
        #-----------------------------------------
        # Update the min and max value in S grid
        #-----------------------------------------
        S_min = np.nanmin( self.S.flat[ IDs ] )
        S_max = np.nanmax( self.S.flat[ IDs ] )
        self.S_min = np.minimum( self.S_min, S_min )
        self.S_max = np.maximum( self.S_max, S_max )

        #-----------------------------------------
        # Update the min and max value in A grid
        #-----------------------------------------

        #-----------------------------------------
        # Update the min and max value in Q grid
        #-----------------------------------------

        #------------------------------------------
        # Update the min and max value in Qs grid
        #------------------------------------------

        #---------------------------------------------
        # Update the min and max value in dz_dt grid
        #---------------------------------------------
        
        #-----------------------------------------
        # Update the min and max value in D grid
        #-----------------------------------------
        D_min = np.nanmin( self.D.flat[ IDs ] )
        D_max = np.nanmax( self.D.flat[ IDs ] )
        if (np.isfinite(D_min)):
            self.D_min = np.minimum( self.D_min, D_min )
        if (np.isfinite(D_max)):
            self.D_max = np.maximum( self.D_max, D_max )
        #----------------------------------------------------
        w = np.where( self.D > 0 )
        self.D_min_pos = np.nanmin( self.D[w] )
                                   
        # Same for self.dz_dt, self.Q, self.Qs, self.S, self.A

        #------------------
        # Optional report
        #------------------
        if (REPORT):
            S_str = str(self.S_min) + ', ' + str(self.S_max)
            print('    min(S), max(S) = ' + S_str + ' [m/m]')

            dz_dt_min = self.dz_dt.min()
            dz_dt_max = self.dz_dt.max()
            print('    min(dz/dt), max(dz/dt) =', dz_dt_min, dz_dt_max)
            print(' ')
            
        
    #   update_mins_and_maxes()
    #-------------------------------------------------------------------
    def update_grid_mins_and_maxes(self, REPORT=False):

        #------------------------------------------------------
        # Note: cProfile shows that nanmin() and nanmax() are
        # fairly costly.  However, this method is called only
        # once at the end of a model run to compute all mins
        # and maxes.  See update_mins_and_maxes( IDs ).
        #
        # Note. ".min()" and ".max()" may be faster than
        # np.nanmin() or even np.min().  Check this.
        #------------------------------------------------------
        
        #-----------------------------------------
        # Update the min and max value in z grid
        #-----------------------------------------
        self.DEM_min = self.DEM.min()
        self.DEM_max = self.DEM.max()
##        self.DEM_min = np.nanmin( self.DEM )
##        self.DEM_max = np.nanmax( self.DEM )
        
        #-----------------------------------------
        # Update the min and max value in S grid
        #-----------------------------------------
        self.S_min = self.S.min()
        self.S_max = self.S.max()
##        self.S_min = np.nanmin( self.S )
##        self.S_max = np.nanmax( self.S)
        #----------------------------------------
        # We've set S=0 on edges, but want to
        # know the smallest positive value.
        #----------------------------------------        
        w = np.where( self.S > 0 )
        self.S_min_pos = np.nanmin( self.S[w] )
        
        #-----------------------------------------
        # Update the min and max value in A grid
        #-----------------------------------------
        self.A_min = self.d8.A.min()
        self.A_max = self.d8.A.max()
##        self.A_min = np.nanmin( self.d8.A )
##        self.A_max = np.nanmax( self.d8.A )
        
        #-----------------------------------------
        # Update the min and max value in Q grid
        #-----------------------------------------
        self.Q_min = self.Q.min()
        self.Q_max = self.Q.max()
##        self.Q_min = np.nanmin( self.Q )
##        self.Q_max = np.nanmax( self.Q )
        #----------------------------------------
        # We've set Q=0 on edges, but want to
        # know the smallest positive value.
        #----------------------------------------        
        w = np.where( self.Q > 0 )
        self.Q_min_pos = np.nanmin( self.Q[w] )
        
        #------------------------------------------
        # Update the min and max value in Qs grid
        #------------------------------------------
        self.Qs_min = self.Qs.min()
        self.Qs_max = self.Qs.max()
##        self.Qs_min = np.nanmin( self.Qs )
##        self.Qs_max = np.nanmax( self.Qs )
        #----------------------------------------
        # We've set Qs=0 on edges, but want to
        # know the smallest positive value.
        #----------------------------------------        
        w = np.where( self.Qs > 0 )
        self.Qs_min_pos = np.nanmin( self.Qs[w] )
        
        #---------------------------------------------
        # Update the min and max value in dz_dt grid
        #---------------------------------------------
        self.dz_dt_min = self.dz_dt.min()
        self.dz_dt_max = self.dz_dt.max()
##        self.dz_dt_min = np.nanmin( self.dz_dt )
##        self.dz_dt_max = np.nanmax( self.dz_dt )
        
        #-----------------------------------------
        # Update the min and max value in D grid
        #-----------------------------------------
        # Note:  At some point, update_D_grid()
        # will compute D[0,0] and get NaN since
        # pID[0,0] = 0 means del_z = 0.
        #-----------------------------------------
        self.D[0,0] = 0   #####
        self.D_min  = self.D.min()
        self.D_max  = self.D.max()
##        self.D_min = np.nanmin( self.D )
##        self.D_max = np.nanmax( self.D )
        #----------------------------------------
        # We've set D=0 on edges, but want to
        # know the smallest positive value.
        #----------------------------------------        
        w = np.where( self.D > 0 )
        self.D_min_pos = np.nanmin( self.D[w] )

        #------------------------------------------
        # Update the min and max value in dt grid
        #------------------------------------------
        self.dt_min  = self.dt_grid.min()
        ## self.dt_max  = self.dt_grid.max()
##        self.dt_min = np.nanmin( self.dt_grid )
##        self.dt_max = np.nanmax( self.dt_grid )
        #-------------------------------------------
        # We've set (dt = dt_limit) for cells with
        # (d8_code == 0), but want to know largest
        # legitimate value.
        #-------------------------------------------        
        w = np.where( self.dt_grid < self.dt_limit )
        self.dt_max = self.dt_grid[w].max()
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            S_str = str(self.S_min) + ', ' + str(self.S_max)
            print('    min(S), max(S) = ' + S_str + ' [m/m]')

            dz_dt_min = self.dz_dt.min()
            dz_dt_max = self.dz_dt.max()
            print('    min(dz/dt), max(dz/dt) =', dz_dt_min, dz_dt_max)
            print(' ')
            
        
    #   update_grid_mins_and_maxes() 
    #-------------------------------------------------------------------
##    def update_dt_grid4(self, IDs, SILENT=True, REPORT=False):
##
##        #---------------------------------------------------------------
##        # Notes: This version is based on the following rule:
##        #
##        #        "No pixel can become lower than its D8 parent in
##        #        one timestep.  A one-pixel pit cannot be raised
##        #        higher than its lowest neighbor in one timestep."
##        #
##        #        Applying this rule everywhere has the consequence
##        #        that no pixel can become higher than any of its
##        #        D8 child pixels in one timestep.
##        #
##        #        z2[ID]  = z[ID]  + z_dot[ID]  * dt
##        #        z2[pID] = z[pID] + z_dot[pID] * dt
##        #
##        #        where pID is the D8 parent pixel for pixel at ID.
##        #        Note that same dt is used for both. If we now
##        #        subtract these 2 equations, we find that:
##        #
##        #        If (z_dot[pID] > z_dot[ID]), then we need:
##        #            dt < (z[ID] - z[pID]) / (z_dot[pID] - z_dot[ID])
##        #        In this case, z[ID] and z[pID] are converging,
##        #        so the slope, S[ID], is getting smaller.
##        
##        #        If (z_dot[pID] = z_dot[ID]), then we need:
##        #            dt < Infinity
##        #        In this case, z[ID] and z[pID] are maintaining
##        #        the same difference (and S[ID]) over time.
##        #
##        #        If (z_dot[pID] < z_dot[ID]), then we need:
##        #            dt > -(z[ID] - z[pID]) / (z_dot[ID] - z_dot[pID])
##        #            dt > (z[ID] - z[pID]) / (z_dot[pID] - z_dot[ID])
##        #            dt > 0  (since both diffs above > 0)
##        #        In this case, z[ID] and z[pID] are diverging,
##        #        so the slope, S[ID], is getting larger.
##        #
##        #--------------------------------------------------------------- 
##        # NB!    If we use CFL_factor = 1, then we'll have
##        #        z2[ID] = z2[pID] instead of z2[ID] > z2[pID].
##        #---------------------------------------------------------------   
##        if not(SILENT):
##            print 'Updating dt_grid...'
##
##        pIDs = self.d8.parent_ID_grid.flat[ IDs ]
##        
##        #------------------------------------
##        # Compute downstream elevation drop
##        #------------------------------------------------
##        # del_z will always be positive as long as
##        # LINK_FLATS = False. (flats have D8 code 0)
##        #------------------------------------------------
##        # NB! del_z will be wrong wherever pID=0 since
##        # it will be computed from DEM[0].  Places with
##        # pID=0 are fixed further down.
##        #------------------------------------------------
##        z_IDs  = self.DEM.flat[ IDs ]
##        z_pIDs = self.DEM.flat[ pIDs ]
##        del_z  = (z_IDs - z_pIDs)
##        
##        #-------------------------------------------------
##        # See Notes above. This can be pos, neg or zero.
##        #---------------------------------------------------
##        # NB! del_dz_dt will be wrong wherever pID=0 since
##        # it will be computed from dz_dt[0].  Places with
##        # pID=0 are fixed further down.
##        #---------------------------------------------------
##        dz_dt_IDs  = self.dz_dt.flat[ IDs  ]
##        dz_dt_pIDs = self.dz_dt.flat[ pIDs ]
##        del_dz_dt  = (dz_dt_IDs - dz_dt_pIDs)
##
##        #-----------------------------------------
##        # Get last update times for IDs and pIDs
##        #-----------------------------------------
##        T_last_IDs  = self.T_last.flat[ IDs ]  # (will equal T_clock now)
##        ## T_last_IDs  = self.T_clock  # (would need to have size(IDs))
##        T_last_pIDs = self.T_last.flat[ pIDs ]
##        term1       = T_last_pIDs * dz_dt_pIDs
##        term2       = T_last_IDs  * dz_dt_IDs
##        TL_term     = (term1 - term2)
##            
##        #------------------------------
##        # Compute stable dt, in years
##        #-----------------------------------------------------
##        # These statements are true regardless of TL term,
##        # because del_dz_dt = d/dt(z[ID] - z[pID]).
##        #-----------------------------------------------------        
##        # If (del_dz_dt < 0), z[ID] & z[pID] are converging.
##        # If (del_dz_dt > 0), z[ID] & z[pID] are diverging.
##        # If (del_dz_dt = 0), z[ID] & z[pID] are neither.
##        #-----------------------------------------------------
##        w1 = np.where( del_dz_dt != 0 )
##        if (w1[0].size > 0):
##            abs_del_dz_dt = np.absolute( del_dz_dt[w1] )
##            T_next_equal = (del_z[w1] + TL_term[w1]) / abs_del_dz_dt
##            dt_equal     = (T_next_equal - self.T_clock)
##            self.dt_grid.flat[ IDs[w1] ] = self.CFL_factor * dt_equal
##
####            abs_del_dz_dt = np.absolute( del_dz_dt[w1] )
####            T_next_equal = (del_z[w1] / abs_del_dz_dt)
####            dt_equal     = (T_next_equal - self.T_clock)
####            dt_equal     = np.absolute(dt_equal)
####            self.dt_grid.flat[ IDs[w1] ] = self.CFL_factor * dt_equal
##            
####            T_next_equal = (del_z[w1] + TL_term[w1]) / del_dz_dt[w1]
####            dt_equal     = (T_next_equal - self.T_clock)
####            ## dt_equal     = np.absolute( dt_equal )   #########
####            self.dt_grid.flat[ IDs[w1] ] = self.CFL_factor * dt_equal
##            
##            #---------------------------------------------------
##            # (10/22/10) If z[ID] & z[pID] are diverging, then
##            # there will never be a future time when they are
##            # equal, but there could be a time in the past.
##            # In this case it seems dt_equal can be negative.
##            #---------------------------------------------------
##            # (10/20/10) This happens at some pits and their
##            # values are set separately below.
##            #---------------------------------------------------
##            w0 = np.where(dt_equal <= 0)
##            if (w0[0].size > 0):
##                bad_IDs = IDs[w1][w0]
##                ## self.dt_grid.flat[ bad_IDs ] = 1e-4       ##############################
##                ## self.dt_grid.flat[ bad_IDs ] = 10000.0    ##############################
##                self.dt_grid.flat[ bad_IDs ] = 5000.0
##                ## self.dt_grid.flat[ bad_IDs ] = np.absolute( self.dt_grid.flat[ bad_IDs ] )
##                
##                #--------------------------------------------------------
##                # (10/25/10) This should be okay since there is not a
##                # constraint on dt to ensure that z[ID] > z[pID].  If
##                # the pixel is reactivated in the future, then, using
##                # dt = (T_clock - T_last) should not cause instability.
##                #--------------------------------------------------------
##                ## self.dt_grid.flat[ bad_IDs ] = self.dt_limit   # (before 10/22/10)              
####                print '############################################'
####                print ' In update_dt_grid4():'
####                print '    dt_equal is NEGATIVE at IDs:'
####                print bad_IDs
####                print 'T_next_equal =', T_next_equal
####                print 'T_clock      =', self.T_clock
####                print 'dt_equal     =', dt_equal
####                print '############################################'
####                print ' '
####                sys.exit()              
##
##        #----------------------------------------------------
##        # (10/21/10) Pits might match some of this, so make
##        # sure it comes before pit handling section.
##        #----------------------------------------------------
##        w2 = np.where( del_dz_dt == 0 )
##        if (w2[0].size > 0):
####            print '################# del_dz_dt = 0 at IDs ='
####            print IDs[w2]
####            print 'and time index =', self.time_index
####            print ' '
##            self.dt_grid.flat[ IDs[w2] ] = self.dt_limit
##    
##        w3 = np.where( del_z <= 0)
##        if (w3[0].size > 0):
####            print '################# del_z <= 0 at IDs ='
####            print IDs[w3]
####            print 'and time index =', self.time_index
####            print ' '
##            self.dt_grid.flat[ IDs[w3] ] = self.dt_limit
##                            
##        #----------------------------------------------------
##        # Deactivate all pixels that have a flow code of 0.
##        # Their z values can then only be changed when they
##        # are synchronized by one of their neighbors.
##        #----------------------------------------------------
##        w4 = np.where( pIDs == 0 )
##        if (w4[0].size > 0):
##            self.dt_grid.flat[ IDs[w4] ] = self.dt_limit
##        ##########################################################
##        #  Above block is for the edge pixels (and flats)
##        ##########################################################
##        
##        #------------------------------------------------------------
##        # For each pit, set its dt-value so that it will fill in
##        # one timestep. It will fill to the level of the pixel that
##        # is its lowest neighbor at that future time. (10/18/10)  
##        #------------------------------------------------------------
##        # NB! Pits are lower than *all* of their neighbors, but
##        #     some of them may not flow into the pit.  We must have
##        #     (del_z > 0)  for the pit's D8 kids, but we could have
##        #     (del_z <= 0) for other neighbors.
##        #------------------------------------------------------------
##        # Two issues: Pit may not become *higher* than any of its
##        # D8 kids (in one step), BUT pit should only be filled to
##        # the point where it is as high as its lowest neighbor.
##        # Are both issues covered as written ?? (10/19/10)
##        #------------------------------------------------------------
##        ## interior_pit_IDs = self.d8.pit_IDs
##        ## interior_pit_IDs = self.get_interior_pit_IDs( FLATS=True )
##        interior_pit_IDs = self.get_interior_pit_IDs( FLATS=False )
##        
##        for k in xrange(len(IDs)):
##            ID  = IDs[k]
##            pID = pIDs[k]
##            if (ID in interior_pit_IDs):
##            ## if (ID in self.d8.pit_IDs):
##                #---------------------------------------------------
##                # (10/20/10) Results seem very similar now whether
##                #  we use D8 kids or all 8 neighbors.  COMPARE.
##                #---------------------------------------------------
##                ## kid_IDs = self.d8.get_kid_IDs( ID )
##                kid_IDs = self.d8.get_neighbors( ID )
##
##                #---------------------------------------------------
##                # (10/20/10) Remove neighbors with flow code of 0.
##                # This has no effect because there can't be any.
##                #---------------------------------------------------                
####                ww = np.where( self.d8.d8_grid.flat[ kid_IDs ] != 0 )
####                if (ww[0].size > 0):
####                    kid_IDs = kid_IDs[ ww ]
##                    
##                if (kid_IDs.size > 0):
##                    z_kids         = self.DEM.flat[ kid_IDs ]
##                    del_z_kids     = (z_kids - z_IDs[k])
##                    dz_dt_kids     = self.dz_dt.flat[ kid_IDs ]
##                    del_dz_dt_kids = (dz_dt_kids - dz_dt_IDs[k])
##                    # Use this if we want to include flats in with pits. (10/20/10)
####                    wcon = np.where( np.logical_and( del_dz_dt_kids < 0,
####                                               del_z_kids > 0 ) )
##
####                    wcon = np.where( np.logical_and( del_dz_dt_kids < 0,
####                                               del_z_kids > 1e-3 ) ) # (converging)
##                    ##                         del_z_kids > self.min_allowed_dz ) ) # (converging)
##                    wcon = np.where( del_dz_dt_kids < 0 ) # (converging)
##                    if (wcon[0].size > 0):
##                        kid_IDs        = kid_IDs[ wcon ]
##                        del_z_kids     = del_z_kids[ wcon ]
##                        dz_dt_kids     = dz_dt_kids[ wcon ]
##                        del_dz_dt_kids = del_dz_dt_kids[ wcon ]
##                        #-----------------------------------------------------
##                        # NB!  In update_T_last_grid(), we set T_last_IDs[k]
##                        #      equal to T_clock, and that is called *before*
##                        #      update_dt_grid4().
##                        #-----------------------------------------------------
##                        T_last_kids  = self.T_last.flat[ kid_IDs ]
##                        term1        = T_last_IDs[k]  * dz_dt_IDs[k]  # (the pit)
##                        term2        = T_last_kids    * dz_dt_kids
##                        TL_term_kids = (term1 - term2)
##                        #---------------------------------------------------        
##                        abs_del_dz_dt = np.absolute(del_dz_dt_kids) 
##                        TN_pit_vals = (del_z_kids + TL_term_kids) / abs_del_dz_dt
##                        dt_equal= (TN_pit_vals - self.T_clock)
##                        #-------------------------------------------------
##                        # (10/21/10) Somehow, dt_equal can be negative,
##                        # so make sure we get a positive value.
##                        #-------------------------------------------------
####                        if (dt_equal.min() < 0):
####                            print '############################################'
####                            print ' ERROR: dt-value is NEGATIVE.  Aborting.'
####                            print '############################################'
####                            print 'del_z_kids =', del_z_kids
####                            print 'dt_equal   =', dt_equal
####                            print ' '
####                            ## sys.exit()
##                        ## wp = np.where( dt_equal > 1e-1 )  ########  DOESN'T WORK (10/22/10)
##                        wp = np.where( dt_equal > 0 )
##
##                        if (wp[0].size > 0):
##                            dt_equal = dt_equal[ wp ]
##                        #---------------------------------------------------
##                        self.dt_grid.flat[ ID ] = dt_equal.min()
##                        if (ID in self.d8.flat_IDs):
##                            self.dt_grid.flat[ ID ] *= self.CFL_factor
##                        ## self.dt_grid.flat[ ID ] = dt_equal.min() * self.CFL_factor
##                        ## self.dt_grid.flat[ ID ] = dt_equal.min() * 1.1
##                        #---------------------------------------------
##                        # (10/22/10) EXPERIMENT.  If one of pit's kids
##                        # is active_ID, then change kid's dt-value.
##                        #---------------------------------------------
####                        if (self.active_ID in kid_IDs):
####                            self.dt_grid.flat[ self.active_ID ] /= self.CFL_factor
##                        #----------------------------------------------------
##                        # (10/22/10) EXPERIMENT. Make sure kid's are
##                        # updated AFTER their pit parent.
##                        #----------------------------------------------------
####                        epsilon = 1e-08
####                        self.dt_grid.flat[ kid_IDs ] = np.maximum( dt_equal.min() + epsilon,
####                                                                      self.dt_grid.flat[ kid_IDs ] )
##                        #----------------------------------------------------
##                        # (10/22/10) EXPERIMENT. Deactivate pits neighbors.
##                        # Seems to be a bad idea.
##                        #----------------------------------------------------
##                        ## self.dt_grid.flat[ kid_IDs ] = self.dt_limit
##                        
##                        #----------------------------------------------------
##                        # (10/20/10) Need a way to deal with "2-pixel pits"
##                        # that goes right about here.
##                        #----------------------------------------------------
####                        wm = np.where( dt_equal == dt_equal.min() )
####                        if (del_z_kids[ wm ] > 1e-6):
####                            self.dt_grid.flat[ ID ] = dt_equal[ wm ]
####                        else:
####                            self.dt_grid.flat[ ID ] = 1.0
##                            ## self.dt_grid.flat[ ID ] = self.dt_limit  # (DEACTIVATE, doesn't work)
##
##        #-----------------------------------------------
##        # (10/22/10) EXPERIMENT to deal with stalling. 
##        #-----------------------------------------------
##        wp = np.where( IDs == self.active_ID )
##        if (wp[0].size > 0):
##            if (self.active_ID_count > 50):
####                ## self.update_entire_DEM()   # (but doesn't update flow grid, etc.)
####                ## self.dt_grid.flat[ self.active_ID ] /= self.CFL_factor
####                ## self.dt_grid.flat[ self.active_ID ] += 50.0
####                print ' '
##                print '###### Progress has stalled at ID =', self.active_ID
####                print '###### Aborting.  Time index =', self.time_index
####                print ' '
####                sys.exit()
##                ## self.dt_grid.flat[ self.active_ID ] /= self.CFL_factor
##                ## self.dt_grid.flat[ self.active_ID ] = 1e-3
##                ## self.dt_grid.flat[ self.active_ID ] = np.maximum( self.dt_grid.flat[ self.active_ID ], 1e-3 )
##                #-----------------------------------------------------------
##                # (10/25/10) The problem with this is that when/if the
##                # pixel gets reactivated (i.e. when synced by a neighbor),
##                # then the new dt = (T_clock - T_last) will be large and
##                # a pit is likely to form.
##                #-----------------------------------------------------------
##                # self.dt_grid.flat[ self.active_ID ] = self.dt_limit
##
##        #-------------------------------------------------------------
##        # (10/22/10) Where dt is very small, divide it by CFL_factor 
##        #-------------------------------------------------------------
####        wp = np.where( IDs == self.active_ID )
####        if (wp[0].size > 0):
####            pID = pIDs[ wp ]  # (parent of active_ID)
####            dt  = self.dt_grid.flat[ self.active_ID ]
####            if (pID in interior_pit_IDs) and (dt < 1e-1):
####                print '########### STALLING AT active_ID =', self.active_ID               
####                self.dt_grid.flat[ self.active_ID ] = self.dt_grid.flat[ pID ] + 1
####                self.dt_grid.flat[ self.active_ID ] = self.dt_grid.flat[ pID ] + 1e-07
####                self.dt_grid.flat[ self.active_ID ] /= self.CFL_factor  ########
####                self.dt_grid.flat[ self.active_ID ] = self.dt_limit     ########
##
####        w5 = np.where( self.dt_grid.flat[ IDs ] < 1e-1 )
####        if (w5[0].size > 0):
######            self.dt_grid.flat[ IDs[w5] ] = np.maximum( self.dt_grid.flat[ IDs[w5] ],
######                                                          self.dt_grid.flat[ pIDs[w5] ] )          
####            self.dt_grid.flat[ IDs[w5] ] = self.dt_grid.flat[ IDs[w5] ] / self.CFL_factor
##            
##        #------------------------------------------
##        # Deactivate pixels where dt is too small
##        #------------------------------------------
####        w5 = np.where( self.dt_grid.flat[ IDs ] < 1e-1 )  ##### EXPERIMENT
####        w5 = np.where( self.dt_grid.flat[ IDs ] < 1e-2 )
####        w5 = np.where( self.dt_grid.flat[ IDs ] < 1e-3 )
####        if (w5[0].size > 0):
####            self.dt_grid.flat[ IDs[w5] ] = self.dt_limit
##                    
##        #--------------------------
##        # Save the min and max dt
##        #------------------------------------------
##        # NB! Don't do min/max over whole grid !!
##        #------------------------------------------
##        dt_vals = self.dt_grid.flat[ IDs ]
##        dt_min  = dt_vals.min()
##        dt_max  = dt_vals.max()
##
####        if (dt_min < self.dt_too_small):
##            
##        self.dt_min = np.minimum( self.dt_min, dt_min )
##        self.dt_max = np.maximum( self.dt_max, dt_max )
##        if not(SILENT):
####            dt_vals = self.dt_grid.flat[ IDs ]
####            dt_min  = dt_vals.min()
####            dt_max  = dt_vals.max()
####            self.dt_min = np.minimum( self.dt_min, dt_min )
####            self.dt_max = np.maximum( self.dt_max, dt_max )
##            print '#########################################'
##            print ' dt_min =', np.around(self.dt_min, 2)
##            print ' in update_dt_grid() (del_z method)'
##            print '#########################################'
##        
##        #-------------------------------
##        # Don't let dt get too small ?
##        #-----------------------------------------------------
##        # dt_min must be less than 1e-4 for the case m=1.5,
##        # n=1.0, K=1.0. Making K smaller allows dt_too_small
##        # to be bigger. Now K is smaller.
##        #-----------------------------------------------------
##        # Fixed units issue in Qs_Grid so should be able to
##        # reduce dt_too_small.
##        #-----------------------------------------------------
##        if (dt_min < self.dt_too_small):
##            print '### WARNING: dt < ' + str(self.dt_too_small)
##            print '        dt = ' + str(dt_min)
##            print ' '
##            if (dt_min <= 0):
##                sys.exit()
####            print '******************************************'
####            print ' Aborting: Stable dt is too small.'
####            print ' Computed dt = ' + str(dt_grid_min)
####            print '******************************************'
####            print ' '
####            sys.exit()
##
##        #------------------
##        # Optional report
##        #------------------
##        ##### REPORT = True  ###################
##        if (REPORT):
##            dt_str     = str(self.dt_min)  + ', ' + str(self.dt_max)
##            del_z_str  = str(del_z.min())  + ', ' + str(del_z.max())
##            del_dz_dt_str = str(del_dz_dt.min()) + ', ' + str(del_dz_dt.max())
##            print '    min(dt),     max(dt)           = ' + dt_str + ' [yrs]'
##            print '    min(del_z),  max(del_z)        = ' + del_z_str
##            print '    min(del_dz_dt), max(del_dz_dt) = ' + del_dz_dt_str
##            print ' '
##            
##    #   update_dt_grid4()
    #-------------------------------------------------------------------
    def update_dt_grid5(self, IDs, SILENT=True, REPORT=False):

        #---------------------------------------------------------------
        # Notes: This version is based on the following stability
        #        condition:
        #
        #        dt < dx * dy * [z(k) - z(p)] / [2 * Qs(k)]
        #
        #        This condition can be derived from the fact that
        #        the rate of sediment transport between cells k and
        #        p, Qs(k), must be less than the rate one gets by
        #        taking half of the sediment volume at cell k that
        #        is higher than cell p and moving it to cell p in
        #        one timestep.
        #
        #        This condition ensures a well-defined timestep
        #        everywhere except at pits and flats.  It may,
        #        however, allow the elevation at cell p to exceed
        #        the elevation at cell k during a (global) timestep
        #        if there is a sediment flux from other neighbors
        #        into cell k and sufficiently low flux out of cell p.
        #
        #        It could be that with global timesteps, z(k) will
        #        be raised by inflow from its neighbors instead of
        #        lowered by the amount used here.  But the result
        #        should still be valid.
        #--------------------------------------------------------------- 
        if not(SILENT):
            print('Updating dt_grid...')

        pIDs = self.d8.parent_ID_grid.flat[ IDs ]

        ## d8_codes = self.d8.d8_grid.flat[ IDs ]  # (9/14/11)
        
        #------------------------------------
        # Compute downstream elevation drop
        #------------------------------------------------
        # del_z will always be positive as long as
        # LINK_FLATS = False. (flats have D8 code 0)
        #------------------------------------------------
        # NB! del_z will be wrong wherever pID=0 since
        # it will be computed from DEM[0].  Places with
        # pID=0 are fixed further down.
        #------------------------------------------------
        z_IDs  = self.DEM.flat[ IDs ]
        z_pIDs = self.DEM.flat[ pIDs ]
        del_z  = (z_IDs - z_pIDs)
        
            
        #------------------------------
        # Compute stable dt, in years
        #------------------------------
        Qs_IDs = self.Qs.flat[ IDs ]
        w1 = np.where( Qs_IDs > 0 )
        if (w1[0].size > 0):
            if (self.da.size == 1):
                da = self.da
            else:
                da = self.da.flat[ IDs[w1] ]
            ## da = (self.rti.xres * self.rti.yres)
            dt_equal = (da / 2) * del_z[w1] / Qs_IDs[w1]
            self.dt_grid.flat[ IDs[w1] ] = self.CFL_factor * dt_equal

##        w3 = np.where( del_z < 0)
##        if (w3[0].size > 0):
####            print '################# del_z <= 0 at IDs ='
####            print IDs[w3]
####            print 'and time index =', self.time_index
####            print ' '
##            self.dt_grid.flat[ IDs[w3] ] = self.dt_limit
                            
        #----------------------------------------------------
        # Deactivate all pixels that have a flow code of 0.
        # Their z values can then only be changed when they
        # are synchronized by one of their neighbors.
        # This affects edge pixels, flats and pits.
        #----------------------------------------------------
        ## w4 = np.where( d8_codes == 0 )  ##### (9/14/11)
        w4 = np.where( pIDs == 0 )
        if (w4[0].size > 0):
            self.dt_grid.flat[ IDs[w4] ] = self.dt_limit

        #-----------------------------------------------
        # (11/18/10) EXPERIMENT to deal with stalling.
        #-----------------------------------------------------------
        # No change in results or speed vs. using the next block
        # with a WHERE call first.
        #-----------------------------------------------------------
        # (10/25/10) The problem with this is that when/if the
        # pixel gets reactivated (i.e. when synced by a neighbor),
        # then the new dt = (T_clock - T_last) will be large and
        # a pit is likely to form.
        #-----------------------------------------------------------
        if (self.active_ID_count > 10):
            print('###### Progress has stalled at ID =', self.active_ID)
            print('###### Deactivating this cell.')
            print(' ')
            self.dt_grid.flat[ self.active_ID ] = self.dt_limit                            ## self.dt_grid.flat[ ID ] = self.dt_limit  # (DEACTIVATE, doesn't work)
                    
        #--------------------------
        # Save the min and max dt
        #------------------------------------------
        # NB! Don't do min/max over whole grid !!
        #------------------------------------------
        dt_vals = self.dt_grid.flat[ IDs ]
        dt_min  = dt_vals.min()
        dt_max  = dt_vals.max()

##        if (dt_min < self.dt_too_small):
            
        self.dt_min = np.minimum( self.dt_min, dt_min )
        self.dt_max = np.maximum( self.dt_max, dt_max )
        if not(SILENT):
##            dt_vals = self.dt_grid.flat[ IDs ]
##            dt_min  = dt_vals.min()
##            dt_max  = dt_vals.max()
##            self.dt_min = np.minimum( self.dt_min, dt_min )
##            self.dt_max = np.maximum( self.dt_max, dt_max )
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
        if (dt_min < self.dt_too_small):
            if (self.DEBUG):
                print('### WARNING: dt < ' + str(self.dt_too_small))
                print('        dt = ' + str(dt_min))
                print(' ')

            #####################################################
            # NOTE: If "Area grid units" in *_d8.cfg file are
            #       set to "km^2" vs. "m^2", code will crash
            #       pretty quickly at this point. (8/31/11)
            #####################################################
            if (dt_min <= 0):
                print('dt_min =', dt_min)
                sys.exit()
##            print '******************************************'
##            print ' Aborting: Stable dt is too small.'
##            print ' Computed dt = ' + str(dt_grid_min)
##            print '******************************************'
##            print ' '
##            sys.exit()

        #------------------
        # Optional report
        #------------------
        ##### REPORT = True  ###################
        if (REPORT):
            dt_str     = str(self.dt_min)  + ', ' + str(self.dt_max)
            del_z_str  = str(del_z.min())  + ', ' + str(del_z.max())
            ## del_dz_dt_str = str(del_dz_dt.min()) + ', ' + str(del_dz_dt.max())
            print('    min(dt),     max(dt)           = ' + dt_str + ' [yrs]')
            print('    min(del_z),  max(del_z)        = ' + del_z_str)
            ## print '    min(del_dz_dt), max(del_dz_dt) = ' + del_dz_dt_str
            print(' ')
            
    #   update_dt_grid5()
    #-------------------------------------------------------------------
    def update_dt_grid6(self, IDs, SILENT=True, REPORT=False):

        #---------------------------------------------------------------
        # Notes: This version is based on the following stability
        #        condition:
        #
        #        dt < dx * dy * [z(k) - z(p)] / [2 * Qs(k)]
        #
        #        This condition can be derived from the fact that
        #        the rate of sediment transport between cells k and
        #        p, Qs(k), must be less than the rate one gets by
        #        taking half of the sediment volume at cell k that
        #        is higher than cell p and moving it to cell p in
        #        one timestep.
        #
        #        This condition ensures a well-defined timestep
        #        everywhere except at pits and flats.  It may,
        #        however, allow the elevation at cell p to exceed
        #        the elevation at cell k during a (global) timestep
        #        if there is a sediment flux from other neighbors
        #        into cell k and sufficiently low flux out of cell p.
        #
        #        It could be that with global timesteps, z(k) will
        #        be raised by inflow from its neighbors instead of
        #        lowered by the amount used here.  But the result
        #        should still be valid.
        #--------------------------------------------------------------- 
        if not(SILENT):
            print('Updating dt_grid...')

        #-----------------------------------------------
        # Should make sure pIDs are correct. (2/28/12)
        #-----------------------------------------------
        pIDs = self.d8.parent_ID_grid.flat[ IDs ]
        d8_codes = self.d8.d8_grid.flat[ IDs ]  # (9/14/11)
        
        #------------------------------------
        # Compute downstream elevation drop
        #------------------------------------------------
        # del_z will always be positive as long as
        # LINK_FLATS = False. (flats have D8 code 0)
        #------------------------------------------------
        # NB! del_z will be wrong wherever pID=0 since
        # it will be computed from DEM[0].  Places with
        # pID=0 are fixed further down.
        #------------------------------------------------
        z_IDs  = self.DEM.flat[ IDs ]
        z_pIDs = self.DEM.flat[ pIDs ]
        del_z  = (z_IDs - z_pIDs)
    
        #----------------------------
        # Get sed. discharge at IDs
        #----------------------------
        Qs_IDs = self.Qs.flat[ IDs ]
        w1 = np.where( Qs_IDs > 0 )
        
        #------------------------------
        # Compute stable dt, in years
        #------------------------------
        if (w1[0].size > 0):
            if (self.da.size == 1):
                da = self.da
            else:
                da = self.da.flat[ IDs[w1] ]
            dt_equal = (da / 2) * del_z[w1] / Qs_IDs[w1]
            self.dt_grid.flat[ IDs[w1] ] = self.CFL_factor * dt_equal

        #----------------------------------------------------
        # For pits and flats, use dz_dt and max_dz_for_pits
        # to set dt.  They are now put on heap. (3/7/12)
        # Edge pixels remain deactivated (from beginning).
        #----------------------------------------------------
        test1 = (d8_codes == 0)
        test2 = self.not_edge_grid.flat[ IDs ]
        w2 = np.where( np.logical_and( test1, test2 ) )
        if (w2[0].size > 0):
            abs_dz_dt = np.absolute( self.dz_dt.flat[ IDs[w2] ] )
            self.dt_grid.flat[ IDs[w2] ] = self.max_dz_for_pits / abs_dz_dt

        
        #----------------------------------------------------
        # Deactivate all pixels that have a flow code of 0.
        # Their z values can then only be changed when they
        # are synchronized by one of their neighbors.
        # This affects edge pixels, flats and pits.
        # (Used before 3/7/12 experiment.)
        #----------------------------------------------------
        # pIDs, d8_codes and Qs should be zero at same IDs.
        #----------------------------------------------------        
##        ## w4 = np.where( pIDs == 0 )      
##        w4 = np.where( d8_codes == 0 )
##        if (w4[0].size > 0):
##            self.dt_grid.flat[ IDs[w4] ] = self.dt_limit

        #------------------------------------------------------
        # (2/29/12) Stalling did occur again at a few pixels
        # when test4 was removed from update_neighbors().
        #------------------------------------------------------
        # (2/28/12) There hasn't been a problem with stalling
        # for a while now.  Eventually remove this test.
        #--------------------------------------------------------
        # (11/18/10) EXPERIMENT to deal with "stalling", where
        # the same active_ID is repeatedly drawn from the heap.
        #-----------------------------------------------------------
        # No change in results or speed vs. using the next block
        # with a WHERE call first.
        #-----------------------------------------------------------
        # (10/25/10) The problem with this is that when/if the
        # pixel gets reactivated (i.e. when synced by a neighbor),
        # then the new dt = (T_clock - T_last) will be large and
        # a pit is likely to form.
        #-----------------------------------------------------------
        if (self.active_ID_count > 10):
            print('###### Progress has stalled at ID =', self.active_ID)
            print('###### Deactivating this cell.')
            print(' ')
            self.dt_grid.flat[ self.active_ID ] = self.dt_limit
                    
        #--------------------------
        # Save the min and max dt
        #------------------------------------------
        # NB! Don't do min/max over whole grid !!
        #------------------------------------------
        dt_vals = self.dt_grid.flat[ IDs ]
        dt_min  = dt_vals.min()
        dt_max  = dt_vals.max()        
        self.dt_min = np.minimum( self.dt_min, dt_min )
        self.dt_max = np.maximum( self.dt_max, dt_max )
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
        if (dt_min < self.dt_too_small):
            if (self.DEBUG):
                print('### WARNING: dt < ' + str(self.dt_too_small))
                print('        dt = ' + str(dt_min))
                print(' ')

            #####################################################
            # NOTE: If "Area grid units" in *_d8.cfg file are
            #       set to "km^2" vs. "m^2", code will crash
            #       pretty quickly at this point. (8/31/11)
            #####################################################
            if (dt_min <= 0):
                print('dt_min =', dt_min)
                sys.exit()
##            print '******************************************'
##            print ' Aborting: Stable dt is too small.'
##            print ' Computed dt = ' + str(dt_grid_min)
##            print '******************************************'
##            print ' '
##            sys.exit()

        #------------------
        # Optional report
        #------------------
        ##### REPORT = True  ###################
        if (REPORT):
            dt_str     = str(self.dt_min)  + ', ' + str(self.dt_max)
            del_z_str  = str(del_z.min())  + ', ' + str(del_z.max())
            ## del_dz_dt_str = str(del_dz_dt.min()) + ', ' + str(del_dz_dt.max())
            print('    min(dt),     max(dt)           = ' + dt_str + ' [yrs]')
            print('    min(del_z),  max(del_z)        = ' + del_z_str)
            ## print '    min(del_dz_dt), max(del_dz_dt) = ' + del_dz_dt_str
            print(' ')
            
    #   update_dt_grid6()
    #-------------------------------------------------------------------
    def update_dz_target_grid2(self, IDs, SILENT=True, REPORT=False):

        #------------------------------------------------------------
        # Note:  dz_target is the change to the DEM that we will
        #        make at a given ID if DEM[ID] is updated at its
        #        scheduled time, T_next[ID].  However, if DEM[ID]
        #        is updated sooner (and using dz_dt[ID]), then
        #        the change to DEM[ID] will be smaller.
        #------------------------------------------------------------
        # Note:  dt_grid was computed with self.CFL_factor.
        #------------------------------------------------------------
        # Note:  dz_target is allowed to keep its sign.
        #------------------------------------------------------------
        # Note:  For IDs that have been "deactivated" by setting
        #        dt_grid[ID] = self.dt_limit, dz_target may be
        #        very large.  (This now happens for edge pixels.)
        #        In this case, the amount that accumulates in
        #        the flux capacitor is unlikely to ever exceed
        #        dz_target.  So they will never be "fully updated"
        #        as "real events" (from the heap) and can only
        #        have their z-value changed during synchronization
        #        with an interior neighbor.  Recall that OK(2006)
        #        paper refers to "active interior cells" (p. 185).
        #------------------------------------------------------------       
        dt    = self.dt_grid.flat[ IDs ]  # (= dt_limit where d8=0)
        dz_dt = self.dz_dt.flat[ IDs ]
        dz    = dt * dz_dt
        self.dz_target_grid.flat[ IDs ] = dz

        #------------------------------------------------------------
        # Set dz_target to 0 for pixels that have been deactivated.
        # This will allow them to get "reactivated" when they have
        # valid flow directions later on.  (Due to flux_cap test.)
        # (10/1/10).  Checked if still helpful on 10/21/10.
        #------------------------------------------------------------
        # Cells with undefined flow directions may have been
        # generated from previous updates.  Now, update_dt_grid5()
        # is called before this and dt_grid is set to dt_limit
        # wherever pID == 0 (or d8_codes == 0).
        #------------------------------------------------------------

        #############################################################
        #############################################################
        # NB! Before 2/28/12, the commented line below was used to
        #     determine cells that should have dz_target = 0) and
        #     these were used in update_neighbors() to determine
        #     which cells to include in "now_IDs".  However, the
        #     new check_zone_info() function then showed that it is
        #     possible to have cells (e.g, when activeID = 2276)
        #     with (dt != dt_limit) at an *interior cell* that has
        #     (d8_code == 0).
        #
        #     In update_dt_grid6(), cells with (d8_code == 0) are
        #     assigned (dt = dt_limit). But the result above shows
        #     that some cells can become pits or flats *after* this
        #     assignment.
        #
        #     The commented line below was setting (dz_target = 0)
        #     for cells with (dt == dt_limit), but was therefore
        #     not including all cells with (d8_code == 0).
        #
        #     In update_neighbors(), the test (dz_target == 0) was
        #     being used to find interior cells with (d8_code == 0),
        #     which because of the above problem was not valid.
        #     Now, update_neighbors() uses a more direct test to
        #     find interior (non-edge) pixels with (d8_code == 0)
        #     to be included in now_IDs.
        #
        #     We do not want to include edge (boundary) cells with
        #     now_IDs, because (1) the OK(2006)paper says not to,
        #     and (2) when a cell one away from the edge (e.g.
        #     ID = 134) gets updated, the set of "now_IDs" will
        #     iteratively grow to include the entire DEM boundary.
        #
        # w = np.where( self.dt_grid.flat[ IDs ] == self.dt_limit )
        #############################################################
        #############################################################

        #-----------------------------------------------------------
        # (3/7/12) Experiment that assigns dt-values to pits/flats
        # and puts them on heap. (Comment out this block.)
        #-----------------------------------------------------------
##        test1 = (self.d8.d8_grid.flat[ IDs ] == 0)
##        test2 = self.not_edge_grid.flat[ IDs ]
##        w = np.where( np.logical_and(test1, test2) )
##        if (w[0].size > 0):
##            self.dz_target_grid.flat[ IDs[w] ] = 0
                
        #-----------------------------------------------------------
        # All cells with (d8_code == 0) are assigned a dt_grid
        # value of dt_limit in update_dt_grid6().  This leads to
        # a very large value of dz_target above, which is set down
        # to zero for the pits and flats but remains large for the
        # edge cells.  On the edges, we will therefore have
        # abs(dz_cap) < abs(dz_target).  See update_neighbors().
        #-----------------------------------------------------------
        
    #   update_dz_target_grid2()
    #-------------------------------------------------------------------
    def update_T_next_grid(self, IDs, SILENT=True, REPORT=False):
             
        #----------------------------------------------------------        
        # Note: These notes are somewhat outdated. (9/14/10)
        #
        #       The 8 neighbors already have a scheduled update
        #       time in T_next, but their state variables have
        #       changed as a result of changing DEM[active_ID].
        #
        #       The change to their z-value that was scheduled to
        #       occur (using previous state variables) has been
        #       saved in the flux_cap_grid.  But this change may
        #       be too large because their update time hasn't
        #       been reached yet.  So we need to know the dz that
        #       would have occurred up to now.
        #
        #       Compute their "current" dt = (T_clock - T_last).
        #       Compute the dz that would have occurred over this
        #           dt using state vars before they were updated.
        #       Maybe just apply this dz and assign each of them
        #           a new T_next.  Then we don't need the flux
        #           capacitor ??
        #       Or save this dz in flux capacitor, assign them a
        #           new T_next based on the just-updated state
        #           and when T_next is reached, apply both the
        #           usual dz update and whatever has been saved
        #           up in the flux capacitor.
        #
        #----------------------------------------------------------
        #       Note that dz_dt integrates all fluxes into and
        #       out of the cell.  Maybe we should be computing
        #       the dz from each neighbor separately using the
        #       corresponding Qs value.
        #----------------------------------------------------------

        #--------------------------------------------------------
        # Update T_next grid for active_ID because it was just
        # updated at time "T_clock".  Note that we rely on
        # T_next grid being up-to-date in order to tell when
        # events have been retracted. See update_active_ID().
        #--------------------------------------------------------       
        dt_vals = self.dt_grid.flat[ IDs ]

        #------------------------------------------------
        # This test was trigged immediately on 2/28/12,
        # so rewrote next section as shown.
        #------------------------------------------------
##        w = np.where( dt_vals >= self.dt_limit )
##        if (w[0].size > 0):
##            print '###########################################'
##            print ' ERROR in update_T_next_grid():'
##            print '       Deactivated pixel put on heap.'
##            print '###########################################'
##            print ' '
##            sys.exit()
            
        #---------------------------------------------------
        # (2/28/12) Don't put things on the heap that have
        # dt >= dt_limit.  This will reduce heap size and
        # perhaps speed things up.
        #---------------------------------------------------
##        w = np.where( dt_vals < self.dt_limit )
##        if (w[0].size == 0):
##            return

        #---------------------------------------------------
        # (2/28/12) Don't put things on the heap that have
        # (d8_code == 0).  The stability condition does not
        # allow us to compute a dt_grid value for them, so
        # they currently get a value of dt_limit, which
        # means they will never be drawn from the heap.
        # This will reduce heap size and perhaps speed
        # things up.
        #---------------------------------------------------
        # NB! Tried this block instead of the one above to
        # resolve the "zmax just increased at time_index
        # 895402" error but got the same exact result.
        # Final heap size also 35635 for both.
        #---------------------------------------------------
        d8_codes = self.d8.d8_grid.flat[ IDs ]

        #---------------------------------------------------
        # (3/5/12) Didn't have this block before today.
        # It may not be needed, but without it T_next
        # may still be set to an old and inaccurate
        # value.  The value is used to check whether an
        # event has been RETRACTED when drawing from heap.
        #---------------------------------------------------
##        w1 = np.where( d8_codes == 0)
##        if (w1[0].size > 0):
##            self.T_next.flat[ IDs[w1] ] = self.dt_limit                   
##        w2 = np.where( d8_codes != 0 )
##        if (w2[0].size == 0):
##            return    
##        dt_vals = dt_vals[ w2 ]
##        new_IDs = IDs[ w2 ]

        #-------------------------------------------
        # (3/7/12) Allow pits & flats to be put on
        # the heap, but not edges.
        #-------------------------------------------
        w1 = np.where( self.not_edge_grid.flat[ IDs ] )
        if (w1[0].size == 0):
            return
        dt_vals = dt_vals[ w1 ]
        new_IDs = IDs[ w1 ]

        #----------------------------------------
        # Put IDs on heap with new update times
        #----------------------------------------
        T_new = (self.T_clock + dt_vals)
        self.T_next.flat[ new_IDs ] = T_new
        self.n_heap += new_IDs.size
        for k in range(new_IDs.size):
            heapq.heappush( self.T_heap, (T_new[k], new_IDs[k]) )

        #---------------
        # For testing
        #---------------
        if (REPORT):
            ID = new_IDs[0]
            print('ID         =', ID)
            print('A[ID]      =', self.d8.A.flat[ ID ])
            print('S[ID]      =', self.S.flat[ ID ])
            print('z[ID]      =', self.DEM.flat[ ID ])
            print('T_next[ID] =', self.T_next.flat[ ID ])
            print('dt[ID]     =', self.dt_grid.flat[ ID ])
            print('dz_dt[ID]  =', self.dz_dt.flat[ ID ])     
            print('T_clock    =', self.T_clock)
            print(' ')
        
    #   update_T_next_grid()
    #-------------------------------------------------------------------
    def save_update_count_grid(self):
        
        count_file = (self.out_directory + self.case_prefix + '_upCount.rtg')
        rtg_files.write_grid( self.update_count, count_file, self.rti,
                              RTG_type='LONG', var_name='update_count' )

        total_updates = np.sum( self.update_count )
        avg_n_updates = (total_updates / self.rti.n_pixels)
        print(' ')
        print('Min number of cell updates   =', self.update_count.min())
        print('Max number of cell updates   =', self.update_count.max())
        print('Total number of cell updates =', total_updates)
        print('Average updates per cell     =', avg_n_updates)
        print(' ')

    #   save_update_count_grid()
    #-------------------------------------------------------------------
    def save_d8_grids(self):

        #----------------------------------
        # Save the D8 flow direction grid
        #----------------------------------
        d8_file = (self.out_directory + self.case_prefix + '_flow.rtg')
        rtg_files.write_grid( self.d8.d8_grid, d8_file, self.rti,
                              RTG_type='BYTE', var_name='d8_code')

        #-------------------------------------
        # Save the D8 contributing area grid
        #-------------------------------------        
        area_file = (self.out_directory + self.case_prefix + '_area.rtg')
        rtg_files.write_grid( self.d8.A, area_file, self.rti,
                              RTG_type='FLOAT', var_name='d8_area')
            
    #   save_d8_grids()
    #-------------------------------------------------------------------
    def save_final_grids(self):

        #---------------------------
        # Save "update count" grid
        #---------------------------
        # self.save_update_count_grid()  # (now called in finalize().)

        #-----------------------------------
        # Save last D8 flow direction grid
        #-----------------------------------
        d8_file = (self.out_directory + self.case_prefix + '_lastflow.rtg')
        rtg_files.write_grid( self.d8.d8_grid, d8_file, self.rti,
                              RTG_type='BYTE', var_name='d8_code')

        #--------------------------------------
        # Save last D8 contributing area grid
        #--------------------------------------        
        area_file = (self.out_directory + self.case_prefix + '_lastA.rtg')
        rtg_files.write_grid( self.d8.A, area_file, self.rti,
                              RTG_type='FLOAT', var_name='d8_area')
        
        #---------------------
        # Save the final DEM
        #---------------------
        DEM_file = (self.out_directory + self.case_prefix + '_lastZ.rtg')
        rtg_files.write_grid( self.DEM, DEM_file, self.rti,
                              RTG_type='FLOAT', var_name='elev')

        #----------------------------
        # Save the final slope grid
        #----------------------------       
        S_file = (self.out_directory + self.case_prefix + '_lastS.rtg')
        rtg_files.write_grid( self.S, S_file, self.rti,
                              RTG_type='FLOAT', var_name='d8_slope')

        #--------------------------------
        # Save the final discharge grid
        #--------------------------------       
        Q_file = (self.out_directory + self.case_prefix + '_lastQ.rtg')
        rtg_files.write_grid( self.Q, Q_file, self.rti,
                              RTG_type='FLOAT', var_name='discharge')

        #-------------------------------------
        # Save the final sed. discharge grid
        #-------------------------------------       
        Qs_file = (self.out_directory + self.case_prefix + '_lastQs.rtg')
        rtg_files.write_grid( self.Qs, Qs_file, self.rti,
                              RTG_type='FLOAT', var_name='sed_discharge')

        #------------------------
        # Save the final D grid
        #------------------------       
        D_file = (self.out_directory + self.case_prefix + '_lastD.rtg')
        rtg_files.write_grid( self.D, D_file, self.rti,
                              RTG_type='FLOAT', var_name='diff_coeff')

        #----------------------------
        # Save the final dz_dt grid
        #----------------------------       
        dz_dt_file = (self.out_directory + self.case_prefix + '_lastdzdt.rtg')
        rtg_files.write_grid( self.dz_dt, dz_dt_file, self.rti,
                              RTG_type='FLOAT', var_name='dz_dt')

        #-------------------------
        # Save the final dt grid
        #-------------------------   
        dt_grid_file = (self.out_directory + self.case_prefix + '_lastdt.rtg')
        rtg_files.write_grid( self.dt_grid, dt_grid_file, self.rti,
                              RTG_type='FLOAT', var_name='dt')
        
        #-----------------------------
        # Save the final T_last grid
        #-----------------------------   
        T_last_file = (self.out_directory + self.case_prefix + '_lastTlast.rtg')
        rtg_files.write_grid( self.T_last, T_last_file, self.rti,
                              RTG_type='FLOAT', var_name='T_last')

        #-----------------------------
        # Save the final T_next grid
        #-----------------------------   
        T_next_file = (self.out_directory + self.case_prefix + '_lastTnext.rtg')
        rtg_files.write_grid( self.T_next, T_next_file, self.rti,
                              RTG_type='FLOAT', var_name='T_next')
        
    #   save_final_grids()
    #-------------------------------------------------------------------

