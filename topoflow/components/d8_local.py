
## Copyright (c) 2001-2013, Scott D. Peckham
## January 2009  (converted from IDL)
## May, July, August, October 2009
## March 2010 (modified to support local timestepping)
## September 2011  (Updated unit tests.)
## February 2012.  (Updated unit_test2(), but link_flats() still
##                  doesn't work correctly for KY_Sub case.)

## Local and Global versions for KY_Sub are now the same
##     including 576 cells near edges that could not be resolved
##     but were not set to flow code of 0.  Same is true for
##     Treynor and Beaver tests.
## Note that linked flats do not have to agree perfectly with
##     those from RT3.

## initialize_computed_vars() reads the DEM from DEM_file,  ################
## but if FILL_PITS_IN_Z0, it actually OVERWRITES that file
## This still needs to be addressed.  (2/23/12)

#---------------------------------------------------------
# Note: This version (March 2010) has been modified to
#       support local timestepping in the Erode model.
#       The main idea is to remove the constraint that
#       calculations are made on grids.  That is, this
#       version allows changes to be made for a subset
#       of the pixels within the domain of the DEM,
#       and then updated in the corresponding grid.
#       This does not require extensive changes due to
#       the way the code was written previously.
#---------------------------------------------------------
#       self.IDs      = pixel IDs to be updated
#                       (rows, cols) from divmod()
#       self.d8_codes = D8 codes at those IDs
#---------------------------------------------------------
#       There is a need to compute an entire flow grid
#       in the initialize() method.  This is done
#       by setting self.IDs to be all of the IDs in the
#       grid and then calling the update_*() methods.
#       See initialize_computed_vars().
#---------------------------------------------------------
## NB!  TF expects d8.codes vs. d8.flow_grid

# from numpy import *    ## ELIMINATE THIS.

import numpy as np

import os        # (for os.chdir, in unit_test())
import os.path
import platform  # (for get_test_info())
import sys
import time

from topoflow.components import d8_base
from topoflow.utils      import fill_pits
from topoflow.utils      import rtg_files

#-------------------------------------------
# For use outside of the TopoFlow package.
#-------------------------------------------
# import d8_base
# import fill_pits
# import rtg_files

#---------------------------------------------------------------------
#
#   class d8_local   (inherits from d8_base.py)
#
#       get_component_name()
#       get_attribute()             # (10/27/11)
#       initialize_computed_vars()
#------------------------------------
#       get_neighbors()
#       is_a_pit()
#       is_an_interior_pit()
#-------------------------------------
#       These next few are old and
#       should be retested. (3/1/12)
#-------------------------------------
#       update_neighborhood_IDs()   # (can do periodic BCs, 4/22/10)
#       update_kid_IDs_BAD()        # (isn't working)
#       update_kid_IDs()            # (can't do periodic BCs, 9/20/10)
#       get_kid_IDs()               # (same as above)
#       add_neighbors_to_IDs()      # (can do periodic BCs, 3/13/10)
#----------------------------------
#       update_flow_grid()
#          start_new_d8_codes()
#          break_flow_grid_ties()
#          link_flats()
#----------------------------------
#       update_parent_ID_grid()     # (can do periodic BCs) 
#       update_parent_IDs()         # (used by erode_d8_local.update_slope_grid())
#       update_flow_from_IDs()
#       update_flow_to_IDs()
#       update_noflow_IDs()         # (untested and unused)
#----------------------------------
#       update_flow_width_grid()
#       update_flow_length_grid()
#       update_area_grid()
#
#-----------------------------------------------------------------------
class d8_component(d8_base.d8_component):

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_D8_Local'

    #   get_component_name() 
    #-----------------------------------------------------------------
    # Note: Do not define an __init__() method here.  It will
    #       override things needed from CSDMS_base.__init__().
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):
    
        map = {'comp_name':          'D8Local',
               'version':            '0.5',
               'model_name':         'D8_Local',
               'model_family':       'Erode',
               'cfg_template_file':  'D8_Local.cfg.in',
               'cfg_extension':      '_d8_local.cfg',
               'cmt_var_prefix':     '/D8Local/Input/Var/',
               'gui_xml_file':       '/home/csdms/cca/erode/0.5/src/share/cmt/gui/D8_Local.xml',
               'dialog_title':       'D8: Local-Timestep Method Parameters',
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
    #-----------------------------------------------------------------
    def initialize_computed_vars(self, DOUBLE=False,
                                 SILENT=True, REPORT=False):
        
        #--------------------------------------------------------
        # Notes:  Need to compute these from d8 flow grid, then
        #         they'll only be changed as necessary during
        #         the "update_*" calls.
        #--------------------------------------------------------
        self.RT3_TEST          = False  ## (default settings)
        self.PERTURB_TEST      = False
        ### self.AREA_GRID_WARNING = False  ## (11/20/11)        
        self.AREA_GRID_WARNING = True   ### (11/20/11)
        
        nx = self.nx  # (Local synonyms)
        ny = self.ny
        self.IDs = None
        if (DOUBLE):
            dtype = 'Float64'
        else:
            dtype = 'Float32'

        #--------------------------------------------------
        # Initialize the grids, due to how d8_local works
        #--------------------------------------------------
        self.d8_grid        = np.zeros([ny, nx], dtype='Int16')
        self.parent_ID_grid = np.zeros([ny, nx], dtype='Int32')
        self.dw = np.zeros([ny, nx], dtype=dtype)
        self.ds = np.zeros([ny, nx], dtype=dtype)
        self.A  = np.zeros([ny, nx], dtype=dtype)

        #------------------------------------------------------
        # These were previously used by get_neighbor_IDs()
        # but fixed bug and now next ones are used. (9/14/11)
        #------------------------------------------------------
        self.row_incs = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1], dtype='Int32')
        self.col_incs = np.array([-1, 0, 1, -1, 0, 1, -1, 0, 1], dtype='Int32')

        #-----------------------------------------------
        # (9/20/10) These get used by update_kid_IDs()
        #-----------------------------------------------
        self.row_incs_CW = np.array([-1,0,1,1,1,0,-1,-1], dtype='Int32')
        self.col_incs_CW = np.array([1,1,1,0,-1,-1,-1,0], dtype='Int32')


        ###############################################
        # Copied remaining code from same method
        # in d8_base.py. (11/20/11)
        ##############################################################
        # BUT, note that erode_base.initialize_DEM() has a READ_FILE
        # option or can call its create_initial_DEM() method.  That
        # happens BEFORE this does, so if DEM_file is found it
        # replaces the DEM from erode_d8_local.py. It is not clean
        # or clear to set DEM_file = '' in the CFG file for this
        # case.  Should probably set a flag of some kind.
        ##############################################################
        if (self.DEM_file == 'NOT_USED'):
            return
        
        #-----------------------------------
        # Read DEM from DEM_file (11/7/11)
        #----------------------------------------------
        # This was removed from start_new_d8_codes().
        #
        # NB! start_new_d8_codes() has a DEM argument
        # that is used by erode_d8_global.py, etc.
        #----------------------------------------------
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
    #-------------------------------------------------------------------
    def get_neighbors(self, IDs, INCLUDE_IDS=False,
                      SILENT=True, REPORT=False):

        #-------------------------------------------------
        # Notes: This doesn't handle periodic BCs yet.
        #-------------------------------------------------
        # (8/23/10, 9/3/10) Original versions.
        #-------------------------------------------------        
        # (9/14/11) Wrote get_neighbors_test() and then
        # found & fixed a serious bug.  Was using only
        # the first 8 of d8.row_incs and d8.col_incs, but
        # there are 9 of each (center is included). Only
        # affected case where IDs.size > 1.
        #-------------------------------------------------
        # (3/1/12) Moved here from erode_d8_local.py and
        # now called by update_flow_grid().
        #-------------------------------------------------        
        rows, cols = divmod( IDs, self.nx )
       
        if (IDs.size == 1):
            #--------------------------------------
            # Includes center pixel & 8 neighbors
            #--------------------------------------
            rows2 = (rows + self.row_incs)
            cols2 = (cols + self.col_incs)
        else:
            #---------------------------------------
            # Includes center pixels & 8 neighbors
            #------------------------------------------------
            # Center row and col (r4 and c4) were excluded
            # before 2/29/12 since IDs were removed at end.
            # But need them for new INCLUDE_IDS option.
            #------------------------------------------------
            r0 = rows + self.row_incs[0]
            r1 = rows + self.row_incs[1]
            r2 = rows + self.row_incs[2]
            r3 = rows + self.row_incs[3]
            r4 = rows + self.row_incs[4]  ## (See note above.)
            r5 = rows + self.row_incs[5]
            r6 = rows + self.row_incs[6]
            r7 = rows + self.row_incs[7]
            r8 = rows + self.row_incs[8]
            #### See note above about r4.
            rows2 = np.concatenate([r0,r1,r2,r3,r4,r5,r6,r7,r8])
            ## rows2 = np.concatenate([r0,r1,r2,r3, r5,r6,r7,r8])
            #---------------------------------------------------------
            c0 = cols + self.col_incs[0]
            c1 = cols + self.col_incs[1]
            c2 = cols + self.col_incs[2]
            c3 = cols + self.col_incs[3]
            c4 = cols + self.col_incs[4]  ## (See note above.)
            c5 = cols + self.col_incs[5]
            c6 = cols + self.col_incs[6]
            c7 = cols + self.col_incs[7]
            c8 = cols + self.col_incs[8]
            #### See note above about c4.
            cols2 = np.concatenate([c0,c1,c2,c3,c4,c5,c6,c7,c8])
            ## cols2 = np.concatenate([c0,c1,c2,c3, c5,c6,c7,c8])
            
##            #-----------------------------------
##            # Includes 8 neighbors in CW order
##            #-----------------------------------
##            r0 = rows + self.row_incs_CW[0]
##            r1 = rows + self.row_incs_CW[1]
##            r2 = rows + self.row_incs_CW[2]
##            r3 = rows + self.row_incs_CW[3]
##            r4 = rows + self.row_incs_CW[4]
##            r5 = rows + self.row_incs_CW[5]
##            r6 = rows + self.row_incs_CW[6]
##            r7 = rows + self.row_incs_CW[7]
##            rows2 = np.concatenate([r0,r1,r2,r3,r4,r5,r6,r7])
##            #------------------------------------------------------
##            c0 = cols + self.col_incs_CW[0]
##            c1 = cols + self.col_incs_CW[1]
##            c2 = cols + self.col_incs_CW[2]
##            c3 = cols + self.col_incs_CW[3]
##            c4 = cols + self.col_incs_CW[4]
##            c5 = cols + self.col_incs_CW[5]
##            c6 = cols + self.col_incs_CW[6]
##            c7 = cols + self.col_incs_CW[7]
##            cols2 = np.concatenate([c0,c1,c2,c3,c4,c5,c6,c7])

        #--------------------------------------
        # Exclude rows and cols out of range
        # better than method below. (9/12/10)
        #--------------------------------------
        w1 = np.where( np.logical_and( rows2 > -1, rows2 < self.ny ) )
        if (w1[0].size > 0):
            rows2 = rows2[ w1 ]
            cols2 = cols2[ w1 ]
        w2 = np.where( np.logical_and( cols2 > -1, cols2 < self.nx ) )
        if (w2[0].size > 0):
            rows2 = rows2[ w2 ]
            cols2 = cols2[ w2 ]
        
        IDs2 = (rows2 * self.nx) + cols2
        
        #------------------------------------
        # Exclude IDs out of range (9/2/10)
        #-------------------------------------------
        # But this gives neighbors of 0 to be
        # [1, 19, 20, 21] instead of [1, 20, 21 ].
        #-------------------------------------------        
##        w1 = np.where( np.logical_and( IDs2 < self.rti.n_pixels, IDs2 >= 0 ) )
##        IDs2 = IDs2[ w1 ]

        #----------------------------------------------------
        # (2/29/12) New option for use in update_d8_vars().
        #----------------------------------------------------
        if (INCLUDE_IDS):
            return np.unique( IDs2 )
        
        #---------------------------
        # Exclude the original IDs
        #--------------------------------------------------------
        # np.setdiff1d removes the IDs in the 2nd argument
        # from those in the 1st argument.  It returns "[]"
        # (which has size 0) if the two arguments are the same.
        # The "assume_unique" keyword is not in all versions.
        #--------------------------------------------------------
        ## nbr_IDs = np.setdiff1d( IDs2, IDs, assume_unique=False )
        ## nbr_IDs = np.setdiff1d( IDs2, IDs )
        if (np.ndim(IDs) == 0):
            nbr_IDs = np.setdiff1d( IDs2, np.array([IDs]) )
        else:
            nbr_IDs = np.setdiff1d( IDs2, IDs )
        
        #-------------------------------------------
        # Exclude the boundary IDs
        # (9/6/10) See update_neighbors() instead.
        #-------------------------------------------
##        if (nbr_IDs.size == 0):
##            return nbr_IDs
##        nbr_IDs = np.setdiff1d( nbr_IDs, self.base_IDs )
        
        return nbr_IDs
        
    #   get_neighbors()
    #-------------------------------------------------------------------
    def is_a_pit(self, ID, SILENT=True, REPORT=False):

        nbr_IDs = self.get_neighbors( ID )
        z_nbrs  = self.DEM.flat[ nbr_IDs ]
        z_ID    = self.DEM.flat[ ID ]
        w       = np.where( z_nbrs <= z_ID )

        return (w[0].size == 0)
        
    #   is_a_pit()
    #-------------------------------------------------------------------
    def is_an_interior_pit(self, ID, SILENT=True, REPORT=False):

        row, col = divmod( ID, self.nx )

        if (row == 0) or (row == self.ny-1) or \
           (col == 0) or (col == self.nx-1):
            return False
        
        nbr_IDs = self.get_neighbors( ID )
        z_nbrs  = self.DEM.flat[ nbr_IDs ]
        z_ID    = self.DEM.flat[ ID ]
        w       = np.where( z_nbrs <= z_ID )

        return (w[0].size == 0)
        
    #   is_an_interior_pit()
    #-------------------------------------------------------------------    
    #-------------------------------------------------------------------
    def update_neighborhood_IDs(self, ID):

        #-------------------------------------
        # Get IDs of 8 neighbor pixels of ID
        #-----------------------------------------------
        # Allow periodic BCs here via use of "%" (mod)
        #-----------------------------------------------    
        nx  = self.nx
        ny  = self.ny
        r0, c0 = divmod( ID, nx)

        #-------------------------------------------------
        # All periodic BCs in the left-right direction ?
        #-------------------------------------------------
        cols = (c0 + self.col_incs) % nx
        if not(self.LR_PERIODIC):
            if (c0 == 0):
                w = np.where( self.col_incs != -1)
                cols = cols[ w ]
            if (c0 == (nx-1)):
                w = np.where( self.col_incs != 1)
                cols = cols[ w ]

        #-------------------------------------------------
        # All periodic BCs in the top-bottom direction ?
        #-------------------------------------------------
        rows = (r0 + self.row_incs) % ny
        if not(self.TB_PERIODIC):
            if (r0 == 0):
                w = np.where( self.row_incs != -1)
                rows = rows[ w ]
            if (r0 == (ny-1)):
                w = np.where( self.row_incs != 1)
                rows = rows[ w ]
                
        #-------------------------------------------------
        # Compute the calendar-style IDs, 3x3 block
        #-------------------------------------------------                    
        IDs  = (rows * nx) + cols
        self.near_IDs = IDs  ##############

    #   update_neighborhood_IDs()
    #-------------------------------------------------------------------
    def update_kid_IDs_BAD(self, ID):
     
        #-----------------        
        # Local synonyms
        #-----------------
        nx = self.nx
        ny = self.ny
        r0, c0 = divmod( ID, nx)

        #-------------------------------------------------
        # All periodic BCs in the left-right direction ?
        #-------------------------------------------------
        cols = (c0 + self.col_incs) % nx
        if not(self.LR_PERIODIC):
            if (c0 == 0):
                w = np.where( self.col_incs != -1)
                cols = cols[ w ]
            if (c0 == (nx-1)):
                w = np.where( self.col_incs != 1)
                cols = cols[ w ]

        #-------------------------------------------------
        # All periodic BCs in the top-bottom direction ?
        #-------------------------------------------------
        rows = (r0 + self.row_incs) % ny
        if not(self.TB_PERIODIC):
            if (r0 == 0):
                w = np.where( self.row_incs != -1)
                rows = rows[ w ]
            if (r0 == (ny-1)):
                w = np.where( self.row_incs != 1)
                rows = rows[ w ]

        #------------------------------------------
        # Find the D8 kids, if any, as the set of
        # neighbor IDs with flow codes towards ID
        #------------------------------------------
        ID_vals = (rows * nx) + cols
        ## d8_vals = self.d8_grid[ rows, cols ]
        d8_vals = self.d8_grid.flat[ ID_vals ]
        ## h_vals  = self.code_opps_3x3.flat
        h_vals  = np.ravel(self.code_opps_3x3)
        print('ID_vals =', ID_vals)
        print('d8_vals =', d8_vals)
        print('h_vals  =', h_vals)
        w       = np.where( d8_vals == h_vals )
        n_kids  = w[0].size
        if (n_kids > 0):
            IDs = ID_vals[ w ]
        else:
            IDs = np.int32([])   # (empty list, size 0)

        self.kid_IDs = IDs
        
    #   update_kid_IDs_BAD()
    #-------------------------------------------------------------------
    def update_kid_IDs(self, ID):
     
        #-----------------        
        # Local synonyms
        #-----------------
        nx = self.nx
        ny = self.ny
        r0, c0 = divmod( ID, nx)

        cols  = (c0 + self.col_incs_CW)   # (NEED CW VERSION)
        rows  = (r0 + self.row_incs_CW)   # (NEED CW VERSION)
        h     = self.code_opps

        #----------------------------------------------
        # Exclude rows and cols out of range.
        # Same as in erode_d8_local8.get_neighbors().
        #----------------------------------------------
        w1 = np.where( np.logical_and( rows > -1, rows < self.ny ) )
        if (w1[0].size > 0):
            rows = rows[ w1 ]
            cols = cols[ w1 ]
            h    = h[ w1 ]
        w2 = np.where( np.logical_and( cols > -1, cols < self.nx ) )
        if (w2[0].size > 0):
            rows = rows[ w2 ]
            cols = cols[ w2 ]
            h    = h[ w2 ]
            
        #------------------------------------------
        # Find the D8 kids, if any, as the set of
        # neighbor IDs with flow codes towards ID
        #------------------------------------------
        ID_vals = (rows * nx) + cols
        d8_vals = self.d8_grid.flat[ ID_vals ]
##        print 'ID_vals =', ID_vals
##        print 'd8_vals =', d8_vals
##        print 'h_vals  =', h
##        print ' '
        w       = np.where( d8_vals == h )   #### DOUBLE-CHECK THIS.
        n_kids  = w[0].size
        if (n_kids > 0):
            self.kid_IDs = ID_vals[ w ]
        else:
            self.kid_IDs = np.int32([])   # (empty list, size 0)
        
    #   update_kid_IDs()
    #-------------------------------------------------------------------
    def get_kid_IDs(self, ID):

        #--------------------------------------------------------
        # Notes:  This is identical to update_kid_IDs(), except
        #         it returns the kid_IDs directly instead of
        #         storing them into self.
        #--------------------------------------------------------
        
        #-----------------        
        # Local synonyms
        #-----------------
        nx = self.nx
        ny = self.ny
        r0, c0 = divmod( ID, nx)

        cols  = (c0 + self.col_incs_CW)   # (NEED CW VERSION)
        rows  = (r0 + self.row_incs_CW)   # (NEED CW VERSION)
        h     = self.code_opps

        #----------------------------------------------
        # Exclude rows and cols out of range.
        # Same as in erode_d8_local8.get_neighbors().
        #----------------------------------------------
        w1 = np.where( np.logical_and( rows > -1, rows < self.ny ) )
        if (w1[0].size > 0):
            rows = rows[ w1 ]
            cols = cols[ w1 ]
            h    = h[ w1 ]
        w2 = np.where( np.logical_and( cols > -1, cols < self.nx ) )
        if (w2[0].size > 0):
            rows = rows[ w2 ]
            cols = cols[ w2 ]
            h    = h[ w2 ]
            
        #------------------------------------------
        # Find the D8 kids, if any, as the set of
        # neighbor IDs with flow codes towards ID
        #------------------------------------------
        ID_vals = (rows * nx) + cols
        d8_vals = self.d8_grid.flat[ ID_vals ]
##        print 'ID_vals =', ID_vals
##        print 'd8_vals =', d8_vals
##        print 'h_vals  =', h
##        print ' '
        w       = np.where( d8_vals == h )
        n_kids  = w[0].size
        if (n_kids > 0):
            kid_IDs = ID_vals[ w ]
        else:
            kid_IDs = np.int32([])   # (empty list, size 0)

        return kid_IDs
    
    #   get_kid_IDs()
    #-------------------------------------------------------------------
    def add_neighbors_to_IDs(self):

        nx  = self.nx
        ny  = self.ny
        IDs = self.IDs
        rows, cols = divmod( IDs, nx )

        #-------------------------------
        # Get IDs of 8 neighbor pixels
        #-----------------------------------------------
        # Allow periodic BCs here via use of "%" (mod)
        #-----------------------------------------------          
        ID0 = ((rows-1) % ny) * nx + ((cols+1) % nx)
        ID1 = ( rows    % ny) * nx + ((cols+1) % nx)
        ID2 = ((rows+1) % ny) * nx + ((cols+1) % nx)
        ID3 = ((rows+1) % ny) * nx + ( cols    % nx)
        ID4 = ((rows+1) % ny) * nx + ((cols-1) % nx)
        ID5 = ( rows    % ny) * nx + ((cols-1) % nx)
        ID6 = ((rows-1) % ny) * nx + ((cols-1) % nx)
        ID7 = ((rows-1) % ny) * nx + ( cols    % nx)

        IDs = np.concatenate((IDs,ID0,ID1,ID2,ID3,ID4,ID5,ID6,ID7))
        IDs = np.unique( IDs )

        self.IDs = IDs
    
    #   add_neighbors_to_IDs()
    #-------------------------------------------------------------------
    def update_flow_grid(self, DEM=None, SILENT=True, REPORT=False):

        #--------------------------------------------------------
        # Note: This over-rides same function in d8_base.py
        #       mainly by using d8_codes instead of d8_grid.
        #--------------------------------------------------------
        
        #--------------------------------------------------------
        # This can be used to test whether "update_area_grid()"
        # is working, even if there is a problem with the
        # "update_flow_grid()" function.
        #--------------------------------------------------------
##        print '### TEST:  Reading existing flow grid instead'
##        print '###        of computing one.'
##        print '################################################'
##        self.read_flow_grid()
##        return
    
        #------------------------------------------------------
        # NOTES:  Direction codes are given by:  |64  128  1|
        #                                        |32   x   2|
        #                                        |16   8   4|

        #         A RiverTools resolve array is returned by:
        #            d8_base.get_resolve_array().
        #------------------------------------------------------
        if not(SILENT):    
            print('Updating D8 flow grid...')

        #------------------------------------------------------
        # Note: self.IDs is initialized to None, but may be
        #       set directly into self by the caller (e.g.
        #       erode_d8_local.py).
        #------------------------------------------------------
        # Note: It is not enough to update the flow codes at
        #       the input IDs because it is then possible for
        #       a neighbor ID to flow towards a cell that is
        #       supposed to flow to it (backflow & loops).
        #       Modified on 3/1/12 to use get_neighbors();
        #       stops error messages from update_area_grid().
        #------------------------------------------------------        
        if (self.IDs == None):
            self.IDs = np.ravel(self.ID_grid)
        else:
            #--------------------------------------------------
            # See erode_d8_local.update_d8_vars() for several
            # other attempts to addresss this (commented),
            # such as "grandparent IDs" that are in IDs.
            #--------------------------------------------------
            # Expand IDs to include all neighbors. (3/1/12)
            #--------------------------------------------------
            # Notice use of new INCLUDE_IDS option.
            #----------------------------------------
            nbr_hood_IDs = self.get_neighbors( self.IDs,
                                               INCLUDE_IDS=True )
            self.IDs = nbr_hood_IDs
            
        #-------------------------------------------------------
        # Save current D8 codes to detect whether they change.
        # Moved here on (9/23/11).
        # D8 codes for IDs are reset to 0 in PERTURB_TEST.
        #-------------------------------------------------------
        self.D8_CODES_CHANGED = False
        self.start_d8_codes = self.d8_grid.flat[ self.IDs ]
        
        #----------------------------------------
        # Assign directions where well-defined,
        # and special codes otherwise.  
        #----------------------------------------
        self.start_new_d8_codes(DEM, SILENT=SILENT, REPORT=REPORT)
        
        #------------------------------------
        # Break ties with "resolve" array ?
        #--------------------------------------------------
        # It seems we will always want BREAK_TIES = True.
        #--------------------------------------------------
        if (self.BREAK_TIES):
            self.break_flow_grid_ties(SILENT=SILENT)
        
        #-------------------
        # Link the flats ?
        #-------------------
        if (self.LINK_FLATS):
            self.link_flats(SILENT=SILENT)

        #--------------------------------------------
        # Have any flow codes changed ?  If not, 
        # we don't need to update A grid. (4/21/10)
        #--------------------------------------------
        w = np.where( self.d8_codes != self.start_d8_codes )
        if (w[0].size > 0):
            self.D8_CODES_CHANGED = True

        #---------------------------------------------
        # Change data type from 2-byte to 1-byte ??
        # No longer seems to be needed. (9/23/11)
        # Note that "valid_code_map" would also set
        # "special codes" to zero.
        #---------------------------------------------
##        final_codes = self.valid_code_map[ self.d8_codes ]
##        self.d8_codes = final_codes
##        self.d8_grid.flat[ self.IDs ] = final_codes

        
        if not(SILENT):
            uniq_codes = np.unique( self.d8_codes )
            print('   unique(codes) =', uniq_codes)
            
            # print '   min(codes), max(codes) =', \
            #       self.d8_codes.min(), self.d8_codes.max()
            # print 'd8_codes ='
            # print self.d8_codes

        #-------------------------------------------------
        # Compare saved flow grid to one just computed
        #-------------------------------------------------
        # Note that unless self.LR_PERIODIC = False and
        # self.TB_PERIODIC = False, the flow grids won't
        # agree on the edges.  This is because we don't
        # set them to zero when using periodic boundary
        # conditions.
        #-------------------------------------------------
        if (self.RT3_TEST):
            code_file = (self.out_directory +
                         self.site_prefix + '_flow.rtg')
            saved_flow_grid = rtg_files.read_grid(code_file, self.rti,
                                                  RTG_type='BYTE')
            w = np.where( saved_flow_grid != np.uint8(self.d8_grid) )
            if (w[0].size == 0):
                print('#### SUCCESS! Flow grids are identical to RT3.')
            else:
                print('##################################################')
                print(' WARNING: Flow grids differ from RT3:')
                print('          Number of pixels =', w[0].size)
                print('          Likely due to how flats are resolved.')
                print(' ')

    #   update_flow_grid()
    #-------------------------------------------------------------------    
    def start_new_d8_codes(self, DEM=None,
                           SILENT=True, REPORT=False):
        
        #--------------------------------------------------------------
        # Notes: If (z is not None), then update D8 codes only for the
        #        pixels with those IDs. (3/2/10)
        #        Should IDs and z be set directly into d8's state
        #        instead of being passed as args?
        #
        #        In caller, modified so that DEM array has
        #        type INTEGER when DEM has type BYTE.  Need a signed
        #        type to compute slopes correctly here.  For example,
        #        (100b - 200b) = 156b.

        #        Use of NumPy's ROLL induces periodic boundaries.
        #--------------------------------------------------------------
        if not(SILENT):
            print('   update_d8_codes(): Initializing grid...')
            ## print '   Starting new D8 codes...'
            ## print '   Initializing D8 codes in update_d8_codes()...'
        
        #--------------------------------------        
        # Compute and save rows,cols from IDs
        #--------------------------------------
        rc_tuple = divmod( self.IDs, self.nx )
        rows = rc_tuple[0]
        cols = rc_tuple[1]
        self.ID_rows = rows    # (reused in other methods)
        self.ID_cols = cols
        
        #-----------------------------
        # Define some local synonyms
        #-----------------------------
        nx   = self.nx
        ny   = self.ny
        #----------------------------------
        # For now, assume that all pixels
        # have the same dimensions.
        #----------------------------------
        dx = self.dx[0]
        dy = self.dy[0]
        dd = self.dd[0]

        #-------------------------------------------------------
        # DEM is read from file by initialize_computed_vars()
        # Depressions should have already been filled, perhaps
        # using the FILL_PITS_IN_Z0 option there.
        #-------------------------------------------------------
        if (DEM == None):
            DEM = self.DEM
            
        #--------------------------
        # Get the z-values at IDs
        #--------------------------------------------------
        # Note: z must be 1D because DEM[rows,cols] is 1D
        #       even if rows and cols are for all pixels.
        #--------------------------------------------------
        z = DEM.flat[ self.IDs ]
        
        #------------------------------
        # Slopes to 8 neighbor pixels
        #-----------------------------------------------
        # Allow periodic BCs here via use of "%" (mod)
        #-----------------------------------------------
##        s1 = (z - DEM[(rows-1) % ny, (cols+1) % nx]) / dd   # (upper-right)
##        s2 = (z - DEM[rows,          (cols+1) % nx]) / dx   # (right)
##        s3 = (z - DEM[(rows+1) % ny, (cols+1) % nx]) / dd   # (lower-right)
##        s4 = (z - DEM[(rows+1) % ny,  cols])         / dy   # (bottom)
##        s5 = (z - DEM[(rows+1) % ny, (cols-1) % nx]) / dd   # (lower-left)
##        s6 = (z - DEM[rows,          (cols-1) % nx]) / dx   # (left)
##        s7 = (z - DEM[(rows-1) % ny, (cols-1) % nx]) / dd   # (upper-left)
##        s8 = (z - DEM[(rows-1) % ny,  cols])         / dy   # (top)

        #------------------------------
        # Slopes to 8 neighbor pixels
        #------------------------------------------------
        # Allow periodic BCs here via use of "%" (mod)
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
        s1 = (z - DEM[rows_U,  cols_R]) / dd   # (upper-right)
        s2 = (z - DEM[rows,    cols_R]) / dx   # (right)
        s3 = (z - DEM[rows_D,  cols_R]) / dd   # (lower-right)
        s4 = (z - DEM[rows_D,  cols])   / dy   # (bottom)
        s5 = (z - DEM[rows_D,  cols_L]) / dd   # (lower-left)
        s6 = (z - DEM[rows,    cols_L]) / dx   # (left)
        s7 = (z - DEM[rows_U,  cols_L]) / dd   # (upper-left)
        s8 = (z - DEM[rows_U,  cols])   / dy   # (top)
        
        #------------------------------
        # Slopes to 8 neighbor pixels
        #-----------------------------------
        # This does not allow periodic BCs
        #-----------------------------------
        ## IDs = self.IDs
        ## s1 = (z - DEM.flat[IDs - nx + 1]) / dd   # (upper-right)
        ## s2 = (z - DEM.flat[IDs + 1])      / dx   # (right)
        ## s3 = (z - DEM.flat[IDs + nx + 1]) / dd   # (lower-right)
        ## s4 = (z - DEM.flat[IDs + nx])     / dy   # (bottom)
        ## s5 = (z - DEM.flat[IDs + nx - 1]) / dd   # (lower-left)
        ## s6 = (z - DEM.flat[IDs - 1])      / dx   # (left)
        ## s7 = (z - DEM.flat[IDs - nx - 1]) / dd   # (upper-left)
        ## s8 = (z - DEM.flat[IDs - nx])     / dy   # (top)
        
        #--------------------------
        # Find the steepest slope
        #--------------------------
##        max_slope = np.maximum(s1, s2)
##        max_slope = np.maximum(max_slope, s3)
##        max_slope = np.maximum(max_slope, s4)
##        max_slope = np.maximum(max_slope, s5)
##        max_slope = np.maximum(max_slope, s6)
##        max_slope = np.maximum(max_slope, s7)
##        max_slope = np.maximum(max_slope, s8)
        #--------------------------------------------
        # Find the steepest slope.
        # Faster way using 3rd, "out" argument ?
        # Works, but isn't any faster. (1/25/12)
        # At least not when used with small arrays.
        #--------------------------------------------
        max_slope = np.maximum( s1, s2 )
        np.maximum( max_slope, s3, max_slope )
        np.maximum( max_slope, s4, max_slope )
        np.maximum( max_slope, s5, max_slope )
        np.maximum( max_slope, s6, max_slope )
        np.maximum( max_slope, s7, max_slope )
        np.maximum( max_slope, s8, max_slope )
      
        #---------------------------------------------------------
        # WARNING: This method of flagging and later breaking
        #          flow direction ties will only work with
        #          power-of-2 flow codes (e.g. Jenson 84 or ARC)
        #---------------------------------------------------------
        # Note:  self.code_list has dtype of 'int16'
        #        max_slope has dtype of float32
        #---------------------------------------------------------
        # print 'max_slope.dtype      =', max_slope.dtype
        # print 'shape(max_slope)     =', shape(max_slope)
        # print 'self.code_list.dtype =', self.code_list.dtype
        #---------------------------------------------------------
        # Note that as a product of boolean and int16, d8_codes
        # has a data type of int16. (g has type int16.)
        #---------------------------------------------------------        
        g        = self.code_list
        d8_codes = (s1 == max_slope) * g[0] + \
                   (s2 == max_slope) * g[1] + \
                   (s3 == max_slope) * g[2] + \
                   (s4 == max_slope) * g[3] + \
                   (s5 == max_slope) * g[4] + \
                   (s6 == max_slope) * g[5] + \
                   (s7 == max_slope) * g[6] + \
                   (s8 == max_slope) * g[7]
        ## print 'd8_codes.dtype =', d8_codes.dtype

        #-----------------------------------------------------
        # (1/25/12) This method seems to be slightly faster
        # but it requires more function calls.
        #-----------------------------------------------------
##        g        = self.code_list
##        d8_codes = np.zeros( s1.size, dtype='int16' )
##        #-------------------------------------------------
##        v0 = np.where( s1 == max_slope )
##        if (v0[0].size != 0):
##            d8_codes[ v0 ] += g[0]
##        #--------------------------------
##        v1 = np.where( s2 == max_slope )
##        if (v1[0].size != 0):
##            d8_codes[ v1 ] += g[1]
##        #--------------------------------
##        v2 = np.where( s3 == max_slope )
##        if (v2[0].size != 0):
##            d8_codes[ v2 ] += g[2]
##        #--------------------------------
##        v3 = np.where( s4 == max_slope )
##        if (v3[0].size != 0):
##            d8_codes[ v3 ] += g[3]
##        #--------------------------------          
##        v4 = np.where( s5 == max_slope )
##        if (v4[0].size != 0):
##            d8_codes[ v4 ] += g[4]
##        #--------------------------------
##        v5 = np.where( s6 == max_slope )
##        if (v5[0].size != 0):
##            d8_codes[ v5 ] += g[5]
##        #--------------------------------
##        v6 = np.where( s7 == max_slope )
##        if (v6[0].size != 0):
##            d8_codes[ v6 ] += g[6]
##        #--------------------------------
##        v7 = np.where( s8 == max_slope )
##        if (v7[0].size != 0):
##            d8_codes[ v7 ] += g[7]

            
        #------------------------------------------------------
        # If we want to fill depressions "naturally" then
        # assign flats and pits a flow code of 0. (3/12/10)
        # If they don't have any outflow, they'll get filled.
        #------------------------------------------------------
        if not(self.LINK_FLATS):
            #---------------------------------------------
            # Set D8 code to zero at pits and store IDs.
            # (9/8/10) Includes "pits" on the edges.
            #---------------------------------------------
            pits = np.where( max_slope < 0 )
            n_pits = pits[0].size
            if (n_pits > 0):
                d8_codes[ pits ] = 0
                self.pit_IDs = self.IDs[ pits ]
            else:
                self.pit_IDs = np.array([])
                
            #----------------------------------------------
            # Set D8 code to zero at flats and store IDs.
            # (9/8/10) Includes "flats" on the edges.
            #----------------------------------------------
            flats = np.where( max_slope == 0 )
            n_flats = flats[0].size
            if (n_flats > 0):
                d8_codes[ flats ] = 0
                self.flat_IDs = self.IDs[ flats ]
            else:
                self.flat_IDs = np.array([])

            #-------------------------------------------
            # This is a bit faster if we don't need to
            # know where pits and flats are located.
            #-------------------------------------------
##            flats_and_pits = where(max_slope <= 0)
##            n_fp = flats_and_pits[0].size
##            if (n_fp != 0):
##                d8_codes[ flats_and_pits ] = 0
                
        else:     
            #------------------------------------------
            # Assign negative codes to flats and pits
            #------------------------------------------
            flats   = np.where(max_slope == 0)
            n_flats = flats[0].size
            if (n_flats != 0):
                d8_codes[ flats ] = -d8_codes[ flats ]
                ### print '########## d8_codes[ flats ] =', d8_codes[ flats ]
                
                #---------------------------------------------------
                # Above is OK, but this gives typecasting error ??
                #---------------------------------------------------
                ## d8_codes[ flats ] = (-1 * d8_codes[ flats ])
            
            #---------------------------------------------
            # There shouldn't be any of these left since
            # they were filled by fill_pits.fill_pits().
            #---------------------------------------------
            pits   = np.where(max_slope < 0)
            n_pits = pits[0].size
            if (n_pits != 0):
                d8_codes[ pits ] = -300
        
        #---------------------------------------------
        # Assign code of zero to NODATA & NaN pixels
        # Don't use NOT wrapped around FINITE. (254)
        #---------------------------------------------
        # Also assign code of zero to pixels
        # that are marked with RT closed-basin code?
        # Streamlines can end at either place.
        #----------------------------------------------
        w = np.where( np.logical_or((z <= self.nodata), (np.isfinite(z) != 1)) )
            #### or (z eq self.closed_basin_code)           
        n_bad = w[0].size
        if (n_bad != 0):
            d8_codes[ w ] = 0
        
        #-----------------------------------------------
        # Set left & right flow grid borders to zero ?
        #-----------------------------------------------
        if not(self.LR_PERIODIC):
            w0 = np.where( np.logical_or(cols == 0, cols == (nx-1)) )
            if (w0[0].size != 0):
                d8_codes[ w0 ] = 0
            
        #-----------------------------------------------
        # Set top & bottom flow grid borders to zero ?
        #-----------------------------------------------
        if not(self.TB_PERIODIC):
            w0 = np.where( np.logical_or(rows == 0, rows == (ny-1)) )
            if (w0[0].size != 0):
                d8_codes[ w0 ] = 0
            
        #----------------------------------------------------
        # Note: Need to copy new D8 codes into self.d8_grid
        #       here, but link_flats() will use them and
        #       make some changes later on.
        #----------------------------------------------------   
        self.d8_codes = d8_codes
        self.d8_grid.flat[ self.IDs ] = d8_codes
                
        if (REPORT):
            dmin   = d8_codes.min()
            dmax   = d8_codes.max()
            print('   --------------------------------------------')
            print('   Data type of flow grid at start =', d8_codes.dtype)
            print('   Number of flats         =', n_flats)
            print('   Number of 1-pixel pits  =', n_pits)
            
##            if not(self.LINK_FLATS):
##                print '   Number of flats & pits  =', n_fp
##            else:
##                print '   Number of flats         =', n_flats
##                print '   Number of 1-pixel pits  =', n_pits
                
            print('   Number of nodata/NaN    =', n_bad)
            print('   min(codes), max(codes)  =', dmin, dmax)
            print('   --------------------------------------------')

    #   start_new_d8_codes()
    #-------------------------------------------------------------------    
    def break_flow_grid_ties(self, SILENT=True):

        #----------------------------------------------------------
        # Notes: This routine resolves all non-flat ties, using
        #        a "tie-breaker array" called "resolve".  Note
        #        that "resolve" always returns one of the eight
        #        D8 flow codes.  Note that "resolve" maps valid
        #        D8 flow codes to themselves.
        #----------------------------------------------------------
       
        if not(SILENT):
            print('   update_d8_codes(): Breaking ties...')
            ## print '   Breaking ties in update_d8_codes():...'

        w = np.where(self.d8_codes > 0)
        if (w[0].size != 0):
            new_codes = self.resolve[ self.d8_codes[w] ]
            self.d8_codes[w] = new_codes
            #-------------------------------------------------
            # (9/23/11) Copy new codes directly into d8_grid
            # since need to know them when linking flats.
            #-------------------------------------------------
            self.d8_grid.flat[ self.IDs[w] ] = new_codes
        
    #   break_flow_grid_ties()
    #-------------------------------------------------------------------    
    def link_flats(self, SILENT=True):

        #--------------------------------------------------------
        # Notes: This procedure uses NumPy array operations to
        #        eliminate all loops and to increase speed.
        #        However, it is not equivalent to LinkUp2,
        #        because the sequential processing of lines
        #        there makes it possible to be linked to a
        #        just-resolved left neighbor.

        #       DIRECTION CODES:      SUBSCRIPTS:
        #       ----------------      -----------
        #          |64 128 1|          |6  7  0|    (i1)
        #          |32  x  2|          |5  x  1|    (i2)
        #          |16  8  4|          |4  3  2|    (i3)

        #       Note that resolve[code] always returns one
        #       of the component directions that made "code";
        #       it never "splits the difference".

        #       Need to make sure edges are not set to zero
        #       yet, but also need to exclude them from
        #       being counted as FLATS.

        #       Can't use the line:
        #           dirs[:,0]=0  &  dirs[:,NS-1]=0

        #       And can't use the line:
        #           FLATS=where((dirs[i2, 1:NS-1] < 0) AND $
        #                       (dirs[i2, 1:NS-1] != -300), NUM)

        #       Current solution is to use a NOTEDGE array.
        #--------------------------------------------------------
        if not(SILENT):
            print('   update_d8_codes(): Linking flats...')
            ## print '   Linking flats in update_d8_codes()...'

        #-----------------        
        # Local synonyms
        #-----------------
        nx = self.nx
        ny = self.ny
        g  = self.code_list   # (RT D8 flow codes)
        h  = self.code_opps   # (opposites of RT D8 flow codes)
        #** g = [1, 2, 4, 8, 16, 32, 64, 128]
        #** h = [16, 32, 64, 128, 1, 2, 4, 8]

        d8_codes = self.d8_codes
        ID_rows  = self.ID_rows
        ID_cols  = self.ID_cols
        ############################
        NOT_EDGE = np.logical_and(ID_rows != 0, ID_rows != (ny-1))
        NOT_EDGE = np.logical_and(NOT_EDGE, ID_cols != 0)
        NOT_EDGE = np.logical_and(NOT_EDGE, ID_cols != (nx-1))
        total_flats = np.int64(0)
        n_reps      = np.int64(0)
        
        while (True):
            n_reps += 1
            STILL_ACTIVE = False
            flats = np.where(np.logical_and(np.logical_and((d8_codes < 0), \
                                                  (d8_codes != -300)), \
                                                  (NOT_EDGE == True)))
            n_flats = flats[0].size
            if (n_flats == 0):
                FINISHED = True
                break

            ## print '>>>> n_flats =', n_flats  #################
            
            #-----------------------------------------------------
            # NB! When using a "D8 subset", flats are no longer
            #     indices into the full D8 flow grid, but into
            #     a smaller, 1D array (d8_codes vs. d8_grid).
            #-----------------------------------------------------
            ## flats = flats[0]  # (not needed. 2/22/12)
            flat_IDs = self.IDs[ flats ]
            rows, cols = divmod( flat_IDs, nx )

            #------------------------------------------------
            # This is noticeably faster because things like
            # (cols+1)%nx are computed once and reused 3x.
            # Results from cProfile show this method's time
            # dropping from 14.6 to 12.1 seconds when used
            # within erode_d8_local.py.
            #------------------------------------------------
            cols_R = (cols + 1) % nx
            cols_L = (cols - 1) % nx
            rows_U = (rows - 1) % ny
            rows_D = (rows + 1) % ny
            
            #---------------------------------------------
            # Flow codes of 8 neighbor pixels (periodic)
            # but only for the "unknown" pixels.
            #---------------------------------------------
            # Allow periodic BCs here via use of "%" (mod)
            # e.g. (-1 % 10) = 9, (10 % 10) = 0.
            # row numbers increase downward.
            #-----------------------------------------------
            d0 = self.d8_grid[rows_U,  cols_R]   # (upper-right)
            d1 = self.d8_grid[rows,    cols_R]   # (right)
            d2 = self.d8_grid[rows_D,  cols_R]   # (lower-right)
            d3 = self.d8_grid[rows_D,  cols]     # (bottom)
            d4 = self.d8_grid[rows_D,  cols_L]   # (lower-left)
            d5 = self.d8_grid[rows,    cols_L]   # (left)
            d6 = self.d8_grid[rows_U,  cols_L]   # (upper-left)
            d7 = self.d8_grid[rows_U,  cols]     # (top)
        
##            #----------------------------------
##            # Flow codes of 8 neighbor pixels
##            # (Doesn't work for edge pixels.)
##            #----------------------------------
##            d0 = self.d8_grid.flat[flats - nx + 1]   # (upper-right)
##            d1 = self.d8_grid.flat[flats + 1]        # (right) 
##            d2 = self.d8_grid.flat[flats + nx + 1]   # (lower-right)
##            d3 = self.d8_grid.flat[flats + nx]       # (bottom)
##            d4 = self.d8_grid.flat[flats + nx - 1]   # (lower-left)
##            d5 = self.d8_grid.flat[flats - 1]        # (left)
##            d6 = self.d8_grid.flat[flats - nx - 1]   # (upper-left)
##            d7 = self.d8_grid.flat[flats - nx]       # (top)
                
            #--------------------------------------------
            # A direction is VALID if:
            # (1) the neighbor in this direction has
            #     a well-defined (gt 0) direction, and
            #
            # (2) this neighbor's direction is not back
            #     towards the center pixel (BACKFLOW)
            #
            #     Later use "READY" idea?
            #     What is the fastest way?
            #--------------------------------------------
            VALID0 = np.logical_and( (d0 > 0), (d0 != h[0]) )
            VALID1 = np.logical_and( (d1 > 0), (d1 != h[1]) )
            VALID2 = np.logical_and( (d2 > 0), (d2 != h[2]) )
            VALID3 = np.logical_and( (d3 > 0), (d3 != h[3]) )
            VALID4 = np.logical_and( (d4 > 0), (d4 != h[4]) )
            VALID5 = np.logical_and( (d5 > 0), (d5 != h[5]) )
            VALID6 = np.logical_and( (d6 > 0), (d6 != h[6]) )
            VALID7 = np.logical_and( (d7 > 0), (d7 != h[7]) )
                
            #--------------------------------------------
            # A direction is allowed only if it was one
            # of the directions that was summed in
            # producing the value of code.
            #------------------------------------------------
            # (9/23/11) Does this only work if g[7] = 128 ?
            #------------------------------------------------
            codes = (-1 * d8_codes[flats])
            #-------------------------------
            ALLOWED7 = (codes >= g[7])
            codes    = (codes - (g[7] * ALLOWED7))
            ALLOWED6 = (codes >= g[6])
            codes    = (codes - (g[6] * ALLOWED6))
            ALLOWED5 = (codes >= g[5])
            codes    = (codes - (g[5] * ALLOWED5)) 
            ALLOWED4 = (codes >= g[4])
            codes    = (codes - (g[4] * ALLOWED4))
            ALLOWED3 = (codes >= g[3])
            codes    = (codes - (g[3] * ALLOWED3)) 
            ALLOWED2 = (codes >= g[2])
            codes    = (codes - (g[2] * ALLOWED2))  
            ALLOWED1 = (codes >= g[1])
            codes    = (codes - (g[1] * ALLOWED1))
            ALLOWED0 = (codes >= g[0])
            #-----------------------------------------------------
            # Next line is not needed and is not in the original
            # original RT3 code. All ALLOWEDs are defined now
            # and codes are no longer needed. (9/23/11)
            #-----------------------------------------------------
            # codes = (codes - (g[0] * ALLOWED0))
               
            ready_code = ((np.logical_and(VALID0, ALLOWED0) * g[0]) + \
                          (np.logical_and(VALID1, ALLOWED1) * g[1]) + \
                          (np.logical_and(VALID2, ALLOWED2) * g[2]) + \
                          (np.logical_and(VALID3, ALLOWED3) * g[3]) + \
                          (np.logical_and(VALID4, ALLOWED4) * g[4]) + \
                          (np.logical_and(VALID5, ALLOWED5) * g[5]) + \
                          (np.logical_and(VALID6, ALLOWED6) * g[6]) + \
                          (np.logical_and(VALID7, ALLOWED7) * g[7]))

            ###########################################################
            # This doesn't actually solve the problem, even though it
            # makes unit_test2() report SUCCESS, because the final
            # area grid in both calls to d8_local is then incomplete.
            # We need to compare final area grid with extension
            # A-after.rtg" to the area grid from D8-global, as done
            # in unit_test().
            ###########################################################
            # Need the next block to fix a subtle bug.  Without
            # it, pixels with IDs = 96093, 96094, 96501 and 96502
            # in unit_test2('KY_SUB','SQUARE') form a self-crossing
            # loop with d8 codes= 4, 32, 1, 32.  This is because
            # iterative flat-linking for an entire grid ensures that
            # there is a "downstream first" ordering.  But this
            # can be violated when flat-linking for a subset of the
            # grid that has flats just outside of the subset.
            ###########################################################            
            # Set ready_code to zero if the pixel that a cell is
            # supposed to flow to is also a flat.  This is meant
            # to resolve a strange "loop" problem that was found
            # using unit_test2('KY_SUB', 'SQUARE').  (2/23/12)
            ###########################################################           
##            FIXABLE        = where(ready_code > 0)  # (returns a tuple)
##            cur_IDs        = flat_IDs[ FIXABLE ]
##            new_d8_codes   = self.resolve[ ready_code[ FIXABLE ] ]
##            new_parent_IDs = cur_IDs + self.inc_map[ new_d8_codes ]
##            p_rows, p_cols = divmod( new_parent_IDs, nx )
##            p_cols_R = (p_cols + 1) % nx
##            p_cols_L = (p_cols - 1) % nx
##            p_rows_U = (p_rows - 1) % ny
##            p_rows_D = (p_rows + 1) % ny
##            #-----------------------------------------------------
##            zc = self.DEM[ p_rows,    p_cols ]     # (center)
##            z0 = self.DEM[ p_rows_U,  p_cols_R ]   # (upper-right)
##            z1 = self.DEM[ p_rows,    p_cols_R ]   # (right)
##            z2 = self.DEM[ p_rows_D,  p_cols_R ]   # (lower-right)
##            z3 = self.DEM[ p_rows_D,  p_cols ]     # (bottom)
##            z4 = self.DEM[ p_rows_D,  p_cols_L ]   # (lower-left)
##            z5 = self.DEM[ p_rows,    p_cols_L ]   # (left)
##            z6 = self.DEM[ p_rows_U,  p_cols_L ]   # (upper-left)
##            z7 = self.DEM[ p_rows_U,  p_cols ]     # (top)
##            #----------------------------------------------
##            # Find the lowest neighbors of new_parent_IDs
##            #----------------------------------------------
##            zmin = np.minimum(z0, z1)
##            np.minimum( zmin, z2, zmin )
##            np.minimum( zmin, z3, zmin )
##            np.minimum( zmin, z4, zmin )
##            np.minimum( zmin, z5, zmin )
##            np.minimum( zmin, z6, zmin )
##            np.minimum( zmin, z7, zmin )
##            #--------------------------------------------------------
##            # If lowest neighbor has same height, this pixel is (or
##            # was) a flat and it may not be safe to link to it.
##            #--------------------------------------------------------            
##            w_flat = where( zmin == zc )
##            if (w_flat[0].size != 0):
##                ready_code[ FIXABLE[0][w_flat] ] = 0  # (WORKS !)
##                ## print '#### Setting some cells as NOT READY.'
##                ## print ' '
##                #------------------------------------------------               
##                # These don't work as expected, as explained
##                # on pages 8-9 in my I2PY 0.2 User's Guide.
##                # First one works for printing but not setting
##                # values, which must be done as shown above.
##                #------------------------------------------------
##                # ready_code[ FIXABLE ][ w_flat ] = 0
##                # ready_code[ FIXABLE[ w_flat ] ] = 0
            ###########################################################
                

            #---------------------------------------------------------
            # Note: ready_code is an array of codes >= 0, many of
            #       which will be zero.  FIXABLE is an array with
            #       the same length as ready_code that contains
            #       True and False.  ready_code[FIXABLE] returns
            #       just the ready_codes that are fixable, but
            #       "forgets" where they go.
            #       FIXABLE is a 2-tuple, with 1st element equal to
            #       the array of 1D subscripts and 2nd element null.
            #---------------------------------------------------------
            FIXABLE = np.where(ready_code > 0)   # (returns a tuple)
            n_fixable = FIXABLE[0].size
            
            if (n_fixable != 0):
                #-----------------------------------
                # Assign flow codes to READY flats
                #-----------------------------------
                total_flats    += n_fixable
                fixable_flats   = flats[0][ FIXABLE ]  #####
                ## old_d8_codes    = d8_codes[ fixable_flats]  # (not used)
                new_d8_codes    = self.resolve[ ready_code[ FIXABLE ] ]
                d8_codes[ fixable_flats ] = new_d8_codes
                #---------------------------------------------------
                # Need to store new d8_codes directly into d8_grid
                # here because of how we get d0,...,d7! (9/23/11)
                #---------------------------------------------------
                self.d8_grid.flat[ flat_IDs[FIXABLE] ] = new_d8_codes
                STILL_ACTIVE = True

                ############################################################
                # When d8_loc.update() is called for all IDs, these two
                # pixels get flow codes of 32.  But when it is called for
                # a square of perturbed IDs, they get codes of 4 and 1.
                ############################################################
                ## print '####### n_fixable =', n_fixable
                ## print '####### size(d8_codes) =', d8_codes.size
##                fixed_IDs = flat_IDs[ FIXABLE ]
##                w1 = where( fixed_IDs == 96093 )
##                if (w1[0].size != 0):
##
##                    print '####### Just linked flat at 96093.'
##                    print '#######    D8[ 96093 ] =', self.d8_grid.flat[ 96093 ]
##                    print '#######    D8[ 96094 ] =', self.d8_grid.flat[ 96094 ]
####                    w3 = where( fixed_IDs == 96094 )
####                    if (w3[0].size != 0):
####                        print '#######    Flat 96094 was also just linked.'
##                w2 = where( fixed_IDs == 96501 )
##                if (w2[0].size != 0):
##                    print '####### Just linked flat at 96501.'
##                    print '#######    D8[ 96501 ] =', self.d8_grid.flat[ 96501 ]
##                    print '#######    D8[ 96502 ] =', self.d8_grid.flat[ 96502 ]
##                    print ' '
##                    print '#######    DEM[ 96093 ] =', self.DEM.flat[96093]
##                    print '#######    DEM[ 96094 ] =', self.DEM.flat[96094]
##                    print '#######    DEM[ 96501 ] =', self.DEM.flat[96501]
##                    print '#######    DEM[ 96502 ] =', self.DEM.flat[96502]
####                    w4 = where( fixed_IDs == 96502 )
####                    if (w4[0].size != 0):
####                        print '#######    Flat 96502 was also just linked.'
####                    print ' '
##                w3 = where( fixed_IDs == 96094 )
##                if (w3[0].size != 0):
##                    print '####### Just linked flat at 96094.'
##                    print '#######    D8[ 96094 ] =', self.d8_grid.flat[ 96094 ]
##                    print '#######    D8[ 96093 ] =', self.d8_grid.flat[ 96093 ]
##                w4 = where( fixed_IDs == 96502 )
##                if (w4[0].size != 0):
##                    print '####### Just linked flat at 96502.'
##                    print '#######    D8[ 96502 ] =', self.d8_grid.flat[ 96502 ]
##                    print '#######    D8[ 96501 ] =', self.d8_grid.flat[ 96501 ]
                ############################################################

      
            #-----------------------------------------------
            # We're finished if we "fixed" everything that
            # was "fixable", even if some flats remain.
            #-----------------------------------------------
            FINISHED = (n_fixable == n_flats)
            if (FINISHED) or not(STILL_ACTIVE):
                break

        #------------------------------------------------
        # Set unresolved code values to zero. (9/23/11)
        #------------------------------------------------
        if not(FINISHED):
            bad = np.where( d8_codes < 0 )
            d8_codes[ bad ] = 0
            self.d8_grid.flat[ self.IDs[ bad ] ] = 0
            if not(SILENT):
                #----------------------------------------------
                # Note: Confirmed that the 576 unlinked flats
                #       for KY_Sub are all near the edges.
                #----------------------------------------------
                print('##################################################')
                print(' WARNING:  Could not link all flats in d8_local.')
                print('           Number of flats remaining =', bad[0].size)
                print('           These are probably near edges.')
                print('##################################################')
                print(' ')
        else:
            if not(SILENT):
                print('   All flats linked. #####')
            
        #---------------------------
        # Finished with while loop
        #---------------------------
        self.total_flats = total_flats
        self.d8_codes = d8_codes    
        
        if not(SILENT):
            print('   Number of iterations =', n_reps, ' (in link_flats())')
        
    #   link_flats()
    #---------------------------------------------------------------------
    def update_parent_ID_grid(self, SILENT=True):

        #--------------------------------------------------------
        # Get a grid which for each grid cell contains the
        # calendar-style index (or ID) of the grid cell that
        # its D8 code says it flows to.
        #--------------------------------------------------------
        # Note: This version can handle periodic boundaries,
        #       as can occur in a landscape evolution model.
        #--------------------------------------------------------
        # Jenson 1984 flow codes:  1, 2, 4, 8, 16, 32, 64, 128
        # ARC/INFO flow codes:     128, 1, 2, 4, 8, 16, 32, 64
        #--------------------------------------------------------
        if not(self.D8_CODES_CHANGED):
            return

##        print '###### CALLING update_parent_ID_grid() #####'
##        print '###### LR_PERIODIC =', self.LR_PERIODIC
##        print '###### TB_PERIODIC =', self.TB_PERIODIC
##        print ' '
        
        if not(SILENT):
            print('Finding parent pixel IDs...')
        
        nx   = self.nx
        ny   = self.ny
        IDs  = self.IDs
        dirs = self.d8_codes
        g    = self.code_list
        
        pIDs = IDs + self.inc_map[ dirs ]
        #---------------------------------------
        # Pixels with invalid flow directions,
        # like edges, get assigned a pID of 0.
        # Note that self.inc_map[0] = 0.
        #---------------------------------------
        wbad = np.where(dirs <= 0)
        nw   = wbad[0].size
        if (nw != 0):    
            pIDs[ wbad ] = 0
            
        #----------------------------------------------
        # Are any of the IDs on the edge of the DEM ?
        #----------------------------------------------
        rows, cols = divmod(IDs, nx)
        T = np.where( rows == 0 )
        B = np.where( rows == ny-1 )
        L = np.where( cols == 0 )
        R = np.where( cols == nx-1 )

        #------------------------------------------------
        # Remap parent IDs for pixels on left and right
        # edges to support periodic boundary conditions
        #------------------------------------------------
        if (self.LR_PERIODIC):
            dR = dirs[R]
            dL = dirs[L]
            #------------------------------------------------------------------------------
            w  = np.where(np.logical_or(np.logical_or((dR == g[0]), (dR == g[1])), (dR == g[2])))
            if (w[0].size != 0):    
                pIDs[R[w]] -= nx  # (subtract nx)
            #------------------------------------------------------------------------------
            w = np.where(np.logical_or(np.logical_or((dL == g[4]), (dL == g[5])), (dL == g[6])))
            if (w[0].size != 0):    
                pIDs[L[w]] += nx
        else:
            if (R[0].size != 0):  pIDs[R] = 0
            if (L[0].size != 0):  pIDs[L] = 0

        #------------------------------------------------
        # Remap parent IDs for pixels on top and bottom
        # edges to support periodic boundary conditions
        #------------------------------------------------
        if (self.TB_PERIODIC):
            dT = dirs[T]
            dB = dirs[B]
            #------------------------------------------------------------------------------
            w = np.where(np.logical_or(np.logical_or((dT == g[0]), (dT == g[6])), (dT == g[7])))
            if (w[0].size != 0):    
                pIDs[T[w]] += self.rti.n_pixels   ## (DOUBLE CHECK THIS)
            #------------------------------------------------------------------------------
            w = np.where(np.logical_or(np.logical_or((dB == g[2]), (dB == g[3])), (dB == g[4])))
            if (w[0].size != 0):    
                pIDs[B[w]] -= self.rti.n_pixels  # (subtract n_pixels)
        else:
            if (T[0].size != 0):  pIDs[T] = 0
            if (B[0].size != 0):  pIDs[B] = 0

        #----------------------------------------
        # Save new parent IDs in parent_ID_grid
        #----------------------------------------
        self.parent_ID_grid.flat[ IDs ] = pIDs
            
    #   update_parent_ID_grid()
    #---------------------------------------------------------------------
    def update_parent_IDs(self):

        #------------------------------------------------------      
        # Notes: See Notes for update_parent_ID_grid() above.
        #------------------------------------------------------
        if not(self.D8_CODES_CHANGED):
            return
        
        #-----------------------------------------
        # Save IDs as a tuple of row indices and
        # column indices, "np.where" style
        #-----------------------------------------
        pIDs = self.parent_ID_grid.flat[ self.IDs ]
        self.parent_IDs = divmod(pIDs, self.nx)

        #------------------------------------------
        # Save IDs as a 1D array of long integers
        # "calendar style"
        #------------------------------------------
        # self.parent_IDs = np.ravel(pIDs)
        
    #   update_parent_IDs()
    #---------------------------------------------------------------------
    def update_flow_from_IDs(self):

        #-----------------------------------------------------------
        # Notes:  This function returns the 4-byte long-integer
        #         array IDs of pixels that flow in a particular
        #         direction.
        #-----------------------------------------------------------
        # Notes:  (3/2/10) Here, self.d8_codes refers to a 1D
        #         array of flow codes (at self.IDs) instead of a
        #         full D8 flow grid.
        #-----------------------------------------------------------
        # Jenson 1984 flow codes:  g = 1, 2, 4, 8, 16, 32, 64, 128
        # ARC/INFO flow codes:     g = 128, 1, 2, 4, 8, 16, 32, 64
        #-----------------------------------------------------------
        if not(self.D8_CODES_CHANGED):
            return
        
        g = self.code_list
        
        w1 = np.where(self.d8_codes == g[0])
        self.n1 = w1[0].size   # (northeast)
        
        w2 = np.where(self.d8_codes == g[1])
        self.n2 = w2[0].size   # (east)
        
        w3 = np.where(self.d8_codes == g[2])
        self.n3 = w3[0].size   # (southeast)
        
        w4 = np.where(self.d8_codes == g[3])
        self.n4 = w4[0].size   # (south)
        
        w5 = np.where(self.d8_codes == g[4])
        self.n5 = w5[0].size   # (southwest)
        
        w6 = np.where(self.d8_codes == g[5])
        self.n6 = w6[0].size   # (west)
        
        w7 = np.where(self.d8_codes == g[6])
        self.n7 = w7[0].size   # (northwest)
        
        w8 = np.where(self.d8_codes == g[7])
        self.n8 = w8[0].size   # (north)

        #------------------------------------------------------
        # Convert indices in d8_codes to indices in d8_grid ?
        #----------------------------------------------------------
        # __builtin__.divmod() returns indices like np.where()
        # but seems no faster than (pIDs / nx, pIDs % nx).
        #----------------------------------------------------------
        self.w1 = divmod( self.IDs[ w1 ], self.nx )
        self.w2 = divmod( self.IDs[ w2 ], self.nx )            
        self.w3 = divmod( self.IDs[ w3 ], self.nx )
        self.w4 = divmod( self.IDs[ w4 ], self.nx )
        self.w5 = divmod( self.IDs[ w5 ], self.nx )
        self.w6 = divmod( self.IDs[ w6 ], self.nx )
        self.w7 = divmod( self.IDs[ w7 ], self.nx )
        self.w8 = divmod( self.IDs[ w8 ], self.nx )
            
    #   update_flow_from_IDs()
    #---------------------------------------------------------------------
    def update_flow_to_IDs(self):

        #---------------------------------------------------- 
        # Note: update_parent_ID_grid() sets parent_ID_grid
        #       and is written to allow periodic BCs.
        #----------------------------------------------------------
        # __builtin__.divmod() returns indices like np.where()
        # but seems no faster than (pIDs / nx, pIDs % nx).
        #----------------------------------------------------------
        # Get IDs of "parent cells" that are downstream
        # of pixels that flow in a given direction.
        #------------------------------------------------
        if not(self.D8_CODES_CHANGED):
            return
        
        if (self.n1 != 0):    # northeast
            p1_IDs  = self.parent_ID_grid[self.w1]
            self.p1 = divmod(p1_IDs, self.nx)
        else:
            self.p1 = None
        #-------------------------------------------------
        if (self.n2 != 0):     # east
            p2_IDs  = self.parent_ID_grid[self.w2]
            self.p2 = divmod(p2_IDs, self.nx)
        else:
            self.p2 = None
        #-------------------------------------------------
        if (self.n3 != 0):     # southeast
            p3_IDs  = self.parent_ID_grid[self.w3]
            self.p3 = divmod(p3_IDs, self.nx)
        else:
            self.p3 = None
        #-------------------------------------------------
        if (self.n4 != 0):     # south 
            p4_IDs  = self.parent_ID_grid[self.w4]
            self.p4 = divmod(p4_IDs, self.nx)
        else:
            self.p4 = None
        #-------------------------------------------------
        if (self.n5 != 0):     # southwest
            p5_IDs  = self.parent_ID_grid[self.w5]
            self.p5 = divmod(p5_IDs, self.nx)
        else:   
            self.p5 = None
        #-------------------------------------------------
        if (self.n6 != 0):     # west 
            p6_IDs  = self.parent_ID_grid[self.w6]
            self.p6 = divmod(p6_IDs, self.nx)
        else:
            self.p6 = None
        #-------------------------------------------------
        if (self.n7 != 0):     # northwest  
            p7_IDs  = self.parent_ID_grid[self.w7]
            self.p7 = divmod(p7_IDs, self.nx)
        else:
            self.p7 = None
        #-------------------------------------------------
        if (self.n8 != 0):     # north 
            p8_IDs  = self.parent_ID_grid[self.w8]
            self.p8 = divmod(p8_IDs, self.nx)
        else:
            self.p8 = None
        #-------------------------------------------------
        ##    print 'p1.shape, p1.size       =', p1.shape, p1.size
        ##    print 'w1[0].shape, w1[0].size =', w1[0].shape, w1[0].size      
        
        #-------------------------------------
        # Some flow directions may not occur
        #-------------------------------------
        self.p1_OK = (self.p1 is not None)
        self.p2_OK = (self.p2 is not None)
        self.p3_OK = (self.p3 is not None)
        self.p4_OK = (self.p4 is not None)
        self.p5_OK = (self.p5 is not None)
        self.p6_OK = (self.p6 is not None)
        self.p7_OK = (self.p7 is not None)
        self.p8_OK = (self.p8 is not None)

    #   update_flow_to_IDs()
    #-------------------------------------------------------------------
    def update_noflow_IDs(self):

        #----------------------------------------------------------
        # Note: These are not used by erode_d8_local.py and this
        #       function is no longer called by d8_base.update().
        #       Also, this version hasn't been tested.  To use
        #       it, self.noflow_IDs should be initialized to
        #       self.edge_IDs in initialize() method. (1/25/12)
        #----------------------------------------------------------
        w0 = np.where( self.d8_codes <= 0 )
        n0 = w0[0].size
        if (n0 == 0): return
        
        new_IDs = self.IDs[ w0 ]
        nf_IDs  = np.concatenate( (self.noflow_IDs, new_IDs) )
        self.noflow_IDs = np.unique( nf_IDs )
        self.n0 = self.noflow_IDs.size

    #   update_noflow_IDs()    
    #-------------------------------------------------------------------    
    def update_flow_width_grid(self, DOUBLE=False, METHOD2=False,
                               SILENT=True, REPORT=False):

        #-------------------------------------------------------------
        # NOTES: This routine returns the flow widths for each
        #        pixel in the DEM.  The extra width of flowing to a
        #        diagonal pixel is taken into account, as well as the
        #        lat/lon-dependence for DEMs with fixed-angle pixels
        #        (Geographic lat/lon).

        #        METHOD2 version ensures that the sum of all flow
        #        widths around a pixel is equal to 2*(dx + dy), but
        #        is incorrect for case of a plane and others.

        #        Flow widths are zero where (flow grid eq 0).
        
        #        Flow widths are for entire pixel and are appropriate
        #        for overland or subsurface flow.
        #------------------------------------------------------------- 
        #        Is this only used by Seepage function now ???
        #-------------------------------------------------------------
        # NB!    np.where returns a 2-tuple, but with "empty"
        #        second part when applied to a 1D array.  This means
        #        that nw = w[0].size still works.
        #-------------------------------------------------------------
        if not(self.D8_CODES_CHANGED):
            return
        
        if not(SILENT):
            print('Updating flow width grid...')

        g  = self.code_list
        fc = self.d8_codes   # (local synonym)

        #----------------
        # Diagonal flow
        #----------------
        wd = np.where(np.logical_or(np.logical_or(np.logical_or((fc == g[0]),  (fc == g[2])), \
                                                    (fc == g[4])), (fc == g[6])))
        nwd = wd[0].size
        if (nwd != 0):
            IDs  = self.IDs[ wd ]
            rows = (IDs / self.nx)
            if not(METHOD2):
                self.dw.flat[IDs] = self.dd[rows]
            else:    
                self.dw.flat[IDs] = (self.dx[rows] + self.dy[rows]) / 4

        #---------------------
        # East and west flow
        #---------------------
        wh  = np.where(np.logical_or((fc == g[1]), (fc == g[5])))
        nwh = wh[0].size
        if (nwh != 0):
            IDs = self.IDs[ wh ]
            rows = (IDs / self.nx)
            self.dw.flat[IDs] = self.dy[ rows ]
            if (METHOD2):    
                self.dw.flat[IDs] = self.dw.flat[IDs] / 2

        #-----------------------
        # North and south flow
        #-----------------------
        wv  = np.where(np.logical_or((fc == g[3]), (fc == g[7])))
        nwv = wv[0].size
        if (nwv != 0):
            IDs  = self.IDs[ wv ]
            rows = (IDs / self.nx)
            self.dw.flat[IDs] = self.dx[ rows ]
            if (METHOD2):    
                self.dw.flat[IDs] = self.dw.flat[IDs] / 2

        #---------------------------
        # Undefined flow direction
        #---------------------------
        wb  = np.where(fc == 0)
        nwb = wb[0].size
        if (nwb != 0):
            #--------------------------------------
            # This prevents divide by zero errors
            #--------------------------------------
            IDs  = self.IDs[ wb ]
            rows = (IDs / self.nx)
            self.dw.flat[IDs] = self.dx[ rows ]

        #------------------
        # Optional report
        #------------------
        if (REPORT):
            dw_vec = self.dw.flat[ self.IDs ]
            dw_str = str(dw_vec.min()) + ', ' + str(dw_vec.max())
            print('    min(dw), max(dw) = ' + dw_str + ' [m]')

    #   update_flow_width_grid()
    #-------------------------------------------------------------------
    def update_flow_length_grid(self, DOUBLE=False,
                                SILENT=True, REPORT=False):

        #-------------------------------------------------------------
        # NOTES: This routine returns the flow lengths for each
        #        pixel in the DEM.  The extra length of flowing to a
        #        diagonal pixel is taken into account, as well as the
        #        latitude-dependence for DEMs with fixed-angle pixels
        #        (Geographic lat/lon).

        #        The only difference between this and the Flow_Widths
        #        function is that roles of dx and dy are switched.

        #        Flow lengths are set to dx[0] for the pixels where
        #       (flow grid eq 0), such as on the edges of the DEM.
        #-------------------------------------------------------------
        if not(self.D8_CODES_CHANGED):
            return
        
        if not(SILENT):
            print('Updating flow length grid...')

        g  = self.code_list
        fc = self.d8_codes  # (local synonym)
              
        #----------------
        # Diagonal flow
        #----------------
        wd = np.where(np.logical_or(np.logical_or(np.logical_or((fc == g[0]),  (fc == g[2])), \
                                                    (fc == g[4])), (fc == g[6])))
        nwd = wd[0].size
        if (nwd != 0):
            IDs  = self.IDs[ wd ]
            rows = (IDs / self.nx)
            self.ds.flat[IDs] = self.dd[ rows ]
        
        #---------------------
        # East and west flow
        #---------------------
        wh  = np.where(np.logical_or((fc == g[1]), (fc == g[5])))
        nwh = wh[0].size
        if (nwh != 0):
            IDs  = self.IDs[ wh ]
            rows = (IDs / self.nx)
            self.ds.flat[IDs] = self.dx[ rows ]
        
        #-----------------------
        # North and south flow
        #-----------------------
        wv  = np.where(np.logical_or((fc == g[3]), (fc == g[7])))
        nwv = wv[0].size
        if (nwv != 0):
            IDs  = self.IDs[ wv ]
            rows = (IDs / self.nx)
            self.ds.flat[IDs] = self.dy[ rows ]

        #---------------------------
        # Undefined flow direction
        #---------------------------
        wb  = np.where(fc == 0)
        nwb = wb[0].size
        if (nwb != 0):
            #--------------------------------------
            # This prevents divide by zero errors
            #--------------------------------------
            IDs  = self.IDs[ wb ]
            rows = (IDs / self.nx)
            self.ds.flat[IDs] = self.dx[ rows ]
            
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            ds_vec = self.ds.flat[ self.IDs ]
            ds_str = str(ds_vec.min()) + ', ' + str(ds_vec.max())
            print('    min(ds), max(ds) = ' + ds_str + ' [m]')

    #   update_flow_length_grid()
    #-------------------------------------------------------------------
    def update_area_grid(self, SILENT=True, REPORT=False):

        #-------------------------------------------------------------
        # Notes: Idea is to find the pixels whose area has
        #        not yet been assigned, and to recursively
        #        assign areas to those for which all children
        #        have been assigned areas.

        #        g = [1, 2, 4, 8, 16, 32, 64, 128]
        #        h = [16, 32, 64, 128, 1, 2, 4, 8]

        #        Notice that the first and last pixel of each
        #        line are assigned an area of zero.

        #        DIRECTION CODES:      SUBSCRIPTS:
        #        ----------------      -----------
        #           |64 128 1|          |6  7  0|    (i1)
        #           |32  x  2|          |5  x  1|    (i2)
        #           |16  8  4|          |4  3  2|    (i3)
        #-------------------------------------------------------------
        #         Wherever the D8 flow code changes, the area grid
        #         of that pixel's parent pixel must be updated, as
        #         well as all of the pixels downstream.
        #
        #         For area grid updates, use the recursive approach
        #         (starting with pixels with no D8 kids and working
        #         downstream). However, we only need to reinitialize
        #         the area grid values for pixels where area could
        #         have changed; any area grid values upstream of
        #         these pixels will still be valid.
        #-------------------------------------------------------------
        if not(self.D8_CODES_CHANGED):
            return
        
        if not(SILENT):
            print('Updating upstream area grid...')

        #-----------------        
        # Local synonyms
        #-----------------
        nx = self.nx
        ny = self.ny
        h  = self.code_opps
        max_reps = 2 * max(nx, ny)   # (used for both loops)
        ## max_reps = 4 * max(nx, ny)   # (used for both loops)
        
        #---------------------------------------
        # For testing:  Save Ap as an RTG file
        # with A=0 at IDs.
        #---------------------------------------
        if (self.PERTURB_TEST):
            self.A.flat[ self.IDs ] = 0
            Ap      = self.A.copy()
            Ap_file = (self.out_directory +
                       self.case_prefix + '_A-perturb.rtg')
            rtg_files.write_grid( Ap, Ap_file, self.rti )
       
        #-----------------------------------------------------
        # (2/22/12) FAILED unit_test2('square') and
        # unit_test2('strip') for KY_Sub, even though TREYNOR
        # and BEAVER passed all 3 cases of unit_test2().
        #-----------------------------------------------------
        IDs = self.IDs.copy()   # (11/25/12, .copy)
        ## IDs = self.IDs       # (11/20/11)
        #-----------------------------------------------------
        # This inequality should be strict because if
        # (IDs.size == self.rti.n_pixels), we don't need
        # to determine which downstream pixels are affected.
        #-----------------------------------------------------
        if (IDs.size < self.rti.n_pixels):
            
            kIDs   = IDs.copy()          # (11/20/11, .copy)
            n_reps = np.int32( 0 )
            
            #-------------------------------------------------
            # Find all pixels downstream of IDs and redefine
            # IDs to include all pixels downstream.
            #-------------------------------------------------
            while (True):
                pIDs = self.parent_ID_grid.flat[ kIDs ]
                pIDs = np.unique( pIDs )
                nk   = kIDs.size  # (2/29/12)
                np   = pIDs.size

                #---------------------------------------------------
                # NOTE: This block is costly as written and kept
                #       getting triggered before (2/29/12).  The
                #       problem was that update_parent_IDs() was
                #       only updating parent IDs for IDs, and it
                #       was still possible for some neighbor cells
                #       to flow toward IDs that flowed to them.
                #       See update_d8_vars() for more info.
                #---------------------------------------------------
                # Other ideas to avoid this costly loop:
                
                # w_backflow = np.in1d( pIDs, kIDs )  # (11/20/11)
                # if (w_backflow.size != 0):
                
                # if (nk == 2) and (np == 2):
                #     if (np.all( kIDs == pIDs )):
                #         break
                #
                # if (nk == np) and (np <= 5):
                #     if (np.all( kIDs == pIDs )):
                #         break    
                #---------------------------------------------------                
                if (n_reps > max_reps):
                    if (self.AREA_GRID_WARNING):
                        print('############################################')
                        print(' WARNING: There is a problem in the')
                        print('          update_area_grid() function.')
                        print(' kIDs =', kIDs)
                        print(' pIDs =', pIDs)
                        print(' d8[pIDs] =', self.d8_grid.flat[ pIDs ])
                        print(' ')
                    break  #######
                #------------------------------------------
                # Try something faster than concatenate ?
                #------------------------------------------
                IDs  = np.concatenate( (IDs, pIDs) )
                IDs  = np.unique( IDs )
                ## kIDs = pIDs         # (11/25/12)
                kIDs    = pIDs.copy()  # (11/20/11, .copy)
                n_reps += 1

                #############################################
                #  BUG FIX (2/16/12)
                #############################################
                ## if (self.PERTURB_TEST):
                ##     print 'pIDs =', pIDs 
                # if (np == 1): break  # (parent_ID[0] = 0)
                #############################################
                if (np == 1) and (pIDs[0] == 0):
                    break
                #############################################

        #---------------------------------------------------------
        # Reset A to 0 for IDs and all pixels downstream of IDs.
        # IDs with d8_code = 0 will not be considered "ready"
        # and will continue to have A=0 at the end.
        #---------------------------------------------------------
        self.A.flat[ IDs ] = 0

        #---------------------------------------
        # These may not be used now. (2/22/12)
        #---------------------------------------
        self.new_A_IDs = IDs
        ######################


        #---------------------------------------
        # For testing:  Save A0 as an RTG file
        # with A=0 at IDs.
        #---------------------------------------
        if (self.PERTURB_TEST):
            A0      = self.A.copy()
            A0_file = (self.out_directory +
                       self.case_prefix + '_A-perturb2.rtg')
            rtg_files.write_grid( A0, A0_file, self.rti )
        
        #--------------------------------------------
        # Convert units for da from m^2 to km^2 ??
        #--------------------------------------------
        # da was stored by self.read_grid_info()
        # Units for da are specified in '_d8.cfg"
        # file as either 'm^2' or 'km^2'
        #-----------------------------------------
        if ('km' in self.A_units.lower()):
            pixel_area = self.da / 1e6
        else:
            pixel_area = self.da

        #------------------
        # Initialize vars
        #------------------------------------------------
        # initialize_computed_vars() initializes self.A
        #------------------------------------------------
        n_reps = np.int32( 0 )  # (reset; used above)
        FINISHED = False
        
        while (True):
        
            STILL_ACTIVE = False
            n_reps += 1

            #--------------------------
            # For debugging (3/11/10)
            #--------------------------
            if (n_reps > max_reps):
                FINISHED = True
                mr_str = str(max_reps)
                print('----------------------------------------')
                print(' ERROR: Number of iterations > ' + mr_str)
                print('        Aborting update_area_grid().')
                print('----------------------------------------')
                break

            #---------------------------------
            # Find IDs that still have A = 0
            #-----------------------------------------
            # Use this if it turns out that we need
            # for A to be defined at pits and flats.
            #-----------------------------------------
##            A_vals   = self.A.flat[ IDs ]
##            w        = np.where( A_vals == 0 )
##            IDs_left = IDs[ w ]
##            n_IDs    = IDs_left.size

            #---------------------------------------------
            # Find IDs that have a valid D8 code & A = 0                
            #-----------------------------------------------------------
            # Note: Pits, flats and edges should have d8_code = 0.
            #       Those in IDs had their area value set to 0 above.
            #       They will not be included in w and not processed,
            #       so they should still have A=0 at the end.
            #       This will make Q and Qs also equal 0 at pits, etc.
            #-----------------------------------------------------------
            # However, it is possible to compute A at pits if all of
            # the D8 kids have known A-values. We are choosing not to.
            #-----------------------------------------------------------
            A_vals   = self.A.flat[ IDs ]
            d8_vals  = self.d8_grid.flat[ IDs ]
            w        = np.where(np.logical_and((A_vals == 0), (d8_vals != 0)))
            IDs_left = IDs[ w ]
            n_IDs    = IDs_left.size
            ## print 'n_reps, n_IDs =', n_reps, n_IDs

            if (n_IDs == 0):
                 FINISHED = True
                 break

            rows, cols = divmod( IDs_left, nx )

            #------------------------------------------------
            # This is noticeably faster because things like
            # (cols+1)%nx are computed once and reused 3x.
            # Results from cProfile show this method's time
            # dropping from 14.6 to 12.1 seconds when used
            # within erode_d8_local.py.
            #------------------------------------------------
            cols_R = (cols + 1) % nx
            cols_L = (cols - 1) % nx
            rows_U = (rows - 1) % ny
            rows_D = (rows + 1) % ny
            
            #-------------------------------------------------
            # Upstream areas of 8 neighbor pixels (periodic)
            # but only for the "unknown" pixels.
            #-------------------------------------------------
            a0 = self.A[rows_U,  cols_R]   # (upper-right)
            a1 = self.A[rows,    cols_R]   # (right)
            a2 = self.A[rows_D,  cols_R]   # (lower-right)
            a3 = self.A[rows_D,  cols]     # (bottom)
            a4 = self.A[rows_D,  cols_L]   # (lower-left)
            a5 = self.A[rows,    cols_L]   # (left)
            a6 = self.A[rows_U,  cols_L]   # (upper-left)
            a7 = self.A[rows_U,  cols]     # (top)
            
            #---------------------------------------------
            # Flow codes of 8 neighbor pixels (periodic)
            # but only for the "unknown" pixels.
            #---------------------------------------------
            d0 = self.d8_grid[rows_U,  cols_R]   # (upper-right)
            d1 = self.d8_grid[rows,    cols_R]   # (right)
            d2 = self.d8_grid[rows_D,  cols_R]   # (lower-right)
            d3 = self.d8_grid[rows_D,  cols]     # (bottom)
            d4 = self.d8_grid[rows_D,  cols_L]   # (lower-left)
            d5 = self.d8_grid[rows,    cols_L]   # (left)
            d6 = self.d8_grid[rows_U,  cols_L]   # (upper-left)
            d7 = self.d8_grid[rows_U,  cols]     # (top)

##            # BEFORE (2/22/12)
##            d0 = self.d8_grid[(rows-1) % ny, (cols+1) % nx]   # (upper-right)
##            d1 = self.d8_grid[rows,          (cols+1) % nx]   # (right)
##            d2 = self.d8_grid[(rows+1) % ny, (cols+1) % nx]   # (lower-right)
##            d3 = self.d8_grid[(rows+1) % ny,  cols]           # (bottom)
##            d4 = self.d8_grid[(rows+1) % ny, (cols-1) % nx]   # (lower-left)
##            d5 = self.d8_grid[rows,          (cols-1) % nx]   # (left)
##            d6 = self.d8_grid[(rows-1) % ny, (cols-1) % nx]   # (upper-left)
##            d7 = self.d8_grid[(rows-1) % ny,  cols]           # (top)

            #-----------------------------------------
            # Initialize all unknown pixels as READY
            #-----------------------------------------
            READY = np.ones([n_IDs], dtype='UInt8')
            
            #----------------------------------
            # A pixel is not READY if any of
            # it's children have unknown area
            #----------------------------------
            # Pixels w/ no children are READY
            # unless they are "edge" pixels
            #----------------------------------
            w0 = np.where(np.logical_and((d0 == h[0]), (a0 == 0)))
            if (w0[0].size != 0):
                READY[w0] = 0
            #---------------------------------------------------
            w1 = np.where(np.logical_and((d1 == h[1]), (a1 == 0)))
            if (w1[0].size != 0):    
                READY[w1] = 0
            #---------------------------------------------------
            w2 = np.where(np.logical_and((d2 == h[2]), (a2 == 0)))
            if (w2[0].size != 0):    
                READY[w2] = 0
            #---------------------------------------------------
            w3 = np.where(np.logical_and((d3 == h[3]), (a3 == 0)))
            if (w3[0].size != 0):    
                READY[w3] = 0
            #---------------------------------------------------
            w4 = np.where(np.logical_and((d4 == h[4]), (a4 == 0)))
            if (w4[0].size != 0):    
                READY[w4] = 0
            #---------------------------------------------------
            w5 = np.where(np.logical_and((d5 == h[5]), (a5 == 0)))
            if (w5[0].size != 0):    
                READY[w5] = 0
            #---------------------------------------------------
            w6 = np.where(np.logical_and((d6 == h[6]), (a6 == 0)))
            if (w6[0].size != 0):    
                READY[w6] = 0
            #---------------------------------------------------
            w7 = np.where(np.logical_and((d7 == h[7]), (a7 == 0)))
            if (w7[0].size != 0):    
                READY[w7] = 0
            
            #----------------------------------------
            # If a pixel is my child and I'm READY,
            # then add it's value to mine.
            #----------------------------------------
            WR      = np.where(READY)
            n_ready = WR[0].size
            WR = WR[0]  ################
            
            if (n_ready != 0):
                #---------------------------
                # This also assigns areas
                # to pixels w/ no children
                #---------------------------
                STILL_ACTIVE = True
                rows = rows[WR]
                cols = cols[WR]
                if (pixel_area.size == 1):
                    self.A[rows, cols] = pixel_area
                else:
                    self.A[rows, cols] = pixel_area[rows, cols]
                #---------------------------------------------                    
                w0 = np.where(d0[WR] == h[0])
                if (w0[0].size != 0):
                    self.A[rows[w0], cols[w0]] += a0[WR[w0]]
                #---------------------------------------------
                w1 = np.where(d1[WR] == h[1])
                if (w1[0].size != 0):
                    self.A[rows[w1], cols[w1]] += a1[WR[w1]]
                #---------------------------------------------
                w2 = np.where(d2[WR] == h[2])
                if (w2[0].size != 0):
                    self.A[rows[w2], cols[w2]] += a2[WR[w2]]
                #---------------------------------------------
                w3 = np.where(d3[WR] == h[3])
                if (w3[0].size != 0):
                    self.A[rows[w3], cols[w3]] += a3[WR[w3]]
                #---------------------------------------------
                w4 = np.where(d4[WR] == h[4])
                if (w4[0].size != 0):
                    self.A[rows[w4], cols[w4]] += a4[WR[w4]]
                #---------------------------------------------
                w5 = np.where(d5[WR] == h[5])
                if (w5[0].size != 0):
                    self.A[rows[w5], cols[w5]] += a5[WR[w5]]
                #---------------------------------------------
                w6 = np.where(d6[WR] == h[6])
                if (w6[0].size != 0):
                    self.A[rows[w6], cols[w6]] += a6[WR[w6]]
                #---------------------------------------------
                w7 = np.where(d7[WR] == h[7])
                if (w7[0].size != 0):
                    self.A[rows[w7], cols[w7]] += a7[WR[w7]]
##            else:
##                print '-------------------------------------------'
##                print ' WARNING: n_ready = 0, n_IDs =', n_IDs
##                print '          n_reps =', n_reps
##                print '          rows =', rows
##                print '          cols =', cols
##                print '          D8[rows, cols] =', self.d8_grid[rows,cols]
##                print '-------------------------------------------'

            #------------------------
            # Are we finished now ?
            #------------------------
            FINISHED = (n_ready == n_IDs)
            
            #------------------------
            # Are we finished now ?
            #------------------------
##            if (n_ready != n_IDs):
##            ### if (n_ready != n_IDs) and (n_ready > 0):
##                FINISHED = False

                  #-------------------------------------------
                  # LATER, may be able to reset A to 0 for
                  # all pixels downstream of IDs as we
                  # iterate to compute A.
##                #-------------------------------------------
##                # Get the IDs of the D8 "parent pixels"
##                #-------------------------------------------
##                # Notes: Must remove duplicates in pIDs !!
##                #-------------------------------------------
##                pIDs  = self.parent_ID_grid.flat[ IDs ]
##                pIDs  = np.unique( pIDs )
##
##                #------------------------------------------
##                # Set A-values of pIDs to 0 to flag
##                # them as needing to be computed/updated.
##                # And so they don't get used by others.
##                #------------------------------------------
##                # self.A.flat[ IDs ]  = 0
##                self.A.flat[ pIDs ] = 0


            
##                #--------------------------------------
##                # Get new IDs, including those that
##                # weren't ready yet in the last batch
##                #--------------------------------------
##                w  = where(logical_not(READY))
##                nw = w[0].size  ## (n_IDs - n_ready)
##                
####                if (nw != 0):
####                    ## print 'Number of IDs to "redo" =', nw
####                    IDs = IDs[ w ]
####                else:
####                    IDs = pIDs
##                    
##                if (nw != 0):
##                    ## print 'Number of IDs to "redo" =', nw
##                    again_IDs = IDs[ w ] 
##                    IDs = np.concatenate( (again_IDs, pIDs) )
##                    IDs = np.unique(IDs)
##                else:
##                    IDs = pIDs
##                    
##                #-------------------
##                # For testing only
##                #-------------------
##                if (IDs.size > self.rti.n_pixels):
##                    FINISHED = True
##                    print '-------------------------------------------'
##                    print ' ERROR: Number of IDs exceeds grid size.'
##                    print '        Aborting update_area_grid().'
##                    print '-------------------------------------------'
##                    print ' '
##                    sys.exit()
##                    # break
##            else:    
##                FINISHED = True
                
            if (FINISHED) or not(STILL_ACTIVE):
                break
        
        ## if not(FINISHED) and not(SILENT):
##        if not(FINISHED):
##            n_left = (n_IDs - n_ready)  ####### DOUBLE CHECK THIS ######
##            print 'Area grid undefined for', n_left, 'pixels.'


        ## print 'n_reps =', n_reps
        
        #-----------------------------------------------
        # Set left & right area grid borders to zero ?
        #-----------------------------------------------
        ID_rows, ID_cols = divmod( self.IDs, nx )
        if not(self.LR_PERIODIC):
##            w0 = np.where( logical_or(ID_cols == 0, ID_cols == (nx-1)) )
##            if (w0.size != 0):
##                self.A.flat[ self.IDs[ w0 ]] = 0
            self.A[:, 0]    = 0
            self.A[:, nx-1] = 0
                
        #-----------------------------------------------
        # Set top & bottom area grid borders to zero ?
        #-----------------------------------------------
        if not(self.TB_PERIODIC):
            w0 = np.where( np.logical_or(ID_rows == 0, ID_rows == (ny-1)) )
##            if (w0.size != 0):
##                self.A.flat[ self.IDs[ w0 ]] = 0
            self.A[0, :]    = 0
            self.A[ny-1, :] = 0          

        #----------------------------------
        # Save area grid as an RTG file ?
        #----------------------------------
        if (self.PERTURB_TEST):
            Af      = self.A.copy()
            Af_file = (self.out_directory +
                       self.case_prefix + '_A-after.rtg')
            rtg_files.write_grid( Af, Af_file, self.rti )
            
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            if ('km' in self.A_units.lower()):
                unit_str = ' [km^2]'
            else:
                unit_str = ' [m^2]'
            A_str = str(self.A.min()) + ', ' + str(self.A.max())
            print('min(A), max(A) = ' + A_str + unit_str)
            print('Number of iterations = ' + str(n_reps))

        #-------------------------------------------------
        # Compare saved area grid to one just computed
        #-------------------------------------------------
        # Note that unless self.LR_PERIODIC = False and
        # self.TB_PERIODIC = False, the area grids won't
        # agree on the edges.  This is because we don't
        # set them to zero when using periodic boundary
        # conditions.
        #-------------------------------------------------
        if (self.RT3_TEST):
            area_file = (self.out_directory +
                         self.site_prefix + '_area.rtg')
            saved_area_grid = rtg_files.read_grid(area_file, self.rti,
                                                  RTG_type='FLOAT')
            w = np.where( saved_area_grid != float32(self.A) )
            if (w[0].size == 0):
                print('#### SUCCESS! Area grids are identical to RT3.')
            else:
                diff = np.absolute(saved_area_grid - self.A)
                print('#################################################')
                print(' WARNING: Area grids differ from RT3:')
                print('          Number of pixels   =', w[0].size)
                print('          Maximum difference =', diff.max())
                print('          Likely due to how flats are resolved')
                print('          or incorrect A_units in CFG file.')
                print(' ')
                
    #   update_area_grid()
    #-------------------------------------------------------------------

