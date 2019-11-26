
## Copyright (c) 2001-2013, Scott D. Peckham
## January 2009  (converted from IDL)
## May, July, August, October 2009
## March 2010 (now uses d8_base.py, as does d8_local.py)
## Februrary 2012 (np.size(a) changed to a.size; faster ??)
##
## NB!  TF expects d8.codes vs. d8.d8_grid

#############################################################
## Note:  channel_base.py calls many of these but then
##        saves the results among its state variables.
#############################################################

## from numpy import *   ## ELIMINATE THIS

import numpy as np

import os        # (for os.chdir(), in unit_test().)
import os.path
import time

from topoflow.components import d8_base
from topoflow.utils      import rtg_files

#-------------------------------------------
# For use outside of the TopoFlow package.
#-------------------------------------------
# import d8_base
# import rtg_files

#---------------------------------------------------------------------
#
#   class d8_global    (inherits from d8_base.py)
#
#       get_component_name()
#       get_attribute()             # (10/27/11)
#       update_flow_grid()
#          start_new_d8_codes()
#          break_flow_grid_ties()
#          link_flats()
#----------------------------------
#       update_parent_ID_grid()     # (can handle periodic BCs) 
#       update_parent_IDs()         # (used by erosion_base.update_slope_grid())
#       update_non_parent_IDs()     # (not working or needed yet)
#       update_flow_from_IDs()
#       update_flow_to_IDs()
#       update_noflow_IDs()
#----------------------------------
#       update_flow_width_grid()
#       update_flow_length_grid()
#       update_area_grid()          # (added on 10/28/09)
#       update_slope_grid()         # (added on 2019-11-16)
#       update_aspect_grid()        # (added on 2019-11-25)
#
#-----------------------------------------------------------------------
class d8_component( d8_base.d8_component ):

    #-----------------------------------------------------------------
    # Note: Do not define an __init__() method here.  It will
    #       override things needed from CSDMS_base.__init__()
    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'D8_Global',
        'version':            '3.5',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #------------------------------------------------------
        'comp_name':          'D8Global',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'd8_global.cfg.in',
        'cfg_extension':      '_d8_global.cfg',
        'cmt_var_prefix':     '/D8Global/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Channels_Kinematic_Wave.xml',
        'dialog_title':       'D8 Global Parameters',
        'time_units':         'seconds' }

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_D8_Global'

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
    def update_flow_grid(self, DEM=None, SILENT=True, REPORT=False):

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
        #            self.get_resolve_array().
        #------------------------------------------------------
        if not(SILENT):    
            print('Updating D8 flow grid...')
            
        #----------------------------------------
        # Assign directions where well-defined,
        # and special codes otherwise
        #----------------------------------------
        self.start_new_d8_codes(DEM, SILENT=SILENT, REPORT=REPORT)
        
        #----------------------------------
        # Break ties with "resolve" array
        #----------------------------------
        if (self.BREAK_TIES):
            self.break_flow_grid_ties(SILENT=SILENT)
        
        #-------------------
        # Link the flats ?
        #-------------------
        if (self.LINK_FLATS): self.link_flats(SILENT=SILENT)
        
        #--------------------------------------------------
        # Change data type from 2-byte to 1-byte (NEW WAY)
        #--------------------------------------------------
        w = np.where( np.logical_or(self.d8_grid < 0, self.d8_grid > 128) )
        if (w[0].size > 0):
            self.d8_grid[w] = 0
        self.d8_grid = self.valid_code_map[ self.d8_grid ]
            
        ## if not(SILENT):
        if (self.DEBUG):
            print('   min(codes), max(codes) =', \
                  self.d8_grid.min(), self.d8_grid.max())
            print('d8_grid =')
            print(self.d8_grid)

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
            code_file = (self.in_directory +
                         self.site_prefix + '_flow.rtg')
            saved_d8_grid = rtg_files.read_grid(code_file, self.rti,
                                                RTG_type='BYTE')
            w = np.where( saved_d8_grid != uint8(self.d8_grid) )
            if (w[0].size == 0):
                print('##### SUCCESS! Flow grids are identical to RT3.')
            else:
                print('#################################################')
                print(' WARNING: Flow grids differ from RT3:')
                print('          Number of pixels   =', w[0].size)
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
  
        #-----------------------------
        # Define some local synonyms
        #-----------------------------
        nx = self.nx
        ny = self.ny
        #----------------------------------
        # For now, assume that all pixels
        # have the same dimensions.
        #----------------------------------
        dx = self.dx[0]
        dy = self.dy[0]
        dd = self.dd[0]
 
        #--------------------------------------------------------------
        # Note: DEM may be passed to this function but may also
        # have been read from DEM_file in initialize_computed_vars().
        #--------------------------------------------------------------
        if (DEM is None):
            DEM = self.DEM

        #------------------------------
        # Slopes to 8 neighbor pixels
        #--------------------------------------------------
        # Allow periodic BCs here via use of np.roll()
        #--------------------------------------------------
        s1 = (DEM - np.roll(np.roll(DEM, 1, axis=0), -1, axis=1))  / dd   # (upper-right)
        s2 = (DEM - np.roll(DEM, -1, axis=1))                         / dx   # (right)
        s3 = (DEM - np.roll(np.roll(DEM, -1, axis=0), -1, axis=1)) / dd   # (lower-right)
        s4 = (DEM - np.roll(DEM, -1, axis=0))                         / dy   # (bottom)
        s5 = (DEM - np.roll(np.roll(DEM, -1, axis=0), 1, axis=1))  / dd   # (lower-left)
        s6 = (DEM - np.roll(DEM, 1, axis=1))                          / dx   # (left)
        s7 = (DEM - np.roll(np.roll(DEM, 1, axis=0), 1, axis=1))   / dd   # (upper-left)
        s8 = (DEM - np.roll(DEM, 1, axis=0))                          / dy   # (top)

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
        #-----------------------------------------
        # Find the steepest slope.
        # Faster way using 3rd, "out" argument ?
        # (1/25/12) Worked in d8_local.py, but
        # wasn't any faster there.  TEST THIS.
        #-----------------------------------------
        max_slope = np.maximum(s1, s2)
        np.maximum(max_slope, s3, max_slope)
        np.maximum(max_slope, s4, max_slope)
        np.maximum(max_slope, s5, max_slope)
        np.maximum(max_slope, s6, max_slope)
        np.maximum(max_slope, s7, max_slope)
        np.maximum(max_slope, s8, max_slope)

        #---------------------------------------------------------
        # WARNING: This method of flagging and later breaking
        #          flow direction ties will only work with
        #          power-of-2 flow codes (e.g. Jenson 84 or ARC)
        #---------------------------------------------------------
        # Note that as a product of boolean and int16, d8_codes
        # has a data type of int16. (g has type int16.)
        #---------------------------------------------------------
        g       = self.code_list
        d8_grid = (s1 == max_slope) * g[0] + \
                  (s2 == max_slope) * g[1] + \
                  (s3 == max_slope) * g[2] + \
                  (s4 == max_slope) * g[3] + \
                  (s5 == max_slope) * g[4] + \
                  (s6 == max_slope) * g[5] + \
                  (s7 == max_slope) * g[6] + \
                  (s8 == max_slope) * g[7]

        #------------------------------------------------------
        # If we want to fill depressions "naturally" then
        # assign flats and pits a flow code of 0. (3/12/10)
        # If they don't have any outflow, they'll get filled.
        #------------------------------------------------------
        if not(self.LINK_FLATS):
            flats_and_pits = np.where(max_slope <= 0)
            n_fp = flats_and_pits[0].size
            if (n_fp != 0):
                d8_grid[ flats_and_pits ] = 0
        else:  
            #------------------------------------------
            # Assign negative codes to flats and pits
            #------------------------------------------
            flats   = np.where(max_slope == 0)
            n_flats = flats[0].size
            if (n_flats != 0):    
                d8_grid[ flats ] = (-1 * d8_grid[ flats ])
            
            #---------------------------------------------
            # There shouldn't be any of these left since
            # they were filled by fill_pits.fill_pits().
            #---------------------------------------------
            pits   = np.where(max_slope < 0)
            n_pits = pits[0].size
            if (n_pits != 0):
                d8_grid[ pits ] = -300
        
        #---------------------------------------------
        # Assign code of zero to NODATA & NaN pixels
        # Don't use NOT wrapped around FINITE. (254)
        #---------------------------------------------
        # Also assign code of zero to pixels
        # that are marked with RT closed-basin code?
        # Streamlines can end at either place.
        #----------------------------------------------
        w = np.where( np.logical_or((DEM <= self.nodata), (np.isfinite(DEM) != 1)) )            
        n_bad = w[0].size
        if (n_bad != 0):
            d8_grid[ w ] = 0
        
        #-----------------------------------------------
        # Set left & right flow grid borders to zero ?
        #-----------------------------------------------
        if not(self.LR_PERIODIC):
            d8_grid[:, 0]      = 0
            d8_grid[:, nx - 1] = 0
            
        #-----------------------------------------------
        # Set top & bottom flow grid borders to zero ?
        #-----------------------------------------------
        if not(self.TB_PERIODIC):
            d8_grid[0, :]      = 0
            d8_grid[ny - 1, :] = 0

        #-------------------------
        # Save D8 grid into self
        #-------------------------
        self.d8_grid = d8_grid
        
        if (REPORT):
            dmin = d8_grid.min()
            dmax = d8_grid.max()
            print('   --------------------------------------------')
            print('   Data type of flow grid at start = ' + str(d8_grid.dtype))
            if not(self.LINK_FLATS):
                print('   Number of flats & pits  = ' + str(n_fp) )
            else:
                print('   Number of flats         = ' + str(n_flats) )
                print('   Number of 1-pixel pits  = ' + str(n_pits) )
            print('   Number of nodata/NaN    = ' + str(n_bad) )
            print('   min(codes), max(codes)  = ' + str(dmin) + ', ' + str(dmax) )
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

        w = np.where(self.d8_grid > 0)
        if (w[0].size != 0):    
            self.d8_grid[w] = self.resolve[ self.d8_grid[w] ]            

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

        d8_grid     = self.d8_grid
        NOT_EDGE    = self.not_edge_grid 
        total_flats = np.int64(0)
        n_reps      = np.int64(0)
        
        while (True):
            n_reps += 1
            STILL_ACTIVE = False
            flats = np.where(np.logical_and(np.logical_and((d8_grid < 0), \
                                                  (d8_grid != -300)), \
                                                  (NOT_EDGE == 1)))

            n_flats = flats[0].size
            if (n_flats == 0):
                ## FINISHED = True  # (not needed)
                break
            
            rows = flats[0]
            cols = flats[1]

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
            d0 = d8_grid[rows_U,  cols_R]   # (upper-right)
            d1 = d8_grid[rows,    cols_R]   # (right)
            d2 = d8_grid[rows_D,  cols_R]   # (lower-right)
            d3 = d8_grid[rows_D,  cols]     # (bottom)
            d4 = d8_grid[rows_D,  cols_L]   # (lower-left)
            d5 = d8_grid[rows,    cols_L]   # (left)
            d6 = d8_grid[rows_U,  cols_L]   # (upper-left)
            d7 = d8_grid[rows_U,  cols]     # (top)
        
##            #----------------------------------
##            # Flow codes of 8 neighbor pixels
##            # (Doesn't work for edge pixels.)
##            #----------------------------------
##            d0 = d8_grid.flat[flats - nx + 1]   # (upper-right)
##            d1 = d8_grid.flat[flats + 1]        # (right) 
##            d2 = d8_grid.flat[flats + nx + 1]   # (lower-right)
##            d3 = d8_grid.flat[flats + nx]       # (bottom)
##            d4 = d8_grid.flat[flats + nx - 1]   # (lower-left)
##            d5 = d8_grid.flat[flats - 1]        # (left)
##            d6 = d8_grid.flat[flats - nx - 1]   # (upper-left)
##            d7 = d8_grid.flat[flats - nx]       # (top)
                
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
            #--------------------------------------------
            codes = (-1 * d8_grid[flats])
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
            
            FIXABLE   = np.where(ready_code > 0)  # (returns a tuple)
            n_fixable = FIXABLE[0].size
            
            if (n_fixable != 0):
                #-----------------------------------
                # Assign flow codes to READY flats
                #-----------------------------------
                total_flats += n_fixable
                fixed_tuple = ( rows[FIXABLE], cols[FIXABLE] )
                new_codes   = self.resolve[ ready_code[ FIXABLE ] ]
                d8_grid[ fixed_tuple] = new_codes
                STILL_ACTIVE = True

            #-----------------------------------------------
            # We're finished if we "fixed" everything that
            # was "fixable", even if some flats remain.
            #-----------------------------------------------
            FINISHED = (n_fixable == n_flats)
            if (FINISHED) or not(STILL_ACTIVE):
                break

        #---------------------------
        # Finished with while loop
        #---------------------------
        self.total_flats = total_flats
        #--------------------------------------
        # This not necessary
        # (9/23/11) But let's make 100% sure.
        #--------------------------------------
        self.d8_grid = d8_grid
        
        if not(SILENT):
            n_rep_str = str( n_reps )
            print('   Number of iterations = ' + n_rep_str + ' (in link_flats())')
        
    #   link_flats()
    #---------------------------------------------------------------------
    def update_parent_ID_grid(self, SILENT=True):

        #-----------------------------------------------------
        # Get a grid which for each grid cell contains the
        # calendar-style index (or ID) of the grid cell that
        # its D8 code says it flows to.
        #-----------------------------------------------------
        # Note: This version can handle periodic boundaries,
        #       as can occur in a landscape evolution model.
        #-----------------------------------------------------
        if not(SILENT):
            print('Finding parent pixel IDs...')
        
        nx = self.nx
        ny = self.ny
        self.parent_ID_grid = self.ID_grid + self.inc_map[self.d8_grid]

        #---------------------------------
        # Iterators for using 1D indices
        #---------------------------------
        dirs = self.d8_grid.flat
        pIDs = self.parent_ID_grid.flat

        #---------------------------------------
        # Get IDs for pixels on the four edges
        #---------------------------------------
        T = np.arange(nx, dtype='Int32')
        B = T + (nx * (ny - 1))
        L = nx * np.arange(ny, dtype='Int32')
        R = L + (nx - 1)

        #------------------------------------------------
        # Remap parent IDs for pixels on left and right
        # edges to support periodic boundary conditions
        #------------------------------------------------
        if (self.LR_PERIODIC):
            w  = np.where(np.logical_or(np.logical_or((dirs[R] == 1), (dirs[R] == 2)), (dirs[R] == 4)))
            if (w[0].size != 0):    
                pIDs[R[w]] -= nx  # (subtract nx)
            #-------------------------------------------------------------------------------------
            w = np.where(np.logical_or(np.logical_or((dirs[L] == 16), (dirs[L] == 32)), (dirs[L] == 64)))
            if (w[0].size != 0):    
                pIDs[L[w]] += nx
        else:
            pIDs[R] = 0
            pIDs[L] = 0

        #------------------------------------------------
        # Remap parent IDs for pixels on top and bottom
        # edges to support periodic boundary conditions
        #------------------------------------------------
        if (self.TB_PERIODIC):
            w = np.where(np.logical_or(np.logical_or((dirs[T] == 1), (dirs[T] == 64)), (dirs[T] == 128)))
            if (w[0].size != 0):    
                pIDs[T[w]] += self.rti.n_pixels   ## (DOUBLE CHECK THIS)
            #-------------------------------------------------------------------------------------
            w = np.where(np.logical_or(np.logical_or((dirs[B] == 4), (dirs[B] == 8)), (dirs[B] == 16)))
            if (w[0].size != 0):    
                pIDs[B[w]] -= self.rti.n_pixels  # (subtract n_pixels)
        else:
            pIDs[T] = 0
            pIDs[B] = 0
            
        #---------------------------------------
        # Pixels with invalid flow directions,
        # like edges, get assigned a pID of 0.
        #---------------------------------------
        wbad = np.where(self.d8_grid <= 0)
        nw   = wbad[0].size
        if (nw != 0):    
            self.parent_ID_grid[wbad] = 0

    #   update_parent_ID_grid()
    #---------------------------------------------------------------------
    def update_parent_IDs(self):

        #---------------------------------------------------------
        # Notes: This version cannot handle periodic boundaries,
        #        and requires that D8 flow codes be set to zero
        #        on the four edges.  This can be done at the end
        #        of the update_flow_grid() function.
        #---------------------------------------------------------
        # NB!  The use of 0's here is important.
        #      If iterating, pID[0]=0.
        #---------------------------------------------------------
       
        #-----------------------------------------------------
        # Get a grid which for each grid cell contains the
        # calendar-style index (or ID) of the grid cell that
        # its D8 code says it flows to.
        #-----------------------------------------------------
        pID_grid = self.parent_ID_grid

        ###################################################
        # Note that "divmod" is builtin and not in numpy.
        ###################################################
        
        #-----------------------------------------
        # Save IDs as a tuple of row indices and
        # column indices, "np.where" style
        #-----------------------------------------
        self.parent_IDs = divmod(pID_grid, self.nx)
        ## self.parent_IDs = (pID_grid / self.nx, pID_grid % self.nx)

        #------------------------------------------
        # Save IDs as a 1D array of long integers
        # "calendar style"
        #------------------------------------------
        # self.parent_IDs = np.ravel(pID_grid)
        
    #   update_parent_IDs()
    #-------------------------------------------------------------------
##    def update_non_parent_IDs(parent_IDs, flow_grid, rti):
##
##        #---------------------------
##        # Get flow grid dimensions
##        #---------------------------
##        nx = rti.ncols
##        ny = rti.nrows
##        
##        #---------------------------------------
##        # Return the IDs of non-parent pixels,
##        # such as ridges, but exclude pixels
##        # with flow code of 0, such as edges
##        # and nodata pixels.
##        #--------------------------------------- 
##        base = zeros([ny, nx], dtype='UInt8')
##        base[parent_IDs] = 1
##
##        wbad = no_flow_IDs(flow_grid, rti)
##        nbad = wbad[0].size
##        if (nbad > 0):
##            base[wbad] = 1   ##########  Should be 1 or 0 ??
##        
##        wnot = where(base == 0)
##        nnot = wnot[0].size
##        if (nnot != 0):    
##            non_parent_IDs = wnot
##        else:    
##            non_parent_IDs = -int32(1)
##
##        return non_parent_IDs
##
##    #   update_non_parent_IDs()
    #---------------------------------------------------------------------
    def update_flow_from_IDs(self):

        #----------------------------------------------------------
        # Notes:  This function returns the 4-byte long-integer
        #         array IDs of pixels that flow in a particular
        #         direction.  Jenson 1984 flow codes are assumed.
        #----------------------------------------------------------
        # Notes:  Later, rename w1 to w_NE, n1 to n_NE, then
        #         use self.code_list[0] vs. 1, etc..  This will
        #         then provide support for ARC flow codes, etc.
        #----------------------------------------------------------
        g = self.code_list
        
        self.w1 = np.where(self.d8_grid == g[0])
        self.n1 = self.w1[0].size   # (northeast)
        
        self.w2 = np.where(self.d8_grid == g[1])
        self.n2 = self.w2[0].size   # (east)
        
        self.w3 = np.where(self.d8_grid == g[2])
        self.n3 = self.w3[0].size   # (southeast)
        
        self.w4 = np.where(self.d8_grid == g[3])
        self.n4 = self.w4[0].size   # (south)
        
        self.w5 = np.where(self.d8_grid == g[4])
        self.n5 = self.w5[0].size   # (southwest)
        
        self.w6 = np.where(self.d8_grid == g[5])
        self.n6 = self.w6[0].size   # (west)
        
        self.w7 = np.where(self.d8_grid == g[6])
        self.n7 = self.w7[0].size   # (northwest)
        
        self.w8 = np.where(self.d8_grid == g[7])
        self.n8 = self.w8[0].size   # (north)

        ########### Same as noflow_IDs  ############
        ## self.w0 = np.where(self.d8_grid <= 0)
        ## self.n0 = self.w0[0].size   #(undefined)
            
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
        ##    print 'p1.shape,    p1.size    =', p1.shape, p1.size
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
        
        #--------------------------------------------------
        # 1/19/07.  Need to set d and u to zero at any ID
        # where flow terminates.  This includes pixels on
        # the edges, those with unresolved flow direction
        # and those where elevation is nodata or NaN.
        # A RiverTools flow grid will have a flow code of
        # zero at all of these places.
        #--------------------------------------------------
        noflow_IDs = np.where(self.d8_grid <= 0)
        num_IDs    = noflow_IDs[0].size
        
        if (num_IDs != 0):
            self.noflow_IDs = noflow_IDs
        else:
            #----------------------------
            # Return IDs of edge pixels
            #----------------------------
            ## self.get_edge_IDs()  # (called by initialize())
            self.noflow_IDs = self.edge_IDs()

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
        if not(SILENT):
            print('Updating flow width grid...')

        fg = self.d8_grid   # (local synonym)

        #----------------
        # Diagonal flow
        #----------------
        wd = np.where(np.logical_or(np.logical_or(np.logical_or((fg == 1), (fg == 4)), \
                                                    (fg == 16)), (fg == 64)))
        nwd = wd[0].size
        if (nwd != 0):
            rows = wd[0]
            if not(METHOD2):
                self.dw[wd] = self.dd[rows]
            else:    
                self.dw[wd] = (self.dx[rows] + self.dy[rows]) / 4

        #---------------------
        # East and west flow
        #---------------------
        wh  = np.where(np.logical_or((fg == 2), (fg == 32)))
        nwh = wh[0].size
        if (nwh != 0):
            self.dw[wh] = self.dy[ wh[0] ] # (wh[0] = rows)
            if (METHOD2):    
                self.dw[wh] = self.dw[wh] / 2

        #-----------------------
        # North and south flow
        #-----------------------
        wv  = np.where(np.logical_or((fg == 8), (fg == 128)))
        nwv = wv[0].size
        if (nwv != 0):
            self.dw[wv] = self.dx[ wv[0] ]  # (wv[0] = rows)
            if (METHOD2):    
                self.dw[wv] = self.dw[wv] / 2

        #---------------------------
        # Undefined flow direction
        #---------------------------
        wb  = np.where(fg == 0)
        nwb = wb[0].size
        if (nwb != 0):
            #--------------------------------------
            # This prevents divide by zero errors
            #--------------------------------------
            self.dw[wb] = self.dx[ wb[0] ] # (wb[0] = rows)

        #------------------
        # Optional report
        #------------------
        if (REPORT):
            dw_str = str(self.dw.min()) + ', ' + str(self.dw.max())
            print('    min(dw), max(dw) = ' + dw_str + ' [m]')
            
##            print '    min(dw) = ' + str(self.dw.min()) + '  [m]'
##            print '    max(dw) = ' + str(self.dw.max()) + '  [m]'

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
        if not(SILENT):
            print('Updating flow length grid...')

        fg = self.d8_grid  # (local synonym)
              
        #----------------
        # Diagonal flow
        #----------------
        wd = np.where(np.logical_or(np.logical_or(np.logical_or((fg == 1), (fg == 4)), \
                                                    (fg == 16)), (fg == 64)))
        nwd = wd[0].size
        if (nwd != 0):
            self.ds[wd] = self.dd[ wd[0] ] # (wd[0] = rows)
        
        #---------------------
        # East and west flow
        #---------------------
        wh  = np.where(np.logical_or((fg == 2), (fg == 32)))
        nwh = wh[0].size
        if (nwh != 0):    
            self.ds[wh] = self.dx[ wh[0] ] # (wh[0] = rows)
        
        #-----------------------
        # North and south flow
        #-----------------------
        wv  = np.where(np.logical_or((fg == 8), (fg == 128)))
        nwv = wv[0].size
        if (nwv != 0):    
            self.ds[wv] = self.dy[ wv[0] ] # (wv[0] = rows)

        #---------------------------
        # Undefined flow direction
        #---------------------------
        wb  = np.where(fg == 0)
        nwb = wb[0].size
        if (nwb != 0):
            #--------------------------------------
            # This prevents divide by zero errors
            #--------------------------------------
            self.ds[wb] = self.dx[ wb[0] ] # (wb[0] = rows)
            
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            ds_str = str(self.ds.min()) + ', ' + str(self.ds.max())
            print('    min(ds), max(ds) = ' + ds_str + ' [m]')
            
##            print '    min(ds) = ' + str(self.ds.min()) + '  [m]'
##            print '    max(ds) = ' + str(self.ds.max()) + '  [m]'

    #   update_flow_length_grid()
    #-------------------------------------------------------------------
    def update_area_grid(self, SILENT=True, REPORT=False):

        #------------------------------------------------------
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
        #------------------------------------------------------
        if not(SILENT):    
            print('Updating upstream area grid...')
        
        #------------------
        # Initialize vars
        #------------------------------------------------
        # initialize_computed_vars() initializes self.A
        # (2019-10-09) Update A grid in-place.
        #------------------------------------------------
        ## self.A = np.minimum(self.A, 0)   # (reset to all zeros)
        ## np.minimum(self.A, 0, self.A)    # (reset to all zeros, in-place)
        self.A[:] = np.minimum(self.A, 0)   # (reset to all zeros, in-place)
        n_reps = np.int32(0)

        #-----------------      
        # Local synonyms
        #-----------------
        nx = self.nx
        ny = self.ny
        h  = self.code_opps
        # g  = self.code_list  # (not used here)

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
    
        while (True):
        
            STILL_ACTIVE = False
            n_reps += 1
            
            #-------------------------------------------------------
            # Added test for (flow_grid != 0), so that edge
            # and nodata pixels won't be considered "READY".
            #-------------------------------------------------------
            # This led to an infinite loop bug because the edge
            # pixels were "ready", then PASS_CHANGE=1b, then their
            # areas were changed back to zero in caller which made
            # them "ready" again.
            #-------------------------------------------------------
            unknown   = np.where(np.logical_and((self.A == 0), (self.d8_grid != 0)))
            rows      = unknown[0]
            cols      = unknown[1]
            n_unknown = rows.size

            #######################################
            # Similar change made in link_flats()
            #######################################
            # if (n_unknown == 0):
            #     # UNFINISHED = False
            #     break
            #######################################
            
            if (n_unknown != 0):    

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

                #--------------------------------------
                # Upstream areas of 8 neighbor pixels
                # for all pixels in the grid
                #--------------------------------------
##                a0 = np.roll(np.roll(self.A, 1, axis=0), -1, axis=1)  # (upper-right)
##                a1 = np.roll(self.A, -1, axis=1)                         # (right)
##                a2 = np.roll(np.roll(self.A, -1, axis=0), -1, axis=1) # (lower-right)
##                a3 = np.roll(self.A, -1, axis=0)                         # (bottom)
##                a4 = np.roll(np.roll(self.A, -1, axis=0), 1, axis=1)  # (lower-left)
##                a5 = np.roll(self.A, 1, axis=1)                          # (left)
##                a6 = np.roll(np.roll(self.A, 1, axis=0), 1, axis=1)   # (upper-left)
##                a7 = np.roll(self.A, 1, axis=0)                          # (top)

                #-----------------------------------------
                # Initialize all unknown pixels as READY
                #-----------------------------------------
                READY = np.ones([n_unknown], dtype='UInt8')
                
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
                
                #------------------------
                # Are we finished now ?
                #------------------------
                if (n_unknown != 0) and (n_ready != n_unknown):   
                    UNFINISHED = True
                else:    
                    UNFINISHED = False
                
            else:    
                UNFINISHED = False

            if not(UNFINISHED) or not(STILL_ACTIVE):
                break
        
        if (UNFINISHED):    
            print('Upstream area not defined for all pixels.')

        #-----------------------------
        # Save area grid for testing
        #-----------------------------
        SAVE_TEST = False
        if (SAVE_TEST):
            file_unit = open('00_AREA_GRID_TEST.rtg', 'wb')
            if (self.rti.SWAP_ENDIAN):
                grid = self.A.copy()
                grid.byteswap(True)
                grid.tofile(file_unit)
            else:
                self.A.tofile( file_unit )
            file_unit.close()
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            if ('km' in self.A_units.lower()):
                unit_str = ' [km^2]'
            else:
                unit_str = ' [m^2]'
            A_str = str(self.A.min()) + ', ' + str(self.A.max())
            print('    min(A), max(A) = ' + A_str + unit_str)
            print('    Number of iterations = ' + str(n_reps))

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
            area_file = (self.in_directory +
                         self.site_prefix + '_area.rtg')
            saved_area_grid = rtg_files.read_grid(area_file, self.rti,
                                                  RTG_type='FLOAT')
            w = np.where( saved_area_grid != np.float32(self.A) )
            if (w[0].size == 0):
                print('##### SUCCESS! Area grids are identical to RT3.')
            else:
                diff = np.absolute(saved_area_grid - self.A)
                print('#################################################')
                print(' WARNING: Area grids differ from RT3:')
                print('          Number of pixels   =', w[0].size)
                print('          Maximum difference =', diff.max())
                print('          Likely due to how flats are resolved.')
                print(' ')
                
    #   update_area_grid()
    #-------------------------------------------------------------------
    def update_slope_grid(self, SILENT=True, REPORT=False):

        pIDs = self.parent_IDs
        self.S = (self.DEM - self.DEM[ pIDs ]) / self.ds

        #----------------------------------------------------    
        # In d8_base.py, noflow_IDs calc is commented out
        # in update().  noflow_IDs includes edge_IDs.
        #----------------------------------------------------
        ## self.S[ self.edge_IDs ] = 0.0
        self.update_noflow_IDs()
        self.S[ self.noflow_IDs ] = 0.0
        ## self.S[ self.S < 0 ] = 0.0

    #   update_slope_grid()
    #-------------------------------------------------------------------
    def update_aspect_grid(self, SILENT=True, REPORT=False):

        #------------------------------------------------------
        # Note: This is a D8 "aspect", really the flow angle,
        #       in radians, counterclockwise from due east.
        #       Map Jenson 1984 D8 codes to angles.
        #------------------------------------------------------
        a = np.linspace( 0.0, 2 * np.pi, 9 )
  
        amap = {1:a[1],   2:a[0],  4:a[7],   8:a[6],
                16:a[5], 32:a[4], 64:a[3], 128:a[2],
                0:a[0]}
                
        d8_grid = self.d8_grid
        aspect = np.zeros( self.d8_grid.shape, dtype='float32')
        aspect[ d8_grid == 0 ]   = 0.0   # (undefined)
        aspect[ d8_grid == 1 ]   = a[1]
        aspect[ d8_grid == 2 ]   = a[0]
        aspect[ d8_grid == 4 ]   = a[7]
        aspect[ d8_grid == 8 ]   = a[6]
        aspect[ d8_grid == 16 ]  = a[5]
        aspect[ d8_grid == 32 ]  = a[4]
        aspect[ d8_grid == 64 ]  = a[3]
        aspect[ d8_grid == 128 ] = a[2]        
        self.aspect = aspect

    #   update_aspect_grid()
    #-------------------------------------------------------------------
    
