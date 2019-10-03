
## Copyright (c) 2001-2010, Scott D. Peckham
## January 2009  (converted from IDL)
## May, July, August 2009
## May 2010 (changes to unit_test() and initialize() )

## NB!  TF expects d8.codes vs. d8.flow_grid

#############################################################
## Note:  channel_base.py calls many of these but then
##        saves the results among its state variables.
#############################################################
## Note:  It might be more clear to refer to flow_grids
##        and flow_codes as d8_grids and d8_codes, etc.
##        (e.g. d8_code_list, d8_code_map, d8_width_grid)
#############################################################

from numpy import *
import numpy

import os, os.path

from . import BMI_base
from . import pixels
from . import rtg_files
from . import tf_utils

from .tf_utils import TF_Print

## from model_output import *

#---------------------------------------------------------------------
#
#   unit_test()
#
#   class d8_base
#     
#       get_attribute()
#       get_status()
#       initialize()
#       read_flow_grid()
#----------------------------------
#       get_flow_code_list()
#       get_flow_code_list_opps()
#       get_parent_inc_map()
#       get_parent_ID_grid()
#       get_parent_IDs()         # (needed for gradients)
#       get_parent_IDs2()        # (not working or needed yet)
#       get_non_parent_IDs()     # (not working or needed yet)
#       get_flow_from_IDs()
#       get_flow_to_IDs()
#       get_edge_IDs()
#       get_noflow_IDs() 
#       get_flow_width_grid()
#       get_flow_length_grid()

#-----------------------------------------------------------------------
def unit_test():

    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory = tf_utils.TF_Test_Directory()
    os.chdir( cfg_directory )
    
    d8 = d8_component()

    cfg_prefix     = 'Case5'
    d8.site_prefix = 'Treynor'
    
    d8.initialize( cfg_prefix=cfg_prefix, mode="driver" )

    print(' ')
    print('size(flow_grid) =', size(d8.flow_grid))
    print('nx              =', d8.nx)
    print('ny              =', d8.ny)
    
#   unit_test()
#---------------------------------------------------------------------
class d8_component(BMI_base.BMI_component):

    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        #----------------
        # New. 10/26/11
        #----------------
        map = {'comp_name':          'TF_D8_component',
               'version':            '3.1',
               'model_name':         'd8_component class',
               'model_family':       'TopoFlow',
               'cfg_template_file':  'None',
               'cfg_extension':      '_d8.cfg',
               'cmt_var_prefix':     'None',
               'gui_xml_file':       'None',
               'dialog_title':       'None',
               'time_step_type':     'fixed',
               'time_units':         'seconds',
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
    def get_status(self):

        #-----------------------------------------------------
        # Notes: Return component status as a string.  The
        #        possible return values are from OpenMI 2.0:
        #
        #           created, initializing, initialized,
        #           updating, updated, finalizing, finalized,
        #           failed (could add "stopped").
        #-----------------------------------------------------
        return self.status

    #   get_status()
    #-----------------------------------------------------------------
    ## def initialize(self, cfg_prefix='Case1', mode="nondriver",
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=True, REPORT=False):

        #--------------------------------------------------------
        # Note:  This function calls functions that compute a
        #        variety of D8 variables and saves them in its
        #        state.  This "d8 state" can be embedded within
        #        another class such as the "channels_base" or
        #        "groundwater_base" class.
        #--------------------------------------------------------
        if not(SILENT):
            print(' ')
            print('D8 component: Initializing...')
        
        self.status = 'initializing'  # (OpenMI 2.0 convention)
        self.mode   = mode

        #------------------------------------------------------
        # Note: A run_model() call or a driver's initialize()
        #       call calling initialize_config_vars() will set
        #       CWD to the location of CFG files.
        #------------------------------------------------------
        # Note: If directories and prefixes are not set in
        #       initialize_config_vars(), then they will
        #       default to CWD and cfg_prefix.
        #------------------------------------------------------
        if (cfg_file == None):
            cfg_extension = self.get_attribute( 'cfg_extension' )
            filename      = self.site_prefix + cfg_extension
            cfg_file      = self.in_directory + filename
            ## self.cfg_file   = os.path.join( os.getcwd(), filename )
            ## self.cfg_prefix = self.site_prefix
        self.cfg_file = cfg_file

        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        ## self.set_constants()
        self.initialize_config_vars() 
        
        ## print '##### tf_d8_base.initialize() is calling read_grid_info'
        self.read_grid_info()
        
        self.read_flow_grid()
        
        self.get_flow_code_list()
        self.get_flow_code_list_opps()
        
        self.get_parent_inc_map()
        self.get_parent_ID_grid()

        self.get_parent_IDs()  # (needed for gradients)
        
        self.get_flow_from_IDs()
        self.get_flow_to_IDs()
        self.get_edge_IDs()
        self.get_noflow_IDs()
        
        self.get_flow_width_grid()
        self.get_flow_length_grid()

        self.status = 'initialized'
        
    #   initialize()
    #-------------------------------------------------------------------
    def read_flow_grid(self):

        #----------------------------------------------------
        # Read a grid of D8 flow codes, same size as DEM.
        #----------------------------------------------------
        TF_Print('Reading D8 flow grid...')
        code_file = (self.in_directory +
                     self.site_prefix + '_flow.rtg')
        self.flow_grid = rtg_files.read_grid(code_file, self.rti,
                                             RTG_type='BYTE')

    #   read_flow_grid()
    #---------------------------------------------------------------------
    def get_flow_code_list(self, ARC=False):

    #-------------------------------------------
    # Notes: RT flow codes  = | 64 128 1 |
    #                         | 32  x  2 |
    #                         | 16  8  4 |

    #        ARC/INFO codes = | 32 64 128 |
    #                         | 16  x   1 |
    #                         |  8  4   2 |
    #-------------------------------------------
        if not(ARC):   
            self.code_list = int16([1, 2, 4, 8, 16, 32, 64, 128])
        else:    
            self.code_list = int16([128, 1, 2, 4, 8, 16, 32, 64])
        
    #   get_flow_code_list()
    #---------------------------------------------------------------------
    def get_flow_code_list_opps(self, ARC=False):

        if not(ARC):    
            self.code_opps = int16([16, 32, 64, 128, 1, 2, 4, 8])
        else:    
            self.code_opps = int16([8, 16, 32, 64, 128, 1, 2, 4])
        
    #   get_flow_code_list_opps()
    #---------------------------------------------------------------------
    def get_parent_inc_map(self):

        #-----------------------------------------
        # Note: parent_ID = ID + incs[flow_code].
        #-----------------------------------------
        nx = self.nx
        incs = int32(array([-nx + 1, 1, nx + 1, nx, nx - 1,
                            -1, -nx - 1, -nx]))
        
        MAP = zeros([129], dtype='Int32')
        MAP[self.code_list] = incs
        
        self.inc_map = MAP
        
    #   get_parent_inc_map()
    #---------------------------------------------------------------------
    def get_parent_ID_grid(self):

        nx = self.nx
        ny = self.ny
        ID_grid = reshape(arange(nx*ny, dtype='Int32'), [ny, nx])
        ## self.ID_grid = ID_grid
        
        #-----------------------------------------------------
        # Get a grid which for each grid cell contains the
        # calendar-style index (or ID) of the grid cell that
        # its D8 code says it flows to.
        #-----------------------------------------------------
        print('Finding parent pixel IDs...')

        self.parent_ID_grid = ID_grid + self.inc_map[self.flow_grid]

    #   get_parent_ID_grid()
    #---------------------------------------------------------------------
    def get_parent_IDs(self):

        #-----------------------------------------
        # NB!  The use of 0's here is important.
        #      If iterating, pID[0]=0.
        #-----------------------------------------
       
        #-----------------------------------------------------
        # Get a grid which for each grid cell contains the
        # calendar-style index (or ID) of the grid cell that
        # its D8 code says it flows to.
        #-----------------------------------------------------
        pID_grid = self.parent_ID_grid

        #---------------------------------------
        # Pixels with invalid flow directions,
        # like edges, get assigned a pID of 0.
        #---------------------------------------
        wbad = where(self.flow_grid <= 0)
        nbad = size(wbad[0])
        if (nbad != 0):
            pID_grid[wbad] = 0

        #-------------------------------------------
        # Save IDs as a tuple of row indices and
        # calendar indices, "numpy.where" style
        #------------------------------------------- 
        self.parent_IDs = (pID_grid / self.nx, pID_grid % self.nx)
  
    #   get_parent_IDs()
    #-------------------------------------------------------------------
##    def get_parent_IDs2(flow_grid, rti):
##
##        #-----------------------------------------------------
##        # Note: This version can handle periodic boundaries,
##        #       as can occur in a landscape evolution model.
##        #-----------------------------------------------------
##        # NB!   This was converted from IDL and is not yet
##        #       working.
##        #-----------------------------------------------------
##        
##        #-----------------------------------------
##        # NB!  The use of 0's here is important.
##        #      If iterating, p[0]=0.
##        #-----------------------------------------
##        dirs = flow_grid   # (synonym reference)
##        
##        #---------------------------
##        # Get flow grid dimensions
##        #---------------------------
##        nx = rti.ncols
##        ny = rti.nrows
##        
##        #-----------------------------------------------------
##        # Get IDs of all "downstream" or "parent" pixels
##        #-----------------------------------------------------
##        # i.e. Get a grid which for each grid cell contains
##        # the calendar-style index (or ID) of the grid cell
##        # that its D8 code says it flows to.
##        #-----------------------------------------------------
##        pID_grid = parent_ID_grid(flow_grid, rti)
##        
##        #-----------------------------
##        # Handle periodic boundaries
##        #-----------------------------
##        T = arange(nx, dtype='Int32')
##        B = T + (nx * (ny - int32(1)))
##        L = nx * arange(ny, dtype='Int32')
##        R = L + (nx - int32(1))
##        #----------------------------------
##        #L = nx * (lindgen(ny-2) + 1L)
##        #R = L + (nx-1L)
##        #-------------------------------
##        # Remap pIDs for edge pixels ?
##        #-------------------------------
##        w = where(logical_or(logical_or((dirs[R] == 1), (dirs[R] == 2)), (dirs[R] == 4)))
##        nw = size(w[0])
##        if (nw != 0):    
##            pIDs[R[w]] = (pIDs[R[w]] - nx)
##        #-------------------------------------------------
##        w = where(logical_or(logical_or((dirs[L] == 16), (dirs[L] == 32)), (dirs[L] == 64)))
##        nw = size(w[0])
##        if (nw != 0):    
##            pIDs[L[w]] += nx
##        #-------------------------------------------------
##        NP = (nx * ny)  #*************************************** DOUBLE CHECK **********
##        w = where(logical_or(logical_or((dirs[T] == 1), (dirs[T] == 64)), (dirs[T] == 128)))
##        nw = size(w[0])
##        if (nw != 0):    
##            pIDs[T[w]] += NP
##        #-------------------------------------------------
##        w = where(logical_or(logical_or((dirs[B] == 4), (dirs[B] == 8)), (dirs[B] == 16)))
##        nw = size(w[0])
##        if (nw != 0):    
##            pIDs[B[w]] = (pIDs[B[w]] - NP)
##        
##        #---------------------------------------
##        # Pixels with invalid flow directions,
##        # like edges, get assigned a pID of 0.
##        #---------------------------------------
##        wbad = where(dirs <= uint8(0))
##        nw = size(wbad[0])
##        if (nw != 0):    
##            pIDs[wbad] = int32(0)
##
##        #-----------------------------
##        # Changed for Python version
##        # Simplifies code elsewhere.
##        #-----------------------------
##        pIDs = (pID_grid / nx, pID_grid % nx)
##        return pIDs
##        
##    #   parent_IDs2()
##    #---------------------------------------------------------------------
##    def non_parent_IDs(parent_IDs, flow_grid, rti):
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
##        nbad = size(wbad[0]])
##        if (nbad > 0):
##            base[wbad] = 1   ##########  Should be 1 or 0 ??
##        
##        wnot = where(base == 0)
##        nnot = size(wnot[0])
##        if (nnot != 0):    
##            non_parent_IDs = wnot
##        else:    
##            non_parent_IDs = -int32(1)
##
##        return non_parent_IDs
##
##    #   non_parent_IDs()
    #---------------------------------------------------------------------
    def get_flow_from_IDs(self):

        #----------------------------------------------------------
        # Notes:  This function returns the 4-byte long-integer
        #         array IDs of pixels that flow in a particular
        #         direction.  RiverTools flow codes are assumed.
        #----------------------------------------------------------
        # Notes:  Later, rename w1 to w_NE, n1 to n_NE, then
        #         use self.code_list[0] vs. 1, etc..  This will
        #         then provide support for ARC flow codes, etc.
        #----------------------------------------------------------   
        self.w1 = where(self.flow_grid == 1)
        self.n1 = size(self.w1[0])   #(northeast)
        
        self.w2 = where(self.flow_grid == 2)
        self.n2 = size(self.w2[0])   #(east)
        
        self.w3 = where(self.flow_grid == 4)
        self.n3 = size(self.w3[0])   #(southeast)
        
        self.w4 = where(self.flow_grid == 8)
        self.n4 = size(self.w4[0])   #(south)
        
        self.w5 = where(self.flow_grid == 16)
        self.n5 = size(self.w5[0])   #(southwest)
        
        self.w6 = where(self.flow_grid == 32)
        self.n6 = size(self.w6[0])   #(west)
        
        self.w7 = where(self.flow_grid == 64)
        self.n7 = size(self.w7[0])   #(northwest)
        
        self.w8 = where(self.flow_grid == 128)
        self.n8 = size(self.w8[0])   #(north)

        ##### Same as noflow_IDs  ####################
##        self.w0 = where(self.flow_grid <= 0)
##        self.n0 = size(self.w0[0])   #(undefined)

        # print 'n1 = ' + str(n1)
                
    #   get_flow_from_IDs()
    #---------------------------------------------------------------------
    def get_flow_to_IDs(self):

        nx = self.nx          
       
        #-------------------------------------------------
        # Get IDs of "parent cells" that are downstream
        # of pixels that flow in a given direction.
        #-------------------------------------------------
        if (self.n1 != 0):    # northeast
            p1_IDs  = self.parent_ID_grid[self.w1]
            self.p1 = (p1_IDs / nx, p1_IDs % nx)
        else:
            self.p1 = None
        #-------------------------------------------------
        if (self.n2 != 0):     # east
            p2_IDs  = self.parent_ID_grid[self.w2]
            self.p2 = (p2_IDs / nx, p2_IDs % nx)
        else:
            self.p2 = None
        #-------------------------------------------------
        if (self.n3 != 0):     # southeast
            p3_IDs  = self.parent_ID_grid[self.w3]
            self.p3 = (p3_IDs / nx, p3_IDs % nx)
        else:
            self.p3 = None
        #-------------------------------------------------
        if (self.n4 != 0):     # south 
            p4_IDs  = self.parent_ID_grid[self.w4]
            self.p4 = (p4_IDs / nx, p4_IDs % nx)
        else:
            self.p4 = None
        #-------------------------------------------------
        if (self.n5 != 0):     # southwest
            p5_IDs  = self.parent_ID_grid[self.w5]
            self.p5 = (p5_IDs / nx, p5_IDs % nx)
        else:   
            self.p5 = None
        #-------------------------------------------------
        if (self.n6 != 0):     # west 
            p6_IDs  = self.parent_ID_grid[self.w6]
            self.p6 = (p6_IDs / nx, p6_IDs % nx)
        else:
            self.p6 = None
        #-------------------------------------------------
        if (self.n7 != 0):     # northwest  
            p7_IDs  = self.parent_ID_grid[self.w7]
            self.p7 = (p7_IDs / nx, p7_IDs % nx)
        else:
            self.p7 = None
        #-------------------------------------------------
        if (self.n8 != 0):     # north 
            p8_IDs  = self.parent_ID_grid[self.w8]
            self.p8 = (p8_IDs / nx, p8_IDs % nx)
        else:
            self.p8 = None
        #-------------------------------------------------
        ##    print 'p1.shape, size(p1)       =', p1.shape, size(p1)
        ##    print 'w1[0].shape, size(w1[0]) =', w1[0].shape, size(w1[0])      
        
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

    #   get_flow_to_IDs()
    #-------------------------------------------------------------------
    def get_edge_IDs(self):

        #---------------------------
        # Get flow grid dimensions
        #---------------------------
        nx = self.rti.ncols
        ny = self.rti.nrows
        
        print('Finding edge pixel IDs...')

        #-------------------------
        # Get IDs of edge pixels
        #-------------------------    
        T_IDs = arange(nx, dtype='Int32')
        B_IDs = T_IDs + (ny - int32(1)) * nx
        L_IDs = (int32(1) + arange(ny - int32(2), dtype='Int32')) * nx
        R_IDs = L_IDs + (nx - int32(1))
        edge_IDs = concatenate([T_IDs, B_IDs, L_IDs, R_IDs])

        #-------------------------------------------
        # Save IDs as a tuple of row indices and
        # calendar indices, "numpy.where" style
        #-------------------------------------------        
        self.edge_IDs = (edge_IDs / nx, edge_IDs % nx)   ##  NB! (row, col)

        #------------------------------------------
        # Save IDs as a 1D array of long-integer,
        # calendar-style indices
        #------------------------------------------
        # self.edge_IDs = edge_IDs
        
    #   get_edge_IDs()
    #-------------------------------------------------------------------
    def get_noflow_IDs(self):
        
        #--------------------------------------------------
        # 1/19/07.  Need to set d and u to zero at any ID
        # where flow terminates.  This includes pixels on
        # the edges, those with unresolved flow direction
        # and those where elevation is nodata or NaN.
        # A RiverTools flow grid will have a flow code of
        # zero at all of these places.
        #--------------------------------------------------
        noflow_IDs = where(self.flow_grid <= 0)
        num_IDs    = size(noflow_IDs[0])
        
        if (num_IDs != 0):
            self.noflow_IDs = noflow_IDs
        else:
            #----------------------------
            # Return IDs of edge pixels
            #----------------------------
            self.get_edge_IDs()
            self.noflow_IDs = self.edge_IDs()

    #   get_noflow_IDs()
    #-------------------------------------------------------------------
    def get_flow_width_grid(self, DOUBLE=False, METHOD2=False):

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
        # NB!    numpy.where returns a 2-tuple, but with "empty"
        #        second part when applied to a 1D array.  This means
        #        that nw = size(w[0]) still works.
        #-------------------------------------------------------------        
        print('Computing flow width grid...')

        #-----------------------
        # Get pixel dimensions
        #-----------------------
        dx, dy, dd = pixels.get_sizes_by_row(self.rti, METERS=True)

        #-------------------------
        # Double or Float type ?
        #-------------------------
        if (DOUBLE):    
            dw = zeros([self.ny, self.nx], dtype='Float64')
        else:    
            dw = zeros([self.ny, self.nx], dtype='Float32')
            dx = float32(dx)
            dy = float32(dy)
            dd = float32(dd)

        #----------------------------------------------
        # Initialize to default value that is used
        # for pixels with flow code of zero.  This is
        # done to avoid "divide by zero" errors.
        #----------------------------------------------
        dw += dx[0]
        ## dw = dw + dx.min()
        
        for row in range(self.ny):
            g = self.flow_grid[row,:]
            
            #----------------
            # Diagonal flow
            #----------------
            wd  = where(logical_or(logical_or(logical_or((g == 1), (g == 4)), \
                                              (g == 16)), (g == 64)))
            nwd = size(wd[0])
            if (nwd != 0):    
                if not(METHOD2):    
                    dw[row, wd] = dd[row]
                else:    
                    dw[row, wd] = (dx[row] + dy[row]) / 4
            
            #---------------------
            # East and west flow
            #---------------------
            wh  = where(logical_or((g == 2), (g == 32)))
            nwh = size(wh[0])
            if (nwh != 0):    
                dw[row, wh] = dy[row]
                if (METHOD2):    
                    dw[row, wh] = dw[row,wh] / 2
            
            #-----------------------
            # North and south flow
            #-----------------------
            wv  = where(logical_or((g == 8), (g == 128)))
            nwv = size(wv[0])
            if (nwv != 0):    
                dw[row, wv] = dx[row]
                if (METHOD2):    
                    dw[row, wv] = dw[row, wv] / 2

        #----------
        # Report
        #----------
        dw_min = dw.min()
        dw_max = dw.max()
        TF_Print('    min(dw) = ' + str(dw_min) + '  [m]')
        TF_Print('    max(dw) = ' + str(dw_max) + '  [m]')

        self.dw = dw

    #   get_flow_width_grid()
    #-------------------------------------------------------------------
    def get_flow_length_grid(self, DOUBLE=False):

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
        TF_Print('Computing flow length grid...')

        #-----------------------
        # Get pixel dimensions
        #-----------------------
        dx, dy, dd = pixels.get_sizes_by_row(self.rti, METERS=True)
            
        #-------------------------
        # Double or Float type ?
        #-------------------------
        if (DOUBLE):    
            ds = zeros([self.ny, self.nx], dtype='Float64')
        else:    
            ds = zeros([self.ny, self.nx], dtype='Float32')
            dx = float32(dx)
            dy = float32(dy)
            dd = float32(dd)
        
        #----------------------------------------------
        # Initialize to default value that is used
        # for pixels with flow code of zero.  This is
        # done to avoid "divide by zero" errors.
        #----------------------------------------------
        ds += dx[0]
        ## ds += dx.min()
        
        for row in range(self.ny):
            g = self.flow_grid[row,:]
            
            #----------------
            # Diagonal flow
            #----------------
            wd  = where(logical_or(logical_or(logical_or((g == 1), (g == 4)), \
                                              (g == 16)), (g == 64)))
            nwd = size(wd[0])
            if (nwd != 0):    
                ds[row, wd] = dd[row]
            
            #---------------------
            # East and west flow
            #---------------------
            wh  = where(logical_or((g == 2), (g == 32)))
            nwh = size(wh[0])
            if (nwh != 0):    
                ds[row, wh] = dx[row]
            
            #-----------------------
            # North and south flow
            #-----------------------
            wv  = where(logical_or((g == 8), (g == 128)))
            nwv = size(wv[0])
            if (nwv != 0):    
                ds[row, wv] = dy[row]
            
        #----------
        # Report
        #----------
        ds_min = ds.min()
        ds_max = ds.max()
        TF_Print('    min(ds) = ' + str(ds_min) + '  [m]')
        TF_Print('    max(ds) = ' + str(ds_max) + '  [m]')

        self.ds = ds

    #   get_flow_length_grid()
    #-------------------------------------------------------------------



