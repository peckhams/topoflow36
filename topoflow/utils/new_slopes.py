#------------------------------------------------------------------------
#  Copyright (c) 2019, Scott D. Peckham
#
#  Oct 2019.  First version to put in utils folder.
#             remove_bad_slopes() is from channels_base.py.

#------------------------------------------------------------------------

#  test1()

#  get_d8_object()
#  get_d8_slope_grid()

#  get_new_slope_grid()

#  remove_bad_slopes()    (old approach)

#------------------------------------------------------------------------

import numpy as np
from . import rtg_files   # (also in utils directory)
from . import rti_files   # (also in utils directory)

from topoflow.components import d8_global

#---------------------------------------------------------------------------
def test1():

    site_prefix = 'Baro_Gam_1min'
    case_prefix = 'Test1'
    cfg_dir = '/Users/peckhams/Desktop/TopoFlow_2019/__Regions/'
    cfg_dir += 'Ethiopia/DEMs/MERIT/Baro_Gam_1min/'
    #----------------------------------------------------
    old_slope_file = (site_prefix + '_old-slope.rtg')
    new_slope_file = (site_prefix + '_new-slope.rtg')
    
    #---------------------------------------------------------- 
    # Compute simple D8 slope grid, with zero slopes in flats
    #----------------------------------------------------------
    get_d8_slope_grid(site_prefix=site_prefix, case_prefix=case_prefix,
                    cfg_dir=cfg_dir, slope_file=old_slope_file)

    #--------------------------------------------------------
    # Compute new, improved slope grid, with no zero slopes
    #--------------------------------------------------------                          
    get_new_slope_grid(site_prefix=site_prefix, case_prefix=case_prefix,
                       cfg_dir=cfg_dir, slope_file=new_slope_file)

#   test1()
#---------------------------------------------------------------------------
def get_d8_object(site_prefix=None, case_prefix=None,
                  cfg_dir=None,
                  SILENT=True, REPORT=False):

    #---------------------------------------------
    # Compute and store a variety of (static) D8
    # flow grid variables.  Embed structure into
    # the "channel_base" component.
    #---------------------------------------------
    d8 = d8_global.d8_component()
    
    #--------------------------------------------------         
    # D8 component builds its cfg filename from these  
    #-------------------------------------------------------------
    # (2/11/2017) The initialize() method in d8_base.py now
    # uses case_prefix (vs. site_prefix) for its CFG file:
    # <site_prefix>_d8_global.cfg.  This is to prevent confusion
    # since this was the only CFG file that used site_prefix.
    #-------------------------------------------------------------    
    d8.site_prefix  = site_prefix
    d8.case_prefix  = case_prefix   # (used in d8_base.py)
    d8.in_directory = cfg_dir
    cfg_file        = cfg_dir + case_prefix + '_d8_global.cfg'
    
    d8.initialize( cfg_file=cfg_file, SILENT=SILENT, \
                   REPORT=REPORT )
    
    #---------------------------------------------------
    # The next 2 "update" calls are needed when we use
    # the new "d8_base.py", but are not needed when
    # using the older "tf_d8_base.py".      
    #---------------------------------------------------
    time = 0
    d8.update(time, SILENT=SILENT, REPORT=REPORT)

    #----------------------------------------------------------- 
    # Note: This is also needed, but is not done by default in
    #       d8.update() because it hurts performance of Erode.
    #----------------------------------------------------------- 
    d8.update_noflow_IDs()
 
    return d8

#   get_d8_object()
#---------------------------------------------------------------------------
def get_d8_slope_grid(site_prefix=None, case_prefix=None,
                    cfg_dir=None, slope_file=None):

    d8 = get_d8_object(site_prefix=site_prefix,
                       case_prefix=case_prefix,
                       cfg_dir=cfg_dir)

    slope = np.zeros([d8.ny, d8.nx], dtype='float64')
     
    z2   = d8.DEM
    pIDs = d8.parent_IDs
    z1   = d8.DEM[ pIDs ]  
    dz   = (z2 - z1)
    wg   = (dz > 0)   # (array of True or False)
    ng   = wg.sum()   # (number of cells with True)
    wb   = np.invert( wg )
    nb   = wb.sum()  # (number of cells with True)
    L    = d8.ds.copy()
    if (ng > 0):
        slope[ wg ] = (dz[wg] / L[wg])

    #-----------------------------------------------------------
    # Set slope to NaN at all "noflow_IDs" (D8 flow code of 0)
    #-----------------------------------------------------------
    slope[ d8.noflow_IDs ] = np.nan
    
    print('Computing simple D8 slope grid...')
    Spmin = slope[ slope > 0 ].min()
    print('   Min slope     = ' + str(slope.min()) )
    print('   Min pos slope = ' + str(Spmin) )
    print('   Max slope     = ' + str(slope.max()) )

    #---------------------------------------
    # Read header_file info for slope_file
    #---------------------------------------
    header_file = (cfg_dir + site_prefix + '.rti')
    grid_info = rti_files.read_info( header_file, REPORT=False)
    
    #------------------------------------
    # Save new slope grid to slope_file
    #------------------------------------             
    slope_file2 = (cfg_dir + slope_file)
    rtg_files.write_grid( slope, slope_file2, grid_info, RTG_type='FLOAT')
    print( 'Finished writing D8 slope grid to file: ')
    print( slope_file2 )
    print()

#   get_d8_slope_grid()                             
#---------------------------------------------------------------------------
def get_new_slope_grid(site_prefix=None, case_prefix=None,
                       cfg_dir=None, slope_file=None):

    #---------------------------------------------------------
    # NB!  The iteration of D8 parent IDs will only terminate
    #      if we have pID[0]=0.   See d8_global.py.
    #---------------------------------------------------------
    # Idea2: If your D8 parent cell (downstream) has S=0,
    #        then iterate to parent of parent until you
    #        reach a cell whose D8 parent has S > 0.
    #        This cell can now be assigned a slope based
    #        on start z and stop z and distance.
    #        Connected cells upstream with S=0 should be
    #        assigned this same slope.
    #        Or, at this point, could compute a new z,
    #        then repeat the process.
    #---------------------------------------------------------
    d8 = get_d8_object(site_prefix=site_prefix,
                       case_prefix=case_prefix,
                       cfg_dir=cfg_dir)

    slope = np.zeros([d8.ny, d8.nx], dtype='float64')
     
    z2   = d8.DEM
    pIDs = d8.parent_IDs
    z1   = d8.DEM[ pIDs ]  
    dz   = (z2 - z1)
    wg   = (dz > 0)   # (array of True or False)
    ng   = wg.sum()   # (number of cells with True)
    wb   = np.invert( wg )
    nb   = wb.sum()  # (number of cells with True)
    L    = d8.ds.copy()
    if (ng > 0):
        slope[ wg ] = (dz[wg] / L[wg])

    print('Computing slope grid...')
    Spmin = slope[ slope > 0 ].min()
    print('   initial min slope     = ' + str(slope.min()) )
    print('   initial min pos slope = ' + str(Spmin) )
    print('   initial max slope     = ' + str(slope.max()) )
    
    #-----------------------------------------------------
    # Get a grid which for each grid cell contains the
    # calendar-style index (or ID) of the grid cell that
    # its D8 code says it flows to.
    #-----------------------------------------------------
    ### pID_grid = d8.parent_ID_grid.copy()
    pID_grid = d8.parent_ID_grid
    ID_grid  = d8.ID_grid

    #-----------------------------------------
    # Save IDs as a tuple of row indices and
    # column indices, "np.where" style
    # divmod() is builtin and not in numpy.
    #-----------------------------------------
    # parent_IDs = divmod(pID_grid, d8.nx)

    #-------------------------------           
    # Get slopes of the D8 parents
    #-------------------------------
    S0 = slope.copy()
    p_slope = slope[ pIDs ]
    w1 = np.logical_and( slope > 0, p_slope == 0 )
    n1 = w1.sum()
    IDs = w1      ##################
    L = L[ IDs ]  #################
    unfinished = np.zeros( n1 ) + 1    ########## 
    n_left = n1          
    max_reps = 500  #######
    n_reps = 0                   
    DONE = (n1 == 0)
    while not(DONE):
        
        #----------------------------------------------------
        # Here, IDs & pIDs are tuples, like from np.where,
        # while ID_vals & pID_vals are long-integer indices.
        #----------------------------------------------------
        pID_vals = pID_grid[ IDs ]
        pIDs = divmod( pID_vals, d8.nx )  # (tuple, like WHERE)

        z1 = d8.DEM[ pIDs ]
        dz = (z2[IDs] - z1)
        ds = d8.ds[ pIDs ]
        L  += ds

        p_slope = S0[ pIDs ]
        ID_vals = d8.ID_grid[ IDs ]  ######
        ## wg = (p_slope > 0)        ###############
        wg = np.logical_and( p_slope > 0, unfinished == 1 )
        ng = wg.sum()
        if (ng > 0):
            ready_ID_vals = ID_vals[ wg ]
            ready_IDs = divmod( ready_ID_vals, d8.nx )
            slope[ ready_IDs ] = dz[ wg ] / L[ wg ]
            unfinished[ wg ] = 0
            n_left = unfinished.sum()
        IDs = pIDs
        
        n_reps += 1
        if (n_reps == max_reps):
            print('Aborting after ' + str(max_reps) + ' reps.') 
        DONE = (n_left == 0) or (pID_vals.max() == 0)   ########    
        ### DONE = (n_reps == max_reps) or (pID_vals.max() == 0)   ########

    #--------------------------------------------------------------    
    # Step 2.  Assign slope to the in-between grid cells with S=0
    #--------------------------------------------------------------
    # If D8 parent was assigned a slope in Step 1, then any of
    # its D8 kids with S=0 should be assigned the same slope.
    # Iterate until all cells with S=0 are assigned a slope.
    #--------------------------------------------------------------
#     pIDs = d8.parent_IDs
#     ## pID_grid = d8.parent_ID_grid.copy()   
#     # w1 = (slope == 0)
#     # n1 = w1.sum()
#     # DONE = (n1 == 0)
#     DONE = False
#     while not(DONE): 
#         pslope = slope[ pIDs ]
#         w2 = np.logical_and( slope == 0, pslope > 0 )
#         n2 = w2.sum()
#         if (n2 > 0):
#             slope[ w2 ] = pslope[ w2 ]
#         DONE = (n2 == 0)
        #-----------------------------------
        # Get IDs for "parents of parents"
        #-----------------------------------
        ## pID_grid = pID_grid[ pIDs ]
        ## pIDs = divmod( pID_grid, d8.nx )

    #-----------------------------------------------------------
    # Set slope to NaN at all "noflow_IDs" (D8 flow code of 0)
    #-----------------------------------------------------------
    slope[ d8.noflow_IDs ] = np.nan
    
    #---------------------------------------
    # Read header_file info for slope_file
    #---------------------------------------
    header_file = (cfg_dir + site_prefix + '.rti')
    grid_info = rti_files.read_info( header_file, REPORT=False)
    
    #------------------------------------
    # Save new slope grid to slope_file
    #------------------------------------             
    slope_file2 = (cfg_dir + slope_file)
    rtg_files.write_grid( slope, slope_file2, grid_info, RTG_type='FLOAT')
    print( 'Finished writing new slope grid to file: ')
    print( slope_file2 )
    print()
    ### print('Finished computing slope grid.')
    Spmin = slope[ slope > 0 ].min()
    print('   n_reps        = ' + str(n_reps) )
    print('   min slope     = ' + str(slope.min()) )
    print('   min pos slope = ' + str(Spmin) )
    print('   max slope     = ' + str(slope.max()) )
    print(' ')
          
#   get_new_slope_grid()     
#---------------------------------------------------------------------------
def remove_bad_slopes(slope, FLOAT=False):

    #------------------------------------------------------------
    # Notes: The main purpose of this routine is to find
    #        pixels that have nonpositive slopes and replace
    #        then with the smallest value that occurs anywhere
    #        in the input slope grid.  For example, pixels on
    #        the edges of the DEM will have a slope of zero.

    #        With the Kinematic Wave option, flow cannot leave
    #        a pixel that has a slope of zero and the depth
    #        increases in an unrealistic manner to create a
    #        spike in the depth grid.

    #        It would be better, of course, if there were
    #        no zero-slope pixels in the DEM.  We could use
    #        an "Imposed gradient DEM" to get slopes or some
    #        method of "profile smoothing".

    #        It is possible for the flow code to be nonzero
    #        at a pixel that has NaN for its slope. For these
    #        pixels, we also set the slope to our min value.

    #        7/18/05. Broke this out into separate procedure.
    #------------------------------------------------------------

    #-----------------------------------
    # Are there any "bad" pixels ?
    # If not, return with no messages.
    #-----------------------------------  
    wb = np.where(np.logical_or((slope <= 0.0), \
                          np.logical_not(np.isfinite(slope))))
    nbad = np.size(wb[0])
    print('size(slope) = ' + str(np.size(slope)) )
    print('size(wb) = ' + str(nbad) )
    
    wg = np.where(np.invert(np.logical_or((slope <= 0.0), \
                                 np.logical_not(np.isfinite(slope)))))
    ngood = np.size(wg[0])
    if (nbad == 0) or (ngood == 0):
        return
    
    #---------------------------------------------
    # Find smallest positive value in slope grid
    # and replace the "bad" values with smin.
    #---------------------------------------------
    print('-------------------------------------------------')
    print('WARNING: Zero or negative slopes found.')
    print('         Replacing them with smallest slope.')
    print('         Use "Profile smoothing tool" instead.')
    S_min = slope[wg].min()
    S_max = slope[wg].max()
    print('         min(S) = ' + str(S_min))
    print('         max(S) = ' + str(S_max))
    print('-------------------------------------------------')
    print(' ')
    slope[wb] = S_min
    
    #--------------------------------
    # Convert data type to double ?
    #--------------------------------
    if (FLOAT):    
        slope = np.float32(slope)
    else:    
        slope = np.float64(slope)
    
    return slope
    
#   remove_bad_slopes()
#---------------------------------------------------------------------------



