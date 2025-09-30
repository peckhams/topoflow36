#------------------------------------------------------------------------
#  Copyright (c) 2019-2023, Scott D. Peckham
#
#  Aug 2023.  Search for 2023-08-24 for extra exit condition. 
#  Oct 2019.  First version to put in utils folder.
#             remove_bad_slopes() is from channels_base.py.

#------------------------------------------------------------------------

#  test1()

#  get_d8_object()
#  get_d8_slope_grid()

#  get_new_slope_grid()
#  get_new_slope_grid0()   # doesn't work

#  smooth_dem()     # experimental

#  remove_bad_slopes()    (old approach)

#------------------------------------------------------------------------

import numpy as np
import copy, os
from . import rtg_files   # (also in utils directory)
from . import rti_files   # (also in utils directory)

from topoflow.components import d8_global

#---------------------------------------------------------------------------
def test1():

    site_prefix = 'Baro_Gam_1min'
    case_prefix = 'Test1'
    cfg_dir = '/Users/peckhams/TF_Tests3/Baro-Gam_60sec/'
    cfg_dir += 'Test1_2015-10_to_2018-10_CHIRPS_Bilinear_cfg'
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
def test2():

    site_prefix = 'Baro_Gam_1min'
    case_prefix = 'Test1'
    cfg_dir = '/Users/peckhams/Desktop/TopoFlow_2019/__Regions/'
    cfg_dir += 'Ethiopia/DEMs/MERIT/Baro_Gam_1min/'
    new_dem_file = site_prefix + '_smooth_DEM.rtg'
    
    smooth_dem(site_prefix=site_prefix, case_prefix=case_prefix, 
               cfg_dir=cfg_dir, dem_file=new_dem_file,
               nsteps=10)
    
#   test2()
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
    # Note that D8 area is also computed as d8.A.     
    #---------------------------------------------------
    time = 0
    d8.update(time, REPORT=REPORT)

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
    topo_dir = d8.topo_directory

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

    #------------------------------------------------
    # Option to return grid, without saving to file
    #------------------------------------------------
    if (slope_file is None):
        return slope

    #---------------------------------------
    # Read header_file info for slope_file
    #---------------------------------------
    header_file = (topo_dir + site_prefix + '.rti')
    grid_info = rti_files.read_info( header_file, REPORT=False)
    
    #------------------------------------
    # Save new slope grid to slope_file
    #------------------------------------
    slope_file2 = slope_file
    if not(os.sep in slope_file):             
        slope_file2 = (topo_dir + slope_file)   
    rtg_files.write_grid( slope, slope_file2, grid_info, RTG_type='FLOAT')
    print( 'Finished writing D8 slope grid to file: ')
    print( slope_file2 )
    print()

#   get_d8_slope_grid()
#---------------------------------------------------------------------------
def get_new_slope_grid(site_prefix=None, case_prefix=None,
                       cfg_dir=None, slope_file='TEST_newslope.rtg'):

    #---------------------------------------------------------
    # NB!  New idea:  2020-09-24.
    #---------------------------------------------------------
    # NB!  The iteration of D8 parent IDs will only terminate
    #      if we have pID[0]=0.   See d8_global.py.
    #---------------------------------------------------------
    # Idea: If your D8 parent cell (downstream) has the same
    #       elevation as you do, then iterate to parent of
    #       parent until you reach a cell with a lower
    #       elevation.
    #       This cell can now be assigned a slope based
    #       on start z and stop z and distance.
    #---------------------------------------------------------
    d8 = get_d8_object(site_prefix=site_prefix,
                       case_prefix=case_prefix,
                       cfg_dir=cfg_dir)
    topo_dir = d8.topo_directory
    
    print('Computing slope grid...')
    print('dtype(d8.ds) =', d8.ds.dtype)  # (float32) 

    #-----------------------------------------------------
    # Get a grid which for each grid cell contains the
    # calendar-style index (or ID) of the grid cell that
    # its D8 code says it flows to.
    #-----------------------------------------------------
    ### pID_grid = d8.parent_ID_grid.copy()
    z2       = d8.DEM
    pIDs     = d8.parent_IDs
    pID_grid = d8.parent_ID_grid
    ID_grid  = d8.ID_grid
    # Change from float32 to float64 to avoid overflow.
    distance = np.float64( d8.ds )   # (initialize distance grid)
    slope    = np.zeros([d8.ny, d8.nx], dtype='float64')
    START    = True
    n_reps   = 0
  
         
    DONE = False
    while not(DONE):   
        z1 = d8.DEM[pIDs]
        dz = (z2 - z1)
        #-------------------------------------------------------
        # Handle both nodata & edge values in DEM (2022-05-02)
        # Otherwise, there are large, erroneous slopes.
        # Issues show up best with a "linear" stretch.
        #-------------------------------------------------------
        w1 = np.logical_or(z1 <= d8.DEM_nodata, d8.d8_grid == 0)
        # w1 = (d8.d8_grid == 0)  # (doesn't work)
        # w1 = np.logical_or(z1 <= d8.DEM_nodata, pIDs == 0)  # (doesn't work)
        dz[ w1 ] = 0

        #-----------------------------------------------        
        # Note that all slopes are initialized to zero
        # and z2 = DEM doesn't change.
        #-----------------------------------------------
        wg = np.logical_and(slope == 0, dz > 0)   # (array of True or False) 
        ng = wg.sum()   # (number of cells with True)
        wb = np.invert( wg )
        nb = wb.sum()   # (number of cells with True)
        if (ng > 0):
            slope[ wg ] = (dz[wg] / distance[wg])

        #-----------------------------
        # Print initial slope values
        #-----------------------------
        if (START):
            Sp_min = slope[ slope > 0 ].min()
            print('   Initial min slope     = ' + str(slope.min()) )
            print('   Initial min pos slope = ' + str(Sp_min) )
            print('   Initial max slope     = ' + str(slope.max()) )
            print()
            START = False

        #-------------------------------------------  
        # Update downstream distance (before pIDs)
        #-------------------------------------------
        distance += d8.ds[ pIDs ]
 
        #--------------------------------------        
        # Update pIDs (D8 parents of parents)
        #--------------------------------------
        pID_vals = pID_grid[ pIDs ]
        pIDs     = divmod( pID_vals, d8.nx)
        pv_max   = pID_vals.max()
        n_reps += 1

        #---------------------------------------
        # Added nreps condition on 2023-08-24.
        # Occurred due to DEM data type issue.
        #---------------------------------------
        DONE = (nb == 0) or (pv_max == 0) or (n_reps > 4*d8.nx)

    if (n_reps > d8.nx):
        print('   ERROR:  Iteration terminated because')
        print('   n_reps > 4 * ncols. Check DEM data type.')
        print()

    #--------------------------
    # Print final slope values
    #--------------------------
    Sp_min = slope[ slope > 0 ].min()
    print('   Final min slope     = ' + str(slope.min()) )
    print('   Final min pos slope = ' + str(Sp_min) )
    print('   Final max slope     = ' + str(slope.max()) )
                
    #-----------------------------------------------------------
    # Set slope to NaN at all "noflow_IDs" (D8 flow code of 0)
    #-----------------------------------------------------------
    # slope[ d8.noflow_IDs ] = np.nan
    # slope[ d8.noflow_IDs ] = 0.0    ### Can change min or max slope
    
    #---------------------------------------
    # Read header_file info for slope_file
    #---------------------------------------
    header_file = (topo_dir + site_prefix + '.rti')
    grid_info = rti_files.read_info( header_file, REPORT=False)
    
    #-----------------------------------------------
    # Save new slope grid to slope_file
    # Don't add topo_dir; may already be full path
    #-----------------------------------------------
    if not(os.sep in slope_file):             
        slope_file2 = (topo_dir + slope_file)
    else:
        slope_file2 = slope_file
    rtg_files.write_grid( slope, slope_file2, grid_info, RTG_type='FLOAT')
    print( 'Finished writing new slope grid to file: ')
    print( slope_file2 )
    print()

    Spmin = slope[ slope > 0 ].min()
    print('   n_reps        = ' + str(n_reps) )
    print('   min slope     = ' + str(np.nanmin(slope)) )
    print('   min pos slope = ' + str(Spmin) )
    print('   max slope     = ' + str(np.nanmax(slope)) )
    print(' ')
          
#   get_new_slope_grid()
#---------------------------------------------------------------------------
def get_new_slope_grid2(site_prefix=None, case_prefix=None,
                        cfg_dir=None, slope_file='TEST_newslope.rtg',
                        REPORT=True):

    #---------------------------------------------------------
    # NB!  There seems to be a bug in get_new_slope_grid()
    #      so this is a new attempt. (2025-02-02)
    #---------------------------------------------------------
    # NB!  The iteration of D8 parent IDs will only terminate
    #      if we have pID[0]=0.   See d8_global.py.
    #---------------------------------------------------------
    # Idea: If your D8 parent cell (downstream) has the same
    #       elevation as you do, then iterate to parent of
    #       parent until you reach a cell with a lower
    #       elevation.
    #       This cell can now be assigned a slope based
    #       on start z and stop z and distance.
    #---------------------------------------------------------
    print('Getting D8 object...')
    d8 = get_d8_object(site_prefix=site_prefix,
                       case_prefix=case_prefix,
                       cfg_dir=cfg_dir)
    topo_dir = d8.topo_directory
    
    if (REPORT):
        print('Computing improved slope grid...')
        # print('  dtype(d8.ds) =', d8.ds.dtype)  # (float32) 

    #-----------------------------------------------------
    # Get a grid which for each grid cell contains the
    # calendar-style index (or ID) of the grid cell that
    # its D8 code says it flows to.
    #-----------------------------------------------------
    # parent_ID is set to zero for: (1) edge cells and
    # (2) cells with D8 flow code of 0 (undefined)
    #-----------------------------------------------------
    # pIDs is likely to contain duplicates, which is OK.
    #-----------------------------------------------------        
    z2       = d8.DEM
    pIDs     = d8.parent_IDs   # a (rows,cols) tuple like np.where
    pID_grid = d8.parent_ID_grid
    ID_grid  = d8.ID_grid
    slope    = np.zeros([d8.ny, d8.nx], dtype='float64')
    n_reps   = 0
   
    ## print('pIDs.type  =', type(pIDs))  # (tuple)

    #------------------
    # Compute dz grid
    #------------------
    z1 = d8.DEM[pIDs]
    dz = (z2 - z1)

    #-------------------------------------------------------
    # Handle both nodata & edge values in DEM (2022-05-02)
    # Otherwise, there are large, erroneous slopes.
    # Issues show up best with a "linear" stretch.
    #-------------------------------------------------------
    ww1 = (z2 <= d8.DEM_nodata)
    nw1 = ww1.sum()
    ww2 = (z1 <= d8.DEM_nodata)
    nw2 = ww2.sum()
    ww3 = (d8.d8_grid == 0)
    nw3 = ww3.sum()
    if (REPORT):
       print('Number of cells in DEM w/ nodata values =', nw1)
       print('Number of cells with D8 flow code of 0  =', nw3)
    #----------------
    dz[ ww1 ] = 0.0    
    dz[ ww2 ] = 0.0
    dz[ ww3 ] = 0.0

    #-------------------------------------------------------    
    # Compute D8 slopes for cells whose D8 parent is lower
    #-------------------------------------------------------
    ds = np.float64(d8.ds)  # Change float32 to float64; avoid overflow
    wg = (dz > 0)   # (array of True or False)
    ng = wg.sum()   # (number of cells with True)
    wb = np.invert( wg )  # where dz == 0.
    nb = wb.sum()  # (number of cells with True)
    if (ng > 0):
        slope[ wg ] = (dz[wg] / d8.ds[wg])
    print('DEM (ncols, nrows) =', d8.nx, ', ', d8.ny)
    print('min(ds), max(ds) =', ds.min(), ', ', ds.max())
    #------------------------------------------------------------------    
    # IMPORTANT!  Cell with ID = 1406120 has (row,col)=(1269,68).
    # This is in the bottom (southernmost) row of the DEM.  Cells in
    # the row above have D8 code of 8 and flow south.
    # The pID for every edge cell is 0, so there could be a dz
    # value > 0 between an edge cell and the upper left corner cell
    # that has ID = 0.
    #------------------------------------------------------------------
    nodata = -999.0
    slope[ ww3 ] = nodata  # set back to zero at end

    #------------------------------------------------    
    # This is used to prevent new slope grid values
    # from affecting the new slope assignments
    # But need to allow overwriting values in slope
    # grid with smaller values.
    #------------------------------------------------
    d8_slope = slope.copy()

    if (REPORT):
        print('Number of cells with nonzero slope =', ng)
        print('Number of cells with zero slope    =', nb)
        print()
    if (nb == 0):
        print('Returning.  There are no cells with dz <= 0.')
        return

    #-----------------------------
    # Print initial slope values
    #-----------------------------
    S_min   = slope[ slope > nodata ].min()
    Sp_min1 = slope[ slope > 0 ].min()
    Sp_min2 = slope[ slope > Sp_min1 ].min()
    Sp_min3 = slope[ slope > Sp_min2 ].min()
    if (REPORT):
        print('Initial min slope = ' + str(S_min) )
        print('Initial max slope = ' + str(slope.max()) )
        print('Initial smallest positive slopes = ')
        print('   ' + str(Sp_min1) )
        print('   ' + str(Sp_min2) )
        print('   ' + str(Sp_min3) )
        print()
    
    #-------------------------------------------
    # Get first cells in the zero-slope chains
    # Note that "ID_slope" = slope grid
    #-------------------------------------------
    # Both slope and pID_slope are grids here
    # IDs & pIDs are tuples: (rows, cols)
    #-----------------------------------------
    pID_slope = slope[ pIDs ]  # this is a grid here
    print('## pID_slope.shape =', pID_slope.shape)
    w0 = (pID_slope == 0)
    n0 = w0.sum()
    normal = np.logical_and(slope > 0, pID_slope == 0)

    #---------------------------------------------------------
    # There are cells that are at start of a zero-slope
    # chain and have slope == 0, but no D8 children and
    # therefore no upstream D8 child with higher elevation.
    # These cells are referred to as "special" vs. "normal".
    #---------------------------------------------------------
    # Can resolve them in the same way, by following the
    # chain downstream until a drop in elevation.  However,
    # there are a few such cells such that their par_slope
    # is nonzero.  We can change (pID_slope == 0) to
    # (slope == 0) to get these as well.
    # The
    #---------------------------------------------------------
    da = np.float32(d8.da / 1e6)   # m2 to km2
    ## special = np.logical_and(da == d8.A, pID_slope == 0)  
    special = np.logical_and(da == d8.A, slope == 0)
    print('## n_normal  =', normal.sum())
    print('## n_special =', special.sum())
#     aa = (da == d8.A)
#     print('## n_da_eq_A =', aa.sum())
#     print('## min(d8.da), max(d8.da) =', da.min(), ', ', da.max())
#     print('## minpos(d8.A) =', d8.A[d8.A > 0].min())
#     print()
    IDs = np.where( np.logical_or( normal, special ))
    ## IDs = np.where( np.logical_and(slope > 0, pID_slope == 0) )
    n1  = IDs[0].size
    #---------------------------------------------------------
#     w1 = np.logical_and(slope > 0, pID_slope == 0)  # both are grids 
#     n1 = w1.sum()

    #--------------------------------------------------
    # Check slopes assigned to "special" flow paths.
    #--------------------------------------------------
    # Many of these flow paths flow into a path that
    # is resolved as a "normal" flow path, so the
    # flat elevation may persist far downstream from
    # the last unresolved 0-slope cell.
    #--------------------------------------------------
    # IDs of cells with no D8 children and slope == 0
    # in the Beaver3 MERIT DEM:
    #    46989 = (249, 114), 6 0-slope cells,
    #       start elevation = 202.722549
    #       end elevation   = 201.598404 (245, 97) many cells further
    #       flows into cell with new slope = 0.00103649
    #       new slope = 0.000586
    #-------------------------------------------------- 
    #    67785 = (135, 165), 4 0-slope cells
    #       start elevation = 211.004166
    #       end elevation   = 210.350250 (149, 137) many cells further
    #       flows into cell with new slope = 0.000188118
    #       new slope = 0.000136
    #-------------------------------------------------- 
    #    49454 = (254, 120), 4 0-slope cells
    #       start elevation =
    #       end elevation   = 
    #       flows into cell with slope = 0.000624501
    #       new slope = 
    #-------------------------------------------------- 
    #    49984 = (374, 121), 4 0-slope cells
    #       start elevation =
    #       end elevation   = 
    #       flows into cell with slope = 0.000597519
    #       new slope = 
    #-------------------------------------------------- 
    #    54080 = (370, 131), 4 0-slope cells
    #       start elevation =
    #       end elevation   = 
    #       flows into cell with slope = 0.00103120
    #       new slope = 
    #--------------------------------------------------  

    if (REPORT):
       print('Number of cells with zero pID_slope =', n0)
       print('Number of start cells =', n1)
       print()
    if (n1 > 0):
        #---------------------------------------------------------
        # Get D8 parent IDs as first ones on zero-slope chains
        #---------------------------------------------------------
        # Should not remove duplicate IDs here because start
        # ID elevations and distances are needed, and we need
        # for all the 1D arrays from here to have same length.
        #---------------------------------------------------------
        # Multiple IDs often have the same D8 parent, so while
        # ID_vals are unique, pID_vals are not.
        # Multiple IDs with same D8 parent can have slope > 0
        # and their D8 parent may have a slope of zero.
        #---------------------------------------------------------
        ID_vals  = ID_grid[ IDs ]
        pID_vals = pID_grid[ IDs ]
        pIDs = divmod(pID_vals, d8.nx)
        #-------------------------------------------------
        # Get unique pID_vals? If so, redefine n1 & pIDs
        #-------------------------------------------------
#         pID_vals = np.unique( pID_vals )
#         n1 = pID_vals.size
#         pIDs = divmod(pID_vals) 
              
        #--------------------------------------------------
        # Get start cell elevations and initial distances
        # from start cells to ends of zero-slope chains
        #--------------------------------------------------
        start_zs  = z2[ IDs ]  # z2 is grid, IDs is tuple
        distances = ds[ IDs ]  # ds is grid
        resolved  = np.zeros( n1, dtype='int8' )        
        maxlen    = np.maximum( d8.nx, d8.ny ) * 2
        chan_IDs  = np.zeros( [n1, maxlen], dtype='int32' )
        print('start_zs.size  =', start_zs.size)   # Beaver3 = 17029
        print('distances.size =', distances.size)  # Beaver3 = 17029
        print()
        #----------------------------------------------
        # Only used to check the distance calculation
        #----------------------------------------------
        chan_ds = np.zeros( [n1, maxlen], dtype='float64')

        #------------------------------------------ 
        # Always include the initial cell in path
        # but only change its slope value if
        # CHANGE_START_SLOPE is True
        #------------------------------------------
        chan_IDs[:, 0] = ID_vals
        chan_IDs[:, 1] = pID_vals  # pID_slope is zero
        chan_ds[:,0]   = ds[ IDs ]
        path_index = 2
        
        #--------------        
        # For testing
        #---------------------------------------------------------
        # As a first test of the algorithm, restricted attention
        # to a single flow path that has many cells with a D8
        # slope of zero.  Did this by setting resolved=1 for all
        # other flowpaths.  Used the Beaver3 DEM from MERIT,
        # that has 3-arcsecond grid cell size. A channel profile
        # extracted with RT4 was used to obtain path distances
        # & all other values were checked with RT4 Value Zoom.
        # This test path has the following attributes:
        #
        # Start cell = (252,129), ID = 53142
        # Elevation of start cell = 205.201416 [m]
        # Path dist of start cell = 0.07366452 [km]
        # D8 slope  of start cell = 0.02680158 [m/m]
        #----------
        # First 0-slope cell = (252, 128), ID = 52732
        # Last 0-slope cell  = (246,  98), ID = 40426
        # Elevation of 0-slope cells = 202.722549 [m].
        # D8 parent of last 0-slope cell = (245, 97) = 40015
        #    Elevation of this cell also = 202.722549 [m].
        #    D8 slope of this cell = 0.01215424 > 0.
        #----------
        # End cell = (245, 96), ID = 39695
        # Elevation of end cell = 201.598404 [m]
        # Path dist of end cell = 3.45685792
        #    (obtained from the saved profile)
        # D8 slope  of end cell =  0.002986
        #----------
        # Elev drop between start and end cells = 
        #    (205.201416 - 201.598404) = 3.60301200 [m]
        #    #### Returned now: 3.603012
        # Flow dist between start and end cells =
        #    (3456.85792 - 73.66452) = 3383.1934 [m]
        #    #### Returned now: 3383.197708
        # New slope = drop/dist = 0.00106497370
        #    #### Returned now: 0.0010649723
        #------------
        # Path IDs (35 cells) =
        # [53142 52732 52322 51912 51502 51092 50682 50272
        # 49862 49452 49042 48633 48223 47813 47403 46992
        # 46581 46170 45760 45350 44940 44530 44120 43710
        # 43301 42891 42481 42071 41660 41659 41248 40837
        # 40426, 40015, 39605]
        #------------
        # Total number of cells in chain =  35
        # Number of zero-slope cells in chain = 32
        # Number of same-elevation cells in chain = 33
        #---------------------------------------------------------
        # All zero-slope cells on path are set to the new slope.
        # If (CHANGE_START_SLOPE) then its slope is reset
        #    even though D8 slope of start cell was nonzero.
        # If (CHANGE_NEAR_END_SLOPE) then its slope is reset 
        #    even though D8 slope of this cell was nonzero.        
        #---------------------------------------------------------
        # Explored option of using a different end cell:
        #    End cell = (245, 97), ID = 40015
        #    Elevation of end cell = 202.722549 [m]  ### same as 0-slope cells
        #    Path dist of end cell = 3.36436796 [km]
        #    D8 slope  of end cell =  0.012154  [m/m/]
        #    Elev drop between start and end cells = 
        #       (205.201416 - 202.722549) = 2.47886700 [m]
        #    Flow dist between start and end cells =
        #       (3364.36796 - 73.66452) = 3290.7034400 [m]
        #    New slope = drop/dist = 0.0007532939        
        #---------------------------------------------------------
        TEST1 = False
        if (TEST1):
            w01 = np.where(ID_vals != 53142)
            resolved[ w01 ] = 1
            # Should only be one now; IDs (not pIDs) are unique
            w02 = np.where(ID_vals == 53142)
            print('start_zs for 53142 =', start_zs[w02])
    
    DONE = (nb == 0) or (n1 == 0)
    n_fixed = 0
    if (REPORT):
        print('Working...')
    while not(DONE):
    
        #--------------------------------------        
        # Update pIDs (D8 parents of parents)
        # We're done when pv_max == 0, since
        # pIDs[0] = 0.
        #----------------------------------------
        IDs      = pIDs   ## these tuples have constant length, n1
        pID_vals = pID_grid[ IDs ]   # parents of parents
        ### pID_vals = pID_grid[ pIDs ]  # parents of parents (SAME?)
        pIDs     = divmod( pID_vals, d8.nx)
        pv_max   = pID_vals.max()
        print('pv_max =', pv_max)

        #-------------------------------------------------
        # Note that pID_slope was initially a 2D grid
        # but here is redefined to be a 1D array
        #-------------------------------------------------
        # Here we use the original D8 slope grid so
        # that results are not affected by new slope
        # assignments.  We allow values in slope grid
        # to be overwritten by smaller values.
        # Shorter side channels with zero slopes that
        # enter longer channels will always get resolved
        # first and tend to get assigned larger slopes.
        #-------------------------------------------------
        pID_slope = d8_slope[ pIDs ]
        
        #-----------------------------------------------
        # Note that pID_slope was initially a 2D grid
        # but here is redefined to be a 1D array
        #-----------------------------------------------
        ## pID_slope = slope[ pIDs ]
        
        ## print('### par_slope.shape =', pID_slope.shape )  # is now 1D
        wg = np.where( pID_slope > 0)   # tuple
        wb = np.where( pID_slope == 0)  # tuple
        wn = np.where( pID_slope < 0)   # tuple  (-999)
        ng = wg[0].size
        nb = wb[0].size
        nn = wn[0].size

        #-----------------------------------   
        # Update distances from start cell
        #--------------------------------------------
        # Need (path_index - 1) here because there
        # is one less distance than number of cells
        #--------------------------------------------
        distances += ds[ IDs ]
        chan_ds[:, path_index-1 ] = ds[IDs]  ### For testing on path       

        if (ng > 0):
            #----------------------------------------------
            # Get elevation of the "end cell" & update dz
            #----------------------------------------------
            # Could go outside if stmt, but then slower.
            # Doesn't use wg.
            # gp stands for "grandparent"
            #----------------------------------------------
            gpID_vals = pID_grid[ pIDs ]
            gpIDs     = divmod( gpID_vals, d8.nx )  #####
            end_zs = z2[ gpIDs ]
            dz = (start_zs - end_zs)

            #-----------------------------------------------------
            # Append the parent cell which has slope > 0 but has
            # same elevation as those in path with slope of 0.
            # Don't increase path_index; it affects the chains
            # which are still in progress.
            #-----------------------------------------------------
            chan_IDs[wg, path_index] = pID_vals[ wg ]

            #-----------------------------------------------------
            # Append the "grandparent cell" which has a lower
            # elevation than those in path with slope of 0.
            # It may be possible that it also has a slope of 0.
            #-----------------------------------------------------
            # Increase distance from parent to grandparent, but
            # don't add ds for grandparent, which is distance to
            # its parent.
            #-----------------------------------------------------            
            chan_IDs[wg, path_index+1] = gpID_vals[ wg ]
            #---------------------------------------------
            prev_IDs = (pIDs[0][wg], pIDs[1][wg])
            distances[wg] += ds[ prev_IDs ]
            chan_ds[wg, path_index ] = ds[prev_IDs]
                              
            #------------------------------------------        
            # Can we eliminate this for loop somehow?
            #------------------------------------------
            ## ready_pIDs = pIDs[ wg ]
            ## new_slopes = (dz[wg[0]] / distances[wg[0]])                      
            for j in range(ng):
                w0 = wg[0][j]  # index for this zero-slope chain
                if (resolved[ w0 ] == 0):
                    new_slope = (dz[w0] / distances[w0])
                    if (TEST1):
                        print('dz.size =', dz.size)
                        print('distances.size =', distances.size)
                        print('start_zs[w0] =', start_zs[w0])
                        print('end_zs[w0]   =', end_zs[w0])
                        print('dz[w0]    =', dz[w0])
                        print('distances[w0] =', distances[w0])
                        print('new_slope =', new_slope)
                        print()
    #                 if (REPORT):
    #                     print('new_slope =', new_slope)
 
                    #-------------------------------------------------                   
                    # Change the slopes of some cells on path
                    #-------------------------------------------------
                    # Recall that D8 slope of start cell and 2nd to
                    # last cell were already greater than zero.
                    # Flags determine whether slopes of these cells
                    # are updated to the new value.
                    #-------------------------------------------------
                    CHANGE_START_SLOPE    = True
                    CHANGE_NEAR_END_SLOPE = True
                    if (CHANGE_START_SLOPE):
                        i1 = 0
                    else:
                        i1 = 1
                    #----------------------------
                    if (CHANGE_NEAR_END_SLOPE):
                        i2 = path_index + 1
                    else:
                        i2 = path_index
                    #------------------------------------------------
                    path_ID_vals  = chan_IDs[w0, :path_index+2]
                    chan_ds_vals  = chan_ds[w0,  :path_index+2 ]
                    ready_ID_vals = chan_IDs[w0, i1:i2]
                    #------------------------------------------------
                    ## TEST2 = (53142 in ready_ID_vals)
                    TEST2 = False
                    if (TEST1 or TEST2):
                        print('path_ID_vals.size =', path_ID_vals.size)
                        print('path_ID_vals =', path_ID_vals) 
                        print()
                        print('ready_ID_vals.size =', ready_ID_vals.size)
                        print('ready_ID_vals =', ready_ID_vals)                 
                        # Note: Last value in chan_ds should be 0.
                        print('chan_ds_vals.size =', chan_ds_vals.size )
                        print('chan_ds_vals =', chan_ds_vals)
                        print()
                    #-----------------------------------------
                    ready_IDs = divmod( ready_ID_vals, d8.nx)
                    n_fixed += ready_ID_vals.size
                    #-----------------------------------------------
                    # We are now using original D8 slope grid to
                    # prevent new short channel slope assignments
                    # from interfering with longer channel slope
                    # assignments that tend to have smaller slopes.
                    # Try to streamline this code
                    #------------------------------------------------
                    nz_slopes  = np.maximum(0, slope[ready_IDs])
                    nz_slopes[ nz_slopes <= 0 ] = 99999.0
                    new_slopes = np.minimum(new_slope, nz_slopes)
                    slope[ ready_IDs ] = new_slopes
                    #----------------------------------
                    ## slope[ ready_IDs ] = new_slope

                    #-----------------------------------------                    
                    # Flag this zero-slope chain as resolved
                    #-----------------------------------------
                    # This should now work because all paths
                    # start with a unique ID (not a pID).
                    #-----------------------------------------
                    resolved[ w0 ] = 1

                    #---------------------------------                    
                    # Flag all zero-slope chains with
                    # this same start_ID as resolved
                    #---------------------------------
#                     start_ID = chan_IDs[w0, 0]
#                     w00 = (chan_IDs[:,0] == start_ID)
#                     resolved[ w00 ] = 1

            ww = (slope == 0)
            n_left = ww.sum() 
            ## n_left = (n1 - n_fixed)
            ## print('  n_fixed     =', n_fixed)
            ## print('  path_length =', k)
            print('  n_left      =', n_left)
            print()
        if (nn > 0):
            #-------------------------------------------------
            # These have reached the edge without resolution
            # Makes progress toward (pv_max == 0)
            #-------------------------------------------------
            resolved[ wn ] = 1   # remove path from consideration           
        if (nb > 0):
            #-----------------------------------------------
            # Append cells to each active zero-slope chain
            # pID_slope is zero here
            # pID_vals is not a grid; computed above
            # Use wb[0] vs. wb here ??   ########
            #-----------------------------------------------
            chan_IDs[wb, path_index] = pID_vals[ wb ]
            path_index += 1

        #---------------------------------------
        # Added nreps condition on 2023-08-24.
        # Occurred due to DEM data type issue.
        #---------------------------------------
        n_reps += 1
        if (REPORT):
            print('  n_reps =', n_reps)

        n_unresolved = (1 - resolved).sum()
        DONE = (n_unresolved == 0) or (pv_max == 0) or (n_reps > maxlen)
        ## DONE = (pv_max == 0) or (n_reps > maxlen)
        ## DONE = (nb == 0) or (pv_max == 0) or (n_reps > maxlen)

    #-------------------- 
    # Why are we DONE ?
    #--------------------
    print() 
    if (n_unresolved == 0):
        print('Finished since n_unresolved = 0.')
    if (pv_max == 0):
        print('Finished since pv_max = 0.')
    if (n_reps > d8.nx):
        print('Finished since n_reps > maxlen.')
        print('   maxlen = 2 * min(nrows,ncols).')
        print('   Check the data type of the DEM.')
    print()

    #----------------------------------
    # How many still have zero slope?
    #----------------------------------
    ww = (slope == 0)
    n_left = ww.sum()
    print('Number of cells with valid D8 flow code that')
    print('  still have zero slope =', n_left)
            
    #-------------------------------------------------
    # Was set to -999 at start to avoid if undefined
    #-------------------------------------------------
    w_noflow = (d8.d8_grid == 0)
    n_noflow = w_noflow.sum()
    slope[ w_noflow ] = 0.0
    print('Number of cells with invalid flow code that')
    print('  have zero slope =', n_noflow)
    print()
    #-----------------------------------------------------------
    # Set slope to NaN at all "noflow_IDs" (D8 flow code of 0)
    #-----------------------------------------------------------
    # slope[ d8.noflow_IDs ] = np.nan
    # slope[ d8.noflow_IDs ] = 0.0    ### Can change min or max slope
    
    #--------------------------
    # Print final slope values
    #--------------------------
    Sp_min1 = slope[ slope > 0 ].min()
    Sp_min2 = slope[ slope > Sp_min1 ].min()
    Sp_min3 = slope[ slope > Sp_min2 ].min()
    if (REPORT):
        print('Final min slope = ' + str(slope.min()) )
        print('Final max slope = ' + str(slope.max()) )
        print('Final smallest positive slopes = ')
        print('   ' + str(Sp_min1) )
        print('   ' + str(Sp_min2) )
        print('   ' + str(Sp_min3) )
        print()

    #---------------------------------------
    # Read header_file info for slope_file
    #---------------------------------------
    header_file = (topo_dir + site_prefix + '.rti')
    grid_info = rti_files.read_info( header_file, REPORT=False)
    
    #-----------------------------------------------
    # Save new slope grid to slope_file
    # Don't add topo_dir; may already be full path
    #-----------------------------------------------
    if not(os.sep in slope_file):             
        slope_file2 = (topo_dir + slope_file)
    else:
        slope_file2 = slope_file
    rtg_files.write_grid( slope, slope_file2, grid_info, RTG_type='FLOAT')
    print( 'Finished writing new slope grid to file: ')
    print( slope_file2 )
    print('   n_reps        = ' + str(n_reps) )
    print()
          
#   get_new_slope_grid2()                              
#---------------------------------------------------------------------------
# def get_new_slope_grid0(site_prefix=None, case_prefix=None,
#                        cfg_dir=None, slope_file=None):
# 
#     #---------------------------------------------------------
#     # NB!  This algorithm doesn't work as intended.
#     #---------------------------------------------------------
#     # NB!  The iteration of D8 parent IDs will only terminate
#     #      if we have pID[0]=0.   See d8_global.py.
#     #---------------------------------------------------------
#     # Idea2: If your D8 parent cell (downstream) has S=0,
#     #        then iterate to parent of parent until you
#     #        reach a cell whose D8 parent has S > 0.
#     #        This cell can now be assigned a slope based
#     #        on start z and stop z and distance.
#     #        Connected cells upstream with S=0 should be
#     #        assigned this same slope.
#     #---------------------------------------------------------
#     d8 = get_d8_object(site_prefix=site_prefix,
#                        case_prefix=case_prefix,
#                        cfg_dir=cfg_dir)
# 
#     slope = np.zeros([d8.ny, d8.nx], dtype='float64')
#      
#     z2   = d8.DEM
#     pIDs = d8.parent_IDs
#     z1   = d8.DEM[ pIDs ]  
#     dz   = (z2 - z1)
#     wg   = (dz > 0)   # (array of True or False)
#     ng   = wg.sum()   # (number of cells with True)
#     wb   = np.invert( wg )
#     nb   = wb.sum()  # (number of cells with True)
#     if (ng > 0):
#         slope[ wg ] = (dz[wg] / d8.ds[wg])
# 
#     print('Computing slope grid...')
#     Spmin = slope[ slope > 0 ].min()
#     print('   initial min slope     = ' + str(slope.min()) )
#     print('   initial min pos slope = ' + str(Spmin) )
#     print('   initial max slope     = ' + str(slope.max()) )
#     
#     #-----------------------------------------------------
#     # Get a grid which for each grid cell contains the
#     # calendar-style index (or ID) of the grid cell that
#     # its D8 code says it flows to.
#     #-----------------------------------------------------
#     ### pID_grid = d8.parent_ID_grid.copy()
#     pID_grid = d8.parent_ID_grid
#     ID_grid  = d8.ID_grid
# 
#     #------------------------------------           
#     # Get slopes of the D8 parent cells
#     #------------------------------------
#     S0 = slope.copy()
#     p_slope = slope[ pIDs ]
#     w1 = np.logical_and( slope > 0, p_slope == 0 )
#     start_ID_vals  = ID_grid[ w1 ]          ##################
#     pID_vals       = pID_grid[ w1 ]
#     n1 = start_ID_vals.size
#     print('Number of initial IDs = ' + str(n1) )
#     start_IDs = divmod( start_ID_vals, d8.nx)
#     pIDs      = divmod( pID_vals, d8.nx)
#     start_z   = d8.DEM[ start_IDs ]
#     L         = d8.ds[ start_IDs ]  # (cumulative stream length)
#     #----------------------------------
#     ### max_reps = 500  #######
#     ### n_reps = 0 
#     n_left = n1 
#     if (n1 > 0):
#         # Next line gives no error even if n1==0.
#         unfinished = np.zeros( n1 ) + 1     
#     DONE = (n1 == 0)
#     while not(DONE):
#         #------------------------------------- 
#         # Fist, get parents of pIDs (not IDs)
#         #----------------------------------------------------
#         # Here, IDs & pIDs are tuples, like from np.where,
#         # while ID_vals & pID_vals are long-integer indices.
#         #----------------------------------------------------
#         # Note: Number of unique pIDs <= number of IDs
#         #----------------------------------------------------
#         ## IDs = pIDs.copy()          # (not an option)
#         ## IDs = copy.copy( pIDs )    # (maybe needed)    
#         IDs = pIDs
#         ID_vals = ID_grid[ IDs ]
#         pID_vals = pID_grid[ IDs ]
#         pIDs = divmod( pID_vals, d8.nx )  # (tuple, like WHERE)
#         L  += d8.ds[ IDs ]
# 
#         #---------------------------------------------------
#         # Allow use of recently assigned slopes
#         # This should prevent overwriting resolved values.
#         #---------------------------------------------------
#         k_slope = slope[ IDs ]
#         p_slope = slope[ pIDs ]
#         ## k_slope = S0[ IDs ]
#         ## p_slope = S0[ pIDs ]
#         
#         test1 = np.logical_and( p_slope > 0, k_slope == 0)
#         wg = np.logical_and( test1, unfinished == 1)
#         ## wg = np.logical_and( p_slope > 0, k_slope == 0 )
#         ## wg = np.logical_and( p_slope > 0, unfinished == 1 )
#         ng = wg.sum()
#         if (ng > 0):
#             ready_ID_vals = ID_vals[ wg ]
#             ready_IDs = divmod( ready_ID_vals, d8.nx )
#             dz = (start_z - d8.DEM[ pIDs ])    # (drop from start to parent)
#             slope[ ready_IDs ] = dz[ wg ] / L[ wg ]
#             unfinished[ wg ] = 0
#             n_left = unfinished.sum()
#             print('n_left = ' + str(n_left) )
#         #---------------------------------------------------
#         if (n_left == 0):
#             print('Step 1 finished: n_left = 0.')
#         #---------------------------------------------------
#         pv_max = pID_vals.max()
#         if (pv_max == 0):
#             print('Step 1 finished: Reached edge of DEM.')
#         
#         DONE = (n_left == 0) or (pv_max == 0) 
#                
# #         n_reps += 1
# #         if (n_reps == max_reps):
# #             print('Aborting after ' + str(max_reps) + ' reps.') 
#         ### DONE = (n_reps == max_reps) or (pID_vals.max() == 0)   ########            
#  
#     #--------------------------------------------------------------    
#     # Step 2.  Assign slope to the in-between grid cells with S=0
#     #--------------------------------------------------------------
#     # If D8 parent was assigned a slope in Step 1, then any of
#     # its D8 kids with S=0 should be assigned the same slope.
#     # Iterate until all cells with S=0 are assigned a slope.
#     #--------------------------------------------------------------
#     pIDs = d8.parent_IDs
#     # DONE = (n1 == 0)
#     DONE = False
#     while not(DONE):
#         #----------------------------------------------
#         # pIDs is not changing but slope is changing.
#         # We're working upstream with the whole grid.
#         #----------------------------------------------------
#         # Note: Since 2 grid cells can have same D8 parent,
#         # all upstream cells with S=0 connected to a given
#         # parent (perhaps along different streams) will be
#         # assigned the same slope.  And the parent itself
#         # may have been assigned 2 or more different slopes
#         # in Step 1 (one overwriting the other).
#         #----------------------------------------------------        
#         pslope = slope[ pIDs ]
#         w2 = np.logical_and( slope == 0, pslope > 0 )
#         n2 = w2.sum()
#         if (n2 > 0):
#             slope[ w2 ] = pslope[ w2 ]
#         DONE = (n2 == 0)
# 
#     #-----------------------------------------------------------
#     # Set slope to NaN at all "noflow_IDs" (D8 flow code of 0)
#     #-----------------------------------------------------------
#     slope[ d8.noflow_IDs ] = np.nan
#     
#     #---------------------------------------
#     # Read header_file info for slope_file
#     #---------------------------------------
#     header_file = (cfg_dir + site_prefix + '.rti')
#     grid_info = rti_files.read_info( header_file, REPORT=False)
#     
#     #------------------------------------
#     # Save new slope grid to slope_file
#     #------------------------------------             
#     slope_file2 = (cfg_dir + slope_file)
#     rtg_files.write_grid( slope, slope_file2, grid_info, RTG_type='FLOAT')
#     print( 'Finished writing new slope grid to file: ')
#     print( slope_file2 )
#     print()
#     ### print('Finished computing slope grid.')
#     Spmin = slope[ slope > 0 ].min()
#     print('   n_reps        = ' + str(n_reps) )
#     print('   min slope     = ' + str(np.nanmin(slope)) )
#     print('   min pos slope = ' + str(Spmin) )
#     print('   max slope     = ' + str(np.nanmax(slope)) )
#     print(' ')
#           
# #   get_new_slope_grid0() 
#---------------------------------------------------------------------------
def smooth_dem(site_prefix=None, case_prefix=None,
               cfg_dir=None, dem_file=None,
               nsteps=100, ncells=3):

    #-------------------------------------------------  
    # Idea:  Use 1D FTCS scheme along D8 flow paths 
    #-------------------------------------------------    
    # For numerical stability, FTCS scheme requires:
    # r = (a * dt / dx^2) < 1/2.
    #-------------------------------------------------
    
    if (dem_file is None):
        dem_file = site_prefix + '_smooth_DEM.rtg'

    d8 = get_d8_object(site_prefix=site_prefix,
                       case_prefix=case_prefix,
                       cfg_dir=cfg_dir)
    topo_dir = d8.topo_directory
    
    z = d8.DEM.copy()
    pIDs = d8.parent_IDs
    ## ID_grid  = d8.ID_grid
    pID_grid = d8.parent_ID_grid
    gpID_grid = pID_grid[ pIDs ]
    gpIDs = divmod( gpID_grid, d8.nx )  # (tuple, like np.where)
    #--------------------------------------
#     d82 = copy.copy(d8)
#     slope = np.zeros([d8.ny, d8.nx], dtype='float64')
    
    # a  = 1
    # dt = 10
    # dx = 1800.0 
    # r = (a * dt) / dx**2
    r = 0.0001
    #------------------
#     DONE = False
#     nsteps = 0
#     while not(DONE):
    #------------------
    for k in range(nsteps):
        #-----------------------------------------------------
        # Next line was a bug; values in flats didn't change
        # But the rest of the profile was smoothed.
        #-----------------------------------------------------
        ## z[:] = z[pIDs] + r * (z - 2*z[pIDs] + z[gpIDs])
        #-----------------------------------------------------
        # Note that max value doesn't change with this method
        # because cells that aren't in pIDs don't change.
        #-----------------------------------------------------
        z[pIDs] = (z + z[gpIDs]) / 2.0
        
#         rhs     = z[pIDs] + r * (z - 2*z[pIDs] + z[gpIDs])
#         z[pIDs] = rhs
        w = np.where( np.logical_or( pID_grid == 0, gpID_grid == 0 ) )
        nw = w[0].size
        if (nw > 0):
            z[w] = d8.DEM[ w ]
            
        ## z[ d8.noflow_IDs ] = d8.DEM[ d8.noflow_IDs ]            
        ### z[ d8.edge_IDs ] = d8.DEM[ d8.edge_IDs ]

        #-------------------------------------------------
        # Check if there are still any D8 slopes of zero
        #-------------------------------------------------
        # It is possible that D8 flow directions changed
        #-------------------------------------------------
#         nsteps += 1
#         d82.update(0, DEM=z)
#         pIDs = d82.parent_IDs
#         gpID_vals = pID_grid[ pIDs ]
#         gpIDs = divmod( gpID_vals, d82.nx )  # (tuple, like WHERE)
#         #-----------------------------------------------------------
#         dz = (z - z[ pIDs ])
#         # dz = (z - z[ pIDs2 ])
#         ## wb = (dz == 0)  # (array of True or False)
#         wb = (dz <= 0)  # (array of True or False)
#         nb = wb.sum()   # (number of cells with True)
#         print('nb = ' + str(nb) )
#         DONE = (nb == 0) or (nsteps > 20)
        
#         wg = (dz > 0)  # (array of True or False)
#         wb   = np.invert( wg )
#         nb   = wb.sum()  # (number of cells with True)
#         ds   = d82.ds.copy()
#         if (ng > 0):
#             slope[ wg ] = (dz[wg] / ds[wg])          

    #-------------------------------------
    # Read header_file info for dem_file
    #-------------------------------------
    header_file = (topo_dir + site_prefix + '.rti')
    grid_info = rti_files.read_info( header_file, REPORT=False)

    #---------------------------
    # Save pID_grid to file
    #---------------------------            
#     pID_grid_file = (cfg_dir + site_prefix + '_pID-grid.rtg')
#     rtg_files.write_grid( pID_grid, pID_grid_file, grid_info, RTG_type='LONG')
        
    #---------------------------
    # Save new DEM to dem_file
    #---------------------------            
    dem_file2 = (topo_dir + dem_file)
    rtg_files.write_grid( z, dem_file2, grid_info, RTG_type='FLOAT')
    print( 'Finished writing new DEM to file: ')
    print( dem_file2 )
    print()
    zpmin = z[ z > 0 ].min()
    print('   r         = ' + str(r) )
    print('   nsteps    = ' + str(nsteps) )
    print('   min z     = ' + str(np.nanmin(z)) )
    print('   min pos z = ' + str(zpmin) )
    print('   max z     = ' + str(np.nanmax(z)) )
    print(' ')
                               
#   smooth_dem()    
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



