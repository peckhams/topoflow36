#------------------------------------------------------------------------
#  Copyright (c) 2019, Scott D. Peckham
#
#  Sep 2019.  First version to put in topoflow utils folder.
#             Started from IDL version:  init_depth.pro.
#------------------------------------------------------------------------

#  test1()

#  remove_bad_slopes()
#  compute_initial_depth()

#------------------------------------------------------------------------

import numpy as np

from . import rtg_files   # (also in utils directory)
from . import rti_files   # (also in utils directory)

# from topoflow.utils import rtg_files
# from topoflow.utils import rti_files

#------------------------------------------------------------------------
def test1():

    compute_initial_depth()

#   test1()
#------------------------------------------------------------------------
def remove_bad_slopes(slope, FLOAT=False, SILENT=False):

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
    #------------------------------------------------------------
    S = slope

    #-----------------------------------
    # Are there any "bad" pixels ?
    # If not, return with no messages.
    #-----------------------------------  
    wb = np.where(np.logical_or((S <= 0.0), \
                          np.logical_not(np.isfinite(S))))
    nbad = np.size(wb[0])
    if not(SILENT):
        print('size(slope) = ' + str(np.size(S)) )
        print('size(wb) = ' + str(nbad) )
    
    wg = np.where(np.invert(np.logical_or((S <= 0.0), \
                                 np.logical_not(np.isfinite(S)))))
    ngood = np.size(wg[0])
    if (nbad == 0) or (ngood == 0):
        return
    
    #---------------------------------------------
    # Find smallest positive value in slope grid
    # and replace the "bad" values with smin.
    #---------------------------------------------
    S_min   = S[ wg ].min()
    S_max   = S[ wg ].max()
    S[ wb ] = S_min  
    if not(SILENT):
        print('-------------------------------------------------')
        print('WARNING: Zero or negative slopes found.')
        print('         Replacing them with smallest slope.')
        print('         Use "Profile smoothing tool" instead.')
        print('         min(S) = ' + str(S_min))
        print('         max(S) = ' + str(S_max))
        print('-------------------------------------------------')
        print(' ')

    #--------------------------------
    # Convert data type to double ?
    #--------------------------------
    if (FLOAT):    
        new_slope = np.float32( S )
    else:    
        new_slope = np.float64( S )
   
    return new_slope
     
#   remove_bad_slopes()
#------------------------------------------------------------------------
def compute_initial_depth( site_prefix=None, cfg_dir=None,
                           baseflow_rate=None, tol=None,
                           SILENT=False,
                           #-----------------------------------
                           area_file=None, slope_file=None,
                           width_file=None, angle_file=None,
                           manning_file=None, d0_file=None):

    #------------------------------------------------------------
    # Note:  This routine uses a grid-based Newton-Raphson
    #        iterative scheme to solve an equation for the
    #        initial depth of water in a channel network
    #        that results from groundwater baseflow.  The
    #        variables involved are:
    #           w = bed bottom width, trapezoid [meters]
    #           A = upstream area [km^2]
    #           S = downstream slope [none]
    #           n = Manning roughness parameter  [s/m^(1/3)]
    #           theta = bank angle [degrees ??]
    #           d = water depth in channel
    #           Ac = wetted cross-section area
    #           P  = wetted cross-section perimeter
    #           Rh = (Ac/P) = hydraulic radius
    #           B = spatially-uniform baseflow volume flux
    #------------------------------------------------------------
    #        The equations used here are:
    #           Q  = v * Ac = B * A    [m3 s-1] (steady-state)
    #           v  = (1/n) R_h^(2/3) * S^(1/2)  [SI units]
    #           Rh = Ac / P [m]
    #           Ac = d * (w + (d * tan(theta)))
    #           P  = w + (2 * d / cos(theta))
    #        If we are given w, n, theta, A, S and B, then
    #        we get an equation for d that cannot be solved
    #        in closed form.  However, we can write the equation
    #        v * Ac = B * A in the form needed to solve for d
    #        (in every grid cell) by Newton's method, i.e.:
    #           F(d) = [v(d) * Ac(d)] - (B * A) = 0
    #------------------------------------------------------------    
    #        Newton's method solves an equation F(d)=0, by
    #        iterating the equation:
    #             d(j+1) = d(j) - [F(d)/F'(d)]
    #        until the difference between d(j+1) and d(j)
    #        drops below the specified tolerance for all
    #        pixels in the grid.  This usually takes about
    #        10 iterations.  Max # of iterations is max_tries.
    #------------------------------------------------------------
    #        Note that NaNs, as can occur on the edges of the
    #        slope grid, are preserved in the depth grid.
    #------------------------------------------------------------
    # NB!    Are multiple roots possible?  If so, does an
    #        alternate root correspond to case of a hydraulic
    #        jump?  Once we have solved for d, we can easily
    #        compute velocity, u, and then use it to select
    #        the appropriate root.
    #------------------------------------------------------------                           
    B = baseflow_rate
    #print('B = ' + str(B) )
    
    #--------------------------------
    # Set tolerance for convergence
    #--------------------------------
    if (tol is None):
        tol = np.float64( 0.001 )

    #----------------------------------------------
    # Read header_file info for all of the files:
    # inputs: area_file, slope_file, width_file,
    #         angle_file, manning_file
    # output:  d0_file (and optional v0_file)
    #----------------------------------------------
    header_file = (cfg_dir + site_prefix + '.rti')
    grid_info   = rti_files.read_info( header_file, REPORT=True)
    byte_order  = grid_info.byte_order
    ncols       = grid_info.ncols
    nrows       = grid_info.nrows
    
    #----------------------------
    # Read the input grid files
    #----------------------------
    A_file = cfg_dir + area_file
    S_file = cfg_dir + slope_file
    w_file = cfg_dir + width_file
    q_file = cfg_dir + angle_file    # q = theta = bank_angle
    n_file = cfg_dir + manning_file
    #---------------------------------------------------------------
    A = rtg_files.read_grid( A_file, grid_info, RTG_type='FLOAT' )
    S = rtg_files.read_grid( S_file, grid_info, RTG_type='FLOAT' )
    w = rtg_files.read_grid( w_file, grid_info, RTG_type='FLOAT' )
    q = rtg_files.read_grid( q_file, grid_info, RTG_type='FLOAT' )
    n = rtg_files.read_grid( n_file, grid_info, RTG_type='FLOAT' )

    #---------------------------------
    # Convert A from [km^2] to [m^2]
    #---------------------------------
    A = A * np.float64(1000000)    # [m^2]

    #------------------------------------
    # Convert theta from [deg] to [rad]
    #------------------------------------
    q = q * (np.pi / np.float64(180))    # [radians]

    #--------------------------------------
    # Remove zero slopes from slope grid
    # or else we'll get NaNs for them.
    #--------------------------------------
    # Better to use profile-smoothed DEM,
    # then zero slopes should be gone.
    #--------------------------------------
    S = remove_bad_slopes( S )
    
    ;------------------------------------------------
    # Option 1:
    # Initialize to very large depth so that we can
    # get largest root (and smallest velocity) if
    # there are multiple roots.
    ;------------------------------------------------
    d0 = np.float64( 500 )    # [meters]
    d  = np.zeros( [nrows, ncols], dtype='float64') + d0

    #------------------------------------------------
    # Option 2:
    # Initialize depth grid, d, for Newton iteration
    # by assuming the case where d << w and R = d,
    # where R = Ac/P = hydraulic radius.
    ;-------------------------------------------------
    ## p1 = np.float64(3) / 5
    ## d = K / w**p1
        
    #-----------------------------------------
    # Construct the constant K, where:
    #    F(d) = Ac^(5/3) - [K * P^(2/3)] = 0
    ;-----------------------------------------
    K = (B * A * n) / np.sqrt(S)

    #------------------
    # Initialize vars
    #------------------
    pow1 = np.float64(5) / 3
    pow2 = np.float64(2) / 3
    pow3 = np.float64(1) / -3
    #--------------------------
    n_tries   = np.int32(0)
    max_tries = np.int32(40)

    #----------------------------------------
    # Iterate entire grid to get d, using a
    # grid-based Newton root-finding method
    #----------------------------------------
    DONE = False
    while not(DONE):
        last_d = d.copy()
        Ac     = d * (w + (d * np.tan(theta)))
        P      = w + (2 * d / np.cos(theta))
        numer  = Ac**pow1 - (K * P**pow2)
        term1  = pow1 * Ac^pow2 * (w + (2*d*np.tan(theta)))
        term2  = pow2 * K * P**pow3 * (2 / np.cos(theta))
        denom  = term1 - term2
        d      = d - (numer/denom)

        #--------------
        # For debugging
        #--------------
        ### print('min(denom) = ' + str(denom.min()) )

        #-------------------------------
        # Compute difference from goal
        # Note: gap = abs(numer/denom)
        #-------------------------------
        gap = np.abs(d - last_d)
        ### gmin = gap.min()  ### HOW TO EXCLUDE NANs?
        ### gmax = gap.max()  ### HOW TO EXCLUDE NANs?
        ### print('     gmin, gmax = ' + str(gmin) + ', ' + str(gmax) )

        #-----------==-------
        # Are we done yet ?
        #--------------------
        n_tries += 1
        wb = np.where( gap > tol )
        nb = np.size( wb[0] )
        DONE = (nb == 0) or (n_tries > max_tries)
        
        #------------------------------
        # Write status to message box
        #------------------------------
        if not(SILENT):
            mstr = 'Pixels left = ' + str(nwg)
            print( mstr )
    
        #----------------------------------
        # Make sure depth is reasonable ?
        #----------------------------------
#         dmin = d.min()   ### HOW TO EXCLUDE NANS?
#         dmax = d.max()   ### HOW TO EXCLUDE NANS?
#         if (dmax > 100.0) or (dmin < 0):
#             DONE = True
#             print('dmin = ' + str(dmin))
#             print('dmax = ' + str(dmax))

    #--------------------------------
    # Failure-to-converge message ?
    #--------------------------------
    if (n_tries > max_tries):
        #------------------------------
        # Write status to message box
        #------------------------------
        if not(SILENT):
            nstr = str(max_tries)
            mstr = 'No convergence after ' + nstr + ' tries'
            print( mstr )

    #----------------------------
    # Save result to depth_file
    #------------------------------------------------
    # Byte swapping may occur in write_grid() based
    # on byte order recorded in grid_info object.
    #------------------------------------------------
    d0_file2 = (cfg_dir + d0_file)
    rtg_files.write_grid( d0, d0_file2, grid_info, RTG_type='FLOAT')
    print( 'Finished writing file: ')
    print( d0_file2 )
    
    #------------------------------------
    # Compute the initial velocity grid
    # using the Manning formula ?
    #------------------------------------
    if (v0_file is not None):
        v0 = (Ac / P)**pow2 * np.sqrt(S) / n

        #---------------------------------
        # Save the initial velocity grid
        #------------------------------------------------
        # Byte swapping may occur in write_grid() based
        # on byte order recorded in grid_info object.
        #------------------------------------------------
        v0_file2 = (cfg_dir + v0_file)
        rtg_files.write_grid( v0, v0_file2, grid_info, RTG_type='FLOAT')
        print( 'Finished writing file: ')
        print( v0_file2 )

    print()

#   compute_initial_depth()                         
#------------------------------------------------------------------------