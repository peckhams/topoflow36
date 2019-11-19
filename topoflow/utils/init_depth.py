#------------------------------------------------------------------------
#  Copyright (c) 2019, Scott D. Peckham
#
#  Sep 2019.  First version to put in topoflow utils folder.
#             Started from IDL version:  init_depth.pro.
#  Nov 2019.  Added get_baseflow_volume_flux().  Edited test2().

#------------------------------------------------------------------------

#  test1()
#  test2()

#  remove_bad_slopes()
#  get_baseflow_volume_flux()   # 2019-11-14
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
def test2():

    site_prefix = 'Baro_Gam_1min'
    cfg_dir = '/Users/peckhams/Desktop/TopoFlow_2019/__Regions/'
    cfg_dir += 'Ethiopia/DEMs/MERIT/Baro_Gam_1min/'
    #--------------------------------------------------------
    area_file    = (site_prefix + '_area.rtg')
    slope_file   = (site_prefix + '_slope.rtg')
    width_file   = (site_prefix + '_chan-w.rtg') 
    manning_file = ( site_prefix + '_chan-n.rtg')
    sinu_file    = (site_prefix + '_sinu.rtg')
    # angle_file = (site_prefix + '_chan-a.rtg')
    #--------------------------------------------------------
    d0_file = (site_prefix + '_d0.rtg')

    #-------------------------------------------------    
    # Estimate bankfull discharge for the Baro River
    #-------------------------------------------------
    A_km2 = np.float64(23567.7)  # [km^2]
    A_m2  = A_km2 * 1e6          # [km^2 -> m^2]
    u_bankfull = 1.5        # [m s-1]
    w_bankfull = 140.0      # [m]
    d_bankfull = 6.0        # [m]
    Q_bankfull = u_bankfull * w_bankfull * d_bankfull  # [m3 s-1]   

    #----------------------------------------------    
    # Estimate baseflow discharge and volume flux
    # See Notes to get_baseflow_volume_flux().
    #----------------------------------------------
    Q_baseflow = Q_bankfull / 8.0 
    B_mps  = (Q_baseflow / A_m2)       # [m s-1]
    B_mmph = B_mps * 1000.0 * 3600.0   # [mm h-1] 
    #----------------------------------------------------------
    print('Baseflow estimate = ' + str(B_mps)  + ' [m s-1]' )
    print('Baseflow estimate = ' + str(B_mmph) + ' [mm h-1]' )

    #----------------------------------------------     
    # Compute the initial channel flow depth grid
    #----------------------------------------------                                          
    compute_initial_depth( site_prefix=site_prefix, cfg_dir=cfg_dir,
            baseflow_rate=B_mps, tol=None, SILENT=False,
            #--------------------------------------------
            area_file=area_file, slope_file=slope_file,
            width_file=width_file, manning_file=manning_file,
            sinu_file=sinu_file, angle_file=None,
            d0_file=d0_file)
 
#   test2()                          
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
        print('         min(slope) = ' + str(S_min))
        print('         max(slope) = ' + str(S_max))
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
def get_baseflow_volume_flux( A_out_km2, Qb_out, REPORT=True,
                              MMPH=False):

    #-----------------------------------------------------------------
    # Notes: A_out = total contributing area (TCA) at outlet [km^2].
    #        Qb_out = baseflow discharge at outlet [m3 s-1]
    #        This could be the minimim value recorded at a gauge.
    #------------------------------------------------------------------
    # Q_out = (u_out * w_out * d_out)   (baseflow or bankfull)
    #
    # As a rough estimate, we could assume that for baseflow Q:
    #    (1) u_baseflow = 1.0  [m s-1]  (roughly u_bankfull / 2)
    #    (2) w_baseflow = (w_bankfull / 2)
    #    (3) d_baseflow = (d_bankfull / 2)  # assume 45 deg bank angle.
    #
    # Expect u_bankfull to be closer to 2 or 3 [m s-1]. (not bigger)
    # With these assumptions, Q_baseflow is about Q_bankfull / 8.
    #------------------------------------------------------------------
    A_out_m2 = A_out_km2 * 1e6           # convert [km^2] to [m^2]
    B_mps    = (Qb_out / A_out_m2)       # [m s-1]
    B_mmph   = B_mps * 1000.0 * 3600.0   # convert [m s-1] to [mm h-1] 

    #------------------
    # Optional report
    #------------------
    if (REPORT):
        print('Baseflow volume flux = ' + str(B_mps)  + ' [m s-1]' )
        print('Baseflow volume flux = ' + str(B_mmph) + ' [mm h-1]' )

    if (MMPH):
        return B_mmph
    else:
        return B_mps

#   get_baseflow_volume_flux()
#------------------------------------------------------------------------
def compute_initial_depth( site_prefix=None, cfg_dir=None,
                           baseflow_rate=None, bank_angle=None,
                           tol=None, SILENT=False,
                           #-----------------------------------
                           area_file=None, slope_file=None,
                           width_file=None, manning_file=None,
                           sinu_file=None, angle_file=None,
                           d0_file=None, v0_file=None):

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
    #           B = spatially-uniform baseflow volume flux [m s-1]
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
    B = baseflow_rate    # baseflow volume flux [m s-1]
    # print('B = ' + str(B) )
    
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
    grid_info   = rti_files.read_info( header_file, REPORT=False)
    byte_order  = grid_info.byte_order
    ncols       = grid_info.ncols
    nrows       = grid_info.nrows
    
    #----------------------------
    # Read the input grid files
    #----------------------------
    A_file = cfg_dir + area_file
    S_file = cfg_dir + slope_file
    w_file = cfg_dir + width_file
    n_file = cfg_dir + manning_file
    #---------------------------------------------------------------
    # Note: S, w and n may have NaNs on edges.  A has zeros.
    #---------------------------------------------------------------    
    A = rtg_files.read_grid( A_file, grid_info, RTG_type='FLOAT' )
    S = rtg_files.read_grid( S_file, grid_info, RTG_type='FLOAT' )
    w = rtg_files.read_grid( w_file, grid_info, RTG_type='FLOAT' )
    n = rtg_files.read_grid( n_file, grid_info, RTG_type='FLOAT' )
    #------------------------------------------------------------------------   
    if (angle_file is not None):
        q_file = cfg_dir + angle_file    # q = theta = bank_angle
        theta = rtg_files.read_grid( q_file, grid_info, RTG_type='FLOAT' )
    else:
        if (bank_angle is None):
            theta = 30.0   # [degrees]
        else:
            theta = bank_angle
    #------------------------------------------------------------------------  
    if (sinu_file is not None):
        sinu_file = cfg_dir + sinu_file  
        sinu = rtg_files.read_grid( sinu_file, grid_info, RTG_type='FLOAT' )
        sinu[ sinu == 0 ] = 1.0  # (starts with zeros on edges)
    else:
        sinu = 1.0
             
    #---------------------------------
    # Convert A from [km^2] to [m^2]
    #---------------------------------
    A = A * np.float64(1000000)    # [m^2]

    #------------------------------------
    # Convert theta from [deg] to [rad]
    #------------------------------------
    theta = theta * (np.pi / np.float64(180))    # [radians]

    #--------------------------------------
    # Remove zero slopes from slope grid
    # or else we'll get NaNs for them.
    #--------------------------------------
    # Better to use profile-smoothed DEM,
    # then zero slopes should be gone.
    #--------------------------------------
    S = remove_bad_slopes( S )

    #---------------------------------    
    # Adjust slopes by the sinuosity
    #---------------------------------
    S = (S / sinu)
    
    #------------------------------------------------
    # Option 1:
    # Initialize to very large, impossible depth so
    # that we can get largest root (and smallest
    # velocity) if there are multiple roots.
    #------------------------------------------------
    # Initialize d, to be saved to d0_file
    #------------------------------------------------    
    d1 = np.float64( 500 )    # [meters]
    d  = np.zeros( [nrows, ncols], dtype='float64') + d1

    #------------------------------------------------
    # Option 2:
    # Initialize depth grid, d, for Newton iteration
    # by assuming the case where d << w and R = d,
    # where R = Ac/P = hydraulic radius.
    #-------------------------------------------------
    ## p1 = np.float64(3) / 5
    ## d = K / w**p1
        
    #-----------------------------------------
    # Construct the constant K, where:
    #    F(d) = Ac^(5/3) - [K * P^(2/3)] = 0
    #-----------------------------------------
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
    print('Iterating...')
    DONE = False
    while not(DONE):
        last_d = d.copy()
        Ac     = d * (w + (d * np.tan(theta)))
        P      = w + (2 * d / np.cos(theta))
        numer  = (Ac**pow1) - (K * P**pow2)
        term1  = pow1 * (Ac**pow2) * (w + (2*d*np.tan(theta)))
        term2  = pow2 * K * (P**pow3) * (2 / np.cos(theta))
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

        #--------------------
        # Are we done yet ?
        #--------------------
        n_tries += 1
        wb = (gap > tol)  # (array of True or False)
        nb = wb.sum()
        ## wb = np.where( gap > tol )
        ## nb = np.size( wb[0] )
        DONE = (nb == 0) or (n_tries > max_tries)
        
        #------------------------------
        # Write status to message box
        #------------------------------
        if not(SILENT):
            mstr = 'Pixels left = ' + str(nb)
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
    rtg_files.write_grid( d, d0_file2, grid_info, RTG_type='FLOAT')
    print( 'Finished writing file: ')
    print( d0_file2 )

    #---------------------------------    
    # Print min and max values of d0
    #---------------------------------
    if not(SILENT):
        print( 'd_min = ' + str(np.nanmin(d)) + ' [m]')
        print( 'd_max = ' + str(np.nanmax(d)) + ' [m]')

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