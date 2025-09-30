#
#  Copyright (c) 2024, Scott D. Peckham
#
#  Mar 2024.  3/22.  Wrote read_dem, set_params,
#             get_d_dx_R, get_d_dx_L, get_d_dy_T, get_d_dy_B,
#             get_qp_gradient, get_qp_slope
#             get_qp_slope_power, get_qp_p_laplacian,
#             get_div_RT, get_div_RB, get_div_LT, get_div_LB.
#
#             3/23.  Wrote get_profile_file, read_long_profile,
#             get_profile_curve, plot_long_profile,
#             get_dem_path, get_histogram_mode,
#             Obtained 1m DEM for Beaver Creek, w/ FLOAT type.
#
#             3/24.  Processed 1m DEM with RiverTools 4.0.
#             Beaver Creek 1m DEM: 27531 cols, 42159 rows.
#             Extracted main channel profile for 1m DEM.
#             Added curve-fitting to plot_long_profile.
#             Wrote get_div_RT, get_div_RB, get_div_LT, and
#             get_div_RB.
#             Wrote get_peckham_profile, and
#             get_whipple_profile.  753 lines of code now.
#
#             3/25. Wrote "profile class" for profile work.
#             Many improvements throughout.
#             Added clip_hilltop to read_profile, plot_profile().
#             Wrote get_narrowest_histogram, get_value.
#             Wrote get_pde_solution and tested.
#  
#             3/26. Wrote save_var_to_rtg.
#             3/27. Downloaded & analyzed Factory Butte DEM.
#             3/28. Wrote fit_curvature_eqn, curvature_eqn.
#                   Wrote get_mask, smooth_dem. 1159 code lines.
#             4/1.  Wrote get_slope_profile, get_kapp_p_profile.
#
#-------------------------------------------------------------------
#  Note: Here, qp stands for "quarter pixel".  If every grid cell
#        or pixel, is divided into 4 parts, then the gradient
#        vector (and continuous angle direction of flow) for each
#        quarter pixel is well-defined because it is determined
#        by exactly 3 elevations that uniquely define a plane.
#        In a river basin, using quarter pixels within grid cells
#        on peaks and ridges also allows rain falling on the grid
#        cell to flow away in multiple directions.  Note that
#        flow direction is ill-defined at such points.
#
#  Note: The ideal, steady-state landscape equation is derived
#        from the fact that the divergence of the 2D vector field
#        of unit-width (depth-integrated) discharge must equal
#        the spatially uniform addition of mass via a uniform
#        rainrate.  That is, for any 2D control volume, there is
#        a net positive amount of fluid leaving that control volume,
#        namely the amount added by rainfall.  A minus sign appears
#        since water flows in the *opposite* direction of the land
#        surface gradient, but this sign is moved to the right-hand
#        side.
#
#  Note: When we compute the PDE operator for a DEM, sometimes the
#        result is positive.  This may be occurring near the tops
#        of peaks and ridges where the slopes are also small (as
#        they are further downstream) and the magnitude of the
#        gradient vector is therefore not actually given by
#        S^gamma, with gamma < 0.  We may need to exclude any of
#        these "hillslope pixels" where the slope decreases in
#        the upstream direction.  Could make a mask of grid cells
#        where computed PDE operator is positive and see if they
#        match up to hillslopes.
#
#  Note: Beside the main channel profile, are other long profiles
#        best-fit by the same value of gamma and Rs?
#
#  Note: Use get_pde_solution() to test this on a known, closed-
#        form, radial solution to the steady-state landform PDE.
#
#  Note: Test on smaller, clean areas in 1m DEM.
#        May need to apply smoothing.
#
#-------------------------------------------------------------------
#  class hamiltonian()
#      __init__()
#      read_dem()
#      get_pde_solution()
#      get_qp_grid()
#      get_d_dx_R()    # "forward"  x-derivative
#      get_d_dx_L()    # "backward" x-derivative
#      get_d_dy_T()    # "forward"  y-derivative
#      get_d_dy_B()    # "backward" y-derivative
#      get_qp_gradient()
#      get_qp_slope()
#      get_qp_slope_power()
#      get_qp_pde_lhs()
#      get_curvatures_RT()
#      fit_curvature_eqn()
#
#      get_value()
#      save_var_to_rtg()
#      get_div_RT()
#      get_div_RB()
#      get_div_LT()
#      get_div_LB()
#      get_histogram_mode()
#
#  class profile()
#      __init__()
#      get_profile_file()
#      read_profile()
#      smooth_profile()
#      clip_profile()
#      get_peckham_profile()
#      get_whipple_profile()
#      get_slope_profile()
#      get_kappa_p_profile()
#      plot_profile()
#
#  kp_from_curvature_eqn()
#  kc_from_curvature_eqn()
#
#-------------------------------------------------------------------
#
#  >>> from topoflow.utils import hamiltonian as ham
#  >>> h = ham.hamiltonian()
#  >>> h.get_qp_pde_lhs()
#  >>> h.get_histogram_mode(bins=80, qp='RT', PLOT=True)
#
#  >>> p = ham.profile()
#  >>> p.plot_profile()
#  >>> print(p.S[0:25])  # (S = slope)
#  >>> p.plot_profile( clip_hilltop=True )
#  >>> p.dataset = 'Beaver_Creek_1m'
#  >>> p.plot_profile()
#  >>> p.plot_profile( clip_hilltop=True )
#
#-------------------------------------------------------------------
 
from topoflow.utils import import_grid as ig
from topoflow.utils import rti_files as rti_files

import numpy as np
import scipy.optimize as spo   # (for curve_fit)
import scipy.signal as sps     # (for convolve2d)
import matplotlib.pyplot as plt

# import matplotlib.dates as mdates
# from matplotlib import cm
# from matplotlib.colors import ListedColormap

#-------------------------------------------------------------------
class hamiltonian():

    #---------------------------------------------------------------
    def __init__( self ):

        #--------------------------    
        # Set the default dataset
        #--------------------------
        self.dataset = 'Factory_Butte_p1'
        ## self.dataset = 'Beaver_Creek_p1'
        ## self.dataset = 'Beaver_Creek_UTM'
        
        # See Peckham papers.
        # But these params are for Beaver_Creek_Geo.
        self.gamma  = -0.701
        self.Rs     = 0.0035
        self.nodata = np.nan
        self.S_bad  = np.nan
        self.SILENT = False

        ### self.nodata = -1
        ### self.S_bad = 1.0  ###############

    #   __init__()
    #---------------------------------------------------------------
    def get_dem_path( self ):

        if (self.dataset == 'Beaver_Creek_UTM'):
            rtg_dir = '/Applications/RiverTools_4.0/data/Beaver_Creek_UTM/'
            ### rtg_file = 'Beaver_DEM.rtg'
            rtg_file = 'Beaver_rawDEM.rtg'
        elif (self.dataset == 'Beaver_Creek_1m'):
            rtg_dir = '/Users/peckhams/DEMs/Beaver_Creek_1m/'
            rtg_file = 'Beaver_Creek_1m_rawDEM.rtg'
        elif (self.dataset == 'Beaver_Creek_p1'):
            rtg_dir = '/Users/peckhams/DEMs/Beaver_Creek_1m/'
            rtg_file = 'Patch1_DEM.rtg'
        elif (self.dataset == 'Beaver_Creek_p1a'):
            rtg_dir = '/Users/peckhams/DEMs/Beaver_Creek_1m/'
            rtg_file = 'Patch1a_DEM.rtg'
        elif (self.dataset == 'BC_small_basin2'):
            rtg_dir = '/Users/peckhams/DEMs/Beaver_Creek_1m/'
            rtg_file = 'BC_small_basin2_rawDEM.rtg'
        elif (self.dataset == 'BC_small_basin2b'):
            rtg_dir = '/Users/peckhams/DEMs/Beaver_Creek_1m/'
            rtg_file = 'BC_small_basin2b_rawDEM.rtg'
        elif (self.dataset == 'Factory_Butte_p1'):
            rtg_dir = '/Users/peckhams/DEMs/Factory_Butte_Utah_1m/'
            rtg_file = 'Factory_Butte_p1_rawDEM.rtg'
        elif (self.dataset == 'Factory_Butte_Loc1a'):
            rtg_dir = '/Users/peckhams/DEMs/Factory_Butte_Utah_1m/'
            rtg_file = 'Factory_Butte_Loc1a_rawDEM.rtg'
        else:
            rtg_dir  = ''
            rtg_file = ''

        self.rtg_dir  = rtg_dir
        self.dem_path = rtg_dir + rtg_file
    
    #   get_dem_path()
    #---------------------------------------------------------------
    def read_dem(self, REPORT=True):

        print('Reading DEM...')
        self.get_dem_path()
           
        self.dem = ig.read_from_rtg( self.dem_path, REPORT=False)

        print('self.dem_path =', self.dem_path)

        self.rti = rti_files.read_info(self.dem_path, SILENT=False, REPORT=REPORT)
        self.nx  = self.rti.ncols
        self.ny  = self.rti.nrows
        self.dx  = self.rti.xres
        self.dy  = self.rti.yres

        ## nx, ny = self.dem.shape

        dtype = str( self.dem.dtype )
        if (dtype.startswith('int')):
            print('### WARNING: This DEM has INTEGER type (not FLOAT).')
            print('###          Changing type to float64.')
            print()
            self.dem = np.float64( self.dem )
            self.rti.data_type = 'DOUBLE'

    #   read_dem()
    #---------------------------------------------------------------
    def smooth_dem(self, REPLACE_DEM=False, n_passes=10):

        #------------------------------------    
        # Gaussian 5x5 kernel approximation
        #------------------------------------
        a = [[1,4,6,4,1], [4,16,24,16,4], [6,24,36,24,6],
             [4,16,24,16,4], [1,4,6,4,1]]
        kernel = np.array(a, dtype='float64') / 256.0
        dem0 = self.dem.copy()
        for k in range(n_passes):
            #---------------------------------------------------                   
            # Note: dem2 will have 4 extra rows & 4 extra cols
            #---------------------------------------------------
            dem2 = sps.convolve2d( dem0, kernel,
                       mode='full', boundary='fill', fillvalue=0)
            dem2 = dem2[2:-2, 2:-2]
            dem0 = dem2.copy()

        print('dem.shape  =', self.dem.shape)
        print('dem2.shape =', dem2.shape)
        self.smooth_DEM = dem2
        self.save_var_to_rtg('smooth_DEM')
        if (REPLACE_DEM):
            self.dem = dem2

    #   smooth_dem()
    #---------------------------------------------------------------
    def get_pde_solution(self, name='radial_log'):
    
        self.dataset = name
        if (name == 'radial_log'):
            # Solution for gamma = -1
            self.gamma = -1.0
            Rs = 0.0035
            a  = 0.0
            b  = 1.0
            nx = 200
            ny = 200
            #---------------
            xmin = -5.0
            xmax = 5.0
            ymin = -5.0
            ymax = 5.0
            x  = np.linspace(xmin, xmax, nx)
            y  = np.linspace(ymin, ymax, ny)
            xx, yy = np.meshgrid(x, y)
            if (Rs > 0):
                z = (-1/Rs)*np.log( a + (xx**2 + yy**2)/b )
            else:
                z = (xx**2 + yy**2)/b

            self.nx  = nx
            self.ny  = ny
            self.dx  = (xmax - xmin)/nx
            self.dy  = (ymax - ymin)/ny
            self.dem = z

    #   get_pde_solution()
    #---------------------------------------------------------------
    def get_qp_grid(self):

        #---------------------------------------    
        # Get empty "quarter pixel" grid, with
        # nx and ny twice as big as the DEM.
        #---------------------------------------
        nx2 = self.nx * 2
        ny2 = self.ny * 2
        qp_grid = np.zeros( (ny2, nx2), dtype='float64' )
        return qp_grid

    #   get_qp_grid()
    #---------------------------------------------------------------
    def get_d_dx_R(self, f):

        #-------------------------------------------------------
        # Remember in 2D numpy arrays, y-axis=0, x-axis=1.
        # Also, shift=-1 means cells with higher x or y-value.
        #-------------------------------------------------------

        #-----------------------------------------------
        # NOTE!  z must have same units as dx and dy!!
        #-----------------------------------------------  
        f_R = np.roll( f, shift=-1, axis=1 )
        f_R[:,-1] = self.nodata
        df_dx_R = (f_R - f) / self.dx
        return df_dx_R

    #   get_d_dx_R()
    #---------------------------------------------------------------
    def get_d_dx_L(self, f):

        #-----------------------------------------------
        # NOTE!  z must have same units as dx and dy!!
        #----------------------------------------------- 
        f_L = np.roll( f, shift=1, axis=1 )
        f_L[:,0] = self.nodata
        df_dx_L = -1 * (f_L - f) / self.dx   #### CHECK -1
        return df_dx_L

    #   get_d_dx_L()
    #---------------------------------------------------------------       
    def get_d_dy_T(self, f):

        #-----------------------------------------------
        # NOTE!  z must have same units as dx and dy!!
        #-----------------------------------------------
        f_T = np.roll( f, shift=1, axis=0 )
        f_T[0,:] = self.nodata
        df_dy_T = (f_T - f) / self.dy
        return df_dy_T

    #   get_d_dy_T()
    #---------------------------------------------------------------
    def get_d_dy_B(self, f):

        #-----------------------------------------------
        # NOTE!  z must have same units as dx and dy!!
        #-----------------------------------------------
        f_B = np.roll( f, shift=-1, axis=0 )
        f_B[-1,:] = self.nodata
        df_dy_B = -1 * (f_B - f) / self.dy    ### CHECK -1
        return df_dy_B

    #   get_d_dy_B()
    #---------------------------------------------------------------
    def get_qp_gradients(self):
      
        print('Computing qp gradients...')
        self.dz_dx_R = self.get_d_dx_R( self.dem )
        self.dz_dx_L = self.get_d_dx_L( self.dem )
        self.dz_dy_T = self.get_d_dy_T( self.dem )
        self.dz_dy_B = self.get_d_dy_B( self.dem )
        
    #   get_qp_gradients()
    #---------------------------------------------------------------
    def get_qp_slopes(self):

        print('Computing slopes...')
        self.S_RT = np.sqrt((self.dz_dx_R ** 2.) + (self.dz_dy_T ** 2.))
        self.S_RB = np.sqrt((self.dz_dx_R ** 2.) + (self.dz_dy_B ** 2.))
        self.S_LT = np.sqrt((self.dz_dx_L ** 2.) + (self.dz_dy_T ** 2.))
        self.S_LB = np.sqrt((self.dz_dx_L ** 2.) + (self.dz_dy_B ** 2.))
        
        S_RT_min = np.nanmin( self.S_RT )
        S_RT_max = np.nanmax( self.S_RT )
        S_RB_min = np.nanmin( self.S_RB )
        S_RB_max = np.nanmax( self.S_RB )
        S_LT_min = np.nanmin( self.S_LT )
        S_LT_max = np.nanmax( self.S_LT )
        S_LB_min = np.nanmin( self.S_LB )
        S_LB_max = np.nanmax( self.S_LB )
        
        if not(self.SILENT):
            print('min, max(S_RT) =', S_RT_min, S_RT_max)
            print('min, max(S_RB) =', S_RB_min, S_RB_max)
            print('min, max(S_LT) =', S_LT_min, S_LT_max)
            print('min, max(S_LB) =', S_LB_min, S_LB_max)
            print()

        #---------------------------------------        
        # Must avoid 0^gamma, since gamma < 1.
        #---------------------------------------
        S_str = str(self.S_bad)
        print('### WARNING:  Changing slope=0 to slope=' + S_str + '.')
        #----------------------------
        w1 = (self.S_RT == 0)
        self.S_RT[ w1 ] = self.S_bad
        w2 = (self.S_RB == 0)
        self.S_RB[ w2 ] = self.S_bad
        w3 = (self.S_LT == 0)
        self.S_LT[ w3 ] = self.S_bad
        w4 = (self.S_LB == 0)
        self.S_LB[ w4 ] = self.S_bad
                        
    #   get_qp_slopes()
    #---------------------------------------------------------------
    def get_qp_slopes_power(self):
    
        print('Computing qp slope powers...')
        self.S_RT_gam = self.S_RT ** (self.gamma - 1)
        self.S_RB_gam = self.S_RB ** (self.gamma - 1)
        self.S_LT_gam = self.S_LT ** (self.gamma - 1)
        self.S_LB_gam = self.S_LB ** (self.gamma - 1)
             
    #   get_qp_slopes_power()
    #---------------------------------------------------------------
    def get_qp_divergence_args(self):

        self.div_arg1_RT = self.S_RT_gam * self.dz_dx_R
        self.div_arg2_RT = self.S_RT_gam * self.dz_dy_T
        #-------------------------------------------------
        self.div_arg1_RB = self.S_RB_gam * self.dz_dx_R
        self.div_arg2_RB = self.S_RB_gam * self.dz_dy_B
        #-------------------------------------------------
        self.div_arg1_LT = self.S_LT_gam * self.dz_dx_L
        self.div_arg2_LT = self.S_LT_gam * self.dz_dy_T
        #-------------------------------------------------
        self.div_arg1_LB = self.S_LB_gam * self.dz_dx_L
        self.div_arg2_LB = self.S_LB_gam * self.dz_dy_B
                         
    #   get_qp_divergence_args()
    #---------------------------------------------------------------
    def get_qp_divergence(self):

        #----------------------------------------------------
        # Compute divergence (p-Laplacian) for Right-Top qp
        #----------------------------------------------------
        print('Computing qp divergence...')
        self.div_RT_x = self.get_d_dx_R( self.div_arg1_RT )
        self.div_RT_y = self.get_d_dy_T( self.div_arg2_RT )
        self.div_RT   = self.div_RT_x + self.div_RT_y
      
        #-------------------------------------------------------
        # Compute divergence (p-Laplacian) for Right-Bottom qp
        #-------------------------------------------------------
        self.div_RB_x = self.get_d_dx_R( self.div_arg1_RB )
        self.div_RB_y = self.get_d_dy_B( self.div_arg2_RB )
        self.div_RB   = self.div_RB_x + self.div_RB_y
          
        #----------------------------------------------------
        # Compute divergence (p-Laplacian) for Left-Top qp
        #----------------------------------------------------
        self.div_LT_x = self.get_d_dx_L( self.div_arg1_LT )
        self.div_LT_y = self.get_d_dy_T( self.div_arg2_LT )
        self.div_LT   = self.div_LT_x + self.div_LT_y
     
        #------------------------------------------------------
        # Compute divergence (p-Laplacian) for Left-Bottom qp
        #------------------------------------------------------
        self.div_LB_x = self.get_d_dx_L( self.div_arg1_LB )
        self.div_LB_y = self.get_d_dy_B( self.div_arg2_LB )
        self.div_LB   = self.div_LB_x + self.div_LB_y
                 
    #   get_qp_divergence()
    #---------------------------------------------------------------
    def get_qp_pde_lhs(self, READ_DEM=True):

        #------------------------------------------------------
        # For quarter pixels, compute:
        # d/dx[ S^(gam-1) * dz/dx] + d/dy[ S^(gam-1) * dz/dy]
        #
        # If ideal, steady-state landscape equation holds,
        # then this should equal -R*, around -0.0035.
        #------------------------------------------------------
        if (READ_DEM):
            self.read_dem()
        ## self.get_qp_grid()        #### better name?
        self.get_qp_gradients()
        self.get_qp_slopes()
        self.get_qp_slopes_power()
        self.get_qp_divergence_args()
        self.get_qp_divergence() 
        self.get_curvatures_RT()

        #-----------------------------------------------
        # Check the mins and maxes of the PDE operator
        #-----------------------------------------------
        div_RT_min = np.nanmin( self.div_RT )
        div_RT_max = np.nanmax( self.div_RT )
        div_RB_min = np.nanmin( self.div_RB )
        div_RB_max = np.nanmax( self.div_RB )
        #---------------------------------------
        div_LT_min = np.nanmin( self.div_LT )
        div_LT_max = np.nanmax( self.div_LT )
        div_LB_min = np.nanmin( self.div_LB )
        div_LB_max = np.nanmax( self.div_LB )
        
        if not(self.SILENT):
            print('min, max: div_RT =', div_RT_min, div_RT_max)
            print('min, max: div_RB =', div_RB_min, div_RB_max)
            print('min, max: div_LT =', div_LT_min, div_LT_max)
            print('min, max: div_LB =', div_LB_min, div_LB_max)
            print()
            print('gamma =', self.gamma)
            print('Finished.')
            print()

    #   get_qp_pde_lhs()
    #---------------------------------------------------------------
    def get_curvatures_RT( self ):
     
        z_x  = self.dz_dx_R
        z_y  = self.dz_dy_T
        S    = self.S_RT
        z_xx = self.get_d_dx_R( z_x )
        z_xy = self.get_d_dy_T( z_x )
        z_yy = self.get_d_dy_T( z_y )

        #----------------------------------------
        # Compute the plan or contour curvature
        #----------------------------------------
        kc = (z_y**2)*(z_xx) - (2*z_x*z_y*z_xy) + (z_x**2)*(z_yy) 
        kc *= -1 / S**3
        self.kappa_c_RT = kc
        
        #--------------------------------
        # Compute the profile curvature
        #--------------------------------
        kp = (z_x**2)*(z_xx) + (2*z_x*z_y*z_xy) + (z_y**2)*(z_yy)
        kp *= -1 / S**2
        self.kappa_p_RT = kp

    #   get_curvatures_RT()
    #---------------------------------------------------------------
    def fit_curvature_eqn(self):

        S  = self.S_RT.flatten()
        kc = self.kappa_c_RT.flatten()
        kp = self.kappa_p_RT.flatten()

        #----------------------        
        # Need to remove NaNs
        #----------------------
        w1 = np.logical_or( np.isnan(S), np.isnan(kc) )
        w2 = np.logical_or( w1, np.isnan(kp) )
        w3 = np.invert(w2)
        S  = S[w3]
        kc = kc[w3]
        kp = kp[w3]
        print('Original Mins and Maxes')
        print('npts =', S.size)
        print('min,max: S  =', S.min(),  ',', S.max()  )
        print('min,max: kc =', kc.min(), ',', kc.max() )
        print('min,max: kp =', kp.min(), ',', kp.max() )
        print()

        #-------------------------------------------------
        # Remove extreme values (main issue is kappa_c?)
        #-------------------------------------------------
#         wa = np.logical_or( kc < -0.3, kc > 0.3 )
#         wb = np.logical_or( kp < -0.3, kp > 0.3 )
#         wa = np.logical_or( kc < -0.1, kc > 0.1 )
#         wb = np.logical_or( kp < -0.1, kp > 0.1 )
#         wa = np.logical_or( kc < -0.05, kc > 0.05 )
#         wb = np.logical_or( kp < -0.05, kp > 0.05 )
#         wa = np.logical_or( kc < -0.3, kc > 0.05 )
#         wb = np.logical_or( kp < -0.3, kp > 0.05 )
#         wa = np.logical_or( kc < -0.3, kc > 0.1 )
#         wb = np.logical_or( kp < -0.3, kp > 0.0 )
#         wa = np.logical_or( kc < -0.3, kc > 0.3 )
#         wb = np.logical_or( kp < -0.3, kp > 0.0 )
        wa = np.logical_or( kc < -0.5, kc > 0.5 )
        wb = np.logical_or( kp < -0.5, kp > 0.0 )
        #---------------------------------------------
#         wa = np.logical_or( kc < -0.2, kc > 0.0 )  # helps
#         wb = np.logical_or( kp < -0.2, kp > 0.0 )
#         wa = np.logical_or( kc < -0.1, kc > 0.0 )  # helps
#         wb = np.logical_or( kp < -0.1, kp > 0.0 )
#         wa = np.logical_or( kc < -0.3, kc > 0.0 )  # helps
#         wb = np.logical_or( kp < -0.3, kp > 0.0 )
        wc = np.invert( np.logical_or( wa, wb) )
        S  = S[wc]
        kc = kc[wc]
        kp = kp[wc]
        print('New Mins and Maxes')
        print('npts =', S.size)
        print('min,max: S  =', S.min(),  ',', S.max()  )
        print('min,max: kc =', kc.min(), ',', kc.max() )
        print('min,max: kp =', kp.min(), ',', kp.max() )
        print()
               
        p0 = np.array([-0.7, 0.03], dtype='float64')  # initial guess
        method = 'lm'
        ## method = 'dogbox'
        # method can be 'lm', 'dogbox', or 'trf'
        
        popt, pcov = spo.curve_fit( curvature_eqn, (S,kc), kp,
                                    p0=p0, method=method )
        gamma = popt[0]
        Rs    = popt[1]
        print('best-fit gamma =', gamma)
        print('best_fit Rs =', Rs)
        print()
                    
    #   fit_curvature_eqn()
    #---------------------------------------------------------------
    def get_mask( self, mask_type='kc_kp_bad1' ): 

        npix = self.nx * self.ny
        ind  = np.arange(npix, dtype='int32')  # flattened array indices

        S    = self.S_RT.flatten()
        kc   = self.kappa_c_RT.flatten()
        kp   = self.kappa_p_RT.flatten()
        
        #-----------------------        
        # Remove NaNs in grids
        #-----------------------
        w1 = np.logical_or( np.isnan(S), np.isnan(kc) )
        w2 = np.logical_or( w1, np.isnan(kp) )
        w3 = np.invert(w2)  # values are not NaN
        S   = S[w3]
        kc  = kc[w3]
        kp  = kp[w3]
        ind = ind[w3]
 
        if (mask_type == 'kc_kp_bad1'):
            wa = np.logical_or( kc < -0.3, kc > 0.0 )
            wb = np.logical_or( kp < -0.3, kp > 0.0 )
            wc = np.logical_or( wa, wb)
        elif (mask_type == 'kc_big_neg'): 
            wc = ( kc < -0.3 )
        elif (mask_type == 'kp_big_neg'): 
            wc = ( kp < -0.3 )
        elif (mask_type == 'kc_big_pos'): 
            wc = ( kc > 0.3 )
        elif (mask_type == 'kp_big_pos'): 
            wc = ( kp > 0.3 )
        elif (mask_type == 'kp_pos'): 
            wc = ( kp > 0.0 )
        elif (mask_type == 'kc_neg'): 
            wc = ( kc < 0.0 )
        else:
            print('ERROR: mask_type not found:', mask_type)

        #------------------------------------------
        # Get the flattened array indices in mask
        #------------------------------------------
        ind = ind[wc]  # still dtype=int32
        n1  = ind.size
        print('Number of cells in grid =', npix)
        print('Number of cells in mask =', n1)
        print('n1/npix =', n1/npix )
        print()

        name = mask_type
        print('ind[0:50] =', ind[0:50])
        indices = np.zeros(n1 + 2, dtype='int32')
        indices[0]    = -1  # Needed for RTM files
        indices[1:-1] = ind
        indices[-1]   = -1
                        
        mask_file = self.dataset + '_' + name + '.rtm'
        mask_path = self.rtg_dir + mask_file
        mask_unit = open( mask_path, 'wb' )
        indices.tofile( mask_unit )
        mask_unit.close()
        
    #   get_mask()
    #---------------------------------------------------------------
    def get_value( self, var_name, replace_nan=None ):
    
        var = getattr(self, var_name)
        
        if (replace_nan is not None):
            w  = np.isnan(var)
            nw = w.sum()
            if (nw > 0):
                var[w] = replace_nan
        return var

    #   get_value()
    #---------------------------------------------------------------
    def save_var_to_rtg( self, var_name):

        grid = self.get_value( var_name )
        #-------------------------------
        # Write grid as binary to file
        #-------------------------------
        prefix   = self.dataset
        rtg_file = prefix + '_' + var_name + '.rtg'
        rtg_path = self.rtg_dir + rtg_file
        rtg_unit = open(rtg_path, 'wb')
        grid.tofile( rtg_unit )    
        rtg_unit.close()
 
        #-------------------------------------       
        # Make an RTI file for this RTG file
        #-------------------------------------
        rti = rti_files.make_info(grid_file=rtg_path,
              ncols=self.nx, nrows=self.ny,
              xres=self.dx, yres=self.dy,
              #---------------------------------
              data_source='hamiltonian.py',
              dtype=grid.dtype,  # will map to RT data type
              ## data_type='DOUBLE',
              ## byte_order='LSB',  # auto-provided
              pixel_geom=1,
              zres=0.001, z_units='METERS',
              box_units='METERS')
#               y_south_edge=0.0,   # auto-provided
#               y_north_edge=None,
#               x_west_edge=0.0,
#               x_east_edge=None,
#               gmin=-9999.0, gmax=-9999.0,
#               UTM_zone='UNKNOWN'
        rti_files.write_info(rtg_path, rti, SILENT=True)
            
    #   save_var_to_rtg()
    #---------------------------------------------------------------
    def get_div_RT(self):
    
        return self.div_RT
        
    #   get_div_RT()
    #---------------------------------------------------------------
    def get_div_RB(self):
    
        return self.div_RB
        
    #   get_div_RB()
    #---------------------------------------------------------------
    def get_div_LT(self):
    
        return self.div_LT
        
    #   get_div_LT()
    #---------------------------------------------------------------
    def get_div_LB(self):
    
        return self.div_LB
        
    #   get_div_LB()
    #---------------------------------------------------------------
    def get_histogram_mode(self, bins='doane', qp='RT',
                           PLOT=False, REMOVE_OUTLIERS=True,
                           vmin=None, vmax=None):

        #-----------------------------------------------
        # See: https://numpy.org/doc/stable/reference/
        #      generated/numpy.histogram.html
        # https://numpy.org/doc/stable/reference/generated/
        #   numpy.histogram_bin_edges.html#numpy.histogram_bin_edges
        #-----------------------------------------------
        # Note: Can set bins to integer, or strings
        #       like, 'auto', 'fd', and 'doane'
        #-----------------------------------------------        
        if (qp == 'RT'):
            a = self.get_div_RT()
        elif (qp == 'RB'):
            a = self.get_div_RB()
        elif (qp == 'LT'):
            a = self.get_div_LT()
        else:
            a = self.get_div_LB()

        w = np.isfinite(a)  # remove the nans
        b = a[w]
        self.npix = b.size   #######

        if (REMOVE_OUTLIERS):
            p = 0.1
            h, edges = np.histogram( b, bins=1000, range=(-20,20) )
            edges_L  = edges[:-1]
            hmax = h.max()
            w1   = (h > p * hmax)
            h    = h[ w1 ]
            e    = edges_L[ w1 ]
#             print('h =', h)
#             print('e =', e)
#             vmin = e[0]
#             vmax = e[-1]
            w2   = (h > 0.5 * hmax)
            e2   = e[ w2 ]
            self.width = (e2[-1] - e2[0])
            if not(self.SILENT):
                print('hmax initially =', hmax)
                print('Width at half maximum =', self.width)
                print()

        if ((vmin is not None) and (vmax is not None)):
            h, edges = np.histogram( b, bins=bins, range=(vmin,vmax) )
        else:
            h, edges = np.histogram( b, bins=bins )
        
        #-----------------------        
        # Find the modal value        
        #-----------------------
        nbins = h.size
        hmax  = h.max()
        w  = np.where( h == hmax)
        w0 = w[0][0]
        ## e0 = edges[ w0-1 ]
        e1 = edges[ w0 ]
        e2 = edges[ w0+1 ]
        
        if not(self.SILENT):
            print('nbins =', nbins)
            print('hmax =', hmax)
            print('w0 =', w0)
            print('left  edge of mode bin =', e1)
            print('right edge of mode bin =', e2)
            Rs = -0.0035
            if (e1 <= Rs) and (e2 > Rs):
                print('Rs =', Rs, 'is in mode bin.')
            else:
                print('Rs =', Rs, 'is not in mode bin.')        
            print()

        self.hmax = hmax
        if not(self.SILENT):
            print('gamma =', self.gamma)
            print('Max histogram height  =', self.hmax)

        #---------------------------------------        
        # Find width at half max (mode height)        
        #---------------------------------------
#         w1 = np.where( np.logical_and(h < hmax/2, edges[:-1] < e1) )
#         i1 = w1[0][-1] 
#         w2 = np.where( np.logical_and(h < hmax/2, edges[:-1] > e2) )
#         i2 = w2[0][0]
#         self.width = edges[i2] - edges[i1]
#         if not(self.SILENT):
#             print('Width at half maximum =', self.width)
#             print()


        if (PLOT):
            x_size = 8.0
            y_size = 4.0
            marker = ','  # pixel
            ## marker = '.'  # point
            figure = plt.figure(1, figsize=(x_size, y_size))
            ## plt.bar( edges[0:-1], h )     
            ## plt.plot( edges[0:-1], h, marker=marker)
            plt.plot( edges[0:-1], h, '.b-')
            plt.show()

    #   get_histogram_mode()
    #---------------------------------------------------------------
    def get_narrowest_histogram( self, bins=30,
                      vmin=-0.5, vmax=0.5,
                      gamma_min=-1.0, gamma_max=-0.5, delta=0.02 ):

        self.SILENT = True
        ### self.read_dem( REPORT=False )
        self.get_pde_solution()
        self.gamma = gamma_min
        while (self.gamma <= gamma_max):
            self.get_qp_p_laplacian( READ_DEM=False )
            self.get_histogram_mode( bins=bins, PLOT=False,
                                     vmin=vmin, vmax=vmax)
            print('gamma =', self.gamma, ', width =', self.width)
            print('hmax =', self.hmax, ', npix =', self.npix)
            print()
            self.gamma += delta

    #   get_narrowest_histogram()
    #---------------------------------------------------------------
    
#-------------------------------------------------------------------    
#-------------------------------------------------------------------
class profile():    

    #---------------------------------------------------------------
    def __init__( self):

        self.n_header = 6  # for RT4 profile files     
        self.dataset  = 'Beaver_Creek_Geo'  # default

    # __init__()
    #---------------------------------------------------------------
    def get_profile_file(self):    

        #----------------------------------------- 
        # This is original file used in papers,
        # extracted for Beaver_Creek_Geo DEM.
        #-----------------------------------------
        if (self.dataset == 'Beaver_Creek_Geo'):
            data_dir  = '/Users/peckhams/Dropbox/Meetings/'
            data_dir += 'Vienna_April_2024_EGU/Talk/'
            data_dir += 'Data_for_Beaver_Creek_profile/'
            filename = 'Beaver2_channelprof1.txt'
        elif (self.dataset == 'Beaver_Creek_1m'):
            data_dir  = '/Users/peckhams/DEMs/'
            data_dir += 'Beaver_Creek_1m/'
            filename  = 'Beaver_Creek_1m_mainprofile.txt'
        elif (self.dataset == 'Beaver_Creek_1m2'):
            data_dir  = '/Users/peckhams/DEMs/'
            data_dir += 'Beaver_Creek_1m/'
            filename  = 'Beaver_Creek_1m_rtsideprofile.txt'
        elif (self.dataset == 'BC_small_basin2'):
            data_dir  = '/Users/peckhams/DEMs/'
            data_dir += 'Beaver_Creek_1m/'
            filename  = 'BC_small_basin2_mainprofile.txt'
        elif (self.dataset == 'BC_small_basin2b'):
            data_dir  = '/Users/peckhams/DEMs/'
            data_dir += 'Beaver_Creek_1m/'
            filename  = 'BC_small_basin2b_mainprofile.txt'
        elif (self.dataset == 'Factory_Butte_Basin1'):
            data_dir  = '/Users/peckhams/DEMs/'
            data_dir += 'Factory_Butte_Utah_1m/'
            filename  = 'Factory_Butte_Basin1_channelprof1.txt'
        elif (self.dataset == 'Factory_Butte_Basin1_p2'):
            data_dir  = '/Users/peckhams/DEMs/'
            data_dir += 'Factory_Butte_Utah_1m/'
            filename  = 'Factory_Butte_Basin1_channelprof2.txt'
                        
        profile_file = data_dir + filename        
        self.profile_file = profile_file
        ### return profile_file

        #------------------------------
        # Count lines in profile_file
        #------------------------------
        n_lines  = 0
        file_unit = open( profile_file, 'r' )
        while (True):
            line = file_unit.readline()
            if (line == ''):
                break
            n_lines += 1     
        file_unit.close()
        
        # nx = 2426 for Beaver_Creek_Geo dataset
        self.nx = (n_lines - self.n_header)

    #   get_profile_file()
    #---------------------------------------------------------------
    def read_profile( self, head_clip=0, tail_clip=0,
                      clip_hilltop=False,
                      REPORT=True):

        #---------------------------------------------------
        # Note: This function reads a longitudinal profile
        #       extracted from a DEM with RiverTools 4.0.
        #---------------------------------------------------
        self.get_profile_file()

        #----------------------------------------    
        # Open file to read & skip header lines
        #----------------------------------------
        file_unit = open( self.profile_file, 'r' )
        for k in range( self.n_header ):
            line = file_unit.readline()
        
        #--------------------------------        
        # Read in the long profile data
        #--------------------------------
        # x units are km; z units are m
        #--------------------------------    
        x = np.zeros( self.nx, dtype='float64')
        z = np.zeros( self.nx, dtype='float64')  
    
        for k in range(self.nx):
            line   = file_unit.readline()
            vals   = line.split()     # split on whitespace
            x_val  = vals[0].strip()  # [km]
            z_val  = vals[1].strip()  # [m]
            x[ k ] = np.float64(x_val) * 1000.0   # [km to m]
            z[ k ] = np.float64(z_val)  # [m]
        file_unit.close()

        #-----------------------------
        # Option to clip the profile
        #-----------------------------
        if (head_clip > 0):
            x = x[head_clip:]
            z = z[head_clip:]
        if (tail_clip > 0):
            i2 = -1 * tail_clip   # i2 < 0
            x = x[:i2]
            z = z[:i2] 

        #-----------------------------
        # Compute the profile slopes
        #-----------------------------
        z_R = np.roll(z, -1)
        x_R = np.roll(x, -1)
        S   = (z - z_R) / (x_R - x)
        S[-1] = 0
        
        #-----------------------------------
        # Option to head clip hilltop;
        # i.e. values before the max slope
        #-----------------------------------
        if (clip_hilltop):
            Smax = S.max()
            w = np.where( S == Smax )
            w0 = w[0][0]
            if (w0 > 0):
                # Max slope is included
                x = x[w0:]
                z = z[w0:]
                S = S[w0:]
            print('## Clipped first', w0, 'profile points.')
                     
        #-------------------------
        # Save, x, z, x0, z0, S0
        #-------------------------
        self.nx = x.size
        self.x  = x
        self.z  = z
        self.S  = S
        self.x0 = x[0]
        self.z0 = z[0]
        self.S0 = (z[0] - z[1]) / (x[1] - x[0])
        
        if (REPORT):
            print('nx =', self.nx)
            print('x0 =', self.x0)
            print('z0 =', self.z0)
            print('S0 =', self.S0)
            print()

    #   read_profile()
    #---------------------------------------------------------------
    def smooth_profile( self, n_passes=1 ):

        #--------------------------------------    
        # Use a Gaussian kernel approximation
        #--------------------------------------
        kern_size = 10
        sigma = 1.0  # (standard deviation)
        kernel = sps.windows.get_window(('gaussian',sigma), kern_size)  
        z1 = self.z.copy()   # profile z-values
        for k in range(n_passes):
            #--------------------------------------------------------                   
            # Note: If mode='full', smoothed profile size will have
            # (n-1) extra elements, where n is the kernel size.
            #-------------------------------------------------------
            z2 = sps.convolve( z1, kernel, mode='same', method='auto')
            ## z2 = z2[2:-2, 2:-2]  # for mode='full'
            z1 = z2.copy()    
        
        self.z_smooth = z1

    #   smooth_profile() 
    #---------------------------------------------------------------   
    def get_peckham_profile( self, x, gamma, Rs):

        #-------------------------------------------------
        # For Beaver_Creek_Geo dataset
        # gamma=-0.701, Rs=0.0035, from paper via DAKOTA
        #-------------------------------------------------
        x0 = self.x0
        z0 = self.z0
        S0 = self.S0
        #----------------------------
        if (gamma != -1):
            p     = (gamma + 1) / gamma
            term1 = S0**(gamma+1)
            term2 = (S0**(gamma) + Rs * (x - x0))**p
            z = z0 + (term1 - term2) / (p * Rs)
        else:
            z = z0 - (1/Rs) * np.log(1 + S0*Rs*(x-x0))
        #-------------------------------------------
        # Also compute slope and profile curvature
        # Does kp = kappa_p = dS/dx?  Sign?
        #-------------------------------------------
        p2 = (1.0 / gamma)  # same sign as gamma
        S  = (S0**gamma + Rs*(x - x0))**p2
        kp = p2 * Rs * (S0**gamma + Rs*(x - x0))**(p2-1)
        self.z_peckham  = z
        self.S_peckham  = S
        self.kp_peckham = kp
        return z

    #   get_peckham_profile()
    #---------------------------------------------------------------
    def get_whipple_profile( self, x, c, p):
     
        x0 = self.x0
        z0 = self.z0
        z  = z0 - (c/p)*(x**p - x0**p)   # p != 0. 
        return z

    #   get_whipple_profile()
    #---------------------------------------------------------------
    def get_slope_profile( self ):

        #----------------------------------    
        # Must call read_profile() first.
        #----------------------------------
        x2 = self.x
        z2 = self.z
        x1 = np.roll(x2, -1)
        z1 = np.roll(z2, -1)
        x1[-1] = x1[-2]  ##### ???
        z1[-1] = z1[-2]  ##### ????
        self.S_profile = (z2 - z1) / (x1 - x2)

    #   get_slope_profile()
    #---------------------------------------------------------------
    def get_kappa_p_profile( self ):
    
        #-------------------------------    
        # Must call read_profile() and
        # get_slope_profile first.
        #-------------------------------
        x2 = self.x
        S2 = self.S_profile
        x1 = np.roll(x2, -1)
        S1 = np.roll(S2, -1)
        x1[-1] = x1[-2]  ##### ???
        S1[-1] = S1[-2]  ##### ????
        self.kp_profile = (S2 - S1) / (x2 - x1)

    #   get_kappa_p_profile()   
    #---------------------------------------------------------------  
    def plot_profile( self, curve='peckham', method='lm',
             var=None,
             head_clip=0, tail_clip=0, clip_hilltop=False, 
             ### x=None, y=None,       
             y2=None, y3=None, y4=None, y5=None,
             y6=None, y7=None, y8=None,
             xmin=None, xmax=None, ymin=None, ymax=None,
             x_name='Distance',  x_units='m', marker=',',
             title='Longitudinal Profile: ',
             y_name='Elevation', y_units='m',
             x_size=8, y_size=4, xticks=None, yticks=None):

        self.get_profile_file()
        title += self.dataset

        #-----------------------------------------        
        # Note: Can use i1 & i2 to clip profile.
        #-----------------------------------------
        self.read_profile(head_clip=head_clip, tail_clip=tail_clip,
                          clip_hilltop=clip_hilltop)

        x = self.x 
        z = self.z
        list1 = ['z_peckham', 'S_peckham', 'kp_peckham', 'kc_peckham_from_eqn']
        if (var in list1):
            curve = 'peckham'

        #-----------------------------------------
        # Get the best-fit curve data for y2
        # method can be 'lm', 'dogbox', or 'trf'
        #-----------------------------------------------------
        # Note that for dataset='Beaver_Creek_Geo', best-fit
        # parameters agree with those in Peckham_et_al_2016.
        #-----------------------------------------------------
        if (curve.lower() == 'peckham'):
            p0 = np.array([-0.701, 0.0035], dtype='float64')  # initial guess
            popt, pcov = spo.curve_fit(self.get_peckham_profile,
                                       x, z, p0=p0, method=method)
            gamma = popt[0]
            Rs    = popt[1]
            y2 = self.get_peckham_profile( x, gamma, Rs)
            print('best-fit gamma =', gamma)
            print('best_fit Rs =', Rs)
            print()
        elif (curve.lower() == 'whipple'):
            p0 = np.array([14.68, 0.133], dtype='float64')  # initial guess
            popt, pcov = spo.curve_fit(self.get_whipple_profile,
                                       x, z, p0=p0, method=method)
            c = popt[0]
            p = popt[1]
            y2 = self.get_whipple_profile( x, c, p)
            print('best-fit c =', c)
            print('best_fit p =', p)
            print()  
        else:
            # Don't do any curve fitting
            pass
             
        #------------------------------------
        # Which variable should be plotted?
        #------------------------------------
        if (var is not None):
            y2 = None
        if (var is None):
            y = z
        elif (var == 'S'):
            self.get_slope_profile()
            y = self.S_profile
            y_name = 'Slope (data)'; y_units = 'm/m'
        elif (var =='kp'):
            self.get_kappa_p_profile()
            y = self.kp_profile
            y_name = 'kappa_p'; y_units = '1/m'   #### check
        elif (var == 'smooth'):
            y = self.z_smooth
            y_name = 'z_smooth'; y_units = 'm'
        elif (var == 'z_peckham'):
            y = self.z_peckham
            y_name = var;  y_units = 'm'
        elif (var == 'S_peckham'):
            y = self.S_peckham
            y_name = var;  y_units = 'm/m'
        elif (var == 'kp_peckham'):
            y = self.kp_peckham
            y_name = var;  y_units = '1/m'
        elif (var == 'kc_peckham_from_eqn'):
            S_kp_tuple = (self.S_peckham, self.kp_peckham)
            y = kc_from_curvature_eqn( S_kp_tuple, gamma, Rs )
            y_name = var;  y_units = '??'
                                    
        #--------------------
        # Create the figure
        #--------------------
        figure = plt.figure(1, figsize=(x_size, y_size))
        # fig, ax = plt.subplots( figsize=(x_size, y_size))

        # Set the plot point marker
        # https://matplotlib.org/3.1.1/api/markers_api.html
        # marker = ','  # pixel
        # marker = '.'  # point (small circle)
        # marker = 'o'  # circle
        # marker = '+'
        # marker = 'x'
      
        plt.plot( x, y, marker=marker)
        if (y2 is not None):
            plt.plot(x, y2, marker=marker)
        if (y3 is not None):
            plt.plot(x, y3, marker=marker)
        if (y4 is not None):
            plt.plot(x, y4, marker=marker)
        if (y5 is not None):
            plt.plot(x, y5, marker=marker)
        if (y6 is not None):
            plt.plot(x, y6, marker=marker)
        if (y7 is not None):
            plt.plot(x, y7, marker=marker)
        if (y8 is not None):
            plt.plot(x, y8, marker=marker)
                        
        plt.xlabel( x_name + ' [' + x_units + ']' )
        plt.ylabel( y_name + ' [' + y_units + ']' )
        if (title is not None):
            plt.title( title )

        plt.ylim( ymin, ymax )
        plt.xlim( xmin, xmax )
        plt.xticks( xticks )
        plt.yticks( yticks )
        #-------------------------------------
        # This may be necessary depending on
        # the data type of ymin, ymax
        #-------------------------------------
        ## plt.ylim( np.array([ymin, ymax]) )
        ## plt.xlim( np.array([xmin, xmax]) )
        plt.show()

    #   plot_profile()
    #---------------------------------------------------------------
#-------------------------------------------------------------------
def kp_from_curvature_eqn( S_kc_tuple, gamma, Rs ):

    # Note: kp
    S, kc = S_kc_tuple
    kp = (S/gamma)*( (Rs/S**gamma) - kc)
    return kp

#   kp_from_curvature_eqn()
#-------------------------------------------------------------------
def kc_from_curvature_eqn( S_kp_tuple, gamma, Rs ):

    #--------------------------------------------------------
    # NOTE:  This is very close to zero when applied to
    #        the "Peckham profile", which is a 1D solution
    #        obtained for kappa_c = 0.
    #-------------------------------------------------------- 
    S, kp = S_kp_tuple
    kc = (Rs/S**gamma) - (kp * gamma / S)
    # kp = (S/gamma)*( (Rs/S**gamma) - kc)  
    return kc

#   kc_from_curvature_eqn()
#-------------------------------------------------------------------



    
    
    
    