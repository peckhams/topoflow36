
## Copyright (c) 2001-2009, Scott D. Peckham
## November 2009 (converted from IDL)

#-----------------------------------------------------------------------
#  Notes: Use the random midpoint displacement method to create
#         a fractal surface/landscape (due to Mandelbrot).
#         This can be used as initial surface for landscape
#         evolution models and is used by Alan Howard's
#         MARSSIM model.
#-----------------------------------------------------------------------
#
#  unit_test()
#  make_fractal_surface()
#
#-----------------------------------------------------------------------
from numpy import *
import numpy

from . import rtg_files
from . import rti_files

#-----------------------------------------------------------------------
def unit_test():

    z = make_fractal_surface(n_levels=8, H=1.5)

    print('min(z), max(z) =', z.min(), z.max())
    print('Finished with unit_test().')
    print(' ')
    
#   unit_test()
#-----------------------------------------------------------------------
def make_fractal_surface(n_levels, H=1.5, RTG_file=None,
                         sigma=float64(1),
                         scale=float64(1),
                         seed=168993,
                         X_WRAP=False, Y_WRAP=False,
                         SILENT=False):

    #---------------------------------------------------------
    # Notes: Can apply SCALE at very end.  A value of about
    #        0.01 should give results similar to Howard's
    #        MATRIX_2D with ERROR=0.02.

    #        H is a fractal exponent of some kind.
    
    #        Set the X_WRAP or Y_WRAP keywords in order to
    #        impose periodic boundaries on the left & right
    #        or top & bottom, respectively.

    #        If n_levels = 7,  nx = 129
    #        If n_levels = 8,  nx = 257
    #        If n_levels = 9,  nx = 513
    #        If n_levels = 10, nx = 1025
    #        If n_levels = 11, nx = 2049
    #----------------------------------------------------------

    if (n_levels > 11):    
        print('********************************************')
        print(' ERROR: Max number of levels is 11,')
        print(' which gives a grid size of 2049 x 2049.')
        print('********************************************')
        print(' ')
        return
    
    if not(SILENT):
        print('Creating a fractal surface...')
    
    #------------------
    # Initialize vars
    #------------------
    factor = float64(1) / sqrt(float64(2) ** H)    #############
    nx     = (int32(2) ** n_levels) + 1
    ny     = nx
    step   = nx - 1
    if not(SILENT):
        print('nx, ny =', nx, ',', ny)
    #----------------------------------------------
    x_vec      = numpy.arange(nx, dtype='Int16')
    y_vec      = numpy.arange(ny, dtype='Int16')
    cols, rows = numpy.meshgrid( x_vec, y_vec )
    ## rows = reshape(repeat(y_vec, nx), (ny, nx))
    ## cols = rot90(rows)  # (will work if nx=ny)
    
    sum_grid = (cols + rows)
    #----------------------------------------------
    DONE = zeros([ny, nx], dtype='UInt8')
    DONE[0,0]           = 1
    DONE[0,nx - 1]      = 1
    DONE[ny - 1,0]      = 1
    DONE[ny - 1,nx - 1] = 1
    #----------------------------------------------
    EDGE = zeros([ny, nx], dtype='UInt8')
    EDGE[:,0]      = 1
    EDGE[:,nx - 1] = 1
    EDGE[0,:]      = 1
    EDGE[ny - 1,:] = 1
    
    #------------------------------
    # Initialize grid of z-values
    #------------------------------
    numpy.random.seed(seed)
    v = random.normal(loc=0.0, scale=1.0, size=(2, 2))
    z = zeros([ny, nx], dtype='Float64')
    z[0,0]           = v[0,0]
    z[0,nx - 1]      = v[0,1]
    z[ny - 1,0]      = v[1,0]
    z[ny - 1,nx - 1] = v[1,1]
    #------------------------------------
    if (X_WRAP):    
        z[0,nx - 1]      = z[0,0]
        z[ny - 1,nx - 1] = z[ny - 1,0]
    if (Y_WRAP):    
        z[ny - 1,0]      = z[0,0]
        z[ny - 1,nx - 1] = z[0,nx - 1]
    #------------------------------------
    zF = z.flat    ## (an iterator to allow 1D indexing)  ##########
    
    for k in range( 1, (n_levels + 1) ):

        if not(SILENT):
            print('Working on level', k)
            
        step = (step / 2)
        
        #---------------------------------------
        # Get midpoint locations of this level
        #---------------------------------------
        w = where(logical_and(logical_and(logical_and(((cols.flat % step) == 0), \
                                                      ((rows.flat % step) == 0)),
                                          logical_not(DONE.flat)), logical_not(EDGE.flat)))
        n_mid = size(w[0])
        #########################
        #     Need this !!
        #########################
        w = w[0]
        
        #-----------------------------------------
        # Break these into two groups, w1 and w2
        #-----------------------------------------
        a1 = where((sum_grid.flat[w] % (2 * step)) == 0)  # (1D array)
        n1 = size(a1[0])
        a2 = where((sum_grid.flat[w] % (2 * step)) != 0)  # (1D array)
        n2 = size(a2[0])
    
        if (n1 != 0):
            w1 = w[a1[0]]
        if (n2 != 0):    
            w2 = w[a2[0]]

        #---------------------------------------------
        # Compute midpoint elevations as the average
        # of the diagonal neighbor elevations plus
        # a rescaled Gaussian random variable
        #---------------------------------------------
        UL = w1 - step * (nx + 1)
        UR = w1 - step * (nx - 1)
        LL = w1 + step * (nx - 1)
        LR = w1 + step * (nx + 1)
        #---------------------------
        ### numpy.random.seed(seed)
        ran = factor * sigma * random.normal(loc=0.0, scale=1.0, size=n1)
        zF[w1] = ((zF[UL] + zF[UR] + zF[LL] + zF[LR]) / float64(4)) + ran
        DONE.flat[w1] = 1
        
        #----------------------------------------------
        # Compute midpoint elevations of remaining
        # pixels at this scale as the average of the
        # nearest neighbor elevations plus a rescaled
        # Gaussian random variable.  n2=0 at start.
        #----------------------------------------------
        if (n2 != 0):    
            T = w2 - (step * nx)
            B = w2 + (step * nx)
            R = w2 + step
            L = w2 - step
            #----------------------------
            ### numpy.random.seed(seed)
            ran = factor * sigma * random.normal(loc=0.0, scale=1.0, size=n2)
            zF[w2] = ((zF[T] + zF[B] + zF[L] + zF[R]) / float64(4)) + ran
            DONE.flat[w2] = 1
        
        #--------------------------------------------
        # Compute elevations of edge pixels at this
        # scale as average of 3 nearest neighbors
        # plus a rescaled Gaussian random variable.
        #--------------------------------------------
        jump = (step * nx)
        #----------------------------
        L = where(logical_and(logical_and((cols.flat == 0), \
                                          ((rows.flat % step) == 0)), \
                              logical_not(DONE.flat)))
        nL = size(L[0])
        T  = L - jump
        B  = L + jump
        R  = L + step
        ### numpy.random.seed(seed)
        ran   = factor * sigma * random.normal(loc=0.0, scale=1.0, size=nL)
        zF[L] = ((zF[T] + zF[B] + zF[R]) / float64(3)) + ran
        DONE.flat[L] = 1
        #-----------------------------------------------------------------------------
        R = where(logical_and(logical_and((cols.flat == (nx - 1)), \
                                          ((rows.flat % step) == 0)), \
                              logical_not(DONE.flat)))
        nR = size(R[0])
        if not(X_WRAP):    
            L = R - step
            T = R - jump
            B = R + jump
            ### numpy.random.seed(seed)
            ran = factor * sigma * random.normal(loc=0.0, scale=1.0, size=nR)
            zF[R] = ((zF[L] + zF[T] + zF[B]) / float64(3)) + ran
        else:    
            zF[R] = zF[L]
        DONE.flat[R] = 1
        #-----------------------------------------------------------------------------
        T = where(logical_and(logical_and((rows.flat == 0), \
                                          ((cols.flat % step) == 0)), \
                              logical_not(DONE.flat)))
        nT = size(T[0])
        L  = T - step
        R  = T + step
        B  = T + jump
        ### numpy.random.seed(seed)
        ran = factor * sigma * random.normal(loc=0.0, scale=1.0, size=nT)
        zF[T] = ((zF[L] + zF[R] + zF[B]) / float64(3)) + ran
        DONE.flat[T] = 1
        #-----------------------------------------------------------------------------
        B  = where(logical_and(logical_and((rows.flat == (ny - 1)), \
                                           ((cols.flat % step) == 0)), \
                               logical_not(DONE.flat)))
        nB = size(B[0])
        if not(Y_WRAP):    
            L = B - step
            R = B + step
            T = B - jump
            ### numpy.random.seed(seed)
            ran = factor * sigma * random.normal(loc=0.0, scale=1.0, size=nB)
            zF[B] = ((zF[L] + zF[R] + zF[T]) / float64(3)) + ran
        else:    
            zF[B] = zF[T]
        DONE.flat[B] = 1
        #-----------------------------------------------------------------------------
    
    #-----------------------
    # Rescale the values ?
    #-----------------------
    if (scale != 1.0):   
        z = (z * scale)
    
    #-------------------------------------------
    # Option to save to RTG file with RTI file
    #-------------------------------------------
    if (RTG_file is not None):
        ###############################
        # CHECK FOR OVERWRITE HERE !
        ###############################
        ###############################
        info = rti_files.make_info( RTG_file, nx, ny,
                                    xres=100.0, yres=100.0,
                                    gmin=z.min(), gmax=z.max(),
                                    data_source="Midpoint displacement method" )
        rti_files.write_info(RTG_file, info)
        RTG_type = 'FLOAT'
        ## RTG_type = 'DOUBLE'
        rtg_files.write_grid( z, RTG_file, info, RTG_type=RTG_type)

    #----------------------
    # Print final message
    #----------------------
    if not(SILENT):
        print('Finished.')
        print(' ')

    return z

#   make_fractal_surface()
#-----------------------------------------------------------------------
