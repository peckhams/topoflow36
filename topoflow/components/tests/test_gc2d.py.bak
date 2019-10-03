## Copyright (c) 2001-2013, Scott D. Peckham

# This is not ready yet.  See gc2d.run_model().

import numpy as np
from topoflow.components import gc2d
## from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
def test_run( t_max=10.0, DEM_file='Animas_200.mat',
                      SILENT=False):

    Toggles.VARIABLE_DT_TOGGLE = 0  # (or change to 1)
    ###################################
    
    print 'Starting GC2D test run...'
    print 'Reading input file...'
    ( H, Zb, Zi, dx, dy ) = load_state(DEM_file=DEM_file,
                                       RESTART_TOGGLE = 0,
                                       INIT_COND_TOGGLE=1 )
    ny, nx = Zb.shape
    
    #------------------
    # Initialize vars
    #------------------
    t           = np.float64(0)
    conserveIce = np.float64(0)  # (total ice mass ??)
    meltrate    = np.zeros( (ny, nx), dtype='Float64' )
  
##    fd_watch = {}
##    fd_watch['thick']  = open( 'thickness_py.bin' , 'wb' )
##    counter = 0

    while (t < t_max):
            
        (dt, t, H, Zi, meltrate, conserveIce) = gc2d.update( t, H, Zb, dx, dy,
                                                     meltrate, conserveIce,
                                                     SILENT=SILENT)
##            COMPRESS_TOGGLE    = Toggles.COMPRESS_TOGGLE,
##            ICEFLOW_TOGGLE     = Toggles.ICEFLOW_TOGGLE,
##            ICESLIDE_TOGGLE    = Toggles.ICESLIDE_TOGGLE,
##            VARIABLE_DT_TOGGLE = Toggles.VARIABLE_DT_TOGGLE,
##            dtDefault=Parameters.dtDefault,
##            dtMax=Parameters.dtMax)       

    #-----------------------
    # Print a short report
    #-----------------------
    print ' '
    print '(nx, ny)       =', nx, ny
    print '(dx, dy)       =', dx, dy
    print '(Hmin, Hmax)   =', H.min(), H.max()
    print '(Zbmin, Zbmax) =', Zb.min(), Zb.max()
    print '(Zimin, Zimax) =', Zi.min(), Zi.max()
    print '(MRmin, MRmax) =', meltrate.min(), meltrate.max()
    print 'conserveIce    =', conserveIce
    print 'Finished.'
    print ' '
    
#   test_run()
#-----------------------------------------------------------------------

