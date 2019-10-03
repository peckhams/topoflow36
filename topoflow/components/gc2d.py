#!/usr/bin/env python

## Copyright (c) 2009, Scott D. Peckham
## Original Matlab version of GC2D, Mark Kessler
## GC2D first converted to Python/NumPy in 2009 by Eric Hutton

################################################################
# NOTE:  TopoFlow can provide "mass balance" for GC2D, but
#        the timescales are very different.  TopoFlow should
#        pass some kind of "net" or cumulative "mass balance"
#        to GC2D at its large timestep.
#
# NOTE:  There is no "load_mask()" function yet, but it is
#        called in a "try" block.
#
# NOTE:  THERMAL_TOGGLE option does not work yet.
#        See notes below regarding undefined vars.
#
# NOTE:  Should carefully test update_vars() due to
#        a bug fix and other changes to the code.
#        Compare to update_vars_OLD().
#
# NOTE:  Check that all "numpy" function calls include "numpy.".
#        Fixed calls to "mean()", "nonzero()", "ravel()",
#        abs() vs. absolute(), max(A,B) vs. maximum(A,B), etc.
#
################################################################

import numpy
import time
import sys
import logging
# import getopt

import scipy    # scipy.signal.convolve, scipy.io.loadmat
from scipy import interpolate
from scipy import signal

# SDP. 10/24/11.  No longer available.  Deprecated?
# from scipy.io.numpyio import fwrite  # used by print_watch_point()

#--------------------------------------------------------------------------------------------------

#   run_model()   # (for testing)

#   ------------------------------
#   Classes (used as structures)
#   ------------------------------
#   MassBalance
#   BoundaryCond
#   Parameters
#   InputParams
#   OutputParams
#   Toggles
#
#   -----------
#   Functions
#   -----------
#   compress_grid()
#   filter2d()
#   add_halo()
#   set_bc()
#   difference_grid()
#   basal_shear_stress()
#   iceflow()
#   ice_sliding()
#   sum_ice_motion()
#   avalanche()
#   calve()
#   mass_balance()
#   mass_conservation()
#   load_dem()
#   load_dem_var()
#   load_mask()             ###### Not written, but called.  #####
#   get_timestep()
#   update_vars()
#   print_watch_point()
#   update()

#   init_valley_glacier()
#   init_ice_sheet()
#   resample_dem()
#   init_ice_surface()
#   load_state()

#   #### load_state_old()
#   #### run_for()

#--------------------------------------------------------------------------------------------------                 
def run_model(t_max=10.0, DEM_file='Animas_200.mat', SILENT=False):

    Toggles.VARIABLE_DT_TOGGLE = 0  # (or change to 1)
    ###################################
    
    print('Starting GC2D test run...')
    print('Reading input file...')
    ( H, Zb, Zi, dx, dy ) = load_state(DEM_file=DEM_file,
                                       RESTART_TOGGLE = 0,
                                       INIT_COND_TOGGLE=1 )
    ny, nx = Zb.shape
    
    #------------------
    # Initialize vars
    #------------------
    t           = numpy.float64(0)
    conserveIce = numpy.float64(0)  # (total ice mass ??)
    meltrate    = numpy.zeros( (ny, nx), dtype='Float64' )
  
##    fd_watch = {}
##    fd_watch['thick']  = open( 'thickness_py.bin' , 'wb' )
##    counter = 0

    while (t < t_max):
            
        (dt, t, H, Zi, meltrate, conserveIce) = update( t, H, Zb, dx, dy,
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
    print(' ')
    print('(nx, ny)       =', nx, ny)
    print('(dx, dy)       =', dx, dy)
    print('(Hmin, Hmax)   =', H.min(), H.max())
    print('(Zbmin, Zbmax) =', Zb.min(), Zb.max())
    print('(Zimin, Zimax) =', Zi.min(), Zi.max())
    print('(MRmin, MRmax) =', meltrate.min(), meltrate.max())
    print('conserveIce    =', conserveIce)
    print('Finished.')
    print(' ')
    
#   run_model()
#-------------------------------------------------------------------------------------------------- 
class MassBalance:    # (enumeration)
   
    ( BAD_VAL , 
     ZERO_BALANCE ,
     CONSTANT_ELA ,
     ELA_LOWERING ,
     ELA_TIME_SERIES ,
     EXTERNAL_FUNC ,
     ELA_LOWERING2 ,
     BALANCE_FILE ,
     D180_TIME_SERIES ) = list(range( 9))

#    class MassBalance
#--------------------------------------------------------------------------------------------------        
class BoundaryCond:   # (enumeration)
   
    ( BAD_VAL ,
     ICE_FREE_BOUND ,
     ZERO_FLUX_BOUND ,
     CONST_FLUX_BOUND ,
     SURF_ELEV_BOUND ,
     SURF_SLOPE_BOUND ) = list(range( 6))

#    class BoundaryCond
#-------------------------------------------------------------------------------------------------- 
class Parameters:    # (structure)

    # Constants
    g        = numpy.float64(9.81)                 # gravitional acceleration  [m/s**2]
    rhoI     = numpy.float64(917)                  # density of ice    [kg/m**3]
    rhoW     = numpy.float64(1000)                 # density of water  [kg/m**3]
    day      = numpy.float64(0.00274)              # length of a day in years  [years]

    # Time
    t         = numpy.float64(0)                   # set time to zero
    tMax      = numpy.float64(100000)              # maximum simulation time in years
    dtMax     = numpy.float64(0.4 * 365*day)       # maximum timestep in years
    dtDefault = dtMax                              # timestep if VARIABLE_DT_TOGGLE==0
    
    sec_per_year = numpy.float64(3600) * 24 * 365  # (SDP, 9/30/09)
    
    # Glacier Properties
    MinGlacThick = numpy.float64(1)

    # Ice Deformation
    glensA = numpy.float64( (6.8e-15)*3.15e7/(1e9) )    # Patterson, 1994; MacGregor, 2000
    ## glensA = numpy.float64( 6.8 * 3.15 * 1e-17)
    
    # Attractor Sliding -- only if ICESLIDE_TOGGLE==1 (generally used)
    UsChar   = numpy.float64(10)
    taubChar = numpy.float64(100000)

    # Standard Sliding -- used if ICESLIDE_TOGGLE==2 (generally not used)
    B                 = numpy.float64(0.0012)     # m/(Pa*yr) -- MacGregor, 2000
    DepthToWaterTable = numpy.float64(20)         # distance from ice surface to water table
    MaxFloatFraction  = numpy.float64(80)         # limits water level in ice
    Hpeff             = numpy.float64(20)         # effective pressure (meters of water)
    
    # Mass Balance
    initELA         = numpy.float64(3350)         # (valley glaciers, try 3500 ice sheets)
    ELAStepSize     = numpy.float64(-50)
    ELAStepInterval = numpy.float64(500)
    gradBz          = numpy.float64(0.01)
    maxBz           = numpy.float64(2)
    tmin            = numpy.float64(200)          # Years, spin-up time
  
    # Avalanching
    angleOfRepose = numpy.float64(30)
    avalanchFreq  = numpy.float64(3)              # average number per year

    # Calving
    seaLevel    = numpy.float64(-100)             # meters
    calvingCoef = numpy.float64(2)                # year^-1

    # Thermal
    c      = numpy.float64(2060)                  # specific heat capacity (J/(kg*K))
    Qg     = numpy.float64(0.05 * 3.15e7)         # Geothermal heat flux (W/m^2)*seconds/year = (J/year)/(m^2)
    gradTz = numpy.float64(-0.0255)               # Geothermal Gradient

    # Only for Ice Sheets ???
    Hbound    = numpy.float64(2000)
    Elev0     = numpy.float64(0)           # reference elevation
    To        = numpy.float64(2.6)         # temperature at Elev0
    lapseRate = numpy.float64(-0.0065)     # degrees per meter
    
#   class Parameters
#-------------------------------------------------------------------------------------------------- 
class InputParams:    # (structure)
   
    CLEAR_FIGURE          = 1
    CONTOUR_INTERVAL      = 50.
    DEBUG_TOGGLE          = 0
    DT_LIMIT              = 0
    ELA_CONTOUR           = 1.
    ICE_CONTOUR           = 1.
    NEW_FIGURE            = 0
    QUIVER_VECS           = 0
    RECONSTRUCT           = 0
    SUBFIGURE             = 0
    THERMAL_CONTOUR       = 0

#   class InputParams
#-------------------------------------------------------------------------------------------------- 
class OutputParams:   # (structure)
    
    plotInterval = 60 * 120       # seconds
    saveInterval = 100            # whole years
    reportInterval = 30           # seconds

    nextPlot = 0                  # initialize to plot on first timestep
    nextSave = 0                  # initialize to save on first timestep
    nextReport = 0                # initialize to report on first timestep

    outputFile = 'savetmp'

#   class OutputParams
#-------------------------------------------------------------------------------------------------- 
class Toggles:   # (structure)

    #------------------------
    # Code behavior toggles
    #-----------------------------------------------------------
    # Toggle or turn on/off segments of the code or select 
    # between multiple possibilities for a given process.
    # Values can be reset in INIT_COND segment.
    # Note that many of these are unused in current version.
    #-----------------------------------------------------------
    GUISTART_TOGGLE     = 0   # started simulation with the gui   (off|on)
        
    SAVE_TOGGLE         = 1   # saving                            (off|on)
    PLOT_TOGGLE         = 1   # plotting                          (off|on)
    REPORT_TOGGLE       = 1   # reporting                         (off|on)
        
    COMPRESS_TOGGLE     = 0   # only simulate area with ice       (off|on)
    VARIABLE_DT_TOGGLE  = 0   # state dependent time step         (off|on)

    INIT_COND_TOGGLE    = 1   # load DEM and climate              (synth|valley|sheet)
    GENERIC_ICE_TOGGLE  = 0   # start with generic ice surface    (off|on)       
        
    ICEFLOW_TOGGLE      = 1   # ice motion by deformation         (off|on)
    ICESLIDE_TOGGLE     = 0   # ice motion by sliding             (off|on|select)        
        
    THERMAL_TOGGLE      = 0   # temp dependance of flow           (off|on)
    FREEZEON_TOGGLE     = 0   # basal ice freeze to bed           (off|on)

    AVALANCHE_TOGGLE    = 0   # avalanche off steep surfaces      (off|on)
    CALVING_TOGGLE      = 0   # calving front                     (off|on)
    ERODE_TOGGLE        = 0   # erode the bed                     (off|on|select)
    
##    CRN_TOGGLE          = 0   # CRN accumulation                  (off|on)

    # MASS_BALANCE_TOGGLE = MassBalance.ELA_LOWERING     # select climate scenerio   (off|on|select)
    MASS_BALANCE_TOGGLE = MassBalance.CONSTANT_ELA     # select climate scenerio   (off|on|select)

    WEST_BC_TOGGLE      = BoundaryCond.ICE_FREE_BOUND  # boundary condition    (no ice|reflect|no flow)
    EAST_BC_TOGGLE      = BoundaryCond.ICE_FREE_BOUND  # boundary condition    (no ice|reflect|no flow)
    SOUTH_BC_TOGGLE     = BoundaryCond.ICE_FREE_BOUND  # boundary condition    (no ice|reflect|no flow)
    NORTH_BC_TOGGLE     = BoundaryCond.ICE_FREE_BOUND  # boundary condition    (no ice|reflect|no flow)

#   class Toggles
#-------------------------------------------------------------------------------------------------- 
#-------------------------------------------------------------------------------------------------- 
def compress_grid( H , Zb , COMPRESS_TOGGLE=False , RESTART_TOGGLE=0,
                   THERMAL_TOGGLE=False ):

    # COMPRESS - ONLY SIMULATE SUB-RECTANGLE THAT CONTAINS ICE
    if (COMPRESS_TOGGLE) and (H.max() > 1) and (RESTART_TOGGLE != 2):
        H_FullSpace  = H.copy()
        Zb_FullSpace = Zb.copy()
        
        if (THERMAL_TOGGLE):
            Ts_FullSpace = Ts.copy()
            Tb_FullSpace = Tb.copy()
            Tm_FullSpace = Tm.copy()

        #[indrw,indcl] = find(H ~= 0);
        indrw, indcl = numpy.where( H != 0 )

        mxrw, mxcl = Zb.shape

        mnrw = max( 0    , min(indrw) - 2 )
        mxrw = min( mxrw , max(indrw) + 2 )
        mncl = max( 0    , min(indcl) - 2 )
        mxcl = min( mxcl , max(indcl) + 2 )

        H  = H [ mnrw:mxrw , mncl:mxcl ]
        Zb = Zb[ mnrw:mxrw , mncl:mxcl ]
        ## Zi = Zb + max( H, 0 )
        ## Zi = Zb + numpy.choose( H<0 , (H,0) )
        Zi = Zb + numpy.maximum(H, 0)

        if (THERMAL_TOGGLE):
            Ts = Ts[ mnrw:mxrw , mncl:mxcl ]
            Tb = Tb[ mnrw:mxrw , mncl:mxcl ]
            Tm = Tm[ mnrw:mxrw , mncl:mxcl ]

        ny, nx          = H.shape
        mx_ny, mx_nx    = Zb_FullSpace.shape
        ny, nx          = Zb.shape
        compression_ratio = (mx_nx * mx_ny) / (nx * ny)
        COMPRESSED_FLAG   = 1
    else:
        ## Zi = Zb + max( H, 0 ) # included for restarts
        ## Zi = Zb + numpy.choose( H<0 , (H,0) )
        Zi = Zb + numpy.maximum(H, 0)
        compression_ratio = 1.
        COMPRESSED_FLAG   = 0

    return ( Zi , compression_ratio , COMPRESSED_FLAG )

#   compress_grid()
#-------------------------------------------------------------------------------------------------- 
def filter2d( b , x , shape='same' ):
   
    return scipy.signal.convolve( b , x , mode=shape )

#   filter2d()
#-------------------------------------------------------------------------------------------------- 
def add_halo( x ):

    x_ext = numpy.concatenate( ( x[:,0,numpy.newaxis] , x     , x[:,-1,numpy.newaxis] ) , axis=1 )
    x_ext = numpy.concatenate( ( [x_ext[0,:]]         , x_ext , [x_ext[-1,:]]         ) )

    return x_ext

#  add_halo()
#-------------------------------------------------------------------------------------------------- 
def set_bc( H , Zb , Zi ,
            THERMAL_TOGGLE  = Toggles.THERMAL_TOGGLE,
            WEST_BC_TOGGLE  = Toggles.WEST_BC_TOGGLE,
            EAST_BC_TOGGLE  = Toggles.EAST_BC_TOGGLE,
            SOUTH_BC_TOGGLE = Toggles.SOUTH_BC_TOGGLE,
            NORTH_BC_TOGGLE = Toggles.NORTH_BC_TOGGLE ):
##            WEST_BC_TOGGLE  = BoundaryCond.ICE_FREE_BOUND ,
##            EAST_BC_TOGGLE  = BoundaryCond.ICE_FREE_BOUND ,
##            SOUTH_BC_TOGGLE = BoundaryCond.ICE_FREE_BOUND ,
##            NORTH_BC_TOGGLE = BoundaryCond.ICE_FREE_BOUND ):
      
    #-------------------------------------------------------
    # MODIFY BOUNDARY CELLS TO ENFORCE BOUNDARY CONDITIONS
    #-------------------------------------------------------
    # DEFAULT BOUNDARY CONDITION IS ZERO FLUX
    #-------------------------------------------------------
    H_ext  = add_halo( H )
    Zb_ext = add_halo( Zb )
    Zi_ext = add_halo( Zi )
        
    if (THERMAL_TOGGLE):
        Ts_ext = add_halo( Ts )
        Tb_ext = add_halo( Tb )
        Tm_ext = add_halo( Tm )

    # WESTERN BOUNDARY CONDITION
    if WEST_BC_TOGGLE == BoundaryCond.SURF_ELEV_BOUND:          # Constant Ice Surface Height
        ZiBound    = numpy.mean(Zb[:,0]) + Hbound
        H_ext[:,0] = ZiBound - Zb_ext[:,0]
    elif WEST_BC_TOGGLE == BoundaryCond.CONST_FLUX_BOUND:       # Constant Ice Flux B.C.
        pass
    elif WEST_BC_TOGGLE == BoundaryCond.SURF_SLOPE_BOUND:       # Constant Ice Surface Slope
        Zi_ext[:,0] = 2*Zi_ext[:,1] - Zi_ext[:,2]
        H_ext [:,0] = Zi_ext[:,0] - Zb_ext[:,0]
        H_ext [:,0] = numpy.maximum( H_ext[:,0], 0 )
    elif WEST_BC_TOGGLE == BoundaryCond.ICE_FREE_BOUND:         # Ice Free Boundary
        H_ext[:,0] = 0

    # EASTERN BOUNDARY CONDITION
    if EAST_BC_TOGGLE == BoundaryCond.SURF_ELEV_BOUND:          # Constant Ice Surface Height
        ZiBound     = numpy.mean(Zb[:,-1]) + Hbound
        H_ext[:,-1] = ZiBound - Zb_ext[:,-1]
    elif EAST_BC_TOGGLE == BoundaryCond.CONST_FLUX_BOUND:       # Constant Ice Flux B.C.
        pass
    elif EAST_BC_TOGGLE == BoundaryCond.SURF_SLOPE_BOUND:       # Constant Ice Surface Slope
        Zi_ext[:,-1] = 2*Zi_ext[:,-2] - Zi_ext[:,-3]
        H_ext [:,-1] = Zi_ext[:,-1] - Zb_ext[:,-1]
        H_ext [:,-1] = numpy.maximum( H_ext[:,-1], 0)
    elif EAST_BC_TOGGLE == BoundaryCond.ICE_FREE_BOUND:         # Ice Free Boundary
        H_ext[:,-1] = 0
            
    # SOUTHERN BOUNDARY CONDITION
    if SOUTH_BC_TOGGLE == BoundaryCond.SURF_ELEV_BOUND:         # Constant Ice Surface Height
        ZiBound    = numpy.mean(Zb[0,:]) + Hbound
        H_ext[0,:] = ZiBound - Zb_ext[0,:]
    elif SOUTH_BC_TOGGLE == BoundaryCond.CONST_FLUX_BOUND:      # Constant Ice Flux B.C.
        pass
    elif SOUTH_BC_TOGGLE == BoundaryCond.SURF_SLOPE_BOUND:      # Constant Ice Surface Slope
        Zi_ext[0,:] = 2*Zi_ext[1,:] - Zi_ext[2,:]
        H_ext [0,:] = Zi_ext[0,:] - Zb_ext[0,:]
        H_ext [0,:] = numpy.maximum( H_ext[0,:], 0 )
    elif SOUTH_BC_TOGGLE == BoundaryCond.ICE_FREE_BOUND:        # Ice Free Boundary
        H_ext[0,:] = 0
            
    # NORTHERN BOUNDARY CONDITION
    if NORTH_BC_TOGGLE == BoundaryCond.SURF_ELEV_BOUND:         # Constant Ice Surface Height
        ZiBound     = numpy.mean(Zb[-1,:]) + Hbound
        H_ext[-1,:] = ZiBound - Zb_ext[-1,:]
    elif NORTH_BC_TOGGLE == BoundaryCond.CONST_FLUX_BOUND:      # Constant Ice Flux B.C.
        pass
    elif NORTH_BC_TOGGLE == BoundaryCond.SURF_SLOPE_BOUND:      # Constant Ice Surface Slope
        Zi_ext[-1,:] = 2*Zi_ext[-2,:] - Zi_ext[-3,:]
        H_ext [-1,:] = Zi_ext[-1,:] - Zb_ext[-1,:]
        H_ext [-1,:] = numpy.maximum( H_ext[-1,:], 0 )
    elif NORTH_BC_TOGGLE == BoundaryCond.ICE_FREE_BOUND:        # Ice Free Boundary
        H_ext[-1,:] = 0
        
    Zi_ext = Zb_ext + H_ext

    return ( H_ext , Zb_ext , Zi_ext )

#   set_bc()
#-------------------------------------------------------------------------------------------------- 
def difference_grid( A , dx , dy ):

    dAdx_ext = ( A[:,1:] - A[:,:-1] ) / dx
    dAdy_ext = ( A[1:,:] - A[:-1,:] ) / dy
    dAdx     = dAdx_ext[1:-1,:]
    dAdy     = dAdy_ext[:,1:-1]

    return ( dAdx , dAdy )

#   difference_grid()
#-------------------------------------------------------------------------------------------------- 
#-------------------------------------------------------------------------------------------------- 
def basal_shear_stress( H_ext , Zi_ext , dx=1. , dy=1. ,
                        g=Parameters.g , rhoI=Parameters.rhoI ):

    #------------------------------------
    # CALCULATE THE BASAL SHEAR STRESS
    #------------------------------------            
    # forward differences  (could use difference_grid())
    dZidxX_ext = ( Zi_ext[:,1:] - Zi_ext[:,:-1] ) / dx
    dZidyY_ext = ( Zi_ext[1:,:] - Zi_ext[:-1,:] ) / dy

    dZidxX     = dZidxX_ext[1:-1,:]
    dZidyY     = dZidyY_ext[:,1:-1]

    HX_ext     = ( H_ext[:,1:] + H_ext[:,:-1] ) / 2.
    HY_ext     = ( H_ext[1:,:] + H_ext[:-1,:] ) / 2.
    HX         = HX_ext[1:-1,:]
    HY         = HY_ext[:,1:-1]
            
    taubxX_ext = -rhoI * g * HX_ext * dZidxX_ext
    taubyY_ext = -rhoI * g * HY_ext * dZidyY_ext

    taubxX     = taubxX_ext[1:-1,:]
    taubyY     = taubyY_ext[:,1:-1]

    taubxY = ( taubxX_ext[:-1,:-1] + taubxX_ext[:-1,1:] + 
               taubxX_ext[1: ,:-1] + taubxX_ext[1: ,1:] ) / 4.
            
    taubyX = ( taubyY_ext[:-1,:-1] + taubyY_ext[:-1,1:] +
               taubyY_ext[1: ,:-1] + taubyY_ext[1: ,1:] ) / 4.
            
    taubX  = numpy.sqrt( taubxX**2 + taubyX**2 )
    taubY  = numpy.sqrt( taubxY**2 + taubyY**2 )

    taubX  = numpy.choose( HX>0 , (0,taubX) )
    taubY  = numpy.choose( HY>0 , (0,taubY) )

    # Fill in zero values with 1 for use in division
    xcmpnt = numpy.choose( numpy.abs(taubX)<1e-5 , ( taubxX / taubX , 0. ) )
    ycmpnt = numpy.choose( numpy.abs(taubY)<1e-5 , ( taubyY / taubY , 0. ) )

    return ( ( xcmpnt , ycmpnt ) , ( taubX , taubY ) , ( HX , HY ) )

#   basal_shear_stress()
#-------------------------------------------------------------------------------------------------- 
def iceflow( taubX , taubY , HX , HY , xcmpnt , ycmpnt ,
             THERMAL_TOGGLE  = Toggles.THERMAL_TOGGLE,
##             THERMAL_TOGGLE=False,
             glensA = Parameters.glensA,
             #----------------------------------------------
             # Remaining values for THERMAL_TOGGLE = True
             #----------------------------------------------
             MinGlacThick = Parameters.MinGlacThick,
             lapseRate    = Parameters.lapseRate ):    # (value for ice sheets ???)
##             MinGlacThick = 1.0,
##             lapseRate = numpy.float64(-0.0065)):    # (value for ice sheets ???)

    #--------------------------------------------
    # CALCULATE ICE VELOCITY DUE TO DEFORMATION
    #--------------------------------------------        
    if (THERMAL_TOGGLE):

        ##################################################################
        # NOTE!  Many of the vars needed by this segment are undefined,
        #        such as: lapseRate (added above), eHs, eTs, eTm, To,
        #        H_ext, Ts_ext and Tm_ext. (SDP, 9/21/09)
        ##################################################################        
        A_ext = numpy.zeros(H_ext.shape , dtype='Float64' )
        ind   = numpy.nonzero( numpy.ravel(H_ext) >= MinGlacThick )

        Ts_ext = To + lapseRate*( Zi_ext - Elev0 )

        
        #A_ext(ind) = interp3( eHs, eTs, eTm, eA, H_ext(ind), Ts_ext(ind), Tm_ext(ind) ) ;
        try:
            numpy.put( A_ext , ind , interpolate.interp3d( eHs , eTs , eTm )( numpy.take(H_ext,ind) , numpy.take(Ts_ext,ind) , numpy.take(Tm_ext,ind) ) )
        except:
            logging.error( "NaN in A, likely H_node exceeds H_glens limits" )
            return -1
            
        AX = ( A_ext[1:-1, :-1] + A_ext[1:-1,1:  ] ) / 2.
        AY = ( A_ext[ :-1,1:-1] + A_ext[1:  ,1:-1] ) / 2.

    else:
        AX = glensA
        AY = glensA

    # Here's the guts of calculating the depth averaged velocity
    UdxX = numpy.abs( .4 * AX * taubX*taubX*taubX * HX ) * xcmpnt
    UdyY = numpy.abs( .4 * AY * taubY*taubY*taubY * HY ) * ycmpnt

    #UdxX = numpy.fix(UdxX*1e6)*1e-6
    #UdyY = numpy.fix(UdyY*1e6)*1e-6

    return ( UdxX , UdyY )

#   iceflow()
#-------------------------------------------------------------------------------------------------- 
def ice_sliding( taubX , taubY , xcmpnt , ycmpnt ,
                 THERMAL_TOGGLE=False,
                 FREEZEON_TOGGLE=False,
                 UsChar=Parameters.UsChar,
                 taubChar=Parameters.taubChar ):

    #------------------------------
    # CALCULATE SLIDING VELOCITY
    #------------------------------        
    # Here's the guts of calculating the sliding velocity 
    UsxX = numpy.choose( numpy.abs(taubX)<1e-5 , ( UsChar * numpy.exp(1 - taubChar / taubX) * xcmpnt ,
                                                   UsChar * numpy.exp(1 - taubChar        ) * xcmpnt ) )
    UsyY = numpy.choose( numpy.abs(taubY)<1e-5 , ( UsChar * numpy.exp(1 - taubChar / taubY) * ycmpnt , 
                                                   UsChar * numpy.exp(1 - taubChar        ) * ycmpnt ) )

    if (THERMAL_TOGGLE and FREEZEON_TOGGLE):
        
        ##################################################################
        # NOTE!  Many of the vars needed by this segment are undefined,
        #        such as:  Tb_ext, Zb_ext, seaLevel.  (SDP, 9/21/09)
        ################################################################## 
        ## notFrozen  = (Tb_ext > -.5) or (Zb_ext < seaLevel)
        notFrozen  = numpy.logical_or( Tb_ext > -0.5, Zb_ext < seaLevel )
        notFrozenX = ( notFrozen[1:-1, :-1] + notFrozen[1:-1,1:  ] ) / 2.
        notFrozenY = ( notFrozen[ :-1,1:-1] + notFrozen[1:  ,1:-1] ) / 2.

        UsxX *= notFrozenX
        UsyY *= notFrozenY

    return ( UsxX , UsyY )

#   ice_sliding()
#-------------------------------------------------------------------------------------------------- 
def sum_ice_motion( UdxX , UdyY , UsxX , UsyY ):

    UxX = (UdxX + UsxX)
    UyY = (UdyY + UsyY)

    return ( UxX , UyY )

#   sum_ice_motion()
#--------------------------------------------------------------------------------------------------             
def avalanche( H , angleOfRepose=Parameters.angleOfRepose ):

    #---------------------------------------
    # AVALANCHE SNOW OFF OF STEEP SURFACES
    #---------------------------------------------------------
    # move ice downslope until the ice surface is everywhere
    # less then or near the angle of repose
    #---------------------------------------------------------            
    ny, nx = Zb.shape
    dHRepose = dx * numpy.tan(angleOfRepose * numpy.pi / 180.)
    Ho       = numpy.maximum( H, 0 )
      
    while True:
        dZidx_down        = numpy.zeros( (ny,nx) , dtype='Float64' )
        dZidx_up          = numpy.zeros( (ny,nx) , dtype='Float64' )
        dZidx_down[:,1:]  = numpy.choose( Zi[:,1:]  < Zi[:,:-1] , ( Zi[:,1:]  - Zi[:,:-1] , 0 ) )
        dZidx_up  [:,:-1] = numpy.choose( Zi[:,:-1] < Zi[:,1:]  , ( Zi[:,:-1] - Zi[:,1:]  , 0 ) )
        dZidx             = numpy.choose( dZidx_up > dZidx_down , ( dZidx_down , dZidx_up ) )

        dZidy_left         = numpy.zeros( (ny,nx) , dtype='Float64' )
        dZidy_right        = numpy.zeros( (ny,nx) , dtype='Float64' )
        dZidy_left [1:,:]  = numpy.choose( Zi[1:,:] < Zi[:-1,:] , ( Zi[1:,:] - Zi[:-1,:] , 0 ) )
        dZidy_right[:-1,:] = numpy.choose( Zi[:-1,:] < Zi[1:,:] , ( Zi[:-1,:] - Zi[1:,:] , 0 ) )
        dZidy              = numpy.choose( dZidy_left > dZidy_right , ( dZidy_right , dZidy_left ) )

        grad  = numpy.sqrt( dZidx**2 + dZidy**2 )
        gradT = dZidy_left + dZidy_right + dZidx_down + dZidx_up
        gradT = numpy.choose( gradT == 0, (gradT,1) )
        grad  = numpy.choose( Ho < 0.1, (grad ,0) )

        mxGrad = grad.max()

        if (mxGrad <= 1.1*dHRepose):
            break

        delH = numpy.choose( grad < dHRepose , ( ( grad - dHRepose)/3. , 0 ) )

        Htmp = Ho.copy()
        Ho   = numpy.choose( Htmp<delH , ( Htmp-delH , 0 ) )
        delH = Htmp - Ho

        delHdn = numpy.zeros( (ny,nx) , dtype='Float64' )
        delHup = numpy.zeros( (ny,nx) , dtype='Float64' )
        delHlt = numpy.zeros( (ny,nx) , dtype='Float64' )
        delHrt = numpy.zeros( (ny,nx) , dtype='Float64' )

        delHup[:,1:  ] = delH[:, :-1] * dZidx_up  [:, :-1]  / gradT[:, :-1]
        delHdn[:, :-1] = delH[:,1:  ] * dZidx_down[:,1:  ]  / gradT[:,1:  ]
        delHrt[1:  ,:] = delH[ :-1,:] * dZidy_right[ :-1,:] / gradT[ :-1,:]
        delHlt[ :-1,:] = delH[1:  ,:] * dZidy_left [1:  ,:] / gradT[1:  ,:]

        Ho = Ho + delHdn + delHup + delHlt + delHrt
        Ho = numpy.maximum( Ho, 0 )

        Zi = Zb + Ho
            
    #H = Ho + (H<0).*H ;
    H = Ho + numpy.choose( H<0 , (0,H) )    ### DOUBLE-CHECK THIS

    return H
        
#   avalanche()
#-------------------------------------------------------------------------------------------------- 
def calve( H , dt , CALVING_TOGGLE=True ):

    if not(CALVING_TOGGLE):
        return
    
    #-------------------------
    # CALVING GLACIER FRONT
    #-----------------------------------------------------------------------
    # one reason this is difficult is that the height of ice in the cell
    # is really just recording the volume of ice, the position of the 
    # margin in the cell not the actual ice height.  Here floation
    # height is assumed (or higher if necessary to account for ice volume)
    #-----------------------------------------------------------------------        
    Hold      = H.copy()
    calvedIce = 0
        
    # Count time backwards with a sshorted timestep until the whole 
    # timestep used during this itteration has been simulated
        
    dtTot = dt
    while (dtTot > 0):
        # Find the calving front, aka the wet glacier margin
        G = H > 1
        W = numpy.logical_and( G==0 , Zb <= seaLevel )
        filt = numpy.array( [[0,1,0],[1,1,1],[0,1,0]] , dtype='Float64' )
        Wfilt = filter2d( filt , W )
        Wfilt[:,(0,-1)] = Wfilt[:,(2,-3)]
        Wfilt[(0,-1),:] = Wfilt[(2,-3),:]
        wetGmargin = Gi * Wfilt > 0
        indWGM = wetGmargin.ravel().nonzero()
 
        # If calving front exists, find water depth, ensure it's positive
        if (indWGM.size > 0):
            ## WDmarg = seaLevel - Zb.flatten()[indWGM]
            WDmarg = seaLevel - Zb.flat[indWGM]
            WDmarg = numpy.maximum( WDmarg, 0 )
            ind    = (WDmarg != 0).nonzero()
            indWGM = numpy.take( indWGM , ind )
            WDmarg = numpy.take( WDmarg , ind )

            #WDmarg = max( 0, seaLevel - Zb(indWGM) ) ;
            #ind = find( WDmarg == 0 ) ;
            #indWGM(ind) = [] ;
            #WDmarg(ind) = [] ;
            
        # If calving front exists, remove some ice
        if (indWGM.size > 0):
            # ice thickness in calving cells
            Hmarg = H.flatten()[indWGM]
            Hmarg = numpy.choose( Hmarg<WDmarg/0.917 , (Hmarg,WDmarg/0.917) )
                
            # A new timestep is calculated such that the calving rate times the 
            # timesstep does not exceed the total contents of any calving cell.
                
            dLinCalvdt   = calvingCoef * WDmarg                         # front migration rate
            dVolCalvdt   = dx * dLinCalvdt * Hmarg                      # rate of volume calved
            volAvailMarg = dx * dx * H.flatten()[indWGM]                # ice volume available
            calvDt = min( dtTot, ( volAvailMarg / dVolCalvdt ).min() )  # calving timestep

            # Remove this calving timestep from total time to calve
            dtTot = dtTot - calvDt
                
            # Convert the volume calved to ice thickness and remove
            calve = dVolCalvdt * calvDt / ( dx * dx )
            H[indWGM] = H[indWGM] - calve
            
            # Record total volume calved for posterity
            calvedIce = calvedIce + calve.sum(asis=0).sum() * dx * dx
                
        else:
            dtTot = 0
            
    # Record ice removal by calving for conservation test
    conserveIce = conserveIce + ( H - Hold ).sum(axis=0).sum()
     
#   calve()
#-------------------------------------------------------------------------------------------------- 
def mass_balance( Zi, t,
                  MASS_BALANCE_TOGGLE=None,
                  initELA=None, ELAStepSize=None, ELAStepInterval=None,
                  tmin=None, gradBz=None, maxBz=None ):                  
##                  MASS_BALANCE_TOGGLE=Toggles.MASS_BALANCE_TOGGLE,
##                  initELA=Parameters.initELA,
##                  tmin=Parameters.tmin ,
##                  ELAStepSize=Parameters.ELAStepSize ,
##                  ELAStepInterval=Parameters.ELAStepInterval ,
##                  gradBz=Parameters.gradBz,
##                  maxBz=Parameters.maxBz ):

    #------------------------------------------------------------
    # (12/4/09) Experiment that worked.  A function in another
    #  Python package can change variables stored in "classes"
    #  like Toggles and Parameters, but if given as defaults
    #  to the arguments of the update() function it will always
    #  use the original values in Toggles and Parameters.
    #------------------------------------------------------------
    if (MASS_BALANCE_TOGGLE == None):
        MASS_BALANCE_TOGGLE = Toggles.MASS_BALANCE_TOGGLE
    if (initELA == None):
        initELA = Parameters.initELA
    if (ELAStepSize == None):
        ELAStepSize = Parameters.ELAStepSize
    if (ELAStepInterval == None):
        ELAStepInterval = Parameters.ELAStepInterval
    if (tmin == None):
        tmin = Parameters.tmin
    if (gradBz == None):
        gradBz = Parameters.gradBz
    if (maxBz == None):
        maxBz = Parameters.maxBz

##    print 'MASS_BALANCE_TOGGLE =', MASS_BALANCE_TOGGLE
##    print 'initELA =', initELA
##    print 'ELAStepSize =', ELAStepSize
##    print 'ELAStepInterval =', ELAStepInterval
##    print 'tmin (spinup) =', tmin
    
    #--------------------------
    # CALCULATE MASS BALANCE
    #---------------------------------------------------------            
    # The imposed mass balance is the imposed climate.
    # There are many possibilities, here are only a few.
    # All must populate the 2D matrix Bxy of size = size(Zb)
    # with values of net precip/melt rate in m/yr.   ###################################
    # Define the scalar, ELA (m), for plotting.
    #---------------------------------------------------------
    if (MASS_BALANCE_TOGGLE == MassBalance.CONSTANT_ELA):
        # Simple ELA, maxBz, gradBz

        ELA = initELA
        #Bxy = min( maxBz , gradBz * ( Zi - ELA ) )
        Bxy = gradBz * ( Zi - ELA )
        Bxy = numpy.choose( Bxy > maxBz , (Bxy, maxBz) )
            
    elif (MASS_BALANCE_TOGGLE == MassBalance.ELA_LOWERING):
        # ELA changing with time experiment
            
        # ELAStepSize = -10 ;       # positive/negative values raise/lower ELA
        # ELAStepInterval = 500 ;
                
        ## ELA = initELA + ELAStepSize * max( 0 , numpy.floor( (t-tmin)/ELAStepInterval ) )
        ELA = initELA + ELAStepSize * numpy.maximum(0, (t-tmin)/ELAStepInterval )  # (SDP, 12/4/09)
        Bxy = gradBz * ( Zi - ELA )
        Bxy = numpy.choose( Bxy > maxBz , (Bxy, maxBz) )

        #----------------
        # For debugging
        #----------------
        DEBUG = False
        if (DEBUG):
            print('t, ELA =', t, ', ', ELA)
            print('min(Bxy), max(Bxy) =', Bxy.min(), ', ', Bxy.max())
        
    elif (MASS_BALANCE_TOGGLE == MassBalance.ELA_LOWERING2):
        # ELA changing with time experiment
            
        tau      = numpy.float64(25)          # intrinsic timescale of ice dynamics 
        tmin     = numpy.float64(0)           # time to begin ELA modification
        initELA  = numpy.float64(4200)        # initial ELA
        stepSize = numpy.float64(-10)         # positive/negative values raise/lower ELA
        dELAdt   = numpy.float64(-0.1) 
                
        ELA = initELA + stepSize * max( 0, numpy.floor( (t-tmin) / (8*tau) ) )
        Bxy = gradBz * ( Zi - ELA )
        Bxy = numpy.choose( Bxy > maxBz , (Bxy, maxBz) )
            
    elif (MASS_BALANCE_TOGGLE == MassBalance.EXTERNAL_FUNC):
        # external mass balance function
        try: Bxy
        except NameError:
            # Mass Balance 2D Must Return Bxy (2d Matrix)
            Bxy = mass_balance_gc2d( t , cellsize , Zi )
            nextGetBxy = t + getBxyInterval
        else:
            if (t >= nextGetBxy):
                Bxy = mass_balance_gc2d( t , cellsize , Zi )
                nextGetBxy = t + getBxyInterval

    elif (MASS_BALANCE_TOGGLE == MassBalance.ELA_TIME_SERIES) or \
         (MASS_BALANCE_TOGGLE == MassBalance.D18O_TIME_SERIES):
        # ELA time series
        ELA = interpolate.interp1d( trecord , ELArecord )( t )
        Bxy = gradBz * ( Zi - ELA )
        Bxy = numpy.choose( Bxy > maxBz , (Bxy, maxBz) )

    elif (MASS_BALANCE_TOGGLE == MassBalance.BALANCE_FILE):
        # external mass balance file
        Bxy = load_dem_var( DEM_file, 'Bxy' )
        ind = numpy.nonzero( numpy.ravel(numpy.abs(Bxy)==min(numpy.abs(Bxy))) )
        ELA = numpy.mean( numpy.take( numpy.ravel(Zi) , ind ) )

    elif (MASS_BALANCE_TOGGLE == MassBalance.ZERO_BALANCE):
        ELA = 0
        Bxy = numpy.zeros( Zb.shape , dtype='Float64' )
                
    else:
        logging.error( "Unrecognized Mass Balance" )
        return -1

    return ( Bxy , ELA )

#   mass_balance()
#--------------------------------------------------------------------------------------------------
def mass_conservation( H_ext , UxX , UyY , HX , HY , dZidxX , dZidyY ,
                       dx=1., dy=1.,
                       MinGlacThick=Parameters.MinGlacThick,
                       BoundaryFlux=0.,    ####  WAS UNDEFINED BEFORE 9/21/09 ####
                       WEST_BC_TOGGLE =Toggles.WEST_BC_TOGGLE,
                       EAST_BC_TOGGLE =Toggles.EAST_BC_TOGGLE,
                       SOUTH_BC_TOGGLE=Toggles.SOUTH_BC_TOGGLE,
                       NORTH_BC_TOGGLE=Toggles.NORTH_BC_TOGGLE ):  
##                       WEST_BC_TOGGLE =BoundaryCond.ICE_FREE_BOUND ,   # (Before 12/4/09)
##                       EAST_BC_TOGGLE =BoundaryCond.ICE_FREE_BOUND ,
##                       SOUTH_BC_TOGGLE=BoundaryCond.ICE_FREE_BOUND ,
##                       NORTH_BC_TOGGLE=BoundaryCond.ICE_FREE_BOUND ):

    #-----------------------------------
    # MASS CONSERVATION -- CONTINUITY
    #--------------------------------------------
    # Ensure that no ice is drawn from the rock
    # CLASS = H_ext >= MinGlacThick
    #--------------------------------------------
    CLASS = numpy.choose( H_ext >= MinGlacThick , (0.,1.) )
            
    DCLASSx = ( CLASS[1:-1,1:  ] - CLASS[1:-1, :-1] ) * numpy.sign( dZidxX )
    DCLASSy = ( CLASS[1:  ,1:-1] - CLASS[ :-1,1:-1] ) * numpy.sign( dZidyY )
            
    UxX = numpy.choose( numpy.abs(DCLASSx+1)<1e-5 , (UxX,0.) )
    UyY = numpy.choose( numpy.abs(DCLASSy+1)<1e-5 , (UyY,0.) )

    # Calculate both components of the ice flux
    qxX = UxX * HX
    qyY = UyY * HY

    #-----------------------------------------------------------
    # Note:  What is appropriate value for BoundaryFlux ??
    #        Was undefined in versions prior to 9/21/09 (SDP).
    #-----------------------------------------------------------
    if (WEST_BC_TOGGLE  == BoundaryCond.CONST_FLUX_BOUND): qxX[: , 0] = BoundaryFlux
    if (EAST_BC_TOGGLE  == BoundaryCond.CONST_FLUX_BOUND): qxX[: ,-1] = BoundaryFlux
    if (SOUTH_BC_TOGGLE == BoundaryCond.CONST_FLUX_BOUND): qyY[0 , :] = BoundaryFlux
    if (NORTH_BC_TOGGLE == BoundaryCond.CONST_FLUX_BOUND): qyY[-1, :] = BoundaryFlux
            
    # Here's the guts of the continuity equation
    dqdxX = ( qxX[ :,1:] - qxX[:  ,:-1] ) / dx
    dqdyY = ( qyY[1:, :] - qyY[:-1,:  ] ) / dy
    dHdt  = -dqdxX - dqdyY

    return ( dHdt , ( qxX , qyY ) )

#   mass_conservation()
#-------------------------------------------------------------------------------------------------- 
#-------------------------------------------------------------------------------------------------- 
def load_dem( DEM_file ):

    # Assume DEM_file is in MatLab format
    vars = scipy.io.loadmat( DEM_file )

    cellsize = numpy.float64(vars['cellsize'])
    easting  = numpy.float64(vars['easting'])
    northing = numpy.float64(vars['northing'])
    topo     = numpy.float64(vars['topo'])

    ny, nx = topo.shape

    logging.info( 'Shape of topo is %d by %d' , ny , nx )
    logging.info( 'Shape of easting is %d'    , easting.size )
    logging.info( 'Shape of northing is %d'   , northing.size )

    if (easting.size != nx):
        sys.exit( 'Easting does not match dimension of topo (%d != %d)' % (easting.size, nx) )
    if (northing.size != ny):
        sys.exit( 'Northing does not match dimension of topo (%d != %d)' % (northing.size, ny) )

    return ( topo , easting , northing , cellsize )

#   load_dem()
#-------------------------------------------------------------------------------------------------- 
def load_dem_var( var_file , val_s ):

    # Assume var_file is in MatLab format,
    # & maybe contains DEM as well.
    vars = scipy.io.loadmat( var_file )

    if (val_s in vars):
        var = vars[val_s]
    else:
        var = None

    return var

#   load_dem_var()
#--------------------------------------------------------------------------------------------------             
def get_timestep( H, Zi_ext, Zi , dHdt, Bxy,
                  dtMax = Parameters.dtMax,
                  dtDefault = Parameters.dtDefault ):

    #---------------------
    # CALCULATE TIMESTEP
    #-----------------------------------------------------------------------            
    # Now that we know the rate of change in ice surface heights due to  
    # ice motion and due to precipitation or melt we need to know over 
    # what period of time we can project forward with these rates and 
    # maintain stability of the ice surface.  The basic idea here is that
    # we don't want to take a timestep any longer then it would take to 
    # reverse the ice surface slope between two cells, such that ice 
    # should be flowing in the other direction.  In fact, let's make our 
    # timestep much less than that.
    #       
    # This calculation sets the timestep such that the change
    # in ice surface elevation nowhere exceeds a set fraction
    # of the local standard deviation in ice surface elevations
    #-----------------------------------------------------------------------
    
    # include ice changes by precip and melt
    dHdtTot = dHdt + Bxy
    adHdt   = numpy.abs(dHdtTot)
            
    # something like standard deviation of 3x3 cell areas around each cell
    filt   = numpy.ones( (3,3) , dtype='Float64' ) / 9.
    ZiMean = filter2d( filt , Zi_ext , 'valid' )
    dHmax  = numpy.sqrt( filter2d( filt, (ZiMean - Zi)**2 ) )
            
    # only consider cells with ice thickness > 10 m
    isGlac = (H > 10.)
            
    # find limiting timestep for each considered cell
    ind = ( numpy.logical_and( numpy.logical_and( adHdt!=0 , dHmax!=0 ) , isGlac!=0 ) ).flatten().nonzero()

    if (ind[0].size > 0):
        dtLimits = dHmax.flatten()[ ind ] / adHdt.flatten()[ ind ]
        dt       = dtLimits.min()
        idt      = ( dtLimits==dt ).nonzero()

        #ind = find( adHdt~=0 & dHmax~=0 & isGlac~=0 ) ;    
        #dtLimits = dHmax(ind)./adHdt(ind) ;
        #[dt, idt] = min( dtLimits ) ;
            
        # locate the x and y position of limiting cell for plotting
        #[rwDT,clDT] = ind2sub( size(adHdt), ind(idt) ) ; 
            
        # limit timestep to dtMax or some fraction of the calculated timestep
        if (dtMax is not None):
            dt = min( dtMax, dt/2. )
            
    else:
        # catch an error, (e.g. if H<10 in all cells )
        #if dt.size==0:
        dt = dtDefault

    #dt = numpy.fix(dt*1e6)*1e-6

    return dt

#   get_timestep()
#-------------------------------------------------------------------------------------------------- 
##def update_vars_OLD( H , Zb , Zi , Bxy , qxX , qyY , dHdt ,
##                     t , dt , conserveIce , dx=1. , dy=1. ):
##
##   t = t + dt
##
##   # numTimeSteps = numTimeSteps + 1 ;
##   # timeSteps(numTimeSteps) = dt ;
##            
##   # increase in ice thicknesses due to precip
##   Bxy_pos  = numpy.choose( Bxy>0 , (0,Bxy) )
##   H       += Bxy_pos * dt
##            
##   # change ice thicknesses due to ice motion
##   H       += dHdt * dt
##            
##   # decrease in ice thicknesses due to melt
##   Bxy_neg  =   numpy.choose( Bxy<0 , (0,Bxy))
##   Bxy_neg  = - numpy.choose( H<-Bxy_neg , (-Bxy_neg,H) )
##   H       += Bxy_neg * dt
##   
##   # record ice addition or removal by climate
##   snowFall    = ( Bxy_neg + Bxy_pos ) * dt
##   conserveIce = conserveIce + snowFall.sum(axis=0).sum()
##            
##   # record ice flux through boundaries
##   qbound = qyY[0,:].sum(axis=0).sum() - qyY[-1,:].sum(axis=0).sum() + \
##            qxX[:,0].sum(axis=0).sum() - qxX[:,-1].sum(axis=0).sum()
##   conserveIce = conserveIce + dt * qbound / dx
##            
##   Zi = Zb + numpy.choose( H<0 , (H,0) )
##        
##   if numpy.isnan(Zi).any():
##      #save workspacedump
##      logging.error( "NaN in ice thickness" )
##      return -1
##
##   return ( t , H , Zi , conserveIce )
##
###  update_vars_OLD()
#-------------------------------------------------------------------------------------------------- 
def update_vars( H , Zb , Zi , Bxy , qxX , qyY , dHdt ,
                 t , dt , conserveIce , dx=1. , dy=1. ):

    #-----------------------------------------------------
    # Note: This function should be checked carefully.
    #       Some parts below may not be correct.
    #       See notes below.  (SDP, 5/12/09)
    #-----------------------------------------------------
    t = t + dt

    # numTimeSteps = numTimeSteps + 1 ;
    # timeSteps(numTimeSteps) = dt ;

    #--------------------------------------------            
    # Increase in ice thicknesses due to precip
    #--------------------------------------------
    Bxy_pos  = numpy.maximum(0, Bxy)  # (SDP)
    H       += Bxy_pos * dt

    #--------------------------------------------           
    # Change ice thicknesses due to ice motion
    # (Should this be after next part ??)
    #--------------------------------------------
    H       += dHdt * dt

    #-------------------------------------------            
    # Decrease in ice thicknesses due to melt,
    # but dH in one timestep can't exceed H.
    # BUG FIX: Added dt in dH part.
    #---------------------------------------------------------
    # Bxy_neg  =   numpy.choose( Bxy<0 , (0,Bxy))
    # Bxy_neg  = - numpy.choose( H<-Bxy_neg , (-Bxy_neg,H) )
    # H       += Bxy_neg * dt
    #---------------------------------------------------------
    Bxy_neg  = numpy.minimum(0, Bxy)
    Bxy_neg  = numpy.maximum(-H/dt, Bxy_neg)
    dH_melt  = Bxy_neg * dt
    H       += dH_melt   # (add a negative quantity)

    #---------------------------------------------
    # Meltrate as water available to runoff,
    # defined here to be a positive contribution
    #---------------------------------------------
    # ice is less dense, so dH_water < dH_ice
    #-----------------------------------------------
    # rhoI = 917     # density of ice,   [kg/m**3]
    # rhoW = 1000    # density of water, [kg/m**3]
    #-----------------------------------------------
    density_ratio = (Parameters.rhoI / Parameters.rhoW)
    meltrate      = -(dH_melt / dt) * density_ratio
    meltrate      = (meltrate / Parameters.sec_per_year)  # [m/yr] -> [m/s]
    meltrate      = numpy.maximum(meltrate, 0, meltrate)
    #######################################################################

    #--------------------------------------------   
    # Record ice addition or removal by climate
    #---------------------------------------------------
    # NOTE: Minus sign added here for the mass-balance
    #       check instead of redefining Bxy_neg as was
    #       done in original version.
    #---------------------------------------------------
    snowFall    = ( -Bxy_neg + Bxy_pos ) * dt
    conserveIce += snowFall.sum(axis=0).sum()
    # conserveIce = conserveIce + snowFall.sum(axis=0).sum()

    #-------------------------------------            
    # Record ice flux through boundaries
    #-------------------------------------  
    qbound = qyY[0,:].sum(axis=0).sum() - qyY[-1,:].sum(axis=0).sum() + \
             qxX[:,0].sum(axis=0).sum() - qxX[:,-1].sum(axis=0).sum()
    conserveIce += (dt * qbound / dx)
    # conserveIce = conserveIce + dt * qbound / dx

    #--------------------------------------------------
    # QUESTION:  Can H be < 0 anywhere at this point?
    #--------------------------------------------------
    dz = numpy.maximum(H, 0)
    Zi = Zb + dz
    ## Zi = Zb + numpy.choose( H<0 , (H,0) )
       
    if numpy.isnan(Zi).any():
      # save workspace dump
      logging.error( "NaN in ice thickness" )
      return -1

    return ( t, H, Zi, meltrate, conserveIce )

#   update_vars()
#-------------------------------------------------------------------------------------------------- 
def print_watch_point( fd , x ):
   
    y = numpy.double( x )
    y.tofile( fd )    # (SDP. 10/24/11.  Untested, unused to avoid fwrite.)
    # fwrite( fd , y.size , y )
    # fd.flush()

#  print_watch_point()
#-------------------------------------------------------------------------------------------------- 
def update( t , H , Zb , dx , dy , meltrate, conserveIce,
            COMPRESS_TOGGLE=None, ICEFLOW_TOGGLE=None,
            ICESLIDE_TOGGLE=None, VARIABLE_DT_TOGGLE=None,
            dtDefault=None, dtMax=None, SILENT=False):
##            COMPRESS_TOGGLE    = Toggles.COMPRESS_TOGGLE,
##            ICEFLOW_TOGGLE     = Toggles.ICEFLOW_TOGGLE,
##            ICESLIDE_TOGGLE    = Toggles.ICESLIDE_TOGGLE,
##            VARIABLE_DT_TOGGLE = Toggles.VARIABLE_DT_TOGGLE,
##            dtDefault=Parameters.dtDefault,
##            dtMax=Parameters.dtMax,
##            SILENT=False):

    #--------------------------------------------------------------
    # (12/4/09) Experiment.  It seems that a function in another
    #  Python package can change variables stored in "classes"
    #  like Toggles and Parameters, but if given as defaults
    #  to the arguments of the update() function will always
    #  use the original values in Toggles and Parameters.
    #--------------------------------------------------------------
    if (COMPRESS_TOGGLE == None):
        COMPRESS_TOGGLE = Toggles.COMPRESS_TOGGLE
    if (ICEFLOW_TOGGLE == None):
        ICEFLOW_TOGGLE = Toggles.ICEFLOW_TOGGLE
    if (ICESLIDE_TOGGLE == None):
        ICESLIDE_TOGGLE = Toggles.ICESLIDE_TOGGLE
    if (VARIABLE_DT_TOGGLE == None):
        VARIABLE_DT_TOGGLE = Toggles.VARIABLE_DT_TOGGLE
    if (dtDefault == None):
        dtDefault = Parameters.dtDefault
    if (dtMax == None):
        dtMax = Parameters.dtMax
        
    #-----------------------------------------------------------
    # COMPRESS - ONLY SIMULATE SUB-RECTANGLE THAT CONTAINS ICE
    #-----------------------------------------------------------
    ( Zi , compression_ratio , COMPRESSED_FLAG ) = compress_grid( H, Zb, COMPRESS_TOGGLE=COMPRESS_TOGGLE )

    #-----------------------------------------------------------
    # MODIFY BOUNDARY CELLS TO ENFORCE BOUNDARY CONDITIONS
    #-----------------------------------------------------------
    ( H_ext , Zb_ext , Zi_ext ) = set_bc( H , Zb , Zi )
    ( dZidxX , dZidyY ) = difference_grid( Zi_ext , dx , dy )

    #-----------------------------------------------------------
    # CALCULATE THE BASAL SHEAR STRESS
    #-----------------------------------------------------------
    ( ( xcmpnt , ycmpnt ) , ( taubX , taubY ) , ( HX , HY ) ) = basal_shear_stress( H_ext , Zi_ext , dx=dx , dy=dy )

    #-----------------------------------------------------------
    # CALCULATE ICE VELOCITY DUE TO DEFORMATION
    #-----------------------------------------------------------
    if (ICEFLOW_TOGGLE):
        ( UdxX , UdyY ) = iceflow( taubX , taubY , HX , HY , xcmpnt , ycmpnt )
    else:
        UdxX = numpy.zeros( xcmpnt.shape , dtype='Float64' )    #### INEFFICIENT TO HAVE INSIDE LOOP
        UdyY = numpy.zeros( ycmpnt.shape , dtype='Float64' )

    #-----------------------------------------------------------
    # CALCULATE SLIDING VELOCITY
    #-----------------------------------------------------------
    if (ICESLIDE_TOGGLE):
        ( UsxX , UsyY ) = ice_sliding( taubX , taubY , xcmpnt , ycmpnt )
    else:
        UsxX = numpy.zeros( xcmpnt.shape , dtype='Float64' )     #### INEFFICIENT TO HAVE INSIDE LOOP
        UsyY = numpy.zeros( ycmpnt.shape , dtype='Float64' )

    #-----------------------------------------------------------
    # Sum all contributions to ice motion
    #-----------------------------------------------------------
    ( UxX , UyY ) = sum_ice_motion( UdxX , UdyY , UsxX , UsyY )

    #-----------------------------------------------------------
    # MASS CONSERVATION -- CONTINUITY
    #-----------------------------------------------------------
    ( dHdt , ( qxX , qyY ) ) = mass_conservation( H_ext , UxX , UyY , HX , HY ,
                                                  dZidxX , dZidyY , dx=dx , dy=dy );

    #-----------------------------------------------------------
    # CALCULATE MASS BALANCE
    #-----------------------------------------------------------
    ( Bxy , ELA ) = mass_balance( Zi , t )
        
    #-----------------------------------------------------------
    # CALCULATE TIMESTEP
    #-----------------------------------------------------------
    if (VARIABLE_DT_TOGGLE):
        dt = get_timestep( H , Zi_ext , Zi , dHdt , Bxy )
    else:
        dt = dtDefault   #### INEFFICIENT TO HAVE INSIDE LOOP

    #-----------------------------------------------------------
    # UPDATE the TIME and ICE THICKNESS
    #-----------------------------------------------------------
    ( t, H, Zi, meltrate, conserveIce ) = update_vars( H , Zb , Zi , Bxy ,
                                                       qxX , qyY , dHdt ,
                                                       t , dt , conserveIce,
                                                       dx=dx , dy=dy )
    #----------------
    # For debugging
    #----------------
    DEBUG = False
    if (DEBUG):
        print('(t, Hmax, MRmax) =', t, ", ", H.max(), ", ", meltrate.max())

    #----------------------------
    # Save H grids at intervals
    #----------------------------
    # need to pass counter
    #----------------------------
##    if ((counter % 1) == 0):
##        print_watch_point( fd_watch['thick'], H )
##    counter = counter + 1

    return (dt, t, H, Zi, meltrate, conserveIce)

#   update()   
#-------------------------------------------------------------------------------------------------- 
def init_valley_glacier( DEM_file='Animas_200.mat', dx_sim=200 ): #, mask_file=None ):

    ( topo , easting , northing , cellsize ) = load_dem( DEM_file )

    #---------------------------------------------------------
    # Resample DEM at new node spacing
    # Now done in "init_valley_glacier" and "init_ice_sheet"
    #---------------------------------------------------------
    dx = numpy.float64(dx_sim)  # set a new dx
    dy = numpy.float64(dx_sim)
    if (cellsize != dx):
        ( topo , easting , northing, cellsize) = resample_dem( topo, dx, dy, cellsize )
        
    #------------------------------------                
    # "watershed_mask" contains AAR and
    #  eroded volume watershed mask
    #------------------------------------
##    try:
##        # load( mask_file );
##        watershed_mask = load_mask( mask_file )     ########  There is no "load_mask()" function here.  #########
##    except:
##        # Use the whole grid if no watershed mask is available
##        watershed_mask = numpy.ones( topo.shape , dtype='Float64' )
##    
##       ## logging.warning( 'No watershed mask found; using the whole grid for AAR and eroded flux calculations.' )
##    
##       # Same message, but not colored red. 
##        print ' '
##        print 'No watershed mask found.  Using the whole grid'
##        print 'for AAR and eroded flux calculations.'
##        print ' '
        
      
##    Mass Balance
##    try:
##        initELA
##    except NameError:
    Parameters.initELA = numpy.float64(3350)
    Parameters.maxBz   = numpy.float64(2)
    Parameters.gradBz  = numpy.float64(1./100.)
        
##        initELA = numpy.float64(3350)   ############### Need to return or store these.
##        maxBz   = numpy.float64(2)
##        gradBz  = numpy.float64(1./100.)

    return ( topo, easting, northing, cellsize, dx, dy )

#   init_valley_glacier()
#-------------------------------------------------------------------------------------------------- 
def init_ice_sheet( DEM_file='Baffin200d', dx_sim=2000 ):

##    DEM_file = 'ValleyNonFjordTopo'

    ( topo , easting , northing, cellsize ) = load_dem( DEM_file )

    #---------------------------------------------------------
    # Resample DEM at new node spacing
    # Now done in "init_valley_glacier" and "init_ice_sheet"
    #---------------------------------------------------------
    dx  = numpy.float64(dx_sim)  # set a new dx
    dy  = numpy.float64(dx_sim)
    if (cellsize != dx):
        ( topo , easting , northing, cellsize) = resample_dem( topo, dx, dy, cellsize )

    ######################################################
    #  NOTE!!  None of these vars are saved or returned
    #          yet, so this is not finished.
    ######################################################
    
#     Parameters.UsChar   = numpy.float64(50)
#     Parameters.taubChar = numpy.float64(50000)
        
    # Calving
    Parameters.seaLevel = -100 ;    # meters
    Parameters.calvingCoef = 10 ;   # 1/year
                
    # Mass Balance
    Parameters.initELA   = numpy.float64(300)
    Parameters.maxBz     = numpy.float64(0.5)
    Parameters.gradBz    = numpy.float64(1./2000)
                
    Parameters.Hbound    = numpy.float64(2000)
 
    Parameters.Elev0     = numpy.float64(0)           # reference elevation
    Parameters.To        = numpy.float64(2.6)         # temperature at Elev0
    Parameters.lapseRate = numpy.float64(-0.0065)     # degrees per meter
                   
    Toggles.COMPRESS_TOGGLE         = 0
    Toggles.GENERIC_ICE_TOGGLE      = 0
    Toggles.MASS_BALANCE_TOGGLE     = MassBalance.CONSTANT_ELA
    Toggles.CALVING_TOGGLE          = 1
    Toggles.ERODE_TOGGLE            = 0
    Toggles.ICESLIDE_TOGGLE         = 0
        
##    Toggles.THERMAL_TOGGLE          = 0
##    Toggles.FREEZEON_TOGGLE         = 0
##    Toggles.HORZTL_ADVECT_TOGGLE    = 0
##    Toggles.GEOTHERMAL_HEAT_TOGGLE  = 0
##    Toggles.STRAIN_HEAT_TOGGLE      = 0
##    Toggles.SLIDING_HEAT_TOGGLE     = 0
##    Toggles.SURFACE_HEAT_FLUX_TOGGLE= 0
##    Toggles.THERMAL_3D_TOGGLE       = 0
        
    Toggles.WEST_BC_TOGGLE  = BoundaryCond.ZERO_FLUX_BOUND
    Toggles.EAST_BC_TOGGLE  = BoundaryCond.ZERO_FLUX_BOUND
    Toggles.SOUTH_BC_TOGGLE = BoundaryCond.ZERO_FLUX_BOUND
    Toggles.NORTH_BC_TOGGLE = BoundaryCond.ZERO_FLUX_BOUND

    return ( topo, easting, northing, cellsize, dx, dy )

#   init_ice_sheet()
#--------------------------------------------------------------------------------------------------                 
def resample_dem( topo, dx, dy, cellsize ):

    if (cellsize == dx):
        print('Error: DEM does not need to be resampled.')
        return -1

    print('Resampling DEM from: cellsize =', cellsize)
    print('    to cellsize =', dx)
    print(' ')

    ny, nx = topo.shape
    xOld = numpy.arange(nx-1) * cellsize
    yOld = numpy.arange(ny-1) * cellsize

    #xOld = (0:nx-1)*cellsize ;
    #yOld = (0:ny-1)*cellsize ;

    XOld,YOld = numpy.meshgrid( xOld , yOld )
        
    #if rem(max(xOld),dx) == 0 and rem(max(yOld),dy) == 0:
    if (max(xOld) % dx == 0) and (max(yOld) % dy == 0):
       nx_New = max(xOld)/dx + 1
       ny_New = max(yOld)/dy + 1
    else:
       nx_New = numpy.ceil( xOld[-1] / dx )
       ny_New = numpy.ceil( yOld[-1] / dy )
            
    x = numpy.arange(nx_New) * dx
    y = numpy.arange(ny_New) * dy
    X,Y = numpy.meshgrid( x , y )
       
    topo     = interpolate.interp2d( XOld , YOld , topo , kind='linear' )( X , Y )
    #topo     = interpolate.interp2d( XOld , YOld , topo, X, Y ) ;

    easting  = interpolate.interp1d( xOld , easting  , kind='linear' )( x )
    northing = interpolate.interp1d( yOld , northing , kind='linear' )( y )
    cellsize = dx

    return (topo, easting, northing, cellsize)

#   resample_dem()
#--------------------------------------------------------------------------------------------------                 
def init_ice_surface(Zb, ZiBound, Hbound, dx): 

    #---------------------------------------------------------             
    # First, rotate the topo such that the ice boundary
    # is on the left side of the simulation.
    # Need to check code; better to rotate DEM prior to use.
    #---------------------------------------------------------   
    ZiBound  = numpy.mean(Zb[:,0]) + Hbound
    H        = numpy.zeros(Zb.shape, dtype='Float64' )
    ny, nx = Zb.shape
    #----------------------------------
    taub      = 200000
    beta      = taub / (Parameters.rhoI * Parameters.g)
    jtermlast = (nx - 2)
    icefree   = False

    #---------------------------------------------------------           
    # For each row, find the cell for which the ice surface 
    # height at the left boundary would be ZiBound if the 
    # terminus that starts in that cell.
    #---------------------------------------------------------
    #for i =1:ny
    for i in range(ny):
    
        mZb   = Zb[i,:]
        slope = -numpy.diff(mZb) / dx

        #-------------------------------------------            
        # Search starts in front of the terminus
        # of the adjacent row that was just found.
        #-------------------------------------------
        jterm = min( jtermlast+1, nx-2 )
        while (jterm > 0):

            #------------------------        
            # Backwater calculation
            #------------------------
            mH = numpy.zeros(mZb.shape, dtype='Float64' )
            for j in range(jterm-1,-1,-1):
                term1  = ( -slope[j]/2. - (mH[j+1]/dx) )**2
                term2  = -(2./dx) * ( slope[j] * mH[j+1] - beta )
                deltaH = -slope[j]*dx/2. - mH[j+1] + dx * numpy.sqrt(term1+term2)
                mH[j]  = mH[j+1] + deltaH

            #------------------------------------------------                
            # The following ensures that the search for
            # the terminus was started beyond the terminus.
            #------------------------------------------------
            mZi = mZb + mH
            if (mZi[0] > ZiBound):
                icefree = True
            elif icefree and (mZi[0] < ZiBound):
                H[i,:] = mH
                jtermlast = jterm
                icefree = False
                break
            else:
                jterm = jterm + 2
                if (jterm >= nx-1):
                    logging.error( "Generic ice overruns boundary" )
                    return -1
                
            jterm = jterm - 1

    ######### CHECK INDENTATION/NESTING HERE #############
    Zi = Zb + H
                   
    ny, nx = Zb.shape
    filt    = numpy.ones( (3,3) , dtype='Float64' ) / 9
    ZiBig   = numpy.zeros( (ny+2,nx+2) , dtype='Float64' )
    ZiBig[1:-1,1:-1] = Zi
    
    for i in range(10):
       ZiBig[(0,-1),:]  = ZiBig[(1,-2),:]
       ZiBig[:,(0,-1)]  = ZiBig[:,(1,-2)]
       ZiBig            = filter2d( filt , ZiBig )
    
    Zi = ZiBig[1:-2,1:-2]
     
    indices = (H == 0)
    Zi[indices] = Zb[indices]

    return (H, Zi)

    # This is done in load_state().
##    conserveIce   = H.sum(axis=0).sum()
##    iceVolumeLast = conserveIce * dx * dy
    
#   init_ice_surface()    
#--------------------------------------------------------------------------------------------------                 
def load_state( DEM_file='Animas_200.mat',
                RESTART_TOGGLE = 0,
                INIT_COND_TOGGLE   = Toggles.INIT_COND_TOGGLE,
                GENERIC_ICE_TOGGLE = Toggles.GENERIC_ICE_TOGGLE ):

    if not((RESTART_TOGGLE == 0) or (RESTART_TOGGLE == 3)):
        # Check what H, Z, dx and dy should be.
        ###  return ( H , Zb , dx , dy )
        return -1

    #--------------------------------------------------------------
    # Note:  Toggles, Parameters, InputParams, OutputParams, etc.
    #        are now defined in "structures" at top of this file.
    #        Will add ability to read new values from a file.
    #--------------------------------------------------------------

    #------------------------------------
    # Reset the "code behavior toggles"
    #-----------------------------------------------------------
    # These toggles turn on/off segments of the code or select 
    # between multiple possibilities for a given process.
    # Values can be reset in INIT_COND segment.
    #-----------------------------------------------------------
    ## toggles = Toggles

    #---------------------------------------
    # Set numerical and physical constants
    #---------------------------------------
    ## params = Parameters

    #-------------------------------
    # Set various input parameters
    #-------------------------------
    ## inparams = InputParams
    
    #----------------------------------------------
    # Set parameters that control output behavior
    #----------------------------------------------
    ## outparams = OutputParams

    #--------------------------
    # Reset GUI parameters ??
    #--------------------------
    #if ( GUISTART_TOGGLE & exist('guiSimParams.mat','file') )
    #   load guiSimParams
    #   delete guiSimParams.mat
    #   clear newInitFile
    #elseif ( ~GUISTART_TOGGLE & exist( './guiPlotParams.mat', 'file' ) )
    #   delete guiPlotParams.mat
    #end

    #----------------------
    # Initialize counters
    #----------------------
    # numTimeSteps = 0
    # timeSteps = zeros(1000000, 1)

    #------------------------------------------------------
    # Initialize bed and ice topography, and climate vars
    #------------------------------------------------------
    # Must define topo, cellsize, dx, and dy
    #-----------------------------------------
    if (Toggles.INIT_COND_TOGGLE == 1):     # Valley glaciers
        ( topo, easting, northing, cellsize, dx, dy ) = init_valley_glacier( DEM_file )
    elif (Toggles.INIT_COND_TOGGLE == 2):   # Ice sheet
        ( topo, easting, northing, cellsize, dx, dy ) = init_ice_sheet( DEM_file )  ## , Bxy_file
##    elif (Toggles.INIT_COND_TOGGLE == 3):   # GUI start
##        ( topo, easting, northing, cellsize ) = load_dem( DEM_file )
##        dx = cellsize
##        dy = dx
##    else:  # Synthetic bedrock topography
##        # (INIT_COND_TOGGLE == 0) ??
##        logging.error( "Must code synthetic initial condition" )
##        return -1
    
    #-------------------------------------------------
    # Make sure "easting" and "northing" are defined
    # (Add this to "load_dem()" function instead.
    #-------------------------------------------------
    ny, nx = topo.shape
    #if !exist('easting') : easting  = numpy.arange( nx )
    #if !exist('northing'): northing = numpy.arange( ny )
    try:              easting
    except NameError: easting  = numpy.arange( nx )
    try:              northing
    except NameError: northing = numpy.arange( ny )

    #---------------------------------------------------------
    # Resample DEM at new node spacing
    # Now done in "init_valley_glacier" and "init_ice_sheet"
    #---------------------------------------------------------
##    if (cellsize != dx):
##        ( topo , easting , northing, cellsize) = resample_dem( topo, dx, dy, cellsize )

    #----------------------------------
    # Set the bed elevation to 'topo'
    #----------------------------------
    Zb     = topo.copy()
    initZb = Zb.copy()

    #--------------------------------------
    # Has ice depth, H, been defined yet?
    #--------------------------------------
    #if !exist('H'): H = numpy.zeros(Zb.shape)
    try: H
    except NameError: H = numpy.zeros( Zb.shape , dtype='Float64' )

    #--------------------------------------
    # Free surface = bed elev + ice depth
    #--------------------------------------
    Zi = (Zb + H)

    #---------------------------------
    # Create a generic ice surface ?
    #---------------------------------
    if (Toggles.GENERIC_ICE_TOGGLE):
        (H, Zi) = init_ice_surface(Zb, ZiBound, Hbound, dx)

    indices = (H == 0)
    Zi[indices] = Zb[indices]  # (won't matter unless we return Zi)
     
    #----------------------------------------------
    # Define X and Y grids (for plotting option?)
    #----------------------------------------------
    ny, nx = Zb.shape
    x   = numpy.arange( nx ) * dx
    y   = numpy.arange( ny ) * dy
    X,Y = numpy.meshgrid( x , y )
    
    #---------------------------------------------
    # Use these to track that ice is conserved ?
    #---------------------------------------------
    conserveIce   = H.sum(axis=0).sum()
    iceVolumeLast = conserveIce * dx * dy

    return ( H, Zb, Zi, dx, dy )

#   load_state()
#--------------------------------------------------------------------------------------------------                 
##def load_state_old( file, RESTART_TOGGLE=0, INIT_COND_TOGGLE=True,
##                    GENERIC_ICE_TOGGLE=False ):
##
##    if not((RESTART_TOGGLE == 0) or (RESTART_TOGGLE == 3)):
##        # Check what H, Z, dx and dy should be.
##        ###  return ( H , Zb , dx , dy )
##        return -1
##
####    if ((RESTART_TOGGLE == 0) or (RESTART_TOGGLE == 3)):
##        
##    #---------------------
##    # Load a saved state
##    #---------------------
##    # CODE BEHAVIOR TOGGLES
##    # toggles turn on/off segments of the code or select 
##    # between multiple possibilities for a given process
##    # values can be reset in INIT_COND segment
##
##    toggles = Toggles
##
##    # OUTPUT BEHAVIOR
##    plotInterval = 60 * 120       # seconds
##    saveInterval = 100            # whole years
##    reportInterval = 30           # seconds
##
##    nextPlot = 0                  # initialize to plot on first timestep
##    nextSave = 0                  # initialize to save on first timestep
##    nextReport = 0                # initialize to report on first timestep
##
##    outputFile = 'savetmp'
##
##    #---------------------------------------
##    # Set numerical and physical constants
##    #---------------------------------------
##    # Defined in Parameters class at top
##    #---------------------------------------
##    params = Parameters
##
##    #-------------------------
##    # RELOAD INPUT ARGUMENTS
##    #-------------------------
##    inputArgs = load_input_args
##
##    #if ( GUISTART_TOGGLE & exist('guiSimParams.mat','file') )
##    #   load guiSimParams
##    #   delete guiSimParams.mat
##    #   clear newInitFile
##    #elseif ( ~GUISTART_TOGGLE & exist( './guiPlotParams.mat', 'file' ) )
##    #   delete guiPlotParams.mat
##    #end
##
##    #----------------------
##    # INITIALIZE COUNTERS
##    #----------------------
##    # numTimeSteps = 0 ;
##    # timeSteps = zeros(1000000,1) ;
##
##    #-------------------------------------
##    # INITIALIZE BED and ICE TOPOGRAPHY,
##    # and CLIMATE VARIABLES
##    #-------------------------------------
##    # Must define topo, cellsize, dx, and dy
##
##    if INIT_COND_TOGGLE:
##
##     ### .mat file contains: 'topo' = matrix of bed elevations and 'cellsize', 
##     ### both in meters. 'easting' and 'northing' are included for plotting
##        
##     if INIT_COND_TOGGLE == 1:    # Valley glaciers
##            
##    #                 DEM_file = 'Yosemite200_rot35_400x650' ;
##    #                 DEM_file = 'Nederland100' ;
##    #                 DEM_file = 'KingsCanyon200Rot256x256shift' ;
##    #                 DEM_file = 'sample200' ;
##    #                 DEM_file = 'animas_200' ;
##    #                 DEM_file = '4J_newDEM_200' ;
##    #                 DEM_file = 'reproj4j_200' ;
##        DEM_file = file
##        DEM_file = 'Animas_200.mat'
##
##        #load( DEM_file ) ;
##        ( topo , easting , northing , cellsize ) = load_dem( DEM_file )
##
##        dx = numpy.float64(200)  # set a new dx
##        dy = numpy.float64(dx)
##            
##        # AAR and eroded volume watershed mask
##        mask_file = 'watershed_mask'
##            
##        try:
##           #load( mask_file );
##           watershed_mask = load_mask( mask_file )
##        except:
##           watershed_mask = numpy.ones( topo.shape , dtype='Float64' )
##           # Use the whole grid if no watershed mask is available
##           msg  = 'No watershed mask found; using the whole grid '
##           msg += 'for AAR and eroded flux calculations.'
##           logging.warning( msg )
##
##           # Same message, but not colored red.               
##           print ' '
##           print 'No watershed mask found; using the whole grid'
##           print 'for AAR and eroded flux calculations.'
##           print ' '
##
##        # Mass Balance
##
##        try: initELA
##        except NameError:
##           initELA = numpy.float64(3350)
##           maxBz   = numpy.float64(2)
##           gradBz  = numpy.float64(1./100.)
##
##     elif INIT_COND_TOGGLE==2:    # Ice sheets
##
##        DEM_file = 'Baffin200d'
##        DEM_file = 'ValleyNonFjordTopo'
##            
##        #load( DEM_file ) ;
##        ( topo , easting , northing ) = load_dem( DEM_file )
##                        
##        dx       = numpy.float64(2000)  # set a new dx
##        dy       = dx 
##            
##        UsChar   = numpy.float64(100)
##        taubChar = numpy.float64(50000)
##            
##        #load( DEM_file, 'Bxy' ) ;
##        Bxy = load_dem_var( DEM_file , 'Bxy' )
##            
##        # Mass Balance
##        initELA   = numpy.float64(3500)
##        maxBz     = numpy.float64(0)
##        gradBz    = numpy.float64(1./100)
##            
##        Hbound    = numpy.float64(2000)
##            
##        Elev0     = numpy.float64(0)           # reference elevation
##        To        = numpy.float64(-30)         # temperature at Elev0
##        lapseRate = numpy.float64(-0.0065)     # degrees per meter
##               
##        COMPRESS_TOGGLE         = 0
##        GENERIC_ICE_TOGGLE      = 0
##        MASS_BALANCE_TOGGLE     = ELA_TIME_SERIES
##        CALVING_TOGGLE          = 1
##        ERODE_TOGGLE            = 0
##
##        THERMAL_TOGGLE          = 0
##        FREEZEON_TOGGLE         = 0
##        HORZTL_ADVECT_TOGGLE    = 0
##        GEOTHERMAL_HEAT_TOGGLE  = 0
##        STRAIN_HEAT_TOGGLE      = 0
##        SLIDING_HEAT_TOGGLE     = 0
##        SURFACE_HEAT_FLUX_TOGGLE= 0
##        THERMAL_3D_TOGGLE       = 0
##
##        WEST_BC_TOGGLE      = ZERO_FLUX_BOUND
##        EAST_BC_TOGGLE      = ZERO_FLUX_BOUND
##        SOUTH_BC_TOGGLE     = ZERO_FLUX_BOUND
##        NORTH_BC_TOGGLE     = ZERO_FLUX_BOUND
##
##     elif INIT_COND_TOGGLE == 3:    # gui_start
##        #load( DEM_file ) ;
##        ( topo , easting , northing ) = load_dem( DEM_file )
##        dy = dx
##         
##     ny, nx = topo.shape
##     #if !exist('easting') : easting  = numpy.arange( nx )
##     #if !exist('northing'): northing = numpy.arange( ny )
##
##     try:              easting
##     except NameError: easting  = numpy.arange( nx )
##     try:              northing
##     except NameError: northing = numpy.arange( ny )
##     
##                
##     # resample DEM at new node spacing
##     if cellsize != dx:
##            
##        ny, nx = topo.shape
##        xOld = numpy.arange(nx-1)*cellsize
##        yOld = numpy.arange(ny-1)*cellsize
##
##        #xOld = (0:nx-1)*cellsize ;
##        #yOld = (0:ny-1)*cellsize ;
##
##        XOld,YOld = numpy.meshgrid( xOld , yOld )
##            
##        #if rem(max(xOld),dx) == 0 and rem(max(yOld),dy) == 0:
##        if max(xOld) % dx == 0 and max(yOld) % dy == 0:
##           nx_New = max(xOld)/dx + 1
##           ny_New = max(yOld)/dy + 1
##        else:
##           nx_New = numpy.ceil( xOld[-1] / dx )
##           ny_New = numpy.ceil( yOld[-1] / dy )
##                
##        x = numpy.arange(nx_New) * dx
##        y = numpy.arange(ny_New) * dy
##
##        X,Y = numpy.meshgrid( x , y )
##           
##        topo     = interpolate.interp2d( XOld , YOld , topo , kind='linear' )( X , Y )
##        #topo     = interpolate.interp2d( XOld , YOld , topo, X, Y ) ;
##
##        easting  = interpolate.interp1d( xOld , easting  , kind='linear' )( x )
##        northing = interpolate.interp1d( yOld , northing , kind='linear' )( y )
##        cellsize = dx
##
##     # Set the bed elevation to 'topo'
##     Zb     = topo.copy()
##     initZb = Zb.copy()
##
##     #if !exist('H'): H = numpy.zeros(Zb.shape)
##     try: H
##     except NameError: H = numpy.zeros( Zb.shape , dtype='Float64' )
##
##     Zi     = H + Zb
##     #clear topo
##        
##     ny, nx = Zb.shape
##     x   = numpy.arange( nx ) * dx
##     y   = numpy.arange( ny ) * dy
##     X,Y = numpy.meshgrid( x , y )
##
##     # Create a generic ice surface
##     if GENERIC_ICE_TOGGLE:
##            
##        # This code segment rotates the topo such that the 
##        # ice boundary is on the left side of the simulation
##        # need to check code; better to rotate DEM prior to use
##            
##        ZiBound = numpy.mean(Zb[:,0]) + Hbound
##        taub    = 200000
##        H       = numpy.zeros(Zb.shape, dtype='Float64' )
##        ny, nx = Zb.shape
##        beta    = taub/(rhoI*g)
##        jtermlast = nx-2
##        icefree = 0
##            
##        # for each row, find the cell for which the ice surface 
##        # height at the left boundary would be ZiBound if the 
##        # terminus that starts in that cell
##        #for i =1:ny
##        for i in range(ny):
##        
##           mZb   = Zb[i,:]
##           slope = -numpy.diff(mZb)/dx
##                
##           # search starts in front of the terminus
##           # of the adjacent row that was just found  
##           jterm = min( jtermlast+1, nx-2 )
##           while jterm > 0:
##            
##              # backwater calculation
##              mH = numpy.zeros(mZb.shape, dtype='Float64' )
##              for j in range(jterm-1,-1,-1):
##        
##                 term1  = ( -slope[j]/2. - (mH[j+1]/dx) )**2
##                 term2  = -(2./dx) * ( slope[j] * mH[j+1] - beta )
##                 deltaH = -slope[j]*dx/2. - mH[j+1] + dx * numpy.sqrt(term1+term2)
##                 mH[j]  = mH[j+1] + deltaH
##                    
##              # the following ensures that the search for
##              # the terminus was started beyond the terminus
##              mZi = mZb + mH
##              if mZi[0] > ZiBound:
##                 icefree = 1
##              elif icefree and mZi[0] < ZiBound:
##                 H[i,:] = mH
##                 jtermlast = jterm
##                 icefree = 0
##                 break
##              else:
##                 jterm = jterm + 2
##                 if jterm >= nx-1:
##                    logging.error( "Generic ice overruns boundary" )
##                    return -1
##                    
##              jterm = jterm - 1
##        
##        Zi = Zb + H
##                       
##        ny,nx = Zb.shape
##        filt    = numpy.ones( (3,3) , dtype='Float64' ) / 9
##        ZiBig   = numpy.zeros( (ny+2,nx+2) , dtype='Float64' )
##        ZiBig[1:-1,1:-1] = Zi
##        
##        for i in range(10):
##           ZiBig[(0,-1),:]  = ZiBig[(1,-2),:]
##           ZiBig[:,(0,-1)]  = ZiBig[:,(1,-2)]
##           ZiBig            = filter2d( filt , ZiBig )
##        
##        Zi = ZiBig[1:-2,1:-2]
##        
##     ind = (H == 0)
##     Zi[ind] = Zb[ind]
##        
##     conserveIce   = H.sum(axis=0).sum()
##     iceVolumeLast = conserveIce * dx * dy
##        
##    else:  # SYNTHETIC BEDROCK TOPOGRAPHY
##        logging.error( "Must code synthetic initial condition" )
##        return -1
##
##    return ( H , Zb , dx , dy )
##
###  load_state_old()
#-------------------------------------------------------------------------------------------------- 
##def run_for( t , t_max , H , Zb , dx , dy , meltrate,
##             ICEFLOW_TOGGLE=True , ICESLIDE_TOGGLE=False ,
##             VARIABLE_DT_TOGGLE=True ,
##             dtDefault=Parameters.dtDefault , dtMax=Parameters.dtMax):
##
##    fd_watch = {}
##    fd_watch['thick']  = open( 'thickness_py.bin' , 'wb' )
##
##    conserveIce = numpy.float64(0.)
##    counter = 0
##    # tic = time.time()
##
##    while (t < t_max):
##
##        #-----------------------------------------------------------
##        # COMPRESS - ONLY SIMULATE SUB-RECTANGLE THAT CONTAINS ICE
##        #-----------------------------------------------------------
##        ( Zi , compression_ratio , COMPRESSED_FLAG ) = compress_grid( H , Zb , COMPRESS_TOGGLE=False )
##
##        #-----------------------------------------------------------
##        # MODIFY BOUNDARY CELLS TO ENFORCE BOUNDARY CONDITIONS
##        #-----------------------------------------------------------
##        ( H_ext , Zb_ext , Zi_ext ) = set_bc( H , Zb , Zi )
##        ( dZidxX , dZidyY ) = difference_grid( Zi_ext , dx , dy )
##
##        #-----------------------------------------------------------
##        # CALCULATE THE BASAL SHEAR STRESS
##        #-----------------------------------------------------------
##        ( ( xcmpnt , ycmpnt ) , ( taubX , taubY ) , ( HX , HY ) ) = basal_shear_stress( H_ext , Zi_ext , dx=dx , dy=dy )
##
##        #-----------------------------------------------------------
##        # CALCULATE ICE VELOCITY DUE TO DEFORMATION
##        #-----------------------------------------------------------
##        if ICEFLOW_TOGGLE:
##            ( UdxX , UdyY ) = iceflow( taubX , taubY , HX , HY , xcmpnt , ycmpnt )
##        else:
##            UdxX = numpy.zeros( xcmpnt.shape , dtype='Float64' )
##            UdyY = numpy.zeros( ycmpnt.shape , dtype='Float64' )
##
##        #-----------------------------------------------------------
##        # CALCULATE SLIDING VELOCITY
##        #-----------------------------------------------------------
##        if ICESLIDE_TOGGLE:
##            ( UsxX , UsyY ) = ice_sliding( taubX , taubY , xcmpnt , ycmpnt )
##        else:
##            UsxX = numpy.zeros( xcmpnt.shape , dtype='Float64' )
##            UsyY = numpy.zeros( ycmpnt.shape , dtype='Float64' )
##
##        #-----------------------------------------------------------
##        # Sum all contributions to ice motion
##        #-----------------------------------------------------------
##        ( UxX , UyY ) = sum_ice_motion( UdxX , UdyY , UsxX , UsyY )
##
##        #-----------------------------------------------------------
##        # MASS CONSERVATION -- CONTINUITY
##        #-----------------------------------------------------------
##        ( dHdt , ( qxX , qyY ) ) = mass_conservation( H_ext , UxX , UyY , HX , HY , dZidxX , dZidyY , dx=dx , dy=dy );
##
##        #-----------------------------------------------------------
##        # CALCULATE MASS BALANCE
##        #-----------------------------------------------------------
##        ( Bxy , ELA ) = mass_balance( Zi , t )
##            
##        #-----------------------------------------------------------
##        # CALCULATE TIMESTEP
##        #-----------------------------------------------------------
##        if VARIABLE_DT_TOGGLE:
##            dt = get_timestep( H , Zi_ext , Zi , dHdt , Bxy , dtMax=dtMax , dtDefault=dtDefault )
##        else:
##            dt = dtDefault
##
##        #-----------------------------------------------------------
##        # UPDATE the TIME and ICE THICKNESS
##        #-----------------------------------------------------------
##        ( t, H, Zi, meltrate, conserveIce ) = update_vars( H , Zb , Zi , Bxy ,
##                                                         qxX , qyY , dHdt ,
##                                                         t , dt , conserveIce,
##                                                         dx=dx , dy=dy )
##
##        if (counter%1 == 0):
##            print_watch_point( fd_watch['thick']  , H )
##        counter = counter + 1
##
###   run_for()
#-------------------------------------------------------------------------------------------------- 
##class Usage( Exception ):
##   
##   def __init__( self , msg ):
##      self.msg = msg
##
#-------------------------------------------------------------------------------------------------- 
