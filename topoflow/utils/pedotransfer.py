
# Copyright (c) 2019, Scott D. Peckham
# August 2019
# See:
# https://csdms.colorado.edu/wiki/Model_help:TopoFlow-Soil_Properties_Page
#-------------------------------------------------------------------

#  read_soil_grid_files()

#  wosten_theta_s()
#  wosten_K_s()
#  wosten_alpha()
#  wosten_n()
#  wosten_l()

#  get_wosten_vars()
#  get_tBC_from_vG_vars()
#  save_soil_hydraulic_vars()

#-------------------------------------------------------------------
#  Set up a "tf4" conda environment (TopoFlow 4.0)
#-------------------------------------------------------------------
#  % conda create --name tf4
#  % conda activate tf4
#  % conda install -c conda-forge gdal   (to read geotiff)
#  % conda install dask
#  % conda install -c conda-forge scipy   (for gamma function)
#         (use conda-forge vs. anaconda; broken?)
#  % conda install -c conda-forge pydap

#-------------------------------------------------------------------
import numpy as np
import gdal, osr  ## ogr
from scipy.special import gamma
import glob

# import os.path

#-------------------------------------------------------------------
def read_soil_grid_files( layer=1 ):

    #-------------------------------------------------------------
    # Read soil property data from ISRIC - SoilGrids files,
    # as needed to compute Wosten (1998) pedotransfer functions.
    # SoilGrids files are in GeoTIFF format.
    #-------------------------------------------------------------
    # Another option is to use the rasterio package:
    #
    # import rasterio
    # with rasterio.open('sample.tif') as r:
    #     ar = r.read()    
    #-------------------------------------------------------------
    # Read percent clay
    # Read percent silt
    # Read percent organic matter
    # Read bulk density
    #-------------------------------------------------------------
    # Here we are using 1km grid cell data for the entire
    # country of South Sudan, downloaded from:
    #  https://soilgrids.org/#!/?layer=ORCDRC_M_sl2_250m&vector=1
    # Go to soilgrids.org.
    # Click on the "Download data" bitmap, on the right.
    # Scroll down to the section with Coverage ID (layer 1km)
    #   and Country droplist.  Choose "South Sudan".
    # Choose layers and click "Download" button.
    # Values are available for 7 different depths:
    #   0.0, 0.05, 0.15, 0.30, 0.60, 1.00, 2.00
    #   sl1, sl2,  sl3,  sl4,  sl5,  sl6,  sl7
    #------------------------------------------------------------- 
    layer_str = str(layer)
    file1 = 'CLYPPT_M_sl' + layer_str + '_1km_South Sudan.tiff'
    file2 = 'SLTPPT_M_sl' + layer_str + '_1km_South Sudan.tiff'
    file3 = 'ORCDRC_M_sl' + layer_str + '_1km_South Sudan.tiff'
    file4 = 'BLDFIE_M_sl' + layer_str + '_1km_South Sudan.tiff'    
    #-------------------------------------------------------------  
    f1 = gdal.Open( file1 )
    # print( f1.RasterCount )
    # print( f1.RasterYSize, f1.RasterXsize )
    C  = f1.ReadAsArray()   # (clay fraction, %, byte type)
    f1 = None  # (close f1)
    #-------------------------------------------------------------    
    f2 = gdal.Open( file2 )
    S  = f2.ReadAsArray()   # (silt fraction, %, byte type)
    f2 = None  # (close f2)
    #-------------------------------------------------------------    
    f3 = gdal.Open( file3 )
    OM  = f3.ReadAsArray()  # (org. matter, g/kg, 
    f3 = None  # (close f3)
    # Convert OM to percent for Wosten 1998 PTF. 
    OM = (OM / 1000.0)
    #-------------------------------------------------------------    
    f4 = gdal.Open( file4 )
    D = f4.ReadAsArray()    # (bulk density, kg/m3)
    f4 = None  # (close f4)
    #------------------------------------------------------------- 
    
    return (C, S, OM, D)
    
#   read_soil_grid_files()
#-------------------------------------------------------------------
def wosten_theta_s( C, S, OM, D, topsoil, subsoil ):

    #--------------------------------    
    # From Wosten (1998). R^2 = 76%
    #-------------------------------------------------------------
    # Equations in Wosten (1998) and Wosten et al. (2001) agree.
    #-------------------------------------------------------------
    p1 = 0.7919 + 0.001691*C - 0.29619*D
    p2 = -0.000001491*S^2 + 0.0000821*OM^2  
    p3 = 0.02427 * (1/C) + 0.01113 * (1/S) + 0.01472 * np.log(S)
    p4 = (-0.0000733 * OM * C) - (0.000619 * D * C)
    p5 = (-0.001183 * D * OM) - (0.0001664 * topsoil * S)

    return (p1 + p2 + p3 + p4 + p5)
    
#   wosten_theta_s()
#-------------------------------------------------------------------
def wosten_K_s( C, S, OM, D, topsoil, subsoil ):

    #---------------------------------   
    # From Wosten (1998). R^2 = 19%.
    # K_s^* = ln(K_s), K_s > 0.
    #---------------------------------
    ####### Are units cm/day ??  ##########
    
    ######################################################
    # NOTE!!  In term p5, the coefficient is given by:
    #         0.02986 in Wosten (1998) and:
    #         0.2986 in Wosten et al. (2001).
    ######################################################
    p1 = 7.755 + 0.0352*S + 0.93*topsoil
    p2 = -0.967*D^2 - 0.000484*C^2 - 0.000322*S^2
    p3 = (0.001 / S) - (0.0748 / OM) - 0.643*np.log(S)
    p4 = (-0.01398 * D * C) - (0.1673 * D * OM)
    #----------------
    # Wosten (1998)
    #----------------
    p5 = (0.02986 * topsoil * C) - 0.03305 * topsoil * S
    #-----------------------    
    # Wosten et al. (2001)
    #-----------------------
    ### p5 = (0.2986 * topsoil * C) - 0.03305 * topsoil * S
    
    Ks =  np.exp( p1 + p2 + p3 + p4 + p5 )
    #-----------------------------------------------
    # Convert units from cm/day to m/sec for TF ??
    #-----------------------------------------------
    # Ks = Ks / (100 * 24.0 * 3600.0)   # [m/sec]
      
    return Ks
    
#   wosten_K_s()
#-------------------------------------------------------------------
def wosten_alpha( C, S, OM, D, topsoil, subsoil):

    #---------------------------------   
    # From Wosten (1998). R^2 = 20%.
    # a^* = ln(a), a > 0.
    #---------------------------------
    p1 = -14.96 + 0.03135*C + 0.0351*S + 0.646*OM + 15.29*D
    p2 = -0.192*topsoil - 4.671*D^2 - 0.000781*C^2 - 0.00687*OM^2
    #----------------
    # Wosten (1998)
    #----------------
    p3 = (0.0449 / OM) + 0.0663*np.log(S) + 0.1482*np.log(OM)
    p4 = (-0.04546 * D * S) + (0.4852 * D * OM)
    #-----------------------    
    # Wosten et al. (2001)
    #-----------------------
    ## p3 = (0.449 / OM) + 0.0663*np.log(S) + 0.1482*np.log(OM)
    ## p4 = (-0.4546 * D * S) - (0.4852 * D * OM)
    p5 = 0.00673 * topsoil * C
   
    return np.exp( p1 + p2 + p3 + p4 + p5)
    
#   wosten_alpha()
#-------------------------------------------------------------------
def wosten_n( C, S, OM, D, topsoil, subsoil ):

    #------------------------------------------   
    # From Wosten (1998). R^2 = 54%.
    # n^* = ln(n-1), n > 1.
    # Wosten (1998) assumes that m = 1 - 1/n.
    #-------------------------------------------------------------
    # Equations in Wosten (1998) and Wosten et al. (2001) agree.
    #-------------------------------------------------------------
    p1 = -25.23 - 0.02195*C + 0.0074*S - 0.1940*OM + 45.5*D
    p2 = -7.24*D^2 + 0.0003658*C^2 + 0.002885*OM^2 - (12.81 / D)
    p3 = (-0.1524 / S) - (0.01958 / OM) - 0.2876 * np.log(S)
    p4 = (-0.0709 * np.log(OM)) - (44.6 * np.log(D))
    p5 = (-0.02264 * D * C) + (0.0896 * D * OM) + (0.00718 * topsoil * C)
    
    return (np.exp(p1 + p2 + p3 + p4 + p5) + 1)

#   wosten_n()
#-------------------------------------------------------------------
def wosten_L( C, S, OM, D, topsoil, subsoil ):

    #-----------------------------------------    
    # From Wosten (1998). R^2 = 12%. 
    # L^* = ln[(L+10)/(10-L)], -10 < L < +10.
    # Mualem (1976) says L should be about 0.5 on average.
    # See:  Wosten et al. (2001) for more about L.
    #-------------------------------------------------------------
    # Equations in Wosten (1998) and Wosten et al. (2001) agree.
    #-------------------------------------------------------------
    p1 = 0.0202 + 0.0006193*C^2  - 0.001136*OM^2
    p2 = -0.2316 * np.log(OM) - (0.03544 * D * C)
    p3 = (0.00283 * D * S) + (0.0488 * D * OM)
    
    s1 = (p1 + p2 + p3)
    
    return 10 * (np.exp(s1) - 1)/(np.exp(s1) + 1)

#   wosten_L()
#-------------------------------------------------------------------
def get_wosten_vars(C, S, OM, D, topsoil, subsoil):

    #----------------------------------------------------------
    # Use the Wosten (1998) pedotransfer functions to compute
    # theta_s, K_s, and van Genuchten parameters, then save
    # them to files.
    #----------------------------------------------------------
    theta_s = wosten_theta_s( C, S, OM, D, topsoil, subsoil )
    K_s     = wosten_K_s( C, S, OM, D, topsoil, subsoil )
    alpha   = wosten_alpha( C, S, OM, D, topsoil, subsoil )
    n       = wosten_n( C, S, OM, D, topsoil, subsoil )
    L       = wosten_L( C, S, OM, D, topsoil, subsoil )
    
    return (theta_s, K_s, alpha, n, L)
    
#   get_wosten_vars()
#-------------------------------------------------------------------
def get_tBC_from_vG_vars( alpha, n, L ):

    #--------------------------------------------------------
    # Convert van Genuchten parameters to those of the
    # transitional Brooks-Corey model.  Although K
    # is computed quite differently for the transitional
    # Brooks-Corey and van Genuchten methods, the equations
    # for ψ are the same if we use the formulas below.
    # For more information, see:
    # https://csdms.colorado.edu/wiki/
    #   Model_help:TopoFlow-Soil_Properties_Page   
    #---------------------------------------------------------
    # NOTE: L is often fixed at 1/2, but Wosten lets it vary.
    # See p. 51 in Wosten (1998).
    # Also see p. 1 in Schaap and van Genuchten (2005).
    #---------------------------------------------------------
    # These equations appear to be general:
    #---------------------------------------------------------
    # (1)  m      = (1 - 1/n)
    # (1') n      = 1 / (1 - m)
    # (1") n-1    = m / (1 - m)
    # (2)  eta    = 2 + (3 * lambda)
    #---------------------------------------------------------
    # These equations come from forcing the transitional
    # Brooks-Corey and van Genuchten equations for pressure
    # head (psi) to match exactly (eta is not involved):
    # tBC params = psi_B, c, lambda  (and eta)
    # vG  params = alpha, m, n, L
    #---------------------------------------------------------
    # NOTE:  We only need alpha and n (not L) to set the
    #        transitional Brooks-Corey parameters.  We
    #        cannot make the functional forms match for K.
    #---------------------------------------------------------
    # (3) psi_B    = 1 / alpha_g
    # (4) c        = n
    # (5) lambda   = m * c = m * n    (Note: c/lambda = 1/m)
    #
    # Using (4) and (5) and (1):
    # (6) lambda = m * c = (1 - 1/n)*c = (1 - 1/n)*n = (n-1)
    #
    # Using (6) and (2)
    # eta    = 2 + (3*lambda) = 2 + 3*(n-1) = 3*n - 1
    # eta    = 3/(1-m) - 1 = [3 - (3-m)]/(1-m) = (2+m)/(1-m)
    #---------------------------------------------------------
    # (n > 1) => 0 < m < 1  
    # (n > 1) => (lambda > 0)
    # (n > 1) => (eta > 2)
    #---------------------------------------------------------  
    psi_B  = (1.0 / alpha)
    c      = n
    lam    = (n - 1)
    eta    = 3*n - 1
    #----------------------------------
    # Compute Green-Ampt parameter, G
    #----------------------------------
    G = -psi_B * gamma(1 + 1/c) * gamma((eta-1)/c) / gamma(eta/c)
    
    return ( psi_B, c, lam, eta, G )

#   get_tBC_from_vG_vars()
#-------------------------------------------------------------------
def save_soil_hydraulic_vars( site_prefix, cfg_directory, layer=1):

    (C, S, OM, D) = read_soil_grid_vars( layer=layer )
    
    topsoil = (layer == 1)
    subsoil = not(topsoil)
    
    (theta_s, K_s, alpha, n, L) = get_wosten_vars(C, S, OM, D, topsoil, subsoil)   
    (psi_B, c, lam, eta, G ) = get_tBC_from_vG_vars(alpha, n, L)
 
    prefix = site_prefix
    # prefix = case_prefix
    
    Ks_file   = prefix + '_2D-Ks.bin'
    qs_file   = prefix + '_2D-qs.bin'
    pB_file   = prefix + '_2D-pB.bin'
    c_file    = prefix + '_2D-c.bin'    
    lam_file  = prefix + '_2D-lam.bin'
    G_file    = prefix + '_2D-G.bin'

    Ks_unit = open(Ks_file, 'wb')
    K_s = np.float32( K_s )  
    K_s.tofile( Ks_unit )
    Ks_unit.close()
    #----------------------------------
    qs_unit = open(qs_file, 'wb')
    theta_s = np.float32( theta_s )
    theta_s.tofile( qs_unit )
    qs_unit.close()
    #----------------------------------
    pB_unit = open(pB_file, 'wb')
    psi_B   = np.float32( psi_B )
    psi_B.tofile( pB_unit)
    pB_unit.close()
    #----------------------------------
    c_unit = open(c_file, 'wb')
    c = np.float32( c )
    c.tofile( c_unit )
    c_unit.close()
    #----------------------------------
    lam_unit = open(lam_file, 'wb')
    lam = np.float32( lam )
    lam.tofile( lam_unit)
    lam_unit.close()
    #----------------------------------
    G_unit = open(G_file, 'wb')
    G = np.float32( G )
    G.tofile( G_unit )
    G_unit.close()
    
#   save_soil_hydraulic_vars()   
#-------------------------------------------------------------------
