

## Copyright (c) 2001-2013, Scott D. Peckham
## January 2009  (converted from IDL)
## May, August 2009

#-----------------------------------------------------------------------

#  Functions:
#     Theta_Field()
#     Theta_Wilting()
#     Theta_Min()
#     Theta_Residual()
#     Theta_Max()
#     --------------------
#     K_of_Theta()
#     Soil_Types()
#     Get_Soil_Params()

#-----------------------------------------------------------------------
def Theta_Field(theta_s, theta_r, psi_B, psi_A, c,
                _lambda, REPORT=None):

    REPORT = (REPORT not in [0,None])
    
    #-------------------------------------
    #Compute "field capacity" theta value
    #-------------------------------------
    psi_f = -float64(340)           #[cm]
    psi_f = (psi_f / float64(100))  # [cm -> meters]
    
    ratio = (psi_f + psi_A) / psi_B    #(should be > 0)
    
    theta_f = (float64(1) + ratio ** c) ** (-_lambda / c)
    theta_f = theta_f * (theta_s - theta_r) + theta_r
    
    #----------------
    #Optional report
    #----------------
    if (REPORT):    
        print('theta_s = ', theta_s)
        print('theta_f = ', theta_f)
        print('theta_r = ', theta_r)
        print(' ')
    
    return theta_f
    
#   Theta_Field
#-----------------------------------------------------------------------
def Theta_Wilting(theta_s, theta_r, psi_B, psi_A, c,
                  _lambda, REPORT=None):

    REPORT = (REPORT not in [0,None])
    
    #----------------------------------------------
    #Compute "permanent wilting point" theta value
    #----------------------------------------------
    psi_w = -float64(15000)          #[cm]
    psi_w = (psi_w / float64(100))   #[cm -> meters]
    
    ratio = (psi_w + psi_A) / psi_B    #(should be > 0)
    
    theta_w = (float64(1) + ratio ** c) ** (-_lambda / c)
    theta_w = theta_w * (theta_s - theta_r) + theta_r
    
    #----------------
    #Optional report
    #----------------
    if (REPORT):    
        print('theta_s = ', theta_s)
        print('theta_w = ', theta_w)
        print('theta_r = ', theta_r)
        print(' ')
    
    return theta_w
    
#   Theta_Wilting
#-----------------------------------------------------------------------
def Theta_Min(theta_s, theta_r, psi_B, psi_A, c,
              _lambda, REPORT=None):

    #---------------------------------------------------------------
    #Notes:  Note that for both B-C and TB-C, psi goes to
    #        -Infinity as theta goes to theta_r (S_eff goes
    #        to zero).  However, natural soils do not have heads
    #        (tensions) less than -31,000 cm.  In this range they
    #        absorb water from the air (H = hygroscopic).  While
    #        initial theta values will always be set to a number
    #        greater than theta_r, evaporation at the surface can
    #        cause theta to drop to values near theta_r.  Here we
    #        use the T-BC equation for theta(psi) to compute a
    #        value theta_H corresponding to psi_H=-31,000 cm.
    #---------------------------------------------------------------
    REPORT = (REPORT not in [0,None])
    
    #---------------------------
    #Compute min value of theta
    #---------------------------
    psi_H = -float64(31000)         #[cm]
    psi_H = (psi_H / float64(100))  #[cm -> meters]
    
    ratio = (psi_H + psi_A) / psi_B    #(should be > 0)
    
    theta_H = (float64(1) + ratio ** c) ** (-_lambda / c)
    theta_H = theta_H * (theta_s - theta_r) + theta_r
    
    #----------------
    #Optional report
    #----------------
    if (REPORT):    
        print('theta_s = ', theta_s)
        print('theta_H = ', theta_H)
        print('theta_r = ', theta_r)
        print(' ')
    
    return theta_H
    
#   Theta_Min
#-----------------------------------------------------------------------
def Theta_Residual(theta_s, psi_B, _lambda, REPORT=None):

    #------------------------------------------------------------
    #Note:  (3/18/08)  This is based on Equation 6-19 (p. 235),
    #       and Figure 6-13 (p. 237) in Dingman (2002), but
    #       uses a pressure head value of -1000000.  See also
    #       Dingman's Table 6-1 (p. 235).

    #       Recall that lambda = (1 / b).

    #       This is called by Get_Soil_Params in GUI_infil.pro.
    #------------------------------------------------------------
    REPORT = (REPORT not in [0,None])
    
    #-------------------------------
    #Compute "residual" theta value
    #-------------------------------
    psi_r = -float64(1000000)        #[cm]
    psi_r = (psi_r / float64(100))   #[cm -> meters]
    
    #-------------------------------------
    #Note:  Both psi's < 0, so ratio > 0
    #-------------------------------------
    theta_r = theta_s * (psi_B / psi_r) ** _lambda
    
    #----------------
    #Optional report
    #----------------
    if (REPORT):    
        print('psi_B    = ', psi_B)
        print('theta_s  = ', theta_s)
        print('theta_r  = ', theta_r)
        print(' ')
    
    return theta_r
    
#   Theta_Residual
#-----------------------------------------------------------------------
def Theta_Max(r0, K_s, eta, _lambda, theta_s, theta_r, \
              theta_i, REPORT=None):

    #-------------------------------------------------------
    #Note:  The limiting value of theta turns out to be the
    #       same for the Brooks-Corey and the transitional
    #       Brooks-Corey relations.  It is computed by
    #       setting K = r0 and solving for theta.
    #-------------------------------------------------------
    REPORT = (REPORT not in [0,None])
    
    #---------------------------
    #Compute max value of theta
    #---------------------------
    if (r0 < K_s):    
        eps = eta / _lambda
        arg = (r0 / K_s) ** (float64(1) / eps)
        theta_u = theta_r + ((theta_s - theta_r) * arg)
    else:    
        theta_u = theta_s
    
    #----------------
    #Optional report
    #----------------
    if (REPORT):    
        print('theta_s = ', theta_s)
        print('theta_u = ', theta_u)
        print('theta_i = ', theta_i)
        print(' ')
    
    return theta_u
    
#   Theta_Max
#-----------------------------------------------------------------------
def K_of_Theta(theta, K_s, theta_s, theta_r, _lambda, REPORT=None):

    #--------------------------------------------------------------
    #Notes:  This function returns the hydraulic conductivity, K,
    #        as a function of the soil moisture, theta, using an
    #        equation that holds for both the "Brooks-Corey" (B-C)
    #        and "transitional Brooks-Corey" (TB-C) cases.

    #        Called by Get_Soil_Params to compute K_i.

    #        lambda = pore size distribution parameter
    #        eta    = "pore-disconnectedness" parameter
    #        eta    = 2d + (3d * lambda)
    #        eps    = eta/lambda

    #        See "Infiltration Theory for Hydrologic Applica-
    #        tions" by R.E. Smith (2002), p. 19-22.
    #--------------------------------------------------------------
    REPORT = (REPORT not in [0,None])
    
    #----------------------------
    # Compute exponent, epsilon
    #----------------------------
    eta = (float64(2) + (float64(3) * _lambda))
    eps = eta / _lambda
    
    #--------------------------------------
    # Compute the "relative conductivity"
    #--------------------------------------
    K_r = ((theta - theta_r) / (theta_s - theta_r)) ** eps
    
    #-----------------------------
    # Compute K from K_s and K_r
    #-----------------------------
    K_r = maximum((minimum(K_r, float64(1))), float64(0))
    K = K_s * K_r
    
    #------------------
    # Optional report
    #------------------
    if (REPORT):    
        print('K = ', K[0:4])
        #print,' '
    
    return K
    
#   K_of_Theta
#-----------------------------------------------------------------------
def Soil_Types(FORM2=None):

    if (FORM2 in [0,None]):    
        types = array(['Sand', 'Loamy sand', 'Sandy loam', \
                       'Silty loam', 'Loam', 'Sandy clay loam', \
                       'Silty clay loam', 'Clay loam', \
                       'Sandy clay', 'Silty clay', 'Clay'])
    else:    
        types = array(['sand', 'loamy_sand', 'sandy_loam', \
                       'silty_loam', 'loam', 'sandy_clay_loam', \
                       'silty_clay_loam', 'clay_loam', \
                       'sandy_clay', 'silty_clay', 'clay'])
    
    return types
    
#   Soil_Types
#-----------------------------------------------------------------------
def Get_Soil_Params(soil_type, REPORT=None):

    #----------------------------------------------------------
    # Notes:  The values here are from Table 6.1 (p. 235) in
    #         Dingman (2002), 2nd. edition.  Data originally
    #         from Clapp and Hornberger (1978).

    #         Values for G were taken from Table 8.1 (p. 136)
    #         in R.E. Smith's monograph.  Typical G for silt
    #         is given as 914.0 mm, but Dingman's table does
    #         not have an entry for silt.
    #----------------------------------------------------------

    REPORT = (REPORT not in [0,None])

    table = array([array([0.395, 1.76E-2, 12.1, 4.05, 82.0]),   \
                   array([0.410, 1.56E-2,  9.0, 4.38, 97.0]),   \
                   array([0.435, 3.47E-3, 21.8, 4.90, 165.0]),  \
                   array([0.485, 7.20E-4, 78.6, 5.30, 724.0]),  \
                   array([0.451, 6.95E-4, 47.8, 5.39, 385.0]),  \
                   array([0.420, 6.30E-4, 29.9, 7.12, 240.0]),  \
                   array([0.477, 1.70E-4, 35.6, 7.75, 1590.0]), \
                   array([0.476, 2.45E-4, 63.0, 8.52, 804.0]),  \
                   array([0.426, 2.17E-4, 15.3, 10.4, 589.0]),  \
                   array([0.492, 1.03E-4, 49.0, 10.4, 3570.0]), \
                   array([0.482, 1.28E-4, 40.5, 11.4, 2230.0])])    # (clay)

    table = float32(table)

##    table = array([array([float32(0.395), float32(1.76E-2), float32(12.1), float32(4.05), float32(82.0)]), \
##                   array([float32(0.410), float32(1.56E-2), float32(9.0), float32(4.38), float32(97.0)]), \
##                   array([float32(0.435), float32(3.47E-3), float32(21.8), float32(4.90), float32(165.0)]), \
##                   array([float32(0.485), float32(7.20E-4), float32(78.6), float32(5.30), float32(724.0)]), \
##                   array([float32(0.451), float32(6.95E-4), float32(47.8), float32(5.39), float32(385.0)]), \
##                   array([float32(0.420), float32(6.30E-4), float32(29.9), float32(7.12), float32(240.0)]), \
##                   array([float32(0.477), float32(1.70E-4), float32(35.6), float32(7.75), float32(1590.0)]), \
##                   array([float32(0.476), float32(2.45E-4), float32(63.0), float32(8.52), float32(804.0)]), \
##                   array([float32(0.426), float32(2.17E-4), float32(15.3), float32(10.4), float32(589.0)]), \
##                   array([float32(0.492), float32(1.03E-4), float32(49.0), float32(10.4), float32(3570.0)]), \
##                   array([float32(0.482), float32(1.28E-4), float32(40.5), float32(11.4), float32(2230.0)])])
    
    types = Soil_Types(FORM2=True)
    ## print 'types =', types
    ## print 'soil_type.lower() =', soil_type.lower()
    
    w  = where(types == soil_type.lower())
    nw = size(w[0])
    if (nw == 0):    
        print('*******************************')
        print(' Sorry, soil type not found. ')
        print('*******************************')
        print(' ')
        return None
    
    #----------------------------
    # Get parameters from table
    #----------------------------
    row = table[w[0][0],:]
    por = row[0]   #[unitless]
    K_s = row[1]   #[cm/s]
    psi_B = row[2]   #[cm]
    b = row[3]   #[unitless]
    G = row[4]   #[mm]
    
    #----------------
    # Convert units
    #----------------
    K_s   = K_s / float64(100)      #([cm/s] -> [m/s])
    psi_B = psi_B / float64(100)    #([cm]   -> [m])
    G     = G / float64(1000)       #([mm]   -> [m])
    
    #------------------------------------
    # psi_B should be negative, right ?
    #------------------------------------
    psi_B = -float64(1) * psi_B
    
    #------------------------
    # Computable parameters
    #------------------------
    _lambda = (float64(1) / b)
    eta = float64(2) + (float64(3) * _lambda)
    
    #---------------------
    # Arbitrary defaults
    #---------------------
    c = float64(1)
    psi_A = float64(0)
    
    #------------------------------------------
    # Set saturated water content to porosity
    #------------------------------------------
    theta_s = por    # (generally a bit less)
    
    #-----------------------------------------------
    # Set residual water content to Theta_Residual
    #-----------------------------------------------
    theta_r = Theta_Residual(theta_s, psi_B, _lambda)
    #** theta_r = 0.05d    ;(before 3/18/08)
    
    #----------------------------------------------
    # Set initial water content to field capacity
    #----------------------------------------------
    theta_i = Theta_Field(theta_s, theta_r, psi_B, psi_A, c, _lambda)
    #** theta_i = 0.1d    ;(before 3/18/08)
    
    #-----------------------------------------
    # Set these values to ensure stability ?
    # NB!  Larger dz allows larger time step
    #-----------------------------------------
    dz = float64(0.03)    #(= 3 cm)
    nz = int32(20)
    
    #------------------------------------------
    # Compute K_i from theta_i, theta_r, etc.
    #------------------------------------------
    K_i = K_of_Theta(theta_i, K_s, theta_s, theta_r, _lambda)
    
    #------------------
    # Optional report
    #------------------
    if (REPORT):    
        print('K_s      = ', K_s, ' [m/s]')
        print('K_i      = ', K_i, ' [m/s]')
        print('porosity = ', por, ' [unitless]')
        print('theta_s  = ', theta_s, ' [unitless]')
        print('theta_i  = ', theta_i, ' [unitless]')
        print('theta_r  = ', theta_r, ' [unitless]')
        print('psi_B    = ', psi_B, ' [m]')
        print('psi_a    = ', psi_a, ' [m]')
        print('b        = ', b, ' [unitless]')
        print('lambda   = ', _lambda, ' [unitless]')
        print('eta      = ', eta, ' [unitless]')
        print('G        = ', G, ' [m]')
        print(' ')
    
    return por, K_s, psi_B, b, _lambda, eta, theta_s, G, theta_i, \
           K_i, theta_r, c, psi_A, dz, nz

#   Get_Soil_Params
#-------------------------------------------------------------------------------------
#  From Rawls et al. (1983)
#  Repeated as Table 2 in Tarboton_Infiltration_Chapter5.pdf, p. 13.
#-------------------------------------------------------------------------------------
# Soil      |   Porosity, n  |  Effective       |  Wetting         |   Hydraulic
# Texture   |                |  porosity, θe    |  front soil      |   conductivity,
#           |                |                  |  suction head    |   Ksat (cm/hr)
#           |                |                  |  |ψf| (cm)       |
#-------------------------------------------------------------------------------------
# Sand            0.437            0.417                4.95                 11.78 
#                 (0.374-0.500)     (0.354-0.480)        (0.97-25.36) 
# #-------------------------------------------------------------------------------------
# Loamy sand        0.437             0.401                6.13                2.99
#                 (0.363-0.506)     (0.329-0.473)        (1.35-27.94)
# #-------------------------------------------------------------------------------------
# Sandy loam        0.453             0.412                11.01                1.09
#                 (0.351-0.555)     (0.283-0.541)        (2.67-45.47)
# #-------------------------------------------------------------------------------------
# Loam            0.463            0.434                8.89                0.34
#                 (0.375-0.551)    (0.334-0.534)        (1.33-59.38)
# #-------------------------------------------------------------------------------------
# Silt loam        0.501            0.486                16.68                0.65
#                 (0.420-0.582)     (0.394-0.578)        (2.92-95.39)
# #-------------------------------------------------------------------------------------
# Sandy            0.398            0.330                21.85                0.15
# clay loam        (0.332-0.464)    (0.235-0.425)        (4.42-108.0)
# #-------------------------------------------------------------------------------------
# Clay loam        0.464            0.309                20.88                0.1
#                 (0.409-0.519)    (0.279-0.501)        (4.79-91.10)
# #-------------------------------------------------------------------------------------
# Silty clay        0.471             0.432                27.30                0.1
# loam            (0.418-0.524)    (0.347-0.517)        (5.67- 131.50)
# #-------------------------------------------------------------------------------------
# Sandy clay        0.430             0.321                23.90                0.06
#                 (0.370-0.490)     (0.207-0.435)        (4.08-140.2)
# #-------------------------------------------------------------------------------------
# Silty clay        0.479            0.423                29.22                0.05
#                 (0.425-0.533)    (0.334-0.512)        (6.13-139.4)
# #-------------------------------------------------------------------------------------
# Clay            0.475            0.385                31.63                0.03
#                 (0.427-0.523)    (0.269-0.501)        (6.39-156.5)
# #-------------------------------------------------------------------------------------











