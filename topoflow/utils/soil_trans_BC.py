
#  Copyright (c) 2019, Scott D. Peckham
#  October 2019  (moved from components/infil_richards_1d.py)

#  Note: This file implements soil retention functions for the
#        transitional Brooks-Corey model (Smith, 2002).
#        soil_vanG.py will be an alternative to this.

#-------------------------------------------------------------------

# get_psi_constants()

# theta_of_psi()      # theta = soil_water__volume_fraction
# K_of_theta()        # K = soil_water__hydraulic_conductivity
# psi_of_theta()

# print_suggested_values()

#-----------------------------------------------------------------------

import numpy as np

#-------------------------------------------------------------------
def get_psi_constants(self):       

    #----------------------------------------------
    # Define some constants in a dictionary
    # Adapted from infil_base.py, set_constants()    
    #----------------------------------------------
    # See Figure 6-13, p. 237 in Dingman (2002)
    #----------------------------------------------------
    # Psi_field_capacity is often used for psi_init.
    # See initialize_theta_i() in infil_richards_1D.py.
    #----------------------------------------------------
    constants = {
    'psi_oven_dry'    :  np.float64(-1e8),      # [m], oven dry
    'psi_air_dry'     :  np.float64(-1e4),      # [m], air dry
    'psi_min'         :  np.float64(-1e4),      # [m], air dry
    'psi_hygro'       :  np.float64(-310),      # [m], hygroscopic
    'psi_wilt'        :  np.float64(-150),      # [m], wilting pt.
    'psi_field'       :  np.float64(-3.4),      # [m], field cap.
    #---------------------------------------------------------------
    'psi_oven_dry_cm' :  np.float64(-1e10),     # [cm], oven dry
    'psi_air_dry_cm'  :  np.float64(-1e6),      # [cm], air dry
    'psi_min_cm'      :  np.float64(-1e6),      # [cm], air dry
    'psi_hygro_cm'    :  np.float64(-31000),    # [cm], hygroscopic
    'psi_wilt_cm'     :  np.float64(-15000),    # [cm], wilting pt.
    'psi_field_cm'    :  np.float64(-340),      # [cm], field cap.
    #---------------------------------------------------------------
    'g'               :  np.float64(9.81) }     # [m s-2], grav. const   
    
    return constants

#   get_psi_constants()
#-----------------------------------------------------------------------
def theta_of_psi(psi, theta_s, theta_r, psi_B, psi_A, c, Lambda,
                 REPORT=False, CM_TO_M=False):

    #---------------------------------------------------------------
    # Notes: This function computes the soil water content, theta
    #        for the give value of pressure head, psi (in cm),
    #        using the soil characteristic relation called
    #        "transitional Brooks-Corey" (TBC).
    #
    #        psi = -1000000 => theta = theta_min (air dry)
    #        psi = -31000   => theta = theta_H (hygroscopic)
    #        psi = -15000   => theta = theta_w (perm. wilting pt.)
    #        psi = -340     => theta = theta_f (field capacity)
    #
    #---------------------------------------------------------------
    # Notes: Note that for both B-C and TB-C, psi goes to
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
    
    #--------------------------------------
    # Convert psi units from cm to meters
    #--------------------------------------
    if (CM_TO_M):
        psi_m = (psi / np.float64(100))    # [cm -> meters]
        ratio = (psi_m + psi_A) / psi_B    # (should be > 0)
    else:
        ratio = (psi + psi_A) / psi_B      # (should be > 0)

#     print('min(ratio)   = ' + str(ratio.min()) )
#     print('max(ratio)   = ' + str(ratio.max()) )
#     print('shape(ratio) = ' + str(ratio.shape) )
#     print()
    #------------------------------------------------
#     print('min(c)     = ' + str(c.min()) )
#     print('max(c)     = ' + str(c.max()) )
#     print('min(lam)   = ' + str(Lambda.min()) )
#     print('max(lam)   = ' + str(Lambda.max()) )
#     print()
         
    cap_theta = (1.0 + (ratio ** c)) ** (-Lambda / c)   # S_eff
    theta = cap_theta * (theta_s - theta_r) + theta_r
    
    #------------------
    # Optional report
    #------------------
    if (REPORT):   
        print('theta_s = ', theta_s)
        print('theta   = ', theta)
        print('theta_r = ', theta_r)
        print(' ')
    
    return theta
    
#   theta_of_psi()
#-----------------------------------------------------------------------
def K_of_theta(theta, K_s, theta_s, theta_r, Lambda,
               REPORT=False):

    #--------------------------------------------------------------
    # Notes: This function returns the hydraulic conductivity, K,
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
    
    #----------------------------
    # Compute exponent, epsilon
    #----------------------------
    eta = (np.float64(2) + (np.float64(3) * Lambda))
    eps = eta / Lambda
    
    #--------------------------------------
    # Compute the "relative conductivity"
    #--------------------------------------
    K_r = ((theta - theta_r) / (theta_s - theta_r)) ** eps
    
    #-----------------------------
    # Compute K from K_s and K_r
    #-----------------------------
    K_r = np.maximum( np.minimum(K_r, 1.0), 0.0 )
    K = K_s * K_r
    
    #------------------
    # Optional report
    #------------------
    if (REPORT):    
        print('K = ', K[0:4])
        # print ' '
    
    return K
    
#   K_of_theta()
#-----------------------------------------------------------------------
def psi_of_theta(theta, K_s, theta_s, theta_r, Lambda, c,
                 psi_B, psi_A, REPORT=False):

    #----------------------------------------
    # S_eff = effective saturation in [0,1]
    #----------------------------------------
    S_eff = (theta - theta_r) / (theta_s - theta_r)
    cpow  = (-c / Lambda)
    arg   = ((S_eff ** cpow )- 1.0) ** (1.0 / c)
    psi   = (psi_B * arg) - psi_A
   
    return psi
    
#   psi_of_theta()         
#-----------------------------------------------------------------------
# def print_suggested_values():
# 
#     ###############################
#     #  NOT READY YET
#     ###############################    
#     if (self.DEBUG):
#         print('Calling print_suggested_values()...')
#         
#     #-------------------------------------------------
#     # theta_r is often set to theta_hygroscopic.
#     # theta_i is often set to theta_field_capacity.
#     #-------------------------------------------------
#     print('=====================================================')
#     for k in range(self.n_layers):
# 
# ##            print 'Ks[k]  =', self.Ks_list[k]
# ##            print 'Ki[k]  =', self.Ki_list[k]    
# ##            print 'qs[k]  =', self.qs_list[k]
# ##            print 'qi[k]  =', self.qi_list[k]
# ##            print 'qr[k]  =', self.qr_list[k]
# ##            print 'pB[k]  =', self.pB_list[k]
# ##            print 'pA[k]  =', self.pA_list[k]
# ##            print 'lam[k] =', self.lam_list[k]
# ##            print 'c[k]   =', self.c_list[k]
# ##            print ' '
# ##            print 'psi_hygro =', self.psi_hygro, ' [m]'
# ##            print 'psi_field =', self.psi_field, ' [m]'
#         
#         #-------------------------------------------------
#         # Compute this by analogy to equations 6-19 and
#         # 6-20 in Dingman (2002), using theta_s instead
#         # of porosity and recalling lambda = (1/b).
#         #-------------------------------------------------
#         # Note that this is not entirely consistent with
#         # the Theta_TBC() function, but that function
#         # requires theta_res as an argument.
#         #-------------------------------------------------
# ##            psi_res   = self.psi_hygro
# ##            theta_sat = self.qs_list[k]
# ##            psi_B     = self.pB_list[k]
# ##            lam       = self.lam_list[k]
# ##            #--------------------------------------
# ##            # Note:  Both psi's < 0, so ratio > 0
# ##            #--------------------------------------
# ##            theta_res = theta_sat * (psi_B / psi_res)**lam
# 
#         #--------------------------------------------
#         # If we trust theta_r, then do this instead
#         #--------------------------------------------
#         theta_res = self.qr_list[k]
#         ###############################################
#         
#         
#         theta_hygro = theta_of_psi( self.psi_hygro,
#                                  self.qs_list[k],
#                                  self.qr_list[k],
#                                  self.pB_list[k],
#                                  self.pA_list[k],
#                                  self.c_list[k],
#                                  self.lam_list[k] )
#         
#         theta_init = theta_of_psi( self.psi_field,
#                                 self.qs_list[k],
#                                 theta_res,         #######
#                                 self.pB_list[k],
#                                 self.pA_list[k],
#                                 self.c_list[k],
#                                 self.lam_list[k] )
# 
#         K_init = K_of_theta( theta_init,       #######
#                              self.Ks_list[k],
#                              self.qs_list[k],
#                              theta_res,        #######
#                              self.lam_list[k] )
# 
#         theta_r = self.qr_list[k]
#         theta_i = self.qi_list[k]
#         K_i     = self.Ki_list[k]
#         print('Suggested initial values for layer', k+1, ':')
#         ## print '   theta_r =', theta_res,  'vs.', theta_r
#         print('   For theta_r =', theta_r)
#         print('   theta_i =', theta_init, '   vs.', theta_i)
#         print('   K_i     =', K_init,     'vs.', K_i)
#         print('   theta_H =', theta_hygro, '  vs.', theta_r, ' (theta_r)')
#         print(' ')
#         
#     print('===========================================================')                                                
#     
# #   print_suggested_values()
#-------------------------------------------------------------------
    
    

