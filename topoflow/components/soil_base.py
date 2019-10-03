
## Copyright (c) 2001-2013, Scott D. Peckham
## January 2009  (converted from IDL)
## May, August 2009
## May 2010 (changes to unit_test())

#------------------------------------------------------------
#  Notes:  This class is used to return values associated
#          with a "closest standard soil type" selection
#          in TopoFlow.

#          These are used by "write_richards_1d_cfg_file()"
#          in "infil_richards_1D.py"
#------------------------------------------------------------
#
#  class soil_base
#
#      __init__()
#      soil_type_list()
#      K_of_theta()
#      -------------------
#      theta_residual()
#      theta_min()
#      theta_wilting()
#      theta_field()
#      theta_max()
#      -------------------
#      initialize()
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.utils import BMI_base

#-----------------------------------------------------------------------
class soil_base( BMI_base.BMI_component ):

    def __init__(self, soil_type='sandy_loam'):

        BMI_base.BMI_component.__init__(self)
        
        self.soil_type = soil_type
  
    #   __init__()
    #-----------------------------------------------------------------------
    def soil_type_list(self, FORM2=False):

        if not(FORM2):    
            types = np.array(['Sand', 'Loamy sand', 'Sandy loam', \
                           'Silt loam', 'Loam', 'Sandy clay loam', \
                           'Silty clay loam', 'Clay loam', \
                           'Sandy clay', 'Silty clay', 'Clay'])
        else:    
            types = np.array(['sand', 'loamy_sand', 'sandy_loam', \
                           'silt_loam', 'loam', 'sandy_clay_loam', \
                           'silty_clay_loam', 'clay_loam', \
                           'sandy_clay', 'silty_clay', 'clay'])
        
        return types
        
    #   soil_type_list()
    #-----------------------------------------------------------------------
    def K_of_theta(self, theta, REPORT=False):

        #--------------------------------------------------------------
        # Notes: This function returns the hydraulic conductivity, K,
        #        as a function of the soil moisture, theta, using an
        #        equation that holds for both the "Brooks-Corey" (B-C)
        #        and "transitional Brooks-Corey" (TB-C) cases.

        #        Called by Get_Soil_Params to compute K_i.

        #        Lambda = pore size distribution parameter
        #        eta    = "pore-disconnectedness" parameter
        #        eta    = 2d + (3d * Lambda)
        #        eps    = eta/Lambda

        #        See "Infiltration Theory for Hydrologic Applica-
        #        tions" by R.E. Smith (2002), p. 19-22.
        #--------------------------------------------------------------
        
        #----------------------------
        # Compute exponent, epsilon
        #----------------------------
        self.eta = (np.float64(2) + (np.float64(3) * self.Lambda))
        eps      = self.eta / self.Lambda
        
        #--------------------------------------
        # Compute the "relative conductivity"
        #--------------------------------------
        K_r = ((theta - self.theta_r) / (self.theta_s - self.theta_r)) ** eps
        
        #-----------------------------
        # Compute K from K_s and K_r
        #-----------------------------
        K_r = np.maximum( np.minimum(K_r, 1), 0)   # (in [0,1])
        K   = self.K_s * K_r
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('K = ', K[0:4])
            #print,' '
        
        return K
        
    #   K_of_theta
    #-----------------------------------------------------------------------
    def theta_residual(self, REPORT=False):

        #------------------------------------------------------------
        #Note:  (3/18/08)  This is based on Equation 6-19 (p. 235),
        #       and Figure 6-13 (p. 237) in Dingman (2002), but
        #       uses a pressure head value of -1000000.  See also
        #       Dingman's Table 6-1 (p. 235).

        #       Recall that Lambda = (1 / b).

        #       This is called by Get_Soil_Params in GUI_infil.pro.
        #------------------------------------------------------------
        
        #---------------------------------
        # Compute "residual" theta value
        #---------------------------------
        psi_r = -np.float64(1000000)        # [cm]
        psi_r = (psi_r / np.float64(100))   # [cm -> meters]
        
        #--------------------------------------
        # Note:  Both psi's < 0, so ratio > 0
        #--------------------------------------
        theta_r = self.theta_s * (self.psi_B / psi_r) ** self.Lambda
 
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('psi_B    = ', self.psi_B)
            print('theta_s  = ', self.theta_s)
            print('theta_r  = ', self.theta_r)
            print(' ')

        self.psi_r   = psi_r
        self.theta_r = theta_r
        
        return theta_r
        
    #   theta_residual()
    #-----------------------------------------------------------------------
    def theta_min(self, REPORT=False):

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
        
        #-----------------------------
        # Compute min value of theta
        #-----------------------------
        psi_H = -np.float64(31000)         #[cm]
        psi_H = (psi_H / np.float64(100))  #[cm -> meters]
        
        ratio = (psi_H + self.psi_A) / self.psi_B    #(should be > 0)
        
        theta_H = (np.float64(1) + ratio ** self.c) ** (-self.Lambda / self.c)
        theta_H = theta_H * (self.theta_s - self.theta_r) + self.theta_r
  
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('theta_s = ', self.theta_s)
            print('theta_H = ', self.theta_H)
            print('theta_r = ', self.theta_r)
            print(' ')

        self.psi_H   = psi_H
        self.theta_H = theta_H
        
        return theta_H
        
    #   theta_min()
    #-----------------------------------------------------------------------
    def theta_wilting(self, REPORT=False):
        
        #------------------------------------------------
        # Compute "permanent wilting point" theta value
        #------------------------------------------------
        psi_w = -np.float64(15000)          #[cm]
        psi_w = (psi_w / np.float64(100))   #[cm -> meters]
        
        ratio = (psi_w + self.psi_A) / self.psi_B    #(should be > 0)
        
        theta_w = (np.float64(1) + ratio ** self.c) ** (-self.Lambda / self.c)
        theta_w = theta_w * (self.theta_s - self.theta_r) + self.theta_r
 
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('theta_s = ', self.theta_s)
            print('theta_w = ', self.theta_w)
            print('theta_r = ', self.theta_r)
            print(' ')

        self.psi_w   = psi_w
        self.theta_w = theta_w
        
        return self.theta_w
    
    #   theta_wilting()
    #-----------------------------------------------------------------------
    def theta_field(self, REPORT=False):
        
        #---------------------------------------
        # Compute "field capacity" theta value
        #---------------------------------------
        psi_f = -np.float64(340)           #[cm]
        psi_f = (psi_f / np.float64(100))  # [cm -> meters]
        
        ratio = (psi_f + self.psi_A) / self.psi_B    #(should be > 0)
        
        theta_f = (np.float64(1) + ratio ** self.c) ** (-self.Lambda / self.c)
        theta_f = theta_f * (self.theta_s - self.theta_r) + self.theta_r

        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('theta_s = ', self.theta_s)
            print('theta_f = ', self.theta_f)
            print('theta_r = ', self.theta_f)
            print(' ')

        self.psi_f   = psi_f
        self.theta_f = theta_f
        
        return theta_f
    
    #   theta_field()
    #-----------------------------------------------------------------------
    def theta_max(self, r0, REPORT=False):

        #---------------------------------------------------------
        # Note:  The limiting value of theta turns out to be the
        #        same for the Brooks-Corey and the transitional
        #        Brooks-Corey relations.  It is computed by
        #        setting K = r0 and solving for theta.
        #---------------------------------------------------------
        
        #-----------------------------
        # Compute max value of theta
        #-----------------------------
        if (r0 < self.K_s):    
            eps = self.eta / self.Lambda
            arg = (r0 / self.K_s) ** (np.float64(1) / eps)
            theta_u = self.theta_r + ((self.theta_s - self.theta_r) * arg)
        else:    
            theta_u = theta_s
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('theta_s = ', theta_s)
            print('theta_u = ', theta_u)
            print('theta_i = ', theta_i)
            print(' ')

        self.theta_u = theta_u
        
        return theta_u
        
    #   theta_max()
    #-----------------------------------------------------------------------
    def initialize(self, cfg_file=None, REPORT=False):

        #-------------------------------------------------------------
#         types = array(['Sand', 'Loamy sand', 'Sandy loam', \
#                        'Silty loam', 'Loam', 'Sandy clay loam', \
#                        'Silty clay loam', 'Clay loam', \
#                        'Sandy clay', 'Silty clay', 'Clay'])
        #-------------------------------------------------------------
        # Notes:  The values here are from Table 6.1 (p. 235) in
        #         Dingman (2002), 2nd. edition.  Data originally
        #         from Clapp and Hornberger (1978).

        #         Values for G were taken from Table 8.1 (p. 136)
        #         in R.E. Smith's monograph.  Typical G for silt
        #         is given as 914.0 mm, but Dingman's table does
        #         not have an entry for silt.
        #----------------------------------------------------------
        table = np.array([np.array([0.395, 1.76E-2, 12.1, 4.05, 82.0]),   \
                       np.array([0.410, 1.56E-2,  9.0, 4.38, 97.0]),   \
                       np.array([0.435, 3.47E-3, 21.8, 4.90, 165.0]),  \
                       np.array([0.485, 7.20E-4, 78.6, 5.30, 724.0]),  \
                       np.array([0.451, 6.95E-4, 47.8, 5.39, 385.0]),  \
                       np.array([0.420, 6.30E-4, 29.9, 7.12, 240.0]),  \
                       np.array([0.477, 1.70E-4, 35.6, 7.75, 1590.0]), \
                       np.array([0.476, 2.45E-4, 63.0, 8.52, 804.0]),  \
                       np.array([0.426, 2.17E-4, 15.3, 10.4, 589.0]),  \
                       np.array([0.492, 1.03E-4, 49.0, 10.4, 3570.0]), \
                       np.array([0.482, 1.28E-4, 40.5, 11.4, 2230.0])])    # (clay)

        table = np.float32(table)

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
        
        types = self.soil_type_list(FORM2=True)
        ## print 'types =', types
        ## print 'soil_type.lower() =', self.soil_type.lower()
        
        w  = np.where( types == self.soil_type.lower() )
        nw = np.size( w[0] )
        if (nw == 0):    
            print('*******************************')
            print(' Sorry, soil type not found. ')
            print('*******************************')
            print(' ')
            return
        
        #----------------------------
        # Get parameters from table
        #----------------------------
        row = table[ w[0][0],: ]
        self.porosity = row[0]   #[unitless]
        self.K_s      = row[1]   #[cm/s]
        self.psi_B    = row[2]   #[cm]
        self.b        = row[3]   #[unitless]
        self.G        = row[4]   #[mm]
        
        #----------------
        # Convert units
        #----------------
        self.K_s   = self.K_s   / np.float64(100)    #([cm/s] -> [m/s])
        self.psi_B = self.psi_B / np.float64(100)    #([cm]   -> [m])
        self.G     = self.G     / np.float64(1000)   #([mm]   -> [m])
        
        #------------------------------------
        # psi_B should be negative, right ?
        #------------------------------------
        self.psi_B = -np.float64(1) * self.psi_B
        
        #------------------------
        # Computable parameters
        #------------------------
        self.Lambda = (np.float64(1) / self.b)
        self.eta    = np.float64(2) + (np.float64(3) * self.Lambda)
        
        #---------------------
        # Arbitrary defaults
        #---------------------
        self.c     = np.float64(1)
        self.psi_A = np.float64(0)

        #------------------------------------------
        # Set saturated water content to porosity
        #------------------------------------------
        self.theta_s = self.porosity   # (generally a bit less)
        
        #-----------------------------------------------
        # Set residual water content to Theta_Residual
        #-----------------------------------------------
        self.theta_r = self.theta_residual()
        #** theta_r = 0.05d    ;(before 3/18/08)
        
        #----------------------------------------------
        # Set initial water content to field capacity
        #----------------------------------------------
        self.theta_i = self.theta_field()
        #** theta_i = 0.1d    ;(before 3/18/08)

        ##########################################
        r0 = self.K_s * 0.95
        self.theta_H = self.theta_min()
        self.theta_w = self.theta_wilting()
        self.theta_f = self.theta_field()
        self.theta_u = self.theta_max(r0)
        
        #-----------------------------------------
        # Set these values to ensure stability ?
        # NB!  Larger dz allows larger time step
        #-----------------------------------------
        self.dz = np.float64(0.03)    #(= 3 cm)
        self.nz = np.int32(20)
        
        #------------------------------------------
        # Compute K_i from theta_i, theta_r, etc.
        #------------------------------------------
        self.K_i = self.K_of_theta(self.theta_i)
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            print('soil_type =', self.soil_type)
            print('K_s       = ', self.K_s,      ' [m/s]')
            print('K_i       = ', self.K_i,      ' [m/s]')
            print('porosity  = ', self.porosity, ' [none]')
            print('theta_s   = ', self.theta_s,  ' [none]')
            print('theta_i   = ', self.theta_i,  ' [none]')
            print('theta_r   = ', self.theta_r,  ' [none]')
            print('psi_B     = ', self.psi_B,    ' [m]')
            print('psi_A     = ', self.psi_A,    ' [m]')
            print('b         = ', self.b,        ' [none]')
            print('lambda    = ', self.Lambda,   ' [none]')
            print('eta       = ', self.eta,      ' [none]')
            print('G         = ', self.G,        ' [m]')
            print(' ')

    #   initialize()
    #-----------------------------------------------------------------------

