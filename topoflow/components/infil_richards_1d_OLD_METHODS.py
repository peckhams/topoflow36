#   Note:  These are old class methods fro infil_richards_1d.py
#          saved here for future reference.  

    #-----------------------------------------------------------------------
    def update_surface_BC_for_theta_old(self, REPORT=False):

        #--------------------------------------------------------
        # Note:  This assumes that psi has already been set at
        #        the surface, and uses TBC relation to compute
        #        theta at the surface.
        #        So must call update_surface_BC_for_psi() 1st.
        #--------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_surface_BC_for_theta()...')

        #------------------------------------
        # Set theta at the surface boundary
        #------------------------------------
        ## print('SINGLE_PROFILE = ' + str(self.SINGLE_PROFILE) )  #######
        if (self.SINGLE_PROFILE):
            psi     = self.p[0]        # [meters]
            theta_s = self.qs[0]
            theta_r = self.qr[0]
            psi_B   = self.pB[0]
            psi_A   = self.pA[0]
            c       = self.c[0]
            Lambda  = self.lam[0]
            self.q[0] = stbc.theta_of_psi(psi, theta_s, theta_r, \
                                  psi_B, psi_A, c, Lambda)
        else:
            #----------------------------------------
            # Now checking if ndim > 1 (2019-10-29)
            #----------------------------------------
            if (self.p.ndim > 1):
                psi = self.p[0,:,:]  # [meters]
            else:
                psi = self.p[0]
            #-----------------------------
            if (self.qs.ndim > 1):
                theta_s = self.qs[0,:,:]
            else:
                theta_s = self.qs[0]
            #-----------------------------
            if (self.qr.ndim > 1):
                theta_r = self.qr[0,:,:]
            else:
                theta_r = self.qr[0]
            #-----------------------------
            if (self.pB.ndim > 1):
                psi_B = self.pB[0,:,:]
            else:
                psi_B = self.pB[0]
            #------------------------------
            if (self.pA.ndim > 1):
                psi_A = self.pA[0,:,:]
            else:
                psi_A = self.pA[0]
            #------------------------------
            if (self.c.ndim > 1):
                c = self.c[0,:,:]
            else:
                c = self.c[0]
            #------------------------------
            if (self.lam.ndim > 1):               
                lam = self.lam[0,:,:]
            else:
                lam = self.lam[0]
            #------------------------------
            self.q[0,:,:] = stbc.theta_of_psi(psi, theta_s, theta_r, \
                                              psi_B, psi_A, c, lam)

        #----------------
        # For debugging
        #----------------
        if (self.DEBUG and self.SINGLE_PROFILE):
            print('In update_surface_BC_for_theta():')
            print('theta[0] =', self.q[0] )
            print('theta[1] =', self.q[1] )

    #   update_surface_BC_for_theta_old()
    #-----------------------------------------------------------------------
    def update_bottom_BC_for_theta_old(self, REPORT=False): 

        if (self.DEBUG):
            print('Calling update_bottom_BC_for_theta()...')

        m = (self.nz - 1)
                            
        #-----------------------------------
        # Set theta at the bottom boundary
        #-----------------------------------
        if (self.SINGLE_PROFILE):
            p   = self.p[m]       # [meters]
            qs  = self.qs[m]
            qr  = self.qr[m]
            pB  = self.pB[m]
            pA  = self.pA[m]
            c   = self.c[m]
            lam = self.lam[m]
            self.q[m] = stbc.theta_of_psi(p, qs, qr, pB, pA, c, lam)
        else:
            #----------------------------------------
            # Now checking if ndim > 1 (2019-10-29)
            #----------------------------------------
            p = self.p[m,:,:]   # [meters]
            #---------------------------------------
            if (self.qs.ndim > 1):
                qs = self.qs[m,:,:]
            else:
                qs = self.qs[m]
            #-----------------------------
            if (self.qr.ndim > 1):
                qr = self.qr[m,:,:]
            else:
                qr = self.qr[m]
            #-----------------------------
            if (self.pB.ndim > 1):
                pB = self.pB[m,:,:]
            else:
                pB = self.pB[m]
            #------------------------------
            if (self.pA.ndim > 1):
                pA = self.pA[m,:,:]
            else:
                pA = self.pA[m]
            #------------------------------
            if (self.c.ndim > 1):
                c = self.c[m,:,:]
            else:
                c = self.c[m]
            #------------------------------
            if (self.lam.ndim > 1):               
                lam = self.lam[m,:,:]
            else:
                lam = self.lam[m]
            #------------------------------            
            self.q[m,:,:] = stbc.theta_of_psi(p, qs, qr, pB, pA, c, lam)

        #----------------
        # For debugging
        #----------------
        if (self.DEBUG and self.SINGLE_PROFILE):
            print('In update_bottom_BC_for_theta():')
            print('theta[m-1] =', self.q[m-1])           
            print('theta[m]   =', self.q[m])
            print(' ')
            
    #   update_bottom_BC_for_theta_old()
    #-----------------------------------------------------------------------
    def update_psi_old(self, REPORT=False):

        #----------------------------------------------------------------
        # Notes: This procedure updates the pressure head, psi, as
        #        a function of the soil moisture, theta, via the
        #        Brooks-Corey (B-C) or "transitional Brooks-Corey"
        #        (TB-C) relation.  The TB-C relation has a continuous
        #        derivative at saturation, unlike the B-C relation.

        #        Psi is < 0 in the unsaturated zone, is 0 at saturation
        #        (e.g. water table) and is > 0 below the water table.

        #        Note that for both B-C and TB-C, psi goes to
        #        -Infinity as theta goes to theta_r (S_eff goes
        #        to zero).  So initial theta values should always
        #        be set to a number greater than theta_r.

        #        For B-C, psi goes to psi_B as theta goes to theta_s,
        #        and psi <= psi_B < 0, or abs(psi) >= abs(psi_B).
        #            pow      = -1d / self.lam
        #            arg      = S_eff^pow
        #            self.psi = self.psiB * arg

        #        For TB-C, psi goes to -psi_a as theta goes to theta_s
        #        and psi <= -psi_a.  If we take psi_a=0, then psi=0 at
        #        saturation (as is commonly assumed).  The hysteresis
        #        effect can be addressed by taking psi_a ne 0 when the
        #        soil is drying (theta decreasing) and perhaps also by
        #        changing the other parameters.

        #        There is a typo in R.E. Smith's AGU monograph in
        #        equation (2.14), where lambda should be -lambda.

        #        See "Infiltration Theory for Hydrologic Applica-
        #        tions" by R.E. Smith (2002), p. 21-22.

        #----------------------------------------------------------------
        # NB!    Due to multiple layers, each input var was set to
        #        a 1D or 3D array by initialize_layer_vars().
        #----------------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_psi()...')
            
        #--------------------------
        # For testing & debugging
        #--------------------------
        #print,'SINGLE_PROFILE  = ', self.SINGLE_PROFILE
        #print,'size(self.q)    = ', self.q.size
        #print,'size(self.p)    = ', self.p.size
        #print,'size(self.qs)   = ', self.qs.size
        #print,'size(self.qr)   = ', self.qr.size
        #print,'size(self.pB)   = ', self.pB.size
        #print,'size(self.pA)   = ', self.pA.size
        #print,'size(self.lam)  = ', self.lam.size
        #print,'size(self.c)    = ', self.c.size
        #print,' '
        
        #---------------------------------------
        # Compute the "effective saturation"
        # Relative saturation = theta/porosity
        #---------------------------------------
        if (self.SINGLE_PROFILE):    
            #--------------------------------
            # All of the vars are 1D arrays
            #--------------------------------------------
            # MIN_VALUE = -150 [m] is wilting point psi
            #--------------------------------------------                                        
            self.p[:] = stbc.psi_of_theta(self.q, self.qs, self.qr,
                                          self.lam, self.c, self.pB, self.pA)
                                          ## MIN_VALUE=-150.0 )
        else:    
            #--------------------------------------
            # Each var is either a 1D or 3D array
            #--------------------------------------
            dim_q   = np.ndim(self.q)
            dim_qs  = np.ndim(self.qs)
            dim_qr  = np.ndim(self.qr)
            dim_pB  = np.ndim(self.pB)
            dim_pA  = np.ndim(self.pA)
            dim_lam = np.ndim(self.lam)
            dim_c   = np.ndim(self.c)
                                                            
            for j in range(self.nz):
                #--------------------------------------------------
                # At a given z, every input var is scalar or grid
                #--------------------------------------------------
                if (dim_q == 3):    
                    q = self.q[j,:,:]
                else:    
                    q = self.q[j]
                if (dim_qs == 3):    
                    qs = self.qs[j,:,:]
                else:    
                    qs = self.qs[j]
                if (dim_qr == 3):    
                    qr = self.qr[j,:,:]
                else:    
                    qr = self.qr[j]
                if (dim_pB == 3):    
                    pB = self.pB[j,:,:]
                else:    
                    pB = self.pB[j]
                if (dim_pA == 3):    
                    pA = self.pA[j,:,:]
                else:    
                    pA = self.pA[j]
                if (dim_lam == 3):    
                    lam = self.lam[j,:,:]
                else:    
                    lam = self.lam[j]
                if (dim_c == 3):    
                    c = self.c[j,:,:]
                else:    
                    c = self.c[j]
                 
                #--------------------------------------------
                # MIN_VALUE = -150 [m] is wilting point psi
                #--------------------------------------------              
                self.p[j,:,:] = stbc.psi_of_theta(q, qs, qr, lam, c, pB, pA ) 
                                                  ### MIN_VALUE=-150.0)
                                               
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('psi   = ', self.p[0:3])
            #print,' '

        if (self.DEBUG):
        ### if (True):
            print('In update_psi():')
            if (self.SINGLE_PROFILE):
                m = (self.nz - 1)
                print('psi[0], theta[0] =', self.p[0], ', ', self.q[0])
                print('psi[m], theta[m] =', self.p[m], ', ', self.q[m])
            #-----------------------------------
            print('min(psi) =', self.p.min() )
            print('max(psi) =', self.p.max() )
            #-----------------------------------
            print('min(q)   =', self.q.min() )
            print('max(q)   =', self.q.max() )
            print('min(qs)  =', self.qs.min() )
            print('max(qs)  =', self.qs.max() )
            print('min(qr)  =', self.qr.min() )
            print('max(qr)  =', self.qr.max() )
#             print('min(Se)  =', S_eff.min()  )
#             print('max(Se)  =', S_eff.max() )
            print('min(pB)  =', self.pB.min() )
            print('max(pB)  =', self.pB.max() )
            print('min(pA)  =', self.pA.min() )
            print('max(pA)  =', self.pA.max() )
            print('min(c)   =', self.c.min() )
            print('max(c )  =', self.c.max() )
            print('min(lam) =', self.lam.min() )
            print('max(lam) =', self.lam.max() )
            print('min(-c/lam) =', np.min( -self.c / self.lam ) )
            print('max(-c/lam) =', np.max( -self.c / self.lam ) )
                      
    #   update_psi_old()
    #-----------------------------------------------------------------------
    def update_K_old(self, REPORT=False):

        #------------------------------------------------------------
        # Notes: This procedure updates the hydraulic conductivity,
        #        K, as a function of the pressure head, psi, via
        #        "transitional Brooks-Corey" (TB-C) relation.
        #        The standard "Brooks-Corey" (B-C) relation is the
        #        special case of c=1, pA=0.

        #        lambda = pore size distribution parameter
        #        eta    = "pore-disconnectedness" parameter
        #        eta    = 2d + (3d * lambda)

        #        There is a typo in R.E. Smith's AGU monograph in
        #        equation (2.14), where eta should be -eta.

        #        See "Infiltration Theory for Hydrologic Applica-
        #        tions" by R.E. Smith (2002), p. 21-22.

        #        For standard Brooks-Corey we would have:
        #            pow = -1d * (eta)
        #            Kr = (p / pB)^pow
        #------------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_K()...')
            
        #-------------------------------------------------
        # Compute K from the "relative conductivity", Kr
        #-------------------------------------------------
        # Use Transitional Brooks-Corey (TB-C)
        # w/ continuous derivative at saturation
        # Note:  q=qs => psi=0, + psiA=0 => K=Ks
        #-----------------------------------------
        if (self.SINGLE_PROFILE):    
            #--------------------------------
            # All of the vars are 1D arrays
            #-----------------------------------------
            # Can there still be multiple layers ???
            #-----------------------------------------            
            self.K[:] = stbc.K_of_psi(self.p, self.pB, self.pA,
                                      self.Ks, self.eta, self.c)
        else:    
            #--------------------------------------
            # Each var is either a 1D or 3D array
            #--------------------------------------
            dim_p   = np.ndim( self.p )
            dim_Ks  = np.ndim( self.Ks )
            dim_pB  = np.ndim( self.pB )
            dim_pA  = np.ndim( self.pA )
            dim_eta = np.ndim( self.eta )
            dim_c   = np.ndim( self.c )

            for j in range(self.nz):
                #--------------------------------------------------
                # At a given z, every input var is scalar or grid
                #--------------------------------------------------
                if (dim_p == 3):    
                    p = self.p[j,:,:]
                else:    
                    p = self.p[j]
                if (dim_Ks == 3):    
                    Ks = self.Ks[j,:,:]
                else:    
                    Ks = self.Ks[j]
                if (dim_pB == 3):    
                    pB = self.pB[j,:,:]
                else:    
                    pB = self.pB[j]
                if (dim_pA == 3):    
                    pA = self.pA[j,:,:]
                else:    
                    pA = self.pA[j]
                if (dim_eta == 3):    
                    eta = self.eta[j,:,:]
                else:    
                    eta = self.eta[j]
                if (dim_c == 3):    
                    c = self.c[j,:,:]
                else:    
                    c = self.c[j]
                #-----------------------------------------------------
                self.K[j,:,:] = stbc.K_of_psi(p, pB, pA, Ks, eta, c)

        if (self.DEBUG):
            print('min(K) =', self.K.min())
            print('max(K) =', self.K.max())
                    
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('K = ', self.K[0:3])
            #print,' '

    #   update_K_old()
