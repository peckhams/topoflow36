
## Copyright (c) 2001-2010, Scott D. Peckham
## January 2009  (converted from IDL)
## May, August 2009
## May 2010  (changes to unit_test() and read_cfg_file()
## June 2010 (Bug fix: Added qH_val and eta_val in
##            set_computed_input_vars(). Unit test. )

#---------------------------------------------------------------------
#  NOTES:  This file defines a Richards 1D infiltration component
#          and related functions.  It inherits from the infiltration
#          "base class" in "infil_base.py".
#---------------------------------------------------------------------
#
#  unit_test()
#
#  class infil_richards_1D
#
#      initialize_layer_vars()
#      ----------------------------
#      get_gui_info()
#      get_cfg_extension()
#      read_cfg_file()
#      set_computed_input_vars()
#      check_input_types()
#      initialize_richards_vars()
#      ----------------------------
#      initialize_theta_r()
#      initialize_theta_i()
#      initialize_K_i()
#      ----------------------------
#      update()
#      update_infil_rate()
#      update_Rg()
#      update_q0()
#      update_theta()
#      check_theta()    # (6/29/10)
#      update_psi()
#      update_K()
#      update_v()
#      update_Zw()
#      ----------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      ------------------------------
#      write_richards_1d_cfg_file()     # (obsolete now)
#      ------------------------------
#      build_layered_var()   ########

#  Functions:
#      Theta_TBC()          (still used)
#      K_of_Theta_TBC()     (used by initialize_K_i())
#      Z_Derivative_1D()    (Mar 2008)
#      Z_Derivative_3D()    (Mar 2007)
#
#-----------------------------------------------------------------------

import numpy as np
# import os

from topoflow.components import infil_base
from topoflow.components import soil_base

from topoflow.utils import model_input
from topoflow.utils import tf_utils

################################
# from save_load import *   # (used by read_cfg_file()
################################

# import matplotlib.pyplot   # (not yet available on beach)

#-----------------------------------------------------------------------
def unit_test():

    ic = infil_component()
    ic.CCA   = False
    ic.DEBUG = True
    ## ic.DEBUG = False

    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    tf_utils.TF_Set_Test_Info( ic )
    
    #-------------------------------
    # Initialize and 1 update call
    #-------------------------------
##    print 'STATUS =', ic.get_status()
##    ic.initialize( mode="driver" )
##    print 'STATUS =', ic.get_status()
##    time_sec = float64(0)
##    ic.update(time_sec)
##    print 'STATUS =', ic.get_status()

    #--------------------------------
    # Run model in stand-alone mode
    #--------------------------------
    ic.run_model( cfg_directory=ic.cfg_directory,
                  cfg_prefix=ic.cfg_prefix )

#   unit_test()
#-----------------------------------------------------------------------
class infil_component(infil_base.infil_component):

    #-------------------------------------------------------------------
    def initialize_layer_vars(self):

        #-------------------------------------------------------
        # Notes: We need to call initialize_layer_vars()
        #        before initialize_config_vars(), which may
        #        call read_cfg_file().  However, this means
        #        we haven't read "n_layers" yet, so just
        #        hardwire it here for now. (5/11/10)
        #-------------------------------------------------------
        n_layers = 3
        # n_layers = self.n_layers
        
        #-------------------------------------------------
        # Get arrays to store soil params for each layer
        #-------------------------------------------------
        self.soil_type = np.zeros(n_layers, dtype='<U200')
        self.dz_val    = np.zeros(n_layers, dtype='Float64')    #### + dz3
        self.nz_val    = np.zeros(n_layers, dtype='Int16')      #### + nz3
        #--------------------------------------------------------
        self.Ks_type   = np.zeros(n_layers, dtype='<U200')
        self.Ki_type   = np.zeros(n_layers, dtype='<U200')
        self.qs_type   = np.zeros(n_layers, dtype='<U200')
        self.qi_type   = np.zeros(n_layers, dtype='<U200')
        self.qr_type   = np.zeros(n_layers, dtype='<U200')
        self.pB_type   = np.zeros(n_layers, dtype='<U200')    
        self.pA_type   = np.zeros(n_layers, dtype='<U200')
        self.lam_type  = np.zeros(n_layers, dtype='<U200')
        self.c_type    = np.zeros(n_layers, dtype='<U200')
        #--------------------------------------------------------        
        self.Ks_file  = np.zeros(n_layers, dtype='<U200')
        self.Ki_file  = np.zeros(n_layers, dtype='<U200')
        self.qs_file  = np.zeros(n_layers, dtype='<U200')
        self.qi_file  = np.zeros(n_layers, dtype='<U200')
        self.qr_file  = np.zeros(n_layers, dtype='<U200')
        self.pB_file  = np.zeros(n_layers, dtype='<U200')
        self.pA_file  = np.zeros(n_layers, dtype='<U200')
        self.lam_file = np.zeros(n_layers, dtype='<U200')
        self.c_file   = np.zeros(n_layers, dtype='<U200')
        #---------------------------------------------------------
        # Note: self.Ks is a Python list.  Initially, each entry
        # is a numpy scalar (type 'np.float64').  However, we
        # can later change any list entry to a scalar or grid
        # (type 'np.ndarray'), according to its "Ks_type".
        #---------------------------------------------------------
        # (5/19/10) Seems we need Ks_val vs Ks here, since
        # we use these to build one big, 3D Ks array.
        #---------------------------------------------------------        
        self.Ks_val  = list(np.zeros(n_layers, dtype='Float64'))
        self.Ki_val  = list(np.zeros(n_layers, dtype='Float64'))
        self.qs_val  = list(np.zeros(n_layers, dtype='Float64'))
        self.qi_val  = list(np.zeros(n_layers, dtype='Float64'))
        self.qr_val  = list(np.zeros(n_layers, dtype='Float64'))
        self.pB_val  = list(np.zeros(n_layers, dtype='Float64'))
        self.pA_val  = list(np.zeros(n_layers, dtype='Float64'))
        self.lam_val = list(np.zeros(n_layers, dtype='Float64'))
        self.c_val   = list(np.zeros(n_layers, dtype='Float64'))
        #------------------------------------------------
        # Note:  These two are computed from the others
        #------------------------------------------------
        self.eta_val = list(np.zeros(n_layers, dtype='Float64'))
        self.qH_val  = list(np.zeros(n_layers, dtype='Float64'))
       
    #   initialize_layer_vars()
    #-------------------------------------------------------------------
    def get_gui_info(self):

        ## num_str = str(tf_utils.TF_Version_Number())
        directory = "/data/progs/topoflow/3.1/gui_info/"
        file_name = "Infil_Richards_1D.cfg"
        self.gui_info_file = (directory + file_name)
        # self.cfg_file = (directory + file_name)
        self.dialog_title  = "Infiltration: Richards-1D Parameters"

    #   get_gui_info()
    #-------------------------------------------------------------------
    def get_cfg_extension(self):

        return '_infil_richards_1d.cfg'
    
    #   get_cfg_extension()
    #-------------------------------------------------------------------
##    def read_cfg_file(self):
##
##        #----------------------------------------------------------
##        # Note: Currently, the configuration files for Green-Ampt
##        #       and Smith-Parlange are very different than those
##        #       for the Richards' 1D method.
##        #----------------------------------------------------------
##        
##        #------------------------------------------
##        # Read parameters from a CFG file that is
##        # in the current working directory.
##        #------------------------------------------
##        print 'Infiltration component (Richards 1D): Reading config file...'
##        file_unit = open(self.cfg_file, 'r')
##
##        #-----------------------------
##        # Skip over the header lines
##        #-----------------------------
##        for k in xrange(1, 5):
##            line = file_unit.readline()
##            
##        #----------------------
##        # Read the infil vars
##        #----------------------
##        method       = Read_Vars(file_unit, data_type='BYTE')
##        method_name  = Read_Vars(file_unit, data_type='STRING')
##        n_layers     = Read_Vars(file_unit, data_type='INTEGER')
##        dt_type, dt  = Read_Vars(file_unit, True)
##        #---------------------------------------------------------
##        self.method      = method
##        self.method_name = method_name
##        self.n_layers    = n_layers
##        self.dt          = float64(dt)   # [seconds]
##        
##        #-------------------------------------
##        # Vars for Richards' equation method
##        #-------------------------------------    
##        for k in xrange(n_layers):
##            Ks_type, Ks   = Read_Vars(file_unit, True)
##            Ki_type, Ki   = Read_Vars(file_unit, True)
##            qs_type, qs   = Read_Vars(file_unit, True)
##            qi_type, qi   = Read_Vars(file_unit, True)
##            qr_type, qr   = Read_Vars(file_unit, True)   #**** added 3/14/08
##            pB_type, pB   = Read_Vars(file_unit, True)
##            pA_type, pA   = Read_Vars(file_unit, True)
##            lam_type, lam = Read_Vars(file_unit, True)
##            c_type, c     = Read_Vars(file_unit, True)
##            dz_type, dz   = Read_Vars(file_unit, True)   ####
##            nz_type, nz   = Read_Vars(file_unit, True)   ####
##            soil_type     = Read_Vars(file_unit, \
##                                      data_type='STRING')
##            #--------------------------------------------------
##            Ks_val, Ks_file  = Load_Var(Ks, Ks_type)
##            self.Ks_val[k]   = Ks_val
##            self.Ks_type[k]  = Type_Code(Ks_type)
##            self.Ks_file[k]  = Ks_file
##            #--------------------------------------------------
##            Ki_val, Ki_file  = Load_Var(Ki, Ki_type)
##            self.Ki_val[k]   = Ki_val
##            self.Ki_type[k]  = Type_Code(Ki_type)
##            self.Ki_file[k]  = Ki_file
##            #--------------------------------------------------
##            qs_val, qs_file  = Load_Var(qs, qs_type)
##            self.qs_val[k]   = qs_val
##            self.qs_type[k]  = Type_Code(qs_type)
##            self.qs_file[k]  = qs_file
##            #--------------------------------------------------
##            qi_val, qi_file  = Load_Var(qi, qi_type)
##            self.qi_val[k]   = qi_val
##            self.qi_type[k]  = Type_Code(qi_type)
##            self.qi_file[k]  = qi_file
##            #--------------------------------------------------
##            qr_val, qr_file  = Load_Var(qr, qr_type)
##            self.qr_val[k]   = qr_val   #**** added 3/14/08
##            self.qr_type[k]  = Type_Code(qr_type)
##            self.qr_file[k]  = qr_file
##            #--------------------------------------------------
##            pB_val, pB_file  = Load_Var(pB, pB_type)
##            self.pB_val[k]   = pB_val
##            self.pB_type[k]  = Type_Code(pB_type)
##            self.pB_file[k]  = pB_file
##            #--------------------------------------------------
##            pA_val, pA_file  = Load_Var(pA, pA_type)
##            self.pA_val[k]   = pA_val
##            self.pA_type[k]  = Type_Code(pA_type)
##            self.pA_file[k]  = pA_file
##            #--------------------------------------------------
##            lam_val, lam_file = Load_Var(lam, lam_type)
##            self.lam_val[k]   = lam_val
##            self.lam_type[k]  = Type_Code(lam_type)
##            self.lam_file[k]  = lam_file
##            #--------------------------------------------------
##            c_val, c_file    = Load_Var(c, c_type)
##            self.c_val[k]    = c_val
##            self.c_type[k]   = Type_Code(c_type)
##            self.c_file[k]   = c_file
##            #--------------------------------------------------
##            #--------------------------------------------------
##            dz_val, dz_file  = Load_Var(dz, dz_type)
##            self.dz_val[k]   = dz_val
##            #--------------------------------------------------
##            nz_val, nz_file  = Load_Var(nz, nz_type)
##            self.nz_val[k]   = int16(nz_val)
##            #--------------------------------------------------
####                self.dz_val[k]    = float64(dz)
####                self.nz_val[k]    = int16(nz)
##            #--------------------------------------------------
##            self.soil_type[k] = soil_type
##     
##        #---------------------
##        # Get output options
##        #---------------------
##        dum_str, save_grid_dt     = Read_Vars(file_unit, True)     
##        save_v0_grids, v0_gs_file = Read_Vars(file_unit, True, data_type='BYTE')
##        save_q0_grids, q0_gs_file = Read_Vars(file_unit, True, data_type='BYTE')
##        save_I_grids, I_gs_file   = Read_Vars(file_unit, True, data_type='BYTE')
##        save_Zw_grids, Zw_gs_file = Read_Vars(file_unit, True, data_type='BYTE')
##        #---------------------------------------------------------------------------
##        self.save_grid_dt  = float64(save_grid_dt)
##        self.SAVE_V0_GRIDS = save_v0_grids
##        self.SAVE_Q0_GRIDS = save_q0_grids
##        self.SAVE_I_GRIDS  = save_I_grids
##        self.SAVE_ZW_GRIDS = save_Zw_grids
##        self.v0_gs_file    = v0_gs_file
##        self.q0_gs_file    = q0_gs_file
##        self.I_gs_file     = I_gs_file
##        self.Zw_gs_file    = Zw_gs_file
##        #---------------------------------------------------------------------------
##        dum_str, save_pixels_dt    = Read_Vars(file_unit, True)
##        save_v0_pixels, v0_ts_file = Read_Vars(file_unit, True, data_type='BYTE')
##        save_q0_pixels, q0_ts_file = Read_Vars(file_unit, True, data_type='BYTE')
##        save_I_pixels, I_ts_file   = Read_Vars(file_unit, True, data_type='BYTE')
##        save_Zw_pixels, Zw_ts_file = Read_Vars(file_unit, True, data_type='BYTE')
##        #---------------------------------------------------------------------------
##        self.save_pixels_dt = float64(save_pixels_dt)
##        self.SAVE_V0_PIXELS = save_v0_pixels
##        self.SAVE_Q0_PIXELS = save_q0_pixels
##        self.SAVE_I_PIXELS  = save_I_pixels
##        self.SAVE_ZW_PIXELS = save_Zw_pixels
##        self.v0_ts_file     = v0_ts_file
##        self.q0_ts_file     = q0_ts_file
##        self.I_ts_file      = I_ts_file
##        self.Zw_ts_file     = Zw_ts_file
##        #---------------------------------------------------------------------------
##        dum_str, save_cube_dt   = Read_Vars(file_unit, True)
##        save_q_cubes, q_cs_file = Read_Vars(file_unit, True, data_type='BYTE')
##        save_p_cubes, p_cs_file = Read_Vars(file_unit, True, data_type='BYTE')
##        save_K_cubes, K_cs_file = Read_Vars(file_unit, True, data_type='BYTE')
##        save_v_cubes, v_cs_file = Read_Vars(file_unit, True, data_type='BYTE')
##        #---------------------------------------------------------------------------
##        self.save_cube_dt = float64(save_cube_dt)
##        self.SAVE_Q_CUBES = save_q_cubes
##        self.SAVE_P_CUBES = save_p_cubes
##        self.SAVE_K_CUBES = save_K_cubes
##        self.SAVE_V_CUBES = save_v_cubes
##        self.q_cs_file  = q_cs_file
##        self.p_cs_file  = p_cs_file
##        self.K_cs_file  = K_cs_file
##        self.v_cs_file  = v_cs_file
##        #------------------------------------------------------------------------------
##        dum_str, save_profile_dt   = Read_Vars(file_unit, True)
##        save_q_profiles, q_ps_file = Read_Vars(file_unit, True, data_type='BYTE')
##        save_p_profiles, p_ps_file = Read_Vars(file_unit, True, data_type='BYTE')
##        save_K_profiles, K_ps_file = Read_Vars(file_unit, True, data_type='BYTE')
##        save_v_profiles, v_ps_file = Read_Vars(file_unit, True, data_type='BYTE')
##        #------------------------------------------------------------------------------
##        self.save_profile_dt = float64(save_profile_dt)
##        self.SAVE_Q_PROFILES = save_q_profiles
##        self.SAVE_P_PROFILES = save_p_profiles
##        self.SAVE_K_PROFILES = save_K_profiles
##        self.SAVE_V_PROFILES = save_v_profiles
##        self.q_ps_file  = q_ps_file
##        self.p_ps_file  = p_ps_file
##        self.K_ps_file  = K_ps_file
##        self.v_ps_file  = v_ps_file
##
##        #-----------------------
##        # Close the config file
##        #-----------------------
##        file_unit.close()
##
##        #---------------------------------------------------------
##        # Make sure that all "save_dts" are larger or equal to
##        # the specified process dt.  There is no point in saving
##        # results more often than they change.
##        # Issue a message to this effect if any are smaller ??
##        #---------------------------------------------------------
##        self.save_grid_dt    = maximum(self.save_grid_dt,    self.dt)
##        self.save_pixels_dt  = maximum(self.save_pixels_dt,  self.dt)
##        self.save_cube_dt    = maximum(self.save_cube_dt,    self.dt)
##        self.save_profile_dt = maximum(self.save_profile_dt, self.dt)
##        
##    #   read_cfg_file()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        self.RICHARDS = True
        self.G_file    = ''   # (still need to be defined)
        self.gam_file  = ''
        self.G_type    = 'Scalar'     # (see below)
        self.gam_type  = 'Scalar'

        #------------------------------------------------------------
        # Compute eta value for each soil layer from lambda values.
        # Depending on lambda, eta values will scalars or grids.
        #------------------------------------------------------------
        for j in range(self.n_layers):
            self.eta_val[j] = float64(2) + (float64(3) * self.lam_val[j])
                             
        #--------------------------------------------------------------
        # Compute a qH value for each soil layer from other values
        # using the Theta_TBC() function.  qH values will be scalars
        # or grids, depending on the args to Theta_TBC().
        #-------------------------------------------------------------
        for j in range(self.n_layers):
            self.qH_val[j] = Theta_TBC( self.psi_hygro_cm, \
                                        self.qs_val[j], self.qr_val[j], \
                                        self.pB_val[j], self.pA_val[j], \
                                        self.c_val[j],  self.lam_val[j] )

        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt    = maximum(self.save_grid_dt,    self.dt)
        self.save_pixels_dt  = maximum(self.save_pixels_dt,  self.dt)
        self.save_profile_dt = maximum(self.save_profile_dt, self.dt)
        self.save_cube_dt    = maximum(self.save_cube_dt,    self.dt)
        
    #   set_computed_input_vars()    
    #-------------------------------------------------------------------
    def check_input_types(self):

        #----------------------------------------------------
        # Notes: ET is often a 2D grid even when the others
        #        are scalars.  See how P_total is defined
        #        in update_surface_influx().
        #----------------------------------------------------
        are_scalars = np.array([
                         self.mp.is_scalar('P'),
                         self.sp.is_scalar('SM'),
                         self.ep.is_scalar('ET'),  #########
                         #-------------------------------
                         self.is_scalar('Ks_val[0]'),
                         self.is_scalar('Ki_val[0]'),
                         self.is_scalar('qs_val[0]'),
                         self.is_scalar('qi_val[0]'),
                         self.is_scalar('qr_val[0]'),
                         self.is_scalar('pB_val[0]'),
                         self.is_scalar('pA_val[0]'),
                         self.is_scalar('c_val[0]'),
                         self.is_scalar('lam_val[0]')])

        self.ALL_SCALARS = np.all(are_scalars)
        
        #----------------------------------------
        # Use the same profile for all pixels ?
        #---------------------------------------------
        # NB! This var only used by Richards' method
        # so it shouldn't appear in "infil_base.py".
        #---------------------------------------------
        self.SINGLE_PROFILE = self.ALL_SCALARS  # (3/19/07)

##        print "#### self.mp.is_scalar('P')  =", self.mp.is_scalar('P')
##        print "#### self.sp.is_scalar('SM') =", self.sp.is_scalar('SM')
##        print "#### self.ep.is_scalar('ET') =", self.ep.is_scalar('ET')    
##        print '#### In check_input_types(), ALL_SCALARS =', self.ALL_SCALARS
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def initialize_richards_vars(self):

        ########################################################
        #  NB! The "P" synonym for "rate" doesn't work here.
        #      Maybe defined in the wrong place ??
        ########################################################
        dtype = 'Float64'
        
        #---------------------------------------
        # Get surface influx to initialize "v"
        #---------------------------------------
        self.update_surface_influx()
        
        #----------------------
        # Compute "total nz"
        #---------------------
        self.nz = sum(self.nz_val)

        #------------------------------------------------
        # Now build a 1D or 3D array for each input var
        #--------------------------------------------------------
        # (3/12/08) Same code should work if (self.n_layers eq 1)
        #--------------------------------------------------------
        self.Ks  = self.build_layered_var(self.Ks_val)
        self.Ki  = self.build_layered_var(self.Ki_val)
        self.qs  = self.build_layered_var(self.qs_val)
        self.qi  = self.build_layered_var(self.qi_val)
        self.qr  = self.build_layered_var(self.qr_val)
        self.pB  = self.build_layered_var(self.pB_val)
        self.pA  = self.build_layered_var(self.pA_val)
        self.lam = self.build_layered_var(self.lam_val)
        self.c   = self.build_layered_var(self.c_val)
        #--------------------------------------------------
        # Note:  eta_val and qH_val are computed from
        #        the others in set_computed_input_vars().      
        #--------------------------------------------------
        self.eta = self.build_layered_var(self.eta_val)
        self.qH  = self.build_layered_var(self.qH_val)

        #--------------
        # For testing
        #--------------
        if (self.DEBUG):
            print('In initialize_richards_vars():')
            print('ALL_SCALARS =', self.ALL_SCALARS)
            print('shape(Ks)   =', np.shape(self.Ks))
            print('shape(Ki)   =', np.shape(self.Ki))
            print('shape(qs)   =', np.shape(self.qs))
            print('shape(qi)   =', np.shape(self.qi))
            print('shape(qr)   =', np.shape(self.qr))
            print('shape(pB)   =', np.shape(self.pB))
            print('shape(pA)   =', np.shape(self.pA))
            print('shape(lam)  =', np.shape(self.lam))
            print('shape(c)    =', np.shape(self.c))
            print('shape(eta)  =', np.shape(self.eta))
            print('shape(qH)   =', np.shape(self.qH))
            print(' ')
            
        #-----------------------------------------------------
        # Compute dz as 1D array from scalars in self.dz_val
        #-----------------------------------------------------
        # NB! Values in self.dz_val are scalars vs. pointers
        # so we can't use the build_layered_var routine.
        #-----------------------------------------------------
        dz_min = self.dz_val.min()
        dz_max = self.dz_val.max()
        if (dz_min == dz_max):
            #----------------------
            # dz is just a scalar
            #----------------------
            self.dz = self.dz_val[0]
        else:
            #-------------------
            # dz is a 1D array
            #-------------------
            self.dz = np.zeros(self.nz, dtype=dtype)

            #--------------------------------------------------
            # Create array of indices.  See build_layered_var
            #--------------------------------------------------
            i = concatenate(([int32(0)], int32(cumsum(self.nz_val))) )
            for j in range(self.n_layers):
                self.dz[ i[j]: i[j+1]-1 ] = self.dz_val[j]

        #----------------------------------------------
        # Compute the z-vector, for plotting profiles
        #----------------------------------------------
        dz = np.repeat(self.dz_val[0], self.nz_val[0])  # (1D ndarray)
        for j in range(1, self.n_layers):
            layer_dz = self.dz_val[j]
            layer_nz = self.nz_val[j]
            dz_j = np.repeat(layer_dz, layer_nz)  # (1D ndarray)
            dz = np.concatenate( (dz, dz_j) )
        ############################################
        # NB! As written (and in IDL version), the
        #     z-vector does not start with 0.
        ############################################
        self.z = np.cumsum(dz)

        #-------------------------------------------------------
        # Note: qi and Ki are created with build_layered_var()
        #-------------------------------------------------------        
        if (self.ALL_SCALARS):
            #----------------------------------
            # Infiltration varies with z only
            #----------------------------------
            self.q = np.zeros(self.nz, dtype=dtype) + self.qi
            self.p = np.zeros(self.nz, dtype=dtype)
            self.K = np.zeros(self.nz, dtype=dtype) + self.Ki
            self.v = np.zeros(self.nz, dtype=dtype)
            #---------------------------------------------------
            self.IN = np.float64(0)   # (infil. rate at surface)
            self.I  = np.float64(0)   # (total infil. depth)
            self.Zw = np.float64(0)   # (wetting front depth)
##            self.I  = float64(1e-6)   # (total infil. depth)
##            self.Zw = float64(1e-6)   # (wetting front depth)

##            if (self.DEBUG):
##                print 'shape(self.v)  =', shape(self.v)
##                print 'shape(P_total) =', shape(self.P_total)
##                print 'type(P_total)  =', type(self.P_total)
                
            #--------------------------------------------
            # Set BC at the surface (done elsewhere ??)
            #--------------------------------------------
            self.v[0] = self.P_total
    
        else:
            #------------------------------------
            # Infiltration varies with x, y & z
            #------------------------------------
            self.q  = zeros((self.nz, self.ny, self.nx), dtype=dtype)
            self.p  = zeros((self.nz, self.ny, self.nx), dtype=dtype) 
            self.K  = zeros((self.nz, self.ny, self.nx), dtype=dtype) 
            self.v  = zeros((self.nz, self.ny, self.nx), dtype=dtype)
            #---------------------------------------------------------------
            self.IN = np.zeros([self.ny, self.nx], dtype=dtype)
            self.I  = np.zeros([self.ny, self.nx], dtype=dtype)
            self.Zw = np.zeros([self.ny, self.nx], dtype=dtype)

            #--------------------------------------
            # Initialize q to qi (qi is 1D or 3D)
            #--------------------------------------
            if (np.size(self.qi) == self.nz):
                for j in range(self.nz):
                    self.q[j,:,:] = self.qi[j]
                # (Can this be done with array operators instead ?)
            else:
                self.q += self.qi

            #--------------------------------------
            # Initialize K to Ki (Ki is 1D or 3D)
            #--------------------------------------
            if (np.size(self.Ki) == self.nz):
                for j in range(self.nz):
                    self.K[j,:,:] = self.Ki[j]
                # (Can this be done with array operators instead ?)
            else:
                self.K += self.Ki

            #--------------------------------------------------
            # If q is now 3D, convert qs to 3D also so we can
            # compute (q - qs) in update_v(). (6/22/10)
            #--------------------------------------------------
            ## if (np.ndim(self.qs) == 1):
            if (np.size(self.qs) == self.nz):
                temp = self.qs.copy()
                self.qs = zeros((self.nz, self.ny, self.nx), dtype=dtype)
                for j in range(self.nz):
                    self.qs[j,:,:] = temp[j]   #######

            #--------------------------------------------------
            # If q is now 3D, convert qH to 3D also so we can
            # compute (q - qH) in ******(). (6/22/10)
            #--------------------------------------------------
            ## if (np.ndim(self.qH) == 1):
            if (np.size(self.qH) == self.nz):
                temp = self.qH.copy()
                self.qH = zeros((self.nz, self.ny, self.nx), dtype=dtype)
                for j in range(self.nz):
                    self.qH[j,:,:] = temp[j]   #######
                    
            #--------------------------------------------
            # Set BC at the surface (done elsewhere ??)
            #--------------------------------------------
            self.v[0,:,:] = self.P_total

        #-------------------------------------------------
        # Print some suggested (i.e. consistent) values
        # for theta_r, theta_i and K_i.  (10/12/10)
        #-------------------------------------------------
        self.print_suggested_values()
        
        ###########################################
        # Override some of the user's settings ??
        ###########################################
##        if (self.SPECIAL_DEFAULTS):
##            self.initialize_theta_r()
##            self.initialize_theta_i()
##            self.initialize_K_i()
            
    #   initialize_richards_vars()
    #-------------------------------------------------------------------
    def initialize_theta_r(self):

        #-------------------------------------------------
        # Note that this is not entirely consistent with
        # the Theta_TBC() function, but that function
        # requires theta_r as an argument.
        #-------------------------------------------------
        # Initialize theta_r to the min allowed value.
        #-------------------------------------------------
        psi_r = self.psi_min_cm
        psi_r = (psi_r / float64(100))   #[cm -> meters]
        
        #--------------------------------------
        # Note:  Both psi's < 0, so ratio > 0
        #--------------------------------------
        self.qr = self.qs * (self.pB / psi_r)**self.lam
    
    #   initialize_theta_r()
    #-------------------------------------------------------------------
    def initialize_theta_i(self):

        #------------------------------------------------
        # Initialize theta_i = qi to the field capacity.
        # Be sure to call initialize_theta_r() first.
        #------------------------------------------------
        self.qi = Theta_TBC( self.psi_field_cm, \
                             self.qs, self.qr, \
                             self.pB, self.pA, \
                             self.c,  self.lam )
        
    #   initialize_theta_i()
    #-------------------------------------------------------------------
    def initialize_K_i(self):

        self.Ki = K_of_Theta_TBC( self.qi, self.Ks, self.qs,
                                  self.qr, self.lam )

    #   initialize_K_i()
    #-------------------------------------------------------------------
    def print_suggested_values(self):

        if (self.DEBUG):
            print('Calling print_suggested_values()...')
            
        #-----------------------------------------------------
        # theta_r is often set to theta_hygroscopic.
        # theta_i is often set to theta_field_capacity.
        #-----------------------------------------------------
        print('=====================================================')
        for k in range(self.n_layers):

##            print 'Ks[k]  =', self.Ks_val[k]
##            print 'Ki[k]  =', self.Ki_val[k]    
##            print 'qs[k]  =', self.qs_val[k]
##            print 'qi[k]  =', self.qi_val[k]
##            print 'qr[k]  =', self.qr_val[k]
##            print 'pB[k]  =', self.pB_val[k]
##            print 'pA[k]  =', self.pA_val[k]
##            print 'lam[k] =', self.lam_val[k]
##            print 'c[k]   =', self.c_val[k]
##            print ' '
##            print 'psi_hygro_cm =', self.psi_hygro_cm
##            print 'psi_field_cm =', self.psi_field_cm
            
            #-------------------------------------------------
            # Compute this by analogy to equations 6-19 and
            # 6-20 in Dingman (2002), using theta_s instead
            # of porosity and recalling lambda = (1/b).
            #-------------------------------------------------
            # Note that this is not entirely consistent with
            # the Theta_TBC() function, but that function
            # requires theta_res as an argument.
            #-------------------------------------------------
##            psi_res   = self.psi_hygro_cm / float64(100)
##            theta_sat = self.qs_val[k]
##            psi_B     = self.pB_val[k]
##            lam       = self.lam_val[k]
##            #--------------------------------------
##            # Note:  Both psi's < 0, so ratio > 0
##            #--------------------------------------
##            theta_res = theta_sat * (psi_B / psi_res)**lam

            #--------------------------------------------
            # If we trust theta_r, then do this instead
            #--------------------------------------------
            theta_res = self.qr_val[k]
            

            theta_hygro = Theta_TBC( self.psi_hygro_cm,
                                     self.qs_val[k],
                                     self.qr_val[k],
                                     self.pB_val[k],
                                     self.pA_val[k],
                                     self.c_val[k],
                                     self.lam_val[k] )
            
            theta_init = Theta_TBC( self.psi_field_cm,
                                    self.qs_val[k],
                                    theta_res,         #######
                                    self.pB_val[k],
                                    self.pA_val[k],
                                    self.c_val[k],
                                    self.lam_val[k] )

            K_init = K_of_Theta_TBC( theta_init,       #######
                                     self.Ks_val[k],
                                     self.qs_val[k],
                                     theta_res,        #######
                                     self.lam_val[k] )

            theta_r = self.qr_val[k]
            theta_i = self.qi_val[k]
            K_i     = self.Ki_val[k]
            print('Suggested initial values for layer', k+1, ':')
            ## print '   theta_r =', theta_res,  'vs.', theta_r
            print('   For theta_r =', theta_r)
            print('   theta_i =', theta_init, '   vs.', theta_i)
            print('   K_i     =', K_init,     'vs.', K_i)
            print('   theta_H =', theta_hygro, '  vs.', theta_r, ' (theta_r)')
            print(' ')
            
        print('=====================================================')                                                
        
    #   print_suggested_values()
    #-------------------------------------------------------------------
    def update(self, time_seconds=None):

        #################################
##        if (self.time_index > 570):
##            self.DEBUG = True
        #################################
            
        #-------------------------------------------------
        # Note: self.IN already set to 0 by initialize()
        #-------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        self.status = 'updating'  # (OpenMI 2.0 convention)
              
        #-------------------------
        # Update computed values 
        #-------------------------
        self.update_surface_influx()  # (= P + SM - ET)

        #------------------------------------
        # Update the Richards eqn variables
        #------------------------------------
        # All layers processed at once
        #------------------------------------
        self.update_theta()
        self.update_psi()
        self.update_K()
        ## self.update_v()
        self.update_v2()
        self.update_Zw()   # (not tested yet ??)

        self.update_infil_rate()
        self.adjust_infil_rate()    # ??????????
        self.update_IN_integral()
        #### self.update_Rg()
        self.update_Rg_integral()
        #### self.update_I()   # (total infiltrated depth)  ############
        self.update_q0()  # (soil moisture at surface)

        #----------------------------------------------
        # Check for NaNs in infiltration (at surface)
        #----------------------------------------------    
        self.check_infiltration()

        #------------------------------------------
        # Read next infil vars from input files ?
        #------------------------------------------
        self.read_input_files()

        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        self.write_output_files()
        ## self.write_output_files(time_seconds)

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        self.update_time()
        self.status = 'updated'  # (OpenMI 2.0 convention)
        if (self.DEBUG):
            print('Completed update:', self.time_index - 1)
            print(' ')

    #   update()
    #-------------------------------------------------------------------
    def update_surface_influx(self):

        if (self.DEBUG):
            print('Calling update_surface_influx()...')
            
        ## P  = self.get_port_data('P',  self.mp, 'METEOROLOGY')
        ## SM = self.get_port_data('SM', self.sp, 'SNOW')
        ## ET = self.get_port_data('ET', self.ep, 'EVAP')

        P  = self.P
        SM = self.SM
        ET = self.ET
        
        ## print 'min(ET), max(ET) =', ET.min(), ET.max()
        
        ### self.P_total = (P + SM)
        self.P_total = (P + SM) - ET

    #   update_surface_influx()
    #-------------------------------------------------------------------
    def update_infil_rate(self):

        #------------------------------------------------------------
        # Notes: IN   = infiltration rate [m/s]
        #        Rg   = groundwater recharge rate [m/s]
        #               (Returned to caller)
        #        Ks   = saturated hydraulic conductivity [m/s]
        #        Ki   = initial hydraulic conductivity [m/s]
        #        qs   = soil moisture content (sat.)  [dimless]
        #        qi   = soil moisture content (init.) [dimless]
        #        qr   = soil residual moisture content [dimless]
        #        pB   = bubbling pressure head [m]
        #        pA   = optional pressure head offset [m]
        #        cid  = cum. infiltration depth (since reset) [m]
        #         P   = precipitation rate [m/s]
        #        SM   = snowmelt rate [m/s]

        #        Note that the infiltration rate has a max possible
        #        value of (P + SM - ET) and asymptotes to Ks as the
        #        total infiltrated depth increases.

        #        Total infiltrated depth is incremented in the
        #        calling function, called Infiltration.
        #------------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_infil_rate()...')
            
        #-------------------------------------------
        # Infiltration rate is flow rate at surface
        #---------------------------------------------
        # Infiltration rate for node just below the
        # surface will & should be different than at
        # the surface and won't compare as well with
        # Green-Ampt, etc.
        #---------------------------------------------
        if (self.SINGLE_PROFILE):    
            self.IN = self.v[0]
            ## self.IN = self.v[1]
        else:
            self.IN = self.v[0,:,:]
            ## self.IN = self.v[1,:,:]
            

##        print 'SINGLE_PROFILE   =', self.SINGLE_PROFILE
##        print 'min(IN), max(IN) =', self.IN.min(), self.IN.max()
        
        #-----------------------------------------------------------
        #Richards' equation is only used in the so-called "upper
        #layers".  There can be between 1 and 3 of these layers.
        #To avoid high computational cost in the less dynamic
        #"lower zone" below, a simplified method is used to route
        #flow through the lower zone to the water table and to
        #estimate a "groundwater recharge rate", Rg.  This is done
        #using the vertical flow rate at the bottom of the set of
        #upper layers and perhaps other information.  If the water
        #table intersects the upper layers, then Rg is computed as
        #the vertical flow rate in the layer just above the water
        #table.
        #-----------------------------------------------------------
        #The simplest method of estimating Rg is to simply set it
        #to the vertical flow rate at the bottom of the upper zone.
        #However, travel time through the lower zone is not taken
        #into account.  If Rg < 0, then water can even be drawn
        #upward from the water table.
        #-----------------------------------------------------------
        #Another simple idea for estimating Rg is to assume that
        #psi varies linearly in the lower zone:   psi = a*z + b.
        #Since psi=0 at the bottom and psi=psi_zb at the top, we
        #can compute a and b as:
        #    a = -psi_zb/(z_wt - z_zb)
        #    b = -a * z_wt
        #The flow rate in the lower zone is then computed as:
        #    v = K(a*z + b) * (1 - a)
        #It follows that the flow rate at the bottom of the lower
        #zone can be written as:
        #    Rg = v(z_wt) = v(z_zb) * (Ks / K(psi_zb)).
        #Since K(psi_zb) <= Ks, we have v(z_zb) < Rg < Ks.
        #If psi_zb < (z_zb - z_wt), then we get Rg < 0 and water
        #can be drawn upward from the water table.
        #-----------------------------------------------------------
        
        #----------------------------------------
        # For testing:  Plot the theta profiles
        #----------------------------------------
        PLOT = False
        if (PLOT):    
            matplotlib.pyplot.figure(1)
            #** wait, 0.005
            
            ymin = self.qi.min()
            ymax = (self.qs + float64(0.05)).max()
            if (self.SINGLE_PROFILE):    
                matplotlib.pyplot.plot(self.z, self.q, marker='+')
                matplotlib.pyplot.xlabel('Depth [meters]')
                matplotlib.pyplot.ylim(array(ymin, ymax))
                matplotlib.pyplot.axis('image')
                matplotlib.pyplot.ylabel('Soil moisture')
                matplotlib.pyplot.show()
            else:    
                matplotlib.pyplot.plot(self.z, self.q[:,2,2], marker='+')
                matplotlib.pyplot.xlabel('Depth [meters]')
                matplotlib.pyplot.ylim(array(ymin, ymax))
                matplotlib.pyplot.axis('image')
                matplotlib.pyplot.ylabel('Soil moisture')
                matplotlib.pyplot.show()
        
        #------------------------------------
        #For testing:  Plot the psi profiles
        #------------------------------------
        if (PLOT):    
            matplotlib.pyplot.figure(2)
            #** wait, 0.005
            yrange = array([-float32(3.0), float32(0.5)])    #(Log case)
            ytitle = '-Log(-Pressure head) [m]'
            #--------------------------------------
            #ytitle = 'Pressure head [m]'
            #yrange = [-20.0, 0.0]  ;(Linear case)
            #--------------------------------------
            if (self.SINGLE_PROFILE):    
                y = -float64(1) * log(absolute(self.p) + 1)
                matplotlib.pyplot.plot(self.z, y, marker='+')
                matplotlib.pyplot.xlabel('Depth [meters]')
                matplotlib.pyplot.ylim(yrange)
                matplotlib.pyplot.axis('image')
                matplotlib.pyplot.ylabel(ytitle)
                matplotlib.pyplot.show()
            else:    
                y = -float64(1) * log(absolute((self.p)[:,2,2]) + float64(1))
                matplotlib.pyplot.plot(self.z, y, marker='+')
                matplotlib.pyplot.xlabel('Depth [meters]')
                matplotlib.pyplot.ylim(yrange)
                matplotlib.pyplot.axis('image')
                matplotlib.pyplot.ylabel(ytitle)
                matplotlib.pyplot.show()
                  
    #   update_infil_rate() 
    #-------------------------------------------------------------------
    def update_Rg(self):
        
        #-----------------------------------------------------
        # Notes:  Override infil_base's method by same name.
        #-----------------------------------------------------
        #  Already updated by update_infil_rate(), but need
        #  this here so it doesn't get overwritten.
        #-----------------------------------------------------
        if (self.DEBUG):
            print('Calling update_Rg()...')
            
        pass
    
    #   update_Rg()
    #-------------------------------------------------------------------
    def update_q0(self):

        if (self.DEBUG):
            print('Callling update_q0()...')
            
        if (self.ALL_SCALARS): 
            self.q0 = self.q[0]
        else:    
            self.q0 = self.q[0,:,:]
            
    #   update_q0()
    #-----------------------------------------------------------------------
    def update_theta(self, REPORT=False):

        #----------------------------------------------------------
        # Notes:  This procedure updates the soil moisture, theta
        #         as a function of the vertical flow rate, v.

        #         q, v, qs and qr are pointers.
        #         dz, nz and dt are scalars.

        #         Theta is called "q" here .
        #----------------------------------------------------------
        # 3/6/08  Bug fix. Was OK for 1D but sign error for 3D.
        #        Introduced Z_Derivative_1D function to clarify
        #        changed sign in dtheta expression and left
        #        Z_Derivative_3D function unchanged.
        #----------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_theta()...')
        
        if (self.SINGLE_PROFILE):    
            #--------------------------------
            # Theta is a 1D array & so is v
            #--------------------------------
            # dz may be scalar or 1D array
            #--------------------------------
            dv_dz = Z_Derivative_1D(self.v, self.dz)
            
            #-------------------------------------------------
            # What should we do at lower boundary ?
            # dv_dz = 0 matches v=const in update_v().
            #-------------------------------------------------
            dv_dz[self.nz - 1] = np.float64(0)
            ## dv_dz[self.nz - 1] = dv_dz[self.nz - 2]
        else:    
            #--------------------------------
            # Theta is a 3D array & so is v
            #--------------------------------
            dv_dz = Z_Derivative_3D(self.v, self.dz)
            
            #-------------------------------------------------
            # What should we do at lower boundary ?
            # dv_dz = 0 matches v=const in Update_Richards_V
            #-------------------------------------------------
            dv_dz[self.nz - 1,:,:] = np.float64(0)
            ## dv_dz[self.nz-1, :, :] = dv_dz[self.nz-2, :, :]
        
        #----------------------------------------------
        # Update soil moisture, theta  (mass balance)
        #----------------------------------------------
        self.q = self.q - (dv_dz * self.dt)

        #-----------------------------
        # Make sure theta <= theta_s
        # and that  theta >= theta_H
        #--------------------------------------------------
        # NB! We don't need this when we check for layers
        # that are filling or losing in update_v().
        #---------------------------------------------------
        # (10/9/10) Next 2 lines lead to error when the
        # time_index reaches 572 for test_plane_csm/plane1
        #---------------------------------------------------
##        self.q = np.minimum( self.q, self.qs )
##        self.q = np.maximum( self.q, self.qH )

        if (self.DEBUG):
        ## if (True):
            print('min(q), max(q) =', self.q.min(), self.q.max())
            ## self.check_theta()
            
        #------------------
        # Optional report
        #------------------
        #if (REPORT):
        #    print 'dv_dz =', dv_dz[0:3]
        #    print 'theta =', self.q[0:3]
        #    # print ' '

    #   update_theta()
    #-----------------------------------------------------------------------
    def check_theta(self):

        w = np.where( logical_or( (self.q < self.qH),
                                     (self.q > self.qs)) )

        if (w[0].size > 0):
            print('############################################')
            print('ERROR: Theta not in [theta_H, theta_s].')
            print('       Aborting model run.')
            print('############################################')
            self.DONE = True
        
    #   check_theta()
    #-----------------------------------------------------------------------
    def update_psi(self, REPORT=False):

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
        #print,'size(self.q)    = ', size(self.q)
        #print,'size(self.p)    = ', size(self.p)
        #print,'size(self.qs)   = ', size(self.qs)
        #print,'size(self.qr)   = ', size(self.qr)
        #print,'size(self.pB)   = ', size(self.pB)
        #print,'size(self.pA)   = ', size(self.pA)
        #print,'size(self.lam)  = ', size(self.lam)
        #print,'size(self.c)    = ', size(self.c)
        #print,' '
        
        #---------------------------------------
        # Compute the "effective saturation"
        # Relative saturation = theta/porosity
        #---------------------------------------
        if (self.SINGLE_PROFILE):    
            #--------------------------------
            # All of the vars are 1D arrays
            #--------------------------------
            S_eff = (self.q - self.qr) / (self.qs - self.qr)
            _pow = -float64(1) * (self.c / self.lam)
            arg = (S_eff ** _pow - 1.0) ** (1.0 / (self.c))
            self.p = (self.pB * arg) - self.pA
        else:    
            #--------------------------------------
            # Each var is either a 1D or 3D array
            #--------------------------------------
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
                #------------------------------------------------
                # NB!  It is OK to raise a grid to a grid power
                #------------------------------------------------
                S_eff = (self.q[j,:,:] - qr) / (qs - qr)     #(grid)
                _pow = -float64(1) * (c / lam)                       #(grid or scalar)
                arg = (S_eff ** _pow - float64(1)) ** (float64(1) / c)             #(grid)
                self.p[j,:,:] = (pB * arg) - pA              #(grid)

        if (self.DEBUG):
            print('min(p), max(p) =', self.p.min(), self.p.max())
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('S_eff = ', S_eff[0:4])
            print('psi   = ', self.p[0:3])
            #print,' '

    #   update_psi()
    #-----------------------------------------------------------------------
    def update_K(self, REPORT=False):

        #------------------------------------------------------------
        # Notes: This procedure updates the hydraulic conductivity,
        #        K, as a function of the pressure head, psi, via
        #        the "Brooks-Corey" (B-C) or "transitional Brooks-
        #        Corey" (TB-C) relation.

        #        lambda = pore size distribution parameter
        #        eta    = "pore-disconnectedness" parameter
        #        eta    = 2d + (3d * lambda)

        #        There is a typo in R.E. Smith's AGU monograph in
        #        equation (2.14), where eta should be -eta.

        #        See "Infiltration Theory for Hydrologic Applica-
        #        tions" by R.E. Smith (2002), p. 21-22.

        #        For standard Brooks-Corey we would have:
        #            pow = -1d * (*eta)
        #            K_r = (*p / *pB)^pow
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
            #--------------------------------
            _pow = -1.0 * (self.eta / self.c)
            Kr = (1.0 + ((self.p + self.pA) / self.pB) ** self.c) ** _pow
            Kr = maximum((minimum(Kr, 1.0)), 0.0)
            self.K = self.Ks * Kr
        else:    
            #--------------------------------------
            # Each var is either a 1D or 3D array
            #--------------------------------------
            dim_Ks  = np.ndim(self.Ks)
            dim_pB  = np.ndim(self.pB)
            dim_pA  = np.ndim(self.pA)
            dim_eta = np.ndim(self.eta)
            dim_c   = np.ndim(self.c)

            for j in range(self.nz):
                #--------------------------------------------------
                # At a given z, every input var is scalar or grid
                #--------------------------------------------------
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
                #------------------------------------------------
                # NB!  It is OK to raise a grid to a grid power
                #------------------------------------------------
                arg = (self.p[j,:,:] + pA) / pB         #(grid)
                _pow = -1.0 * (eta / c)          #(grid or scalar)
                Kr = (1.0 + arg ** c) ** _pow    #(grid)
                Kr = maximum((minimum(Kr, 1.0)), 0.0)   #(grid)
                self.K[j,:,:] = (Ks * Kr)               #(grid)

        if (self.DEBUG):
            print('min(K), max(K) =', self.K.min(), self.K.max())
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('K = ', self.K[0:3])
            #print,' '

    #   update_K()
    #-----------------------------------------------------------------------
    def update_v(self, REPORT=False):

        #-----------------------------------------------------------
        # Notes: This procedure updates the vertical flow rate
        #        at each level, v, as a function of psi, K & theta.
        #        q, v, psi, K, qs and qr are pointers.
        #        dz, nz, and dt are scalars.

        #        If R is a grid or grid sequence, then the infil
        #        vars are initialized as 3D vs. 1D by the routine
        #        initialize_layer_vars(), and no action here is
        #        required.

        #        P_total = (P + SM - ET) for current timestep
        #              (could be either grid or scalar)
        #              and is computed & passed by caller.
        #        vB  = flow rate [m/s] at bottom of a cell
        #        vT  = flow rate [m/s] at top of a cell
        #              (vT[k] = vB[k-1], except for k=0, see fig.)
        #        dz  = z-distance between nodes
        #        nn  = number of nodes on z-axis

        ###########################################################
        #        K_bar is a "mean value" of K.
        #        Using K_bar = K doesn't work for the case of
        #        redistribution due to evaporation, but seems
        #        to work OK for many other cases.
        ###########################################################
        #
        #        K, psi and theta are assumed to be uniform
        #        within any given soil layer, while the flow
        #        rates are for the boundaries between layers.

        #        If one is a scalar or grid, they all are.

        #        The first derivative of psi is computed using
        #        psi values on either side of a boundary and
        #        the z-distance between the layer centers.
        #-----------------------------------------------------------
        #
        #               surface
        #               =======  vT[0]
        #
        #        vB[0]  =======  vT[1]
        #
        #        vB[1]  =======  vT[2]
        #        
        #-----------------------------------------------------------        

        #----------------
        # For debugging
        #----------------
        if (self.DEBUG):
            print('Calling update_v()...')
            print('   SINGLE_PROFILE =', self.SINGLE_PROFILE)
##        print 'type(self.p)   =', type(self.p)
##        print 'type(self.K)   =', type(self.K)
##        print 'nz             =', self.nz
##        print 'size(self.dz)  =', np.size(self.dz)
##        print ' '
        
        #-----------------------------------------------
        # NB! There are better ways to compute K_bar
        #     than the mean value used here.  Another
        #     method is discussed by R.E. Smith, at
        #     the top of page 192 in his AGU monograph.
        #-----------------------------------------------
        if (self.SINGLE_PROFILE):    
            #--------------------------------
            # All of the vars are 1D arrays
            #--------------------------------
            # dp_dz = (p_below - p) / dz
            #--------------------------------            
            dp_dz   = Z_Derivative_1D(self.p, self.dz)
            K_below = np.roll(self.K, -1, axis=0)
            K_bar   = (self.K + K_below) / 2.0
            vB = K_bar * (1.0 - dp_dz)      # (bottom of cell)
            vT = np.roll(vB, 1, axis=0)  # (top of cell)
            
            #--------------------------------------
            # If R is not a scalar, then the next
            # line will generate an error
            #--------------------------------------------------
            # Initialize_Infil_Vars now sets self.all_scalars
            # (used to define SINGLE_PROFILE) based on
            # infil, precip, snowmelt and ET. (3/19/07)
            #--------------------------------------------------
##            if (self.q[0] < self.qs[0]):
##                vT[0] = self.P_total    # (not saturated yet)
##            else:
##                #-----------------------------------------------
##                # Recall that d/dt(theta) = -d/dz(v).  So once
##                # surface has reached saturation, we must have
##                # the same v-value in the top 2 cells.
##                #-----------------------------------------------
##                print '### REACHED SATURATION at time_index =', self.time_index
##                vT[0] = vB[0]
##                print 'theta[0], vT[0] =', self.q[0], vT[0]
##                print 'p[0], p[1]  =', self.p[0], self.p[1]
##                print 'K[0], Ks[0] =', self.K[0], self.Ks[0]
##                print ' '
                
            vT[0] = self.P_total   # (before 11/9/10)

##            print 'theta[0], vT[0] =', self.q[0], vT[0]
##            print 'p[0], p[1]  =', self.p[0], self.p[1]
##            print 'K[0], Ks[0] =', self.K[0], self.Ks[0]
##            print ' '
            
##            print 'psi[0]   =', self.p[0]
##            print 'dp_dz[0] =', dp_dz[0]
##            print 'K[0]     =', self.K[0]
##            print 'K[1]     =', self.K[1]
##            print 'K_bar[0] =', K_bar[0]
##            print 'vB[0]    =', vB[0]
##            print 'vT[0]    =', vT[0]
##            print '----------------------------------'
                
            #-------------------------------------------------
            # "Gravity drainage" bottom BC (dp_dz=0 => v=Ki)   ;************
            # Assumes wetting front doesn't reach bottom ?
            #-------------------------------------------------
            #** vB[self.nz - 1] = *Ki
            ## vB[self.nz - 1] = self.Ki[self.nz - 1]  # (11/2/10)
            
            #-------------------------------
            # Simple bottom BC (dv_dz = 0)                     ;************
            #-------------------------------
            vB[self.nz - 1] = vB[self.nz - 2]
        else:
            #-------------------------------------
            # Compute pressure gradient, dpsi_dz
            #-------------------------------------
            # print 'Computing dp_dz...'
            dp_dz = Z_Derivative_3D(self.p, self.dz)     #(pointer args)
            
            #----------------------------------------
            # Theta, K, psi and v are all 3D arrays
            #----------------------------------------
            # print 'Computing K_bar...'
            K_below = np.roll(self.K, -1, axis=0)
            K_bar   = (self.K + K_below) / 2.0
            # print 'Computing vB...'
            vB = K_bar * (1.0 - dp_dz)
            vT = np.roll(vB, 1, axis=0)
            # print 'Setting vT[0,:,:]...'
            vT[0,:,:] = self.P_total        #(works if R is scalar or grid)
            
            #-------------------------------------------------
            # "Gravity drainage" bottom BC (dp_dz=0 => v=Ki)   ;************
            # Assumes wetting front doesn't reach bottom ?
            #-------------------------------------------------
            #** vB[*,*,self.nz-1L] = *self.Ki
            
            #-------------------------------
            # Simple bottom BC (dv/dz = 0)                     ;************
            #-------------------------------
            # print 'Setting simple bottom BC...'
            vB[self.nz - 1,:,:] = vB[self.nz - 2,:,:]
        
        #----------------------------------------------------
        # Flow rate into any layer must be less than the
        # "amount of space available", while flow rate
        # out must be less than "amount of water available"
        #-------------------------------------------------------
        # With this, it seems to be unnecessary to force theta
        # to stay within limits in update_theta() function.
        #-------------------------------------------------------
        filling   = where((vT - vB) >= 0)
        n_filling = filling[0].size
        losing    = where((vT - vB)   < 0)
        n_losing  = losing[0].size

        #----------------
        # For debugging
        #----------------
##        print 'At n_filling test...'
##        print 'shape(qs) =', np.shape(self.qs)
##        print 'shape(q)  =', np.shape(self.q)
##        print 'shape(vT) =', np.shape(vT)
##        print 'shape(vB) =', np.shape(vB)
##        print 'shape(v)  =', np.shape(self.v)

        #-----------------------------------------------------
        # Note that qs may need to be converted from 1D to
        # 3D in initialize_richards_vars() for this to work.
        #-----------------------------------------------------
        if (n_filling != 0):
            space_avail = (self.dz / self.dt) * (self.qs - self.q)
            # print 'shape(space_avail) =', np.shape(space_avail)
            vT[filling] = (minimum(vT[filling], (space_avail[filling] + vB[filling])))

        # print 'At n_losing test...'
        if (n_losing != 0):    
            #-----------------------------------------------------------
            # Note that for both B-C and TB-C, psi goes to -Infinity
            # as theta goes to theta_r (S_eff goes to zero).  However,
            # natural soils do not have heads (tensions) less than
            # -31,000 cm.  In this range they absorb water from the air
            # (hygroscopic).  The function Theta_Min finds the smallest
            # theta-value, theta_H, corresponding to this limiting
            # psi-value, and uses it to determine water available for
            # exfiltration. It is incorrect to use theta_i or theta_r.
            # During evaporation, the surface flow rate equals the
            # evaporation rate (which is < 0) until the soil moisture
            # at the surface drops to theta_H.  It then gradually
            # approaches a value of 0 (from values < 0) because though
            # water is drawn from below it cannot keep up with the
            # demand at the surface.  This can be seen in the plot of
            # infiltration rate when the EVAP keyword is set, in
            # richards0.pro.
            #-----------------------------------------------------------
            water_avail = (self.dz / self.dt) * (self.q - self.qH)
            vT[losing] = (maximum(vT[losing], (vB[losing] - water_avail[losing])))
        
        #------------------------
        # Update heap var for v
        #------------------------
        # print 'Setting self.v to vT...'
        self.v = vT

        ## if (self.q[0] >= self.qs[0]):
        if (self.time_index > 580):
            print('v[0], v[1] =', self.v[0], ', ', self.v[1])  ################

            
        #-----------------------------------
        # Return flow rate in bottom layer
        #-----------------------------------
        if (self.SINGLE_PROFILE):    
            self.Rg = self.v[self.nz - 1]
        else:    
            self.Rg = self.v[self.nz - 1,:,:]

        if (self.DEBUG):
            print('min(v), max(v) =', self.v.min(), self.v.max())
            
        #------------------
        # Optional report
        #------------------
        #if (REPORT) then begin
        #    print,'dpsi_dz = ', dp_dz[0:3]
        #    print,'vT      = ', vT[0:3]
        #    print,'vB      = ', vB[0:3]
        #endif
        
    #   update_v()
    #-----------------------------------------------------------------------
    def update_v2(self, REPORT=False):

        #-----------------------------------------------------------
        # Notes: This procedure updates the vertical flow rate
        #        at each level, v, as a function of psi, K & theta.
        #        q, v, psi, K, qs and qr are pointers.
        #        dz, nz, and dt are scalars.

        #        If R is a grid or grid sequence, then the infil
        #        vars are initialized as 3D vs. 1D by the routine
        #        initialize_layer_vars(), and no action here is
        #        required.

        #        P_total = (P + SM - ET) for current timestep
        #              (could be either grid or scalar)
        #              and is computed & passed by caller.
        #        v   = flow rate [m/s] at bottom of a cell
        #        dz  = z-distance between nodes
        #        nn  = number of nodes on z-axis

        ###########################################################
        #        K_bar is a "mean value" of K.
        #        Using K_bar = K doesn't work for the case of
        #        redistribution due to evaporation, but seems
        #        to work OK for many other cases.
        ###########################################################
        #
        #        K, psi and theta are assumed to be uniform
        #        within any given soil layer, while the flow
        #        rates are for the boundaries between layers.

        #        If one is a scalar or grid, they all are.

        #        The first derivative of psi is computed using
        #        psi values on either side of a boundary and
        #        the z-distance between the layer centers.
        #-----------------------------------------------------------        

        #----------------
        # For debugging
        #----------------
        if (self.DEBUG):
            print('Calling update_v()...')
            print('   SINGLE_PROFILE =', self.SINGLE_PROFILE)
##        print 'type(self.p)   =', type(self.p)
##        print 'type(self.K)   =', type(self.K)
##        print 'nz             =', self.nz
##        print 'size(self.dz)  =', np.size(self.dz)
##        print ' '
        
        #-----------------------------------------------
        # NB! There are better ways to compute K_bar
        #     than the mean value used here.  Another
        #     method is discussed by R.E. Smith, at
        #     the top of page 192 in his AGU monograph.
        #-----------------------------------------------
        if (self.SINGLE_PROFILE):    
            #--------------------------------
            # All of the vars are 1D arrays
            #--------------------------------
            # dp_dz = (p_below - p) / dz
            #--------------------------------            
            dp_dz   = Z_Derivative_1D(self.p, self.dz)
            K_below = np.roll(self.K, -1, axis=0)
            K_bar   = (self.K + K_below) / 2.0
            v       = K_bar * (1.0 - dp_dz)      # (bottom of cell)
            
            #--------------------------------------
            # If R is not a scalar, then the next
            # line will generate an error
            #--------------------------------------------------
            # Initialize_Infil_Vars now sets self.all_scalars
            # (used to define SINGLE_PROFILE) based on
            # infil, precip, snowmelt and ET. (3/19/07)
            #--------------------------------------------------
            if (self.q[0] < self.qs[0]):
                v[0] = self.P_total    # (not saturated yet)
##            else:
##                #-----------------------------------------------
##                # Recall that d/dt(theta) = -d/dz(v).  So once
##                # surface has reached saturation, we must have
##                # the same v-value in the top 2 cells.
##                #-----------------------------------------------
##                print '### REACHED SATURATION at time_index =', self.time_index
##                v[0] = v[1]
                
##                print 'theta[0], vT[0] =', self.q[0], vT[0]
##                print 'p[0], p[1]  =', self.p[0], self.p[1]
##                print 'K[0], Ks[0] =', self.K[0], self.Ks[0]
##                print ' '
                
            ## v[0] = self.P_total   # (before 11/9/10)

##            print 'theta[0], v[0] =', self.q[0], v[0]
##            print 'p[0], p[1]  =', self.p[0], self.p[1]
##            print 'K[0], Ks[0] =', self.K[0], self.Ks[0]
##            print ' '
            
##            print 'psi[0]   =', self.p[0]
##            print 'dp_dz[0] =', dp_dz[0]
##            print 'K[0]     =', self.K[0]
##            print 'K[1]     =', self.K[1]
##            print 'K_bar[0] =', K_bar[0]
##            print 'v[0]      =', v[0]
##            print '----------------------------------'
                
            #-------------------------------------------------
            # "Gravity drainage" bottom BC (dp_dz=0 => v=Ki)   ;************
            # Assumes wetting front doesn't reach bottom ?
            #-------------------------------------------------
            ## v[self.nz - 1] = self.Ki[self.nz - 1]  # (11/2/10)
            
            #-------------------------------
            # Simple bottom BC (dv_dz = 0)                     ;************
            #-------------------------------
            v[self.nz - 1] = v[self.nz - 2]
        else:
            #-------------------------------------
            # Compute pressure gradient, dpsi_dz
            #-------------------------------------
            # print 'Computing dp_dz...'
            dp_dz = Z_Derivative_3D(self.p, self.dz)     #(pointer args)
            
            #----------------------------------------
            # Theta, K, psi and v are all 3D arrays
            #----------------------------------------
            # print 'Computing K_bar...'
            K_below = np.roll(self.K, -1, axis=0)
            K_bar   = (self.K + K_below) / 2.0
            # print 'Computing vB...'
            v = K_bar * (1.0 - dp_dz)
            # print 'Setting vT[0,:,:]...'
            v[0,:,:] = self.P_total        #(works if R is scalar or grid)
            
            #-------------------------------------------------
            # "Gravity drainage" bottom BC (dp_dz=0 => v=Ki)   ;************
            # Assumes wetting front doesn't reach bottom ?
            #-------------------------------------------------
            #** v[*,*,self.nz-1L] = *self.Ki
            
            #-------------------------------
            # Simple bottom BC (dv/dz = 0)                     ;************
            #-------------------------------
            # print 'Setting simple bottom BC...'
            v[self.nz - 1,:,:] = v[self.nz - 2,:,:]

        #---------------------
        # Save the flow rates
        #---------------------
        self.v = v
        
        #----------------
        # For debugging
        #----------------
##        print 'At n_filling test...'
##        print 'shape(qs) =', np.shape(self.qs)
##        print 'shape(q)  =', np.shape(self.q)
##        print 'shape(vT) =', np.shape(vT)
##        print 'shape(vB) =', np.shape(vB)
##        print 'shape(v)  =', np.shape(self.v)

        ## if (self.q[0] >= self.qs[0]):
        ## if (self.time_index > 580) and (self.time_index < 700):
        print('v[0], v[1] =', self.v[0], ', ', self.v[1])  ################
 
        #-----------------------------------
        # Return flow rate in bottom layer
        #-----------------------------------
        if (self.SINGLE_PROFILE):    
            self.Rg = self.v[self.nz - 1]
        else:    
            self.Rg = self.v[self.nz - 1,:,:]

        if (self.DEBUG):
            print('min(v), max(v) =', self.v.min(), self.v.max())
            
        #------------------
        # Optional report
        #------------------
        #if (REPORT) then begin
        #    print,'dpsi_dz = ', dp_dz[0:3]
        #endif
        
    #   update_v2()
    #-----------------------------------------------------------------------
    def update_Zw(self, REPORT=False):

        #------------------------------------------------------------
        # Note: This procedure attempts to identify the depth of the
        #       wetting front from examination of the theta values.

        #       Notice that it is not assumed that the soil moisture
        #       profile is a decreasing function from surface down.
        #       If soil moisture profile starts to increase as we
        #       approach the water table, this should still work.

        #       Notice also that limiting theta value may not be
        #       equal to theta_s.  For example, it approaches a
        #       smaller value for a sustained (R lt K_s).
        #------------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_Zw()...')
            
        if (self.SINGLE_PROFILE):    
            q_below = np.roll(self.q, -1, axis=0)
            diff = (self.q - q_below)
            diff[self.nz - 1] = float64(0)
            indices = where(diff > 0)        # (must be > not >=)
            nd = indices[0].size    
            
            if (nd == 0):    
                #----------------------------------
                # Wetting front is at the surface
                #----------------------------------------------
                # This can happen for an equilibrium profile
                # that is monotonically increasing with depth
                # or for case where theta = theta_i for all z
                #----------------------------------------------
                self.Zw = float64(0)
            else:    
                imax = indices[0][nd - 1]  ########################
                
                #----------------------------------
                # This is one way to define Z_wet
                #----------------------------------
                #;*self.Zw = imax * (*self.dz)
                
                #----------------------------------------------
                # Get min and max theta of decreasing section
                #----------------------------------------------
                frac = float64(0.2)
                tmax = nanmax(self.q[indices])
                tmin = nanmin(self.q[indices])
                tmid = tmin + frac * (tmax - tmin)
                w    = where(self.q[0: (imax+1)] > tmid)
                nw2  = size(w[0])
                if (nw2 > 0):    
                    imax2 = w[0][nw2 - 1]   ###################
                    self.Zw = imax2 * self.dz
                else:    
                    self.Zw = float32(0)         #*******  IS THIS RIGHT ? *******
        else:    
            
            q_below = np.roll(self.q, -1, axis=0)
            diff = (self.q - q_below)
            diff[self.nz - int32(1),:,:] = float64(0)
            n_dz = size(self.dz)     #(either 1 or self.nz)
            
            #-----------------------------
            #Zero out the 2D array *self.Zw
            #---------------------------------------
            #Note: *self.Zw should have already been
            #set to be a 2D array of correct size
            #---------------------------------------
            self.Zw = minimum((maximum(self.Zw, float64(0))), float64(0))
            
            #-------------------------------------
            #Loop over the z-levels to find Z_wet
            #-------------------------------------
            for j in range(self.nz - 1):           #(nz-1) vs. nz
                diff_j = diff[j,:,:]
                next_diff_j = diff[j + 1,:,:]     #*******
                if (n_dz == 1):    
                    dz = self.dz
                else:    
                    dz = self.dz[j]
                
                #------------------------------------------------------
                #NB!  We are looking for a local min in theta profile.
                #------------------------------------------------------
                #NB!  If theta is same at all levels, then we will
                #never get (diff_j GT 0) so Z_wet will remain at 0,
                #even if (theta eq theta_s) at all levels !!!!    ****************
                #How can we get Z_wet = Z_bot = (nz-1) * dz ???  ****************
                #------------------------------------------------------
                #NB!  If theta increases at all levels, which can
                #happen for an equilibrium profile with water table,
                #then we will never get (diff_j GT 0.0) and Z_wet will
                #remain at 0.
                #------------------------------------------------------
                #NB!  For (j eq (nz-2)), we have (next_diff_j EQ 0.0).
                #------------------------------------------------------
                IDs = where(logical_and((diff_j > 0), \
                            (next_diff_j <= 0)))
                n_IDs = size(IDs[0])
                
                if (n_IDs != 0):    
                    self.Zw[IDs] = (j * dz)

    #   update_Zw()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        #-----------------------------------------------------
        # Notes:  Override infil_base's method by same name.
        #-----------------------------------------------------
        self.Ks_unit  = []  # (empty lists to hold file objects)
        self.Ki_unit  = []
        self.qs_unit  = []
        self.qi_unit  = []
        self.qr_unit  = []
        self.pB_unit  = []
        self.pA_unit  = []
        self.lam_unit = []
        self.c_unit   = []
        
        for j in range(self.n_layers):
            self.Ks_unit.append(  model_input.open_file(self.Ks_type[j],  self.Ks_file[j]) )
            self.Ki_unit.append(  model_input.open_file(self.Ki_type[j],  self.Ki_file[j]) )
            self.qs_unit.append(  model_input.open_file(self.qs_type[j],  self.qs_file[j]) )
            self.qi_unit.append(  model_input.open_file(self.qi_type[j],  self.qi_file[j]) )
            self.qr_unit.append(  model_input.open_file(self.qr_type[j],  self.qr_file[j]) )
            self.pB_unit.append(  model_input.open_file(self.pB_type[j],  self.pB_file[j]) )
            self.pA_unit.append(  model_input.open_file(self.pA_type[j],  self.pA_file[j]) )
            self.lam_unit.append( model_input.open_file(self.lam_type[j], self.lam_file[j]) )
            self.c_unit.append(   model_input.open_file(self.c_type[j],   self.c_file[j]) )

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        #-----------------------------------------------------
        # Notes:  Override infil_base's method by same name.
        #-----------------------------------------------------
        rti = self.rti

        for j in range(self.n_layers):        
            Ks_val = model_input.read_next(self.Ks_unit[j], self.Ks_type[j], rti)
            if (Ks_val is not None): self.Ks_val[j] = Ks_val

            Ki_val = model_input.read_next(self.Ki_unit[j], self.Ki_type[j], rti)
            if (Ki_val is not None): self.Ki_val[j]  = Ki_val

            qs_val = model_input.read_next(self.qs_unit[j], self.qs_type[j], rti)
            if (qs_val is not None): self.qs_val[j]  = qs_val

            qi_val = model_input.read_next(self.qi_unit[j], self.qi_type[j], rti)
            if (qi_val is not None): self.qi_val[j]  = qi_val
            
            qr_val = model_input.read_next(self.qr_unit[j], self.qr_type[j], rti)
            if (qr_val is not None): self.qr_val[j]  = qr_val

            pB_val = model_input.read_next(self.pB_unit[j], self.pB_type[j], rti)
            if (pB_val is not None): self.pB_val[j]  = pB_val

            pA_val = model_input.read_next(self.pA_unit[j], self.pA_type[j], rti)
            if (pA_val is not None): self.pA_val[j]  = pA_val

            lam_val = model_input.read_next(self.lam_unit[j], self.lam_type[j], rti)
            if (lam_val is not None): self.lam_val[j]  = lam_val

            c_val = model_input.read_next(self.c_unit[j], self.c_type[j], rti)
            if (c_val is not None): self.c_val[j]  = c_val

            #---------------------------------------------------------
            # If we read a lambda value from a file, then we need to
            # compute and save corresponding eta = [2 + (3*lambda)]
            #---------------------------------------------------------
            #### if not(self.lam_unit[j].closed):  ############
            if (self.lam_type[j] == 1) or (self.lam_type[j] == 3):
                self.eta_val[j] = float64(2) + (float64(3) * self.lam_val[j])
                                 
            #-----------------------------------------
            # Update qH, given by Theta_TBC function
            #-----------------------------------------
            self.qH_val[j] = Theta_TBC( self.psi_hygro_cm, \
                                        self.qs_val[j], self.qr_val[j], \
                                        self.pB_val[j], self.pA_val[j], \
                                        self.c_val[j],  self.lam_val[j] )

    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):

        #-----------------------------------------------------
        # Notes:  Override infil_base's method by same name.
        #-----------------------------------------------------
        for j in range(self.n_layers):
            if (self.Ks_type[j]  != 'Scalar'): self.Ks_unit[j].close()        
            if (self.Ki_type[j]  != 'Scalar'): self.Ki_unit[j].close()
            if (self.qs_type[j]  != 'Scalar'): self.qs_unit[j].close()
            if (self.qi_type[j]  != 'Scalar'): self.qi_unit[j].close()
            if (self.qr_type[j]  != 'Scalar'): self.qr_unit[j].close()
            if (self.pB_type[j]  != 'Scalar'): self.pB_unit[j].close()
            if (self.pA_type[j]  != 'Scalar'): self.pA_unit[j].close()        
            if (self.lam_type[j] != 'Scalar'): self.lam_unit[j].close()
            if (self.c_type[j]   != 'Scalar'): self.c_unit[j].close()
            #------------------------------------------------------------
##            if (self.Ks_file[j]  != ''): self.Ks_unit[j].close()        
##            if (self.Ki_file[j]  != ''): self.Ki_unit[j].close()
##            if (self.qs_file[j]  != ''): self.qs_unit[j].close()
##            if (self.qi_file[j]  != ''): self.qi_unit[j].close()
##            if (self.qr_file[j]  != ''): self.qr_unit[j].close()
##            if (self.pB_file[j]  != ''): self.pB_unit[j].close()
##            if (self.pA_file[j]  != ''): self.pA_unit[j].close()        
##            if (self.lam_file[j] != ''): self.lam_unit[j].close()
##            if (self.c_file[j]   != ''): self.c_unit[j].close()
            
    #   close_input_files()
    #------------------------------------------------------------------- 
    def write_richards_1d_cfg_file(self, in_directory=None,
                                   case_prefix=None,
                                   n_layers=3,
                                   dt=float64(60),
                                   soil_types = ['silty_loam',
                                                 'loam',
                                                 'clay_loam']):

        if (in_directory is None):
            in_directory = self.in_directory
        if (case_prefix is None):
            case_prefix = self.case_prefix      

        n_soils = len(soil_types)

        #----------------------
        # Open a new CFG file
        #----------------------
        cfg_file = (in_directory + case_prefix + '_infil_R1D.cfg')
        cfg_unit = open(cfg_file, 'w')
        
        #------------------------------
        # Write method, timestep, etc.
        #------------------------------
        format = '%-23s %-14s %-23s %-14s\n'
        top = [
            '\n',
            '--------------------------------------------------\n',
            '  Infiltration Process Variables \n',
            '--------------------------------------------------\n',
            (format % ('Method code:', '4', ' ', ' ')),
            (format % ('Method name:', 'Richards_1D', ' ', ' ')),
            (format % ('Number of layers:', str(n_layers), ' ', ' ')),
            (format % ('Time step:',       'Scalar', str(dt), '[sec]'))  ]
        cfg_unit.writelines(top)
        
        #------------------------------------------------
        # Write (scalar) soil parameters for each layer
        #------------------------------------------------
        for k in range(n_layers):
            #-------------------------------------------
            # If (n_layers > n_soils), repeat last one
            #-------------------------------------------
            if (k > (n_soils-1)):  k = (n_layers-1)
            soil = soil_base.soil_base(soil_types[k])
            soil.initialize()
            middle = [
                (format % ('Ks:',               'Scalar', str(soil.K_s),     '[m/s]')),
                (format % ('Ki:',               'Scalar', str(soil.K_i),     '[m/s]')),
                (format % ('qs:',               'Scalar', str(soil.theta_s), '[none]')),
                (format % ('qi:',               'Scalar', str(soil.theta_i), '[none]')),
                (format % ('qr:',               'Scalar', str(soil.theta_r), '[none]')),
                (format % ('pB:',               'Scalar', str(soil.psi_B),   '[m]')),
                (format % ('pA:',               'Scalar', str(soil.psi_A),   '[m]')),
                (format % ('lambda:',           'Scalar', str(soil.Lambda),  '[none]')),
                (format % ('c:',                'Scalar', str(soil.c),       '[none]')),
                (format % ('dz:',               'Scalar', str(soil.dz),      '[m]')),
                (format % ('nz:',               'Scalar', str(soil.nz),      '[none]')),
                (format % ('Closest soil_type:', soil_types[k], ' ', ' '))  ]           
            cfg_unit.writelines(middle)
        
        #-------------------------------
        # Write default output options
        #-------------------------------
        cp = case_prefix
        bottom = [
            (format % ('Save grid timestep:',    'Scalar', '60.00000000',       '[sec]')),
            (format % ('Save v0 grids:',         '0',       cp + '_2D-v0.nc',  '[m/s]')),
            (format % ('Save q0 grids:',         '0',       cp + '_2D-q0.nc',  '[none]')),
            (format % ('Save I  grids:',         '0',       cp + '_2D-I.nc',   '[m]')),
            (format % ('Save Zw grids:',         '0',       cp + '_2D-Zw,nc',  '[m]')),
            (format % ('Save pixels timestep:',  'Scalar',  '60.00000000',      '[sec]')),
            (format % ('Save v0 pixels:',        '0',       cp + '_0D-v0.nc',  '[m/s]')),
            (format % ('Save q0 pixels:',        '0',       cp + '_0D-q0.nc',  '[none]')),
            (format % ('Save I  pixels:',        '0',       cp + '_0D-I.nc',   '[m]')),
            (format % ('Save Zw pixels:',        '0',       cp + '_0D-Zw.nc',  '[m]')),
            (format % ('Save cube timestep:',    'Scalar',  '60.00000000',      '[sec]')),
            (format % ('Save q cubes:',          '0',       cp + '_3D-q.nc',   '[none]')),
            (format % ('Save p cubes:',          '0',       cp + '_3D-p.nc',   '[m]')),
            (format % ('Save K cubes:',          '0',       cp + '_3D-K.nc',   '[m/s]')),
            (format % ('Save v cubes:',          '0',       cp + '_3D-v.nc',   '[m/s]')),
            (format % ('Save profile timestep:', 'Scalar',  '60.00000000',      '[sec]')),
            (format % ('Save q profiles:',       '0',       cp + '_1D-q.nc',   '[none]')),
            (format % ('Save p profiles:',       '0',       cp + '_1D_p.nc',   '[m]')),
            (format % ('Save K profiles:',       '0',       cp + '_1D_K.nc',   '[m/s]')),
            (format % ('Save v profiles:',       '0',       cp + '_1D_v.nc',   '[m/s]'))  ]
        cfg_unit.writelines(bottom)

        #-----------------
        # Close the file
        #-----------------
        cfg_unit.close()
        print('Finished writing new Richard_1D CFG file to:')
        print('   ' + cfg_file)
        print(' ')
        
    #   write_richards_1d_cfg_file()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def build_layered_var(self, v_by_layer):

        #-----------------------------------------------------------
        # Notes:  This routine examines user-selected parameters
        #         for each soil layer.  If all layers have a
        #         scalar value, then a 1D array (a z-profile) is
        #         constructed for this variable.  Due to IDL's
        #         dynamic data typing, it will get used correctly
        #         by the "Update_Richards" routines.  If any layer
        #         has a 2D value, then a 3D array is constructed
        #         for this variable.

        #         Note that self.nz was previously set to the sum:
        #            long(total(self.nz_val))
        #-----------------------------------------------------------

        #-----------------------------
        # Create an array of indices
        #-----------------------------
        #i[0] = 0
        #i[1] = self.nz_val[0]
        #i[2] = self.nz_val[0] + self.nz_val[1]
        #etc.
        #----------------------------------------------
        i = concatenate(([int32(0)], int32(cumsum(self.nz_val))) )
        
        #----------------------------------------
        # Do all layers have a scalar parameter
        # value for this particular variable ??
        #----------------------------------------
        nmax = int16(1)
        for j in range(self.n_layers):
            nj = size(v_by_layer[j])
            nmax = maximum(nmax, nj)
        ALL_SCALARS = (nmax == 1)
        
        #-------------------------------------------
        # Build a "data cube" from layer variables
        #-------------------------------------------
        if (ALL_SCALARS):
            #----------------------------------------------
            # All layers have a scalar value for this var
            #----------------------------------------------
            var = zeros([self.nz], dtype='Float64')
            for j in range(self.n_layers):
                var[i[j]: i[j + 1]] = v_by_layer[j]
        else:    
            #--------------------------------------------------------
            # Note that all nz "levels" in a given layer can be
            # initialized to a grid, but not yet to different grids
            #--------------------------------------------------------
            var = zeros([self.nz, self.ny, self.nx], dtype='Float64')
            for j in range(self.n_layers):
                for k in range(i[j], i[j + 1]):
                    var[k,:,:] = v_by_layer[j]
                
                #----------------------------------------------------
                # Next line doesn't work if v_by_layer[j] is a grid
                #----------------------------------------------------
                #  var[*, *, i[j]:i[j+1]-1 ] = v_by_layer[j]

        return var

    #   build_layered_var
    #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def Theta_TBC(psi_cm, theta_s, theta_r, psi_B, psi_A, c, Lambda,
              REPORT=False):

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
    psi_m = (psi_cm / float64(100))    # [cm -> meters]
    
    ratio = (psi_m + psi_A) / psi_B    # (should be > 0)
    
    theta = (float64(1) + ratio ** c) ** (-Lambda / c)
    theta = theta * (theta_s - theta_r) + theta_r
    
    #------------------
    # Optional report
    #------------------
    if (REPORT):    
        print('theta_s = ', theta_s)
        print('theta   = ', theta)
        print('theta_r = ', theta_r)
        print(' ')
    
    return theta
    
#   Theta_TBC()
#-----------------------------------------------------------------------
def K_of_Theta_TBC(theta, K_s, theta_s, theta_r, Lambda,
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
    eta = (float64(2) + (float64(3) * Lambda))
    eps = eta / Lambda
    
    #--------------------------------------
    # Compute the "relative conductivity"
    #--------------------------------------
    K_r = ((theta - theta_r) / (theta_s - theta_r)) ** eps
    
    #-----------------------------
    # Compute K from K_s and K_r
    #-----------------------------
    K_r = maximum((minimum(K_r, 1.0)), 0.0)
    K = K_s * K_r
    
    #------------------
    # Optional report
    #------------------
    if (REPORT):    
        print('K = ', K[0:4])
        # print ' '
    
    return K
    
#   K_of_Theta_TBC()
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def Z_Derivative_1D(v, dz):

    #----------------------------------------------------------
    # Notes:  v is a 1D array and dz is a scalar or 1D array.
    #         The result is a 1D array, same size as v.

    #        This function does not worry about the wrap
    #        around affect of ROLL at bottom.  This must
    #        be handled by the caller.
    #----------------------------------------------------------
    v_below = np.roll(v, -1, axis=0) 
    dv_dz   = (v_below - v) / dz
    
    return dv_dz
    
#   Z_Derivative_1D()
#-----------------------------------------------------------------------
def Z_Derivative_3D(v, dz):

    #------------------------------------------------------------
    # Notes:  v is a 3D array (or data cube) and dz is a scalar
    #         or 1D array.  The result is a 3D array, same size
    #         as v.

    #         This function does not worry about the wrap
    #         around affect of ROLL at bottom.  This must
    #         be handled by the caller.
    #------------------------------------------------------------
    n_dz = size(dz)
    v_below = np.roll(v, -1, axis=0)
    
    if (n_dz == 1):    
        #-----------------
        # dz is a scalar
        #-----------------
        dv_dz = (v_below - v) / dz
    else:    
        dv_dz = (v_below - v)
        for j in range(n_dz):
            dv_dz[j,:,:] = dv_dz[j,:,:] / dz[j]
    
    return dv_dz
    
#   Z_Derivative_3D()
#-----------------------------------------------------------------------

