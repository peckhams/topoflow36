"""
# This class defines a hydrologic infiltration component that numerically
solves the 1D version of the Richards Equation.  That is, while each grid
cell can have its own profile, the horizontal/lateral components of the
Darcy velocity field are assumed to be negligible.

This class inherits from the infiltration "base class" in "infil_base.py".

Note:  The numerical scheme used by this component cannot handle sharp
changes (e.g. discontinuities) in hydraulic conductivity.

See: Smith, R.E. (2002) Infiltration Theory for Hydrologic Applications,
Water Resources Monograph 15, AGU.
"""
## Copyright (c) 2001-2019, Scott D. Peckham
##
## January 2013   (Revised handling of input/output names).
## October 2012   (CSDMS Standard Names and BMI)
## January 2009  (converted from IDL)
## May, August 2009
## May 2010  (changes to unit_test() and read_cfg_file()
## June 2010 (Bug fix: Added qH_list and eta_list in
##            set_computed_input_vars(). Unit test. )
## November 2010 (New approach to BCs and update_theta().)

#---------------------------------------------------------------------
#
#  unit_test()
#
#  class infil_component         # (inherits from infil_base.py)
#
#      get_component_name()
#      get_attribute()           # (10/26/11)
#      get_input_var_names()     # (10/25/12)
#      get_output_var_names()    # (10/25/12)
#      get_var_name()            # (10/25/12)
#      get_var_units()           # (10/25/12)
#      ----------------------------
#      initialize_layer_vars()
#      set_computed_input_vars()
#      check_input_types()
#      initialize_computed_vars()
#      ----------------------------
#      initialize_theta_r()
#      initialize_theta_i()
#      initialize_K_i()
#      ----------------------------
#      update()
#      update_surface_BC()
#      update_bottom_BC()
#      update_theta()
#      check_theta()    # (6/29/10)
#      update_psi()
#      update_K()
#      update_v()
#      update_Zw()
#      ----------------------
#      update_infil_rate()
#      update_Rg()
#      update_q0()
#      ----------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      ------------------------------
#      build_layered_var()    # (Moved into infil_base.py)

#  Functions:
#      Theta_TBC()          (still used)
#      K_of_Theta_TBC()     (used by initialize_K_i())
#      Z_Derivative_1D()    (Mar 2008)
#      Z_Derivative_3D()    (Mar 2007)
#      Z_Forward_Average()
#      Z_Backward_Average()
#
#-----------------------------------------------------------------------

import numpy as np
import os, sys

from topoflow.components import infil_base
from topoflow.components import soil_base

from topoflow.utils import model_input
from topoflow.utils import tf_utils  ## (for unit_test only)
from topoflow.utils import soil_trans_BC as stbc

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
    _att_map = {
        'model_name':         'TopoFlow_Infiltration_Richards_1D',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'InfilRichards1D',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Infil_Richards_1D.cfg.in',
        'cfg_extension':      '_infil_richards_1d.cfg',
        'cmt_var_prefix':     '/InfilRichards1D/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Infil_Richards_1D.xml',
        'dialog_title':       'Infiltration: Richards 1D Parameters',
        'time_units':         'seconds' }
    
    _input_var_names = [
        'atmosphere_water__rainfall_volume_flux',           # (P_rain)  
        'glacier_ice__melt_volume_flux',                    # (MR)
        'land_surface__elevation',                          # (elev)
        'land_surface_water__evaporation_volume_flux',      # (ET)
        'snowpack__melt_volume_flux',                       # (SM)
        'soil_water_sat-zone_top_surface__elevation' ]      # (h_table)               

    #-----------------------------------------------------------
    # We use "bubbling_pressure_head" vs. "air_entry_pressure"
    # because base quantity is "head" with units of length.
    # We could use "air_entry_pressure_head", though.
    #-----------------------------------------------------------
    # We may want to make some of these available by layer.
    #-------------------------------------------------------------------
    # lambda = brooks-corey pore-size distribution parameter
    # b      = brooks-corey pore-size distribution index" = 1 / lambda
    # c      = brooks-corey-smith pore connectedness index
    #        = (eta / lambda) = (2*b + 3)
    #-------------------------------------------------------------------
    # See infil_base.set_constants() for how these constants are set:
    #   psi_oven_dry, psi_air_dry, psi_min, psi_hygro,
    #   psi_wilt and psi_field.
    #-------------------------------------------------------------------
    # porosity is set in soil_base.py.  #####################
    #-------------------------------------------------------------------           
    _output_var_names = [
        'model__time_step',                                # dt
        # 'model_grid_cell__area',                         # da
        # 'soil__porosity',                                # phi
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux', # vol_IN
        'soil_surface_water__infiltration_volume_flux',    # IN
        'soil_surface_water__volume_fraction',             # q0
        # 'soil_water__brooks-corey_b_parameter',          # b
        'soil_water__brooks-corey_eta_parameter',          # eta
        'soil_water__brooks-corey_lambda_parameter',       # lam
        'soil_water__brooks-corey-smith_c_parameter',      # c
        'soil_water__brooks-corey-smith_pressure_head_offset_parameter',  # pA
        'soil_water__bubbling_pressure_head',              # pB
        # 'soil_water__field-capacity_volume_fraction',    # qf  ######## CHECK
        'soil_water__hydraulic_conductivity',              # K
        'soil_water__hygroscopic_volume_fraction',         # qH
        'soil_water__initial_hydraulic_conductivity',      # Ki
        'soil_water__initial_volume_fraction',             # qi
        # 'soil_water__normalized_volume_fraction',        # S_eff
        'soil_water__pressure_head',                       # p
        # 'soil_water__relative_hydraulic_conductivity',   # K_rel
        'soil_water__residual_volume_fraction',            # qr
        'soil_water__saturated_hydraulic_conductivity',    # Ks
        'soil_water__saturated_volume_fraction',           # qs
        'soil_water__volume_fraction',                     # q
        # 'soil_water__wilting-point_volume_fraction',     # qw
        'soil_water_flow__z_component_of_darcy_velocity',  # v
        'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux',  # vol_Rg
        'soil_water_sat-zone_top__recharge_volume_flux',   # Rg
        'soil_water_wetting-front__depth' ]                # Zw

    _var_name_map = {
        'atmosphere_water__rainfall_volume_flux':          'P_rain',
        'glacier_ice__melt_volume_flux':                   'MR',
        'land_surface__elevation':                         'elev',
        'land_surface_water__evaporation_volume_flux':     'ET',
        'snowpack__melt_volume_flux':                      'SM',
        'soil_water_sat-zone_top_surface__elevation':      'h_table',
        #--------------------------------------------------------------
        'model__time_step':                                'dt',
        ## 'model_grid_cell__area':                        'da', 
        # 'soil__porosity':                                'phi',
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux': 'vol_IN',
        'soil_surface_water__infiltration_volume_flux':    'IN',
        'soil_surface_water__volume_fraction':             'q0',
        # 'soil_water__brooks-corey_b_parameter':          'b',
        'soil_water__brooks-corey_eta_parameter':          'eta',
        'soil_water__brooks-corey_lambda_parameter':       'lam',
        'soil_water__brooks-corey-smith_c_parameter':      'c',
        'soil_water__brooks-corey-smith_pressure_head_offset_parameter': 'pA',      
        'soil_water__bubbling_pressure_head':              'pB',
        # 'soil_water__field-capacity_volume_fraction':    'qf',  ######  CHECK
        'soil_water__hydraulic_conductivity':              'K',
        'soil_water__hygroscopic_volume_fraction':         'qH',
        'soil_water__initial_hydraulic_conductivity':      'Ki',
        'soil_water__initial_volume_fraction':             'qi',
        # 'soil_water__normalized_volume_fraction',        'S_eff',
        'soil_water__pressure_head':                       'p',
        # 'soil_water__relative_hydraulic_conductivity':   'K_rel',
        'soil_water__residual_volume_fraction':            'qr',
        'soil_water__saturated_hydraulic_conductivity':    'Ks',
        'soil_water__saturated_volume_fraction':           'qs',
        'soil_water__volume_fraction':                     'q',
        # 'soil_water__wilting-point_volume_fraction':     'qw',
        'soil_water_flow__z_component_of_darcy_velocity':  'v',
        'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux': 'vol_Rg',
        'soil_water_sat-zone_top__recharge_volume_flux':   'Rg',
        'soil_water_wetting-front__depth':                 'Zw' }

    _var_units_map = {
        'atmosphere_water__rainfall_volume_flux':          'm s-1',   
        'glacier_ice__melt_volume_flux':                   'm s-1',
        'land_surface__elevation':                         'm',
        'land_surface_water__evaporation_volume_flux':     'm s-1',
        'snowpack__melt_volume_flux':                      'm s-1',
        'soil_water_sat-zone_top_surface__elevation':      'm',
        #------------------------------------------------------------
        'model__time_step': 's',     ############## CHECK
        # 'model_grid_cell__area': 'm2',
        # 'soil__porosity': '1',
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux': 'm3',
        'soil_surface_water__infiltration_volume_flux': 'm s-1',
        'soil_surface_water__volume_fraction': '1',
        # 'soil_water__brooks-corey_b_parameter': '1', 
        'soil_water__brooks-corey_eta_parameter': '1',
        'soil_water__brooks-corey_lambda_parameter': '1',
        'soil_water__brooks-corey-smith_c_parameter': '1',
        'soil_water__brooks-corey-smith_pressure_head_offset_parameter': 'm',     
        'soil_water__bubbling_pressure_head': 'm',
        # 'soil_water__field-capacity_volume_fraction': '1',
        'soil_water__hydraulic_conductivity': 'm s-1',
        'soil_water__initial_hydraulic_conductivity': 'm s-1',
        'soil_water__initial_volume_fraction': '1',
        'soil_water__pressure_head': 'm',
        'soil_water__hygroscopic_volume_fraction': '1',
        # 'soil_water__normalized_volume_fraction': '1',
        # 'soil_water__relative_hydraulic_conductivity': '1',
        'soil_water__residual_volume_fraction': '1',
        'soil_water__saturated_hydraulic_conductivity': 'm s-1',
        'soil_water__saturated_volume_fraction': '1',
        'soil_water__volume_fraction': '1',
        # 'soil_water__wilting-point_volume_fraction': '1',
        'soil_water_flow__z_component_of_darcy_velocity': 'm s-1',
        'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux': 'm3',
        'soil_water_sat-zone_top__recharge_volume_flux': 'm s-1',
        'soil_water_wetting-front__depth': 'm' }
        
    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Infiltration_Richards_1D'

    #   get_component_name() 
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print('###################################################')
            print(' ERROR: Could not find attribute: ' + att_name)
            print('###################################################')
            print(' ')

    #   get_attribute()
    #-------------------------------------------------------------------
    def get_input_var_names(self):

        #--------------------------------------------------------
        # Note: These are currently variables needed from other
        #       components vs. those read from files or GUI.
        #--------------------------------------------------------   
        return self._input_var_names
    
    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):
 
        return self._output_var_names
    
    #   get_output_var_names()
    #-------------------------------------------------------------------
    def get_var_name(self, long_var_name):
            
        return self._var_name_map[ long_var_name ]

    #   get_var_name()
    #-------------------------------------------------------------------
    def get_var_units(self, long_var_name):

        return self._var_units_map[ long_var_name ]
   
    #   get_var_units()
    #-------------------------------------------------------------------
##    def get_var_type(self, long_var_name):
##
##        #---------------------------------------
##        # So far, all vars have type "double",
##        # but use the one in BMI_base instead.
##        #---------------------------------------
##        return 'float64'
##    
##    #   get_var_type()
    #-------------------------------------------------------------------
    def initialize_layer_vars(self):

        #---------------------------------------------------
        # Notes:  The sequence of function calls is:
        #
        # initialize()                  # in BMI_base.py
        #    initialize_config_vars()   # in BMI_base.py
        #       read_config_file()      # in BMI_base.py
        #           initialize_layer_vars()   # here
        #           (after n_layers was read)
        #       set_computed_input_vars()     # below
        #----------------------------------------------------
        n_layers = self.n_layers
        
        #-------------------------------------------------
        # Get arrays to store soil params for each layer
        #-------------------------------------------------
        self.soil_type = np.zeros(n_layers, dtype='<U200')
        self.dz_val    = np.zeros(n_layers, dtype='float64')    #### + dz3
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
        # Actual variable arrays built with build_layered_var.
        #---------------------------------------------------------
        # (5/19/10) We need Ks_list vs Ks here, since we use
        # these to build one big, 3D Ks array.
        #---------------------------------------------------------     
        self.Ks_list  = list( np.zeros(n_layers, dtype='float64') )
        self.Ki_list  = list( np.zeros(n_layers, dtype='float64') )
        self.qs_list  = list( np.zeros(n_layers, dtype='float64') )
        self.qi_list  = list( np.zeros(n_layers, dtype='float64') )
        self.qr_list  = list( np.zeros(n_layers, dtype='float64') )
        self.pB_list  = list( np.zeros(n_layers, dtype='float64') )
        self.pA_list  = list( np.zeros(n_layers, dtype='float64') )
        self.lam_list = list( np.zeros(n_layers, dtype='float64') )
        self.c_list   = list( np.zeros(n_layers, dtype='float64') )
        #------------------------------------------------
        # Note:  These two are computed from the others
        #------------------------------------------------
        self.eta_list = list( np.zeros(n_layers, dtype='float64') )
        self.qH_list  = list( np.zeros(n_layers, dtype='float64') )
       
    #   initialize_layer_vars()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        self.RICHARDS = True
        self.G_file    = ''   # (still need to be defined)
        self.gam_file  = ''
        self.G_type    = 'Scalar'     # (see below)
        self.gam_type  = 'Scalar'

        #------------------------------------------------------------
        # Compute eta value for each soil layer from lambda values.
        # Depending on lambda, eta values will be scalars or grids.
        #------------------------------------------------------------
        for j in range(self.n_layers):
#             print('type(lam_list[j]) = ' + str(type(self.lam_list[j])) )
#             print('lam_list[j].dtype = ' + str(self.lam_list[j].dtype) )
#             print('lam_list[j].shape = ' + str(self.lam_list[j].shape) )
            self.eta_list[j] = np.float64(2) + (np.float64(3) * self.lam_list[j])
                             
        #-------------------------------------------------------------
        # Compute a qH value for each soil layer from other values
        # using the theta_of_psi() function.  qH values will be
        # scalars or grids, depending on the args to theta_of_psi().
        #-------------------------------------------------------------
        for j in range(self.n_layers):
            self.qH_list[j] = stbc.theta_of_psi( self.psi_hygro, \
                                   self.qs_list[j], self.qr_list[j], \
                                   self.pB_list[j], self.pA_list[j], \
                                   self.c_list[j],  self.lam_list[j] )
#             self.qH_list[j] = Theta_TBC( self.psi_hygro, \
#                                          self.qs_list[j], self.qr_list[j], \
#                                          self.pB_list[j], self.pA_list[j], \
#                                          self.c_list[j],  self.lam_list[j] )

        #------------------------------------------------------------
        # 2019-10-30  #############################
        # Compute a Ki value for each soil layer from other values
        # using the K_of_theta() function.  Ki values will be
        # scalars or grids, depending on the args to K_of_theta().
        #------------------------------------------------------------
        # This ensures that K_i is consistent with theta_i.
        # But need to make sure K[0] = Ki.  ##################
        #------------------------------------------------------------        
        for j in range(self.n_layers):
            self.Ki_list[j] = stbc.K_of_theta( self.qi_list[j], \
                                   self.Ks_list[j], self.qs_list[j], \
                                   self.qr_list[j], self.lam_list[j] )
                 
            if (self.DEBUG):
                print('min (Ki_list[j]) =', self.Ki_list[j].min() )
                print('max (Ki_list[j]) =', self.Ki_list[j].max() )
                                                                  
        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt    = np.maximum(self.save_grid_dt,    self.dt)
        self.save_pixels_dt  = np.maximum(self.save_pixels_dt,  self.dt)
        self.save_profile_dt = np.maximum(self.save_profile_dt, self.dt)
        self.save_cube_dt    = np.maximum(self.save_cube_dt,    self.dt)
        
    #   set_computed_input_vars()    
    #-------------------------------------------------------------------
    def check_input_types(self):

        #----------------------------------------------------
        # Notes: ET is often a 2D grid even when the others
        #        are scalars.  See how P_total is defined
        #        in update_surface_influx().
        #----------------------------------------------------
        are_scalars = np.array([
                         self.is_scalar('P_rain'),
                         self.is_scalar('SM'),
                         self.is_scalar('ET'),  #########
                         #----------------------------------
                         self.is_scalar('Ks_list[0]'),
                         self.is_scalar('Ki_list[0]'),
                         self.is_scalar('qs_list[0]'),
                         self.is_scalar('qi_list[0]'),
                         self.is_scalar('qr_list[0]'),
                         self.is_scalar('pB_list[0]'),
                         self.is_scalar('pA_list[0]'),
                         self.is_scalar('c_list[0]'),
                         self.is_scalar('lam_list[0]')])

        self.ALL_SCALARS = np.all(are_scalars)
        
        #----------------------------------------
        # Use the same profile for all pixels ?
        #---------------------------------------------
        # NB! This var only used by Richards' method
        # so it shouldn't appear in "infil_base.py".
        #---------------------------------------------
        self.SINGLE_PROFILE = self.ALL_SCALARS  # (3/19/07)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):
     
        dtype = 'float64'
        self.vol_IN = self.initialize_scalar( 0, dtype=dtype)
        self.vol_Rg = self.initialize_scalar( 0, dtype=dtype)

        #---------------------------------------
        # Get surface influx to initialize "v"
        #---------------------------------------
        self.update_surface_influx()
        
        #-----------------------------------------------------
        # Compute dz as 1D array from scalars in self.dz_val      
        #-----------------------------------------------------
        # Compute the z-vector, for plotting profiles
        #----------------------------------------------
        ### self.nz = np.sum( self.nz_val )
        self.build_layer_z_vector()

        #------------------------------------------------
        # Now build a 1D or 3D array for each input var
        #--------------------------------------------------------
        # (3/12/08) Same code should work if (self.n_layers eq 1)
        #--------------------------------------------------------
        self.Ks  = self.build_layered_var( self.Ks_list )
        self.Ki  = self.build_layered_var( self.Ki_list )
        self.qs  = self.build_layered_var( self.qs_list )
        self.qi  = self.build_layered_var (self.qi_list )
        self.qr  = self.build_layered_var( self.qr_list )
        self.pB  = self.build_layered_var( self.pB_list )
        self.pA  = self.build_layered_var( self.pA_list )
        self.lam = self.build_layered_var( self.lam_list )
        self.c   = self.build_layered_var( self.c_list )
        #--------------------------------------------------
        # Note:  eta_list and qH_list are computed from
        #        the others in set_computed_input_vars().      
        #--------------------------------------------------
        self.eta = self.build_layered_var( self.eta_list )
        self.qH  = self.build_layered_var( self.qH_list )

        #--------------
        # For testing
        #--------------
        if (self.DEBUG):
            print('In initialize_computed_vars():')
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
 
        #-------------------------------------------------------
        # Note: qi and Ki are created with build_layered_var()
        #-------------------------------------------------------
        # Note: Don't need to initialize Rg here.
        #-------------------------------------------------------          
        if (self.ALL_SCALARS):
            #----------------------------------
            # Infiltration varies with z only
            #----------------------------------
            self.q = np.zeros(self.nz, dtype=dtype) + self.qi
            self.p = np.zeros(self.nz, dtype=dtype)
            self.K = np.zeros(self.nz, dtype=dtype) + self.Ki
            self.v = np.zeros(self.nz, dtype=dtype)
            #---------------------------------------------------------
            self.IN = np.float64(0)   # (infil. rate at surface)
            self.I  = np.float64(0)   # (total infil. depth)
            self.Zw = np.float64(0)   # (wetting front depth)
            # self.Rg = np.float64(0)   # (infil. rate at water table)
            #---------------------------------------------------------
            # Initialize I to 1e-6 to avoid divide by zero at start?
            #---------------------------------------------------------
##            self.I  = np.float64(1e-6)   # (total infil. depth)
##            self.Zw = np.float64(1e-6)   # (wetting front depth)

##            if (self.DEBUG):
##                print 'shape(self.v)  =', np.shape(self.v)
##                print 'shape(P_total) =', np.shape(self.P_total)
##                print 'type(P_total)  =', type(self.P_total)
                
            #--------------------------------------------
            # Set BC at the surface (done elsewhere ??)
            #--------------------------------------------
            self.v[0] = self.P_total
    
        else:
            #------------------------------------
            # Infiltration varies with x, y & z
            #------------------------------------
            self.q  = np.zeros((self.nz, self.ny, self.nx), dtype=dtype)
            self.p  = np.zeros((self.nz, self.ny, self.nx), dtype=dtype) 
            self.K  = np.zeros((self.nz, self.ny, self.nx), dtype=dtype) 
            self.v  = np.zeros((self.nz, self.ny, self.nx), dtype=dtype)
            #---------------------------------------------------------------
            self.IN = np.zeros([self.ny, self.nx], dtype=dtype)
            self.I  = np.zeros([self.ny, self.nx], dtype=dtype)
            self.Zw = np.zeros([self.ny, self.nx], dtype=dtype)
            # self.Rg = np.zeros([self.ny, self.nx], dtype=dtype)

            #--------------------------------------
            # Initialize q to qi (qi is 1D or 3D)
            #--------------------------------------

            if (np.size(self.qi) == self.nz):
                for j in range(self.nz):
                    self.q[j,:,:] = self.qi[j]
                # (Can this be done with array operators instead ?)
            else:
                self.q += self.qi
            
            if (self.DEBUG):
                print('Initialized theta to theta_i.')
                print('   min(qi) = ', self.qi.min() )
                print('   max(qi) = ', self.qi.max() )
                print('   min(q)  = ', self.q.min()  )
                print('   max(q)  = ', self.q.max()  )

            #--------------------------------------
            # Initialize K to Ki (Ki is 1D or 3D)
            #--------------------------------------

            if (np.size(self.Ki) == self.nz):
                for j in range(self.nz):
                    self.K[j,:,:] = self.Ki[j]
                # (Can this be done with array operators instead ?)
            else:
                self.K += self.Ki
               
            if (self.DEBUG): 
                print('Initialized K to K_i.')
                print('   min(Ki) = ', self.Ki.min() )
                print('   max(Ki) = ', self.Ki.max() )
                print('   min(K)  = ', self.K.min()  )
                print('   max(K)  = ', self.K.max()  )
        
            #--------------------------------------------------
            # If q is now 3D, convert qs to 3D also so we can
            # compute (q - qs) in update_v(). (6/22/10)
            #--------------------------------------------------
            ## if (np.ndim(self.qs) == 1):
            if (np.size(self.qs) == self.nz):
                temp = self.qs.copy()
                self.qs = np.zeros((self.nz, self.ny, self.nx), dtype=dtype)
                for j in range(self.nz):
                    self.qs[j,:,:] = temp[j]   #######

            #--------------------------------------------------
            # If q is now 3D, convert qH to 3D also so we can
            # compute (q - qH) in ******(). (6/22/10)
            #--------------------------------------------------
            ## if (np.ndim(self.qH) == 1):
            if (np.size(self.qH) == self.nz):
                temp = self.qH.copy()
                self.qH = np.zeros((self.nz, self.ny, self.nx), dtype=dtype)
                for j in range(self.nz):
                    self.qH[j,:,:] = temp[j]   #######
                    
            #--------------------------------------------
            # Set BC at the surface (done elsewhere ??)
            #--------------------------------------------
            self.v[0,:,:] = self.P_total

        ##############################################
        # Set initial values of psi and v (11/12/10)
        ##############################################
        self.update_psi()
        self.update_v()
        
        #-------------------------------------------------
        # Print some suggested (i.e. consistent) values
        # for theta_r, theta_i and K_i.  (10/12/10)
        #-------------------------------------------------
        ## self.print_suggested_values()
        
        ###########################################
        # Override some of the user's settings ??
        ###########################################
##        if (self.SPECIAL_DEFAULTS):
##            self.initialize_theta_r()
##            self.initialize_theta_i()
##            self.initialize_K_i()
            
    #   initialize_computed_vars()
    #-------------------------------------------------------------------
#     def initialize_theta_r(self):
# 
#         #-------------------------------------------------
#         # Note that this is not entirely consistent with
#         # the Theta_TBC() function, but that function
#         # requires theta_r as an argument.
#         #-------------------------------------------------
#         # Initialize theta_r to the min allowed value.
#         #-------------------------------------------------
#         psi_r = self.psi_min
#         
#         #--------------------------------------
#         # Note:  Both psi's < 0, so ratio > 0
#         #--------------------------------------
#         self.qr = self.qs * (self.pB / psi_r)**self.lam
#     
#     #   initialize_theta_r()
#     #-------------------------------------------------------------------
#     def initialize_theta_i(self):
# 
#         #------------------------------------------------
#         # Initialize theta_i = qi to the field capacity.
#         # Be sure to call initialize_theta_r() first.
#         #------------------------------------------------
#         self.qi = stbc.theta_of_psi( self.psi_field, \
#                              self.qs, self.qr, \
#                              self.pB, self.pA, \
#                              self.c,  self.lam )
# #         self.qi = Theta_TBC( self.psi_field, \
# #                              self.qs, self.qr, \
# #                              self.pB, self.pA, \
# #                              self.c,  self.lam )
#                                                 
#     #   initialize_theta_i()
    #-------------------------------------------------------------------
    def initialize_K_i(self):

        self.Ki = stbc.K_of_theta( self.qi, self.Ks, self.qs,
                                   self.qr, self.lam )
                                  
#         self.Ki = K_of_Theta_TBC( self.qi, self.Ks, self.qs,
#                                   self.qr, self.lam )

    #   initialize_K_i()
    #-------------------------------------------------------------------
    def print_suggested_values(self):

        if (self.DEBUG):
            print('Calling print_suggested_values()...')
            
        #-------------------------------------------------
        # theta_r is often set to theta_hygroscopic.
        # theta_i is often set to theta_field_capacity.
        #-------------------------------------------------
        print('=====================================================')
        for k in range(self.n_layers):

##            print 'Ks[k]  =', self.Ks_list[k]
##            print 'Ki[k]  =', self.Ki_list[k]    
##            print 'qs[k]  =', self.qs_list[k]
##            print 'qi[k]  =', self.qi_list[k]
##            print 'qr[k]  =', self.qr_list[k]
##            print 'pB[k]  =', self.pB_list[k]
##            print 'pA[k]  =', self.pA_list[k]
##            print 'lam[k] =', self.lam_list[k]
##            print 'c[k]   =', self.c_list[k]
##            print ' '
##            print 'psi_hygro =', self.psi_hygro, ' [m]'
##            print 'psi_field =', self.psi_field, ' [m]'
            
            #-------------------------------------------------
            # Compute this by analogy to equations 6-19 and
            # 6-20 in Dingman (2002), using theta_s instead
            # of porosity and recalling lambda = (1/b).
            #-------------------------------------------------
            # Note that this is not entirely consistent with
            # the Theta_TBC() function, but that function
            # requires theta_res as an argument.
            #-------------------------------------------------
##            psi_res   = self.psi_hygro
##            theta_sat = self.qs_list[k]
##            psi_B     = self.pB_list[k]
##            lam       = self.lam_list[k]
##            #--------------------------------------
##            # Note:  Both psi's < 0, so ratio > 0
##            #--------------------------------------
##            theta_res = theta_sat * (psi_B / psi_res)**lam

            #--------------------------------------------
            # If we trust theta_r, then do this instead
            #--------------------------------------------
            theta_res = self.qr_list[k]

            theta_hygro = stbc.theta_of_psi( self.psi_hygro,
                                     self.qs_list[k],
                                     self.qr_list[k],
                                     self.pB_list[k],
                                     self.pA_list[k],
                                     self.c_list[k],
                                     self.lam_list[k] )
            
            theta_init = stbc.theta_of_psi( self.psi_field,
                                    self.qs_list[k],
                                    theta_res,         #######
                                    self.pB_list[k],
                                    self.pA_list[k],
                                    self.c_list[k],
                                    self.lam_list[k] )

            K_init = stbc.K_of_theta( theta_init,       #######
                                      self.Ks_list[k],
                                      self.qs_list[k],
                                      theta_res,        #######
                                      self.lam_list[k] )
                                     
#             theta_hygro = Theta_TBC( self.psi_hygro,
#                                      self.qs_list[k],
#                                      self.qr_list[k],
#                                      self.pB_list[k],
#                                      self.pA_list[k],
#                                      self.c_list[k],
#                                      self.lam_list[k] )
#             
#             theta_init = Theta_TBC( self.psi_field,
#                                     self.qs_list[k],
#                                     theta_res,         #######
#                                     self.pB_list[k],
#                                     self.pA_list[k],
#                                     self.c_list[k],
#                                     self.lam_list[k] )
# 
#             K_init = K_of_Theta_TBC( theta_init,       #######
#                                      self.Ks_list[k],
#                                      self.qs_list[k],
#                                      theta_res,        #######
#                                      self.lam_list[k] )

            theta_r = self.qr_list[k]
            theta_i = self.qi_list[k]
            K_i     = self.Ki_list[k]
            print('Suggested initial values for layer', k+1, ':')
            ## print '   theta_r =', theta_res,  'vs.', theta_r
            print('   For theta_r =', theta_r)
            print('   theta_i =', theta_init, '   vs.', theta_i)
            print('   K_i     =', K_init,     'vs.', K_i)
            print('   theta_H =', theta_hygro, '  vs.', theta_r, ' (theta_r)')
            print(' ')
            
        print('===========================================================')                                                
        
    #   print_suggested_values()
    #-------------------------------------------------------------------
    def update(self, time_seconds=None):

        #################################
##        if (self.time_index > 570):
##            self.DEBUG = True
        #################################
        ## self.DEBUG = True
        
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
        self.update_surface_BC()
        self.update_bottom_BC()
        self.update_theta()
        self.update_psi()
        self.update_K()
        self.update_v()
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
        self.write_output_files()   # (Bug fix: 2019-10-27)
        ## self.write_output_files(time_seconds)

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        self.update_time()
        self.status = 'updated'  # (OpenMI 2.0 convention)
        if (self.DEBUG):
            print('Completed update() for time index =', self.time_index - 1)
            print(' ')

    #   update()
    #-------------------------------------------------------------------
    def update_surface_influx(self):

        if (self.DEBUG):
            print('Calling update_surface_influx()...')

        #------------------------------------        
        # These are now embedded references
        #------------------------------------
        P_rain = self.P_rain
        SM     = self.SM
        ET     = self.ET
        ## print 'min(ET), max(ET) =', ET.min(), ET.max()

        ###############################        
        # What if this is negative ?
        ###############################
        self.P_total = (P_rain + SM) - ET
        ### self.P_total = (P_rain + SM)

    #   update_surface_influx()
    #-----------------------------------------------------------------------
    def update_surface_BC(self, REPORT=False):

        #----------------------------------------------------------
        # Notes: Boundary conditions at the surface and at the
        #        bottom of the domain must be specified.  The
        #        approach used here is to specify values of psi,
        #        and to use theta(psi) (from TBC), to specify
        #        corresponding values of theta (and maybe K?).
        #----------------------------------------------------------
        # At the surface, the so-called "flux boundary condition"
        # is used prior to ponding (i.e. surface saturation).
        # We solve the following for psi[0]:
        #     v[0] = Kbar[0] * {1 - (psi[1] - psi[0])/dz} = r.
        # So we get:
        #
        #     psi[0] = {(r/Kbar[0]) - 1} * dz + psi[1].
        #
        # After ponding, we have psi[0] = 0.
        #----------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_surface_BC()...')

        if (self.dz.size == 1):
            dz = self.dz
        else:
            dz = self.dz[0]
        
        #----------------------------------
        # Set psi at the surface boundary
        #----------------------------------
        if (self.SINGLE_PROFILE):
            r = self.P_total
            if (self.q[0] < self.qs[0]):
                #------------------------------
                # Top layer is not saturated.
                #------------------------------
                Kbar      = (self.K[0] + self.K[1]) / 2.0  ################
                self.p[0] = ((r / Kbar) - 1) * dz + self.p[1]
            else:
                #--------------------------
                # Top layer is saturated.
                #--------------------------
                self.p[0] = 0.0
                Kbar      = r   # (just for report at end)
        else:
            p0   = self.p[0,:,:]
            p1   = self.p[1,:,:]
            Kbar = (self.K[0,:,:] + self.K[1,:,:]) / 2.0  ##############

            if (self.DEBUG):
                print('###  min(p0) = ', p0.min() )
                print('###  max(p0) = ', p0.max() )
                print('###  min(p1) = ', p1.min() )
                print('###  max(p1) = ', p1.max() )
                print('###  min(Kbar) = ', Kbar.min() )
                print('###  max(Kbar) = ', Kbar.max() )
                print('###  min(P_total) = ', self.P_total.min() )
                print('###  max(P_total) = ', self.P_total.max() )
                print()
                             
            #-------------------------------------
            # Where is top layer NOT saturated ?
            #-------------------------------------
            ## w1 = np.where( self.p[0,:,:] <  0)
            #-----------------------------------------------------
            # This makes w1 an array of True or False and should
            # be faster.  Don't need to check if w1 is empty.
            #-----------------------------------------------------
            w1 = ( self.q[0,:,:] <  self.qs[0,:,:] )
            if (self.P_total.size > 1):   ## BUG FIX: 2019-10-29
                r = self.P_total[ w1 ]
            else:
                r = self.P_total
            p0[w1] = ((r / Kbar[w1]) - 1) * dz + p1[w1]            
            #----------------------------------- 
            # This uses WHERE in the usual way
            #-----------------------------------           
#             w1 = np.where( self.q[0,:,:] <  self.qs[0,:,:] )
#             n1 = w1[0].size     
#             if (n1 != 0):
#                 r      = self.P_total[ w1 ]  ########
#                 p0[w1] = ((r / Kbar[w1]) - 1) * dz + p1[w1]

            #---------------------------------
            # Where is top layer saturated ?
            #---------------------------------
            ## w2 = np.where( self.p[0,:,:] >= 0 ) 
            #-------------------------------------------------- 
            ### w2 = ( self.q[0,:,:] >= self.qs[0,:,:] )
            w2 = np.invert( w1 )
            p0[w2] = 0.0
            #----------------------------------- 
            # This uses WHERE in the usual way
            #-----------------------------------                 
#             w2 = np.where( self.q[0,:,:] >= self.qs[0,:,:] )
#             n2 = w2[0].size
#             ## n2 = (self.rti.n_pixels - n1)                
#             if (n2 != 0):
#                 p0[w2] = 0.0

            #----------------------------------                
            # Set pressure head for top layer
            #----------------------------------
            self.p[0,:,:] = p0
            if (self.DEBUG):
                print('###  min(p0) = ', p0.min(), '(final)' )
                print('###  max(p0) = ', p0.max(), '(final)' )
            
        #------------------------
        # Check for stability ?
        #------------------------
        if (self.CHECK_STABILITY):
            if (self.SINGLE_PROFILE):
                wpos = (self.p[0] > 0)
            else:
                wpos = (self.p[0,:,:] > 0)
            npos = wpos.sum()
            if (npos > 0):
                print('############################################')
                print('ERROR in update_surface_BC():')
                print('   Some pressure head values are positive.')
                print('   Trying reducing timestep dt in CFG')
                print('   and check dz and nz for each layer.')
                print('############################################')
                print()
                self.DONE = True
                sys.exit()
                                                 
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
#             self.q[0] = Theta_TBC(psi, theta_s, theta_r, \
#                                   psi_B, psi_A, c, Lambda)
        else:
            #----------------------------------------
            # Now checking if ndim > 1 (2019-10-29)
            #----------------------------------------
            psi = self.p[0,:,:]    # [meters]
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
                Lambda = self.lam[0,:,:]
            else:
                Lambda = self.lam[0]
            #------------------------------
            self.q[0,:,:] = stbc.theta_of_psi(psi, theta_s, theta_r, \
                                      psi_B, psi_A, c, Lambda)
#             self.q[0,:,:] = Theta_TBC(psi, theta_s, theta_r, \
#                                       psi_B, psi_A, c, Lambda)

        #----------------
        # For debugging
        #----------------
        if (self.DEBUG and self.SINGLE_PROFILE):
            print('In update_surface_BC():')
            print('theta[0] =', self.q[0] )
            print('theta[1] =', self.q[1] )
            print('psi[0]   =', self.p[0])
            print('psi[1]   =', self.p[1])
            print('r        =', r )
            print('Kbar     =', Kbar )
            print('(r/Kbr - 1) =', ((r/Kbar) - 1))

    #   update_surface_BC()
    #-----------------------------------------------------------------------
    ## def update_bottom_BC(self, REPORT=False, BC='WATER_TABLE'):
    def update_bottom_BC(self, REPORT=False, BC='NO_FLOW'):
        
        #-----------------------------------------------------------
        # Notes: Boundary conditions at the surface and at the
        #        bottom of the domain must be specified.  The
        #        approach used here is to specify values of psi,
        #        and to use theta(psi) (from TBC), to specify
        #        corresponding values of theta (and maybe K?).
        #-----------------------------------------------------------
        # At the bottom, one of 3 BCs can be used:
        #
        # (1) gravity drainage:  no gradient in capillary pressure
        #                        (recall that H = (psi - z))
        #                           =>  d/dz(psi) = 0
        #                           =>  psi[n-1] = psi[n-2].
        #
        # (2) no flow:           v=0 => d/dz(psi) = 1
        #                           =>  psi[n-1] = psi[n-2] + dz
        #
        # (3) water table:       psi[n-1] = 0
        #
        # (4) fixed:             psi[n-1] = (initial value always)
        #
        #-----------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_bottom_BC()...')

        m = (self.nz - 1)
        
        #----------------------------------
        # Set psi at the bottom boundary
        #----------------------------------
        if (BC == 'WATER_TABLE'):
            if (self.SINGLE_PROFILE):
                self.p[m] = 0.0
            else:
                self.p[m,:,:] = 0.0
        elif (BC == 'GRAVITY_DRAINAGE'):
            if (self.SINGLE_PROFILE):
                self.p[m] = self.p[m-1]
            else:
                self.p[m,:,:] = self.p[m-1,:,:]
        elif (BC == 'NO_FLOW'):
            if (self.dz.size == 1):
                dz = self.dz
            else:
                dz = self.dz[m]
            #----------------------------------------------
            # Should we use self.dz[m] or self.dz[m-1] ??
            #----------------------------------------------
            if (self.SINGLE_PROFILE):
                self.p[m] = self.p[m-1] + dz
            else:
                self.p[m,:,:] = self.p[m-1,:,:] + dz
##        elif (BC == 'FIXED'):
##            if (self.SINGLE_PROFILE):
##                self.p[m] = value at time zero
##            else:
##                self.p[m,:,:] = value at time zero
                
        #-----------------------------------
        # Set theta at the bottom boundary
        #-----------------------------------
        if (self.SINGLE_PROFILE):
            psi     = self.p[m]       # [meters]
            theta_s = self.qs[m]
            theta_r = self.qr[m]
            psi_B   = self.pB[m]
            psi_A   = self.pA[m]
            c       = self.c[m]
            Lambda  = self.lam[m]
            self.q[m] = stbc.theta_of_psi(psi, theta_s, theta_r, \
                                  psi_B, psi_A, c, Lambda)
#             self.q[m] = Theta_TBC(psi, theta_s, theta_r, \
#                                   psi_B, psi_A, c, Lambda)
        else:
            #----------------------------------------
            # Now checking if ndim > 1 (2019-10-29)
            #----------------------------------------
            psi = self.p[m,:,:]   # [meters]
            #---------------------------------------
            if (self.qs.ndim > 1):
                theta_s = self.qs[m,:,:]
            else:
                theta_s = self.qs[m]
            #-----------------------------
            if (self.qr.ndim > 1):
                theta_r = self.qr[m,:,:]
            else:
                theta_r = self.qr[m]
            #-----------------------------
            if (self.pB.ndim > 1):
                psi_B = self.pB[m,:,:]
            else:
                psi_B = self.pB[m]
            #------------------------------
            if (self.pA.ndim > 1):
                psi_A = self.pA[m,:,:]
            else:
                psi_A = self.pA[m]
            #------------------------------
            if (self.c.ndim > 1):
                c = self.c[m,:,:]
            else:
                c = self.c[m]
            #------------------------------
            if (self.lam.ndim > 1):               
                Lambda = self.lam[m,:,:]
            else:
                Lambda = self.lam[m]
            #------------------------------            
            self.q[m,:,:] = stbc.theta_of_psi(psi, theta_s, theta_r, \
                                      psi_B, psi_A, c, Lambda)
#             self.q[m,:,:] = Theta_TBC(psi, theta_s, theta_r, \
#                                       psi_B, psi_A, c, Lambda)

        #----------------
        # For debugging
        #----------------
        ## if (self.SINGLE_PROFILE):
        if (self.DEBUG and self.SINGLE_PROFILE):
            print('In update_bottom_BC():')
            print('theta[m-1] =', self.q[m-1])           
            print('theta[m]   =', self.q[m])
            print('psi[m-1]   =', self.p[m-1])
            print('psi[m]     =', self.p[m] )
            print(' ')
            
    #   update_bottom_BC()
    #-----------------------------------------------------------------------
    def update_theta(self, REPORT=False):

        #----------------------------------------------------------
        # Notes:  This procedure updates the soil moisture, theta
        #         Theta is called "q" here .
        #----------------------------------------------------------
        # NB!  Forward derivative has last value "wrong",
        #      while backward derivative has first value "wrong".
        #      This is OK, as long as we explicitly set the
        #      first (surface) and last (bottom) values of
        #      theta using the boundary conditions.
        #----------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_theta()...')

        #-----------------------------------------------
        # Compute z[i+1] - z[i-1] (same for 1D and 3D)
        #-----------------------------------------------
        # First & last values of zdiff will be "wrong"
        #-----------------------------------------------
        z_below = np.roll( self.z, -1, axis=0 )
        z_above = np.roll( self.z,  1, axis=0 )
        z_diff  = (z_below - z_above)
        n_dz    = z_diff.size  # (should equal self.nz)
        #---------------------------------------------------
#         print('### min(z_diff) = ', z_diff.min() )
#         print('### max(z_diff) = ', z_diff.max() )
#         print('### n_dz        = ', n_dz )
#         print('### self.nz     = ', self.nz )
#         print('### z =')
#         print( self.z )
#         print('### z_diff =')
#         print( z_diff )

        #-----------------------------------------------
        # This should also work
        #------------------------
##        if (self.dz.size == 1):
##            z_diff = 2.0 * self.dz
##        else:
##           dz_below = np.roll( self.dz, -1, axis=0 )
##           dz_above = np.roll( self.dz,  1, axis=0 )
##           z_diff   = (dz_below + dz_above)
##        n_dz = z_diff.size

#         print('#### min(p)  = ', self.p.min() )
#         print('#### max(p)  = ', self.p.max() )
#         print('#### min(K)  = ', self.K.min() )
#         print('#### max(K)  = ', self.K.max() )
#         print('#### min(dz) = ', self.dz.min() )
#         print('#### max(dz) = ', self.dz.max() )
#         print('#### dz = ')
#         print( self.dz )
                               
        if (self.SINGLE_PROFILE): 
            #------------------------------------
            # Theta, psi, K and v are 1D arrays
            #------------------------------------
            # dz may be scalar or 1D array
            #------------------------------------
            dp_dz_1 = Z_Derivative_1D( self.p, self.dz )
            dp_dz_2 = Z_Derivative_1D( self.p, self.dz, BACKWARD=True )
            K_bar_1 = Z_Forward_Average( self.K )
            K_bar_2 = Z_Backward_Average( self.K )
            term1   = K_bar_2 * (dp_dz_2 - 1.0)
            term2   = K_bar_1 * (dp_dz_1 - 1.0)
            #------------------------------------------------------
            d_theta = (-2.0 * self.dt / z_diff) * (term1 - term2)

            #-----------------------------------------------
            # Surface and bottom values of theta were
            # already set in update_boundary_conditions().
            # dtheta has "wrong" values at those places.
            #-----------------------------------------------
            d_theta[0]         = 0.0 
            d_theta[self.nz-1] = 0.0

#             print('min(d_theta) =', d_theta.min() )
#             print('max(d_theta) =', d_theta.max() )
            
        else:    
            #------------------------------------
            # Theta, psi, K and v are 3D arrays
            #------------------------------------
            dp_dz_1 = Z_Derivative_3D( self.p, self.dz)
            dp_dz_2 = Z_Derivative_3D( self.p, self.dz, BACKWARD=True )
            K_bar_1 = Z_Forward_Average( self.K )
            K_bar_2 = Z_Backward_Average( self.K )
            term1   = K_bar_2 * (dp_dz_2 - 1.0)
            term2   = K_bar_1 * (dp_dz_1 - 1.0)
            #----------------------------------------------
            d_theta = (-2.0 * self.dt) * (term1 - term2)
            
            if (n_dz == 1):    
                d_theta = d_theta / z_diff   # (z_diff is a scalar)
            else:    
                for j in range(n_dz):
                    d_theta[j,:,:] = d_theta[j,:,:] / z_diff[j]
   
            #-----------------------------------------------
            # Surface and bottom values of theta were
            # already set in update_boundary_conditions().
            # dtheta has "wrong" values at those places.
            #-----------------------------------------------
            d_theta[0,:,:]         = 0.0 
            d_theta[self.nz-1,:,:] = 0.0
        
        #----------------------------------------------
        # Update soil moisture, theta  (mass balance)
        #----------------------------------------------
        self.q += d_theta

        #-----------------------------------------------
        # Option to check stability, in the sense that
        # theta values are still in range.
        #-----------------------------------------------
        # Recall that: S_eff  = (q - qr) / (qs - qr)
        # and that S_eff must be in [0, 1].
        #-----------------------------------------------
        if (self.CHECK_STABILITY):
            self.check_theta()

        if (self.DEBUG):
            self.check_theta()

    #   update_theta()
    #-----------------------------------------------------------------------
    def check_theta(self):

        ## bad1 = (self.q < self.qH)  ####### (was triggered, theta=0.217)
        ## bad1 = (self.q < self.qr)  #### (wrong shape)
        tol  = 0.01
        bad1 = (self.q < 0.01)
        bad2 = (self.q > self.qs + tol)
        bad3 = np.logical_not( np.isfinite( self.q ) )
        nb1  = bad1.sum()
        nb2  = bad2.sum()
        nb3  = bad3.sum()
        if (nb1 > 0) or (nb2 > 0) or (nb3 > 0):
            print('################################################')
            print('ERROR detected in update_theta():')
            print('  theta is the soil water content.')
            if (nb1 > 0):
                print('  theta < 0.01 (a residual value)')
                ## print('  theta < theta_H (hygroscopic lower limit)')
                qmin = self.q.min()
                print('  min(theta) = ' + str(qmin) )
#                 if (self.SINGLE_PROFILE):
#                     print('theta_H    = ' + str(self.qH) )
            if (nb2 > 0):
                print('  theta > theta_s (saturated upper limit)')
                qmax = self.q.max()
                print('  max(theta) = ' + str(qmax) )
#                 if (self.SINGLE_PROFILE):
#                     print('theta_s    = ' + str(self.qs) )
            print('Try reducing infil. timestep, dt, in CFG file.')
            print('################################################')
            print()
            #-----------------------------------------------
            # Next line will not stop the model run unless
            # this component is the driver component.
            #-----------------------------------------------
            self.DONE = True
            sys.exit()
        else:
            if (self.DEBUG):
                print('All theta values are in range.')
                print('min(theta) =', self.q.min() )
                print('max(theta) =', self.q.max() )
                print()

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
            #--------------------------------
            S_eff  = (self.q - self.qr) / (self.qs - self.qr)
            cpow   = (-self.c / self.lam)
            arg    = ((S_eff ** cpow )- 1.0) ** (1.0 / self.c)
            self.p = (self.pB * arg) - self.pA
            
            #------------------------------------------------------
            # S_eff = effective saturation or normalized vol frac
            # Could be shared, but requires more storage.
            # b = brooks-corey_b_parameter = 1 / lambda.
            #------------------------------------------------------
            ## self.S_eff = S_eff
            ## self.b = 1 / self.lam
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
                S_eff = (self.q[j,:,:] - qr) / (qs - qr)      #(grid)
                cpow  = (-c / lam)                            #(grid or scalar)
                arg   = ((S_eff ** cpow) - 1.0) ** (1.0 / c)  #(grid)
                self.p[j,:,:] = (pB * arg) - pA               #(grid)

                #------------------------------------------------------
                # S_eff = effective saturation or normalized vol frac
                # Could be shared, but requires more storage.
                # b = brooks-corey_b_parameter = 1 / lambda.
                #------------------------------------------------------
                ## self.S_eff = S_eff
                ## self.b = 1 / lam

        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            print('S_eff = ', S_eff[0:4])
            print('psi   = ', self.p[0:3])
            #print,' '

        if (self.DEBUG):
            print('In update_psi():')
            if (self.SINGLE_PROFILE):
                m = (self.nz - 1)
                print('psi[0], theta[0] =', self.p[0], ', ', self.q[0])
                print('psi[m], theta[m] =', self.p[m], ', ', self.q[m])
            #-----------------------------------
            print('min(q)   =', self.q.min() )
            print('max(q)   =', self.q.max() )
            print('min(qs)  =', self.qs.min() )
            print('max(qs)  =', self.qs.max() )
            print('min(qr)  =', self.qr.min() )
            print('max(qr)  =', self.qr.max() )
            print('min(Se)  =', S_eff.min()  )
            print('max(Se)  =', S_eff.max() )
            print('min(pB)  =', self.pB.min() )
            print('max(pB)  =', self.pB.max() )
            print('min(pA)  =', self.pA.min() )
            print('max(pA)  =', self.pA.max() )
            print('min(psi) =', self.p.min() )
            print('max(psi) =', self.p.max() )

                        
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
            epow = (-self.eta / self.c)
            Kr = (1.0 + ((self.p + self.pA) / self.pB) ** self.c) ** epow
            Kr = np.maximum( np.minimum(Kr, 1.0), 0.0)
            self.K = self.Ks * Kr
        else:    
            #--------------------------------------
            # Each var is either a 1D or 3D array
            #--------------------------------------
            dim_Ks  = np.ndim( self.Ks )
            dim_pB  = np.ndim( self.pB )
            dim_pA  = np.ndim( self.pA )
            dim_eta = np.ndim( self.eta )
            dim_c   = np.ndim( self.c )

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
                arg  = (self.p[j,:,:] + pA) / pB        #(grid)
                epow = (-eta / c)                       #(grid or scalar)
                Kr = (1.0 + arg ** c) ** epow           #(grid)
                Kr = np.maximum( np.minimum(Kr, 1.0), 0.0)   #(grid)
                self.K[j,:,:] = (Ks * Kr)               #(grid)

                #--------------------------------------------------                
                # Kr = (K/Ks) = relative hydraulic conductivity
                # This allows sharing, but requires more storage.
                #--------------------------------------------------
                ## self.Krel = Kr

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
        #        dp_dz is undefined for the last element, so
        #        we need to set v there.  Exactly how we do
        #        it doesn't really matter since v is not used
        #        directly in the computations.
        #-----------------------------------------------------------
            
        #----------------
        # For debugging
        #----------------
        if (self.DEBUG):
            print('Calling update_v()...')
            print('   SINGLE_PROFILE =', self.SINGLE_PROFILE)
        
        #-----------------------------------------------
        # NB! There are better ways to compute K_bar
        #     than the mean value used here.  Another
        #     method is discussed by R.E. Smith, at
        #     the top of page 192 in his AGU monograph.
        #-----------------------------------------------
        if (self.SINGLE_PROFILE):    
            #----------------------------------------
            # Theta, psi, K and v are all 1D arrays
            #----------------------------------------
            # dp_dz = (p_below - p) / dz
            #----------------------------------------         
            dp_dz  = Z_Derivative_1D( self.p, self.dz )
            K_bar  = Z_Forward_Average( self.K )
            self.v = K_bar * (1.0 - dp_dz)   # (bottom of cell)

            #-------------------
            # See Notes above.
            #-------------------
            m = self.nz - 1
            self.v[ m ] = self.v[ m-1 ]
        else:
            #----------------------------------------
            # Theta, psi, K and v are all 3D arrays
            #----------------------------------------
            # dp_dz = (p_below - p) / dz
            #----------------------------------------
            dp_dz  = Z_Derivative_3D( self.p, self.dz )
            K_bar  = Z_Forward_Average( self.K )
            self.v = K_bar * (1.0 - dp_dz)   # (bottom of cell)
            
            #-------------------
            # See Notes above.
            #-------------------
            m = self.nz - 1
            self.v[ m,:,: ] = self.v[ m-1,:,: ]
        
        #----------------
        # For debugging
        #----------------
##        print 'At n_filling test...'
##        print 'shape(qs) =', np.shape(self.qs)
##        print 'shape(q)  =', np.shape(self.q)
##        print 'shape(vT) =', np.shape(vT)
##        print 'shape(vB) =', np.shape(vB)
##        print 'shape(v)  =', np.shape(self.v)

        #-----------------------------------
        # Return flow rate in bottom layer
        #-----------------------------------
        if (self.SINGLE_PROFILE):    
            self.Rg = self.v[self.nz - 1]
        else:    
            self.Rg = self.v[self.nz - 1,:,:]

        if (self.DEBUG):
            print('min(v), max(v) =', self.v.min(), self.v.max())
        
    #   update_v()
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
            diff[self.nz - 1] = np.float64(0)
            indices = np.where(diff > 0)        # (must be > not >=)
            nd = indices[0].size    
            
            if (nd == 0):    
                #----------------------------------
                # Wetting front is at the surface
                #----------------------------------------------
                # This can happen for an equilibrium profile
                # that is monotonically increasing with depth
                # or for case where theta = theta_i for all z
                #----------------------------------------------
                self.Zw = np.float64(0)
            else:    
                imax = indices[0][nd - 1]  ########################
                
                #----------------------------------------------
                # Get min and max theta of decreasing section
                #----------------------------------------------
                frac = np.float64(0.2)
                tmax = np.nanmax( self.q[indices] )
                tmin = np.nanmin( self.q[indices] )
                tmid = tmin + frac * (tmax - tmin)
                w    = np.where( self.q[0: (imax+1)] > tmid )
                nw2  = np.size( w[0] )
                if (nw2 > 0):    
                    imax2 = w[0][nw2 - 1]   ###################
                    self.Zw = imax2 * self.dz
                else:    
                    self.Zw = np.float64(0)         #####  IS THIS RIGHT ?
        else:    
            
            q_below = np.roll( self.q, -1, axis=0 )
            diff = (self.q - q_below)
            diff[self.nz - np.int32(1),:,:] = np.float64(0)
            n_dz = self.dz.size     #(either 1 or self.nz)
            
            #--------------------------------
            # Zero out the 2D array self.Zw
            #------------------------------------------
            # Note: *self.Zw should have already been
            # set to be a 2D array of correct size
            #------------------------------------------
            self.Zw = np.minimum( np.maximum(self.Zw, np.float64(0)), np.float64(0) )
            
            #---------------------------------------
            # Loop over the z-levels to find Z_wet
            #---------------------------------------
            for j in range(self.nz - 1):           #(nz-1) vs. nz
                diff_j = diff[j,:,:]
                next_diff_j = diff[j + 1,:,:]     #*******
                if (n_dz == 1):    
                    dz = self.dz
                else:    
                    dz = self.dz[j]
                
                #--------------------------------------------------------
                # NB!  We are looking for a local min in theta profile.
                #--------------------------------------------------------
                # NB!  If theta is same at all levels, then we will
                # never get (diff_j GT 0) so Z_wet will remain at 0,
                # even if (theta eq theta_s) at all levels !!!!    **********
                # How can we get Z_wet = Z_bot = (nz-1) * dz ???   **********
                #--------------------------------------------------------
                # NB!  If theta increases at all levels, which can
                # happen for an equilibrium profile with water table,
                # then we will never get (diff_j GT 0.0) and Z_wet will
                # remain at 0.
                #--------------------------------------------------------
                # NB!  For (j eq (nz-2)), we have (next_diff_j EQ 0.0).
                #--------------------------------------------------------
                IDs = np.where( np.logical_and((diff_j > 0), \
                                               (next_diff_j <= 0)))
                n_IDs = IDs[0].size
                
                if (n_IDs != 0):    
                    self.Zw[IDs] = (j * dz)

    #   update_Zw()
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
        
        #-------------------------------------------------------------
        # Richards' equation is only used in the so-called "upper
        # layers".  There can be between 1 and 3 of these layers.
        # To avoid high computational cost in the less dynamic
        # "lower zone" below, a simplified method is used to route
        # flow through the lower zone to the water table and to
        # estimate a "groundwater recharge rate", Rg.  This is done
        # using the vertical flow rate at the bottom of the set of
        # upper layers and perhaps other information.  If the water
        # table intersects the upper layers, then Rg is computed as
        # the vertical flow rate in the layer just above the water
        # table.
        #-------------------------------------------------------------
        # The simplest method of estimating Rg is to simply set it
        # to the vertical flow rate at the bottom of the upper zone.
        # However, travel time through the lower zone is not taken
        # into account.  If Rg < 0, then water can even be drawn
        # upward from the water table.
        #-------------------------------------------------------------
        # Another simple idea for estimating Rg is to assume that
        # psi varies linearly in the lower zone:   psi = a*z + b.
        # Since psi=0 at the bottom and psi=psi_zb at the top, we
        # can compute a and b as:
        #     a = -psi_zb/(z_wt - z_zb)
        #     b = -a * z_wt
        # The flow rate in the lower zone is then computed as:
        #     v = K(a*z + b) * (1 - a)
        # It follows that the flow rate at the bottom of the lower
        # zone can be written as:
        #     Rg = v(z_wt) = v(z_zb) * (Ks / K(psi_zb)).
        # Since K(psi_zb) <= Ks, we have v(z_zb) < Rg < Ks.
        # If psi_zb < (z_zb - z_wt), then we get Rg < 0 and water
        # can be drawn upward from the water table.
        #-------------------------------------------------------------
        
        #----------------------------------------
        # For testing:  Plot the theta profiles
        #----------------------------------------
        PLOT = False
        if (PLOT):    
            matplotlib.pyplot.figure(1)
            #** wait, 0.005
            
            ymin = self.qi.min()
            ymax = (self.qs + np.float64(0.05)).max()
            if (self.SINGLE_PROFILE):    
                matplotlib.pyplot.plot(self.z, self.q, marker='+')
                matplotlib.pyplot.xlabel('Depth [meters]')
                matplotlib.pyplot.ylim(np.array(ymin, ymax))
                matplotlib.pyplot.axis('image')
                matplotlib.pyplot.ylabel('Soil moisture')
                matplotlib.pyplot.show()
            else:    
                matplotlib.pyplot.plot(self.z, self.q[:,2,2], marker='+')
                matplotlib.pyplot.xlabel('Depth [meters]')
                matplotlib.pyplot.ylim(np.array(ymin, ymax))
                matplotlib.pyplot.axis('image')
                matplotlib.pyplot.ylabel('Soil moisture')
                matplotlib.pyplot.show()
        
        #--------------------------------------
        # For testing:  Plot the psi profiles
        #--------------------------------------
        if (PLOT):    
            matplotlib.pyplot.figure(2)
            ### wait, 0.005
            yrange = np.array([-np.float32(3.0), np.float32(0.5)])    #(Log case)
            ytitle = '-Log(-Pressure head) [m]'
            #--------------------------------------
            # ytitle = 'Pressure head [m]'
            # yrange = [-20.0, 0.0]  ;(Linear case)
            #--------------------------------------
            if (self.SINGLE_PROFILE):    
                y = -np.float64(1) * np.log(np.absolute(self.p) + 1)
                matplotlib.pyplot.plot(self.z, y, marker='+')
                matplotlib.pyplot.xlabel('Depth [meters]')
                matplotlib.pyplot.ylim(yrange)
                matplotlib.pyplot.axis('image')
                matplotlib.pyplot.ylabel(ytitle)
                matplotlib.pyplot.show()
            else:    
                y = -np.float64(1) * np.log(np.absolute((self.p)[:,2,2]) + 1)
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
            print('Calling update_q0()...')
            
        if (self.ALL_SCALARS): 
            self.q0 = self.q[0]
        else:    
            self.q0 = self.q[0,:,:]
            
    #   update_q0()
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
        
        for k in range(self.n_layers):
            self.Ks_file[k]  = self.in_directory + self.Ks_file[k]
            self.Ki_file[k]  = self.in_directory + self.Ki_file[k]
            self.qs_file[k]  = self.in_directory + self.qs_file[k]
            self.qi_file[k]  = self.in_directory + self.qi_file[k]
            self.qr_file[k]  = self.in_directory + self.qr_file[k]
            self.pB_file[k]  = self.in_directory + self.pB_file[k]
            self.pA_file[k]  = self.in_directory + self.pA_file[k]
            self.lam_file[k] = self.in_directory + self.lam_file[k]
            self.c_file[k]   = self.in_directory + self.c_file[k]

            self.Ks_unit.append(  model_input.open_file(self.Ks_type[k],  self.Ks_file[k]) )
            self.Ki_unit.append(  model_input.open_file(self.Ki_type[k],  self.Ki_file[k]) )
            self.qs_unit.append(  model_input.open_file(self.qs_type[k],  self.qs_file[k]) )
            self.qi_unit.append(  model_input.open_file(self.qi_type[k],  self.qi_file[k]) )
            self.qr_unit.append(  model_input.open_file(self.qr_type[k],  self.qr_file[k]) )
            self.pB_unit.append(  model_input.open_file(self.pB_type[k],  self.pB_file[k]) )
            self.pA_unit.append(  model_input.open_file(self.pA_type[k],  self.pA_file[k]) )
            self.lam_unit.append( model_input.open_file(self.lam_type[k], self.lam_file[k]) )
            self.c_unit.append(   model_input.open_file(self.c_type[k],   self.c_file[k]) )

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        #-----------------------------------------------------
        # Notes:  Override infil_base's method by same name.
        #-----------------------------------------------------
        if (self.DEBUG):
            print('Calling read_input_files()...')

        rti = self.rti

        for j in range(self.n_layers):        
            Ks_val = model_input.read_next(self.Ks_unit[j], self.Ks_type[j], rti)
            if (Ks_val is not None): self.Ks_list[j] = Ks_val

            Ki_val = model_input.read_next(self.Ki_unit[j], self.Ki_type[j], rti)
            if (Ki_val is not None): self.Ki_list[j]  = Ki_val

            qs_val = model_input.read_next(self.qs_unit[j], self.qs_type[j], rti)
            if (qs_val is not None): self.qs_list[j]  = qs_val

            qi_val = model_input.read_next(self.qi_unit[j], self.qi_type[j], rti)
            if (qi_val is not None): self.qi_list[j]  = qi_val
            
            qr_val = model_input.read_next(self.qr_unit[j], self.qr_type[j], rti)
            if (qr_val is not None): self.qr_list[j]  = qr_val

            pB_val = model_input.read_next(self.pB_unit[j], self.pB_type[j], rti)
            if (pB_val is not None): self.pB_list[j]  = pB_val

            pA_val = model_input.read_next(self.pA_unit[j], self.pA_type[j], rti)
            if (pA_val is not None): self.pA_list[j]  = pA_val

            lam_val = model_input.read_next(self.lam_unit[j], self.lam_type[j], rti)
            if (lam_val is not None):
                #---------------------------------------------------------
                # If we read a lambda value from a file, then we need to
                # compute and save corresponding eta = [2 + (3*lambda)]
                # BUG FIX:  2019-10-29
                #---------------------------------------------------------
                self.lam_list[j] = lam_val
                self.eta_list[j] = np.float64(2) + (np.float64(3) * lam_val)
                
            c_val = model_input.read_next(self.c_unit[j], self.c_type[j], rti)
            if (c_val is not None): self.c_list[j]  = c_val
 
            #----------------------------------------------
            # Update qH, given by theta_of_psi() function
            #----------------------------------------------
            self.qH_list[j] = stbc.theta_of_psi( self.psi_hygro, \
                                         self.qs_list[j], self.qr_list[j], \
                                         self.pB_list[j], self.pA_list[j], \
                                         self.c_list[j],  self.lam_list[j] )
#             self.qH_list[j] = Theta_TBC( self.psi_hygro, \
#                                          self.qs_list[j], self.qr_list[j], \
#                                          self.pB_list[j], self.pA_list[j], \
#                                          self.c_list[j],  self.lam_list[j] )

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

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def Theta_TBC(psi, theta_s, theta_r, psi_B, psi_A, c, Lambda,
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
        psi_m = (psi/ np.float64(100))        # [cm -> meters]
        ratio = (psi_m + psi_A) / psi_B    # (should be > 0)
    else:
        ratio = (psi + psi_A) / psi_B      # (should be > 0)
        
    theta = (1.0 + (ratio ** c)) ** (-Lambda / c)
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
    
#   K_of_Theta_TBC()
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def Z_Derivative_1D(v, dz, BACKWARD=False):

    #----------------------------------------------------------
    # Notes:  v is a 1D array and dz is a scalar or 1D array.
    #         The result is a 1D array, same size as v.

    #        This function does not worry about the wrap
    #        around affect of ROLL at bottom.  This must
    #        be handled by the caller.

    #        (11/11/10) Added BACKWARD keyword.
    #----------------------------------------------------------
    if not(BACKWARD):
        v_below = np.roll(v, -1, axis=0) 
        dv_dz   = (v_below - v) / dz
        ## dv_dz[self.nz - 1] = ????
    else:
        v_above = np.roll(v, 1, axis=0)
        dv_dz   = (v - v_above) / dz
        ## dv_dz[0] = ????
        
    return dv_dz
    
#   Z_Derivative_1D()
#-----------------------------------------------------------------------
def Z_Derivative_3D(v, dz, BACKWARD=False):

    #------------------------------------------------------------
    # Note:  v is a 3D array (or data cube) and dz is a scalar
    #        or 1D array.  The result is a 3D array, same size
    #        as v.

    #        This function does not worry about the wrap
    #        around affect of ROLL at bottom.  This must
    #        be handled by the caller.

    #        (11/11/10) Added BACKWARD keyword.
    #-------------------------------------------------------------
    # Note:  If dz_val[j] is the same for all layers (j), so
    #        that dz_val.min() == dz_val.max(), then self.dz
    #        is set to a scalar value.  self.dz and self.nz are
    #        both set in infil_base.py (build_layer_z_vector()).
    #-------------------------------------------------------------    
    n_dz = dz.size

    if not(BACKWARD):
        v_below = np.roll(v, -1, axis=0)
        
        if (n_dz == 1):    
            dv_dz = (v_below - v) / dz  # (dz is a scalar)
        else:    
            dv_dz = (v_below - v)
            for j in range(n_dz):
                dv_dz[j,:,:] = dv_dz[j,:,:] / dz[j]
        ## dv_dz[self.nz - 1] = ????
    else:
        v_above = np.roll(v, 1, axis=0)

        # v_above has same min & max
#         print('min(v) = ' + str(v.min()) )
#         print('max(v) = ' + str(v.max()) )
        
        if (n_dz == 1):
            dv_dz = (v - v_above) / dz  # (dz is a scalar)
        else:    
            dv_dz = (v - v_above)
            for j in range(n_dz):
                dv_dz[j,:,:] = dv_dz[j,:,:] / dz[j]
        ## dv_dz[0] = ????
        
    return dv_dz
    
#   Z_Derivative_3D()
#-----------------------------------------------------------------------
def Z_Forward_Average( v ):

    #----------------------------------------------
    # Notes: This should work for both 1D and 3D.
    #        For 3D, "axis=0" is the z-axis.
    #----------------------------------------------
    v_below = np.roll(v, -1, axis=0)
    v_avg   = (v_below + v) / 2.0
    ## v_avg[self.nz - 1] = ?????
    
    return v_avg

#   Z_Forward_Average()
#-----------------------------------------------------------------------
def Z_Backward_Average( v ):

    #----------------------------------------------
    # Notes: This should work for both 1D and 3D.
    #        For 3D, "axis=0" is the z-axis.
    #----------------------------------------------
    v_above = np.roll(v, 1, axis=0)
    v_avg   = (v + v_above) / 2.0
    ## v_avg[0] = ??????
    
    return v_avg

#   Z_Backward_Average()
#-----------------------------------------------------------------------





