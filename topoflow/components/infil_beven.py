"""
This class defines a hydrologic infiltration component that implements the model
introduced by Beven (1984) where saturated hydraulic conductivity, Ks, decreases
exponentially with subsurface depth, z.

This class inherits from the infiltration "base class" in "infil_base.py".

Note:  THIS HAS NOT BEEN TESTED YET.  STILL NEED TO CREATE A CFG FILE.

See: Beven, K. (1984) "Infiltration into a class of vertically non-uniform soils",
Hydrol. Sciences J., 29(4), 425-434.
"""
## Copyright (c) 2009-2016, Scott D. Peckham
## May 2009
## May 2010 (changes to unit_test())
## Nov 2016.  Completed.

#-----------------------------------------------------------------------
#
#  class infil_component         # (inherits from infil_base.py)
#
#      get_component_name()
#      get_attribute()           # (10/26/11)
#      get_input_var_names()     # (5/15/12)
#      get_output_var_names()    # (5/15/12)
#      get_var_name()            # (5/15/12)
#      get_var_units()           # (5/15/12)
#      -------------------------
#      initialize_layer_vars()     # (5/11/10)
#      set_computed_input_vars()
#      check_input_types()
#      update_infil_rate()
#      -------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      ----------------------
#      update_outfile_names()

#  Functions:
#      Beven_Exp_K_Infil_Rate_v1  (Nov. 2016)
#      Beven_Exp_K_Infil_Rate_1D  (Not written)
#      Beven_Exp_K_Infil_Rate_3D  (Not written)
#      Beven_Exp_K_Infil_Rate_v2  (Not written)

#-----------------------------------------------------------------------

import numpy as np

# import os

from topoflow.components import infil_base

from topoflow.utils import model_input
from topoflow.utils.tf_utils import TF_Print, TF_String

#-----------------------------------------------------------------------
class infil_component( infil_base.infil_component ):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Infiltration_Beven',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'InfilBeven',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Infil_Beven.cfg.in',
        'cfg_extension':      '_infil_beven.cfg',
        'cmt_var_prefix':     '/InfilBeven/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Infil_Beven.xml',
        'dialog_title':       'Infiltration: Beven Parameters',
        'time_units':         'seconds' }

    _input_var_names = [
        'atmosphere_water__rainfall_volume_flux',           # (P_rain)  
        'glacier_ice__melt_volume_flux',                    # (MR)
        'land_surface__elevation',                          # (elev)
        'land_surface_water__evaporation_volume_flux',      # (ET)
        'snowpack__melt_volume_flux',                       # (SM)
        'soil_water_sat-zone_top_surface__elevation' ]      # (h_table)

        ##  'land_surface_water__baseflow_volume_flux',     # (GW)

    #-----------------------------------------------------------------
    # These are input vars provided by the user, but can also
    # be retrieved by other components as "output_vars".
    #-----------------------------------------------------------------
    #    'soil_water__saturated_hydraulic_conductivity', # Ks
    #    'soil_water__saturated_volume_fraction',        # qs
    #    'soil_water__initial_volume_fraction',          # qi
    #
    #-----------------------------------------------------------------
    _output_var_names = [
        'model__time_step',                   # dt
        # 'model_grid_cell__area',           # da
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux', # vol_IN
        'soil_surface_water__infiltration_volume_flux',                       # IN
        'soil_surface_water__time_integral_of_infiltration_volume_flux',      # I
        'soil_water__green-ampt_capillary_length',         # G
        'soil_water__potential_infiltration_volume_flux',  # fc
        'soil_water__initial_hydraulic_conductivity',      # Ki
        'soil_water__initial_volume_fraction',             # qi
        'soil_water__saturated_hydraulic_conductivity',    # Ks
        'soil_water__saturated_volume_fraction',           # qs
        'soil_water_flow__z_component_of_darcy_velocity',  # v
        'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux',  # vol_Rg
        'soil_water_sat-zone_top__recharge_volume_flux' ]  # Rg

    _var_name_map = {
        'atmosphere_water__rainfall_volume_flux':          'P_rain',   
        'glacier_ice__melt_volume_flux':                   'MR',
        'land_surface__elevation':                         'elev',
        'land_surface_water__evaporation_volume_flux':     'ET',
        'snowpack__melt_volume_flux':                      'SM',
        'soil_water_sat-zone_top_surface__elevation':      'h_table',
        #--------------------------------------------------------------
        'model__time_step': 'dt',
        # 'model_grid_cell__area': 'da',
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux': 'vol_IN',
        'soil_surface_water__infiltration_volume_flux':    'IN',
        'soil_surface_water__time_integral_of_infiltration_volume_flux': 'I',
        'soil_water__green-ampt_capillary_length':         'G',
        'soil_water__initial_hydraulic_conductivity':      'Ki',
        'soil_water__initial_volume_fraction':             'qi',
        'soil_water__potential_infiltration_volume_flux':  'fc',
        'soil_water__saturated_hydraulic_conductivity':    'Ks',
        'soil_water__saturated_volume_fraction':           'qs',
        'soil_water_flow__z_component_of_darcy_velocity':  'v',
        'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux': 'vol_Rg',
        'soil_water_sat-zone_top__recharge_volume_flux':    'Rg' }

    _var_units_map = {
        'atmosphere_water__rainfall_volume_flux':          'm s-1',   
        'glacier_ice__melt_volume_flux':                   'm s-1',
        'land_surface__elevation':                         'm',
        'land_surface_water__evaporation_volume_flux':     'm s-1',
        'snowpack__melt_volume_flux':                      'm s-1',
        'soil_water_sat-zone_top_surface__elevation':      'm',
        #-----------------------------------------------------------
        'model__time_step': 's',
        ## 'model_grid_cell__area': 'm2', 
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux': 'm3',
        'soil_surface_water__infiltration_volume_flux': 'm s-1',
        'soil_surface_water__time_integral_of_infiltration_volume_flux': 'm',
        'soil_water__green-ampt_capillary_length': 'm',
        'soil_water__initial_hydraulic_conductivity': 'm s-1',
        'soil_water__initial_volume_fraction': '1',
        'soil_water__potential_infiltration_volume_flux': 'm s-1',
        'soil_water__saturated_hydraulic_conductivity': 'm s-1',
        'soil_water__saturated_volume_fraction': '1',
        'soil_water_flow__z_component_of_darcy_velocity': 'm s-1',
        'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux': 'm3',
        'soil_water_sat-zone_top__recharge_volume_flux': 'm s-1' }

    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Infiltration_Beven'

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
    def initialize_layer_vars(self):

        #----------------------------------------------------------
        # Notes: initialize_config_vars() calls read_cfg_file().
        #        If read_cfg_file() finds a variable "n_layers",
        #        then it calls initialize_layer_vars() so that
        #        subsequent layer variables - as indicated by a
        #        subscript in the CFG file - can be read directly
        #        into a list or array.
        #----------------------------------------------------------
        n_layers = self.n_layers
        
        #-------------------------------------------------
        # Get arrays to store soil params for each layer
        #-------------------------------------------------
        self.soil_type = np.zeros([n_layers], dtype='<U200')
##        self.dz_val    = np.zeros([n_layers], dtype='Float64')    #### + dz3
##        self.nz_val    = np.zeros([n_layers], dtype='Int16')      #### + nz3
        #----------------------------------------------------------
        self.Ks_type  = np.zeros(n_layers, dtype='<U200')
        self.qs_type  = np.zeros(n_layers, dtype='<U200')
        self.qi_type  = np.zeros(n_layers, dtype='<U200')
        self.f_type   = np.zeros(n_layers, dtype='<U200')
        self.C_type   = np.zeros(n_layers, dtype='<U200')
        #--------------------------------------------------------        
        self.Ks_file  = np.zeros(n_layers, dtype='<U200')
        self.qs_file  = np.zeros(n_layers, dtype='<U200')
        self.qi_file  = np.zeros(n_layers, dtype='<U200')
        self.f_file   = np.zeros(n_layers, dtype='<U200')
        self.C_file   = np.zeros(n_layers, dtype='<U200')
        #---------------------------------------------------------
        # Note: self.Ks is a Python list.  Initially, each entry
        # is a numpy scalar (type 'np.float64').  However, we
        # can later change any list entry to a scalar or grid
        # (type 'np.ndarray'), according to its "Ks_type".
        #---------------------------------------------------------
        # While CFG file for Richards 1D uses "Ks_val[0]", the
        # CFG file for Green-Ampt, etc. just uses "Ks[0]".
        # Too late to change it to be consistent.
        #---------------------------------------------------------   
        # NOTE!  In "initialize_computed_vars()", these lists
        # will be used to build ndarrays *with the same names*.
        #---------------------------------------------------------   
        self.Ks  = list( np.zeros(n_layers, dtype='Float64') )
        self.qs  = list( np.zeros(n_layers, dtype='Float64') )
        self.qi  = list( np.zeros(n_layers, dtype='Float64') )
        self.f   = list( np.zeros(n_layers, dtype='Float64') )
        self.C   = list( np.zeros(n_layers, dtype='Float64') )
        #--------------------------------------------------------- 
#         self.Ks_list  = list( np.zeros(n_layers, dtype='Float64') )
#         self.qs_list  = list( np.zeros(n_layers, dtype='Float64') )
#         self.qi_list  = list( np.zeros(n_layers, dtype='Float64') )
#         self.f_list   = list( np.zeros(n_layers, dtype='Float64') )
#         self.C_list   = list( np.zeros(n_layers, dtype='Float64') )
        
    #   initialize_layer_vars()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        #--------------------------------------------------------
        # Define these here, so all components can use the same
        # output file functions, like "open_output_files()".
        #--------------------------------------------------------
        self.RICHARDS       = False
        self.SAVE_Q0_GRIDS  = False
        self.SAVE_ZW_GRIDS  = False
        self.SAVE_Q0_PIXELS = False
        self.SAVE_ZW_PIXELS = False

        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt    = np.maximum(self.save_grid_dt,    self.dt)
        self.save_pixels_dt  = np.maximum(self.save_pixels_dt,  self.dt)
        
    #   set_computed_input_vars()
    #-------------------------------------------------------------------
    def check_input_types(self):

        #------------------------------------------------------
        # Notes: Usually this will be overridden by a given
        #        method of computing ET.  But this one should
        #        work for Green-Ampt and Smith-Parlange.
        #------------------------------------------------------
        are_scalars = np.array([
                         self.is_scalar('P_rain'),
                         self.is_scalar('SM'),
                         self.is_scalar('h_table'),
                         #----------------------------
                         self.is_scalar('Ks[0]'),
                         self.is_scalar('qs[0]'),
                         self.is_scalar('qi[0]'),
                         self.is_scalar('f[0]'),
                         self.is_scalar('C[0]')  ])

        self.ALL_SCALARS = np.all(are_scalars)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):
  
        #-----------------------------------------------------------
        # Note: h  = water table elevation [m]
        #       z  = land surface elevation [m]
        #
        #       Currently, h and z are always grids,
        #       so IN will be a grid.
        #       z, h and IN must be compatible.
        #-----------------------------------------------------------
        self.vol_IN = self.initialize_scalar( 0, dtype='float64')
        self.vol_Rg = self.initialize_scalar( 0, dtype='float64')
        
        if (self.ALL_SCALARS):
            #-----------------------------------------------------
            # Note: "I" is initialized to 1e-6 to avoid a divide
            #       by zero when first computing fc, which does
            #       have a singularity at the origin.
            #-----------------------------------------------------
            self.IN     = self.initialize_scalar( 0,    dtype='float64')
            self.Rg     = self.initialize_scalar( 0,    dtype='float64') 
            self.I      = self.initialize_scalar( 1e-6, dtype='float64')
            self.tp     = self.initialize_scalar( -1,   dtype='float64')
            self.fp     = self.initialize_scalar( 0,    dtype='float64')
            # self.r_last = self.initialize_scalar( 0,  dtype='float64') # (P+SM at prev step)
        else:
            self.IN     = self.initialize_grid( 0,    dtype='float64')
            self.Rg     = self.initialize_grid( 0,    dtype='float64')
            self.I      = self.initialize_grid( 1e-6, dtype='float64')
            self.tp     = self.initialize_grid( -1,   dtype='float64')
            self.fp     = self.initialize_grid( 0,    dtype='float64')
            # self.r_last = self.initialize_grid( 0,   dtype='float64')

        #-----------------------------------------------------
        # Compute dz as 1D array from scalars in self.dz_val      
        #-----------------------------------------------------
        # Compute the z-vector, for plotting profiles
        #----------------------------------------------
        self.build_layer_z_vector()

        #------------------------------------------------
        # Now build a 1D or 3D array for each input var
        #----------------------------------------------------------
        # (3/12/08) Same code should work if (self.n_layers eq 1)
        #----------------------------------------------------------
        # Convert from lists to arrays; same name. (11/15/16)
        #----------------------------------------------------------
        self.Ks  = self.build_layered_var( self.Ks )
        self.qs  = self.build_layered_var( self.qs )
        self.qi  = self.build_layered_var (self.qi )
        self.f   = self.build_layered_var( self.f)
        self.C   = self.build_layered_var( self.C)
     
    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def update_infil_rate(self):

        Beven_Exp_K_Infil_Rate_v1(self)
                
    #   update_infil_rate()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        #-------------------------------------------------------
        # NB! Beven, Green-Ampt & Smith-Parlange currently only
        #     support ONE layer (n_layers == 1).
        #-------------------------------------------------------
        self.Ks_unit  = []  # (empty lists to hold file objects)
        self.qs_unit  = []
        self.qi_unit  = []
        self.f_unit   = []
        self.C_unit   = []

        for k in range(self.n_layers):
            self.Ks_file[k] = self.in_directory + self.Ks_file[k]
            self.qs_file[k] = self.in_directory + self.qs_file[k]
            self.qi_file[k] = self.in_directory + self.qi_file[k]
            self.f_file[k]  = self.in_directory + self.f_file[k]
            self.C_file[k]  = self.in_directory + self.C_file[k]

            self.Ks_unit.append( model_input.open_file(self.Ks_type[k], self.Ks_file[k]) )
            self.qs_unit.append( model_input.open_file(self.qs_type[k], self.qs_file[k]) )
            self.qi_unit.append( model_input.open_file(self.qi_type[k], self.qi_file[k]) )
            self.f_unit.append(  model_input.open_file(self.f_type[k],  self.f_file[k])  )
            self.C_unit.append(  model_input.open_file(self.C_type[k],  self.C_file[k]) )

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        rti = self.rti

        #-------------------------------------------------------
        # All grids are assumed to have data type of Float32.
        #-------------------------------------------------------
        # NB! Beven, Green-Ampt & Smith-Parlange currently only
        #     support ONE layer (n_layers == 1).
        #------------------------------------------------------- 
        for k in range(self.n_layers):
            Ks = model_input.read_next(self.Ks_unit[k], self.Ks_type[k], rti)
            if (Ks is not None): self.Ks[k] = Ks

            qs = model_input.read_next(self.qs_unit[k], self.qs_type[k], rti)
            if (qs is not None): self.qs[k] = qs

            qi = model_input.read_next(self.qi_unit[k], self.qi_type[k], rti)
            if (qi is not None): self.qi[k] = qi
   
            f  = model_input.read_next(self.f_unit[k], self.f_type[k], rti)
            if (f is not None): self.f[k] = f

            C = model_input.read_next(self.C_unit[k], self.C_type[k], rti)
            if (C is not None): self.C[k] = C
           
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):

        #--------------------------------------------------------
        # NB! Beven, Green-Ampt & Smith-Parlange currently only
        #     support ONE layer (n_layers == 1).
        #--------------------------------------------------------        
        for k in range(self.n_layers):
            if (self.Ks_type[k]  != 'Scalar'): self.Ks_unit[k].close()        
            if (self.qs_type[k]  != 'Scalar'): self.qs_unit[k].close()
            if (self.qi_type[k]  != 'Scalar'): self.qi_unit[k].close()
            if (self.f_type[k]   != 'Scalar'): self.f_unit[k].close()
            if (self.C_type[k]   != 'Scalar'): self.C_unit[k].close()
          
    #   close_input_files()
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.v0_gs_file = (self.out_directory + self.v0_gs_file)
        self.I_gs_file  = (self.out_directory + self.I_gs_file)
        #---------------------------------------------------------
        self.v0_ts_file = (self.out_directory + self.v0_ts_file)
        self.I_ts_file  = (self.out_directory + self.I_ts_file)

##        self.v0_gs_file = (self.case_prefix + '_2D-v0.rts')
##        self.I_gs_file  = (self.case_prefix + '_2D-I.rts')
##        self.q0_gs_file = (self.case_prefix + '_2D-q0.rts')
##        self.Zw_gs_file = (self.case_prefix + '_2D-Zw.rts')
##        #---------------------------------------------------------
##        self.v0_ts_file = (self.case_prefix + '_0D-v0.txt')
##        self.I_ts_file  = (self.case_prefix + '_0D-I.txt')
##        self.q0_ts_file = (self.case_prefix + '_0D-q0.txt')
##        self.Zw_ts_file = (self.case_prefix + '_0D-Zw.txt')

    #   update_outfile_names() 
    #-------------------------------------------------------------------

#-----------------------------------------------------------------------
def Beven_Exp_K_Infil_Rate_v1(self):

    #------------------------------------------------------------
    #Notes:  This comes from: Beven, K. (1984) "Infiltration
    #        into a class of vertically non-uniform soils",
    #        Hydrol. Sciences J., 29(4), 425-434.
    #
    #        IN  = infiltration rate [m/s]
    #        Ks  = saturated hydraulic conductivity [m/s]
    #            = Ks0 * exp(f * z)   # Beven, f < 0.
    #        Ki  = initial hydraulic conductivity [m/s]
    #              (not applicable in this method)
    #        qs  = soil moisture content [dimless]
    #        qi  = soil initial moisture content [dimless]
    #        C   = "storage suction factor" [m]
    #            = (delta_psi * delta_theta) ~ constant > 0
    #        f   = parameter in K*(z) = K0 * exp(f*z) [1/m]
    #              K* < Ks, is "eff. K behind wetting front"
    #              (f < 0, roughly between -1 and -13)
    #              Is f roughly equal to (-gamma/G) in the
    #              Smith-Parlange method?  (Compare denoms.)
    #              f is sometimes written as (del_theta/m).
    #        K*  = effective conductivity behind wetting front
    #            = K0 * exp( f * z ) < Ks
    #        K0  = K*(0) = coeff. in previous equation [m/s]
    #              (Here, we assume that K0 = Ks[0])
    #         I  = cum. infiltration depth (since reset) [m]
    #         P  = precipitation rate [m/s]
    #        SM  = snowmelt rate [m/s]
    #         r  = (P + SM)  [m/s]

    #        Note that the infiltration rate has a max possible
    #        value of (P + SM).  It seems that this method does
    #        not asymptote to Ks as I increases.  Can we simply
    #        add Ks?  We should reset I between "events", but
    #        haven't decided how to do this yet.

    #        Total infiltrated depth, I, is incremented in the
    #        infiltration base class.
    #------------------------------------------------------------

    #---------------------------------------------
    # Note that Ks is used to store K* and C is
    # different than c used for Richards' method.
    # Need to add f to the set of infil vars.
    # Also need to add GUI to collect all vars.
    # Note:  t1<0, t3<0 & t1/t3 > 0; C>0 & t2>0.
    #---------------------------------------------
    dq = (self.qs[0] - self.qi[0])
    K0 = self.Ks[0]   ##### ????   ## K0 is a model parameter < Ks0.
    t1 = (K0 * self.f) / dq        ## f is a model parameter, in [-13, -1]?
    t2 = (self.C + self.I)         ## C is a model parameter in [0, 0.1] m (Beven)
    t3 = (np.float64(1) - np.exp(-np.float64(1) * self.f * self.I / dq))
    fc = t1 * t2 / t3
    
    #-----------------------------------------------
    # Initially, IN = r, and all of the incoming
    # water infiltrates.  IN cannot exceed r.
    # Ponding time, Tp, is time until (IN lt r).
    #-----------------------------------------------
    self.fc = fc
    self.IN = np.minimum(fc, self.P_total)
    
    #-------------------------------------
    # Is P_total less than Ks anywhere ?
    # If so, set IN = P_total there.
    #-------------------------------------
    # self.check_low_rainrate()  # (now done in infil_base.py)
       
#   Beven_Exp_K_Infil_Rate_v1
#-----------------------------------------------------------------------



