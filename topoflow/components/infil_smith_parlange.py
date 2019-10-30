"""
# This class defines a hydrologic infiltration component that implements
the well-known Smith-Parlange 3-parameter model (with gamma).

This class inherits from the infiltration "base class" in "infil_base.py".

See: Smith, R.E. (2002) Infiltration Theory for Hydrologic Applications,
Water Resources Monograph 15, AGU.
"""
## Copyright (c) 2009-2016, Scott D. Peckham
##
## January 2013   (Revised handling of input/output names).
## October 2012   (CSDMS Standard Names and BMI)
## January 2009  (converted from IDL)
## May, August 2009
## May 2010 (changes to initialize(), read_cfg_file(), etc.)

#-----------------------------------------------------------------------
# 
#  class infil_component         # (inherits from infil_base.py)
#
#      get_component_name()
#      get_attribute()           # (10/26/11)
#      get_input_var_names()     # (10/23/12)
#      get_output_var_names()    # (10/23/12)
#      get_var_name()            # (10/23/12)
#      get_var_units()           # (10/23/12)
#      -------------------------
#      initialize_layer_vars()
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
#      Smith_Parlange_Infil_Rate_v1
#      Smith_Parlange_Infil_Rate_1D
#      Smith_Parlange_Infil_Rate_3D
#      Smith_Parlange_Infil_Rate_v2  (uses previous 2)

#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.components import infil_base

from topoflow.utils import model_input
from topoflow.utils.tf_utils import TF_Print, TF_String

#-----------------------------------------------------------------------
class infil_component(infil_base.infil_component):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Infiltration_Smith_Parlange',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'InfilSmithParlange',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Infil_Smith_Parlange.cfg.in',
        'cfg_extension':      '_infil_smith_parlange.cfg',
        'cmt_var_prefix':     '/InfilSmithParlange/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Infil_Smith_Parlange.xml',
        'dialog_title':       'Infiltration: Smith-Parlange Parameters',
        'time_units':         'seconds' }

    _input_var_names = [
        'atmosphere_water__rainfall_volume_flux',           # (P_rain)  
        'glacier_ice__melt_volume_flux',                    # (MR)
        'land_surface__elevation',                          # (elev)  ##### IS THIS USED?
        'land_surface_water__evaporation_volume_flux',      # (ET)
        'snowpack__melt_volume_flux',                       # (SM)
        'soil_water_sat-zone_top_surface__elevation' ]      # (h_table)

        ##  'land_surface_water__baseflow_volume_flux',     # (GW)

    #-----------------------------------------------------------------
    # These are input vars provided by the user, but can also
    # be retrieved by other components as "output_vars".
    #-----------------------------------------------------------------
    #    'soil_water__green-ampt_capillary_length',      # G
    #    'soil_water__initial_hydraulic_conductivity',   # Ki
    #    'soil_water__initial_volume_fraction',          # qi
    #    'soil_water__saturated_hydraulic_conductivity', # Ks
    #    'soil_water__saturated_volume_fraction',        # qs
    #
    #-----------------------------------------------------------------
    # 'land_water__domain_time_integral_of__infiltration_rate'
    # vs. 'basin_cumulative_infiltrated_water_volume'.
    #-----------------------------------------------------------------
    # Does "ground water" connote the "saturated zone" ?
    #-----------------------------------------------------------------
    # "vertical" vs. "downward" doesn't establish a sign convention.
    #-----------------------------------------------------------------
    # "soil_water" vs. "subsurface_water" or "ground_water".
    #-----------------------------------------------------------------
    # The parameter G appears in Green-Ampt but seems to be much
    # older, according to Smith's book, p. 69.
    #-----------------------------------------------------------------
    _output_var_names = [
        'model__time_step',                   # dt
         # 'model_grid_cell__area',          # da
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux', # vol_IN
        'soil_surface_water__infiltration_volume_flux',                   # IN
        'soil_surface_water__time_integral_of_infiltration_volume_flux',  # I
        'soil_water__green-ampt_capillary_length',         # G (Also used for S-P.)
        'soil_water__potential_infiltration_volume_flux',  # fc
        'soil_water__initial_hydraulic_conductivity',      # Ki
        'soil_water__initial_volume_fraction',             # qi
        'soil_water__saturated_hydraulic_conductivity',    # Ks
        'soil_water__saturated_volume_fraction',           # qs
        'soil_water__smith-parlange_gamma_parameter',      # gam  #### NOT IN G-A.
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
        'soil_surface_water__infiltration_volume_flux':   'IN',
        'soil_surface_water__time_integral_of_infiltration_volume_flux': 'I',
        'soil_water__green-ampt_capillary_length':        'G',  # (Also used for S-P.)
        'soil_water__initial_hydraulic_conductivity':     'Ki',
        'soil_water__initial_volume_fraction':            'qi',
        'soil_water__potential_infiltration_volume_flux': 'fc',
        'soil_water__saturated_hydraulic_conductivity':   'Ks',
        'soil_water__saturated_volume_fraction':          'qs',
        'soil_water__smith-parlange_gamma_parameter':     'gam', ##### NOT IN G-A.
        'soil_water_flow__z_component_of_darcy_velocity': 'v',
        'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux': 'vol_Rg',
        'soil_water_sat-zone_top__recharge_volume_flux':  'Rg' }

    _var_units_map = {
        'atmosphere_water__rainfall_volume_flux':          'm s-1',   
        'glacier_ice__melt_volume_flux':                   'm s-1',
        'land_surface__elevation':                         'm',
        'land_surface_water__evaporation_volume_flux':     'm s-1',
        'snowpack__melt_volume_flux':                      'm s-1',
        'soil_water_sat-zone_top_surface__elevation':      'm',
        #------------------------------------------------------------
        'model__time_step': 's',
        # 'model_grid_cell__area': 'm2', 
        'soil_surface_water__time_integral_of_infiltration_volume_flux': 'm',
        'soil_surface_water__infiltration_volume_flux': 'm s-1',
        'soil_surface_water__domain_time_integral_of_infiltration_volume_flux': 'm3',
        'soil_water__green-ampt_capillary_length': 'm',  # (Also used for S-P.)
        'soil_water__initial_hydraulic_conductivity': 'm s-1',
        'soil_water__initial_volume_fraction': '1',
        'soil_water__potential_infiltration_volume_flux': 'm s-1',
        'soil_water__saturated_hydraulic_conductivity': 'm s-1',
        'soil_water__saturated_volume_fraction': '1',
        'soil_water__smith-parlange_gamma_parameter': '1',  ##### NOT IN G-A.
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
  
        return 'TopoFlow_Infiltration_Smith-Parlange'

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

        #----------------------------------------------------------
        # Notes: initialize_config_vars() calls read_cfg_file().
        #        If read_cfg_file() finds a variable "n_layers",
        #        then it calls initialize_layer_vars() so that
        #        subsequent layer variables - as indicated by a
        #        subscript in the CFG file - can be read directly
        #        into a list or array.
        #----------------------------------------------------------
        ## n_layers = 1  (before 11/15/16)
        n_layers = self.n_layers
        
        #-------------------------------------------------
        # Get arrays to store soil params for each layer
        #-------------------------------------------------
        self.soil_type = np.zeros([n_layers], dtype='<U200')
##        self.dz_val    = np.zeros([n_layers], dtype='Float64')    #### + dz3
##        self.nz_val    = np.zeros([n_layers], dtype='Int16')      #### + nz3
        #----------------------------------------------------------
        self.Ks_type  = np.zeros(n_layers, dtype='<U200')
        self.Ki_type  = np.zeros(n_layers, dtype='<U200')
        self.qs_type  = np.zeros(n_layers, dtype='<U200')
        self.qi_type  = np.zeros(n_layers, dtype='<U200')
        self.G_type   = np.zeros(n_layers, dtype='<U200')
        self.gam_type = np.zeros(n_layers, dtype='<U200')
        #--------------------------------------------------------        
        self.Ks_file  = np.zeros(n_layers, dtype='<U200')
        self.Ki_file  = np.zeros(n_layers, dtype='<U200')
        self.qs_file  = np.zeros(n_layers, dtype='<U200')
        self.qi_file  = np.zeros(n_layers, dtype='<U200')
        self.G_file   = np.zeros(n_layers, dtype='<U200')
        self.gam_file = np.zeros(n_layers, dtype='<U200')
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
        self.Ki  = list( np.zeros(n_layers, dtype='Float64') )
        self.qs  = list( np.zeros(n_layers, dtype='Float64') )
        self.qi  = list( np.zeros(n_layers, dtype='Float64') )
        self.G   = list( np.zeros(n_layers, dtype='Float64') )
        self.gam = list( np.zeros(n_layers, dtype='Float64') )
        #-------------------------------------------------------------        
#         self.Ks_list  = list( np.zeros(n_layers, dtype='Float64') )
#         self.Ki_list  = list( np.zeros(n_layers, dtype='Float64') )
#         self.qs_list  = list( np.zeros(n_layers, dtype='Float64') )
#         self.qi_list  = list( np.zeros(n_layers, dtype='Float64') )
#         self.G_list   = list( np.zeros(n_layers, dtype='Float64') )
#         self.gam_list = list( np.zeros(n_layers, dtype='Float64') )

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
                         #-----------------------------
                         self.is_scalar('Ks[0]'),
                         self.is_scalar('Ki[0]'),
                         self.is_scalar('qs[0]'),
                         self.is_scalar('qi[0]'),
                         self.is_scalar('G[0]'),
                         self.is_scalar('gam[0]')  ])


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
        #--------------------------------------------------------
        # (3/12/08) Same code should work if (self.n_layers eq 1)
        #--------------------------------------------------------
        # Convert from lists to arrays; same name. (11/15/16)
        #--------------------------------------------------------
        self.Ks  = self.build_layered_var( self.Ks )
        self.Ki  = self.build_layered_var( self.Ki )
        self.qs  = self.build_layered_var( self.qs)
        self.qi  = self.build_layered_var (self.qi )
        self.G   = self.build_layered_var( self.G )
        self.gam = self.build_layered_var( self.gam )  #######
     
    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def update_infil_rate(self):

        #-----------------------------------------------------
        # This function is not totally correct but is stable
        # and gives results similar to the correct method.
        #-----------------------------------------------------
        Smith_Parlange_Infil_Rate_v1(self)

        #-----------------------------------------------------
        # Numerically more correct method (??) but unstable
        #-----------------------------------------------------
        ## r = self.P_total
        ## r_last = ???
        ## n = self.time_index
        ## Smith_Parlange_Infil_Rate_v2(self, r, r_last, n)
       
    #   update_infil_rate()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        #-------------------------------------------------------
        # This method works for Green-Ampt and Smith-Parlange
        # but must be overridden for Richards 1D.
        #-------------------------------------------------------
        # NB! Green-Ampt and Smith-Parlange currently only
        #     support ONE layer (n_layers == 1).
        #-------------------------------------------------------
        self.Ks_unit  = []  # (empty lists to hold file objects)
        self.Ki_unit  = []
        self.qs_unit  = []
        self.qi_unit  = []
        self.G_unit   = []
        self.gam_unit = []

        for k in range(self.n_layers):
            self.Ks_file[k] = self.in_directory + self.Ks_file[k]
            self.Ki_file[k] = self.in_directory + self.Ki_file[k]
            self.qs_file[k] = self.in_directory + self.qs_file[k]
            self.qi_file[k] = self.in_directory + self.qi_file[k]
            self.G_file[k]  = self.in_directory + self.G_file[k]
            self.gam_file[k] = self.in_directory + self.gam_file[k]

            self.Ks_unit.append(  model_input.open_file(self.Ks_type[k],  self.Ks_file[k]) )
            self.Ki_unit.append(  model_input.open_file(self.Ki_type[k],  self.Ki_file[k]) )
            self.qs_unit.append(  model_input.open_file(self.qs_type[k],  self.qs_file[k]) )
            self.qi_unit.append(  model_input.open_file(self.qi_type[k],  self.qi_file[k]) )
            self.G_unit.append(   model_input.open_file(self.G_type[k],   self.G_file[k])  )
            self.gam_unit.append( model_input.open_file(self.gam_type[k], self.gam_file[k]) )

    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        rti = self.rti

        #-------------------------------------------------------
        # All grids are assumed to have data type of Float32.
        #-------------------------------------------------------
        # This method works for Green-Ampt and Smith-Parlange
        # but must be overridden for Richards 1D.
        #-------------------------------------------------------
        # NB! Green-Ampt and Smith-Parlange currently only
        #     support ONE layer (n_layers == 1).
        #------------------------------------------------------- 
        for k in range(self.n_layers):
            Ks = model_input.read_next(self.Ks_unit[k], self.Ks_type[k], rti)
            if (Ks is not None): self.Ks_list[k] = Ks

            Ki = model_input.read_next(self.Ki_unit[k], self.Ki_type[k], rti)
            if (Ki is not None): self.Ki_list[k] = Ki

            qs = model_input.read_next(self.qs_unit[k], self.qs_type[k], rti)
            if (qs is not None): self.qs_list[k] = qs

            qi = model_input.read_next(self.qi_unit[k], self.qi_type[k], rti)
            if (qi is not None): self.qi_list[k] = qi
            
            G  = model_input.read_next(self.G_unit[k], self.G_type[k], rti)
            if (G is not None): self.G_list[k] = G

            gam = model_input.read_next(self.gam_unit[k], self.gam_type[k], rti)
            if (gam is not None): self.gam_list[k] = gam
          
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):

        #-------------------------------------------------------
        # This method works for Green-Ampt and Smith-Parlange
        # but must be overridden for Richards 1D.
        #-------------------------------------------------------
        # NB! Green-Ampt and Smith-Parlange currently only
        #     support ONE layer (n_layers == 1).
        #-------------------------------------------------------        
        for k in range(self.n_layers):
            if (self.Ks_type[k]  != 'Scalar'): self.Ks_unit[k].close()        
            if (self.Ki_type[k]  != 'Scalar'): self.Ki_unit[k].close()
            if (self.qs_type[k]  != 'Scalar'): self.qs_unit[k].close()
            if (self.qi_type[k]  != 'Scalar'): self.qi_unit[k].close()
            if (self.G_type[k]   != 'Scalar'): self.G_unit[k].close()
            if (self.gam_type[k] != 'Scalar'): self.gam_unit[k].close()
            #------------------------------------------------------------
##            if (self.Ks_file[k]  != ''): self.Ks_unit[k].close()        
##            if (self.Ki_file[k]  != ''): self.Ki_unit[k].close()
##            if (self.qs_file[k]  != ''): self.qs_unit[k].close()
##            if (self.qi_file[k]  != ''): self.qi_unit[k].close()
##            if (self.G_file[k]   != ''): self.G_unit[k].close()
##            if (self.gam_file[k] != ''): self.gam_unit[k].close()
          
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
def Smith_Parlange_Infil_Rate_v1(self):

    #------------------------------------------------------------
    #Notes:  IN  = infiltration rate [m/s]
    #        Ks  = saturated hydraulic conductivity [m/s]
    #        Ki  = initial hydraulic conductivity [m/s]
    #        qs  = soil moisture content [dimless]
    #        qi  = soil initial moisture content [dimless]
    #         G  = capillary length scale [m]
    #         I  = cum. infiltration depth (since reset) [m]
    #         P  = precipitation rate [m/s]

    #        SM  = snowmelt rate [m/s]
    #         r  = (P + SM)  [m/s]

    #        Note that the infiltration rate has a max possible
    #        value of (P + SM) and asymptotes to Ks as the
    #        total infiltrated depth increases.  Need to reset
    #        this depth between "events", but haven't decided
    #        how to do this yet.

    #        Total infiltrated depth is incremented in the
    #        calling function, called Infiltration.
    #------------------------------------------------------------
    Ks  = self.Ks[0]   # (synonyms)
    Ki  = self.Ki[0]
    qs  = self.qs[0]
    qi  = self.qi[0]
    G   = self.G[0]
    gam = self.gam[0]
    
    dK = (Ks - Ki)
    dq = (qs - qi)
    t2 = np.exp((gam * self.I) / (G * dq)) - np.float64(1)
    IN = (gam * dK / t2) + Ks
    
    #---------------------------------------------------
    # Initially, IN = P_total, and all of the incoming
    # water infiltrates.  IN cannot exceed P_total.
    # Ponding time, Tp, is time until (IN < P_total).
    #---------------------------------------------------
    self.IN = np.minimum(IN, self.P_total)
    
    #-----------------------------------
    #Is P_total less than Ks anywhere ?
    #If so, set IN = P_total there.
    #-----------------------------------
    ## self.check_low_rainrate()  # (now done in infil_base.py)
    
#   Smith_Parlange_Infil_Rate_v1
#-----------------------------------------------------------------------
def Smith_Parlange_Infil_Rate_1D(self, r, r_last, n):

    #------------------------------------------------------------------
    #Notes:  This version was written on 5/8-11/06 using the method
    #        for dealing with variable rainfall that is described
    #        on p. 114-116 in Smith, R.E. (2002) Infiltration Theory
    #        for Hydrologic Applications, Water Resources Monograph
    #        15, AGU.  The method used previously gives very similar
    #        results, at least for sufficiently small time steps, but
    #        is technically incorrect.  Due to its simplicity, it may
    #        still be useful as a faster approximate method if we can
    #        estimate error for various parameter settings.

    #        The main difference with this version is that ponding
    #        times may occur between timesteps and two different
    #        cases are considered - ponding in the middle of a time
    #        step and ponding during a rise in rainrate at the very
    #        end of a timestep (MID_PONDING vs. END_PONDING). Also
    #        ponding time was not specifically computed before.

    #        In order to test the second case (END_PONDING) uniform

    #        precip was used with rates of 5e-5 m/s and 1e-4 m/s with
    #        durations of 8 minutes and 5 minutes.  The soil type of
    #        "silty clay" was selected to set infiltration parameters,
    #        the infiltration rate was set to 0.1 minutes (6 seconds,
    #        same as default for channel process).  The model was then
    #        run for 200 time steps and the infiltration rate was
    #        plotted (TEST_0D-v0.txt).  These settings triggered an
    #        END_PONDING event and results appeared to be correct.
    #        However, with these settings a virtually identical result
    #        was obtained by the much simpler incorrect method used by
    #        Green_AmptINfil_Rate_v1.  However, when infiltration
    #        time step was increased to 0.3 minutes (18 seconds),
    #        there was a fairly significant difference and had to
    #        change initial f from (r/2) to (0.4 * r) to avoid a jump
    #        in the infiltration curve.
    #------------------------------------------------------------------
    Ks  = self.Ks[0]   # (synonyms)
    Ki  = self.Ki[0]
    qs  = self.qs[0]
    qi  = self.qi[0]
    G   = self.G[0]
    gam = self.gam[0]
    
    #-----------------------------
    # If r <= Ks, ponding cannot
    # occur, so return f = r
    #-----------------------------
    if (r <= Ks):    
        return r
    
    #----------------------------------
    #Time at start and end of timestep
    #----------------------------------
    t_start = n * self.dt
    t_end = (n + 1) * self.dt
    
    #-----------------------
    #Other common variables
    #-----------------------
    dq  = (qs - qi)
    dK  = (Ks - Ki)
    fac = (G * dq / gam)
    
    PONDED = (self.tp != -float64(1))
    R_GT_KS = (r > Ks)
    if (not(PONDED) and R_GT_KS):    
        #--------------------------------------------
        #Compute test values of Ip  (Smith eqn 5.47)
        #--------------------------------------------
        #If (r_last lt Ks) then we'll have Ip1 < 0,
        #so ponding at end of interval can't occur.
        #But it seems it should be able to if r is
        #big enough or time step is long enough.
        #So set Ip1 large in this case ??
        #--------------------------------------------
        #*** if (r_last le *self.Ks) AND (n eq 1) then begin
        if (r_last <= Ks):    
            Ip1 = np.float64(999999)
        else:
            Ip1 = fac * np.log(float64(1) + (gam * dK / (r_last - Ks)))
        Ip2 = fac * np.log(float64(1) + (gam * dK / (r - Ks)))
        dIp = (Ip2 - self.I)
        
        #--------------------------------------------
        # Does ponding occur within this timestep ?
        #--------------------------------------------
        #   Ia(n) = sum_{k=0}^n (r_k * dt)
        #   Ia(n-1) < Ip(r_n) < Ia(n)
        #   0       < [Ip(r_n) - Ia(n-1)] < r_n*dt
        #   I = Ia(n-1), not updated yet by caller
        #------------------------------------------
        MID_PONDING = np.logical_and((dIp > 0), (dIp < (r * self.dt)))
        if (MID_PONDING):    
            print('MID_PONDING occurred.')
            #------------------------
            #For this case: Ip = Ip2
            #------------------------
            self.tp = t_start + (dIp / r)    #(r gt Ks, above)
            self.fp = r
        
        #------------------------------------------
        # Does ponding occur at end of timestep ?
        #------------------------------------------
        #   Ip(r_n) < Ia(n-1) < Ip(r_{n-1})
        #   [Ip(r_n) - Ia(n-1)] < 0  AND
        #    Ia(n-1) < Ip(r_{n-1})
        #----------------------------------------
        END_PONDING = np.logical_and(np.logical_and((dIp < 0), (self.I < Ip1)), np.logical_not(MID_PONDING))
        if (END_PONDING):    
            print('END_PONDING occurred.')
            #-----------------------------
            # For this case: Ip = self.I
            #-----------------------------
            self.tp = t_end
            self.fp = ((gam * dK) / (np.exp(self.I / fac) - np.float64(1))) + Ks
        
        #--------------------------
        # For debugging & testing
        #--------------------------
        if (logical_or(MID_PONDING, END_PONDING)):    
            TF_Print('tp = ' + TF_String(self.tp))
            TF_Print(' ')
        
    
    #--------------------------------------------------
    #If ponding still hasn't occurred, all "rain" will
    #infiltrate; infil rate = r.  RETURN to caller and
    #increment the cumulative depth, self.I by (r * dt)
    #--------------------------------------------------
    PONDED = (self.tp != -float64(1))
    if logical_not(PONDED):    
        return r
    
    #---------------------------------------------------
    #If ponding from the start or in previous 2 tests
    #then forge ahead and compute infiltration rate
    #via Newton-Raphson iteration using fp and tp.
    #---------------------------------------------------
    #First, set initial guess for iteration, f.
    #*self.fp was either set above or in a previous call.
    #---------------------------------------------------
    #   f = r / 2d    ;(converged usually,
    #                   but jumped to f < 0 sometimes)
    #   f = r         ;(converges, but jumps to f < 0)
    #   f  = 0.9 * r  ;(didn't converge)
    #   f = r / 100d  ;(didn't converge)
    #   f = 2d * r    ;(didn't converge)
    #   f = 1000d     ;(didn't converge)
    #---------------------------------------------------
    fp = self.fp
    f  = np.float64(0.41) * r
    #** f  = 0.5d * r    ;(converged for complex rain, but 0.4 and 0.6 didn't)
    #**                  ;(for simple rain, jumped to f < 0)
    #-----------------
    CONVERGED = False
    tol = np.float32(1E-8)
    n_tries = np.int16(0)
    n_max = np.int16(20)
    while (n_tries < n_max) and not(CONVERGED):
    #------------------------------------------
    # We're using eqn. (6.27) in Smith's book
    # for t(f) and solving it for f(t_end).
    #------------------------------------------
        term1 = (t_end - self.tp) * dK * (np.float64(1) - gam) / (G * dq)
        term2 = np.log(np.float64(1) + (gam * dK) / (f - Ks)) / gam
        term3 = np.log((f - Ki) / (f - Ks))
        term4 = np.log(np.float64(1) + (gam * dK) / (fp - Ks)) / gam
        term5 = np.log((fp - Ki) / (fp - Ks))
        h = term1 - term2 + term3 + term4 - term5   #(eq 0)
        #---------------------------------------------
        termA = dK / (f - Ks)
        termB = np.float64(1) / ((f - Ks) + (gam * dK))
        termC = np.float64(1) / (f - Ki)
        dh_df = termA * (termB - termC)
        del_f = -(h / dh_df)
        f = f + del_f
        CONVERGED = (np.absolute(del_f) <= tol)
        n_tries = n_tries + 1
    
    #--------------------------------------
    #Did Newton-Raphson fail to converge ?
    #--------------------------------------
    if not(CONVERGED):    
        TF_Print('********************************************')
        TF_Print('  ERROR:  Failure to converge.')
        TF_Print('  Max number of iterations exceeded')
        TF_Print('  while computing infiltration rate.')
        TF_Print('********************************************')
        TF_Print(' ')
        return np.float64(0)
    
    #----------------------------------------
    #Return infiltration rate at time, t_end
    #----------------------------------------
    self.IN = np.minimum(f, r)
    ## return np.minimum(f, r)
    
#  Smith_Parlange_Infil_Rate_1D
#-----------------------------------------------------------------------
def Smith_Parlange_Infil_Rate_3D(self, r, r_last, n):

    Ks  = self.Ks[0]   # (synonyms)
    Ki  = self.Ki[0]
    qs  = self.qs[0]
    qi  = self.qi[0]
    G   = self.G[0]
    gam = self.gam[0]
    
    #------------------------------------
    # Time at start and end of timestep
    #------------------------------------
    t_start = n * self.dt
    t_end = (n + 1) * self.dt
    
    #-------------------------
    # Other common variables
    #-------------------------
    dq  = (qs - qi)
    dK  = (Ks - Ki)
    fac = (G * dq / gam)
    
    #-------------------------------------------------------
    # NB!  wf is where ponding *could* occur during this
    # time step.  Ponding can't occur where (r le *self.Ks)     *** CHECK ***
    # and Ip2 is undefined in this case.
    #-------------------------------------------------------
    wf = np.where( np.logical_and((self.tp == -np.float64(1)), (r > Ks)) )
    n_flux = np.size( wf[0] )
    if (n_flux > 0):    
        #----------------------------------------------
        # Compute test values of Ip  (Smith eqn 5.47)
        #----------------------------------------------
        # If (r_last lt Ks) then we'll have Ip1 < 0,
        # so ponding at end of interval can't occur.
        # But it seems it should be able to if r is
        # big enough or time step is long enough.
        # So set Ip1 large in this case ??
        #----------------------------------------------
        if (r_last <= Ks):    
            Ip1 = np.float64(999999)
        else:    
            Ip1 = fac * np.log(np.float64(1) + (gam * dK / (r_last - Ks)))
        Ip2 = fac * np.log(np.float64(1) + (gam * dK / (r - Ks)))
        dIp = (Ip2 - self.I)
        
        #----------------------------------------
        # Does ponding occur for any grid cells
        # within this timestep ?
        #----------------------------------------
        dIp = (Ip2[wf] - (self.I)[wf])
        wm = np.where( np.logical_and((dIp > 0), (dIp < (r[wf] * self.dt))) )
        n_mid = np.size( wm[0] )   ########
        if (n_mid != 0):    
            #----------------------
            # Note that: Ip = Ip2
            #----------------------
            self.tp[ wf[wm] ] = t_start + (dIp[wm] / r[wf[wm]])
            self.fp[ wf[wm] ] = r[wf[wm]]
        
        #----------------------------------------
        # Does ponding occur for any grid cells
        # at the end of this timestep ?
        #----------------------------------------
        we = np.where( np.logical_and((dIp < 0), (self.I < Ip1)) )
        n_end = np.size( we[0] )
        if (n_end != 0):    
            #-------------------------
            # Note that: Ip = self.I
            #-------------------------
            self.tp[ wf[we] ] = t_end
            fp_array = ((gam * dK) / (np.exp(self.I / fac) - np.float64(1))) + Ks  #***
            self.fp[ wf[we] ] = fp_array[ wf[we] ]
        
    
    #--------------------------------------------
    # Initialize the grid of infiltration rates
    #--------------------------------------------
    ss = self.tp.shape
    ny = ss[0]
    nx = ss[1]
    #-------------------------------------------------
    ## ss = idl_func.size(self.tp, dimensions=True)
    ## nx = ss[0]
    ## ny = ss[1]
    f = np.zeros([ny, nx], dtype='Float64')
    
    #---------------------------------------------------
    # For grid cells where ponding has not yet occurred,
    # all "rain" should infiltrate and we set f=r.
    # This will include cases where r < Ks.
    #---------------------------------------------------------
    # Grid cells that *have* ponded but that now have r < Ks        **********
    # are special.  They should also be assigned f = r,
    # since we always have fc > Ks, but are not part of wf.
    # This case is handled at the end when we return (f < R)
    # (the lesser of f and R), since r < Ks < fc.
    # (self.I will be incremented by (f * dt) in caller.)
    #---------------------------------------------------------
    # wp is where ponding has occurred, just now or earlier
    # wf is where ponding has not occurred.
    #---------------------------------------------------------
    wp = np.where( self.tp != -np.float64(1) )
    n_ponded = np.size( wp[0] )
    wf = np.where( self.tp == -np.float64(1) )
    n_flux = np.size( wf[0] )
    if (n_flux != 0):    
        f[wf] = r[wf]
    if (n_ponded == 0):    
        return f
    
    #-----------------------------------------------------
    # For grid cells where ponding has occurred, we must
    # compute the infiltration rate via grid-based
    # Newton-Raphson iteration using fp and tp arrays.
    #-----------------------------------------------------
    # First, we must initialize fp for those grid cells.
    #-----------------------------------------------------
    f[wp] = np.float64(0.41) * r[wp]
    #** f[wp] = 0.5d * r[wp]
    #-----------------------
    CONVERGED = False             #(all grid cells)
    tol = np.float32(1E-8)            #(tolerance)
    n_tries = np.int16(0)
    n_max = np.int16(20)
    while (n_tries < n_max) and not(CONVERGED):
    #------------------------------------------
    # We're using eqn. (6.27) in Smith's book
    # for t(f) and solving it for f(t_end).
    #------------------------------------------
    # Some vars here may be scalars and others
    # may be 2D arrays, but h must be 2D ??
    # We compute h as a 2D array here because        ;**********
    # of this issue, which is not optimal.
    #-------------------------------------------
        term1 = (t_end - self.tp) * dK * (np.float64(1) - gam) / (G * dq)
        term2 = np.log(np.float64(1) + (gam * dK) / (f - Ks)) / gam
        term3 = np.log((f - Ki) / (f - Ks))
        term4 = np.log(np.float64(1) + (gam * dK) / (fp - Ks)) / gam
        term5 = np.log((fp - Ki) / (fp - Ks))
        h = term1 - term2 + term3 + term4 - term5   #(eq 0)
        #---------------------------------------------
        termA = dK / (f - Ks)
        termB = np.float64(1) / ((f - Ks) + (gam * dK))
        termC = np.float64(1) / (f - Ki)
        dh_df = termA * (termB - termC)
        #---------------------------------------------
        del_f = -(h / dh_df)
        f[wp] = f[wp] + del_f[wp]
        CONVERGED = (np.absolute(del_f[wp]) <= tol)
        n_tries = n_tries + 1
    
    #----------------------------------------
    # Did Newton-Raphson fail to converge ?
    #----------------------------------------
    if not(CONVERGED):    
        TF_Print('********************************************')
        TF_Print('  ERROR:  Failure to converge.')
        TF_Print('  Max number of iterations exceeded')
        TF_Print('  while computing infiltration rate.')
        TF_Print('********************************************')
        TF_Print(' ')
        return np.float64(0)
    
    #-------------------------------
    # Is r less than Ks anywhere ?
    # If so, set f = r there.
    #-------------------------------
    #***  Check_Low_Rainrate, self, f, r   ;(NO LONGER NEEDED, SEE ABOVE)
    
    #-------------------------------------------
    # Return infiltration rates at time, t_end
    #-------------------------------------------
    self.IN = np.minimum(f, r)
    ## return np.minimum(f, r)
    
#   Smith_Parlange_Infil_Rate_3D
#-----------------------------------------------------------------------
def Smith_Parlange_Infil_Rate_v2(self, r, r_last, n):

    #------------------------------------------------------------
    # Notes: This was written on 5/11/06 using information
    #        from Smith, R.E. (2002) Infiltration Theory for
    #        Hydrologic Applications, Water Resources Monograph
    #        15, AGU, pages 114-116.

    #        IN  = infiltration rate [m/s]
    #        Ks  = saturated hydraulic conductivity [m/s]
    #        Ki  = initial hydraulic conductivity [m/s]
    #        qs  = soil moisture content [dimless]
    #        qi  = soil initial moisture content [dimless]
    #         G  = capillary length scale [m]
    #         I  = cum. infiltration depth (since reset) [m]
    #         P  = precipitation rate [m/s]
    #        SM  = snowmelt rate [m/s]
    #         n  = time step (for computing t_start & t_end)

    #        Note that the infiltration rate has a max possible
    #        value of (P + SM) and asymptotes to Ks as the
    #        total infiltrated depth increases.  Need to reset
    #        this depth between "events", but haven't decided
    #        how to do this yet.

    #        Total infiltrated depth is incremented in the
    #        calling function, called Infiltration.
    #------------------------------------------------------------
    SINGLE_PROFILE = self.all_scalars  #(3/19/07)
    
    if (SINGLE_PROFILE):    
        Smith_Parlange_Infil_Rate_1D(self, r, r_last, n)
    else:    
        Smith_Parlange_Infil_Rate_3D(self, r, r_last, n)
    
#  Smith_Parlange_Infil_Rate_v2
#-----------------------------------------------------------------------
