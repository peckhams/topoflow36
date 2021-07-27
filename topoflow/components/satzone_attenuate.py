#
#  Copyright (c) 2021, Scott D. Peckham
#
#  July 2021. Started with satzone_darcy_layers.py.
#             This component takes Rg = P_rain * fraction.
#
#-----------------------------------------------------------------------
#  NOTES:  This file defines a groundwater component that treats
#          the water table recharge rate, Rg, to be a constant
#          fraction of P_rain (from meteorology component.)
#          This will route that fraction through a much slower
#          groundwater flowpath so that it re-emerges at a later
#          time as baseflow.  This should attenuate "flashy"
#          hydrographs as observed for Ethiopia.
#          This inherits from the groundwater "base class" in
#          "satzone_base.py".
#-----------------------------------------------------------------------
#
#  class satzone_component
#
#      get_component_name()
#      get_attribute()
#      get_input_var_names()
#      get_output_var_names()
#      get_var_name()
#      get_var_units() 
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.components import satzone_base

#-----------------------------------------------------------------------
class satzone_component( satzone_base.satzone_component ):

    _att_map = {
        'model_name':         'TopoFlow_Saturated_Zone_Attenuate',
        'version':            '3.6',        
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        'time_units':         'seconds',
        #-------------------------------------------------------------
        'comp_name':          'SatZoneAttenuate',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Satzone_Attenuate.cfg.in',
        'cfg_extension':      '_satzone_attenuate.cfg'}
        # 'cmt_var_prefix':     '/SatZoneAttenuate/Input/Var/',
        # 'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Satzone_Attenuate.xml',
        # 'dialog_title':       'Saturated Zone: Attenuate Parameters',

      
    #----------------------------------------------
    # What about ET? (Taking water off the ground
    # water surface???  (Bolton, 5/16/2012)
    #----------------------------------------------
    _input_var_names = [
        'atmosphere_water__rainfall_volume_flux',           # (P_rain@met)
        'channel_water_x-section__mean_depth',              # (d@channels)
        'land_surface_water__evaporation_volume_flux',      # (ET@evap)
        'soil_water_sat-zone_top__recharge_volume_flux'  ]  # (Rg@infil)     

    _output_var_names = [
        'land_surface__elevation',                          # elev
        'land_surface_water__baseflow_volume_flux',         # GW
        'land_surface_water__domain_time_integral_of_baseflow_volume_flux',  # vol_GW
        'model__time_step',                                 # dt
#         'model_soil_layer-0__porosity',                     # qs[0]
#         'model_soil_layer-0__saturated_thickness',          # y[0,:,:]
#         'model_soil_layer-0__thickness',                    # th[0,:,:]
#         'model_soil_layer-1__porosity',                     # qs[1]
#         'model_soil_layer-1__saturated_thickness',          # y[1,:,:]
#         'model_soil_layer-1__thickness',                    # th[1,:,:]
#         'model_soil_layer-2__porosity',                     # qs[2]
#         'model_soil_layer-2__saturated_thickness',          # y[2,:,:]
#         'model_soil_layer-2__thickness',                    # th[2,:,:]
        #----------------------------------------------
        # These are for *all* soil layers (not used).
        #----------------------------------------------
        # 'model_soil_layer__porosity',                       # qs
        # 'model_soil_layer__saturated_thickness',            # y
        # 'model_soil_layer__thickness',                      # th
        #----------------------------------------
        # The "top_layer" is same as "layer_0".
        #----------------------------------------
        'soil_water_sat-zone_top_surface__elevation',   # h_table  #############
        'soil_top-layer__porosity',                     # qs[0,:,:]
        'soil_top-layer__saturated_thickness',          # y[0,:,:]
        'soil_top-layer__thickness' ]                   # th[0,:,:]

        #-------------------------------------------
        # These are read from GUI/file, but could
        # still be returned.
        #-------------------------------------------
        # 'soil_water_sat-zone_top_surface__initial_elevation' ]   # h0_table

    #-------------------------------------------------------------------
    # Note: The variables qs, th and y are ndarrays.  If we define
    #       another variable as a slice or subset of these, such as
    #       qs_top = qs[0], or y_top = y[0,:,:], then they will
    #       also change whenever the main ndarray changes.
    #       To see this, try:
    #           >>> a = np.ones((3,3))
    #           >>> b = a[0,:]
    #           >>> print a
    #           >>> print b
    #           >>> a[0,:] = 2
    #           >>> print a
    #           >>> print b
    #       With this trick, we can avoid slices and subscripts in
    #       the var_name_map, which getattr and setattr don't support.
    #-------------------------------------------------------------------
    _var_name_map = {
        'atmosphere_water__rainfall_volume_flux':        'P_rain',
        'channel_water_x-section__mean_depth':           'd',      # channels comp
        'soil_water_sat-zone_top__recharge_volume_flux': 'Rg',
        #------------------------------------------------------------------------
        'land_surface__elevation': 'elev',
        'land_surface_water__baseflow_volume_flux': 'GW',  
        'land_surface_water__domain_time_integral_of_baseflow_volume_flux': 'vol_GW',
        'land_surface_water__evaporation_volume_flux': 'ET',
        'model__time_step': 'dt',
        #----------------------------------------------------------------
        # These are defined in satzone_base.py.  (9/22/14)
#         'model_soil_layer-0__porosity':            'qs_layer_0', ## 'qs[0]',
#         'model_soil_layer-0__saturated_thickness': 'y_layer_0',  ## 'y[0,:,:]',
#         'model_soil_layer-0__thickness':           'th_layer_0', ## 'th[0,:,:]',
#         'model_soil_layer-1__porosity':            'qs_layer_1',
#         'model_soil_layer-1__saturated_thickness': 'y_layer_1',
#         'model_soil_layer-1__thickness':           'th_layer_1',
#         'model_soil_layer-2__porosity':            'qs_layer_2',
#         'model_soil_layer-2__saturated_thickness': 'y_layer_2',
#         'model_soil_layer-2__thickness':           'th_layer_2',

        #----------------------------------------------
        # These are for *all* soil layers (not used).
        #----------------------------------------------
        # 'model_soil_layers__porosity':            'qs',
        # 'model_soil_layers__saturated_thickness': 'y',
        # 'model_soil_layers__thickness':           'th',
        #----------------------------------------
        # The "top_layer" is same as "layer_0".
        #----------------------------------------
        'soil_water_sat-zone_top_surface__elevation': 'h_table',
        'soil_top-layer__porosity':                   'qs_layer_0',  ## 'qs[0]',
        'soil_top-layer__saturated_thickness':        'y_layer_0',   ## 'y[0,:,:]',  
        'soil_top-layer__thickness':                  'th_layer_0' } ## 'th[0],
        
    _var_units_map = {
        'atmosphere_water__rainfall_volume_flux':              'm s-1',
        'channel_water_x-section__mean_depth': 'm',      # channels comp
        'soil_water_sat-zone_top__recharge_volume_flux': 'm s-1',
        #----------------------------------------------------------------
        'land_surface__elevation': 'm',
        'land_surface_water__baseflow_volume_flux': 'm s-1',
        'land_surface_water__domain_time_integral_of_baseflow_volume_flux': 'm3',
        'land_surface_water__evaporation_volume_flux': 'm s-1',
        'model__time_step': 's',         ############# CHECK UNITS
#         'model_soil_layer-0__porosity': '1',
#         'model_soil_layer-0__saturated_thickness': 'm',
#         'model_soil_layer-0__thickness':'m',
#         'model_soil_layer-1__porosity': '1',
#         'model_soil_layer-1__saturated_thickness': 'm',
#         'model_soil_layer-1__thickness': 'm',
#         'model_soil_layer-2__porosity': '1',
#         'model_soil_layer-2__saturated_thickness': 'm',
#         'model_soil_layer-2__thickness': 'm',
        #----------------------------------------------
        # These are for *all* soil layers (not used).
        #----------------------------------------------
        # 'model_soil_layers__porosity': '1',
        # 'model_soil_layers__saturated_thickness': 'm',
        # 'model_soil_layers__thickness': 'm',
        #----------------------------------------
        # The "top_layer" is same as "layer_0".
        #----------------------------------------
        'soil_water_sat-zone_top_surface__elevation': 'm',
        'soil_top-layer__porosity': '1',
        'soil_top-layer__saturated_thickness': 'm',
        'soil_top-layer__thickness': 'm' }

    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Satzone_Darcy_Layers'

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
    def update_water_volume_change(self):

        #-----------------------------------------------------------
        # Notes: Q_gw = total subsurface flux [m^3/s] (horizontal)
        #        Rg = rate at which water from the surface
        #             arrives at the water table [m/s]
        #        da = pixel area [m^2]
        #        dt = GW timestep [sec]
        #        w1 = IDs of pixels that flow in direction 1
        #        p1 = IDs of parent pixels for "w1 pixels" 
        #-----------------------------------------------------------
        fraction = 1.0
        self.Rg = self.P_rain * fraction  ###########
        
        #-------------------------------------------------------
        # Compute dV = total amount of water to be added to or
        # removed from the soil column during the subsurface
        # flow timestep.  Initialize with vertical recharge.
        #-------------------------------------------------------
        dt  = self.dt
        self.dV = (self.Rg * self.da) * dt
    
        #-----------------------------------------
        # Add contributions from neighbor pixels
        #-------------------------------------------------------------
        # Each grid cell passes flow to *one* downstream neighbor.
        # Note that multiple grid cells can flow toward a given grid
        # cell, so a grid cell ID may occur in d8.p1 and d8.p2, etc.
        #-------------------------------------------------------------
        if (self.d8.p1_OK):    
            self.dV[ self.d8.p1 ] += (dt * self.Q_gw[self.d8.w1])
        if (self.d8.p2_OK):    
            self.dV[ self.d8.p2 ] += (dt * self.Q_gw[self.d8.w2])
        if (self.d8.p3_OK):    
            self.dV[ self.d8.p3 ] += (dt * self.Q_gw[self.d8.w3])
        if (self.d8.p4_OK):    
            self.dV[ self.d8.p4 ] += (dt * self.Q_gw[self.d8.w4])
        if (self.d8.p5_OK):    
            self.dV[ self.d8.p5 ] += (dt * self.Q_gw[self.d8.w5])
        if (self.d8.p6_OK):    
            self.dV[ self.d8.p6 ] += (dt * self.Q_gw[self.d8.w6])
        if (self.d8.p7_OK):    
            self.dV[ self.d8.p7 ] += (dt * self.Q_gw[self.d8.w7])
        if (self.d8.p8_OK):    
            self.dV[ self.d8.p8 ] += (dt * self.Q_gw[self.d8.w8])

        #----------------------------------------------------
        # Subtract the amount that flows out to D8 neighbor
        #----------------------------------------------------
        self.dV -= (self.Q_gw * dt)  # (in place)

        #------------------------------
        # Convert dV to a water depth
        #------------------------------
        self.dzw = ( self.dV / self.da)

        #----------------
        # For debugging
        #----------------          
#         dzw_min = self.dzw.min()
#         dzw_max = self.dzw.max()
#         print('   dzw_min = ' + str(dzw_min))
#         print('   dzw_max = ' + str(dzw_max))
#         print ' '
                
    #   update_water_volume_change()
    #-------------------------------------------------------------------

