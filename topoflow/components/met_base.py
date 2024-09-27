"""
This file defines a "base class" for meteorology components as
well as any functions used by most or all meteorology methods.
That is, all meteorology components inherit methods from this class.
(This is currently the only one.)  The methods of this class should
be over-ridden as necessary for different methods of modeling
meteorology. This class, in turn, inherits from the "BMI base class"
in BMI_base.py.
"""
#  To DEBUG, set self.DEBUG = True in initialize().
#-----------------------------------------------------------------------
#
#  Copyright (c) 2001-2023, Scott D. Peckham
#  Feb 2024.  Add surface temperature module based on dew point
#             temperature
#  Sep 2023.  Fixed bug where PRECIP_ONLY was No but evaluated
#               to True when initialize() called read_input_files()
#               because it wasn't yet mapped to True or False.
#               Modified read_config_file() in BMI_base.py as a
#               general solution.  Using 0 or 1 should also work.
#             Added update_P_snow_integral() which must be
#               called after update_P_snow() in update().
#             Added update_P_rain_integral() which must be
#               called after update_P_rain() in update().
#               No longer calling update_P_integral().
#             Removed unneeded unit conversion to mbar in
#               update_vapor_pressure(). Get units from e_sat.
#             Removed unneeded snowpack vars from _input_var_names.
#  Aug 2023.  Fixed "CONVERT_K_TO_C" bug in read_input_files().
#             Fixed another possible bug.  Search for 2023-08-28.
#             Ability to set T_rain_snow temperature threshold.
#             Should update_P_integral() only use P_rain?  ################
#  Jul 2021.  Added start_year to CFG file and also in:
#             set_missing_cfg_options().  Now passed to functions
#             in solar_funcs.py that have a "year" keyword.
#  Jul 2020.  Separate initialize_input_file_vars().
#             Updated update_bulk_aero_conductance().
#             Updated read_input_files().
#  May 2020.  Separate set_rain_to_zero() method (readability)
#             Better PRECIP_ONLY support; set_computed_input_vars()
#             Bug fix: problem with read_input_files() call.
#             Separate: set_slope_angle(), set_aspect_angle().
#  Sep 2014.  Fixed sign error in update_bulk_richardson_number().
#             Ability to compute separate P_snow and P_rain.
#  Aug 2014.  New CSDMS Standard Names and clean up.
#  Nov 2013.  Converted TopoFlow to a Python package.
#
#  Jan 2013. Revised handling of input/output names.
#  Oct 2012. CSDMS Standard Names (version 0.7.9) and BMI.
#  May 2012. P is now a 1D array with one element and mutable,
#            so any comp with ref to it can see it change.
#  Jun 2010. update_net_shortwave_radiation(), etc.  
#  May 2010. Changes to initialize() and read_cfg_file().
#  Aug 2009. Improvements
#  Jan 2009. Converted from IDL.
#
#-----------------------------------------------------------------------
# Notes: Do we ever need to distinguish between a surface
#        temperature and snow temperature (in the snow) ?
#        Recall that a separate T_soil_x variable is used
#        to compute Qc.
#
#        Cp_snow is from NCAR CSM Flux Coupler web page

#        Does "land_surface_air__latent_heat_flux" make sense? (2/5/13)   
#-----------------------------------------------------------------------    
#
#  class met_component     (inherits from BMI_base.py)
#
#      get_component_name()
#      get_attribute()           # (10/26/11)
#      get_input_var_names()     # (5/15/12)
#      get_output_var_names()    # (5/15/12)
#      get_var_name()            # (5/15/12)
#      get_var_units()           # (5/15/12)
#      ---------------------
#      set_constants()
#      initialize()
#      update()
#      finalize()
#      ----------------------------
#      set_missing_cfg_options()
#      set_computed_input_vars()
#      set_slope_angle()
#      set_aspect_angle()
#      initialize_input_file_vars()   # (7/3/20)
#      initialize_computed_vars()
#      ----------------------------
#      update_P_integral()
#      update_P_max()
#      update_P_rain()    # (9/14/14, new method)
#      update_P_snow()    # (9/14/14, new method)
#      update_P_rain_integral()  # 9/23
#      update_P_snow_integral()  # 9/23
#      ------------------------------------
#      update_bulk_richardson_number()
#      update_bulk_aero_conductance()
#      update_sensible_heat_flux()
#      update_saturation_vapor_pressure()
#      update_vapor_pressure()
#      update_relative_humidity() # (8/29/23, unused, need spec. humidity)
#      update_dew_point()                    # (7/6/10)
#      update_precipitable_water_content()   # (7/6/10)
#      ------------------------------------
#      update_latent_heat_flux()
#      update_conduction_heat_flux()
#      update_advection_heat_flux()
#      ------------------------------------
#      update_julian_day()                   # (7/1/10)
#      update_albedo()                       # (3/1/24)
#      update_net_shortwave_radiation()      # (7/1/10)
#      update_em_air()                       # (7/1/10)
#      update_net_longwave_radiation()       # (7/1/10)
#      update_net_energy_flux()              # ("Q_sum")
#      ------------------------------------
#      open_input_files()
#      read_input_files()
#      set_rain_to_zero()     # (2020-05-05)
#      close_input_files()
#      ------------------------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
#  Functions:
#      compare_em_air_methods()
#
#-----------------------------------------------------------------------

import numpy as np
import os

from topoflow.components import solar_funcs as solar

from topoflow.utils import BMI_base
from topoflow.utils import model_input
from topoflow.utils import model_output
from topoflow.utils import rtg_files
from topoflow.utils import time_utils

#-----------------------------------------------------------------------
class met_component( BMI_base.BMI_component ):

    #-------------------------------------------------------------------
    _att_map = {
        'model_name':         'TopoFlow_Meteorology',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham, Lauren A. Bolotin',
        'grid_type':          'uniform',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'Meteorology',
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Meteorology.cfg.in',
        'cfg_extension':      '_meteorology.cfg',
        'cmt_var_prefix':     '/Meteorology/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Meteorology.xml',
        'dialog_title':       'Meteorology: Method 1 Parameters',
        'time_units':         'seconds' }

    #---------------------------------------------------------
    # Note that SWE = "snow water equivalent", but it really
    # just means "liquid_equivalent".
    #---------------------------------------------------------   
    _input_var_names = [ 'snowpack__depth', # h_snow
                        'glacier_ice__thickness', # h_ice
                         'snowpack__z_mean_of_mass-per-volume_density']  # rho_snow
#         'snowpack__depth' ]                            # h_snow
#         'snowpack__liquid-equivalent_depth',           # h_swe
#         'snowpack__melt_volume_flux' ]                 # SM   (MR used for ice?)

    #-----------------------------------------------------------
    # albedo, emissivity and transmittance are dimensionless.
    #-----------------------------------------------------------
    # "atmosphere_aerosol_dust__reduction_of_transmittance" vs.
    # This TF parameter comes from Dingman, App. E, p. 604.
    #-----------------------------------------------------------
    # There is an Optical_Air_Mass function in solar_funcs.py.
    # However, this quantity is not saved in comp state.
    #
    # "optical_path_length_ratio" vs. "optical_air_mass" OR
    # "airmass_factor" OR "relative_airmass" OR
    # "relative_optical_path_length"
    #-----------------------------------------------------------
    # Our term "liquid_equivalent_precipitation" is widely
    # used on the Internet, with 374,000 Google hits.
    #--------------------------------------------------------------
    # Note: "bulk exchange coefficient" has 2460 Google hits.
    #       It is closely related to a "transfer coefficient"
    #       for mass, momentum or heat.  There are no CF
    #       Standard Names with "bulk", "exchange" or "transfer".
    #
    # Zhang et al. (2000) use "bulk exchange coefficient" in a
    # nonstandard way, with units of velocity vs. unitless.
    #
    # Dn = bulk exchange coeff for the conditions of
    #      neutral atmospheric stability [m/s]
    # Dh = bulk exchange coeff for heat  [m/s]
    # De = bulk exchange coeff for vapor [m/s]
    #---------------------------------------------------------------
    # Now this component uses T_air to break the liquid-equivalent
    # precip rate into separate P_rain and P_snow components.
    # P_rain is used by channel_base.update_R()
    # P_snow is used by snow_base.update_depth()
    #---------------------------------------------------------------
    _output_var_names = [
        # 'atmosphere__optical_path_length_ratio',                           # M_opt [1]  (in solar_funcs.py)
        # 'atmosphere__von_karman_constant',                                 # kappa
        'atmosphere_aerosol_dust__reduction_of_transmittance',               # dust_atten  ##### (from GUI)
        'atmosphere_air-column_water-vapor__liquid-equivalent_depth',        # W_p ("precipitable depth")
        'atmosphere_bottom_air__brutsaert_emissivity_canopy_factor',         # canopy_factor
        'atmosphere_bottom_air__brutsaert_emissivity_cloud_factor',          # cloud_factor
        'atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance',   # De [m s-1], latent
        'atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance', # Dh [m s-1], sensible
        'atmosphere_bottom_air__emissivity',                                 # em_air
        'atmosphere_bottom_air__mass-per-volume_density',                    # rho_air
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity',       # Cp_air
        'atmosphere_bottom_air__neutral_bulk_aerodynamic_conductance',       # Dn [m s-1], neutral
        'atmosphere_bottom_air__pressure',                                   # p0
        'atmosphere_bottom_air__temperature',                                # T_air
        'atmosphere_bottom_air_flow__bulk_richardson_number',                # Ri [1]
        'atmosphere_bottom_air_flow__log_law_roughness_length',              # z0_air
        'atmosphere_bottom_air_flow__reference-height_speed',                # uz
        'atmosphere_bottom_air_flow__speed_reference_height',                # z
        'atmosphere_bottom_air_land_net-latent-heat__energy_flux',           # Qe [W m-2]
        'atmosphere_bottom_air_land_net-sensible-heat__energy_flux',         # Qh [W m-2]
        'atmosphere_bottom_air_water-vapor__dew_point_temperature',          # T_dew
        'atmosphere_bottom_air_water-vapor__partial_pressure',               # e_air # (insert "reference_height" ??)
        'atmosphere_bottom_air_water-vapor__relative_saturation',            # RH
        'atmosphere_bottom_air_water-vapor__saturated_partial_pressure',     # e_sat_air         
        'atmosphere_water__domain_time_integral_of_precipitation_leq-volume_flux',  # vol_P
        'atmosphere_water__domain_time_integral_of_rainfall_volume_flux',           # vol_PR
        'atmosphere_water__domain_time_integral_of_snowfall_leq-volume_flux',       # vol_PS
        'atmosphere_water__domain_time_max_of_precipitation_leq-volume_flux',     # P_max
        'atmosphere_water__precipitation_leq-volume_flux',                        # P [m s-1]
        'atmosphere_water__rainfall_volume_flux',            # P_rain [m s-1] (liquid)      
        'atmosphere_water__snowfall_leq-volume_flux',        # P_snow [m s-1]
        'earth__standard_gravity_constant',                  # g   [m s-2]
        'land_surface__albedo',                              # albedo
        'land_surface__aspect_angle',                        # alpha  (from GUI)
        'land_surface__emissivity',                          # em_surf
        'land_surface__latitude',                            # lat_deg [degrees]
        'land_surface__longitude',                           # lon_deg [degrees]
        'land_surface__slope_angle',                         # beta  (from GUI)
        'land_surface__temperature',                         # T_surf   ### OR JUST "land__temperature"?
        # 'land_surface_air__temperature',                          # T_air
        'land_surface_air_water-vapor__partial_pressure',           # e_surf # (insert "reference_height" ??)
        'land_surface_air_water-vapor__saturated_partial_pressure', # e_sat_surf
        'land_surface_net-longwave-radiation__energy_flux',  # Qn_LW [W m-2]
        'land_surface_net-shortwave-radiation__energy_flux', # Qn_SW [W m-2]
        'land_surface_net-total-energy__energy_flux',        # Q_sum [W w-2]
        'model__time_step',                                  # dt
        'physics__stefan_boltzmann_constant',                # sigma  [W m-2 K-4]
        'physics__von_karman_constant',                      # kappa  [1]
        'water__mass-specific_latent_fusion_heat',           # Lf     [J kg-1]
        'water__mass-specific_latent_vaporization_heat',     # Lv     [J kg-1]
        'water-liquid__mass-per-volume_density' ]            # rho_H2O
        
    #-----------------------------------------
    # These are used only in solar_funcs.py
    # Later, create a Radiation component.
    #---------------------------------------------
    # Should we allow "day" as a base quantity ?
    # "day_length" is confusing.  Think about "date" also.
    # Maybe something like:
    #
    #    "earth__mean_solar_rotation_period"
    #    "earth__sidereal_rotation_period"
    #    "earth__stellar_rotation_period"   (relative to "fixed stars")
    #         maybe:  "earth__complete_rotation_period" ??
    #
    #    OR:
    #    "earth_mean_solar_day__duration"
    #    "earth_sidereal_day__duration"
    #    "earth_stellar_day__duration"
    #
    #    OR perhaps:
    #    "earth_mean_solar_day__rotation_period"
    #    "earth_sidereal_day__rotation_period"
    #    "earth_stellar_day__rotation_period"
    #
    #    "stellar rotation period" gives 84,500 Google hits.
    #    "solar_rotation_period" gives 41,100 Google hits.
    #    "sidereal_roation_period" gives 86,000 Google hits.
    #    "stellar day" gives 136,000 Google hits (but many unrelated).
    #
    #    NB! "stellar_rotation_period" is ambiguous since it is also
    #         used for the rotation period of a star.
    #
    #    "earth_mean_solar_day__hour_count"  ("standard_day" ?)
    #    "earth_sidereal_day__hour_count"
    #    "earth_sidereal_day__duration"
    #    "earth__rotation_period"   = "sidereal_day"
    #
    #    "earth_stellar_day__period"  ??
    #    "earth_stellar_day__duration" ??
    #
    #------------------------------------------------------------------
    # For "earth__rotation_rate", it seems this should be based on
    # the sidereal day (23.93 hours) instead of the mean solar day.
    #------------------------------------------------------------------
    # There are at least a few online sources that use both terms:
    # "equivalent latitude" and "equivalent longitude".  See:
    # "The Construction and Application of a Martian Snowpack Model".
    #------------------------------------------------------------------
    # Adopt the little-used term:  "topographic_sunrise" ?
    # Or maybe "illuminated_topography", or "local_sunrise" ??
    #------------------------------------------------------------------
    # For angle relations between the earth and the sun, should we
    # just use the adjective "solar" in the quantity name or include
    # sun in the object name?  We could also use terms like:
    #     earth_to_sun__declination_angle
    #     earth_to_sun__right_ascension_angle
    #
    #------------------------------------------------------------------
    # The adjective "local" in "earth_local_apparent_noon__time"
    # may be helpful in other contexts such as:
    # 'earth__local_longitude' and 'land_surface__local_elevation'.
    #------------------------------------------------------------------    
    # 'earth__autumnal_equinox_date',
    # 'earth__autumnal_equinox_time',
    # 'earth_axis__ecliptic_tilt_angle',  # tilt_angle
    # 'earth__julian_day_number',         ########
    # 'earth__julian_day_angle',
    # 'earth__local_apparent_noon_time'
    # 'earth__mean_radius', 
    # 'earth__mean_solar_day_duration',   # (exactly 24 hours)
    # 'earth_orbit__eccentricity',
    # 'earth_orbit__period',     # (one year)
    # 'earth__perihelion_julian_day',  ######
    # 'earth__rotation_period',        ######
    # 'earth__rotation_rate', # Omega       ###### What about Angular Velocity ?
    # 'earth__sidereal_day_duration',     # (one rotation = 23.934470 hours)
    # 'earth__solar_declination_angle',
    # 'earth__solar_hour_angle',
    # 'earth__solar_irradiation_constant',  ## (or "insolation_constant" ??)
    # 'earth__solar_right_ascension_angle',
    # 'earth__solar_vertical_angle',   (complement of zenith angle)
    # 'earth__solar_zenith_angle',
    # 'earth__stellar_day_duration',  # (relative to the "fixed stars")
    # 'earth__summer_solstice_date',
    # 'earth__summer_solstice_time',
    # 'earth__topographic_sunrise_equivalent_latitude',
    # 'earth__topographic_sunrise_equivalent_longitude',  (flat_lon + offset) 
    # 'earth__topographic_sunrise_equivalent_longitude_offset',
    # 'earth__topographic_sunrise_time',
    # 'earth__topographic_sunset_time',
    # 'earth_true_solar_noon___time',  #####
    #     'earth_clock__true_solar_noon_time'
    # 'earth__vernal_equinox_date',
    # 'earth__vernal_equinox_time',
    # 'earth__winter_solstice_date',
    # 'earth__winter_solstice_time',
    #
    # What about a "slope_corrected" or "topographic" version of K_dir ?
    #
    # 'land_surface__backscattered_shortwave_irradiation_flux', # K_bs
    # 'land_surface__diffuse_shortwave_irradiation_flux',       # K_dif
    # 'land_surface__direct_shortwave_irradiation_flux',        # K_dir
    # 'land_surface__global_shortwave_irradiation_flux',        # K_glob = K_dif + K_dir
    #------------------------------------------------------------------


    #------------------------------------------------------------------   
    # Maybe we should rename "z" to "z_ref" and "uz" to "uz_ref" ?
    #------------------------------------------------------------------   
    _var_name_map = {
        'snowpack__depth': 'h_snow',
        'glacier_ice__thickness': 'h_ice',
        'snowpack__z_mean_of_mass-per-volume_density': 'rho_snow',
        # 'snowpack__liquid-equivalent_depth': 'h_swe',
        # 'snowpack__melt_volume_flux': 'SM',
        #-----------------------------------------------------------------
        # 'atmosphere__optical_path_length_ratio': 'M_opt',   # (in solar_funcs.py)
        # 'atmosphere__von_karman_constant': 'kappa',
        'atmosphere_aerosol_dust__reduction_of_transmittance': 'dust_atten',
        'atmosphere_air-column_water-vapor__liquid-equivalent_depth': 'W_p',   #########
        'atmosphere_bottom_air__brutsaert_emissivity_canopy_factor': 'canopy_factor',
        'atmosphere_bottom_air__brutsaert_emissivity_cloud_factor':  'cloud_factor',
        'atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance':   'De',
        'atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance': 'Dh',
        'atmosphere_bottom_air__emissivity': 'em_air',               
        'atmosphere_bottom_air__mass-per-volume_density': 'rho_air',
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'Cp_air',
        'atmosphere_bottom_air__neutral_bulk_heat_aerodynamic_conductance':  'Dn',
        'atmosphere_bottom_air__pressure':                           'p0',
        'atmosphere_bottom_air__temperature':                        'T_air',
        'atmosphere_bottom_air_flow__bulk_richardson_number':        'Ri',        
        'atmosphere_bottom_air_flow__log_law_roughness_length':      'z0_air', ## (not "z0")
        'atmosphere_bottom_air_flow__reference-height_speed':        'uz',
        'atmosphere_bottom_air_flow__speed_reference_height':        'z',
        'atmosphere_bottom_air_land_net-latent-heat__energy_flux':   'Qe',
        'atmosphere_bottom_air_land_net-sensible-heat__energy_flux': 'Qh',
        'atmosphere_bottom_air_water-vapor__dew_point_temperature':  'T_dew',  
        'atmosphere_bottom_air_water-vapor__partial_pressure':       'e_air',
        'atmosphere_bottom_air_water-vapor__relative_saturation':    'RH',
        'atmosphere_bottom_air_water-vapor__saturated_partial_pressure': 'e_sat_air',
        'atmosphere_water__domain_time_integral_of_precipitation_leq-volume_flux': 'vol_P',
        'atmosphere_water__domain_time_integral_of_rainfall_volume_flux':          'vol_PR',
        'atmosphere_water__domain_time_integral_of_snowfall_leq-volume_flux':      'vol_PS', 
        'atmosphere_water__domain_time_max_of_precipitation_leq-volume_flux': 'P_max',       
        'atmosphere_water__precipitation_leq-volume_flux': 'P',
        'atmosphere_water__rainfall_volume_flux':          'P_rain',    
        'atmosphere_water__snowfall_leq-volume_flux':      'P_snow',
        'earth__standard_gravity_constant':                'g',
        'land_surface__albedo':                            'albedo',
        'land_surface__aspect_angle':                      'alpha',
        'land_surface__emissivity':                        'em_surf',
        'land_surface__latitude':                          'lat_deg',
        'land_surface__longitude':                         'lon_deg',
        'land_surface__slope_angle':                       'beta',
        'land_surface__temperature':                       'T_surf',
         # 'land_surface_air__temperature': 'T_surf',
        'land_surface_air_water-vapor__partial_pressure':           'e_surf',
        'land_surface_air_water-vapor__saturated_partial_pressure': 'e_sat_surf',        
        'land_surface_net-longwave-radiation__energy_flux':         'Qn_LW',
        'land_surface_net-shortwave-radiation__energy_flux':        'Qn_SW',
        'land_surface_net-total-energy__energy_flux':               'Q_sum',
        'model__time_step':                                         'dt',
        'physics__stefan_boltzmann_constant':                       'sigma',
        'physics__von_karman_constant':                             'kappa',
        'water__mass-specific_latent_fusion_heat':                  'Lf',
        'water__mass-specific_latent_vaporization_heat':            'Lv',
        'water-liquid__mass-per-volume_density': 'rho_H2O' }

    #-----------------------------------------------------------------
    # Note: The "update()" function calls several functions with the
    #       MBAR keyword set to get units of "mbar" vs. "kPa".
    #-----------------------------------------------------------------
    # Note: We need to be careful with whether units are C or K,
    #       for all "thermal" quantities (e.g. thermal_capacity).
    #-----------------------------------------------------------------
    # Note: ARHYTHM had 3 "bulk exchange coefficients" that are all
    #       equal and therefore have the same units of [m s-1].
    #       Double-check that this is what is intended.   ##########
    #-----------------------------------------------------------------
    # Note: "atmosphere_column_water__liquid_equivalent_depth" has
    #       units of "cm", as in Dingman's book.  Make sure it gets
    #       used correctly in equations.
    #-----------------------------------------------------------------
    # Note: slope_angle and aspect_angle have units of RADIANS.
    #       aspect_angle is measured CW from north.
    #       RT files ending in "_mf-angle.rtg" and "fd-aspect.rtg"
    #       contain aspect values.  The former are in [0, 2 Pi]
    #       while the latter are in [-Pi, Pi] and both measure
    #       CCW from due east.  They are converted for use here.
    #-----------------------------------------------------------------    
    _var_units_map = {
        'snowpack__depth':                   'm',
        'glacier_ice__thickness':           'm',
        'snowpack__z_mean_of_mass-per-volume_density':     'kg m-3',
#         'snowpack__liquid-equivalent_depth': 'm',
#         'snowpack__melt_volume_flux':        'm s-1',
        #-------------------------------------------------------------
        # 'atmosphere__optical_path_length_ratio': '1',
        # 'atmosphere__von_karman_constant': '1',
        'atmosphere_aerosol_dust__reduction_of_transmittance':        '1',
        'atmosphere_air-column_water-vapor__liquid-equivalent_depth': 'cm',           # (see Notes above)
        'atmosphere_bottom_air__brutsaert_emissivity_canopy_factor':  '1',
        'atmosphere_bottom_air__brutsaert_emissivity_cloud_factor':   '1',
        'atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance': 'm s-1',   # (see Notes above)
        'atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance': 'm s-1', # (see Notes above)  
        'atmosphere_bottom_air__emissivity': '1',                            
        'atmosphere_bottom_air__mass-per-volume_density': 'kg m-3',
        'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity': 'J kg-1 K-1',   # (see Notes above)
        'atmosphere_bottom_air__neutral_bulk_heat_aerodynamic_conductance': 'm s-1',   # (see Notes above)
        'atmosphere_bottom_air__pressure':                               'mbar',
        'atmosphere_bottom_air__temperature':                            'deg_C',      # (see Notes above)
        'atmosphere_bottom_air_flow__bulk_richardson_number':            '1',
        'atmosphere_bottom_air_flow__log_law_roughness_length':          'm',
        'atmosphere_bottom_air_flow__reference-height_speed':            'm s-1',
        'atmosphere_bottom_air_flow__speed_reference_height':            'm',
        'atmosphere_bottom_air_land_net-latent-heat__energy_flux':       'W m-2',
        'atmosphere_bottom_air_land_net-sensible-heat__energy_flux':     'W m-2',
        'atmosphere_bottom_air_water-vapor__dew_point_temperature':      'deg_C',
        'atmosphere_bottom_air_water-vapor__partial_pressure':           'mbar', # (see Notes above)
        'atmosphere_bottom_air_water-vapor__relative_saturation':        '1',
        'atmosphere_bottom_air_water-vapor__saturated_partial_pressure': 'mbar',     # (see Notes above)
        'atmosphere_water__domain_time_integral_of_precipitation_leq-volume_flux': 'm3',
        'atmosphere_water__domain_time_integral_of_rainfall_volume_flux':          'm3',
        'atmosphere_water__domain_time_integral_of_snowfall_leq-volume_flux':      'm3',
        'atmosphere_water__domain_time_max_of_precipitation_leq-volume_flux': 'm s-1',        
        'atmosphere_water__precipitation_leq-volume_flux': 'm s-1',
        'atmosphere_water__rainfall_volume_flux': 'm s-1',      # (see Notes above)  
        'atmosphere_water__snowfall_leq-volume_flux': 'm s-1',  # (see Notes above)
        'earth__standard_gravity_constant': 'm s-2',
        'land_surface__albedo': '1',
        'land_surface__aspect_angle': 'radians',                      # (see Notes above)
        'land_surface__emissivity': '1',
        'land_surface__latitude': 'degrees',
        'land_surface__longitude': 'degrees',
        'land_surface__slope_angle': 'radians',
        'land_surface__temperature': 'deg_C',
        # 'land_surface_air__temperature': 'deg_C',
        'land_surface_air_water-vapor__partial_pressure':           'mbar',
        'land_surface_air_water-vapor__saturated_partial_pressure': 'mbar',
        'land_surface_net-longwave-radiation__energy_flux': 'W m-2',
        'land_surface_net-shortwave-radiation__energy_flux': 'W m-2',
        'land_surface_net-total-energy__energy_flux': 'W m-2',
        'model__time_step': 's',
        'physics__stefan_boltzmann_constant': 'W m-2 K-4',
        'physics__von_karman_constant': '1',
        'water__mass-specific_latent_fusion_heat': 'J kg-1',
        'water__mass-specific_latent_vaporization_heat': 'J kg-1',
        'water-liquid__mass-per-volume_density': 'kg m-3' }

    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Meteorology'

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
        return np.array( self._input_var_names )
    
    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):
 
        return np.array( self._output_var_names )
    
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
    def set_constants(self):

        #---------------------------------
        # Define some physical constants
        #---------------------------------
        self.g        = np.float64(9.81)    # [m s-2, gravity]
        self.kappa    = np.float64(0.408)   # [1]  (von Karman)
        self.rho_H2O  = np.float64(1000)    # [kg m-3]
        self.rho_air  = np.float64(1.2614)  # [kg m-3]
        self.Cp_air   = np.float64(1005.7)  # [J kg-1 K-1]
        self.Lv       = np.float64(2500000) # [J kg-1] Latent heat of vaporiz.
        self.Lf       = np.float64(334000)  # [J kg-1 = W s kg-1], Latent heat of fusion
        self.sigma    = np.float64(5.67E-8) # [W m-2 K-4]  (Stefan-Boltzman constant)
        self.C_to_K   = np.float64(273.15)  # (add to convert deg C to K)

        self.twopi         = np.float64(2) * np.pi
        self.one_seventh   = np.float64(1) / 7
        self.hours_per_day = np.float64(24)
        self.secs_per_day  = np.float64(3600) * self.hours_per_day

        #---------------------------
        # See update_latent_heat()
        #-----------------------------------------------------------        
        # According to Dingman (2002, p. 273), constant should
        # be 0.622 instead of 0.662 (Zhang et al., 2000, p. 1002).
        # Is this constant actually the dimensionless ratio of
        # the molecular weight of water to that of dry air ?
        #-----------------------------------------------------------
        self.latent_heat_constant = np.float64(0.622)
        # self.latent_heat_constant = np.float64(0.662)
        
        #-------------------------------------
        # Unit conversion factors for precip 
        #-------------------------------------
        self.mmph_to_mps = (np.float64(1) / np.float64(3600000))
        self.mmps_to_mps = (np.float64(1) / np.float64(1000))
        self.mph_to_mps  = (np.float64(1) / np.float64(3600))
        self.mps_to_mmph = np.float64(3600000)
        self.forever     = np.float64(999999999)  # [minutes]
        
        #------------------------------------------------
        # Only needed for method 1, where all rates and
        # durations are read as 1D arrays from GUI.
        # Method 1 may be removed in a future version.
        #------------------------------------------------
##        self.method1_rates     = None
##        self.method1_durations = None
##        self.method1_n_rates   = 0
        
    #   set_constants()           
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        ## self.DEBUG  = True   # override value in BMI_base.
        self.SILENT = SILENT
        if not(self.SILENT):
            print(' ')
            print('Meteorology component: Initializing...')
            
        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
                
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()
        self.initialize_config_vars()
        self.set_missing_cfg_options()
        ## print('self.P_file = ' + self.P_file )

        ### self.read_grid_info() # NOW IN initialize_config_vars()
        
        ## print '    Calling initialize_basin_vars()...'
        self.initialize_basin_vars()  # (5/14/10)
        ## print '    Calling initialize_time_vars()...'
        self.initialize_time_vars()

        #---------------------------------------------------
        # See "initialize_computed_vars()" method below.
        #---------------------------------------------------
        # Note: read_input_files() is called further down,
        #       and uses self.is_scalar(P), which uses
        #       self.P.ndim to determine if scalar.  So
        #       scalar/grid needs to be set at start.
        #       Check other input variables for this issue.    #######
        #       Needs to be before "Disabled" test.
        #-------------------------------------------------------
        # Data read from float32 files is converted to float64.
        #-------------------------------------------------------
        self.RAIN_OVER = False
        P_type = self.P_type.lower()
        ## dtype = 'float32'   ## TEST: 2020-05-04
        dtype = 'float64'

        if (P_type == 'scalar'):
            #-------------------------------------------------
            # (5/19/12) This makes P "mutable", which allows
            # its updated values to be seen by any component
            # that has a reference to it.
            #-------------------------------------------------
            # If P_type is 'scalar', then P was initialized
            # in initialize_config_vars().  But we still
            # need to initialize P_rain and P_snow here.
            #-------------------------------------------------
            self.P_rain = self.initialize_scalar(0, dtype=dtype)
            self.P_snow = self.initialize_scalar(0, dtype=dtype)

        #------------------------------------------------------
        # NB! "Sample steps" must be defined before we return
        #     Check all other process modules.
        #------------------------------------------------------
        if (self.comp_status.lower() == 'disabled'):
            if not(self.SILENT):
                print('Meteorology component: Disabled in CFG file.')
            self.e_air    = self.initialize_scalar(0, dtype=dtype)
            self.e_surf   = self.initialize_scalar(0, dtype=dtype)
            self.em_air   = self.initialize_scalar(0, dtype=dtype)
            self.Qn_SW    = self.initialize_scalar(0, dtype=dtype)
            self.Qn_LW    = self.initialize_scalar(0, dtype=dtype)
            self.Q_sum    = self.initialize_scalar(0, dtype=dtype)
            self.Qc       = self.initialize_scalar(0, dtype=dtype)
            self.Qa       = self.initialize_scalar(0, dtype=dtype)
            self.DONE     = True
            self.status   = 'initialized'
            return

        #----------------------------------------
        # Initialize vars to be read from files
        # Also see read_input_files() below.
        #----------------------------------------
        self.initialize_input_file_vars()

        #-----------------------------------------------
        # Read from files as needed to initialize vars 
        #-----------------------------------------------
        self.open_input_files()
        self.read_input_files()  # (initializes P)
        
        ## self.check_input_types()  # (not needed so far)

        #--------------------------------------------
        # Set any input variables that are computed
        #--------------------------------------------------
        # NOTE:  Must be called AFTER read_input_files().
        #--------------------------------------------------
        ## print('#### CALLING set_computed_input_vars() in met_base.py...')
        self.set_computed_input_vars()
                
        #-----------------------
        # Initialize variables
        #-----------------------
        ## print '    Calling initialize_computed_vars()...'
        self.initialize_computed_vars() # (after read_input_files)
                
        if not(self.PRECIP_ONLY):
            self.open_output_files() 
        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, dt=-1.0):

        #----------------------------------------------------------
        # Note: The read_input_files() method is first called by
        #       the initialize() method.  Then, the update()
        #       method is called one or more times, and it calls
        #       other update_*() methods to compute additional
        #       variables using input data that was last read.
        #       Based on this pattern, read_input_files() should
        #       be called at end of update() method as done here.
        #       If the input files don't contain any additional
        #       data, the last data read persists by default.
        #----------------------------------------------------------
        
        #--------------------------------
        # Has component been disabled ?
        #--------------------------------
        if (self.comp_status.lower() == 'disabled'):
            # Note: self.status should be 'initialized'.
            return

        #----------------------------------------
        # Read next met vars from input files ?
        #-----------------------------------------------------       
        # Note: read_input_files() is called by initialize()
        # and those values must be used for the "update"
        # calls before reading new ones.
        #-----------------------------------------------------
        self.status = 'updating'
        if (self.time_index > 0):
            self.read_input_files()
 
        #-------------------------------------------------           
        # Set precip to zero on edges of DEM for the sake
        # of mass balance checking (2023-09-11) ?
        # This doesn't seem to be helpful or needed.
        #-------------------------------------------------
#         if (self.P_type.lower() in ['grid', 'grid_sequence']):
#             self.P[ self.edge_IDs ] = 0.0
        
        #-------------------------------------------
        # Update computed values related to precip
        # P is initialized in initialize()
        #-------------------------------------------
        self.update_P_integral()       # update vol_P (leq)
        self.update_P_max()
        self.update_P_rain()
        self.update_P_snow()
        self.update_P_rain_integral()  # update vol_PR
        self.update_P_snow_integral()  # update vol_PS (leq)
                
        #-------------------------
        # Update computed values
        #-------------------------
        if not(self.PRECIP_ONLY):
            self.update_saturation_vapor_pressure(MBAR=True)
            self.update_vapor_pressure()
            self.update_dew_point() ###
            self.update_T_surf()
            self.update_saturation_vapor_pressure(MBAR=True, SURFACE=True)  ########
            self.update_bulk_richardson_number()
            self.update_bulk_aero_conductance()
            self.update_sensible_heat_flux()
            self.update_precipitable_water_content()  ###
            self.update_vapor_pressure(SURFACE=True)   ########
            self.update_latent_heat_flux()      # (uses e_air and e_surf)
            self.update_conduction_heat_flux()
            self.update_advection_heat_flux()
            self.update_julian_day()
            self.update_albedo(method='aging')
            self.update_net_shortwave_radiation()
            self.update_em_air()
            self.update_net_longwave_radiation()
            self.update_net_energy_flux()  # (at the end)
            
        #----------------------------------------
        # If placed here near bottom, then the
        # first values (e.g. P) are used twice.
        #----------------------------------------        
        ##  if (self.time_index > 0):
        ##      self.read_input_files()
      
        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # Components use own self.time_sec by default.
        #-----------------------------------------------
        if not(self.PRECIP_ONLY):
            self.write_output_files()
            ## self.write_output_files( time_seconds )

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        self.update_time( dt )
        self.status = 'updated'  # (OpenMI)

    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        #--------------------------------
        # Has component been disabled ?
        #--------------------------------
        if (self.comp_status.lower() == 'disabled'):
            # Note: self.status should be 'initialized'.
            return

        self.status = 'finalizing' 
        self.close_input_files()   ##  TopoFlow input "data streams"
        if not(self.PRECIP_ONLY):
            self.close_output_files()
        self.status = 'finalized'  # (OpenMI)

        if not(self.SILENT):
            self.print_final_report(comp_name='Meteorology component')
     
    #   finalize()
    #-------------------------------------------------------------------
    def set_missing_cfg_options(self):

        #-------------------------------------------------------  
        # Note: This is called in initialize() AFTER calling
        #       initialize_config_vars().  It is used to set
        #       newer toggles, etc. that may not have been
        #       set in the CFG file.
        #-------------------------------------------------------
        if not(hasattr(self, 'PRECIP_ONLY')):
            # See: set_computed_input_vars()
            self.PRECIP_ONLY = False  # don't set to 'No' here.
        if not(hasattr(self, 'P_factor')):
            self.P_factor = 1.0    # (2022-02-15)

        #----------------------------------------------------
        # Option to set the 50-50 rain-snow temp. threshold
        # See Jennings et al. (2018) for more info.
        #----------------------------------------------------
        if not(hasattr(self, 'T_rain_snow_type')):
            self.T_rain_snow_type = 'Scalar'
        if not(hasattr(self, 'T_rain_snow')):
            self.T_rain_snow = 1.0  # [deg C]
                              
        #----------------------------------------------------
        # Toggle to use SATTERLUND or BRUTSAERT methods
        # for computing e_air and em_air. (Not in GUI yet.)
        #----------------------------------------------------
        if not(hasattr(self, 'SATTERLUND')):
            self.SATTERLUND = False

        #------------------------------------------------
        # New option in model_input.read_next() to read
        # a precip time series from a NextGen CSV file
        #------------------------------------------------
        if not(hasattr(self, 'NGEN_CSV')):
            self.NGEN_CSV = False

        #--------------------------------------------
        # (2021-07-28) Added start_year to CFG file
        #--------------------------------------------
        if not(hasattr(self, 'start_year')):
            self.start_year = solar.Current_Year()
    
        #----------------------------------------------------
        # Options to save forcing/calculated variable outputs
        #----------------------------------------------------
        if not(hasattr(self, 'SAVE_QSW_GRIDS')):
            self.SAVE_QSW_GRIDS = False
        if not(hasattr(self, 'SAVE_QLW_GRIDS')):
            self.SAVE_QLW_GRIDS = False
        if not(hasattr(self, 'SAVE_TSURF_GRIDS')):
            self.SAVE_TSURF_GRIDS = False
        if not(hasattr(self, 'SAVE_ALB_GRIDS')):
            self.SAVE_ALB_GRIDS = False
        #-------------------------------------------
        if not(hasattr(self, 'SAVE_QSW_PIXELS')):
            self.SAVE_QSW_PIXELS = False
        if not(hasattr(self, 'SAVE_QLW_PIXELS')):
            self.SAVE_QLW_PIXELS = False
        if not(hasattr(self, 'SAVE_TSURF_PIXELS')):
            self.SAVE_TSURF_PIXELS = False
        if not(hasattr(self, 'SAVE_ALB_PIXELS')):
            self.SAVE_ALB_PIXELS = False
                                        
    #   set_missing_cfg_options()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        #-----------------------------------------------
        # Convert precip rate units from mm/h to m/s ?
        #-----------------------------------------------
        # NB! read_input_files() does this for files.
        #-----------------------------------------------
        if (self.P_type.lower() == 'scalar'):
            #------------------------------------------------
            # (2/7/13) Use "*=" here to preserve reference.
            #------------------------------------------------
            print('Scalar rainrate set to:', self.P, ' [mmph]')
            self.P *= self.mmph_to_mps  # (after print stmt)

        #------------------------------------------------
        # Convert PRECIP_ONLY string to True or False
        # 2023-09-01. This is now done in a general way
        # in read_config_file() in BMI_base.py.
        #------------------------------------------------
#         on_list = ['yes', 'true', 'on', '1']
#         setting = self.PRECIP_ONLY.lower()
#         self.PRECIP_ONLY = (setting in on_list)

        #---------------------------------------
        # Print info message about PRECIP_ONLY
        # See set_missing_cfg_options().
        #---------------------------------------
        if (self.PRECIP_ONLY and not(self.SILENT)):
            print('------------------------------------------')
            print(' NOTE: Since PRECIP_ONLY = True, output')
            print('       variables for met component will')
            print('       not be computed or saved to files.')
            print('       And evap component may not work.')
            print('------------------------------------------')
            print()

        #---------------------------------------------
        # Convert GMT_offset from string to int
        # because GUI can't use ints in droplist yet
        #---------------------------------------------
        self.GMT_offset = np.int16( self.GMT_offset )

        #------------------------------------------------
        # Convert start_month from string to integer
        # January should be 1.  See solar.Julian_Day().
        #------------------------------------------------
        month_list = ['January', 'February', 'March', 'April',
                      'May', 'June', 'July', 'August', 'September',
                      'October', 'November', 'December']
        self.start_month = month_list.index( self.start_month ) + 1

        #--------------------------------
        # See set_missing_cfg_options()
        #--------------------------------
  
        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt   = np.maximum(self.save_grid_dt,    self.dt)
        self.save_pixels_dt = np.maximum(self.save_pixels_dt,  self.dt)
        
    #   set_computed_input_vars()
    #-------------------------------------------------------------------
    def set_slope_angle(self):
    
        #-------------------------------------------------
        # Read slope grid & convert to slope angle, beta
        # NB!  RT slope grids have NaNs on edges.
        #-------------------------------------------------
        slopes = rtg_files.read_grid( self.slope_grid_file, self.rti,
                                      RTG_type='FLOAT' )
        self.slopes = slopes
        beta   = np.arctan( slopes )
        beta   = (self.twopi + beta) % self.twopi
        #---------------------------------------------
        w_nan = np.where( np.logical_not(np.isfinite(beta)) )
        n_nan = np.size(w_nan[0])
        if (n_nan != 0):    
            beta[ w_nan ] = np.float64(0)
        #------------------------------------------------------------------
        w_bad = np.where( np.logical_or( (beta < 0), (beta > np.pi / 2) ) )
        n_bad = np.size(w_bad[0])
        if (n_bad != 0):    
            print('ERROR: In met_base.py, some slope angles are out')
            print('       of range.  Returning without setting beta.')
            print()
            return

        self.beta = beta
            
    #   set_slope_angle()
    #-------------------------------------------------------------------
    def set_aspect_angle(self):
        #------------------------------------------------------
        # Read aspect grid.  Alpha must be CW from north.
        # NB!  RT aspect grids have NaNs on edges.
        #---------------------------------------------------------
        # RT files ending in "_mf-angle.rtg" and "fd-aspect.rtg"
        # contain aspect values.  The former are in [0, 2 Pi]
        # while the latter are in [-Pi, Pi] and both measure
        # CCW from due east.
        #---------------------------------------------------------
        aspects = rtg_files.read_grid( self.aspect_grid_file, self.rti,
                                       RTG_type='FLOAT' )
        alpha   = (np.pi / 2) - aspects
        alpha   = (self.twopi + alpha) % self.twopi
        #-----------------------------------------------
        w_nan = np.where( np.logical_not( np.isfinite(alpha) ) )
        n_nan = np.size( w_nan[0] )
        if (n_nan != 0):    
            alpha[ w_nan ] = np.float64(0)
            
        self.alpha = alpha
            
    #   set_aspect_angle()
    #-------------------------------------------------------------------
    def initialize_input_file_vars(self):
             
        #----------------------------------------
        # Initialize vars to be read from files
        #--------------------------------------------------------
        # initialize_var() initializes as scalar or grid.
        # Need this in order to use "update_var()", which makes
        # variables "mutable", so that its updated values can
        # be seen by any component that has a reference to it.
        #--------------------------------------------------------
        # If a variable's type (called [var]_type) is 'Scalar',
        # then it is initialized in initialize_config_vars().
        # For the 'Time_Series', 'Grid', and 'Grid_Sequence'
        # types, read_config_file() initializes them to '0.0',
        # so that self has the attribute.
        #--------------------------------------------------------
        dtype = 'float64'
        if (self.P_type.lower() != 'scalar'):
            self.P = self.initialize_var(self.P_type, dtype=dtype)
            self.P_rain = self.initialize_var(self.P_type, dtype=dtype)
            self.P_snow = self.initialize_var(self.P_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.T_air_type.lower() != 'scalar'):
            self.T_air = self.initialize_var(self.T_air_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.T_surf_type.lower() != 'scalar'):
            self.T_surf = self.initialize_var(self.T_surf_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.RH_type.lower() != 'scalar'):
            self.RH = self.initialize_var(self.RH_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.p0_type.lower() != 'scalar'):
            self.p0 = self.initialize_var(self.p0_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.uz_type.lower() != 'scalar'):
            self.uz = self.initialize_var(self.uz_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.z_type.lower() != 'scalar'):
            self.z = self.initialize_var(self.z_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.z0_air_type.lower() != 'scalar'):
            self.z0_air = self.initialize_var(self.z0_air_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.albedo_type.lower() != 'scalar'):
            self.albedo = self.initialize_var(self.albedo_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.em_surf_type.lower() != 'scalar'):
            self.em_surf = self.initialize_var(self.em_surf_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.dust_atten_type.lower() != 'scalar'):
            self.dust_atten = self.initialize_var(self.dust_atten_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.cloud_factor_type.lower() != 'scalar'):
            self.cloud_factor = self.initialize_var(self.cloud_factor_type, dtype=dtype)
        #--------------------------------------------------------------------
        if (self.canopy_factor_type.lower() != 'scalar'):
            self.canopy_factor = self.initialize_var(self.canopy_factor_type, dtype=dtype)
            
    #   initialize_input_file_vars()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        #------------------------------------------------------
        # Note: Some of these require "self.rti", which is
        #       only stored by read_grid_info() after the
        #       set_computed_input_vars() function is called.
        #       So these parts can't go there.
        #------------------------------------------------------
        
        #---------------------------------------
        # Add self.topo_directory to:
        #   slope_grid_file & aspect_grid_file
        #------------------------------------------------
        # Note: Maybe put these into read_input_files()
        #------------------------------------------------
        self.slope_grid_file  = (self.topo_directory + self.slope_grid_file) 
        self.aspect_grid_file = (self.topo_directory + self.aspect_grid_file)  
        self.set_slope_angle()
        self.set_aspect_angle()

        #------------------------------------------
        # For albedo calculations:
        # Start 'number of days since major snowfall' at 0
        self.n = 0    
        # How many time steps are in a three day period? How many days are in each timestep?:
        dt_per_3days = np.int64(259200/self.dt)
        self.days_per_dt = self.dt/86400         
        # Initialize array of rolling 3 day total snowfall    
        P_snow_3day = np.zeros(dt_per_3days)
        # Initialize grid for rolling 3 day total snowfall
        self.P_snow_3day_grid = self.initialize_grid(0)
        # Combine the grid and array so you have [# of time steps, # of rows, # of columns]
        self.P_snow_3day_grid = np.zeros((P_snow_3day.size, self.P_snow_3day_grid.shape[0], self.P_snow_3day_grid.shape[1])) 
        #------------------------------------------

    
        #---------------------------        
        # Create lon and lat grids
        #---------------------------
        if (self.rti.pixel_geom == 0):
            self.lon_deg = solar.Longitude_Grid( self.rti )
            self.lat_deg = solar.Latitude_Grid( self.rti )

##            print 'Lon grid ='
##            print self.lon_deg
##            print 'Lat grid ='
##            print self.lat_deg
            
            #-----------------------------
            # Write grids to RTG files ?
            #-----------------------------
##            lon_file = (self.out_directory + self.site_prefix + '_lons.bin')
##            rtg_files.write_grid( self.lon_deg, lon_file, self.rti )
##            lat_file = (self.out_directory + self.site_prefix + '_lats.bin')
##            rtg_files.write_grid( self.lat_deg, lat_file, self.rti )
        else:
            if not(self.SILENT):    
                print('SORRY: met_base.py cannot yet create lon and lat')
                print('       grids for this DEM because it uses UTM')
                print('       coordinates.  Using lat/lon for Denver, CO.')
                print()
            #--------------------------------------------
            # For now, use scalar values for Denver, CO
            #--------------------------------------------
            self.lon_deg = np.float64( -104.9841667 )
            self.lat_deg = np.float64( 39.7391667 )
            ## return

        #-------------------------------------------------
        # Initialize max precip rate with the first rate
        #------------------------------------------------
        # Note: Need this here because rate may be
        #       zero at the end of update_precip_rate()
        #------------------------------------------------
        # vol_P is used for mass balance check.
        #------------------------------------------------
        dtype = 'float64'
        P_max = self.P.max()       # (after read_input_files)
        ## self.P_max = self.P.max()
        self.P_max  = self.initialize_scalar( P_max, dtype=dtype)
        self.vol_P  = self.initialize_scalar( 0, dtype=dtype)
        self.vol_PR = self.initialize_scalar( 0, dtype=dtype)
        self.vol_PS = self.initialize_scalar( 0, dtype=dtype)
        
        #----------------------------------------------------
        # Initialize vars that are computed from input vars
        #----------------------------------------------------
#         self.update_bulk_richardson_number()
#         # self.update_bulk_aero_conductance()  ## Needs h_snow; not avail. yet.
#         self.update_sensible_heat_flux()
#         self.update_saturation_vapor_pressure(MBAR=True)
#         self.update_saturation_vapor_pressure(MBAR=True, SURFACE=True)  ########
#         self.update_vapor_pressure()
#         self.update_dew_point() ###
#         self.update_precipitable_water_content()  ###
#         self.update_vapor_pressure(SURFACE=True)   ########
#         self.update_latent_heat_flux()      # (uses e_air and e_surf)
#         self.update_conduction_heat_flux()
#         self.update_advection_heat_flux()
#         self.update_julian_day()
#         self.update_net_shortwave_radiation()
#         self.update_em_air()
#         self.update_net_longwave_radiation()
#         self.update_net_energy_flux()  # (at the end)
                   
        #----------------------------------------------------------
        # For using new framework which embeds references from
        # meteorology to snow, etc., these need to be defined
        # in the initialize step.  However, they will most likely
        # change from scalar to grid during update, so we need to
        # check that the reference isn't broken when the dtype
        # changes. (5/17/12)
        #---------------------------------------------------------- 
        # These depend on grids alpha and beta, so will be grids.
        #----------------------------------------------------------                
        self.Qn_SW  = self.initialize_grid( 0, dtype=dtype)
        self.Qn_LW  = self.initialize_grid( 0, dtype=dtype)
        self.Qn_tot = self.initialize_grid( 0, dtype=dtype)
        self.Q_sum  = self.initialize_grid( 0, dtype=dtype)
        #---------------------------------------------------------- 
        # These may be scalars or grids.
        #---------------------------------         
        self.Qe     = self.initialize_scalar( 0, dtype='float64')
        self.e_air  = self.initialize_scalar( 0, dtype='float64')
        self.e_surf = self.initialize_scalar( 0, dtype='float64')
        self.em_air = self.initialize_scalar( 0, dtype='float64')
        self.Qc     = self.initialize_scalar( 0, dtype='float64')
        self.Qa     = self.initialize_scalar( 0, dtype='float64')
         
        #------------------------------------
        # Initialize the decimal Julian day
        #---------------------------------------
        # (2021-07-28)  Added year=start_year.
        #---------------------------------------
        self.year       = self.start_year.copy()
        self.julian_day = solar.Julian_Day( self.start_month,
                                            self.start_day,
                                            self.start_hour,
                                            year=self.start_year)
        self.start_datetime = time_utils.get_datetime_str(
                                 self.start_year, self.start_month,
                                 self.start_day,  self.start_hour,
                                 0, 0)
        
    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def update_P_integral(self):

        #---------------------------------------------------
        # Notes: This can be used for mass balance checks,
        #        such as now done by print_final_report()
        #        in topoflow.py.
        #---------------------------------------------------     
        if (self.DEBUG):
            print('Calling update_P_integral()...')
     
        #-------------------------------------------------
        # Update mass total for P, sum over all pixels
        #-------------------------------------------------
        # We need to include total precip here, that is,
        # P = P_rain + P_snow (liquid equivalent), not
        # just P_rain, for use in a mass balance check.
        # P_rain and da are both either scalar or grid.
        #-------------------------------------------------   
        volume = np.double(self.P * self.da * self.dt)  # [m^3]
        if (np.size(volume) == 1):
            self.vol_P += (volume * self.rti.n_pixels)
        else:
            self.vol_P += np.sum(volume)
   
        # For testing
#         print('In update_P_integral:')
#         print('vol_P =', self.vol_P)
#         print('time  =', self.time)
#         print()
        
    #   update_P_integral()
    #-------------------------------------------------------------------
    def update_P_max(self):

        if (self.DEBUG):
            print('Calling update_P_max()...')
            
        #-----------------------------------------
        # Save the maximum precip. rate in [m/s]
        #-------------------------------------------
        # Must use "fill()" to preserve reference.
        #-------------------------------------------
        self.P_max.fill( np.maximum(self.P_max, self.P.max()) )
        ## self.P_max = np.maximum(self.P_max, self.P.max())

        #---------------
        # For testing
        #--------------
        ## print '##### P =', self.P * 1000 * 3600  # (mmph)
        ## print '##### P_max =', self.P_max * 1000 * 3600  # (mmph)

    #   update_P_max()  
    #-------------------------------------------------------------------
    def update_P_rain(self):

        #-----------------------------------------------------------
        # Note:  This routine is written so that it doesn't matter
        #        whether P and T_air are grids or scalars.
        #        For scalars: 1.5 * True = 1.5, 1.5 * False = 0.
        #        Here are the possible combinations for checking.
        #-----------------------------------------------------------
        # P       T_air     P_rain
        #----------------------------
        # scalar  scalar    scalar    
        # scalar  grid      grid
        # grid    scalar    grid
        # grid    grid      grid
        #----------------------------
        if (self.DEBUG):
            print('Calling update_P_rain()...')

        #-------------------------------------------------
        # P_rain is the precip that falls as liquid that
        # can contribute to runoff production.
        #-------------------------------------------------
        # P_rain is used by channel_base.update_R.
        #-------------------------------------------------
        P_rain = self.P * (self.T_air > self.T_rain_snow)
        
        if ((np.ndim( self.P_rain ) == 0) & (self.T_air_type.lower() == 'scalar')): 
            self.P_rain.fill( P_rain )   #### (mutable scalar)
        else:
            self.P_rain = P_rain
  
        if (self.DEBUG):
            if (self.P_rain.max() > 0):
                print('   >> Rain is falling...')

#         print('shape(P_rain) = ', self.P_rain.shape )
#         print('min(P_rain)   = ', self.P_rain.min() )
#         print('max(P_rain)   = ', self.P_rain.max() )
#         print()
               
        #--------------
        # For testing
        #--------------
        ## print 'shape(P)      =', shape(self.P)
        ## print 'shape(T_air)  =', shape(self.T_air)
        ## print 'shape(P_rain) =', shape(self.P_rain)
        ## print 'T_air         =', self.T_air

        #########################################
        #### Old note, to remember for later.
        #--------------------------------------------------
        # (2/7/13) We must use "*=" to preserve reference
        # if P is a "mutable scalar".
        #--------------------------------------------------
                             
    #   update_P_rain()
    #-------------------------------------------------------------------
    def update_P_snow(self):

        #----------------------------------------------------
        # Notes:  Rain and snow may fall simultaneously at
        #         different grid cells in the model domain.
        #----------------------------------------------------
        if (self.DEBUG):
            print('Calling update_P_snow()...')
            
        #-------------------------------------------------
        # P_snow is the precip that falls as snow or ice
        # that contributes to the snow depth.  This snow
        # may melt to contribute to runoff later on.
        #-------------------------------------------------
        # P_snow is a "water equivalent" volume flux
        # that was determined from a total volume flux
        # and a rain-snow temperature threshold.
        #-------------------------------------------------
        # P_snow is used by snow_base.update_depth.
        #-------------------------------------------------
        P_snow = self.P * (self.T_air <= self.T_rain_snow)
        
        if ((np.ndim( self.P_snow ) == 0) & (self.T_air_type.lower() == 'scalar')):
            self.P_snow.fill( P_snow )   #### (mutable scalar)
        else:
            self.P_snow = P_snow

        if (self.DEBUG):
            if (self.P_snow.max() > 0):
                print('   >> Snow is falling...')
                     
    #   update_P_snow()
    #-------------------------------------------------------------------
    def update_P_rain_integral(self):

        #---------------------------------------------------
        # Notes: This can be used for mass balance checks,
        #        such as now done by print_final_report()
        #        in topoflow.py.
        #        Must be called AFTER update_P_rain().
        #---------------------------------------------------     
        if (self.DEBUG):
            print('Calling update_P_rain_integral()...')
     
        #------------------------------------------------
        # Update mass total for P, sum over all pixels
        #------------------------------------------------
        # 2023-08-31. This one only uses P_rain.
        # P_rain and da are both either scalar or grid.
        #------------------------------------------------   
        volume = np.double(self.P_rain * self.da * self.dt)  # [m^3]
        if (np.size(volume) == 1):
            self.vol_PR += (volume * self.rti.n_pixels)
        else:
            self.vol_PR += np.sum(volume)
   
        if (self.DEBUG):
            print('  time, vol_PR =', self.time, ', ', self.vol_PR)
            print()
        
    #   update_P_rain_integral()
    #-------------------------------------------------------------------
    def update_P_snow_integral(self):

        #---------------------------------------------------
        # Notes: This can be used for mass balance checks,
        #        such as now done by print_final_report()
        #        in topoflow.py.
        #        Must be called AFTER update_P_snow().
        #---------------------------------------------------     
        if (self.DEBUG):
            print('Calling update_P_rain_integral()...')
     
        #----------------------------------------------------
        # Update mass total for P_snow, sum over all pixels
        #----------------------------------------------------
        # 2023-09-11. This one only uses P_snow.
        # P_snow and da are both either scalar or grid.
        #------------------------------------------------   
        volume = np.double(self.P_snow * self.da * self.dt)  # [m^3]
        if (np.size(volume) == 1):
            self.vol_PS += (volume * self.rti.n_pixels)
        else:
            self.vol_PS += np.sum(volume)
   
        if (self.DEBUG):
            print('  time, vol_PS =', self.time, ', ', self.vol_PS)
            print()
        
    #   update_P_snow_integral()
    #-------------------------------------------------------------------
    def update_bulk_richardson_number(self):

        if (self.DEBUG):
            print('Calling update_bulk_richardson_number()...')

        #---------------------------------------------------------------
        # (9/6/14)  Found a typo in the Zhang et al. (2000) paper,
        # in the definition of Ri.  Also see Price and Dunne (1976).
        # We should have (Ri > 0) and (T_surf > T_air) when STABLE.
        # This also removes problems/singularities in the corrections
        # for the stable and unstable cases in the next function.      
        #---------------------------------------------------------------
        # Notes: Other definitions are possible, such as the one given
        #        by Dingman (2002, p. 599).  However, this one is the
        #        one given by Zhang et al. (2000) and is meant for use
        #        with the stability criterion also given there.
        #---------------------------------------------------------------
        top     = self.g * self.z * (self.T_air - self.T_surf)  
        bot     = (self.uz)**2.0 * (self.T_air + np.float64(273.15))
        self.Ri = (top / bot)

        if (self.DEBUG):
            Ri_min = self.Ri.min()
            Ri_max = self.Ri.max()
            print('  Ri_min, Ri_max =', Ri_min, ', ', Ri_max)
            print()
            
    #   update_bulk_richardson_number()
    #-------------------------------------------------------------------        
    def update_bulk_aero_conductance(self):

        if (self.DEBUG):
            print('Calling update_bulk_aero_conductance()...')

        #----------------------------------------------------------------
        # Notes: Dn       = bulk exchange coeff for the conditions of
        #                   neutral atmospheric stability [m/s]
        #        Dh       = bulk exchange coeff for heat  [m/s]
        #        De       = bulk exchange coeff for vapor [m/s]
        #        h_snow   = snow depth [m]
        #        z0_air   = surface roughness length scale [m]
        #                   (includes vegetation not covered by snow)
        #        z        = height that has wind speed uz [m]
        #        uz       = wind speed at height z [m/s]
        #        kappa    = 0.408 = von Karman's constant [unitless]
        #        RI       = Richardson's number (see function)
        #----------------------------------------------------------------
        h_snow = self.h_snow  # (ref from new framework)
        
        #---------------------------------------------------
        # Compute bulk exchange coeffs (neutral stability)
        # using the logarithm "law of the wall".
        #-----------------------------------------------------
        # Note that "arg" = the drag coefficient (unitless).
        #-----------------------------------------------------
        # Dn will be a grid if any of the variables:
        #   z, h_snow, z0_air, or uz is a grid.
        #-----------------------------------------------------                
        arg = self.kappa / np.log((self.z - h_snow) / self.z0_air)
        Dn  = self.uz * (arg)**2.0
        
        T_AIR_SCALAR  = self.is_scalar('T_air')
        T_SURF_SCALAR = self.is_scalar('T_surf')

        if (T_AIR_SCALAR and T_SURF_SCALAR):
            if (self.T_air == self.T_surf):
               nw = 0
            else:
               nw = 1
        else:
            w  = (self.T_air != self.T_surf)   # (boolean array)
            nw = w.sum()
        
        if (nw == 0):
            #--------------------------------------------
            # All pixels are neutral. Set Dh = De = Dn.
            #--------------------------------------------
            self.Dn = Dn
            self.Dh = Dn
            self.De = Dn
            return
        
        #-------------------------------------
        # One or more pixels are not neutral
        # so make a correction using RI
        #---------------------------------------------
        # NB!  RI could be a grid when Dn is a
        # scalar, and this will change Dn to a grid.
        #---------------------------------------------
        # Ri = Richardson_Number(z, uz, T_air, T_surf)
        #--------------------------------------------
        # Before 12/21/07.  Has bug if RI is a grid
        #--------------------------------------------
        # w_stable = where(*T_air gt *T_surf, n_stable)
        # if (n_stable ne 0) then begin
        #     Dn[w_stable] = Dn[w_stable]/(1d + (10d * RI))
        # endif
        # w_unstable = where(*T_air lt *T_surf, n_unstable)
        # if (n_unstable ne 0) then begin
        #----------------------------------------------
        # Multiplication and substraction vs. opposites
        # for the stable case.  Zhang et al. (2000)
        # Hopefully not just a typo.
        #----------------------------------------------
        #    Dn[w_unstable] = Dn[w_unstable]*(1d - (10d * self.Ri))
        # endif
        
        #-----------------
        # After 12/21/07
        #------------------------------------------------------------
        # If T_air, T_surf or uz is a grid, then Ri will be a grid.
        # This version makes only one call to WHERE, so its faster.
        #------------------------------------------------------------
        # Multiplication and substraction vs. opposites for the
        # stable case (Zhang et al., 2000); hopefully not a typo.
        # It plots as a smooth curve through Ri=0.
        #------------------------------------------------------------
        # (9/7/14)  Modified so that Dn is saved, but Dh = De.
        #------------------------------------------------------------        
        Dh = Dn.copy()   ### (9/7/14.  Save Dn also.)
        nD = Dh.size
        nR = self.Ri.size
        if (nR > 1):    
            #--------------------------
            # Case where RI is a grid
            #--------------------------
            ws = (self.Ri > 0)     # where stable
            ns = ws.sum()
            wu = np.invert( ws )   # where unstable
            nu = wu.sum()
        
            if (nD == 1):
                #-----------------------------------
                # Convert Dh to a grid, same as Ri
                #-----------------------------------
                Dh = Dh + np.zeros( self.Ri.shape, dtype='float64' )

            #----------------------------------------------------------
            # If (Ri > 0), or (T_surf > T_air), then STABLE. (9/6/14)
            #---------------------------------------------------------- 
            # When ws and wu are boolean arrays, don't
            # need to check whether any are True. 
            #------------------------------------------- 
            # Dh[ws] = Dh[ws] / (np.float64(1) + (np.float64(10) * self.Ri[ws]))  
            # Dh[wu] = Dh[wu] * (np.float64(1) - (np.float64(10) * self.Ri[wu]))
            #-----------------------------------------------------------------------        
            if (ns != 0):  
                Dh[ws] = Dh[ws] / (np.float64(1) + (np.float64(10) * self.Ri[ws]))
            if (nu != 0):    
                Dh[wu] = Dh[wu] * (np.float64(1) - (np.float64(10) * self.Ri[wu]))
        else:    
            #----------------------------
            # Case where Ri is a scalar
            #--------------------------------
            # Works if Dh is grid or scalar
            #--------------------------------
            if (self.Ri > 0):    
                Dh = Dh / (np.float64(1) + (np.float64(10) * self.Ri))
            else:    
                Dh = Dh * (np.float64(1) - (np.float64(10) * self.Ri))

        #----------------------------------------------------
        # NB! We currently assume that these are all equal.
        #----------------------------------------------------
        self.Dn = Dn
        self.Dh = Dh
        self.De = Dh   ## (assumed equal)

        if (self.DEBUG):
            Dn_min = self.Dn.min()
            Dn_max = self.Dn.max()
            print('  Dn_min, Dn_max =', Dn_min, ', ', Dn_max)
            Dh_min = self.Dh.min()
            Dh_max = self.Dh.max()
            print('  Dh_min, Dh_max =', Dh_min, ', ', Dh_max)
            print()
       
    #   update_bulk_aero_conductance()
    #-------------------------------------------------------------------
    def update_sensible_heat_flux(self):

        #--------------------------------------------------------
        # Notes: All the Q's have units of W/m^2 = J/(m^2 s).
        #        Dh is returned by Bulk_Exchange_Coeff function
        #        and is not a pointer.
        #--------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_sensible_heat_flux()...')
            
        #---------------------
        # Physical constants
        #---------------------
        # rho_air = 1.225d   ;[kg m-3, at sea-level]
        # Cp_air  = 1005.7   ;[J kg-1 K-1]
        
        #-----------------------------
        # Compute sensible heat flux
        #-----------------------------
        delta_T = (self.T_air - self.T_surf)
        self.Qh = (self.rho_air * self.Cp_air) * self.Dh * delta_T

        if (self.DEBUG):
            Qh_min = self.Qh.min()
            Qh_max = self.Qh.max()
            print('  Qh_min, Qh_max =', Qh_min, ', ', Qh_max)
            print()

    #   update_sensible_heat_flux()
    #-------------------------------------------------------------------
    def update_saturation_vapor_pressure(self, MBAR=False,
                                         SURFACE=False):

        if (self.DEBUG):
            print('Calling update_saturation_vapor_pressure()...')

        #----------------------------------------------------------------
        # Notes: Saturation vapor pressure is a function of temperature.
        #        T is temperature in Celsius.  By default, the method
        #        of Brutsaert (1975) is used, but if the SATTERLUND
        #        keyword is set then the method of Satterlund (1979) is
        #        used.  When plotted, they look almost identical.  See
        #        the compare_em_air_method routine in this file.
        #        Dingman (2002) uses the Brutsaert method.
        #        Liston (1995, EnBal) uses the Satterlund method.

        #        By default, the result is returned with units of kPa.
        #        Set the MBAR keyword for units of millibars.
        #        100 kPa = 1 bar = 1000 mbars
        #                => 1 kPa = 10 mbars
        #----------------------------------------------------------------
        # NB!    Here, 237.3 is correct, and not a misprint of 273.2.
        #        See footnote on p. 586 in Dingman (Appendix D).
        #----------------------------------------------------------------
        # Also see: topoflow.utils.met_utils.py   #################
        #----------------------------------------------------------------
        # NOTE: If the temperature, T_air or T_surf, is constant in
        #       time, so that T_air_type or T_surf_type is in
        #       ['Scalar', 'Grid'], and if it has been initialized
        #       correctly, then there is no need to recompute e_sat.
        #----------------------------------------------------------------                        
        if (SURFACE):
#             HAVE_VAR   = hasattr(self, 'e_sat_surf'))
#             T_CONSTANT = (self.T_surf_type in ['Scalar', 'Grid'])
#             if (HAVE_VAR and T_CONSTANT): return
            T = self.T_surf
        else:
#             HAVE_VAR   = hasattr(self, 'e_sat_air')
#             T_CONSTANT = (self.T_air_type in ['Scalar', 'Grid'])
#             if (HAVE_VAR and T_CONSTANT): return
            T = self.T_air
        
        if not(self.SATTERLUND):    
            #------------------------------
            # Use Brutsaert (1975) method
            #------------------------------
            term1 = (np.float64(17.3) * T) / (T + np.float64(237.3))
            e_sat = np.float64(0.611) * np.exp(term1)  # [kPa]
        else:    
            #-------------------------------
            # Use Satterlund (1979) method     #### DOUBLE CHECK THIS (7/26/13)
            #-------------------------------
            term1 = np.float64(2353) / (T + np.float64(273.15))
            e_sat = np.float64(10) ** (np.float64(11.4) - term1)   # [Pa]
            e_sat = (e_sat / np.float64(1000))  # [kPa]

        #-----------------------------------
        # Convert units from kPa to mbars?
        #-----------------------------------
        if (MBAR):    
            e_sat = (e_sat * np.float64(10))   # [mbar]

        if (SURFACE):
            self.e_sat_surf = e_sat
        else:
            self.e_sat_air  = e_sat

        #----------------
        # For debugging
        #----------------
        if (self.DEBUG):
            if (SURFACE):
                e_sat_surf_min = self.e_sat_surf.min()
                e_sat_surf_max = self.e_sat_surf.max()
                print('  e_sat_surf_min, e_sat_surf_max =', 
                       e_sat_surf_min, ', ', e_sat_surf_max)
                print()
            else:
                e_sat_air_min = self.e_sat_air.min()
                e_sat_air_max = self.e_sat_air.max()
                print('  e_sat_air_min, e_sat_air_max =', 
                       e_sat_air_min, ', ', e_sat_air_max)
                print()
  
    #   update_saturation_vapor_pressure()
    #------------------------------------------------------------------- 
    def update_vapor_pressure(self, SURFACE=False):

        if (self.DEBUG):
            print('Calling update_vapor_pressure()...')

        #---------------------------------------------------
        # Notes: T is temperature in Celsius
        #        RH = relative humidity, in [0,1]
        #             by definition, it equals (e / e_sat)
        #        e has units of kPa.
        #---------------------------------------------------
        if (SURFACE):
            e_sat = self.e_sat_surf
        else:
            e_sat = self.e_sat_air

        #-------------------------------------------------         
        # RH is in [0,1], so e gets same units as e_sat.
        # So we never need to convert units of e.
        # 2023-09-01.  Removed unit conversion to mbar.
        #-------------------------------------------------    
        e = (self.RH * e_sat)

        if (SURFACE):
            self.e_surf = e
        else:
            self.e_air  = e

        #----------------
        # For debugging
        #----------------
        if (self.DEBUG):
            if (SURFACE):
                e_surf_min = self.e_surf.min()
                e_surf_max = self.e_surf.max()
                print('  e_surf_min, e_surf_max =', 
                       e_surf_min, ', ', e_surf_max)
                print()
            else:
                e_air_min = self.e_air.min()
                e_air_max = self.e_air.max()
                print('  e_air_min, e_air_max =', 
                       e_air_min, ', ', e_air_max)
                print()
                 
    #   update_vapor_pressure()
    #-------------------------------------------------------------------
    def update_dew_point(self):

        if (self.DEBUG):
            print('Calling update_dew_point()...')

        #-----------------------------------------------------------
        # Notes:  The dew point is a temperature in degrees C and
        #         is a function of the vapor pressure, e_air.
        #         Vapor pressure is a function of air temperature,
        #         T_air, and relative humidity, RH.
        #-----------------------------------------------------------

        #-------------------------------------------         
        # This formula needs e_air in kPa units.
        # See: Dingman (2002, Appendix D, p. 587).
        # 2023-09-01.  But it may contain a bug.
        #-------------------------------------------        
#         e_air_kPa = self.e_air / np.float64(10) # [mbar -> kPa]
#         log_vp    = np.log( e_air_kPa )
#         top = log_vp + np.float64(0.4926)
#         bot = np.float64(0.0708) - (np.float64(0.00421) * log_vp)
#         self.T_dew = (top / bot)    # [degrees C]
         #-------------------------------------------         
        # This formula needs e_air in Pa units.
        # See: Dingman (2015, 3.2.5, p. 114).
        #-------------------------------------------        
#         e_air_Pa = self.e_air * 100 # [mbar -> Pa]
#         self.T_dew = (top / bot)    # [degrees C]
        #-----------------------------------------------         
        # This formula needs e_air in mbar units.
        # See: https://en.wikipedia.org/wiki/Dew_point
        #-----------------------------------------------
        a = 6.1121   # [mbar]
        b = 18.678
        c = 257.14   # [deg C]
        # d = 234.5    # [deg C]
        log_term   = np.log( self.e_air / a)
        self.T_dew =  c * log_term / (b - log_term)  # [deg C]

        if (self.DEBUG):
            Td_min = self.T_dew.min()
            Td_max = self.T_dew.max()
            print('  T_dew_min, T_dew_max =', Td_min, ', ', Td_max, '[C]')
            print()
             
    #   update_dew_point()
    #-------------------------------------------------------------------
    def update_T_surf(self):
        # -------------------------------------------------
        # Estimate T_surf using T_dew (Raleigh et al. 2013).
        # Only run this function if T_surf is provided 
        # as a scalar or grid so that it still varies in time
        # -------------------------------------------------
        if ((self.T_surf_type.lower() == 'scalar') or (self.T_surf_type.lower() == 'grid')):
        # -------------------------------------------------
        # If snow and/or ice are present,  T_surf cannot
        # exceed 0 deg C
        # -------------------------------------------------
            T_surf = np.where(((self.h_snow > 0) | (self.h_ice > 0)), # where snow or ice exists
                              np.minimum(self.T_dew, np.float64(0)), # T_surf is either T_dew or 0, whichever is lower 
                              self.T_dew) # everywhere else, T_surf = T_dew
            self.T_surf = T_surf

    # update_T_surf()
    #-------------------------------------------------------------------
    def update_precipitable_water_content(self):

        if (self.DEBUG):
            print('Calling update_precipitable_water_content()...')

        #------------------------------------------------------------
        # Notes:  W_p is precipitable water content in centimeters,
        #         which depends on air temp and relative humidity.
        #------------------------------------------------------------
        arg      = np.float64( 0.0614 * self.T_dew )
        self.W_p = np.float64(1.12) * np.exp( arg )  # [cm]

        if (self.DEBUG):
            Wp_min = self.W_p.min()
            Wp_max = self.W_p.max()
            print('  W_p_min, W_p_max =', Wp_min, ', ', Wp_max, '[cm]')
            print()
            
    #   update_precipitable_water_content()
    #-------------------------------------------------------------------
    def update_latent_heat_flux(self):

        if (self.DEBUG):
            print('Calling update_latent_heat_flux()...')

        #--------------------------------------------------------
        # Notes:  Pressure units cancel out because e_air and
        #         e_surf (in numer) have same units (mbar) as
        #         p0 (in denom).
        #--------------------------------------------------------        
        # According to Dingman (2002, p. 273), constant should
        # be 0.622 instead of 0.662 (Zhang et al., 2000).
        #--------------------------------------------------------
        const   = self.latent_heat_constant
        factor  = (self.rho_air * self.Lv * self.De)
        delta_e = (self.e_air - self.e_surf)
        self.Qe = factor * delta_e * (const / self.p0)

        if (self.DEBUG):
            Qe_min = self.Qe.min()
            Qe_max = self.Qe.max()
            print('  T_air, T_surf =', self.T_air, ', ', self.T_surf, ' [deg C]')
            print('  RH = ', self.RH, '[0-1]')
            print('  h_snow = ', self.h_snow[34, 46], ' [m]')
            print('  uz = ', self.uz,'[m/s]')
            print('  Qe_min, Qe_max =', Qe_min, ', ', Qe_max, '[W m-2]')
            print()
            
    #   update_latent_heat_flux()
    #-------------------------------------------------------------------
    def update_conduction_heat_flux(self):

        #-----------------------------------------------------------------
        # Notes: The conduction heat flux from snow to soil for
        #        computing snowmelt energy, Qm, is close to zero.
        #        Currently, self.Qc = 0 in initialize().

        #        However, the conduction heat flux from surface and sub-
        #        surface for computing Qet is given by Fourier's Law,
        #        namely Qc = Ks(Tx - Ts)/x.

        #        All the Q's have units of W/m^2 = J/(m^2 s).
        #-----------------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_conduction_heat_flux()...')
            print('   Qc is always equal to 0 now.')
            print()

    #   update_conduction_heat_flux()
    #-------------------------------------------------------------------
    def update_advection_heat_flux(self):

        #------------------------------------------------------
        # Notes: Currently, self.Qa = 0 in initialize().
        #        All the Q's have units of W/m^2 = J/(m^2 s).
        #------------------------------------------------------        
        if (self.DEBUG):
            print('Calling update_advection_heat_flux()...')
            print('   Qa is always equal to 0 now.')
            print()
        
    #   update_advection_heat_flux()
    #-------------------------------------------------------------------
    def update_julian_day(self):

        if (self.DEBUG):
            print('Calling update_julian_day()...')

        #-----------------------------------   
        # Update the julian_day and year ?
        #-----------------------------------
        datetime = time_utils.get_current_datetime(
                              self.start_datetime,
                              self.time_min, time_units='minutes')
        (y,m1,d,h,m2,s) = time_utils.split_datetime_str( datetime, ALL=True )
        self.year = y
        ### self.year.fill(y)  # (if year is 0D ndarray, mutable)
 
        #----------------------------------
        # Update the *decimal* Julian day
        #----------------------------------
        self.julian_day = solar.Julian_Day( m1, d, hour_num=h, year=y)
        # print('Julian Day =', self.julian_day)
     
        #----------------------------------
        # Update the *decimal* Julian day
        #--------------------------------------------------
        # Before 2021-07-29, but doesn't stay in [1,365].
        #--------------------------------------------------
        ## self.julian_day += (self.dt / self.secs_per_day) # [days]
  
        #------------------------------------------
        # Compute the offset from True Solar Noon
        # clock_hour is in 24-hour military time
        # but it can have a decimal part.
        #------------------------------------------
        dec_part   = self.julian_day - np.int16(self.julian_day)
        clock_hour = dec_part * self.hours_per_day
        ## print '    Computing solar_noon...'
        solar_noon = solar.True_Solar_Noon( self.julian_day,
                                            self.lon_deg,
                                            self.GMT_offset,
                                            DST_offset=None,  #####
                                            year=self.year)
        ## print '    Computing TSN_offset...'
        self.TSN_offset = (clock_hour - solar_noon)    # [hours]

        if (self.DEBUG):
            print('  Julian day =', self.julian_day)
            print('  clock hour =', clock_hour)
            sn_min = solar_noon.min()
            sn_max = solar_noon.max()
            ts_min = self.TSN_offset.min()
            ts_max = self.TSN_offset.max()
            print('  solar_noon min, max =', sn_min, ', ', sn_max)
            print('  TSN_offset min, max =', ts_min, ', ', ts_max)
            print() 
            
    #   update_julian_day()
    #-------------------------------------------------------------------
    def update_albedo(self, method = 'aging'):
        # Only use this routine if time varying albedo is not supplied as an input:
        if ((self.albedo_type.lower() == 'scalar') or (self.albedo_type.lower() == 'grid')):
            if method == 'aging':
                albedo = self.albedo
                #------------------------------------------------
                # Dynamic albedo accounting for aging snow
                #------------------------------------------------
                # (Rohrer and Braun 1994): alpha = alpha0 + K * e^(-nr)
                # alpha = albedo
                # alpha0 = minimum snowpack albedo (~0.4)
                # K = constant (~0.44)
                # n = number of days since last major snowfall, at least 3 cm over 3 days
                # r = recession coefficient = 0.05 for temperatures < than 
                # 0 deg C, 0.12 for temperatures > 0 deg C
                #------------------------------------------------
                r = np.where((self.T_air > 0), 0.12, 0.05)
                K = 0.44
                alpha0 = 0.4

                self.P_snow_3day_grid = np.roll(self.P_snow_3day_grid, -1, axis=0) # you can roll on different axes (time axis), shape of the DEM and time axis and roll on the time axis
                ws_density_ratio = (self.rho_H2O/self.rho_snow)
                self.P_snow_3day_grid[np.size(self.P_snow_3day_grid, axis = 0)  - 1] = self.P_snow*self.dt*ws_density_ratio

                P_snow_3day_grid_total = np.sum(self.P_snow_3day_grid, axis=0) # maybe multipy by timestep here # also make sure to only sum over time axis
                # self.P_snow_3day_grid_total = P_snow_3day_grid_total # if you want to output and make sure it's working properly

                self.n = np.where((P_snow_3day_grid_total >= 0.03), 0, self.n)
                self.n = np.where((P_snow_3day_grid_total < 0.03), self.n+self.days_per_dt, self.n)
                snow_albedo = alpha0 + K * np.exp(-self.n*r)

                albedo = np.where((self.h_snow > 0), # where snow exists
                                snow_albedo, albedo)
                albedo = np.where(((self.h_snow == 0) & (self.h_ice > 0)), # where ice exists without snow
                                np.float64(0.3), albedo)
                albedo = np.where(((self.h_snow == 0) & (self.h_ice == 0)), # where there is no snow or ice (tundra)
                                np.float64(0.15), albedo)
                self.albedo = albedo
        #------------------------------------------------
        # Simple dynamic albedo depending on ice vs. snow vs. bare ground (tundra) using values from Dingman
        #------------------------------------------------
        if method == 'simple':
            albedo = self.albedo
            albedo = np.where((self.h_snow > 0), # where snow exists
                              np.float64(0.75), albedo)
            albedo = np.where(((self.h_snow == 0) & (self.h_ice > 0)), # where ice exists without snow
                              np.float64(0.3), albedo)
            albedo = np.where(((self.h_snow == 0) & (self.h_ice == 0)), # where there is no snow or ice (tundra)
                              np.float64(0.15), albedo)
            self.albedo = albedo
    #   update_albedo()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def update_net_shortwave_radiation(self):

        #---------------------------------------------------------
        # Notes:  If time is before local sunrise or after local
        #         sunset then Qn_SW should be zero.
        #---------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_net_shortwave_radiation()...')

        #--------------------------------
        # Compute Qn_SW for this time
        #--------------------------------
        K_cs = solar.Clear_Sky_Radiation( self.lat_deg,
                                           self.julian_day,
                                           self.W_p,
                                           self.TSN_offset,
                                           self.alpha,
                                           self.beta,
                                           self.albedo,
                                           self.dust_atten )
        
        #-------------------------------------------
        # 2024-03-06: Fix missing account for albedo 
        # in net shortwave radiation calcs
        # Dingman 3rd Edition 2015 Eqn. 6B1.1:
        # net shortwave = Kin * (1-albedo)
        #-------------------------------------------
        Qn_SW = K_cs * (1-self.albedo)
        
        if (np.ndim( self.Qn_SW ) == 0):
            self.Qn_SW.fill( Qn_SW )   #### (mutable scalar)
        else:
            self.Qn_SW[:] = Qn_SW  # [W m-2]

        if (self.DEBUG):
            Qsw_min = self.Qn_SW.min()
            Qsw_max = self.Qn_SW.max()
            print('  Qn_SW_min, Qn_SW_max =', Qsw_min, ', ', Qsw_max, '[W m-2]')
            print()
        
    #   update_net_shortwave_radiation()
    #-------------------------------------------------------------------
    def update_em_air(self):

        if (self.DEBUG):
            print('Calling update_em_air()...')

        #---------------------------------------------------------
        # NB!  The Brutsaert and Satterlund formulas for air
        #      emissivity as a function of air temperature are in
        #      close agreement; see compare_em_air_methods().
        #      However, we must pay close attention to whether
        #      equations require units of kPa, Pa, or mbar.
        #
        #             100 kPa = 1 bar = 1000 mbars
        #                => 1 kPa = 10 mbars
        #---------------------------------------------------------
        # NB!  Temperatures are assumed to be given with units
        #      of degrees Celsius and are converted to Kelvin
        #      wherever necessary by adding C_to_K = 273.15.
        #
        #      RH = relative humidity [unitless]
        #---------------------------------------------------------
        # NB!  I'm not sure about how F is added at end because
        #      of how the equation is printed in Dingman (2002).
        #      But it reduces to other formulas as it should.
        #---------------------------------------------------------
        T_air_K = self.T_air + self.C_to_K
        
        if not(self.SATTERLUND):
            #-----------------------------------------------------
            # Brutsaert (1975) method for computing emissivity
            # of the air, em_air.  This formula uses e_air with
            # units of kPa. (From Dingman (2002, p. 196).)
            # See notes for update_vapor_pressure().
            #-----------------------------------------------------
            e_air_kPa = self.e_air / np.float64(10)  # [kPa]
            F       = self.canopy_factor
            C       = self.cloud_factor
            term1   = (1.0 - F) * 1.72 * (e_air_kPa / T_air_K) ** self.one_seventh
            term2   = (1.0 + (0.22 * C ** 2.0))
            em_air = (term1 * term2) + F
        else:
            #--------------------------------------------------------
            # Satterlund (1979) method for computing the emissivity
            # of the air, em_air, that is intended to "correct
            # apparent deficiencies in this formulation at air
            # temperatures below 0 degrees C" (see G. Liston)
            # Liston cites Aase and Idso(1978), Satterlund (1979)
            #--------------------------------------------------------
            e_air_mbar = self.e_air
            eterm  = np.exp(-1 * (e_air_mbar)**(T_air_K / 2016) )
            em_air = 1.08 * (1.0 - eterm)

        #-------------------------- 
        # Update em_air, in-place
        #---------------------------------------------------------
        # NB! Currently, em_air is always initialized as scalar,
        #     but could change to grid after assignment.  Must
        #     determine if scalar or grid in initialize().
        #---------------------------------------------------------
        self.em_air = em_air
#         if (np.ndim( self.em_air ) == 0):
#             self.em_air.fill( em_air )   #### (mutable scalar)
#         else:
#             self.em_air[:] = em_air

        if (self.DEBUG):
            ema_min = self.em_air.min()
            ema_max = self.em_air.max()
            print('  em_air_min, em_air_max =', ema_min, ', ', ema_max)
            print()
                       
    #   update_em_air()    
    #-------------------------------------------------------------------
    def update_net_longwave_radiation(self):

        #----------------------------------------------------------------
        # Notes: Net longwave radiation is computed using the
        #        Stefan-Boltzman law.  All four data types
        #        should be allowed (scalar, time series, grid or
        #        grid stack).
        #
        #        Qn_LW = (LW_in - LW_out)
        #        LW_in   = em_air  * sigma * (T_air  + 273.15)^4
        #        LW_out  = em_surf * sigma * (T_surf + 273.15)^4
        #
        #        Temperatures in [deg_C] must be converted to
        #        [K].  Recall that absolute zero occurs at
        #        0 [deg_K] or -273.15 [deg_C].
        #
        #----------------------------------------------------------------
        # First, e_air is computed as:
        #   e_air = RH * 0.611 * exp[(17.3 * T_air) / (T_air + 237.3)]
        # Then, em_air is computed as:
        #   em_air = (1 - F) * 1.72 * [e_air / (T_air + 273.15)]^(1/7) *
        #             (1 + 0.22 * C^2) + F
        #----------------------------------------------------------------
        if (self.DEBUG):
            print('Calling update_net_longwave_radiation()...')

        #--------------------------------
        # Compute Qn_LW for this time
        #--------------------------------
        T_air_K  = self.T_air  + self.C_to_K
        T_surf_K = self.T_surf + self.C_to_K
        LW_in    = self.em_air  * self.sigma * (T_air_K)** 4.0
        LW_out   = self.em_surf * self.sigma * (T_surf_K)** 4.0

        #----------------------------------------------------        
        # 2023-08-29.  The next line was here before today,
        # and accounts for the amount of longwave radiation
        # from the air that is reflected from the surface.
        # See: https://daac.ornl.gov/FIFE/guides/
        #        Longwave_Radiation_UNL.html
        # It reduces the net longwave radiation.
        #----------------------------------------------------
        LW_out += ((1.0 - self.em_surf) * LW_in)

        self.Qn_LW[:] = (LW_in - LW_out)   # [W m-2]

        #--------------------------------------------------------------  
        # Can't do this yet.  Qn_LW is always initialized grid now
        # but will often be created above as a scalar. (9/23/14)
        #--------------------------------------------------------------
#         if (np.ndim( self.Qn_LW ) == 0):
#             self.Qn_LW.fill( Qn_LW )   #### (mutable scalar)
#         else:
#             self.Qn_LW[:] = Qn_LW  # [W m-2]

        if (self.DEBUG):
            Qlw_min = self.Qn_LW.min()
            Qlw_max = self.Qn_LW.max()
            print('  Qn_LW_min, Qn_LW_max =', Qlw_min, ', ', Qlw_max, '[W m-2]')
            print()
                    
    #   update_net_longwave_radiation()
    #-------------------------------------------------------------------
#     def update_net_total_radiation(self):
# 
#         #-------------------------------------------------------
#         # Notes: Added this on 9/11/14.  Not used;  see Q_sum.
#         #------------------------------------------------------------
#         #        Qn_SW = net shortwave radiation flux (solar)
#         #        Qn_LW = net longwave radiation flux (air, surface)
#         #------------------------------------------------------------       
#         if (self.DEBUG):
#             print('Calling update_net_total_radiation()...')
#             
#         Qn_tot = self.Qn_SW + self.Qn_LW   # [W m-2]
# 
#         if (np.ndim( self.Qn_tot ) == 0):
#             self.Qn_tot.fill( Qn_tot )   #### (mutable scalar)
#         else:
#             self.Qn_tot[:] = Qn_tot  # [W m-2]
#                        
#     #   update_net_total_radiation()
    #-------------------------------------------------------------------
    def update_net_energy_flux(self):

        if (self.DEBUG):
            print('Calling update_net_energy_flux()...')

        #------------------------------------------------------
        # Notes: Q_sum is used by "snow_energy_balance.py".
        #------------------------------------------------------
        #        Qm    = energy used to melt snowpack (if > 0)
        #        Qn_SW = net shortwave radiation flux (solar)
        #        Qn_LW = net longwave radiation flux (air, surface)
        #        Qh    = sensible heat flux from turbulent convection
        #                between snow surface and air
        #        Qe    = latent heat flux from evaporation, sublimation,
        #                and condensation
        #        Qa    = energy advected by moving water (i.e. rainfall)
        #                (ARHYTHM assumes this to be negligible; Qa=0.)
        #        Qc    = energy flux via conduction from snow to soil
        #                (ARHYTHM assumes this to be negligible; Qc=0.)
        #        Ecc   = cold content of snowpack = amount of energy
        #                needed before snow can begin to melt [J m-2]

        #        All Q's here have units of [W m-2].
        #        Are they all treated as positive quantities ?

        #        rho_air  = density of air [kg m-3]
        #        rho_snow = density of snow [kg m-3]
        #        Cp_air   = specific heat of air [J kg-1 K-1]
        #        Cp_snow  = heat capacity of snow [J kg-1 K-1]
        #                 = ???????? = specific heat of snow
        #        Kh       = eddy diffusivity for heat [m2 s-1]
        #        Ke       = eddy diffusivity for water vapor [m2 s-1]
        #        Lv       = latent heat of vaporization [J kg-1]
        #        Lf       = latent heat of fusion [J kg-1]
        #        ------------------------------------------------------
        #        Dn       = bulk exchange coeff for the conditions of
        #                   neutral atmospheric stability [m/s]
        #        Dh       = bulk exchange coeff for heat
        #        De       = bulk exchange coeff for vapor
        #        ------------------------------------------------------
        #        T_air    = air temperature [deg_C]
        #        T_surf   = surface temperature [deg_C]
        #        T_snow   = average snow temperature [deg_C]
        #        RH       = relative humidity [unitless] (in [0,1])
        #        e_air    = air vapor pressure at height z [mbar]
        #        e_surf   = surface vapor pressure [mbar]
        #        ------------------------------------------------------
        #        h_snow   = snow depth [m]
        #        z        = height where wind speed is uz [m]
        #        uz       = wind speed at height z [m/s]
        #        p0       = atmospheric pressure [mbar]
        #        T0       = snow temperature when isothermal [deg_C]
        #                   (This is usually 0.)
        #        z0_air   = surface roughness length scale [m]
        #                   (includes vegetation not covered by snow)
        #                   (Values from page 1033: 0.0013, 0.02 [m])
        #        kappa    = von Karman's constant [unitless] = 0.41
        #        dt       = snowmelt timestep [seconds]
        #----------------------------------------------------------------
        Q_sum = self.Qn_SW + self.Qn_LW + self.Qh + \
                self.Qe + self.Qa + self.Qc    # [W m-2]

        if (np.ndim( self.Q_sum) == 0):
            self.Q_sum.fill( Q_sum )   #### (mutable scalar)
        else:
            self.Q_sum[:] = Q_sum  # [W m-2]

        if (self.DEBUG):
            Qn_min = self.Q_sum.min()
            Qn_max = self.Q_sum.max()
            print('  Q_sum_min, Q_sum_max =', Qn_min, ', ', Qn_max, '[W m-2]')
            print()
                        
    #   update_net_energy_flux()   
    #-------------------------------------------------------------------  
    def open_input_files(self):

        #------------------------------------------------------------
        # Note: 2020-05-03. Changed in_directory to met_directory,
        #       which now defaults to cfg_directory.
        #       See set_directories() in BMI_base.py.
        #------------------------------------------------------------
        if (self.DEBUG):
            print('Calling open_input_files()...')
 
        #-------------------------------------------------------      
        # Note: slope and aspect files were read already from
        #       topo_directory, in initialize_computed_vars().
        #-------------------------------------------------------
        self.P_file      = self.met_directory + self.P_file
        self.T_air_file  = self.met_directory + self.T_air_file
        self.T_surf_file = self.met_directory + self.T_surf_file
        self.RH_file     = self.met_directory + self.RH_file
        self.p0_file     = self.met_directory + self.p0_file
        self.uz_file     = self.met_directory + self.uz_file
        self.z_file      = self.met_directory + self.z_file
        self.z0_air_file = self.met_directory + self.z0_air_file

        self.albedo_file        = self.met_directory + self.albedo_file
        self.em_surf_file       = self.met_directory + self.em_surf_file
        self.dust_atten_file    = self.met_directory + self.dust_atten_file
        self.cloud_factor_file  = self.met_directory + self.cloud_factor_file
        self.canopy_factor_file = self.met_directory + self.canopy_factor_file

        self.P_unit      = model_input.open_file(self.P_type,      self.P_file,
                                 NGEN_CSV=self.NGEN_CSV)
        self.T_air_unit  = model_input.open_file(self.T_air_type,  self.T_air_file)
        self.T_surf_unit = model_input.open_file(self.T_surf_type, self.T_surf_file)
        self.RH_unit     = model_input.open_file(self.RH_type,     self.RH_file)
        self.p0_unit     = model_input.open_file(self.p0_type,     self.p0_file)
        self.uz_unit     = model_input.open_file(self.uz_type,     self.uz_file)
        self.z_unit      = model_input.open_file(self.z_type,      self.z_file)
        self.z0_air_unit = model_input.open_file(self.z0_air_type, self.z0_air_file)
               
        #-----------------------------------------------
        # These are needed to compute Qn_SW and Qn_LW.
        #-----------------------------------------------
        self.albedo_unit        = model_input.open_file(self.albedo_type,
                                                        self.albedo_file)
        self.em_surf_unit       = model_input.open_file(self.em_surf_type,
                                                        self.em_surf_file)
        self.dust_atten_unit    = model_input.open_file(self.dust_atten_type,
                                                        self.dust_atten_file)
        self.cloud_factor_unit  = model_input.open_file(self.cloud_factor_type,
                                                        self.cloud_factor_file)
        self.canopy_factor_unit = model_input.open_file(self.canopy_factor_type,
                                                        self.canopy_factor_file)
        #----------------------------------------------------------------------------
        # Note: GMT_offset plus slope and aspect grids will be read separately.
        #----------------------------------------------------------------------------

##        self.Qn_SW_unit  = model_input.open_file(self.Qn_SW_type,  self.Qn_SW_file)
##        self.Qn_LW_unit  = model_input.open_file(self.Qn_LW_type,  self.Qn_LW_file)
        
    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        if (self.DEBUG):
            print('Calling read_input_files()...')

        rti = self.rti
           
        #--------------------------------------------------------
        # All grids are assumed to have a data type of float32.
        #--------------------------------------------------------
        # NB! read_next() returns None if TYPE arg is "Scalar".
        #--------------------------------------------------------
        # 2022-11-18.  Added NGEN_CSV option in model_input.py.
        #--------------------------------------------------------
        if (self.NGEN_CSV):
            # units_factor = self.mph_to_mps
            P_units_factor = self.mmps_to_mps
        else:
            P_units_factor = self.mmph_to_mps
        P = model_input.read_next(self.P_unit, self.P_type, rti,
                                  units_factor=P_units_factor,
                                  NGEN_CSV=self.NGEN_CSV)
                                                  
        # print('MET: (time,P) =', self.time, P)
        # print('shape(P) =', P.shape)
                    
        if (P is None):
            ## print('read_next() returned P is None.')
            self.set_rain_to_zero()  # (2020-05-05)
        else:
            #---------------------------------------------
            # Added 3rd "factor" argument on 2022-02-15.
            #---------------------------------------------
            ## print('read_next() returned P is NOT None.')
            ## print('P.min, P.max =', P.min(), P.max() )   
            self.update_var( 'P', P, self.P_factor )  # (In BMI_base.py)

            if not(self.SILENT):
                if (self.DEBUG or (self.time_index == 0)):
                    print('In met_base read_input_files():')
                    print('   time = ' + str(self.time) )
                    Pmin_str = str( P.min() * self.mps_to_mmph )
                    Pmax_str = str( P.max() * self.mps_to_mmph )
                    print('   min(P) = ' + Pmin_str + ' [mmph]')
                    print('   max(P) = ' + Pmax_str + ' [mmph]')
                    print(' ')

        ###############################################################
        # If any of these are scalars (read from a time series file)
        # then we'll need to use "fill()" method to prevent breaking
        # the reference to the "mutable scalar". (2/7/13)
        ###############################################################
        T_air = model_input.read_next(self.T_air_unit, self.T_air_type, rti)
        ## print('## T_air_type =', self.T_air_type )
        if (T_air is not None):
            # For testing
            # print('min(T_air) =', T_air.min() )
            # print('max(T_air) =', T_air.max() )
            # print('T_air.dtype =', T_air.dtype )
            # print('T_air.shape =', T_air.shape )
            # print('Assuming air temperature units are Kelvin.')
            #------------------------------------------------
            # 2023-08-28.  This is a bug.  Already Celsius.
            #------------------------------------------------
#             CONVERT_K_TO_C = True   ###############
#             if (CONVERT_K_TO_C):
#                 T_air -= 273.15
            self.update_var( 'T_air', T_air )

        #----------------------------------------------------------
        # 2023-09-01.  Note that in CFG files, boolean flags are
        # set to "Yes" or "No".  However, both strings evaluate
        # to True in Python, so need to make sure it is correctly
        # remapped to True of False.  Modified read_config_file()
        # in BMI_base.py to always do this for any CFG var whose
        # value is set to Yes/True/On or No/False/Off.
        #----------------------------------------------------------
        # 2023-09-05.  Moved this down so we also read T_air for
        # the degree-day snow component.
        #----------------------------------------------------------         
#         if (self.DEBUG):
#            print('In read_input_files(), self.PRECIP_ONLY =', self.PRECIP_ONLY)         
        if (self.PRECIP_ONLY):  # 2022-11-28, for speed.
            return
        ###########################################

        T_surf = model_input.read_next(self.T_surf_unit, self.T_surf_type, rti)
        if (T_surf is not None):
            self.update_var( 'T_surf', T_surf )

        RH = model_input.read_next(self.RH_unit, self.RH_type, rti)
        if (RH is not None):
            self.update_var( 'RH', RH )

        p0 = model_input.read_next(self.p0_unit, self.p0_type, rti)
        if (p0 is not None):
            self.update_var( 'p0', p0 )

        uz = model_input.read_next(self.uz_unit, self.uz_type, rti)
        if (uz is not None):
            self.update_var( 'uz', uz )

        z = model_input.read_next(self.z_unit, self.z_type, rti)
        if (z is not None):
            self.update_var( 'z', z )

        z0_air = model_input.read_next(self.z0_air_unit, self.z0_air_type, rti)
        if (z0_air is not None):
            self.update_var( 'z0_air', z0_air )

        #----------------------------------------------------------------------------
        # These are needed to compute Qn_SW and Qn_LW.
        #----------------------------------------------------------------------------
        # Note: We could later write a version of read_next() that takes "self"
        #       and "var_name" as args and that uses "exec()".
        #----------------------------------------------------------------------------
        albedo = model_input.read_next(self.albedo_unit, self.albedo_type, rti)
        if (albedo is not None):
            self.update_var( 'albedo', albedo )
            ## self.albedo = albedo

        em_surf = model_input.read_next(self.em_surf_unit, self.em_surf_type, rti)
        if (em_surf is not None):
            self.update_var( 'em_surf', em_surf )
            ## self.em_surf = em_surf

        dust_atten = model_input.read_next(self.dust_atten_unit, self.dust_atten_type, rti)
        if (dust_atten is not None):
            self.update_var( 'dust_atten', dust_atten )
            ## self.dust_atten = dust_atten

        cloud_factor = model_input.read_next(self.cloud_factor_unit, self.cloud_factor_type, rti)
        if (cloud_factor is not None):
            self.update_var( 'cloud_factor', cloud_factor )
            ## self.cloud_factor = cloud_factor

        canopy_factor = model_input.read_next(self.canopy_factor_unit, self.canopy_factor_type, rti)
        if (canopy_factor is not None):
            self.update_var( 'canopy_factor', canopy_factor )
            ## self.canopy_factor = canopy_factor

        #------------------------------------------------------------
        # EXPERIMENTAL, NOT FINISHED
        # Read variables from files into scalars or grids while
        # again making sure to preserve references (in-place),
        # but all in one function call. (11/15/16)
        #------------------------------------------------------------
#         model_input.read_next2(self, 'T_air',   rti)
#         model_input.read_next2(self, 'T_surf',  rti)
#         model_input.read_next2(self, 'RH',      rti)
#         model_input.read_next2(self, 'p0',      rti)
#         model_input.read_next2(self, 'uz',      rti)
#         model_input.read_next2(self, 'z',       rti)
#         model_input.read_next2(self, 'z0_air',  rti)
#         #----------------------------------------------------
#         model_input.read_next2(self, 'albedo',  rti)
#         model_input.read_next2(self, 'em_surf', rti)
#         model_input.read_next2(self, 'dust_atten',    rti)
#         model_input.read_next2(self, 'cloud_factor',  rti)
#         model_input.read_next2(self, 'canopy_factor', rti)

        #-------------------------------------------------------------
        # Compute Qsw_prefactor from cloud_factor and canopy factor.
        #-------------------------------------------------------------
        ## self.Qsw_prefactor = 
        
        #-----------------------------------------------------------
        # These are now computed by two functions in this file:
        # "update_net_shortwave_radiation()" and
        # "update_net_longwave_radiation()", called from update().
        # This could be used to read them from a file.
        #-----------------------------------------------------------        
##        Qn_SW = model_input.read_next(self.Qn_SW_unit, self.Qn_SW_type, rti)
##        if (Qn_SW is not None):
##            self.update_var( 'Qn_SW', Qn_SW )
##
##        Qn_LW = model_input.read_next(self.Qn_LW_unit, self.Qn_LW_type, rti)
##        if (Qn_LW is not None):
##            self.update_var( 'Qn_LW', Qn_LW )
         
    #   read_input_files()
    #-------------------------------------------------------------------  
    def set_rain_to_zero(self):

        #-------------------------------------------------------     
        # Note: As a variable, rain rate, P, is special.
        #       When there are no more data values left in an
        #       input file, we can assume that P=0 after that.
        #       Other variables may persist the last value
        #       that was read, but not P.
        #-------------------------------------------------------
        #       Uniform rain for a given duration can be
        #       specified by setting dt to duration in the
        #       meteorology CFG file, and then setting the
        #       P value to a single scalar or grid.
        #-------------------------------------------------------
        #       Using fill() method writes values in-place and
        #       does not "break the reference".  Note that
        #       scalar P is saved as a 0D numpy array.
        #-------------------------------------------------------
        #       2/7/13. Note that we don't change P from grid
        #       to scalar since that could cause trouble for
        #       other comps that use P, so we just zero it out.
        #--------------------------------------------------
        if (self.RAIN_OVER):
            return  # (do nothing, already set to zero)

        # Note: read_next() returns None for scalar type.
        SCALAR_P = (self.P_type.lower() == 'scalar')
        # if (SCALAR_P and (self.time == 0)):  ## EQUIVALENT HERE
        if (SCALAR_P and (self.time_sec < self.dt)):
            return 

        self.RAIN_OVER = True
        self.P.fill( 0 )  # (works for scalar or grid)
        
        #-----------------------------------------------
        # Either self.P_type is "Scalar" or we've read
        # all of the data in the rain_rates file.
        #-----------------------------------------------
        REPORT = True
        if not(REPORT):
            return
      
        if not(self.SILENT):
            if (SCALAR_P):
                print('Reached end of scalar rainfall duration.')
                print('  P set to 0 by read_input_files().')
            else:
                print('Reached end of precip file:')
                ## print( self.P_file )
                print('  P set to 0 by met_base.read_input_files().')
                print('  time =', self.time_min, 'minutes.')
            
        # print '######### In met_base.read_input_files() #######'
        # print 'self.P_type =', self.P_type
        # print 'self.P      =', self.P

    #   set_rain_to_zero()
    #-------------------------------------------------------------------  
    def close_input_files(self):

        if (self.DEBUG):
            print('Calling close_input_files()...')

        if (self.P_type      != 'Scalar'): self.P_unit.close()
        if (self.T_air_type  != 'Scalar'): self.T_air_unit.close()
        if (self.T_surf_type != 'Scalar'): self.T_surf_unit.close()
        if (self.RH_type     != 'Scalar'): self.RH_unit.close()
        if (self.p0_type     != 'Scalar'): self.p0_unit.close()
        if (self.uz_type     != 'Scalar'): self.uz_unit.close()
        if (self.z_type      != 'Scalar'): self.z_unit.close()
        if (self.z0_air_type != 'Scalar'): self.z0_air_unit.close()

        #---------------------------------------------------
        # These are needed to compute Qn_SW and Qn_LW.
        #---------------------------------------------------
        if (self.albedo_type        != 'Scalar'): self.albedo_unit.close()
        if (self.em_surf_type       != 'Scalar'): self.em_surf_unit.close()
        if (self.dust_atten_type    != 'Scalar'): self.dust_atten_unit.close()
        if (self.cloud_factor_type  != 'Scalar'): self.cloud_factor_unit.close()
        if (self.canopy_factor_type != 'Scalar'): self.canopy_factor_unit.close()
        
##        if (self.Qn_SW_type  != 'Scalar'): self.Qn_SW_unit.close()        
##        if (self.Qn_LW_type  != 'Scalar'): self.Qn_LW_unit.close()
        
##        if (self.P_file      != ''): self.P_unit.close()
##        if (self.T_air_file  != ''): self.T_air_unit.close()
##        if (self.T_surf_file != ''): self.T_surf_unit.close()
##        if (self.RH_file     != ''): self.RH_unit.close()
##        if (self.p0_file     != ''): self.p0_unit.close()
##        if (self.uz_file     != ''): self.uz_unit.close()
##        if (self.z_file      != ''): self.z_unit.close()
##        if (self.z0_air_file != ''): self.z0_air_unit.close()
##        #--------------------------------------------------------
##        if (self.Qn_SW_file  != ''): self.Qn_SW_unit.close()        
##        if (self.Qn_LW_file  != ''): self.Qn_LW_unit.close()
        
    #   close_input_files()
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        if (self.DEBUG):
            print('Calling update_outfile_names()...')
        
        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
        self.ea_gs_file  = (self.out_directory + self.ea_gs_file  )
        self.es_gs_file  = (self.out_directory + self.es_gs_file  )
        self.Qsw_gs_file = (self.out_directory + self.Qsw_gs_file )
        self.Qlw_gs_file = (self.out_directory + self.Qlw_gs_file )
        self.ema_gs_file = (self.out_directory + self.ema_gs_file )
        self.tsurf_gs_file = (self.out_directory + self.tsurf_gs_file )
        self.alb_gs_file = (self.out_directory + self.alb_gs_file)
        #------------------------------------------------------------
        self.ea_ts_file  = (self.out_directory + self.ea_ts_file  )
        self.es_ts_file  = (self.out_directory + self.es_ts_file  )
        self.Qsw_ts_file = (self.out_directory + self.Qsw_ts_file )
        self.Qlw_ts_file = (self.out_directory + self.Qlw_ts_file )
        self.ema_ts_file = (self.out_directory + self.ema_ts_file )
        self.tsurf_ts_file = (self.out_directory + self.tsurf_ts_file)
        self.alb_ts_file = (self.out_directory + self.alb_ts_file)
        
##        self.ea_gs_file = (self.case_prefix + '_2D-ea.rts')
##        self.es_gs_file = (self.case_prefix + '_2D-es.rts')
##        #-----------------------------------------------------
##        self.ea_ts_file = (self.case_prefix + '_0D-ea.txt')
##        self.es_ts_file = (self.case_prefix + '_0D-es.txt')

    #   update_outfile_names()   
    #-------------------------------------------------------------------  
    def open_output_files(self):

        if (self.DEBUG):
            print('Calling open_output_files()...')
        model_output.check_netcdf( SILENT=self.SILENT )
        self.update_outfile_names()
        
        #--------------------------------------
        # Open new files to write grid stacks
        #--------------------------------------
        if (self.SAVE_EA_GRIDS):
            model_output.open_new_gs_file( self, self.ea_gs_file, self.rti,
                                           ## var_name='e_air',
                                           var_name='ea',
                                           long_name='vapor_pressure_in_air',
                                           units_name='mbar')
            
        if (self.SAVE_ES_GRIDS):
            model_output.open_new_gs_file( self, self.es_gs_file, self.rti,
                                           ## var_name='e_surf',
                                           var_name='es',
                                           long_name='vapor_pressure_at_surface',
                                           units_name='mbar')
        if (self.SAVE_QSW_GRIDS):
            model_output.open_new_gs_file( self, self.Qsw_gs_file, self.rti,
                                           var_name='Qsw',
                                           ## var_name='Qn_SW',
                                           long_name='net_shortwave_radiation',
                                           units_name='W/m^2')
            
        if (self.SAVE_QLW_GRIDS):
            model_output.open_new_gs_file( self, self.Qlw_gs_file, self.rti,
                                           var_name='Qlw',
                                           ## var_name='Qn_LW',
                                           long_name='net_longwave_radiation',
                                           units_name='W/m^2')
        if (self.SAVE_EMA_GRIDS):
            model_output.open_new_gs_file( self, self.ema_gs_file, self.rti,
                                           var_name='ema',
                                           long_name='air_emissivity',
                                           units_name='none')
            
        if (self.SAVE_TSURF_GRIDS):
            model_output.open_new_gs_file( self, self.tsurf_gs_file, self.rti,
                                            var_name='tsurf',
                                            long_name='surface_temperature',
                                            units_name='C')  
            
        if (self.SAVE_ALB_GRIDS):
            model_output.open_new_gs_file( self, self.alb_gs_file, self.rti,
                                            var_name='alb',
                                            long_name='albedo',
                                            units_name='none') 
        #--------------------------------------
        # Open new files to write time series
        #--------------------------------------
        IDs = self.outlet_IDs 
        if (self.SAVE_EA_PIXELS):
            model_output.open_new_ts_file( self, self.ea_ts_file, IDs,
                                           ## var_name='e_air',
                                           var_name='ea',
                                           long_name='vapor_pressure_in_air',
                                           units_name='mbar')

        if (self.SAVE_ES_PIXELS):
            model_output.open_new_ts_file( self, self.es_ts_file, IDs,
                                           ## var_name='e_surf',
                                           var_name='es',
                                           long_name='vapor_pressure_at_surface',
                                           units_name='mbar')

        if (self.SAVE_QSW_PIXELS):
            model_output.open_new_ts_file( self, self.Qsw_ts_file, IDs,
                                           var_name='Qsw',
                                           long_name='net_shortwave_radiation',
                                           units_name='W/m^2')

        if (self.SAVE_QLW_PIXELS):
            model_output.open_new_ts_file( self, self.Qlw_ts_file, IDs,
                                           var_name='Qlw',
                                           long_name='net_longwave_radiation',
                                           units_name='W/m^2')
            
        if (self.SAVE_EMA_PIXELS):
            model_output.open_new_ts_file( self, self.ema_ts_file, IDs,
                                           var_name='ema',
                                           long_name='air_emissivity',
                                           units_name='none')
            
        if (self.SAVE_TSURF_PIXELS):
            model_output.open_new_ts_file( self, self.tsurf_ts_file, IDs,
                                            var_name='tsurf',
                                            long_name='surface_temperature',
                                            units_name='C')      
        if (self.SAVE_ALB_PIXELS):
            model_output.open_new_ts_file( self, self.alb_ts_file, IDs,
                                            var_name='alb',
                                            long_name='albedo',
                                            units_name='none')   

    #   open_output_files()
    #-------------------------------------------------------------------
    def write_output_files(self, time_seconds=None):

        if (self.DEBUG):
            print('Calling write_output_files()...')

        #-----------------------------------------
        # Allows time to be passed from a caller
        #-----------------------------------------
        if (time_seconds is None):
            time_seconds = self.time_sec
        model_time = int(time_seconds)
        
        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
        if (model_time % int(self.save_grid_dt) == 0):
            self.save_grids()
        if (model_time % int(self.save_pixels_dt) == 0):
            self.save_pixel_values()

    #  write_output_files()
    #-------------------------------------------------------------------
    def close_output_files(self):
    
        if (self.SAVE_EA_GRIDS):   model_output.close_gs_file( self, 'ea')
        if (self.SAVE_ES_GRIDS):   model_output.close_gs_file( self, 'es')
        if (self.SAVE_QSW_GRIDS):  model_output.close_gs_file( self, 'Qsw')
        if (self.SAVE_QLW_GRIDS):  model_output.close_gs_file( self, 'Qlw')
        if (self.SAVE_EMA_GRIDS):  model_output.close_gs_file( self, 'ema')
        if (self.SAVE_TSURF_GRIDS): model_output.close_gs_file(self, 'tsurf')
        if (self.SAVE_ALB_GRIDS): model_output.close_gs_file(self, 'alb')
        #-------------------------------------------------------------------
        if (self.SAVE_EA_PIXELS):  model_output.close_ts_file( self, 'ea') 
        if (self.SAVE_ES_PIXELS):  model_output.close_ts_file( self, 'es') 
        if (self.SAVE_QSW_PIXELS): model_output.close_ts_file( self, 'Qsw') 
        if (self.SAVE_QLW_PIXELS): model_output.close_ts_file( self, 'Qlw')
        if (self.SAVE_EMA_PIXELS): model_output.close_ts_file( self, 'ema')
        if (self.SAVE_TSURF_PIXELS): model_output.close_ts_file(self, 'tsurf')
        if (self.SAVE_ALB_PIXELS): model_output.close_ts_file(self, 'alb')


    #   close_output_files()        
    #-------------------------------------------------------------------  
    def save_grids(self):
       
        if (self.SAVE_EA_GRIDS):
            model_output.add_grid( self, self.e_air,  'ea', self.time_min )
            
        if (self.SAVE_ES_GRIDS):
            model_output.add_grid( self, self.e_surf, 'es', self.time_min )

        if (self.SAVE_QSW_GRIDS):
            model_output.add_grid( self, self.Qn_SW, 'Qsw', self.time_min )
            
        if (self.SAVE_QLW_GRIDS):
            model_output.add_grid( self, self.Qn_LW, 'Qlw', self.time_min )

        if (self.SAVE_EMA_GRIDS):
            model_output.add_grid( self, self.em_air, 'ema', self.time_min )
        
        if (self.SAVE_TSURF_GRIDS):
            model_output.add_grid( self, self.T_surf, 'tsurf', self.time_min )

        if (self.SAVE_ALB_GRIDS):
            model_output.add_grid( self, self.albedo, 'alb', self.time_min )
            
    #   save_grids()            
    #-------------------------------------------------------------------  
    def save_pixel_values(self):

        IDs  = self.outlet_IDs
        time = self.time_min   ######
        
        if (self.SAVE_EA_PIXELS):
            model_output.add_values_at_IDs( self, time, self.e_air,  'ea', IDs )
            
        if (self.SAVE_ES_PIXELS):
            model_output.add_values_at_IDs( self, time, self.e_surf, 'es', IDs )

        if (self.SAVE_QSW_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Qn_SW, 'Qsw', IDs )
            
        if (self.SAVE_QLW_PIXELS):
            model_output.add_values_at_IDs( self, time, self.Qn_LW, 'Qlw', IDs )

        if (self.SAVE_EMA_PIXELS):
            model_output.add_values_at_IDs( self, time, self.em_air, 'ema', IDs )

        if (self.SAVE_TSURF_PIXELS):
            model_output.add_values_at_IDs( self, time, self.T_surf, 'tsurf', IDs )

        if (self.SAVE_ALB_PIXELS):
            model_output.add_values_at_IDs( self, time, self.albedo, 'alb', IDs )

            
    #   save_pixel_values()
    #-------------------------------------------------------------------
#---------------------------------------------------------------------------------
def compare_em_air_methods():

    #--------------------------------------------------------------
    # Notes:  There are two different methods that are commonly
    #         used to compute the vapor pressure of air, e_air,
    #         and then the emissivity of air, em_air, for use in
    #         longwave radiation calculations.  This routine
    #         compares them graphically.
    #
    # NB!     This hasn't been tested since conversion from IDL.
    #-------------------------------------------------------------
    import matplotlib.pyplot
    
    T_air = np.arange(80, dtype='float32') - np.float64(40)   #[Celsius]  (-40 to 40)
    RH  = np.float64(1.0)
    C2K = np.float64(273.15)
    
    #--------------------------
    # Brutsaert (1975) method
    #--------------------------
    term1   = (np.float64(17.3) * T_air) / (T_air + np.float64(237.3))   ######### DOUBLE CHECK THIS (7/26/13)
    e_air1  = RH * np.float64(0.611) * np.exp( term1 )  # [kPa]
    em_air1 = np.float64(1.72) * (e_air1 / (T_air + C2K)) ** (np.float64(1) / 7)
    
    #---------------------------
    # Satterlund (1979) method
    #----------------------------
    # NB! e_air has units of Pa
    #----------------------------
    term2   = np.float64(2353) / (T_air + C2K)
    e_air2  = RH * np.float64(10) ** (np.float64(11.40) - term2)   # [Pa]
    eterm   = np.exp(-np.float64(1) * (e_air2 / np.float64(100)) ** ((T_air + C2K) / np.float64(2016)))
    em_air2 = np.float64(1.08) * (np.float64(1) - eterm)
    
    #----------------------------
    # Plot the two e_air curves
    #--------------------------------
    # These two agree quite closely
    #--------------------------------
    matplotlib.pyplot.figure(figsize=(8, 6), dpi=80)
    matplotlib.pyplot.show()
    matplotlib.pyplot.plot(T_air, e_air1)
    matplotlib.pyplot.show()
    ## oplot(T_air, (e_air2 / np.float64(1000)), psym=-3)   # [Pa -> kPa]
    
    #-----------------------------
    # Plot the two em_air curves
    #--------------------------------------------------
    # These two don't agree very well for some reason
    #--------------------------------------------------
    matplotlib.pyplot.figure(figsize=(8, 6), dpi=80)
    matplotlib.pyplot.show()
    matplotlib.pyplot.plot(T_air, em_air1)
    matplotlib.pyplot.show()
    ## oplot(T_air, em_air2, psym=-3)
    
#   compare_em_air_Methods
#---------------------------------------------------------------------------------
    
        
