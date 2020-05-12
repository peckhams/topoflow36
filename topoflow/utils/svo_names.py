"""
Cross-component (global) SVO name handling functions
"""
# Copyright (c) 2020, Scott D. Peckham

# May 2020.  Added:  SVO name for 'vol_swe'.
#            Added:  SVO name for 'vol_soil'.

# import <something> here, if needed

#---------------------------------------------------------------------
#
#   test1()
#   get_svo_name()
#   get_short_name()
#
#   get_short_name_map()     # local to TopoFlow
#   get_svo_name_map()
#
#---------------------------------------------------------------------
def test1( short_name = 'T_air' ):
  
    svo_name = get_svo_name( short_name )
    print('Short TF name =', short_name)
    print('SVO name =', svo_name)

#   test1()
#---------------------------------------------------------------------
def get_svo_name( short_name ):

    svo_name_map  = get_svo_name_map()
    svo_name      = svo_name_map[ short_name ]
    return svo_name

#   get_svo_name()
#---------------------------------------------------------------------
def get_short_name( svo_name ):

    short_name_map = get_short_name_map()
    short_name     = short_name_map[ svo_name ]
    return short_name

#   get_short_name()
#---------------------------------------------------------------------
def get_svo_name_map():

    #------------------------------------------------
    # Note:  How to create an "inverse map":
    # inv_map = dict(zip(map.values(), map.keys()))
    #------------------------------------------------
    short_name_map = get_short_name_map()
    svo_name_map   = dict( zip( short_name_map.values(),
                                short_name_map.keys() ) )

    return svo_name_map
    
#   get_svo_name_map()
#---------------------------------------------------------------------
def get_short_name_map():

    #---------------------------------------------------------------
    # Note:  Think of this as a "list of symbols" for the entire
    #        TopoFlow application (all components).
    #---------------------------------------------------------------
    # Note:  "land_surface__elevation" & "land_surface__slope"
    #        have different short names in some components, and
    #        this has not yet been taken into account. 
    #        In satzone component, "elev" is used, but "DEM" is
    #        used in the D8, Erode, and smooth_dem components.
    #        Mapping to "elev" and "S_bed" for now. (2020-01-09)
    #---------------------------------------------------------------
    # Note:  Names starting with "model_" are not official SVO
    #        standard names, but some might be used by EMELI.
    #---------------------------------------------------------------    
    short_name_map = {
    'atmosphere__optical_path_length_ratio' : 'M_opt',
    'atmosphere__von_karman_constant' : 'kappa',
    'atmosphere_aerosol_dust__reduction_of_transmittance' : 'dust_atten',
    'atmosphere_air-column_water-vapor__liquid-equivalent_depth' : 'W_p',
    'atmosphere_bottom_air__brutsaert_emissivity_canopy_factor' : 'canopy_factor',
    'atmosphere_bottom_air__brutsaert_emissivity_cloud_factor' : 'cloud_factor',
    'atmosphere_bottom_air__bulk_latent_heat_aerodynamic_conductance' : 'De',
    'atmosphere_bottom_air__bulk_sensible_heat_aerodynamic_conductance' : 'Dh',
    'atmosphere_bottom_air__emissivity' : 'em_air',
    'atmosphere_bottom_air__mass-per-volume_density' : 'rho_air',
    'atmosphere_bottom_air__mass-specific_isobaric_heat_capacity' : 'Cp_air',
    'atmosphere_bottom_air__neutral_bulk_heat_aerodynamic_conductance' : 'Dn',
    'atmosphere_bottom_air__pressure' : 'p0',
    'atmosphere_bottom_air__temperature' : 'T_air',
    'atmosphere_bottom_air_flow__bulk_richardson_number' : 'Ri',
    'atmosphere_bottom_air_flow__log_law_roughness_length' : 'z0_air',
    'atmosphere_bottom_air_flow__reference-height_speed' : 'uz',
    'atmosphere_bottom_air_flow__speed_reference_height' : 'z',
    'atmosphere_bottom_air_land_net-latent-heat__energy_flux' : 'Qe',
    'atmosphere_bottom_air_land_net-sensible-heat__energy_flux' : 'Qh',
    'atmosphere_bottom_air_water-vapor__dew_point_temperature' : 'T_dew',
    'atmosphere_bottom_air_water-vapor__partial_pressure' : 'e_air',
    'atmosphere_bottom_air_water-vapor__relative_saturation' : 'RH',
    'atmosphere_bottom_air_water-vapor__saturated_partial_pressure' : 'e_sat_air',
    'atmosphere_water__domain_time_integral_of_precipitation_leq-volume_flux' : 'vol_P',
    'atmosphere_water__domain_time_max_of_precipitation_leq-volume_flux' : 'P_max',
    'atmosphere_water__geomorphic_precipitation_leq-volume_flux' : 'R',
    'atmosphere_water__precipitation_leq-volume_flux' : 'P',
    'atmosphere_water__rainfall_volume_flux' : 'P_rain',
    'atmosphere_water__snowfall_leq-volume_flux' : 'P_snow',
    'basin_outlet_water_flow__half_of_fanning_friction_factor' : 'f_outlet',
    'basin_outlet_water_x-section__mean_depth' : 'd_outlet',
    'basin_outlet_water_x-section__peak_time_of_depth' : 'Td_peak',
    'basin_outlet_water_x-section__peak_time_of_volume_flow_rate' : 'T_peak',
    'basin_outlet_water_x-section__peak_time_of_volume_flux' : 'Tu_peak',
    'basin_outlet_water_x-section__time_integral_of_volume_flow_rate' : 'vol_Q',
    'basin_outlet_water_x-section__time_max_of_mean_depth' : 'd_peak',
    'basin_outlet_water_x-section__time_max_of_volume_flow_rate' : 'Q_peak',
    'basin_outlet_water_x-section__time_max_of_volume_flux' : 'u_peak',
    'basin_outlet_water_x-section__volume_flow_rate' : 'Q_outlet',
    'basin_outlet_water_x-section__volume_flux' : 'u_outlet',
    'bedrock__uplift_rate' : 'U',
    'canals__count' : 'n_canals',
    'canals_entrance__x_coordinate' : 'canals_in_x',
    'canals_entrance__y_coordinate' : 'canals_in_y',
    'canals_entrance_water__volume_flow_rate' : 'Q_canals_in',
    'canals_entrance_water__volume_fraction' : 'Q_canals_fraction',
    'canals_exit__x_coordinate' : 'canals_out_x',
    'canals_exit__y_coordinate' : 'canals_out_y',
    'canals_exit_water__volume_flow_rate' : 'Q_canals_out',
    'channel_bottom_surface__slope' : 'S_bed',
    'channel_bottom_water_flow__domain_max_of_log_law_roughness_length' : 'z0val_max',
    'channel_bottom_water_flow__domain_min_of_log_law_roughness_length' : 'z0val_min',
    'channel_bottom_water_flow__log_law_roughness_length' : 'z0val',
    'channel_bottom_water_flow__magnitude_of_shear_stress' : 'tau',
    'channel_bottom_water_flow__shear_speed' : 'u_star',
    'channel_centerline__sinuosity' : 'sinu',
    'channel_water__volume' : 'vol',
    'channel_water_flow__domain_max_of_manning_n_parameter' : 'nval_max',
    'channel_water_flow__domain_min_of_manning_n_parameter' : 'nval_min',
    'channel_water_flow__froude_number' : 'froude',
    'channel_water_flow__half_of_fanning_friction_factor' : 'f',
    'channel_water_flow__manning_n_parameter' : 'nval',
    'channel_water_surface__slope' : 'S_free',
    'channel_water_total-sediment__volume_flow_rate' : 'Qs',
    'channel_water_total-sediment__volume_flow_rate_law_area_exponent' : 'm',
    'channel_water_total-sediment__volume_flow_rate_law_coefficient' : 'K',
    'channel_water_total-sediment__volume_flow_rate_law_slope_exponent' : 'n',
    'channel_water_x-section__domain_max_of_mean_depth' : 'd_max',
    'channel_water_x-section__domain_max_of_volume_flow_rate' : 'Q_max',
    'channel_water_x-section__domain_max_of_volume_flux' : 'u_max',
    'channel_water_x-section__domain_min_of_mean_depth' : 'd_min',
    'channel_water_x-section__domain_min_of_volume_flow_rate' : 'Q_min',
    'channel_water_x-section__domain_min_of_volume_flux' : 'u_min',
    'channel_water_x-section__hydraulic_radius' : 'Rh',
    'channel_water_x-section__initial_mean_depth' : 'd0',
    'channel_water_x-section__mean_depth' : 'd',
    'channel_water_x-section__volume_flow_rate' : 'Q',
    'channel_water_x-section__volume_flow_rate_law_area_exponent' : 'p',
    'channel_water_x-section__volume_flux' : 'u',
    'channel_water_x-section__wetted_area' : 'A_wet',
    'channel_water_x-section__wetted_perimeter' : 'P_wet',
    'channel_x-section_trapezoid_bottom__width' : 'width',
    'channel_x-section_trapezoid_side__flare_angle' : 'angle',
    'earth__standard_gravity_constant' : 'g',
    'glacier_ice__domain_time_integral_of_melt_volume_flux' : 'vol_MR',
    'glacier_ice__melt_volume_flux' : 'MR',
    'glacier_ice__thickness' : 'H',
    'glacier_top_surface__elevation' : 'Zi',
    'land_surface__albedo' : 'albedo',
    'land_surface__aspect_angle' : 'alpha',
    'land_surface__domain_max_of_increment_of_elevation' : 'dz_max',
    #'land_surface__elevation' : 'DEM',   # used by d8, erode, smooth_dem
    #'land_surface__elevation' : 'Zb',    # used by gc2d, ice_base
    'land_surface__elevation' : 'elev',   # used by infil, satzone
    'land_surface__emissivity' : 'em_surf',
    'land_surface__increment_of_elevation' : 'dz',
    'land_surface__initial_elevation' : 'z0',
    'land_surface__latitude' : 'lat_deg',
    'land_surface__longitude' : 'lon_deg',
    'land_surface__slope' : 'S_bed',
    #'land_surface__slope' : 'S',   # used by d8, erode, etc.
    'land_surface__slope_angle' : 'beta',
    'land_surface__temperature' : 'T_surf',
    'land_surface__time_derivative_of_elevation' : 'dz_dt',
    'land_surface_air_water-vapor__partial_pressure' : 'e_surf',
    'land_surface_air_water-vapor__saturated_partial_pressure' : 'e_sat_surf',
    'land_surface_contour-segment__total_contributing_area' : 'A',
    'land_surface_net-longwave-radiation__energy_flux' : 'Qn_LW',
    'land_surface_net-shortwave-radiation__energy_flux' : 'Qn_SW',
    'land_surface_net-total-energy__energy_flux' : 'Q_sum',
    'land_surface_soil__conduction_heat_flux' : 'Qc',
    'land_surface_water__area_integral_of_depth' : 'vol_land',
    'land_surface_water__baseflow_volume_flux' : 'GW',    
    'land_surface_water__depth' : 'd_flood',
    'land_surface_water__domain_time_integral_of_baseflow_volume_flux' : 'vol_GW',
    'land_surface_water__domain_time_integral_of_evaporation_volume_flux' : 'vol_ET',
    'land_surface_water__domain_time_integral_of_runoff_volume_flux' : 'vol_R',
    'land_surface_water__evaporation_volume_flux' : 'ET',
    'land_surface_water__priestley-taylor_alpha_coefficient' : 'alpha',
    'land_surface_water__runoff_volume_flux' : 'R',
    'model__time_step' : 'dt',
    'model_domain_boundary__lowering_rate' : 'BLR',
    'model_grid_cell__area' : 'da',
    'model_grid_cell__d8_flow_length' : 'ds',
    'model_grid_cell__d8_flow_width' : 'dw',
    'model_grid_cell__diameter' : 'dd',
    'model_grid_cell__x_length' : 'dx',
    'model_grid_cell__y_length' : 'dy',
    'model_soil_layer-0__porosity' : 'qs_layer_0',
    'model_soil_layer-0__saturated_thickness' : 'y_layer_0',
    'model_soil_layer-0__thickness' : 'th_layer_0',
    'model_soil_layer-1__porosity' : 'qs_layer_1',
    'model_soil_layer-1__saturated_thickness' : 'y_layer_1',
    'model_soil_layer-1__thickness' : 'th_layer_1',
    'model_soil_layer-2__porosity' : 'qs_layer_2',
    'model_soil_layer-2__saturated_thickness' : 'y_layer_2',
    'model_soil_layer-2__thickness' : 'th_layer_2',
    'network_channel_water__volume' : 'vol_chan',
    'physics__stefan_boltzmann_constant' : 'sigma',
    'physics__von_karman_constant' : 'kappa',
    'sinks__count' : 'n_sinks',
    'sinks__x_coordinate' : 'sinks_x',
    'sinks__y_coordinate' : 'sinks_y',
    'sinks_water__volume_flow_rate' : 'Q_sinks',
    'snowpack__degree-day_coefficient' : 'c0',
    'snowpack__degree-day_threshold_temperature' : 'T0',
    'snowpack__depth' : 'h_snow',
    'snowpack__domain_time_integral_of_liquid-equivalent_depth' : 'vol_swe',
    'snowpack__domain_time_integral_of_melt_volume_flux' : 'vol_SM',
    'snowpack__energy-per-area_cold_content' : 'Ecc',
    'snowpack__initial_depth' : 'h0_snow',
    'snowpack__initial_liquid-equivalent_depth' : 'h0_swe',
    'snowpack__liquid-equivalent_depth' : 'h_swe',
    'snowpack__melt_volume_flux' : 'SM',
    'snowpack__z_mean_of_mass-per-volume_density' : 'rho_snow',
    'snowpack__z_mean_of_mass-specific_isobaric_heat_capacity' : 'Cp_snow',
    'soil__porosity' : 'phi',
    'soil__reference_depth_temperature' : 'T_soil_x',   # x -> z, more clear?
    'soil__temperature_reference_depth' : 'soil_x',
    'soil__thermal_conductivity' : 'K_soil',
    'soil_surface__temperature' : 'T_surf',
    'soil_surface_water__domain_time_integral_of_infiltration_volume_flux' : 'vol_IN',
    'soil_surface_water__infiltration_volume_flux' : 'IN',
    'soil_surface_water__time_integral_of_infiltration_volume_flux' : 'I',
    'soil_surface_water__volume_fraction' : 'q0',
    'soil_top-layer__porosity' : 'qs_layer_0',
    'soil_top-layer__saturated_thickness' : 'y_layer_0',
    'soil_top-layer__thickness' : 'th_layer_0',
    'soil_water__brooks-corey_b_parameter' : 'b',
    'soil_water__brooks-corey_eta_parameter' : 'eta',
    'soil_water__brooks-corey_lambda_parameter' : 'lam',
    'soil_water__brooks-corey-smith_c_parameter' : 'c',
    'soil_water__brooks-corey-smith_pressure_head_offset_parameter' : 'pA',
    'soil_water__bubbling_pressure_head' : 'pB',
    'soil_water__domain_time_integral_of_volume_fraction' : 'vol_soil',
    'soil_water__field-capacity_volume_fraction' : 'qf',
    'soil_water__green-ampt_capillary_length' : 'G',
    'soil_water__hydraulic_conductivity' : 'K',
    'soil_water__hygroscopic_volume_fraction' : 'qH',
    'soil_water__initial_hydraulic_conductivity' : 'Ki',
    'soil_water__initial_volume_fraction' : 'qi',
    'soil_water__normalized_volume_fraction' : 'S_eff',
    'soil_water__potential_infiltration_volume_flux' : 'fc',
    'soil_water__pressure_head' : 'p',
    'soil_water__relative_hydraulic_conductivity' : 'K_rel',
    'soil_water__residual_volume_fraction' : 'qr',
    'soil_water__saturated_hydraulic_conductivity' : 'Ks',
    'soil_water__saturated_volume_fraction' : 'qs',
    'soil_water__smith-parlange_gamma_parameter' : 'gam',
    'soil_water__volume_fraction' : 'q',
    'soil_water__wilting-point_volume_fraction' : 'qw',
    'soil_water_flow__z_component_of_darcy_velocity' : 'v',
    'soil_water_sat-zone_top__domain_time_integral_of_recharge_volume_flux' : 'vol_Rg',
    'soil_water_sat-zone_top__recharge_volume_flux' : 'Rg',
    'soil_water_sat-zone_top_surface__elevation' : 'h_table',
    'soil_water_wetting-front__depth' : 'Zw',
    'sources__count' : 'n_sources',
    'sources__x_coordinate' : 'sources_x',
    'sources__y_coordinate' : 'sources_y',
    'sources_water__volume_flow_rate' : 'Q_sources',
    'water__mass-specific_latent_fusion_heat' : 'Lf',
    'water__mass-specific_latent_vaporization_heat' : 'Lv',
    'water-liquid__mass-per-volume_density' : 'rho_H2O' }

    return short_name_map

#   get_short_name_map()
#---------------------------------------------------------------------




