
S.D. Peckham
May 13, 2020

Compare outputs in logfile and log window for these:
    1. __Richards_1layer_30mmph_75min_rain_scalar
    2. __Richards_1layer_30mmph_75min_rain_series
    3. __Richards_1layer_30mmph_75min_rain_grid
    4. __Richards_1layer_30mmph_75min_rain_stack

Channel Flow:  Kinematic Wave method
Meteorology:   Meteorology (main)
Infiltration:  Richards 1D method
Evaporation:   None
Snowmelt:      None
Icemelt:       None
Diversions:    None

Key TopoFlow Driver Settings:
    dt = 6.0 sec
    stop_method  = Until_model_time
    T_stop_model = 100 (minutes)

Key Channel Flow Settings:
    dt             = 6.0 sec
    d0_type        = Scalar
    qs_list[0]     = 0.485    (saturated water content)
    qi_list[0]     = 0.4      (initial water content, relatively high)
    d0             = 0.0
    FLOOD_OPTION   = 0
    save_grid_dt   = 60.0 sec
    save_pixesl_dt = 60.0 sec
    SAVE_Q_GRIDS   = Yes
    SAVE_Q_PIXELS  = Yes

Key Meteorology Settings:
    dt = 60.0 sec    (must match time interval in time series file)
    PRECIP_ONLY = Yes
    P_type = Time_Series
    P = [case_prefix]_rain_uniform_30.txt    #######

Key Infiltration Settings:
    dt = 0.5 sec
    n_layers = 1  (parameters are for a "silt loam" soil type)
    save_grid_dt   = 60.0 sec
    save_pixesl_dt = 60.0 sec
    SAVE_Q_PROFILES  = Yes  (save soil moisture z profiles for monitored pixels)
    q_ps_file        = [case_prefix]_1D-q.nc
    SAVE_P_PROFILES  = Yes  (save pressure head z-profiles)
    p_ps_file        = [case_prefix]_1D-p.nc
    SAVE_K_PROFILES  = Yes  (save hydraul. conductivity z-profiles)
    K_ps_file        = [case_prefix]_1D-K.nc
    SAVE_V_PROFILES  = Yes  (save Darcy flow velocity z-profiles)
    v_ps_file        = case_prefix]_1D-v.nc


  