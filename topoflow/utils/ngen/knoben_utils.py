#
#  Copyright, 2025, Scott D. Peckham
#
#  knoben_utils.py
#
#  These functions allow Knoben indicators to be computed from
#  USGS GAGES-II hydroclimate data for 9067 GAGES-II CONUS basins.
#  GAGES-II has montly temperature and precip data, but does not
#  provide monthly PET data needed to compute Knoben indicators.
#  Here, the Harmon method is used to compute monthly PET from
#  the other data in GAGES-II.
#
#  calc_knoben_indicators() is called in gages2_utils.py.
#
#  See Berghuijs et al. (2014), Knoben et al. (2018), Mai et al. (2022).
#
#------------------------------------------------------------------------
#  To get started:
#
#  % conda activate tf36
#  % cd <path-to-knoben_utils.py>
#  % python
#  >>> import knoben_utils as ku
#
#------------------------------------------------------------------------
#
#  get_day_angle()
#  get_declination()
#  get_total_PET_for_day()   ## via Hamon 1963 method
#
#  get_days_per_month()
#  get_total_PET_for_month()          # sum daily pet for all days in month
#  get_harmon_monthly_PET_array()
#
#  get_knoben_monthly_moisture_index()
#  get_knoben_monthly_MI_array()
#  get_annual_aridity()          # standard definition, as a check
#  get_knoben_annual_aridity()
#  get_knoben_annual_seasonality()
#  get_knoben_annual_snow_fraction()
#  calc_knoben_indicators()
#
#  plot_monthly_precip_values()  # not ready
#  plot_monthly_temp_values()    # not ready
#
#------------------------------------------------------------------------

import numpy as np
## import matplotlib.pyplot as plt

from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import gages2_utils as g2

#------------------------------------------------------------------------
def get_day_angle( Julian_day, DEGREES=False ):

    #---------------------------------------------------------
    # Notes:  The Julian day does not need to be an integer;
    #         decimal values can be used for more precision.
    #---------------------------------------------------------

    #-------------------------------------    
    # Use this if Julian Day starts at 1
    #-------------------------------------
    ## angle = (2 * np.pi) * (Julian_day - np.float64(1)) / np.float64(365)

    #-------------------------------------    
    # Use this if Julian Day starts at 0
    #-------------------------------------
    # Don't use Days_Per_Year() here.
    #-----------------------------------------------------------
    # We should be using 366 vs. 365 for leap years, but would
    # then need to pass year to every Day_Angle() call.
    #-----------------------------------------------------------
    angle = (2 * np.pi) * Julian_day / np.float64(365)
            
    if (DEGREES):    
        angle = angle * (np.float64(180) / np.pi)
    
    return angle
    
#   get_day_angle()
#------------------------------------------------------------------------
def get_declination( day_angle, DEGREES=False, DMS=False ):

    ########################################################
    # NB! Make sure that DEGREES and DMS default to False.
    ########################################################

    #-----------------------------------------------------------
    # Note:  The declination reaches its lowest value of -23.5
    #        degrees on the Winter Solstice (Dec. 21/22) and
    #        reaches its highest value of 23.5 degrees on the
    #        Summer Solstice (June 21/22).  It is zero for
    #        both the Vernal Equinox (Mar. 20/21) and the
    #        Autumnal Equinox (Sept. 22/23).  The value of
    #        23.4397 degrees is the fixed tilt angle of the
    #        Earth's axis from from the plane of the ecliptic.
    #-----------------------------------------------------------  
    delta = np.float64(0.006918) - \
            (np.float64(0.399912) * np.cos(day_angle)) + \
            (np.float64(0.070257) * np.sin(day_angle)) - \
            (np.float64(0.006758) * np.cos(np.float64(2) * day_angle)) + \
            (np.float64(0.000907) * np.sin(np.float64(2) * day_angle)) - \
            (np.float64(0.002697) * np.cos(np.float64(3) * day_angle)) + \
            (np.float64(0.001480) * np.sin(np.float64(3) * day_angle))
    
    #------------------------------------
    # Convert from radians to degrees ?
    #------------------------------------
    if (DEGREES):    
        delta = delta * (np.float64(180) / np.pi)
    
    #----------------------------------------
    # Convert from radians to "decimal DMS"
    #----------------------------------------
    if (DMS):    
        delta = delta * (np.float64(180) / np.pi)  # [decimal degrees]
        deg = np.int16(delta)
        min = np.int16((delta - deg) * np.float64(60))
        sec = np.int16(((delta - deg) * np.float64(60) - min) * np.float64(60))
        delta = deg + (min / np.float64(100)) + (sec / np.float64(10000))   # [decimal DMS, DD.MMSS]
    
    return delta
    
#   get_declination()
#---------------------------------------------------------------------
def get_total_PET_for_day( julian_day, T_avg, lat_deg ):

    #-----------------------------------------------------------------
    # Note: Compute the daily-average PET (potential evaporation)
    #       using the Hamon 1963 method.  
    #           PET = c * (n_daylight_hours/12) * P_sat
    # 
    #       See:
    #       https://www.hec.usace.army.mil/confluence/hmsdocs/
    #       hmstrm/evaporation-and-transpiration/hamon-method
    #
    #       e_sat = saturation vapor pressure at T_avg
    #       P_sat = saturation vapor density at T_avg
    #
    #       Source of formulas:
    #       PET (Hamon, 1963)
    #       solar_dec (solar declination)
    #       n_daylight_hours (Allen et al., 1998)
    #       sunset_hour_angle (Allen et al., 1998)
    #       e_sat (Allen et al., 1998)
    #       P_sat (Wiederhold, 1997)
    #-----------------------------------------------------------------
    # Allen, R.G., L.S. Pereira, D. Raes, M. Smith (1998)
    # Crop evapotranspiration-guidelines for computing crop water 
    # requirements - FAO irrigation and drainage
    #-----------------------------------------------------------------
    day_angle = get_day_angle(julian_day)
    solar_dec = get_declination(day_angle)
    lat_rad   = lat_deg * (np.pi / 180)
    tan_product = np.tan(lat_rad) * np.tan(solar_dec)
    sunset_hour_angle = np.arccos(-1 * tan_product)
    n_daylight_hours  = (24 / np.pi) * sunset_hour_angle
    #-----------------------------------------------------
    # This formula for e_sat is due to Brutsaert and the
    # units are kPa.  For more info, see the function:
    # Saturation_Vapor_Pressure() in solar_funcs.py
    #-----------------------------------------------------
    # Pa = Pascals = unit of pressure = N / m2.
    # N  = Newton  = unit of force    = kg * m / s2
    #-----------------------------------------------------
    # The denominator in P_sat formula is a conversion
    # from degrees Celsius to Kelvin.  So units of P_sat
    # are (kPa/K).
    #-----------------------------------------------------
    # A typical value of e_sat at 25 degC = 3.165 kPa
    # A typical value of P_sat at 20 degC = 17.2 g/m3
    
    #  kPa        Pa           kg * m
    #  ---- = ---------- =  --------------
    #   K     1000 * K      s2 * 1000 * K 
    #-----------------------------------------------------    
    e_sat = 0.6108 * np.exp( 17.27 * T_avg / (T_avg + 237.3)) # [kPa]

    #--------------------------------------------------    
    # Units of P_sat should be either g/m3 or kg/m3.
    # The factor of 216.7 must also cancel the Kelvin 
    # units in the denominator.
    #--------------------------------------------------
    P_sat = 216.7 * e_sat / (T_avg + 273.16)

    #-----------------------------------------------------
    # Hamon constant should be calibrated and validated.
    # The default value of the Hamon constant used in
    # HEC HMS (see URL above) are given as:
    #    0.0065 [in/(g/m3)]
    #    0.1651 [mm/(g/m3)] = 0.01651 [cm/(g/m3)]
    # Confirmed that:  0.0065 inches = 0.1651 mm.    
    #-----------------------------------------------------
    hamon_constant = 0.1651
    PET = hamon_constant * (n_daylight_hours / 12) * P_sat

    #-------------------------------------------------------
    # Note: Using 0.1651 gave results that are consistent
    #       with the annual total PET value from GAGES-II
    #       and the aridity values.
    #       So actual units of PET here must be cm.
    #-------------------------------------------------------
    return PET

#   get_total_PET_for_day()
#---------------------------------------------------------------------
def get_days_per_month( month_num ):

    day_map = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    day_map = np.array(day_map)
    n_days  = day_map[ month_num ]
    return n_days

#   get_days_per_month()
#---------------------------------------------------------------------
def get_total_PET_for_month( month_num, T_avg, lat_deg ):

    #------------------------------------------------
    # Note: Get estimate of average PET for a given
    #       month_num that ranges from 1 to 12.
    #       PET units here are centimeters.
    #------------------------------------------------
    day_map  = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    day_map  = np.array(day_map)
    start_days = np.zeros(13, dtype='int16')
    start_days[1:] = day_map.cumsum()
    
    PET_sum = np.float32(0)
    day1 = start_days[month_num - 1]
    dayz = start_days[month_num]
    for julian_day in range(day1, dayz):
        PET_sum += get_total_PET_for_day( julian_day, T_avg, lat_deg )
    PET_month = PET_sum  # [cm]

    #------------------------------    
    # This would give the average
    #------------------------------
    ### PET_month = (PET_sum / day_map[month_num-1])

    return PET_month

#   get_total_PET_for_month()
#----------------------------------------------------------------------
def get_harmon_monthly_PET_array( T_monthly_arr, lat_deg):

    PET_monthly_arr = np.zeros(12, dtype='float32')

    for month_num in range(1,13):
        T_mon   = T_monthly_arr[ month_num-1 ]
        PET_mon = get_total_PET_for_month( month_num, T_mon, lat_deg )
        PET_monthly_arr[ month_num-1] = PET_mon
        
    return PET_monthly_arr

#   get_harmon_monthly_PET_array()
 #---------------------------------------------------------------------
def get_knoben_monthly_moisture_index(P_month, PET_month):

    #------------------------------------------------------------
    # Note: MI is a version of Thornthwaite's "moisture index",
    #       computed for a given month and latitude.
    #       See Knoben et al. (2018)
    #------------------------------------------------------------
    ## print('## P_month, PET_month [cm] =', P_month, PET_month)
    ## print('## standard aridity = ', PET_month / P_month)
    if (PET_month < P_month):
        MI = 1 - (PET_month / P_month)
    elif (PET_month == P_month):
        MI = 0
    else:
        MI = (P_month / PET_month) - 1
    return MI

#   get_knoben_monthly_moisture_index()
#---------------------------------------------------------------------
def get_knoben_monthly_MI_array( P_monthly_arr, PET_monthly_arr ):

    #---------------------------------------------    
    # Compute and save average MI for each month
    #---------------------------------------------
    MI_monthly_arr = np.zeros(12, dtype='float32')
    for k in range(12):
        P_month   = P_monthly_arr[k]
        PET_month = PET_monthly_arr[k]
        #-------------------------------
        MI = get_knoben_monthly_moisture_index(P_month, PET_month)
        MI_monthly_arr[k] = MI
    
    return MI_monthly_arr

#   get_knoben_MI_monthly_array()
#---------------------------------------------------------------------
def get_annual_aridity(P_monthly_arr, PET_monthly_arr):

    #----------------------------------------------------------
    # Note: Both precip and PET must have same units, which
    #       will typically be mm or cm.  In GAGES-II, precip
    #       has units of cm and the annual PET it provides
    #       has units of mm.  Values are in [0, 5.8].
    #---------------------------------------------------------- 
    P_sum   = P_monthly_arr.sum()    # [cm]
    PET_sum = PET_monthly_arr.sum()  # [cm]
    aridity = (PET_sum / P_sum)
    return aridity

#   get_annual_aridity()
#---------------------------------------------------------------------
def get_knoben_annual_aridity( MI_monthly_arr ):

    #------------------------------------------------------------
    # Note: I_m is the annual average aridity as defined by
    #       Knoben et al. (2018), which is computed from a
    #       version of Thornthwaite's monthly "moisture index".
    #       See Knoben et al. (2018)
    #------------------------------------------------------------
    MI_sum = np.sum( MI_monthly_arr )
    I_m = (MI_sum / 12)
    
    return I_m

#   get_knoben_annual_aridity()
#---------------------------------------------------------------------
def get_knoben_annual_seasonality( MI_monthly_arr ):

    #------------------------------------------------------------
    # Note: I_mr is a measure of the "seasonality of aridity",
    #       as defined by Knoben et al. (2018). It is simply
    #       the range of the monthly aridity values.
    #       See Knoben et al. (2018)
    #------------------------------------------------------------
    I_mr = MI_monthly_arr.max() - MI_monthly_arr.min()
    return I_mr

#   get_knoben_annual_seasonality()
#---------------------------------------------------------------------
def get_knoben_annual_snow_fraction( P_monthly_arr, T_monthly_arr ):

    T0 = 0  # (rain-to-snow temp threshold from Knoben 2018)
    w = np.where(T_monthly_arr < T0)
    fs = (P_monthly_arr[w].sum() / P_monthly_arr.sum())
    return fs

#   get_knoben_annual_snow_fraction()
#---------------------------------------------------------------------
def calc_knoben_indicators( lat_deg, gages2_val_dict ):

    #---------------------------------------------
    # Note:  lat = latitude for the USGS site ID
    #---------------------------------------------
    P_monthly_arr = g2.get_gages2_monthly_precip_array( gages2_val_dict)
    T_monthly_arr = g2.get_gages2_monthly_temp_array( gages2_val_dict)

    PET_monthly_arr = get_harmon_monthly_PET_array( T_monthly_arr, lat_deg )
    #-----------------------------------------------------
    # Compute standard aridity as "ARIDITY2", as a check
    #-----------------------------------------------------
    aridity = get_annual_aridity( P_monthly_arr, PET_monthly_arr)
       
    MI_monthly_arr = get_knoben_monthly_MI_array( P_monthly_arr, PET_monthly_arr )
    ### get_monthly_MI_array( gages2_val_dict )

    K_aridity = get_knoben_annual_aridity( MI_monthly_arr )
    K_seasonality = get_knoben_annual_seasonality( MI_monthly_arr )
    K_snow_frac = get_knoben_annual_snow_fraction( P_monthly_arr, T_monthly_arr )
    
    return aridity, K_aridity, K_seasonality, K_snow_frac
   
#   calc_knoben_indicators()
#---------------------------------------------------------------------
#---------------------------------------------------------------------
def plot_monthly_precip_values( gages2_val_dict, site_id ):

    # Plot and fit a sine curve.
    # P_monthly_arr = g2.get_gages2_monthly_precip_array( gages2_val_dict)
    pass

#   plot_monthly_precip_values()
#---------------------------------------------------------------------
def plot_monthly_temp_values( gages2_val_dict, site_id ):

    # Plot and fit a sine curve.
    # T_monthly_arr = g2.get_gages2_monthly_temp_array( gages2_val_dict)
    pass

#   plot_monthly_temp_values()
#---------------------------------------------------------------------

## site_id not in dictionary = 01372058
## site_id not in dictionary = 01484085
## site_id not in dictionary = 02300021
## site_id not in dictionary = 02300042
## site_id not in dictionary = 02310747
## site_id not in dictionary = 02322698
## site_id not in dictionary = 02326550
## site_id not in dictionary = 02343801
## site_id not in dictionary = 02403500
## site_id not in dictionary = 07250085
## site_id not in dictionary = 08041749
## site_id not in dictionary = 08041780
## site_id not in dictionary = 08074000
## site_id not in dictionary = 08075500
## site_id not in dictionary = 10172860
## site_id not in dictionary = 10309101
## site_id not in dictionary = 14211720











