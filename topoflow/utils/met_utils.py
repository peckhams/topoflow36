
# Copyright (c) 2023, Scott D. Peckham
#
# Aug 2023.  Created to compute RH from spec. humidity.
#
#-------------------------------------------------------------------
#
#  Functions:
#
#  saturation_vapor_pressure()
#  relative_humidity()
#
#-------------------------------------------------------------------

import numpy as np
# from topoflow.utils import met_base

#------------------------------------------------------------------- 
def saturation_vapor_pressure( T, method=None, MBAR=True ):

    #-------------------------------------------------
    # Notes: See met_base.py:
    #           update_saturation_vapor_pressure() &
    #           compare_em_air_methods()
    #
    #        T = T_air = air temperature, Celsius
    #-------------------------------------------------
    T = np.float32(T)  # This sets the data type.
    if (method == None):
        method = 'BRUTSAERT'

    if (method == 'BRUTSAERT'):    
        #------------------------------
        # Use Brutsaert (1975) method
        #------------------------------
        term1 = (17.3 * T) / (T + 237.3)
        e_sat = 0.611 * np.exp( term1 )   # [kPa]
    elif (method == 'SATTERLUND'):    
        #-------------------------------
        # Use Satterlund (1979) method
        #-------------------------------
        term1 = 2353 / (T + 273.15)
        e_sat = 10 ** (11.4 - term1)  # [Pa]
        e_sat = (e_sat / 1000)        # [kPa]
    elif (method == 'BOLTON'):
        #-------------------------------------------------
        # Use Bolton (1980) method, similar to Brutsaert
        #-------------------------------------------------
        #   The computation of Equiv. Potential Temp. 
        #   http://www.eol.ucar.edu/projects/ceop/dm/
        #          documents/refdata_report/eqns.html
        #-------------------------------------------------
        term1 = (17.67 * T) / (T + 243.5)
        e_sat =  0.6112 * np.exp( term1 )  # [kPa]

    #-----------------------------------
    # Convert units from kPa to mbars?
    #-----------------------------------
    if (MBAR):    
        e_sat = (e_sat * 10)   # [mbar]
        
    return e_sat

#   saturation_vapor_pressure()    
#-------------------------------------------------------------------   
def relative_humidity(q_air, T_air, P_surf, method=None):

    #-------------------------------------------------------------- 
    # Notes: Convert specific humidity to relative humidity.
    #        GLDAS has spec. humidity but not RH.
    #--------------------------------------------------------------   
    # See: https://earthscience.stackexchange.com/questions/2360/
    #   how-do-i-convert-specific-humidity-to-relative-humidity
    # Based on R version by David LeBauer.
    #
    # q_air  = specific humidity, dimensionless (kg/kg)
    #        = ratio of water mass / total air mass in air
    # T_air  = air temperature, degrees C
    # P_surf = surface air pressure, mbar
    # RH     = rel. humidity
    #        = actual water mix. ratio / saturation mix. ratio
    #-------------------------------------------------------------- 
    e_sat_air = saturation_vapor_pressure( T_air, method=method)
    e_air     =  q_air * P_surf / ((0.378 * q_air) + 0.622)
    RH        = (e_air / e_sat_air)
 
    #---------------------------------   
    # Make sure RH is in range [0,1]
    #---------------------------------
    if (isinstance(RH, np.ndarray)):
        np.minimum(RH, 1, RH)  # Assume RH is an array
        np.maximum(RH, 0, RH)
    else:
        RH = min(RH, 1)
        RH = max(RH, 0)

    return RH

#   relative_humidity() 
#-------------------------------------------------------------------




    
