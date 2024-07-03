
# Copyright (c) 2023-2024, Scott D. Peckham
#
# Nov 2023. Added Mercator_XY().
#           Added Lambert_Azimuthal_Equal_Area_XY()
# Oct 2023. Moved Albers projection code from hydrofab_tools.py.
#           Add note how Albers_XY differs from Albers_XY2.
#           Moved Lambert_AEA code from hlr_tools.py. 
#
#---------------------------------------------------------------------
#
#  % conda activate tf36
#  % python
#  >>> from topoflow.utils import projections as proj
#  >>> x,y = proj.Lambert_Azimuthal_Equal_Area_XY( lon, lat )
#
#---------------------------------------------------------------------
#
#  Albers_XY()
#  Albers_q()
#  Albers_XY2()
#  Lambert_Azimuthal_Equal_Area_XY()
#  Mercator_XY()
#
#---------------------------------------------------------------------

import numpy as np

#---------------------------------------------------------------------
def Albers_XY( lon_deg, lat_deg, REPORT=True ):

    #----------------------------------------------------
    # Note:  Albers_XY(-96, 23) = (0,0).
    # Bounds of conterminous US are:
    #     maxlat = 49.382808
    #     minlat = 24.521208
    #     maxlon = -66.945392
    #     minlon = -124.736342
    # Albers_XY( minlon, maxlat) =
    # Albers_XY( -124.736342, 49.382808 )
    #----------------------------------------------------
    # Note: lon,lat of Bellingham, Washington =
    #       -122.485886, 48.769768 
    #       lon,lat of Neah Bay, Washington =
    #       (-124d, 36m 33.59s), (48d, 21m, 33.59s)
    #       -124.6250, 48.3681  
    #----------------------------------------------------
    d2r  = (np.pi / 180.)  # (degrees -> radians)
    lat  = lat_deg * d2r
    lon  = lon_deg * d2r
    
    #---------------------------------------------
    # Note: Formulas from Snyder's Book, p. 100.
    #---------------------------------------------
    R    = 6378137      # (radius of Earth, meters)  (using semi-major axis)
    lat1 = 29.5 * d2r   # (standard parallel 1)
    lat2 = 45.5 * d2r   # (standard parallel 2)
    lon0 = -96.0 * d2r  # (center lon)
    lat0 = 23.0 * d2r   # (center lat)
    
    n     = (np.sin(lat1) + np.sin(lat2)) / 2
    C     = (np.cos(lat1)**2) + (2 * n * np.sin(lat1) )
    rho   = (R / n) * np.sqrt(C - (2 * n * np.sin(lat)))
    rho0  = (R / n) * np.sqrt(C - (2 * n * np.sin(lat0))) 
    theta = n * (lon - lon0)

    x = rho * np.sin( theta )
    y = rho0 - (rho * np.cos( theta ))
    
    if (REPORT):
        print('lon, lat =', lon_deg, ',', lat_deg)
        print('x,   y   =', x, ',', y)

    return (x,y)
        
#   Albers_XY()
#------------------------------------------------------------------------
def Albers_q( lat, e ):

    term1 = np.sin(lat) / (1 - (e * np.sin(lat))**2)
    term2 = np.log( (1 - e*np.sin(lat)) / (1 + e*np.sin(lat)) )
    term2 = term2 / (2*e)
    return (1 - e**2) * (term1 - term2)

#   Albers_q
#------------------------------------------------------------------------
def Albers_XY2( lon_deg, lat_deg, REPORT=True ):

    #----------------------------------------------------
    # Note:  Albers_XY2(-96, 23) = (0,0).
    # Bounds of conterminous US are:
    #     maxlat = 49.382808
    #     minlat = 24.521208
    #     maxlon = -66.945392
    #     minlon = -124.736342
    # Geographic center of the USA:
    #     http://www.kansastravel.org/geographicalcenter.htm
    #     lon, lat = -98d 35m, 39d 50m,
    #     lon, lat = -98.58333, 39.833333
    # HUC01 cat-29 (minlon, minlat):
    #     lon, lat = -69.2950005, 46.77957759999999
    #------------------------------------------
    # Albers_XY2( minlon, maxlat) =
    # Albers_XY2( -124.736342, 49.382808 )
    # Also try:
    # Albers_XY2( -131.5, 51.5)
    # Albers_XY2( -98.58333, 38.833333)
    #----------------------------------------------------
    # Note: lon,lat of Bellingham, Washington =
    #       -122.485886, 48.769768 
    #       lon,lat of Neah Bay, Washington =
    #       (-124d, 36m 33.59s), (48d, 21m, 33.59s)
    #       -124.6250, 48.3681  
    #----------------------------------------------------
    d2r = (np.pi / 180.)  # (degrees -> radians)
    lat = lat_deg * d2r
    lon = lon_deg * d2r
    
    #---------------------------------------------
    # Note: Formulas from Snyder's Book, p. 100.
    #---------------------------------------------
    a    = 6378137      # (Earth ellipsoid semi-major axis)
    f    = (1 / 298.257222101004)  # (Earth flattening)
    e    = np.sqrt( f * (2 - f) )
    lat1 = 29.5 * d2r   # (standard parallel 1)
    lat2 = 45.5 * d2r   # (standard parallel 2)
    lon0 = -96.0 * d2r  # (center lon)
    lat0 = 23.0 * d2r   # (center lat)
    
    m1    = np.cos(lat1) / np.sqrt(1 - (e * np.sin(lat1))**2)
    m2    = np.cos(lat2) / np.sqrt(1 - (e * np.sin(lat2))**2)
    q     = Albers_q( lat, e )
    q0    = Albers_q( lat0, e)
    q1    = Albers_q( lat1, e)
    q2    = Albers_q( lat2, e)
    n     = (m1**2 - m2**2) / (q2 - q1)    
    C     = (m1**2 + (n * q1))
    rho   = (a / n) * np.sqrt(C - (n * q))
    rho0  = (a / n) * np.sqrt(C - (n * q0)) 
    theta = n * (lon - lon0)

    x = rho * np.sin( theta )
    y = rho0 - (rho * np.cos( theta ))
  
    if (REPORT):
        print('lon, lat =', lon_deg, ',', lat_deg)
        print('x,   y   =', x, ',', y)
        print()

    return x,y
        
#   Albers_XY2()
#-----------------------------------------------------------------------
def Lambert_Azimuthal_Equal_Area_XY( lon, lat,
            lon0=-100.0, lat1=45.0):

    #---------------------------------------------
    # Note: Formulas from Snyder's Book, p. 185.
    #---------------------------------------------
    # R = radius of Clarke 1866 Authalic Sphere
    # Central_Meridian   = lon0 = -100.0
    # Latitude of Origin = lat0 = 45.0
    #-------------------------------------------------
    # Equivalent PROJ4 command string:
    # +proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0
    #    +a=6370997 +b=6370997 +units=m +no_defs
    #-------------------------------------------------
    # Was able to test using this PROJ4 string at:
    #    https://mygeodata.cloud/cs2cs/
    # and got perfect agreement.
    #-------------------------------------------------
    # EPSG = 4326 is Geographic coords (WGS 84)
    # HLR data uses: Lambert_Azimuthal_Equal_Area
    #   with lon0 = -100, lat1 = 45.0, and both east
    #   easting (+x_0)  and northing (+y_0) = 0.
    #-------------------------------------------------
    # EPSG = 2163 is a Lambert Azimuthal Equal Area
    # projection w/ almost identical info to our PRJ
    # https://spatialreference.org/ref/epsg/2163/
    # https://epsg.io/2163
    # Download PRJ file and compare.
    # Also very similar to EPSG = 9311, which is a
    #  "non-deprecated replacement" for EPSG 2163?
    #-------------------------------------------------     
    R   = 6370997.0   # [meters]
    d2r = (np.pi / 180.0)
    lon_rad  = lon  * d2r   # psi in Snyder
    lon0_rad = lon0 * d2r
    lon_diff_rad = (lon_rad - lon0_rad)
    lat_rad  = lat  * d2r   # phi in Snyder
    lat1_rad = lat1 * d2r

    #----------------------------------------------------------
    # Note: Origin of map projection is at (lon0, lat1).
    #       x is positive to the east of the map origin
    #       and increases from west to east as it should.
    #       (-100, 45)   -> (0.0, 0.0)
    #       (-100.5, 45) -> (-39313.012957, 121.2945992)
    #       (-99.5, 45)  -> ( 39313.012957, 121.2945992)
    #----------------------------------------------------------
    # Note: y is positive to the north of the map origin
    #       and increases from south to north as it should.
    #       (-100, 45.5) -> (0.0,  55597.26072638)
    #       (-100, 44.5) -> (0.0, -55597.26072638)
    #----------------------------------------------------------   
    kterm1 = np.sin(lat1_rad) * np.sin(lat_rad)
    kterm2 = np.cos(lat1_rad) * np.cos(lat_rad) * np.cos(lon_diff_rad)
    kp     = np.sqrt(2 / (1 + kterm1 + kterm2))  # k_prime in Snyder
    #----------------------------------------------------------
    x      = R * kp * np.cos(lat_rad) * np.sin(lon_diff_rad)
    #----------------------------------------------------------
    yterm1 = np.cos(lat1_rad) * np.sin(lat_rad)
    yterm2 = np.sin(lat1_rad) * np.cos(lat_rad) * np.cos(lon_diff_rad)
    y      = R * kp * (yterm1 - yterm2)
    return x,y

#   Lambert_Azimuthal_Equal_Area_XY()
#------------------------------------------------------------------------
def Mercator_XY( lon, lat, lon0=0.0, SPHERE=False):

    #---------------------------------------------
    # Note: Formulas from Snyder's Book, p. 41.
    #---------------------------------------------
    # R = radius of Clarke 1866 Authalic Sphere
    # Central_Meridian   = lon0 = -100.0
    #-------------------------------------------------
    # Equivalent PROJ4 command string:
    # +proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 
    #    +datum=WGS84 +units=m +no_defs
    #-------------------------------------------------
    # Was able to test using this PROJ4 string at:
    #    https://mygeodata.cloud/cs2cs/
    # and got perfect agreement.
    #-------------------------------------------------
    # EPSG code = 3395, for WGS 84 / World Mercator
    # (World between 80 deg South and 84 deg North.)
    #-------------------------------------------------
    ## R = 6371000.7900    # [meters]  sphere used for WGS 84 ?
    ## R = 6366707.0195  # [meters]  sphere used for WGS 84 ?
    R   = 6370997.0   # [meters]  Clarke 1866 Authalic Sphere 
    d2r = (np.pi / 180.0)
    lon_rad  = lon  * d2r   # lambda in Snyder
    lon0_rad = lon0 * d2r
    lon_diff_rad = (lon_rad - lon0_rad)
    lat_rad  = lat  * d2r   # phi in Snyder
    ### lat1_rad = lat1 * d2r

    #----------------------------------------------------------
    # Note: Origin of map projection is at (lon0, 0.0).
    #       x is positive to the east of the map origin
    #       and increases from west to east as it should.
    #       (-100, 45)   -> (-11131949.079327356, 5591295.918553136)
    #----------------------------------------------------------
    if (SPHERE): 
        x = R * lon_diff_rad
        y = R * np.log( np.tan( np.pi/4 + lat_rad/2))
    else:
        a = 6378137.0       # [meters] semi-major axis
        b = 6356752.314245  # [meters] semi-minor axis
        e = np.sqrt(1 - (b/a)**2.0)
        ## inv_f = 298.257223563
        factor = (1 - e*np.sin(lat_rad))/(1 + e*np.sin(lat_rad))
        factor = factor**(e/2)
        x = a * lon_diff_rad
        y = a * np.log( np.tan( np.pi/4 + lat_rad/2) * factor)
    return x,y

#   Mercator_XY()
#------------------------------------------------------------------------


