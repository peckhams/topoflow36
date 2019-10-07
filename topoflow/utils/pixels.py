
#########################################################
# 1/16/09.  Using the factor "dtor" here looks strange
#           in meters_per_degree_lon() and
#           in meters_per_degree_lat()
#           but seems to be correct.
#########################################################
        
## Copyright (c) 2009-2010, Scott D. Peckham
## January 12-16, 2009
## May 2010 (changes to unit_test())

#-------------------------------------------------------------------
#
#  Functions:
#
#  unit_test()
#
#  get_sizes_by_row()
#  get_da()
#  meters_per_degree_lon()
#  meters_per_degree_lat()
#
#-------------------------------------------------------------------

import numpy as np

from .tf_utils import TF_Print

#-------------------------------------------------------------------
def unit_test():

    from . import rti_files
    
    print('Testing meters_per_degree_lon()...')
    MPD_lon = meters_per_degree_lon(0)
    print('meters_per_degree_lon(0)   =', MPD_lon)
    print('meters_per_3_arcsec_lon(0) =', MPD_lon / 1200.0)
    print(' ')
    #------------------------------------------------------------
    print('Testing meters_per_degree_lat()...')
    MPD_lat = meters_per_degree_lat(0)
    print('meters_per_degree_lat(0)   =', MPD_lat)
    print('meters_per_3_arcsec_lat(0) =', MPD_lat / 1200.0)
    print(' ')
    #------------------------------------------------------------    
    print('Testing get_da() with KY_Sub (fixed-angle pixels)...')
    in_directory  = '/Applications/RIVIX/RiverTools_3.0/basins/KY_Sub/'
    site_prefix = 'KY_Sub'
    RTI_file    = (in_directory + site_prefix + '.rti')
    info = rti_files.read_info( RTI_file, REPORT=True )
    da   = get_da(info, REPORT=True, VERBOSE=True)
    print('da = ')
    print(da)
    print('shape(da) =', shape(da))
    print(' ')
    #------------------------------------------------------------    
    print('Testing get_da() with Beaver... (fixed-length pixels)')
    in_directory  = '/Applications/RIVIX/RiverTools_3.0/basins/Beaver_Creek_KY/'
    site_prefix = 'Beaver'
    RTI_file    = (in_directory + site_prefix + '.rti')
    info = rti_files.read_info( RTI_file, REPORT=True )
    da   = get_da(info, REPORT=True, VERBOSE=True)
    print('da = ')
    print(da)
    print('shape(da) =', shape(da))
    print(' ')
    
#   unit_test()                    
#-------------------------------------------------------------------
def get_sizes_by_row(rti, REPORT=False, METERS=False):

    #----------------------------------------------------------
    # NOTES: This routine returns the xsize, ysize, diagonal
    #        size and area of pixels, based on the pixel
    #        geometry associated with the current DEM, in
    #        kilometers.  RTI_file contains georef rti.

    #        For fixed-angle pixels, the lat/lon-dependence
    #        of both xsizes and ysizes (on an ellipsoid)
    #        is taken into account.

    #        The calculations are done efficiently using
    #        IDL array operations, so dx, dy, dd, and da
    #        are returned as 1D arrays of length NROWS.

    #        The vars dx,dy,dd,da are returned as DOUBLES.
    #        This is necessary.

    #        Note that:
    #            (minlat, maxlat-yresdeg) -> y=(ny-1, 0) and
    #            (minlon, maxlon-xresdeg) -> x=(nx-1, 0).
    #        This is related to the fact that !order=1 in RT.
    #        Lats are for bottom edge of pixel.
    #        Note that 3600 arcsecs = 1 degree.

    #        NOTE that dd=np.sqrt( dx^2 + dy^2) and da=(dx * dy)
    #        are very good approximations as long as dx and
    #        dy are not too large.
    #----------------------------------------------------------
    
    #-------------------------
    # Kilometers or Meters ?
    #-------------------------
    if (METERS):  
        ufactor = np.float64(1)
    else:    
        ufactor = np.float64(1000)
    
    if (rti.pixel_geom == 0):    
        #---------------
        # Compute lats
        #---------------
        DTORD   = (np.pi / np.float64(180))
        ycoords = np.arange(rti.nrows, dtype='int16')
        yresdeg = (rti.yres / np.float64(3600))    #(arcsecs -> degrees)
        lats = (rti.y_north_edge - (yresdeg * (ycoords + np.float64(1))))
        
        #----------------------------------------
        # Compute pixel sizes using lats & lons
        #----------------------------------------
        # NB!  dy becomes a vector !!
        #      Also for (pixel_geom eq 1) below !
        #------------------------------------------
        MPD_LON = meters_per_degree_lon(lats)         #(vector)
        MPD_LAT = meters_per_degree_lat(lats)         #(vector)
        dx = (rti.xres / np.float64(3600) * MPD_LON / ufactor)  #(vector)
        dy = (rti.yres / np.float64(3600) * MPD_LAT / ufactor)  #(vector)
        dd = np.sqrt(dx ** 2 + dy ** 2)
        da = (dx * dy)   
    else:    
        dx = (rti.xres / ufactor)   #(meters or km)
        dy = (rti.yres / ufactor)
        dd = np.sqrt(dx ** 2 + dy ** 2)
        da = (dx * dy)
        #-------------------------------
        # Return dx, dy, etc. as grids
        #-------------------------------
        dx = np.zeros([rti.nrows], dtype='float64') + dx
        dy = np.zeros([rti.nrows], dtype='float64') + dy
        dd = np.zeros([rti.nrows], dtype='float64') + dd
        da = np.zeros([rti.nrows], dtype='float64') + da
    
    #------------------
    # Optional report
    #------------------
    if (REPORT):    
        if (rti.pixel_geom == 0):    
            RADEG = (np.float64(180) / np.pi)
            #---------------------
            print('Pixel geometry = Fixed-angle ')
            print(' ')
            print('Actual south edge lat   = ' + str(rti.y_south_edge))
            print('Computed south edge lat = ' + str(lats[rti.nrows - 1] * RADEG))
            print('Actual north edge lat   = ' + str(rti.y_north_edge))
            print('Computed north edge lat = ' + str((lats[0] + yresdeg) * RADEG))
        else:    
            print('Pixel geometry = Fixed-length ')
            print(' ')
            print('Actual south edge y   = ' + str(rti.y_south_edge))
            print('Actual north edge y   = ' + str(rti.y_north_edge))
        print(' ')
        print('Min(dx), Max(dx) = ' + str(dx.min()) + str(dx.max()))
        print('Min(dy), Max(dy) = ' + str(dy.min()) + str(dy.max()))
        print('Min(dd), Max(dd) = ' + str(dd.min()) + str(dd.max()))
        print('Min(da), Max(da) = ' + str(da.min()) + str(da.max()))
        print(' ')

    return (dx, dy, dd)

#   get_sizes_by_row()
#-------------------------------------------------------------------
def get_da(rti, METERS=False, REPORT=False, VERBOSE=False):
    
    #------------------------------------
    # Pixel dimensions; convert km to m
    # These are planform dimensions.
    #---------------------------------------
    dx, dy, dd = get_sizes_by_row(rti, METERS=True)
    
    #-------------------------------------------
    # 7/13/06.  Allow da to be scalar or grid.
    # For speed;  was always grid before.
    #-------------------------------------------
    if (rti.pixel_geom == 1):
        da = (dx[0] * dy[0])
        
        if (VERBOSE):
            TF_Print('dx = ' + str(dx[0]) + '  [m]')
            TF_Print('dy = ' + str(dy[0]) + '  [m]')
            TF_Print('da = ' + str(da) + '  [m^2]')
    else:
        #---------------------------------
        # Convert da from 1D to 2D array
        # Then subscript with the wk's.
        #---------------------------------
        TF_Print('Computing pixel area grid...')
        nx = rti.ncols
        ny = rti.nrows

        ######################################
        # DOUBLE CHECK THIS.  SEEMS CORRECT
        # IS THERE A MORE EFFICIENT WAY ??
        ######################################
        da_by_row = (dx * dy)
        da = np.reshape(np.repeat(da_by_row, nx), (ny, nx))

        #------------------------------------------------------
        ## This resulted from I2PY and doesn't seem right.
        ## matrixmultiply() was depracated in favor of dot().
        #------------------------------------------------------
        ## self.da = (transpose(matrixmultiply(transpose(repeat(1,nx)), \
        ##                                     transpose(self.da_by_row))))

        if (VERBOSE):
            da_min = da.min()
            da_max = da.max()
            TF_Print('    min(da) = ' + str(da_min) + '  [m^2]')
            TF_Print('    max(da) = ' + str(da_max) + '  [m^2]')
            dx_min = dx.min()
            dx_max = dx.max()
            TF_Print('    min(dx) = ' + str(dx_min) + '  [m]')
            TF_Print('    max(dx) = ' + str(dx_max) + '  [m]')
            dy_min = dy.min()
            dy_max = dy.max()
            TF_Print('    min(dy) = ' + str(dy_min) + '  [m]')
            TF_Print('    max(dy) = ' + str(dy_max) + '  [m]')
            TF_Print(' ')

    return da

#   get_da()
#-------------------------------------------------------------------
def meters_per_degree_lon(lat_deg,
                          a_radius=np.float64(6378137.0),  # [meters]
                          b_radius=np.float64(6356752.3),  # [meters]
                          mean_elev=np.float64(0)):

    #----------------------------------------------------------
    # NOTES: This formula comes from both:
    #        Snyder, J.P. (1987) Map projections - A working
    #        manual, USGS Prof. Paper 1395, p. 25.
    #        and

    #        Ewing, C.E., and Mitchell, M.M. (1970)
    #        Introduction to Geodesy, American Elsevier
    #        Publishing Company, New York, New York, pp. 8-26

    #        The above references give a formula for the
    #        radius of curvature of a parallel of latitude,
    #        Rp, on page 19.  Arclength between two longitudes
    #        lam1 and lam2 is then:  Rp(lat) * (lam2 - lam1).
    #        Note that the argument, lat_deg, is the geodetic
    #        latitude in decimal degrees, and not the
    #        geocentric latitude; otherwise the square root
    #        term wouldn't be needed.  Every parallel of
    #        latitude is a circle.

    #        Based on the diagram on page 14, it can be seen
    #        that if N is extended by the mean_elevation, then
    #        we must add (mean_elev * np.cos(lat)) to Rp.

    #        Checked with an example from an ESRI manual.
    #        For Clarke 1866 ellipsoid, one degree of latitude
    #        at Equator is 111.321 km, and at 60 degrees north
    #        it is 55.802 km.

    #        a_radius = equatorial radius of ellipsoid (meters)
    #        b_radius = polar radius of ellipsoid (meters)

    #        For WGS_1984 ellipsoid:
    #            a_radius = 6378.1370, b_radius = 6356.7523
    #        For CLARKE_1866 ellipsoid:
    #            a_radius = 6378.2064, b_radius = 6356.5838
    #        Default arguments are for WGS_1984.
    #-----------------------------------------------------------
    
    #-------------------------------
    # Compute flattening ratio, f,
    # and eccentricity, e
    #-------------------------------
    f = (a_radius - b_radius) / a_radius
    e = np.sqrt((np.float64(2) * f) - f ** np.float64(2))
    
    #--------------------
    # Compute the value
    #--------------------
    dtor = np.pi / np.float64(180)
    lat_rad = (lat_deg * dtor)
    p1 = mean_elev * np.cos(lat_rad)
    p2 = a_radius * np.cos(lat_rad) / np.sqrt(np.float64(1) - \
                                        (e * np.sin(lat_rad)) ** np.float64(2))

    return (dtor * (p1 + p2))
    
#   meters_per_degree_lon()
#-------------------------------------------------------------------
def meters_per_degree_lat(lat_deg,
                          a_radius=np.float64(6378137.0),  # [meters]
                          b_radius=np.float64(6356752.3),  # [meters]
                          mean_elev=np.float64(0)):

    #----------------------------------------------------------
    # NOTES: This formula comes from both:
    #        Snyder, J.P. (1987) Map projections - A working
    #        manual, USGS Prof. Paper 1395, p. 25.
    #        and

    #        Ewing, C.E., and Mitchell, M.M. (1970)
    #        Introduction to Geodesy, American Elsevier
    #        Publishing Company, New York, New York, pp. 8-26.

    #        The above references give a formula for the
    #        radius of curvature of a meridian of longitude,
    #        R, on page 19.  Arclength along a meridian is
    #        given by the integral of R(phi) * d_phi between
    #        two geodetic latitudes, phi1 and phi2.
    #        Note that the argument, lat_deg, is the geodetic
    #        latitude in decimal degrees, and not the
    #        geocentric latitude. Every meridian of longitude
    #        is the same ellipse.

    #        Based on the diagram on page 18, it can be seen
    #        that if the mean elevation of a region is nonzero,
    #        then R should be replaced by (R + mean_elev).

    #        Default arguments are for WGS_1984.
    #----------------------------------------------------------
    
    #-------------------------------
    # Compute flattening ratio, f,
    # and eccentricity, e
    #-------------------------------
    f = (a_radius - b_radius) / a_radius
    e = np.sqrt((np.float64(2) * f) - f ** np.float64(2))
    
    #--------------------
    # Compute the value
    #--------------------
    dtor    = np.pi / np.float64(180)
    lat_rad = (lat_deg * dtor)
    p = a_radius * (np.float64(1) - e ** np.float64(2)) / (np.float64(1) - (e * np.sin(lat_rad)) ** np.float64(2)) ** np.float64(1.5)
       
    return (dtor * (mean_elev + p))
    
#   meters_per_degree_lat()
#-------------------------------------------------------------------    

