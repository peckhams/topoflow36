"""
This module contains some scripts for extracting information from
a NextGen "hydrofabric" consisting of 2 GeoJSON files:
   catchment-data.geojson and nexus-data.geojson
such as geographic bounding box and outlet coordinates.
The get_dem() method subsets a GDAL VRT DEM to create a DEM for
a specified catchment.
"""
#------------------------------------------------------------------------
#
#  Copyright (C) 2022.  Scott D. Peckham
#
#  Example use:
#
#  >>> import hydrofab_tools as hft
#  >>> c = hft.catchment()
#  >>> c.print_info('cat-29')
#  >>> dem = c.get_dem('cat-29')
#
#  Note: To use this code, you should first create a conda environment
#        called "hydrofabric" (or similar) and install geopandas
#        and gdal in that environment via:
#        % conda create --name hydrofabric
#        % conda activate hydrofabric
#        % conda install gdal (OR conda install -c conda-forge gdal)
#
#  Note: In order to use the get_dem() method as written you need to
#        first install the AWS CLI for your OS, and then run the
#        command "aws configure".  After that, GDAL will be able to
#        find the required AWD configuration options in the files:
#        ~/.aws/credentials and  ~/.aws/config.
#       
#------------------------------------------------------------------------

from osgeo import gdal, osr
# import glob

# Note: This version uses json (built-in) instead of geopandas
import json
import numpy as np

#------------------------------------------------------------------------
"""
This module contains some scripts for extracting information from
a NextGen "hydrofabric" consisting of 2 GeoJSON files:
   catchment-data.geojson and nexus-data.geojson
such as geographic bounding box and outlet coordinates.
The get_dem() method subsets a GDAL VRT DEM to create a DEM for
a specified catchment.
"""
#------------------------------------------------------------------------
#
#  Copyright (C) 2022.  Scott D. Peckham
#
#  Example use:
#
#  >>> import hydrofab_tools as hft
#  >>> c = hft.catchment()
#  >>> c.print_info('cat-29')
#  >>> dem = c.get_dem('cat-29')
#
#  Note: To use this code, you should first create a conda environment
#        called "hydrofabric" (or similar) and install gdal in that
#        environment via:
#        % conda create --name hydrofabric
#        % conda activate hydrofabric
#        % conda install gdal (OR conda install -c conda-forge gdal)
#
#  Note: In order to use the get_dem() method as written you need to
#        first install the AWS CLI for your OS, and then run the
#        command "aws configure".  After that, GDAL will be able to
#        find the required AWD configuration options in the files:
#        ~/.aws/credentials and  ~/.aws/config.
#       
#------------------------------------------------------------------------

from osgeo import gdal, osr
# import glob

# Note: This version uses json (built-in) instead of geopandas
import json
import numpy as np

#------------------------------------------------------------------------
#
#  class catchment()
#      __init__()
#      read_data()
#      get_bounding_box()
#      get_raster_cellsize()   #####
#      get_raster_bounds()     #####
#      get_dem()
#      get_outlet_coords()
#      get_headwater_catchments()  # 11/30/22
#      get_area_sqkm()
#      ocean_outlet_count()
#      print_info()
#
#  get_hydrofab_dem_info()
#  Albers_XY()  # (function)
#  Albers_q()
#  Albers_XY2()
#
#------------------------------------------------------------------------
class catchment:
    #--------------------------------------------------------------------
    def __init__(self, catchment_data_file=None, nexus_data_file=None,
                       spatial_dir=None):

        # Rename data_file to json_file ??
        if (spatial_dir is None):
            # Relative to the top level of ngen repo folder.
            # spatial_dir = './data/topoflow/spatial/'
            spatial_dir = '/Users/peckhams/Dropbox/GitHub/ngen/data/topoflow/spatial/'
            ## spatial_dir = './'
        if (catchment_data_file is None):
            catchment_data_file = 'catchment_data_HUC01.geojson'
        if (nexus_data_file is None):
            nexus_data_file     = 'nexus_data_HUC01.geojson'
        self.catchment_data_file = spatial_dir + catchment_data_file
        self.nexus_data_file     = spatial_dir + nexus_data_file
        #------------------
        self.read_data()
    
    #   __init__()    
    #--------------------------------------------------------------------
    def read_data(self):

        #--------------------------------------------------------------     
        # Note: This version reads the GeoJSON files directly without
        #       adding geopandas as a dependency.  It reads the data
        #       into a "simulated dataframe".
        #--------------------------------------------------------------
        print('Reading catchment info files...')
        cat_unit = open( self.catchment_data_file, 'r')
        nex_unit = open( self.nexus_data_file, 'r')
        self.catchment_data = json.load( cat_unit )
        self.nexus_data     = json.load( nex_unit )
        cat_unit.close()
        nex_unit.close()
        nc = len(self.catchment_data['features'])
        print('Read info for', nc, 'catchments.')
        print()
                       
    #   read_data()   
    #--------------------------------------------------------------------
    def get_bounding_box(self, cat_id_str, buffer_sec=0, REPORT=False):

        #--------------------------------------------------------------
        # Notes:  Get bounding box coords for a given catchment from
        #         its id string, e.g. id_str = 'cat-29'.
        #         Polygon coordinates are stored in:
        #            catchment-data.geojson (POLYGON type geometry)
        #         Coordinates can be given as Geographic or Albers.
        #--------------------------------------------------------------   
        df = self.catchment_data['features']  # returns list of dictionaries
        crs_name = self.catchment_data['crs']['properties']['name']
        record = next((item for item in df if item["id"] == cat_id_str), False)
        if not(record):
            print('SORRY, cat_id not found:', cat_id_str)
            print('Returning.')
            print()
            return -1

        rgc = record['geometry']['coordinates'][0]
        pairs   = np.array( rgc, dtype='float64')
        xcoords = pairs[:,0]
        ycoords = pairs[:,1]
        self.minlon = xcoords.min()
        self.minlat = ycoords.min()
        self.maxlon = xcoords.max()
        self.maxlat = ycoords.max()
        bounds = [self.minlon, self.minlat, self.maxlon, self.maxlat] 
                
        # print('## record["id"] =', record['id'])
        # print('## record["geometry"]["coordinates"][0] =', rgc)
        # print('## type(rgc) =', type(rgc) )  # list
        # print('## rgc[0] =', rgc[0])
        # print('## xcoords[0:2] =', xcoords[0:2])
        # print('## ycoords[0:2] =', ycoords[0:2])
        # print()

        #------------------------------------------------
        # Option to specify a buffer in arcseconds
        #------------------------------------------------
        # (3 arcsecs) * (1 deg / 3600 arcsecs) =
        # 0.000833333 degrees  ~ 92.6 meters
        #------------------------------------------------
        # (6 arcsecs) * (1 deg / 3600 arcsecs) =
        # 0.001666666 degrees  ~ 185.2 meters        
        #------------------------------------------------
        if (buffer_sec > 0):
           buffer_deg = (buffer_sec / 3600.0)
           self.minlon -= buffer_deg
           self.minlat -= buffer_deg
           self.maxlon += buffer_deg
           self.maxlat += buffer_deg
           bounds = [self.minlon, self.minlat, self.maxlon, self.maxlat]

        #---------------------------------------------------------
        # Option to specify a buffer in decimal degrees, as a %.
        #---------------------------------------------------------
#         if (buffer_pct > 0):
#             lon_range = (self.maxlon - self.minlon)
#             lat_range = (self.maxlat - self.minlat)
#             dx = buffer_pct * lon_range
#             dy = buffer_pct * lat_range
#             self.minlon -= dx
#             self.minlat -= dy
#             self.maxlon += dx
#             self.maxlat += dy
#             bounds = [self.minlon, self.minlat, self.maxlon, self.maxlat]
     
        if (REPORT):
            print('Bounding box for: ' + cat_id_str)
            print('  minlon, maxlon =', bounds[0], ', ', bounds[2])
            print('  minlat, maxlat =', bounds[1], ', ', bounds[3])
            print('  crs_name       =', crs_name )
            print()

        self.bounds = bounds   
        return bounds
        
    #   get_bounding_box()
    #-------------------------------------------------------------------
    def get_raster_cellsize( self, gdal_unit ):

        # Note:  This is from regrid.py.
        geotransform = gdal_unit.GetGeoTransform()
        # ulx  = geotransform[0]
        xres = geotransform[1]
        # xrtn = geotransform[2]
        #-----------------------
        # uly  = geotransform[3]
        # yrtn = geotransform[4]  # (not yres !!)
        yres = geotransform[5]  # (not yrtn !!)
    
        return (xres, abs(yres))

    #   get_raster_cellsize()
    #------------------------------------------------------------------- 
    def get_raster_bounds( self, ds, VERBOSE=False):

        #-------------------------------------------------------------
        # Note:  The bounds depend on the map projection and are not
        # necessarily a Geographic bounding box of lons and lats.
        # ds = dataset = "gdal_unit"   ########  
        #-------------------------------------------------------------
        # See:
        # https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html
        # and search on "geotransform".  An example of gdal.SetGeoTransform
        # gives: [xmin, pixel_size, 0, ymax, 0, -pixel_size].
        # Also says args are:
        # [ulx, xDist, rtnX, uly, yDist, rtnY]
        # This is consistent with information below.
        #-------------------------------------------------------------    
        # ulx = upper left x  = xmin
        # uly = upper left y  = ymax
        # lrx = lower right x = xmax
        # lry = lower right y = ymin
        #-----------------------------

        #----------------------------------------------------------  
        # Notice the strange order of parameters here is CORRECT.
        # It is not:  ulx, xres, xskew, uly, yres, yskew
        #----------------------------------------------------------
        ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
        lrx = ulx + (ds.RasterXSize * xres)
        lry = uly + (ds.RasterYSize * yres)

        if (VERBOSE):
            print('ulx, uly   =', ulx, uly)
            print('lrx, lry   =', lrx, lry)
            print('xres, yres = ', xres, yres)
            print('xskew, yskew =', xskew, yskew)
            print('----------------------------------')

        return [ulx, lry, lrx, uly]  # [xmin, ymin, xmax, ymax]

    #   get_raster_bounds()
    #--------------------------------------------------------------------
    def get_dem( self, cat_id_str='cat-29', RESAMPLE_ALGO='bilinear',
                 CLIP_ALBERS=False, CLIP_GEO=False,
                 zfactor=100, REPORT=True,
                 SAVE_DEM=True, DEM_out_file=None):

        #------------------------------------------------------------
        # Note:  Use GDAL VRT format to get DEMs for catchments.
        #------------------------------------------------------------
        # Note:  The CLIP_GEO option will return a DEM that uses
        #        Geographic coordinates and has fixed-angle grid
        #        cell sizes (e.g. 2 arcseconds).
        #        The CLIP_ALBERS option will return a DEM that uses
        #        Albers Conic Equal Area coordinates and has fixed-
        #        length grid cell sizes (e.g. 100 meters).
        #        The original hydrofabric DEM uses this same Albers
        #        projection, but the catchment polygons described
        #        in the catchment & nexus GeoJSON files use
        #        Geographic coordinates.
        #------------------------------------------------------------             
        if not(CLIP_ALBERS or CLIP_GEO):
            CLIP_ALBERS=True
        if (SAVE_DEM and (DEM_out_file is None)):
            DEM_out_file = cat_id_str + '_rawDEM.tif'
            # DEM_out_file = cat_id_str + '_DEM.tif'
           
        #--------------------------------------------------------  
        # This shows some of the other keywords to gdal.Warp.
        #--------------------------------------------------------      
        # WarpOptions(options=[], format=None, outputBounds=None,
        # outputBoundsSRS=None, xRes=None, yRes=None, targetAlignedPixels=False,
        # width=0, height=0, srcSRS=None, dstSRS=None, srcAlpha=False,
        # dstAlpha=False, warpOptions=None, errorThreshold=None,
        # warpMemoryLimit=None, creationOptions=None, outputType=GDT_Unknown,
        # workingType=GDT_Unknown, resampleAlg=None, srcNodata=None,
        # dstNodata=None, multithread=False, tps=False, rpc=False,
        # geoloc=False, polynomialOrder=None, transformerOptions=None,
        # cutlineDSName=None, cutlineLayer=None, cutlineWhere=None,
        # cutlineSQL=None, cutlineBlend=None, cropToCutline=False,
        # copyMetadata=True, metadataConflictValue=None,
        # setColorInterpretation=False, callback=None, callback_data=None)
        # 
        # Create a WarpOptions() object that can be passed to gdal.Warp()
        # Keyword arguments are : options --- can be be an array of strings,
        # a string or let empty and filled from other keywords.
                    
        #-----------------------------------   
        # Specify the resampling algorithm
        #-----------------------------------
        algo_dict = {
        'nearest'     : gdal.GRA_NearestNeighbour,
        'bilinear'    : gdal.GRA_Bilinear,
        'cubic'       : gdal.GRA_Cubic,
        'cubicspline' : gdal.GRA_CubicSpline,
        'lanczos'     : gdal.GRA_Lanczos,
        'average'     : gdal.GRA_Average,
        'min'         : gdal.GRA_Min,
        'max'         : gdal.GRA_Max,
        'mode'        : gdal.GRA_Mode,
        'med'         : gdal.GRA_Med }
    
        resample_algo = algo_dict[ RESAMPLE_ALGO ]

        #----------------------------------------------------    
        # Set config options for accessing an AWS S3 bucket?
        #-----------------------------------------------------------------
        # If you have installed the AWS CLI for your OS, and have run
        # the command "aws configure", then these do not need to be set.
        # In that case, GDAL will find them in ~/.aws/credentials and
        # ~/.aws/config.
        #-----------------------------------------------------------------
        # print('Setting GDAL config options...')
        # gdal.SetConfigOption('AWS_REGION', 'us-east-1')
        # gdal.SetConfigOption('AWS_SECRET_ACCESS_KEY', 'MY_SECRET_ACCESS_KEY')
        # gdal.SetConfigOption('AWS_ACCESS_KEY_ID', 'MY_ACCESS_KEY_ID')

        #-----------------------------------------------------------
        # vsis3 = vsi s3'. VSI = Virtual Server Infrastructure
        #-----------------------------------------------------------
        # VRT = Virtual Raster Table (a GDAL format)
        # A VRT file is a type of XML file.
        # From Mike: vrt_file = "s3://nextgen-hydrofabric/dem.vrt"
        # I think this is considered a "public" S3 bucket.
        # The projection is: Albers Conic Equal Area and spans the
        # entire US with 30 meter grid cell size.
        # Mike says that elevation units are "cm", but this is not
        #   specified in the VRT metadata. So must divide by 100.
        #-----------------------------------------------------------
        print('Opening VRT file in S3 bucket...')
        bucket_name = 'nextgen-hydrofabric/'   # (AWS S3 bucket name)
        dem_file = 'dem.vrt'
        vrt_path = '/vsis3/' + bucket_name + dem_file
        vrt_unit = gdal.Open( vrt_path )

        #--------------------------------------------------------
        # Get info about the full DEM without downloading it
        # (Should be able to get anything in the VRT XML file.)
        #--------------------------------------------------------
        (dx, dy) = self.get_raster_cellsize( vrt_unit )  #########
        vrt_xres = dx
        vrt_yres = dy
        # This block assumes xres is Degrees (it isn't here)
        ## vrt_xres_deg = dx
        ## vrt_yres_deg = dy
        ## vrt_xres_sec = (vrt_xres_deg * 3600.0)
        ## vrt_yres_sec = (vrt_yres_deg * 3600.0)
        vrt_ncols    = vrt_unit.RasterXSize
        vrt_nrows    = vrt_unit.RasterYSize
        # geotransform = ulx, xres, xskew, uly, yskew, yres
        vrt_gt       = vrt_unit.GetGeoTransform()   # (tuple)
        vrt_ulx      = vrt_gt[0]
        vrt_uly      = vrt_gt[3]
        vrt_lry      = vrt_uly - (vrt_nrows * vrt_yres)   #########
        vrt_bounds   = self.get_raster_bounds( vrt_unit )
        vrt_proj     = vrt_unit.GetProjection()
        vrt_meta     = vrt_unit.GetMetadata()   # (can be empty dict)
        vrt_srs      = osr.SpatialReference(wkt=vrt_proj)
        vrt_geog_cs  = vrt_srs.GetAttrValue('geogcs')

        if (vrt_srs.IsProjected):
            vrt_proj_cs = vrt_srs.GetAttrValue('projcs')
        else:
            vrt_proj_cs = 'None'
        #-----------------------------------------------------
        # WARNING: This may cause full DEM to be downloaded.
        #-----------------------------------------------------
        # vrt_band   = vrt_unit.GetRasterBand(1)
        # vrt_nodata = vrt_band.GetNoDataValue()
        # zmin,zmax,zmean,zstdv = vrt_band.GetStatistics(True, True)
        
        ## vrt_dem    = vrt_unit.ReadAsArray()
        ## dem_dtype  = vrt_dem.dtype
        ## zmin       = vrt_dem.min()
        ## zmax       = vrt_dem.max()        
        #----------------------------------------------------     
        if (REPORT):
            print('VRT DEM Information:')
            print('  ncols    =', vrt_ncols)
            print('  nrows    =', vrt_nrows)
            print('  xres     =', vrt_xres)
            print('  yres     =', vrt_yres)
            print('  xy_units  =', 'meters')
            print('  DEM units =', 'cm')     # not saved in metadata!
            print('  DEM dtype =', 'int32')  # long, 4-byte int
            print('  bounds   =', vrt_bounds)
            print('  geog cs  =', vrt_geog_cs)
            print('  proj cs  =', vrt_proj_cs)
            print('  metadata =', vrt_meta)
            ##### print('  proj     =', vrt_proj)    # WKT, wordy
            ## print('zmin    =', zmin)
            ## print('zmax    =', zmax)
            ## print('VRT DEM dtype =', dem_dtype)
            ## print('VRT DEM min, max  =', dem_min, ', ', dem_max)
            print()

        #---------------------------------------------
        # Check the VRT bounding box info (all good)
        #---------------------------------------------
#         vrt_xrange = vrt_bounds[2] - vrt_bounds[0]
#         vrt_yrange = vrt_bounds[3] - vrt_bounds[1]
#         print('vrt_xrange =', vrt_xrange)
#         print('  ncols * xres =', vrt_ncols * vrt_xres)
#         print('vrt_yrange =', vrt_yrange)
#         print('  nrows * yres =', vrt_nrows * vrt_yres)
#         print()
        
        #------------------------------------------
        #  GDAL bounds  = [ulx, lry, lrx, uly] 
        #               = [xmin, ymin, xmax, ymax]
        #---------------------------------------------------------
        # Note:  get_bounding_box() uses the catchment & nexus
        #        GeoJSON files and returns a geographic bounding
        #        box: minlon, minlat, maxlon, maxlat
        #---------------------------------------------------------           
        out_bounds_geo = self.get_bounding_box( cat_id_str, REPORT=True )

        #--------------------------------------------------------
        # Set the output Spatial Reference System via EPSG code
        #--------------------------------------------------------
        # The NextGen Hydrofabric SRC appears to be unique and
        # does not have an EPSG code.  However, this one seems
        # to be a somewhat close match (EPSG: 5070)        
        #
        # Also see:
        # https://desktop.arcgis.com/en/arcmap/10.3/guide-books/
        #         map-projections/albers-equal-area-conic.htm
        # https://spatialreference.org/ref/esri/102003/html/
        # https://en.wikipedia.org/wiki/Albers_projection
        #----------------------------------------------------
        # EPSG Codes via https://epsg.io/:
        #
        # 4326: WGS 84 - WGS84 - World Geodetic System 1984
        # 5070: NAD83 / Conus Albers
        # 7019: GRS 1980
        # 6269: MAGNA-SIRGAS / San Jose del Guaviare urban grid
        # 9122: degree (supplier to define representation)
        # 4269: NAD83
        # 9001: metre
        #----------------------------------------------------
        # Convert coords from lat/lon to Hydrofabric Albers
        #----------------------------------------------------
        wgs84 = osr.SpatialReference()
        wgs84.ImportFromEPSG( 4326 )
        wkt_albers = vrt_unit.GetProjection()
        # wkt = wkt.replace('unnamed', 'Albers')  ####
        albers = osr.SpatialReference()
        albers.ImportFromWkt( wkt_albers )
        ### albers.ImportFromEPSG( 5070 )   #####
       
        #-------------------------------------------------
        # Option to check the WKT (well known text) info
        #-------------------------------------------------
        # See: https://www.geoapi.org/3.0/javadoc/org/
        #      opengis/referencing/doc-files/WKT.html
        #-------------------------------------------------
#         wkt1 = wgs84.ExportToWkt()        
#         wkt1 = wgs84.ExportToPrettyWkt()
#         print('##### wkt for WGS 84 =')
#         print( wkt1 )
#         print()
#         wkt2 = albers.ExportToPrettyWkt()
#         print('##### wkt2 for Albers =')
#         print( wkt2 )
#         print()

        #--------------------------------------------------
        # Option to clip using original Albers projection
        # or the similar EPSG:5070 projection
        #--------------------------------------------------
        if (CLIP_ALBERS):  
            #-----------------------------------------------------
            # Transform geographic bounding box to Albers extent
            #---------------------------------------------------------------
            albers2latlon = osr.CoordinateTransformation( albers, wgs84 )
            latlon2albers = osr.CoordinateTransformation( wgs84, albers )
            #---------------------------------------------------------------
            minlon = out_bounds_geo[0]
            minlat = out_bounds_geo[1]
            maxlon = out_bounds_geo[2]
            maxlat = out_bounds_geo[3]
            #-----------------------------------------------------------
            # NOTE:  When using the GDAL Python package version 3.*
            #        the order of arguments to TransformPoint must be
            #        flipped (i.e. y,x vs. x,y) as compared to version
            #        2.*.  This means that most online examples (that
            #        use version 2.*) are wrong.
            #        See:  https://github.com/OSGeo/gdal/issues/1546
            #-----------------------------------------------------------
            albers_ul_xyz = latlon2albers.TransformPoint( maxlat, minlon )
            albers_ll_xyz = latlon2albers.TransformPoint( minlat, minlon )
            albers_ur_xyz = latlon2albers.TransformPoint( maxlat, maxlon )
            albers_lr_xyz = latlon2albers.TransformPoint( minlat, maxlon )
            print('Albers UL xyz:', albers_ul_xyz)
            print('Albers LL xyz:', albers_ll_xyz)
            print('Albers UR xyz:', albers_ur_xyz)
            print('Albers LR xyz:', albers_lr_xyz)
            print()
            xvals = [albers_ul_xyz[0], albers_ll_xyz[0],
                     albers_ur_xyz[0], albers_lr_xyz[0]]
            yvals = [albers_ul_xyz[1], albers_ll_xyz[1],
                     albers_ur_xyz[1], albers_lr_xyz[1]]
            xmin = np.floor( np.min( xvals ) )
            xmax = np.ceil( np.max( xvals ) )
            ymin = np.floor( np.min( yvals ) )
            ymax = np.ceil( np.max( yvals ) )
            #----------------------------------------------
            out_srs = None
            # out_srs = "EPSG:5070"
            out_bounds_albers = [xmin, ymin, xmax, ymax]
            out_bounds   = out_bounds_albers
            out_xy_units = 'meters'
            #-----------------------------------------------------
            # For hydrofabric DEM (dem.vrt), xres = yres = 30 m,
            # which is roughly 1 arcsecond at the equator.
            #-----------------------------------------------------
            out_xres = 90.0    # [meters]   ########
            out_yres = 90.0    # [meters]   ########

            #------------------------------------------
            # Option to compute new DEM cols and rows
            #   relative to original VRT DEM
            #------------------------------------------
#             col_min = np.int32((xmin - vrt_ulx) / vrt_xres)
#             col_max = np.int32((xmax - vrt_ulx) / vrt_xres)
#             row_min = np.int32((ymin - vrt_lry) / vrt_yres)
#             row_max = np.int32((ymax - vrt_lry) / vrt_yres)
#     #         row_min = np.int32((ymin - vrt_uly) / vrt_yres)
#     #         row_max = np.int32((ymax - vrt_uly) / vrt_yres)
#             print('col_min, col_max =', col_min, col_max)
#             print('row_min, row_max =', row_min, row_max)
#             print()
     
        #--------------------------------------------------------
        # Option to clip using Geographic bounding box (WGS 84)
        # In gdal.Warp, must set: outputBounds=out_bounds
        # and check other keywords like: xRes, yRes.
        #--------------------------------------------------------
        if (CLIP_GEO):      
            out_srs = "EPSG:4326"
            out_bounds = out_bounds_geo
            out_xy_units = 'degrees'
            # If dstSRS = WGS84, out_xres units must be degrees
            out_xres = 3.0/3600  # 3 arcseconds -> degrees
            out_yres = 3.0/3600

        #----------------------------------------
        # Change/increase the grid cell size ??
        #----------------------------------------
        # out_xres = vrt_xres   # no change
        # out_yres = vrt_yres
        #-------------------------
        # out_xres = (2 * vrt_xres)
        # out_yres = (2 * vrt_yres)
        #----------------------------------------------------------
        # For hydrofabric DEM (dem.vrt), xres = yres = 30 meters,
        # which is roughly 1 arcsecond at the equator.
        # If dstSRS = WGS84, out_xres units must be degrees
        #----------------------------------------------------------        
        # out_xres = 2./3600  # (2 arcseconds -> degrees)
        # out_yres = 2./3600  # (2 arcseconds -> degrees)

        #------------------------------------------------  
        # Resample & clip and write new grid to GeoTIFF
        #------------------------------------------------
        # See notes near "METHOD 2" below.
        #------------------------------------------------
        # Only set xRes and yRes if they are Degrees
        #------------------------------------------------        
        print('Creating DEM for ' + cat_id_str + '...')     
        out_unit = gdal.Warp( DEM_out_file, vrt_unit,
            format = None,  # (don't save to file yet)
            #### format = 'GTiff',  # (output format string)
            #### outputType=gdal.GDT_Float32,
            outputBounds=out_bounds,  #######
            dstSRS=out_srs,  #######
            ### srcSRS="EPSG:####",  # don't need to specify
            xRes=out_xres, yRes=out_yres,  #################
            resampleAlg = resample_algo )

        #----------------------
        # Get info on new DEM
        #----------------------
        ## print('Getting new DEM info...')
        out_ncols  = out_unit.RasterXSize
        out_nrows  = out_unit.RasterYSize
        (out_xres, out_yres) = self.get_raster_cellsize( out_unit )  #####
        ## out_bounds = self.get_raster_bounds( out_unit )
        out_band   = out_unit.GetRasterBand(1)
        VRT_nodata = out_band.GetNoDataValue()
        out_gt     = out_unit.GetGeoTransform()   # (tuple)
        dem        = out_unit.ReadAsArray()

        #-------------------------------------------
        # Convert GDAL DataType to Numpy data type
        #-------------------------------------------
#         gdal_type_map = {1: "int8", 2: "uint16", 3: "int16",
#                          4: "uint32", 5: "int32", 6: "float32",
#                          7: "float64", 10: "complex64", 11: "complex128"}
#         print('GDAL data type  =', out_band.DataType)
#         print('Numpy data type =', gdal_type_map[ out_band.DataType ])

        #------------------------------------------
        # Print out all the non-method attributes
        # to find:  DataType, XSize, & YSize)
        #------------------------------------------
#         for att in dir(out_band):
#             if not(callable(getattr(out_band, att))):
#                 print( att, getattr(out_band, att) )        
#         print()

        #------------------------------
        # Change the DEM nodata value
        #------------------------------
        ## print('###### Nodata after warp =', VRT_nodata)
        w1 = (dem == VRT_nodata)   # boolean array
        w2 = np.invert( w1 )
        out_nodata = -9999.0
        # out_nodata = np.nan
        #----------------------------
        dem       = (dem / zfactor)   # (convert cm to m)
        dem[ w1 ] = out_nodata
        dem = np.float32( dem )   #############
        dem_dtype = dem.dtype  # (orig: int32, then: float64, then: float32)
        #------------------------------
        dem_min   = np.min( dem[w2] )
        dem_max   = np.max( dem[w2] )

        #-----------------------------------------
        # Write rescaled DEM to new GeoTIFF file
        #-----------------------------------------
        if (SAVE_DEM):
            ncols  = out_ncols
            nrows  = out_nrows
            nbands = 1
            data_type = gdal.GDT_Float32
            #------------------------------------------  
            driver    = gdal.GetDriverByName('GTiff')
            outRaster = driver.Create(DEM_out_file, ncols, nrows, nbands, data_type)
            outRaster.SetGeoTransform( out_gt )
            outband = outRaster.GetRasterBand(1)
            outband.WriteArray( dem )
            outband.SetNoDataValue( out_nodata )
            outband.SetUnitType( 'm' )   ######
            outRasterSRS = osr.SpatialReference()
            if (CLIP_ALBERS):
                outRasterSRS.ImportFromWkt( wkt_albers )
            else:
                epsg_code = 4326  # (WGS 84)
                outRasterSRS.ImportFromEPSG( epsg_code )
            outRaster.SetProjection( outRasterSRS.ExportToWkt() )
            outband.FlushCache()
       
        #------------------------
        # Save DEM info in self
        #------------------------
        self.DEM_bounds = out_bounds
        self.DEM_ncols  = out_ncols
        self.DEM_nrows  = out_nrows
        self.DEM_xres   = out_xres
        self.DEM_yres   = out_yres
        self.DEM_min    = dem_min
        self.DEM_max    = dem_max
         
        if (REPORT):
            print('Info for ' + cat_id_str + ' DEM:')
            print('  ncols  =', out_ncols)
            print('  nrows  =', out_nrows)
            print('  xres   =', out_xres)
            print('  yres   =', out_yres)
            print('  xy_units  =', out_xy_units)
            print('  DEM units =', 'meters')
            print('  DEM dtype =', dem_dtype)
            print('  nodata    =', out_nodata)
            print('  min, max  =', dem_min, ',', dem_max)
            print()

        #----------------------------------        
        # Close both in_file and out_file
        #----------------------------------
        print('Closing files...')
        vrt_unit = None   # Close vrt_file
        out_unit = None   # Close out_file
        print('Finished.')
        print() 
        return dem
                           
    #   get_dem()
    #--------------------------------------------------------------------
    def get_outlet_coords(self, cat_id_str, REPORT=False):

        #--------------------------------------------------------------
        # Notes:  Get outlet coordinates for a given catchment from
        #         its id string, e.g. id_str = 'cat-29'.
        #         First, get value of "toid" string stored in:
        #            catchment-data.geojson, e.g. 'nex-26' for cat-29
        #         Then, get coords from 'geometry' (POINT type)   
        #--------------------------------------------------------------
        crs_name = self.catchment_data['crs']['properties']['name']
        
        df1 = self.catchment_data['features']    # list of dict 
        record1 = next((item for item in df1 if item["id"] == cat_id_str), False)
        if not(record1):
            print('SORRY, cat_id not found:', cat_id_str)
            print('Returning.')
            print()
            return -1
        nexus_id_str = record1['properties']['toid']
        self.nexus_id_str = nexus_id_str

        df2 = self.nexus_data['features']        # list of dict
        record2 = next((item for item in df2 if item["id"] == nexus_id_str), False)
        if not(record2):
            print('SORRY, nexus_id not found:', nexus_id_str)
            print('Returning.')
            print()
            return -1

        nexus_type = record2['properties']['nexus_type']
        if (nexus_type != 'junction'):
            print( 'WARNING: nexus_type =', nexus_type )

        rgc = record2['geometry']['coordinates']
        coords = np.array( rgc , dtype='float64' )
        self.outlet_x = coords[0]
        self.outlet_y = coords[1]

        if (crs_name.endswith('CRS84')):
            self.outlet_lon = self.outlet_x
            self.outlet_lat = self.outlet_y
            crs_type = 'Geographic'
        elif (crs_name.endswith('5070')):
            crs_type = 'Albers Equal Area Conic'
        else:
            crs_type = 'Unknown'

        if (REPORT):
            print('Outlet coords for: ' + cat_id_str)
            print('  x, y =', self.outlet_x, ',', self.outlet_y )
            print('  crs_name =', crs_name )
            print('  crs_type =', crs_type)
            print()
            
        return coords
        
    #   get_outlet_coords()
    #--------------------------------------------------------------------
    def get_headwater_catchments(self, REPORT=False):

        #----------------------    
        # Get list of "toids"
        #----------------------
        # features -> properties -> toid
        df = self.catchment_data['features']    # list of dict
        toid_num_list = list() 
        for record in df:
            toid_str = record['properties']['toid']
            num_str = toid_str[4:]  # remove "nex-"
            toid_num_list.append( num_str )

        #-----------------------------------------------       
        # Create list of all "ids" not in "toids" list
        #-----------------------------------------------
        head_ids = list()
        for record in df:
            id_str  = record['id']
            num_str = id_str[4:]  # remove "cat-"
            if (num_str not in toid_num_list):
                head_id_str = 'cat-' + num_str
                head_ids.append( head_id_str )

        #---------------------------------------------
        # Option to sort these by their "area_sqkm",
        # and perhaps remove the smallest ones
        #---------------------------------------------
        n_heads = len(head_ids)
        areas = np.zeros( n_heads, dtype='float64' )
        k = 0
        for cat_id_str in head_ids:
            areas[ k ] = self.get_area_sqkm( cat_id_str )
            k += 1
        s = np.argsort( -1 * areas )   # -1 => large to small
#         amin = a.min()
#         amax = a.max()
        head_ids_arr = np.array( head_ids )
        head_ids_arr = head_ids_arr[ s ]
        head_ids     = head_ids_arr.tolist()
        
        if (REPORT):
            for cat_id_str in head_ids:
                print(cat_id_str, self.get_area_sqkm(cat_id_str))

        self.head_ids = head_ids
        return head_ids
        
    #   get_headwater_catchments()
    #--------------------------------------------------------------------    
    def get_area_sqkm(self, cat_id_str, REPORT=False):    

        crs_name = self.catchment_data['crs']['properties']['name']
        df = self.catchment_data['features']
        record = next((item for item in df if item["id"] == cat_id_str), False)
        if not(record):
            print('SORRY, cat_id not found:', cat_id_str)
            print('Returning.')
            print()
            return -1       

        self.area_sqkm = np.float64( record['properties']['area_sqkm'] )

        if (crs_name.endswith('CRS84')):
            crs_type = 'Geographic'
        elif (crs_name.endswith('5070')):
            crs_type = 'Albers Equal Area Conic'
        else:
            crs_type = 'Unknown'
            
        if (REPORT):
            print('Area for: ' + cat_id_str + ' [km2] =')
            print('  ', self.area_sqkm )
            print('  crs_name       =', crs_name )
            print('  crs_type       =', crs_type )
            print()
                       
        return self.area_sqkm
             
    #   get_area_sqkm()
    #--------------------------------------------------------------------
    def ocean_outlet_count(self, REPORT=False):
 
        df    = self.catchment_data['features']     
        recs  = [item for item in df if item['properties']['toid'] == 'nex-0']
        count = len(recs)
        self.ocean_outlet_count = count
        if (REPORT):
            print('Number of catchments draining to ocean =', count)
            print()  
        return count

    #   ocean_outlet_count()
    #--------------------------------------------------------------------     
    def print_info(self, cat_id_str='cat-29'):
   
        print('Extracting catchment info...') 
        area_sqkm = self.get_area_sqkm( cat_id_str )
        if (area_sqkm == -1):
            print('Catchment not found:', cat_id_str)
            return

        outlet_coords = self.get_outlet_coords( cat_id_str )
        bounds        = self.get_bounding_box( cat_id_str )
         
        print('Info for catchment: ' + cat_id_str)
        print('   catchment ID   = ', cat_id_str )
        print('   outlet lon,lat = ', outlet_coords )
        print('   area_sqkm =', area_sqkm, '[km2]' )
        print('   minlon    =', bounds[0] )
        print('   minlat    =', bounds[1] )
        print('   maxlon    =', bounds[2] )
        print('   maxlat    =', bounds[3] )
        print('   nexus_ID  =', self.nexus_id_str )
        print()
        
    #   print_info()
    #--------------------------------------------------------------------
    
#---------------------------------------------------------------------
def get_hydrofab_dem_info():

    #------------------------------------
    # bounds = [ulx, lry, lrx, uly]
    #        = [xmin, ymin, xmax, ymax]
    #------------------------------------
    dem_info = dict()
    dem_info['ncols'] = 163008
    dem_info['nrows'] = 114503
    dem_info['dx']    = 30.0
    dem_info['dy']    = 30.0
    dem_info['xmin']  = -2470965.0   # ulx
    dem_info['xmax']  =  2419275.0   # lrx
    dem_info['ymin']  =  186285.0    # lry     
    dem_info['ymax']  =  3621375.0   # uly 
    return dem_info
    
#   get_hydrofab_dem_info()
#------------------------------------------------------------------------
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

    ulx = -2470965  # -2.4709650000000005e+06
    uly = 3621375   # 3.6213750000000019e+06
    
    if (REPORT):
        print('lon, lat =', lon_deg, ',', lat_deg)
        print('x,   y   =', x, ',', y)
        print()
        print('x - ulx =', x - ulx)
        print('uly - y =', uly - y)

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
    d2r  = (np.pi / 180.)  # (degrees -> radians)
    lat  = lat_deg * d2r
    lon  = lon_deg * d2r
    
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

    dem_info = get_hydrofab_dem_info()
    ncols = dem_info['ncols']
    nrows = dem_info['nrows']
    xres  = dem_info['dx']
    yres  = dem_info['dy']
    ulx   = dem_info['xmin']
    uly   = dem_info['ymax']
    lry   = dem_info['ymin']

    col = (x - ulx) / xres
    row = (y - lry) / yres
    ## row = (y - uly) / yres
    # row = (nrows - 1) - (y - uly)/yres
    
    if (REPORT):
        print('lon, lat =', lon_deg, ',', lat_deg)
        print('x,   y   =', x, ',', y)
        print('ulx, uly =', ulx, ',', uly)
        print('lry =', lry)
        print('ncols, nrows =', ncols, nrows)
        print('col, row     =', col, row)
        print()
        print('x - ulx =', x - ulx)
        print('uly - y =', uly - y)

    return (x,y)
        
#   Albers_XY2()
#------------------------------------------------------------------------

      
       
       
       
       
#  class catchment()
#      __init__()
#      read_data()
#      get_bounding_box()
#      get_raster_cellsize()   #####
#      get_raster_bounds()     #####
#      get_dem()
#      get_outlet_coords()
#      get_headwater_catchments()  # 11/30/22
#      get_area_sqkm()
#      ocean_outlet_count()
#      print_info()
#
#  get_hydrofab_dem_info()
#  Albers_XY()  # (function)
#  Albers_q()
#  Albers_XY2()
#
#------------------------------------------------------------------------
class catchment:
    #--------------------------------------------------------------------
    def __init__(self, catchment_data_file=None, nexus_data_file=None,
                       spatial_dir=None):

        # Rename data_file to json_file ??
        if (spatial_dir is None):
            # Relative to the top level of ngen repo folder.
            # spatial_dir = './data/topoflow/spatial/'
            spatial_dir = '/Users/peckhams/Dropbox/GitHub/ngen/data/topoflow/spatial/'
            ## spatial_dir = './'
        if (catchment_data_file is None):
            catchment_data_file = 'catchment_data_HUC01.geojson'
        if (nexus_data_file is None):
            nexus_data_file     = 'nexus_data_HUC01.geojson'
        self.catchment_data_file = spatial_dir + catchment_data_file
        self.nexus_data_file     = spatial_dir + nexus_data_file
        #------------------
        self.read_data()
    
    #   __init__()    
    #--------------------------------------------------------------------
    def read_data(self):

        #--------------------------------------------------------------     
        # Note: This version reads the GeoJSON files directly without
        #       adding geopandas as a dependency.  It reads the data
        #       into a "simulated dataframe".
        #--------------------------------------------------------------
        print('Reading catchment info files...')
        cat_unit = open( self.catchment_data_file, 'r')
        nex_unit = open( self.nexus_data_file, 'r')
        self.catchment_data = json.load( cat_unit )
        self.nexus_data     = json.load( nex_unit )
        cat_unit.close()
        nex_unit.close()
        nc = len(self.catchment_data['features'])
        print('Read info for', nc, 'catchments.')
        print()
                       
    #   read_data()   
    #--------------------------------------------------------------------
    def get_bounding_box(self, cat_id_str, buffer_sec=0, REPORT=False):

        #--------------------------------------------------------------
        # Notes:  Get bounding box coords for a given catchment from
        #         its id string, e.g. id_str = 'cat-29'.
        #         Polygon coordinates are stored in:
        #            catchment-data.geojson (POLYGON type geometry)
        #         Coordinates can be given as Geographic or Albers.
        #--------------------------------------------------------------   
        df = self.catchment_data['features']  # returns list of dictionaries
        crs_name = self.catchment_data['crs']['properties']['name']
        record = next((item for item in df if item["id"] == cat_id_str), False)
        if not(record):
            print('SORRY, cat_id not found:', cat_id_str)
            print('Returning.')
            print()
            return -1

        rgc = record['geometry']['coordinates'][0]
        pairs   = np.array( rgc, dtype='float64')
        xcoords = pairs[:,0]
        ycoords = pairs[:,1]
        self.minlon = xcoords.min()
        self.minlat = ycoords.min()
        self.maxlon = xcoords.max()
        self.maxlat = ycoords.max()
        bounds = [self.minlon, self.minlat, self.maxlon, self.maxlat] 
                
        # print('## record["id"] =', record['id'])
        # print('## record["geometry"]["coordinates"][0] =', rgc)
        # print('## type(rgc) =', type(rgc) )  # list
        # print('## rgc[0] =', rgc[0])
        # print('## xcoords[0:2] =', xcoords[0:2])
        # print('## ycoords[0:2] =', ycoords[0:2])
        # print()

        #------------------------------------------------
        # Option to specify a buffer in arcseconds
        #------------------------------------------------
        # (3 arcsecs) * (1 deg / 3600 arcsecs) =
        # 0.000833333 degrees  ~ 92.6 meters
        #------------------------------------------------
        # (6 arcsecs) * (1 deg / 3600 arcsecs) =
        # 0.001666666 degrees  ~ 185.2 meters        
        #------------------------------------------------
        if (buffer_sec > 0):
           buffer_deg = (buffer_sec / 3600.0)
           self.minlon -= buffer_deg
           self.minlat -= buffer_deg
           self.maxlon += buffer_deg
           self.maxlat += buffer_deg
           bounds = [self.minlon, self.minlat, self.maxlon, self.maxlat]

        #---------------------------------------------------------
        # Option to specify a buffer in decimal degrees, as a %.
        #---------------------------------------------------------
#         if (buffer_pct > 0):
#             lon_range = (self.maxlon - self.minlon)
#             lat_range = (self.maxlat - self.minlat)
#             dx = buffer_pct * lon_range
#             dy = buffer_pct * lat_range
#             self.minlon -= dx
#             self.minlat -= dy
#             self.maxlon += dx
#             self.maxlat += dy
#             bounds = [self.minlon, self.minlat, self.maxlon, self.maxlat]
     
        if (REPORT):
            print('Bounding box for: ' + cat_id_str)
            print('  minlon, maxlon =', bounds[0], ', ', bounds[2])
            print('  minlat, maxlat =', bounds[1], ', ', bounds[3])
            print('  crs_name       =', crs_name )
            print()

        self.bounds = bounds   
        return bounds
        
    #   get_bounding_box()
    #-------------------------------------------------------------------
    def get_raster_cellsize( self, gdal_unit ):

        # Note:  This is from regrid.py.
        geotransform = gdal_unit.GetGeoTransform()
        # ulx  = geotransform[0]
        xres = geotransform[1]
        # xrtn = geotransform[2]
        #-----------------------
        # uly  = geotransform[3]
        # yrtn = geotransform[4]  # (not yres !!)
        yres = geotransform[5]  # (not yrtn !!)
    
        return (xres, abs(yres))

    #   get_raster_cellsize()
    #------------------------------------------------------------------- 
    def get_raster_bounds( self, ds, VERBOSE=False):

        #-------------------------------------------------------------
        # Note:  The bounds depend on the map projection and are not
        # necessarily a Geographic bounding box of lons and lats.
        # ds = dataset = "gdal_unit"   ########  
        #-------------------------------------------------------------
        # See:
        # https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html
        # and search on "geotransform".  An example of gdal.SetGeoTransform
        # gives: [xmin, pixel_size, 0, ymax, 0, -pixel_size].
        # Also says args are:
        # [ulx, xDist, rtnX, uly, yDist, rtnY]
        # This is consistent with information below.
        #-------------------------------------------------------------    
        # ulx = upper left x  = xmin
        # uly = upper left y  = ymax
        # lrx = lower right x = xmax
        # lry = lower right y = ymin
        #-----------------------------

        #----------------------------------------------------------  
        # Notice the strange order of parameters here is CORRECT.
        # It is not:  ulx, xres, xskew, uly, yres, yskew
        #----------------------------------------------------------
        ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
        lrx = ulx + (ds.RasterXSize * xres)
        lry = uly + (ds.RasterYSize * yres)

        if (VERBOSE):
            print('ulx, uly   =', ulx, uly)
            print('lrx, lry   =', lrx, lry)
            print('xres, yres = ', xres, yres)
            print('xskew, yskew =', xskew, yskew)
            print('----------------------------------')

        return [ulx, lry, lrx, uly]  # [xmin, ymin, xmax, ymax]

    #   get_raster_bounds()
    #--------------------------------------------------------------------
    def get_dem( self, cat_id_str='cat-29', RESAMPLE_ALGO='bilinear',
                 CLIP_ALBERS=False, CLIP_GEO=False,
                 zfactor=100, REPORT=True,
                 SAVE_DEM=True, DEM_out_file=None):

        #------------------------------------------------------------
        # Note:  Use GDAL VRT format to get DEMs for catchments.
        #------------------------------------------------------------
        # Note:  The CLIP_GEO option will return a DEM that uses
        #        Geographic coordinates and has fixed-angle grid
        #        cell sizes (e.g. 2 arcseconds).
        #        The CLIP_ALBERS option will return a DEM that uses
        #        Albers Conic Equal Area coordinates and has fixed-
        #        length grid cell sizes (e.g. 100 meters).
        #        The original hydrofabric DEM uses this same Albers
        #        projection, but the catchment polygons described
        #        in the catchment & nexus GeoJSON files use
        #        Geographic coordinates.
        #------------------------------------------------------------             
        if not(CLIP_ALBERS or CLIP_GEO):
            CLIP_ALBERS=True
        if (SAVE_DEM and (DEM_out_file is None)):
            DEM_out_file = cat_id_str + '_rawDEM.tif'
            # DEM_out_file = cat_id_str + '_DEM.tif'
           
        #--------------------------------------------------------  
        # This shows some of the other keywords to gdal.Warp.
        #--------------------------------------------------------      
        # WarpOptions(options=[], format=None, outputBounds=None,
        # outputBoundsSRS=None, xRes=None, yRes=None, targetAlignedPixels=False,
        # width=0, height=0, srcSRS=None, dstSRS=None, srcAlpha=False,
        # dstAlpha=False, warpOptions=None, errorThreshold=None,
        # warpMemoryLimit=None, creationOptions=None, outputType=GDT_Unknown,
        # workingType=GDT_Unknown, resampleAlg=None, srcNodata=None,
        # dstNodata=None, multithread=False, tps=False, rpc=False,
        # geoloc=False, polynomialOrder=None, transformerOptions=None,
        # cutlineDSName=None, cutlineLayer=None, cutlineWhere=None,
        # cutlineSQL=None, cutlineBlend=None, cropToCutline=False,
        # copyMetadata=True, metadataConflictValue=None,
        # setColorInterpretation=False, callback=None, callback_data=None)
        # 
        # Create a WarpOptions() object that can be passed to gdal.Warp()
        # Keyword arguments are : options --- can be be an array of strings,
        # a string or let empty and filled from other keywords.
                    
        #-----------------------------------   
        # Specify the resampling algorithm
        #-----------------------------------
        algo_dict = {
        'nearest'     : gdal.GRA_NearestNeighbour,
        'bilinear'    : gdal.GRA_Bilinear,
        'cubic'       : gdal.GRA_Cubic,
        'cubicspline' : gdal.GRA_CubicSpline,
        'lanczos'     : gdal.GRA_Lanczos,
        'average'     : gdal.GRA_Average,
        'min'         : gdal.GRA_Min,
        'max'         : gdal.GRA_Max,
        'mode'        : gdal.GRA_Mode,
        'med'         : gdal.GRA_Med }
    
        resample_algo = algo_dict[ RESAMPLE_ALGO ]

        #----------------------------------------------------    
        # Set config options for accessing an AWS S3 bucket?
        #-----------------------------------------------------------------
        # If you have installed the AWS CLI for your OS, and have run
        # the command "aws configure", then these do not need to be set.
        # In that case, GDAL will find them in ~/.aws/credentials and
        # ~/.aws/config.
        #-----------------------------------------------------------------
        # print('Setting GDAL config options...')
        # gdal.SetConfigOption('AWS_REGION', 'us-east-1')
        # gdal.SetConfigOption('AWS_SECRET_ACCESS_KEY', 'MY_SECRET_ACCESS_KEY')
        # gdal.SetConfigOption('AWS_ACCESS_KEY_ID', 'MY_ACCESS_KEY_ID')

        #-----------------------------------------------------------
        # vsis3 = vsi s3'. VSI = Virtual Server Infrastructure
        #-----------------------------------------------------------
        # VRT = Virtual Raster Table (a GDAL format)
        # A VRT file is a type of XML file.
        # From Mike: vrt_file = "s3://nextgen-hydrofabric/dem.vrt"
        # I think this is considered a "public" S3 bucket.
        # The projection is: Albers Conic Equal Area and spans the
        # entire US with 30 meter grid cell size.
        # Mike says that elevation units are "cm", but this is not
        #   specified in the VRT metadata. So must divide by 100.
        #-----------------------------------------------------------
        print('Opening VRT file in S3 bucket...')
        bucket_name = 'nextgen-hydrofabric/'   # (AWS S3 bucket name)
        dem_file = 'dem.vrt'
        vrt_path = '/vsis3/' + bucket_name + dem_file
        vrt_unit = gdal.Open( vrt_path )

        #--------------------------------------------------------
        # Get info about the full DEM without downloading it
        # (Should be able to get anything in the VRT XML file.)
        #--------------------------------------------------------
        (dx, dy) = self.get_raster_cellsize( vrt_unit )  #########
        vrt_xres = dx
        vrt_yres = dy
        # This block assumes xres is Degrees (it isn't here)
        ## vrt_xres_deg = dx
        ## vrt_yres_deg = dy
        ## vrt_xres_sec = (vrt_xres_deg * 3600.0)
        ## vrt_yres_sec = (vrt_yres_deg * 3600.0)
        vrt_ncols    = vrt_unit.RasterXSize
        vrt_nrows    = vrt_unit.RasterYSize
        # geotransform = ulx, xres, xskew, uly, yskew, yres
        vrt_gt       = vrt_unit.GetGeoTransform()   # (tuple)
        vrt_ulx      = vrt_gt[0]
        vrt_uly      = vrt_gt[3]
        vrt_lry      = vrt_uly - (vrt_nrows * vrt_yres)   #########
        vrt_bounds   = self.get_raster_bounds( vrt_unit )
        vrt_proj     = vrt_unit.GetProjection()
        vrt_meta     = vrt_unit.GetMetadata()   # (can be empty dict)
        vrt_srs      = osr.SpatialReference(wkt=vrt_proj)
        vrt_geog_cs  = vrt_srs.GetAttrValue('geogcs')

        if (vrt_srs.IsProjected):
            vrt_proj_cs = vrt_srs.GetAttrValue('projcs')
        else:
            vrt_proj_cs = 'None'
        #-----------------------------------------------------
        # WARNING: This may cause full DEM to be downloaded.
        #-----------------------------------------------------
        # vrt_band   = vrt_unit.GetRasterBand(1)
        # vrt_nodata = vrt_band.GetNoDataValue()
        # zmin,zmax,zmean,zstdv = vrt_band.GetStatistics(True, True)
        
        ## vrt_dem    = vrt_unit.ReadAsArray()
        ## dem_dtype  = vrt_dem.dtype
        ## zmin       = vrt_dem.min()
        ## zmax       = vrt_dem.max()        
        #----------------------------------------------------     
        if (REPORT):
            print('VRT DEM Information:')
            print('  ncols    =', vrt_ncols)
            print('  nrows    =', vrt_nrows)
            print('  xres     =', vrt_xres)
            print('  yres     =', vrt_yres)
            print('  xy_units  =', 'meters')
            print('  DEM units =', 'cm')     # not saved in metadata!
            print('  DEM dtype =', 'int32')  # long, 4-byte int
            print('  bounds   =', vrt_bounds)
            print('  geog cs  =', vrt_geog_cs)
            print('  proj cs  =', vrt_proj_cs)
            print('  metadata =', vrt_meta)
            ##### print('  proj     =', vrt_proj)    # WKT, wordy
            ## print('zmin    =', zmin)
            ## print('zmax    =', zmax)
            ## print('VRT DEM dtype =', dem_dtype)
            ## print('VRT DEM min, max  =', dem_min, ', ', dem_max)
            print()

        #---------------------------------------------
        # Check the VRT bounding box info (all good)
        #---------------------------------------------
#         vrt_xrange = vrt_bounds[2] - vrt_bounds[0]
#         vrt_yrange = vrt_bounds[3] - vrt_bounds[1]
#         print('vrt_xrange =', vrt_xrange)
#         print('  ncols * xres =', vrt_ncols * vrt_xres)
#         print('vrt_yrange =', vrt_yrange)
#         print('  nrows * yres =', vrt_nrows * vrt_yres)
#         print()
        
        #------------------------------------------
        #  GDAL bounds  = [ulx, lry, lrx, uly] 
        #               = [xmin, ymin, xmax, ymax]
        #---------------------------------------------------------
        # Note:  get_bounding_box() uses the catchment & nexus
        #        GeoJSON files and returns a geographic bounding
        #        box: minlon, minlat, maxlon, maxlat
        #---------------------------------------------------------           
        out_bounds_geo = self.get_bounding_box( cat_id_str, REPORT=True )

        #--------------------------------------------------------
        # Set the output Spatial Reference System via EPSG code
        #--------------------------------------------------------
        # The NextGen Hydrofabric SRC appears to be unique and
        # does not have an EPSG code.  However, this one seems
        # to be a somewhat close match (EPSG: 5070)        
        #
        # Also see:
        # https://desktop.arcgis.com/en/arcmap/10.3/guide-books/
        #         map-projections/albers-equal-area-conic.htm
        # https://spatialreference.org/ref/esri/102003/html/
        # https://en.wikipedia.org/wiki/Albers_projection
        #----------------------------------------------------
        # EPSG Codes via https://epsg.io/:
        #
        # 4326: WGS 84 - WGS84 - World Geodetic System 1984
        # 5070: NAD83 / Conus Albers
        # 7019: GRS 1980
        # 6269: MAGNA-SIRGAS / San Jose del Guaviare urban grid
        # 9122: degree (supplier to define representation)
        # 4269: NAD83
        # 9001: metre
        #----------------------------------------------------
        # Convert coords from lat/lon to Hydrofabric Albers
        #----------------------------------------------------
        wgs84 = osr.SpatialReference()
        wgs84.ImportFromEPSG( 4326 )
        wkt_albers = vrt_unit.GetProjection()
        # wkt = wkt.replace('unnamed', 'Albers')  ####
        albers = osr.SpatialReference()
        albers.ImportFromWkt( wkt_albers )
        ### albers.ImportFromEPSG( 5070 )   #####
       
        #-------------------------------------------------
        # Option to check the WKT (well known text) info
        #-------------------------------------------------
        # See: https://www.geoapi.org/3.0/javadoc/org/
        #      opengis/referencing/doc-files/WKT.html
        #-------------------------------------------------
#         wkt1 = wgs84.ExportToWkt()        
#         wkt1 = wgs84.ExportToPrettyWkt()
#         print('##### wkt for WGS 84 =')
#         print( wkt1 )
#         print()
#         wkt2 = albers.ExportToPrettyWkt()
#         print('##### wkt2 for Albers =')
#         print( wkt2 )
#         print()

        #--------------------------------------------------
        # Option to clip using original Albers projection
        # or the similar EPSG:5070 projection
        #--------------------------------------------------
        if (CLIP_ALBERS):  
            #-----------------------------------------------------
            # Transform geographic bounding box to Albers extent
            #---------------------------------------------------------------
            albers2latlon = osr.CoordinateTransformation( albers, wgs84 )
            latlon2albers = osr.CoordinateTransformation( wgs84, albers )
            #---------------------------------------------------------------
            minlon = out_bounds_geo[0]
            minlat = out_bounds_geo[1]
            maxlon = out_bounds_geo[2]
            maxlat = out_bounds_geo[3]
            #-----------------------------------------------------------
            # NOTE:  When using the GDAL Python package version 3.*
            #        the order of arguments to TransformPoint must be
            #        flipped (i.e. y,x vs. x,y) as compared to version
            #        2.*.  This means that most online examples (that
            #        use version 2.*) are wrong.
            #        See:  https://github.com/OSGeo/gdal/issues/1546
            #-----------------------------------------------------------
            albers_ul_xyz = latlon2albers.TransformPoint( maxlat, minlon )
            albers_ll_xyz = latlon2albers.TransformPoint( minlat, minlon )
            albers_ur_xyz = latlon2albers.TransformPoint( maxlat, maxlon )
            albers_lr_xyz = latlon2albers.TransformPoint( minlat, maxlon )
            print('Albers UL xyz:', albers_ul_xyz)
            print('Albers LL xyz:', albers_ll_xyz)
            print('Albers UR xyz:', albers_ur_xyz)
            print('Albers LR xyz:', albers_lr_xyz)
            print()
            xvals = [albers_ul_xyz[0], albers_ll_xyz[0],
                     albers_ur_xyz[0], albers_lr_xyz[0]]
            yvals = [albers_ul_xyz[1], albers_ll_xyz[1],
                     albers_ur_xyz[1], albers_lr_xyz[1]]
            xmin = np.floor( np.min( xvals ) )
            xmax = np.ceil( np.max( xvals ) )
            ymin = np.floor( np.min( yvals ) )
            ymax = np.ceil( np.max( yvals ) )
            #----------------------------------------------
            out_srs = None
            # out_srs = "EPSG:5070"
            out_bounds_albers = [xmin, ymin, xmax, ymax]
            out_bounds   = out_bounds_albers
            out_xy_units = 'meters'
            #-----------------------------------------------------
            # For hydrofabric DEM (dem.vrt), xres = yres = 30 m,
            # which is roughly 1 arcsecond at the equator.
            #-----------------------------------------------------
            out_xres = 90.0    # [meters]   ########
            out_yres = 90.0    # [meters]   ########

            #------------------------------------------
            # Option to compute new DEM cols and rows
            #   relative to original VRT DEM
            #------------------------------------------
#             col_min = np.int32((xmin - vrt_ulx) / vrt_xres)
#             col_max = np.int32((xmax - vrt_ulx) / vrt_xres)
#             row_min = np.int32((ymin - vrt_lry) / vrt_yres)
#             row_max = np.int32((ymax - vrt_lry) / vrt_yres)
#     #         row_min = np.int32((ymin - vrt_uly) / vrt_yres)
#     #         row_max = np.int32((ymax - vrt_uly) / vrt_yres)
#             print('col_min, col_max =', col_min, col_max)
#             print('row_min, row_max =', row_min, row_max)
#             print()
     
        #--------------------------------------------------------
        # Option to clip using Geographic bounding box (WGS 84)
        # In gdal.Warp, must set: outputBounds=out_bounds
        # and check other keywords like: xRes, yRes.
        #--------------------------------------------------------
        if (CLIP_GEO):      
            out_srs = "EPSG:4326"
            out_bounds = out_bounds_geo
            out_xy_units = 'degrees'
            # If dstSRS = WGS84, out_xres units must be degrees
            out_xres = 3.0/3600  # 3 arcseconds -> degrees
            out_yres = 3.0/3600

        #----------------------------------------
        # Change/increase the grid cell size ??
        #----------------------------------------
        # out_xres = vrt_xres   # no change
        # out_yres = vrt_yres
        #-------------------------
        # out_xres = (2 * vrt_xres)
        # out_yres = (2 * vrt_yres)
        #----------------------------------------------------------
        # For hydrofabric DEM (dem.vrt), xres = yres = 30 meters,
        # which is roughly 1 arcsecond at the equator.
        # If dstSRS = WGS84, out_xres units must be degrees
        #----------------------------------------------------------        
        # out_xres = 2./3600  # (2 arcseconds -> degrees)
        # out_yres = 2./3600  # (2 arcseconds -> degrees)

        #------------------------------------------------  
        # Resample & clip and write new grid to GeoTIFF
        #------------------------------------------------
        # See notes near "METHOD 2" below.
        #------------------------------------------------
        # Only set xRes and yRes if they are Degrees
        #------------------------------------------------        
        print('Creating DEM for ' + cat_id_str + '...')     
        out_unit = gdal.Warp( DEM_out_file, vrt_unit,
            format = None,  # (don't save to file yet)
            #### format = 'GTiff',  # (output format string)
            #### outputType=gdal.GDT_Float32,
            outputBounds=out_bounds,  #######
            dstSRS=out_srs,  #######
            ### srcSRS="EPSG:####",  # don't need to specify
            xRes=out_xres, yRes=out_yres,  #################
            resampleAlg = resample_algo )

        #----------------------
        # Get info on new DEM
        #----------------------
        ## print('Getting new DEM info...')
        out_ncols  = out_unit.RasterXSize
        out_nrows  = out_unit.RasterYSize
        (out_xres, out_yres) = self.get_raster_cellsize( out_unit )  #####
        ## out_bounds = self.get_raster_bounds( out_unit )
        out_band   = out_unit.GetRasterBand(1)
        VRT_nodata = out_band.GetNoDataValue()
        out_gt     = out_unit.GetGeoTransform()   # (tuple)
        dem        = out_unit.ReadAsArray()

        #-------------------------------------------
        # Convert GDAL DataType to Numpy data type
        #-------------------------------------------
#         gdal_type_map = {1: "int8", 2: "uint16", 3: "int16",
#                          4: "uint32", 5: "int32", 6: "float32",
#                          7: "float64", 10: "complex64", 11: "complex128"}
#         print('GDAL data type  =', out_band.DataType)
#         print('Numpy data type =', gdal_type_map[ out_band.DataType ])

        #------------------------------------------
        # Print out all the non-method attributes
        # to find:  DataType, XSize, & YSize)
        #------------------------------------------
#         for att in dir(out_band):
#             if not(callable(getattr(out_band, att))):
#                 print( att, getattr(out_band, att) )        
#         print()

        #------------------------------
        # Change the DEM nodata value
        #------------------------------
        ## print('###### Nodata after warp =', VRT_nodata)
        w1 = (dem == VRT_nodata)   # boolean array
        w2 = np.invert( w1 )
        out_nodata = -9999.0
        # out_nodata = np.nan
        #----------------------------
        dem       = (dem / zfactor)   # (convert cm to m)
        dem[ w1 ] = out_nodata
        dem = np.float32( dem )   #############
        dem_dtype = dem.dtype  # (orig: int32, then: float64, then: float32)
        #------------------------------
        dem_min   = np.min( dem[w2] )
        dem_max   = np.max( dem[w2] )

        #-----------------------------------------
        # Write rescaled DEM to new GeoTIFF file
        #-----------------------------------------
        if (SAVE_DEM):
            ncols  = out_ncols
            nrows  = out_nrows
            nbands = 1
            data_type = gdal.GDT_Float32
            #------------------------------------------  
            driver    = gdal.GetDriverByName('GTiff')
            outRaster = driver.Create(DEM_out_file, ncols, nrows, nbands, data_type)
            outRaster.SetGeoTransform( out_gt )
            outband = outRaster.GetRasterBand(1)
            outband.WriteArray( dem )
            outband.SetNoDataValue( out_nodata )
            outband.SetUnitType( 'm' )   ######
            outRasterSRS = osr.SpatialReference()
            if (CLIP_ALBERS):
                outRasterSRS.ImportFromWkt( wkt_albers )
            else:
                epsg_code = 4326  # (WGS 84)
                outRasterSRS.ImportFromEPSG( epsg_code )
            outRaster.SetProjection( outRasterSRS.ExportToWkt() )
            outband.FlushCache()
       
        #------------------------
        # Save DEM info in self
        #------------------------
        self.DEM_bounds = out_bounds
        self.DEM_ncols  = out_ncols
        self.DEM_nrows  = out_nrows
        self.DEM_xres   = out_xres
        self.DEM_yres   = out_yres
        self.DEM_min    = dem_min
        self.DEM_max    = dem_max
         
        if (REPORT):
            print('Info for ' + cat_id_str + ' DEM:')
            print('  ncols  =', out_ncols)
            print('  nrows  =', out_nrows)
            print('  xres   =', out_xres)
            print('  yres   =', out_yres)
            print('  xy_units  =', out_xy_units)
            print('  DEM units =', 'meters')
            print('  DEM dtype =', dem_dtype)
            print('  nodata    =', out_nodata)
            print('  min, max  =', dem_min, ',', dem_max)
            print()

        #----------------------------------        
        # Close both in_file and out_file
        #----------------------------------
        print('Closing files...')
        vrt_unit = None   # Close vrt_file
        out_unit = None   # Close out_file
        print('Finished.')
        print() 
        return dem
                           
    #   get_dem()
    #--------------------------------------------------------------------
    def get_outlet_coords(self, cat_id_str, REPORT=False):

        #--------------------------------------------------------------
        # Notes:  Get outlet coordinates for a given catchment from
        #         its id string, e.g. id_str = 'cat-29'.
        #         First, get value of "toid" string stored in:
        #            catchment-data.geojson, e.g. 'nex-26' for cat-29
        #         Then, get coords from 'geometry' (POINT type)   
        #--------------------------------------------------------------
        crs_name = self.catchment_data['crs']['properties']['name']
        
        df1 = self.catchment_data['features']    # list of dict 
        record1 = next((item for item in df1 if item["id"] == cat_id_str), False)
        if not(record1):
            print('SORRY, cat_id not found:', cat_id_str)
            print('Returning.')
            print()
            return -1
        nexus_id_str = record1['properties']['toid']
        self.nexus_id_str = nexus_id_str

        df2 = self.nexus_data['features']        # list of dict
        record2 = next((item for item in df2 if item["id"] == nexus_id_str), False)
        if not(record2):
            print('SORRY, nexus_id not found:', nexus_id_str)
            print('Returning.')
            print()
            return -1

        nexus_type = record2['properties']['nexus_type']
        if (nexus_type != 'junction'):
            print( 'WARNING: nexus_type =', nexus_type )

        rgc = record2['geometry']['coordinates']
        coords = np.array( rgc , dtype='float64' )
        self.outlet_x = coords[0]
        self.outlet_y = coords[1]

        if (crs_name.endswith('CRS84')):
            self.outlet_lon = self.outlet_x
            self.outlet_lat = self.outlet_y
            crs_type = 'Geographic'
        elif (crs_name.endswith('5070')):
            crs_type = 'Albers Equal Area Conic'
        else:
            crs_type = 'Unknown'

        if (REPORT):
            print('Outlet coords for: ' + cat_id_str)
            print('  x, y =', self.outlet_x, ',', self.outlet_y )
            print('  crs_name =', crs_name )
            print('  crs_type =', crs_type)
            print()
            
        return coords
        
    #   get_outlet_coords()
    #--------------------------------------------------------------------
    def get_headwater_catchments(self, REPORT=False):

        #----------------------    
        # Get list of "toids"
        #----------------------
        # features -> properties -> toid
        df = self.catchment_data['features']    # list of dict
        toid_num_list = list() 
        for record in df:
            toid_str = record['properties']['toid']
            num_str = toid_str[4:]  # remove "nex-"
            toid_num_list.append( num_str )

        #-----------------------------------------------       
        # Create list of all "ids" not in "toids" list
        #-----------------------------------------------
        head_ids = list()
        for record in df:
            id_str  = record['id']
            num_str = id_str[4:]  # remove "cat-"
            if (num_str not in toid_num_list):
                head_id_str = 'cat-' + num_str
                head_ids.append( head_id_str )

        #---------------------------------------------
        # Option to sort these by their "area_sqkm",
        # and perhaps remove the smallest ones
        #---------------------------------------------
        n_heads = len(head_ids)
        areas = np.zeros( n_heads, dtype='float64' )
        k = 0
        for cat_id_str in head_ids:
            areas[ k ] = self.get_area_sqkm( cat_id_str )
            k += 1
        s = np.argsort( -1 * areas )   # -1 => large to small
#         amin = a.min()
#         amax = a.max()
        head_ids_arr = np.array( head_ids )
        head_ids_arr = head_ids_arr[ s ]
        head_ids     = head_ids_arr.tolist()
        
        if (REPORT):
            for cat_id_str in head_ids:
                print(cat_id_str, self.get_area_sqkm(cat_id_str))

        self.head_ids = head_ids
        return head_ids
        
    #   get_headwater_catchments()
    #--------------------------------------------------------------------    
    def get_area_sqkm(self, cat_id_str, REPORT=False):    

        crs_name = self.catchment_data['crs']['properties']['name']
        df = self.catchment_data['features']
        record = next((item for item in df if item["id"] == cat_id_str), False)
        if not(record):
            print('SORRY, cat_id not found:', cat_id_str)
            print('Returning.')
            print()
            return -1       

        self.area_sqkm = np.float64( record['properties']['area_sqkm'] )

        if (crs_name.endswith('CRS84')):
            crs_type = 'Geographic'
        elif (crs_name.endswith('5070')):
            crs_type = 'Albers Equal Area Conic'
        else:
            crs_type = 'Unknown'
            
        if (REPORT):
            print('Area for: ' + cat_id_str + ' [km2] =')
            print('  ', self.area_sqkm )
            print('  crs_name       =', crs_name )
            print('  crs_type       =', crs_type )
            print()
                       
        return self.area_sqkm
             
    #   get_area_sqkm()
    #--------------------------------------------------------------------
    def ocean_outlet_count(self, REPORT=False):
 
        df    = self.catchment_data['features']     
        recs  = [item for item in df if item['properties']['toid'] == 'nex-0']
        count = len(recs)
        self.ocean_outlet_count = count
        if (REPORT):
            print('Number of catchments draining to ocean =', count)
            print()  
        return count

    #   ocean_outlet_count()
    #--------------------------------------------------------------------     
    def print_info(self, cat_id_str='cat-29'):
   
        print('Extracting catchment info...') 
        area_sqkm = self.get_area_sqkm( cat_id_str )
        if (area_sqkm == -1):
            print('Catchment not found:', cat_id_str)
            return

        outlet_coords = self.get_outlet_coords( cat_id_str )
        bounds        = self.get_bounding_box( cat_id_str )
         
        print('Info for catchment: ' + cat_id_str)
        print('   catchment ID   = ', cat_id_str )
        print('   outlet lon,lat = ', outlet_coords )
        print('   area_sqkm =', area_sqkm, '[km2]' )
        print('   minlon    =', bounds[0] )
        print('   minlat    =', bounds[1] )
        print('   maxlon    =', bounds[2] )
        print('   maxlat    =', bounds[3] )
        print('   nexus_ID  =', self.nexus_id_str )
        print()
        
    #   print_info()
    #--------------------------------------------------------------------
    
#---------------------------------------------------------------------
def get_hydrofab_dem_info():

    #------------------------------------
    # bounds = [ulx, lry, lrx, uly]
    #        = [xmin, ymin, xmax, ymax]
    #------------------------------------
    dem_info = dict()
    dem_info['ncols'] = 163008
    dem_info['nrows'] = 114503
    dem_info['dx']    = 30.0
    dem_info['dy']    = 30.0
    dem_info['xmin']  = -2470965.0   # ulx
    dem_info['xmax']  =  2419275.0   # lrx
    dem_info['ymin']  =  186285.0    # lry     
    dem_info['ymax']  =  3621375.0   # uly 
    return dem_info
    
#   get_hydrofab_dem_info()
#------------------------------------------------------------------------
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

    ulx = -2470965  # -2.4709650000000005e+06
    uly = 3621375   # 3.6213750000000019e+06
    
    if (REPORT):
        print('lon, lat =', lon_deg, ',', lat_deg)
        print('x,   y   =', x, ',', y)
        print()
        print('x - ulx =', x - ulx)
        print('uly - y =', uly - y)

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
    d2r  = (np.pi / 180.)  # (degrees -> radians)
    lat  = lat_deg * d2r
    lon  = lon_deg * d2r
    
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

    dem_info = get_hydrofab_dem_info()
    ncols = dem_info['ncols']
    nrows = dem_info['nrows']
    xres  = dem_info['dx']
    yres  = dem_info['dy']
    ulx   = dem_info['xmin']
    uly   = dem_info['ymax']
    lry   = dem_info['ymin']

    col = (x - ulx) / xres
    row = (y - lry) / yres
    ## row = (y - uly) / yres
    # row = (nrows - 1) - (y - uly)/yres
    
    if (REPORT):
        print('lon, lat =', lon_deg, ',', lat_deg)
        print('x,   y   =', x, ',', y)
        print('ulx, uly =', ulx, ',', uly)
        print('lry =', lry)
        print('ncols, nrows =', ncols, nrows)
        print('col, row     =', col, row)
        print()
        print('x - ulx =', x - ulx)
        print('uly - y =', uly - y)

    return (x,y)
        
#   Albers_XY2()
#------------------------------------------------------------------------

      
       
       
       
       