
# Copyright (c) 2023, Scott D. Peckham
#
# Jun 2023. Started from gages2_tools.py.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils import mopex_tools as mt
#  >>> mt.creat_tsv( nf_max=500, SWAP_XY=False )
#
#---------------------------------------------------------------------
#
#  get_basin_repo_dir()
#  get_mopex_data_dir()
#  create_tsv()
#
#---------------------------------------------------------------------

import numpy as np
from osgeo import ogr, osr
import json, sys, time

from topoflow.utils import shape_utils as su

#---------------------------------------------------------------------
def get_basin_repo_dir():

    #-----------------------------------
    # Modify this directory as needed.
    #-----------------------------------
    repo_dir  = '/Users/peckhams/Dropbox/NOAA_NextGen/'
    repo_dir += '__NextGen_Example_Basin_Repo/'
    return repo_dir

#   get_basin_repo_dir()
#---------------------------------------------------------------------
def get_mopex_data_dir():

    #-----------------------------------
    # Modify this directory as needed.
    #-----------------------------------
    repo_dir  = get_basin_repo_dir()
    data_dir  = repo_dir + 'MOPEX/'
    data_dir += 'Hydrologic Synthesis Project 2009_18GB/'
    data_dir += 'spatialdata/MOPEX431_basins/'
    return data_dir
       
#   get_mopex_data_dir()
#---------------------------------------------------------------------
def create_tsv( data_dir=None, nf_max=50,
                shape_file='MOPEX431_basins.shp',                    
                prj_file  ='MOPEX431_basins.prj',
                tsv_file='new_mopex431.tsv',
                SWAP_XY=False,  # False works for MOPEX.
                REPORT=True):

    if (data_dir is None):
        data_dir = get_mopex_data_dir()

    insert_key = 'Longitude'

    su.read_shapefile(data_dir=data_dir, shape_file=shape_file,
                      prj_file=prj_file, tsv_file=tsv_file,
                      nf_max=nf_max,
                      SWAP_XY=SWAP_XY, REPORT=REPORT,
                      ADD_BOUNDING_BOX=True, insert_key=insert_key,
                      filter_key=None, filter_value=None)

#   create_tsv()
#---------------------------------------------------------------------


