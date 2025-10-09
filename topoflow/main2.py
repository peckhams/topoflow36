#! /usr/bin/env python
#
# Copyright (c) 2022-2024, Scott D. Peckham
#
# class tf36_bmi()
#
# run_model()
#
#-----------------------------------------------------------------------
# Example usage 1:
#
# conda activate tf36
# % python
# >>> from topoflow import main2
# >>> cat_id_str  = 'cat-209'
# >>> start_dir   = '/Users/peckhams/Dropbox/GitHub/'
# >>> ngen_dir    = start_dir + 'ngen/'
# >>> forcing_dir = ngen_dir  + 'data/topoflow/forcing/huc01/'
# >>> spatial_dir = ngen_dir  + 'data/topoflow/spatial/huc01/'
# >>> catchment_file = 'catchment_data_HUC01.geojson'
# >>> nexus_file = 'nexus_data_HUC01.geojson'
# >>> main2.run_model( cat_id_str=cat_id_str, ngen_dir=ngen_dir,
#           spatial_dir=spatial_dir, forcing_dir=forcing_dir,
#           catchment_file=catchment_file,  nexus_file=nexus_file )
#
#-----------------------------------------------------------------------
# Example usage 2:
#
# conda activate tf36
# % python
# >>> from topoflow import main2
# >>> cat_id_str  = 'cat-11223',
# >>> start_dir   = '/Users/peckhams/Dropbox/GitHub/'
# >>> ngen_dir    = start_dir + 'ngen/'
# >>> forcing_dir = ngen_dir  + 'data/topoflow/forcing/lauren_test/'
# >>> spatial_dir = ngen_dir  + 'data/topoflow/spatial/lauren_test/'
# >>> catchment_file = 'gauge_01073000_catchment_data.geojson'
# >>> nexus_file = 'gauge_01073000_nexus_data.geojson'
# >>> main2.run_model( cat_id_str=cat_id_str, ngen_dir=ngen_dir,
#           spatial_dir=spatial_dir, forcing_dir=forcing_dir,
#           catchment_file=catchment_file, nexus_file=nexus_file )
#
#-----------------------------------------------------------------------
from topoflow.framework import multi_bmi

#-----------------------------------------------------------------------
class tf36_bmi( multi_bmi.multi_bmi ):
    # Note: This just renames the class.
    pass

#-----------------------------------------------------------------------
def run_model( cfg_file=None, cat_id_str='cat-209',
               ngen_dir=None, spatial_dir=None, forcing_dir=None,
               catchment_file=None, nexus_file=None, SILENT=False ):

    #--------------------------------------------------------------
    # Note: The purpose of this function is to run TopoFlow using
    #       the multi-BMI mechanism (one BMI exposed to NextGen)
    #       while using the same directories and files as will be
    #       used when TopoFlow is run from within NextGen.
    #       This was first developed and used for Fall AGU 2022.
    #--------------------------------------------------------------
    if (ngen_dir is None):
        print('========================================================')
        print('ERROR: You must specify the full path to the directory')
        print('where you have installed NextGen, using the ngen_dir')
        print('keyword:  e.g., ngen_dir="~/ngen".')
        print('========================================================')
        print()
        return
    if (spatial_dir is None):
        print('========================================================')
        print('ERROR: You must specify the full path to the directory')
        print('that contains the hydrofabric GeoJSON files for NextGen')
        print('using the spatial_dir keyword.  For example:')
        print('   ~/ngen/data/topoflow/spatial/huc01/')
        print('========================================================')
        print()
        return
    if (forcing_dir is None):
        print('========================================================')
        print('ERROR: You must specify the full path to the directory')
        print('that contains the forcing (or met) data for NextGen')
        print('using the forcing_dir keyword.  For example:')
        print('   ~/ngen/data/topoflow/forcing/huc01/')
        print('========================================================')
        print()
        return
    if (catchment_file is None):
        print('===========================================================')
        print('ERROR: You must specify the filename of the GeoJSON')
        print('or GPKG file that contains the hydrofabric catchment data.')
        print('This file must be in the spatial_dir, e.g.')
        print('   ~/ngen/data/topoflow/spatial/')
        print('For example: catchment_data_HUC01.geojson') 
        print('===========================================================')
        print()
        return
    if (nexus_file is None):
        print('===========================================================')
        print('ERROR: You must specify the filename of the GeoJSON')
        print('or GPKG file that contains the hydrofabric nexus data.')
        print('This file must be in the spatial_dir:')
        print('   e.g. ~/ngen/data/topoflow/spatial/')
        print('For example: nexus_data_HUC01.geojson') 
        print('===========================================================')
        print()
        return
    #------------------------------------
    # If missing, add a slash character
    #------------------------------------
    if (ngen_dir[-1] != '/'):
        ngen_dir += '/'
    if (spatial_dir[-1] != '/'):
        spatial_dir += '/'
    if (forcing_dir[-1] != '/'):
        forcing_dir += '/' 

    #-------------------------------------------------        
    # Create instance of tf36_bmi and set parameters       
    #-------------------------------------------------
    tf = tf36_bmi( SILENT=SILENT )
    tf.ngen_dir       = ngen_dir
    tf.forcing_dir    = forcing_dir
    tf.spatial_dir    = spatial_dir
    tf.catchment_file = catchment_file
    tf.nexus_file     = nexus_file   
    #-----------------------------------------------------
    tf.site_prefix = cat_id_str
    tf.case_prefix = 'Test1'
    
    if (cfg_file is None):
        #----------------------------
        # Create a default CFG file
        #----------------------------
        if (cat_id_str is not None):
            cfg_file = tf.get_test_cfg_file( test=cat_id_str)
        else:
            cfg_file = tf.get_test_cfg_file( test='Treynor_June_07_67' )
        # cfg_file = tf.get_test_cfg_file( test='Treynor_June_20_67' )
    #------------------------------------------------------------   
    # Note: tf.run_model() calls initialize(), and initialize()
    #       calls prepare_tf_inputs(). This creates an instance
    #       of get_inputs() class, defined in prepare_inputs.py.
    #       Input params from here are passed there.
    #------------------------------------------------------------     
    tf.run_model( cfg_file=cfg_file )

#   run_model()
#-------------------------------------------------------------------
