#! /usr/bin/env python
#
# Copyright (c) 2022-2024, Scott D. Peckham
#
# class tf36_bmi()
#
# run_model()
#
#-----------------------------------------------------------------------
# Example usage 1:  Create nextgen CFG file with all info.
#
# conda activate tf36
# % python
# >>> from topoflow import main_ngen as mn
# >>> ngen_dir   = '/Users/peckhams/Dropbox/GitHub/ngen/'
# >>> basin_dir  = ngen_dir + 'data/topoflow/input_files/cat-11223/'
# >>> cfg_file   = basin_dir + 'Test1_cfg/Test1_nextgen.cfg'
# >>> mn.run_model( cfg_file=cfg_file )
#
#-----------------------------------------------------------------------
from topoflow.framework import multi_bmi

#-----------------------------------------------------------------------
class tf36_bmi( multi_bmi.multi_bmi ):
    # Note: This just renames the class.
    pass

#-----------------------------------------------------------------------
def run_model( cfg_file=None, SILENT=False ):

    #--------------------------------------------------------------
    # Note: The purpose of this function is to run TopoFlow using
    #       the multi-BMI mechanism (one BMI exposed to NextGen)
    #       while using the same directories and files as will be
    #       used when TopoFlow is run from within NextGen.
    #       If not already available, it automatically downloads
    #       a DEM for the specified NextGen catchment and then
    #       creates all required input files, including a set of
    #       default CFG files.
    #       This was first developed and used for Fall AGU 2022.
    #--------------------------------------------------------------
    if (cfg_file is None):
        print('SORRY, you must specify the full path to a')
        print('TopoFlow/NextGen CFG file.  If "ngen_dir" is')
        print('the NextGen directory, then the "basin_dir" will')
        print('typically have the form:')
        print('   basin_dir = ngen_dir + "data/topoflow/input_files/cat-84/"')
        print('and the cfg_file will then have the form:')
        print('   cfg_file = basin_dir + "Test1_cfg/Test1_nextgen.cfg".')
        print() 
        return

    #-------------------------------------------        
    # Create an instance of the tf36_bmi class       
    #-------------------------------------------
    tf = tf36_bmi( SILENT=SILENT )
    
#     if (cfg_file is None):
#         #----------------------------
#         # Create a default CFG file
#         #----------------------------
#         if (cat_id_str is not None):
#             cfg_file = tf.get_test_cfg_file( test=cat_id_str)
#         else:
#             cfg_file = tf.get_test_cfg_file( test='Treynor_June_07_67' )
#         # cfg_file = tf.get_test_cfg_file( test='Treynor_June_20_67' )

    #------------------------------------------------------------   
    # Note: tf.run_model() calls initialize(), and initialize()
    #       calls prepare_inputs_for_ngen_basin(). This creates
    #       an instance of get_inputs() class, defined in
    #       prepare_inputs.py.
    #------------------------------------------------------------     
    tf.run_model( cfg_file=cfg_file )

#   run_model()
#-------------------------------------------------------------------
