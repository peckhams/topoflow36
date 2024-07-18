#! /usr/bin/env python
#
# Copyright (c) 2022, Scott D. Peckham
#
# class tf36_bmi()
#
# run_model()
#
#-----------------------------------------------------------------------

from topoflow.framework import multi_bmi


#-----------------------------------------------------------------------
class tf36_bmi( multi_bmi.multi_bmi ):
    # Note: This just renames the class.
    pass

#-----------------------------------------------------------------------
def run_model( cfg_file=None, cat_id_str='cat-209',
               ngen_dir=None, SILENT=False ):

    if (ngen_dir is None):
        print('========================================================')
        print('ERROR: You must specify the full path to the directory')
        print('where you have installed NextGen, using the ngen_dir')
        print('keyword:  e.g., ngen_dir="~/ngen".')
        print('========================================================')
        print()
        return
    tf = tf36_bmi( SILENT=SILENT )
    tf.ngen_dir    = ngen_dir
    tf.site_prefix = cat_id_str
    tf.case_prefix = 'Test1'
    
    if (cfg_file is None):
        if (cat_id_str is not None):
            cfg_file = tf.get_test_cfg_file( test=cat_id_str)
        else:
            cfg_file = tf.get_test_cfg_file( test='Treynor_June_07_67' )
        # cfg_file = tf.get_test_cfg_file( test='Treynor_June_20_67' )
        
    tf.run_model( cfg_file=cfg_file )

#   run_model()
#-------------------------------------------------------------------
