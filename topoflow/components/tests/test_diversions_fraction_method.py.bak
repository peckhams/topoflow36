## Copyright (c) 2001-2013, Scott D. Peckham

from topoflow.components import diversions_fraction_method
from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
def test_instantiate():

    d = diversions_fraction_method.diversions_component()
    d.CCA   = False
    d.DEBUG = False
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory   = tf_utils.TF_Test_Directory()
    os.chdir( cfg_directory )
    cfg_prefix      = tf_utils.TF_Test_Case_Prefix()
    d.site_prefix   = 'Treynor'  #############  
    
    ## d.run_model( cfg_directory=cfg_directory,
    ##              cfg_prefix=cfg_prefix )
                 
    d.initialize( cfg_prefix=cfg_prefix, mode="driver" )

    print 'use_sources =', d.use_sources
    print 'use_sinks   =', d.use_sinks
    print 'use_canals  =', d.use_canals
    print ' '
    print 'source_file =', d.source_file
    print 'sink_file   =', d.sink_file
    print 'canal_file  =', d.canal_file
    
    print 'Finished with unit_test().'
    
#   test_instantiate()
#-----------------------------------------------------------------------

