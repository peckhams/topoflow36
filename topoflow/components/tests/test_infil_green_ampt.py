## Copyright (c) 2001-2013, Scott D. Peckham

from topoflow.components import infil_green_ampt
from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
def test_instantiate():

    c = infil_green_ampt.infil_component
    c.CCA   = False
    c.DEBUG = False
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory = tf_utils.TF_Test_Directory()
    cfg_prefix    = tf_utils.TF_Test_Case_Prefix()
    c.site_prefix = 'Treynor'  #############  
  
    c.run_model( cfg_directory=cfg_directory,
                 cfg_prefix=cfg_prefix )
    
##    c.initialize( cfg_prefix=cfg_prefix, mode="driver" )
    
#   test_instantiate()
#-----------------------------------------------------------------------

