## Copyright (c) 2001-2013, Scott D. Peckham

from topoflow.components import evap_priestley_taylor
from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
def test_instantiate():

    ec = evap_priestley_taylor.evap_component()
    ec.CCA   = False
    ec.DEBUG = False
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory   = tf_utils.TF_Test_Directory()
    os.chdir( cfg_directory )
    cfg_prefix      = tf_utils.TF_Test_Case_Prefix()
    ec.site_prefix  = 'Treynor'  #############  
  
    ec.run_model( cfg_directory=cfg_directory,
                  cfg_prefix=cfg_prefix )

##    ec.initialize( cfg_prefix=cfg_prefix, mode="driver" )
 
##    print 'nx =', ec.nx
##    print 'ny =', ec.ny
    
#   test_instantiate()
#-----------------------------------------------------------------------

