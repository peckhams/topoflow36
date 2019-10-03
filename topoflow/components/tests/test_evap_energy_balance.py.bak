## Copyright (c) 2001-2013, Scott D. Peckham

from topoflow.components import evap_energy_balance
from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
def test_instantiate():

    g = evap_energy_balance.evap_component()
    g.CCA   = False
    g.DEBUG = False
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory = tf_utils.TF_Test_Directory()
    cfg_prefix    = tf_utils.TF_Test_Case_Prefix()
    g.site_prefix = 'Treynor'  #############  
  
    g.run_model( cfg_directory=cfg_directory,
                 cfg_prefix=cfg_prefix )
    
##    g.initialize( cfg_prefix=cfg_prefix, mode="driver" )
    
    print 'nx =', g.nx
    print 'ny =', g.ny
    
#   test_instantiate()
#-----------------------------------------------------------------------

