## Copyright (c) 2001-2013, Scott D. Peckham

from topoflow.components import ice_base
from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
def test_instantiate():

    c = ice_base.ice_component
    c.CCA   = False
    c.DEBUG = False
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory = tf_utils.TF_Test_Directory()
    cfg_prefix    = tf_utils.TF_Test_Case_Prefix()
    c.site_prefix = 'Treynor'  #############  

    #-------------------------------
    # Initialize and 1 update call
    #-------------------------------
##    print 'STATUS =', c.get_status()
##    c.initialize( mode="driver" )
##    print 'STATUS =', c.get_status()
##    time_sec = np.float64(0)
##    c.update(time_sec)
##    print 'STATUS =', ic.get_status()
    
    c.run_model( cfg_directory=cfg_directory,
                 cfg_prefix=cfg_prefix )
    
##    c.initialize( cfg_prefix=cfg_prefix, mode="driver" )

    #-----------------------------------------
    # Save Zi for Animas as flat binary grid
    #-----------------------------------------   
##    c.save_matlab_dem_as_rtg() # (mat_file = 'Animas_200.mat')
    
#   test_instantiate()
#-----------------------------------------------------------------------

