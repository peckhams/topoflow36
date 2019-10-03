## Copyright (c) 2001-2013, Scott D. Peckham

from topoflow.components import topoflow_driver
from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
def test_instantiate():

    c = topoflow_driver.topoflow_driver()
    c.CCA   = False
    c.DEBUG = False
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory   = tf_utils.TF_Test_Directory()
    cfg_prefix      = tf_utils.TF_Test_Case_Prefix()
    c.site_prefix = 'Treynor'  #############    

    #-------------------------------------------
    # Call initialize() and call update() once
    #-------------------------------------------
##    print 'STATUS =', c.get_status()
##    c.initialize( mode="driver" )
##    print 'STATUS =', c.get_status()
##    time_sec = np.float64(0)
##    c.update(time_sec)
##    print 'STATUS =', c.get_status()

    c.run_model( cfg_directory=cfg_directory,
                 cfg_prefix=cfg_prefix )

#   test_instantiate()
#-----------------------------------------------------------------------



