## Copyright (c) 2001-2013, Scott D. Peckham

from topoflow.components import channels_diffusive_wave
from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
def test_instantiate():

    c = channels_diffusive_wave.channels_component()
    c.CCA   = False
    c.DEBUG = False
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory = tf_utils.TF_Test_Directory()
    cfg_prefix    = tf_utils.TF_Test_Case_Prefix()
    c.site_prefix = 'Treynor'  ###########      

    #-------------------------------------------
    # Call initialize() and call update() once
    #-------------------------------------------
##    print 'STATUS =', c.get_status()
##    c.initialize(mode="driver")
##    print 'STATUS =', c.get_status()
##    time_sec = float64(0)
##    c.update(time_sec)
##    print 'STATUS =', c.get_status()

    c.run_model( cfg_directory=cfg_directory,
                 cfg_prefix=cfg_prefix )
    
##    print 'KINEMATIC WAVE =', c.KINEMATIC_WAVE
##    print 'DIFFUSIVE WAVE =', c.DIFFUSIVE_WAVE
##    print 'DYNAMIC WAVE   =', c.DYNAMIC_WAVE
##    print ' '
##    print 'nx          =', c.nx
##    print 'ny          =', c.ny
    
    
#   test_instantiate()
#-----------------------------------------------------------------------

