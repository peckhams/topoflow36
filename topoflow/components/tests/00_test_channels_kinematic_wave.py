## Copyright (c) 2001-2013, Scott D. Peckham

from topoflow.components import channels_kinematic_wave
from topoflow.framework import framework
# from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
#
# test1()   # Instantiate and call run_model() method.
#
#-----------------------------------------------------------------------
def test1():

    c = channels_kinematic_wave.channels_component()
    c.CCA   = False
    c.DEBUG = False

    #----------------------------------------------
    # Get path to the current file (framework.py)
    # At top need: "#! /usr/bin/env python" ??
    #----------------------------------------------
    paths = framework.get_package_paths()
    framework_dir = paths['framework']
    examples_dir  = paths['examples']
    #-------------------------------------------
    # Set cfg_prefix and cfg_directory for test
    #-------------------------------------------
    cfg_prefix    = 'June_20_67'
    cfg_directory = examples_dir + 'Treynor_Iowa/'
    # c.site_prefix = 'Treynor'  #### Still need this ?
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
##    cfg_directory = tf_utils.TF_Test_Directory()
##    cfg_prefix    = tf_utils.TF_Test_Case_Prefix()
##    c.site_prefix = 'Treynor'  ###########      

    ## cfg_directory = '/data/progs/topoflow/3.0/data/test_plane_canal/'
    #cfg_directory = '/Applications/TopoFlow/Data/Test_Plane_Canal/'
    #cfg_prefix    = 'Case5'
    #c.site_prefix = 'plane' ########
    
    #-------------------------------------------
    # Call initialize() and call update() once
    #-------------------------------------------
##    print 'STATUS =', c.get_status()
##    c.initialize(mode="driver")
##    print 'STATUS =', c.get_status()
##    time_sec = float64(0)
##    c.update(time_sec)
##    print 'STATUS =', c.get_status()

    #-------------------------------------------------
    # Can't run this without the framework, because
    # it provides refs to vars from other components.
    # This fails when self.P is needed (from met).
    #-------------------------------------------------
    c.run_model( cfg_directory=cfg_directory,
                 cfg_prefix=cfg_prefix )
    
##    print 'KINEMATIC WAVE =', c.KINEMATIC_WAVE
##    print 'DIFFUSIVE WAVE =', c.DIFFUSIVE_WAVE
##    print 'DYNAMIC WAVE   =', c.DYNAMIC_WAVE
##    print ' '
##    print 'nx          =', c.nx
##    print 'ny          =', c.ny
    
#   test1()
#-----------------------------------------------------------------------

