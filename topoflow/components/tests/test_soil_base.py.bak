## Copyright (c) 2001-2013, Scott D. Peckham

from topoflow.components import soil_base
from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
def test_instantiate( soil_type='sandy_loam' ):

    soil = soil_base.soil_base( soil_type=soil_type )
    soil.CCA   = False
    soil.DEBUG = False
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory    = tf_utils.TF_Test_Directory()
    os.chdir( cfg_directory )
    cfg_prefix       = tf_utils.TF_Test_Case_Prefix()
    soil.site_prefix = 'Treynor'  #############
    
    soil.initialize( cfg_prefix=cfg_prefix, REPORT=True )

    print 'soil_type =', c.soil_type
    print 'Note:  t_r < t_H < t_w < t_f < t_i < t_u < t_s'
    print ' '
    print 'theta_r =', soil.theta_r   # ("absolute" min)
    print 'theta_H =', soil.theta_H   # (achievable min)
    print 'theta_w =', soil.theta_w
    print 'theta_f =', soil.theta_f
    print 'theta_i =', soil.theta_i
    print 'theta_u =', soil.theta_u   # (max possible)
    print 'theta_s =', soil.theta_s
    print ' '
    print 'K_s        =', soil.K_s
    print 'K(theta_s) =', soil.K_of_theta(soil.theta_s)
    print ' '
    print 'K_i        =', soil.K_i
    print 'K(theta_i) =', soil.K_of_theta(soil.theta_i)
    print ' '
    
#   test_instantiate()
#-----------------------------------------------------------------------

