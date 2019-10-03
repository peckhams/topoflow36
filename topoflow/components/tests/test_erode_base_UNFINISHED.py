
# Copyright (c) 2001-2013, Scott D. Peckham
#
#-----------------------------------------------------------------------
#
#  unit_test()
#
#-----------------------------------------------------------------------
def unit_test( n_steps=10 ):

    c = erosion_component()
    c.CCA = False
    c.DEBUG = True
##    c.SILENT = False
##    c.REPORT = True
    
    #-------------------------------------------------
    # NOTE: The Treynor_Iowa DEM has no depressions!
    #-------------------------------------------------  
##    directory   = tf_utils.TF_Test_Directory()
##    site_prefix = tf_utils.TF_Test_Site_Prefix()
##    case_prefix = tf_utils.TF_Test_Case_Prefix()

    #-------------------
    # Use this on a Mac
    #-------------------
    # cfg_directory = '/Applications/Erode/Data/Test1/'

    #--------------------
    # Use this on beach
    #--------------------
    cfg_directory = '/home/csdms/models/erode/0.5/share/data/Test_100x100/'
    cfg_prefix = 'Test1'

    #--------------------------------------------------------
    # Note: CSDMS_base.run_model() changes to cfg_directory
    #--------------------------------------------------------
    c.run_model( cfg_directory=cfg_directory,
                 cfg_prefix=cfg_prefix,
                 n_steps=n_steps )

#   unit_test()
#-----------------------------------------------------------------------


