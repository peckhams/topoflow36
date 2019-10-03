
# Copyright (c) 2001-2013, Scott D. Peckham
#
#-----------------------------------------------------------------------
#
#  unit_test()
#
#-----------------------------------------------------------------------
def unit_test( n_steps=10 ):

    c = erosion_component()
    c.CCA    = False
    c.DEBUG  = False
    c.SILENT = True    # (now the default in CSDMS_base.py)

    
##    c.DEBUG = True
##    c.SILENT = False

##    c.REPORT = True
    
    #-------------------------------------------------
    # NOTE: The Treynor_Iowa DEM has no depressions!
    #-------------------------------------------------  
##    directory   = tf_utils.TF_Test_Directory()
##    site_prefix = tf_utils.TF_Test_Site_Prefix()
##    case_prefix = tf_utils.TF_Test_Case_Prefix()

    #-----------------------------------
    # Use these when running on Mac.
    #-----------------------------------
##    cfg_directory   = '/Applications/Erode/Data/Test_50x50/'
##    cfg_prefix = 'Test'

##    cfg_directory   = '/Applications/Erode/Data/Test_100x100/'
##    cfg_prefix = 'Test'

##    cfg_directory   = '/Applications/Erode/Data/Test_200x200/'
##    cfg_prefix = 'Test'

##    cfg_directory   = '/Applications/Erode/Data/Test_400x400/'
##    cfg_prefix = 'Test'

    # These next few are obsolete now.
##    cfg_directory   = '/Applications/Erode/Data/Test2/'
##    ## cfg_directory   = '/data/sims/erode/3.1/Test2/'
##    cfg_prefix = 'Test2'

##    cfg_directory   = '/Applications/Erode/Data/Test4/'     # (80x80, 11/8/10)
##    cfg_prefix = 'Test4'

##    cfg_directory   = '/Applications/Erode/Data/Test2B/'  # (50x50, 11/8/10)
##    cfg_prefix = 'Test2B'

##    cfg_directory   = '/Applications/Erode/Data/Test3/'  # (20x20)
##    ## cfg_directory   = '/data/sims/erode/3.1/Test3/'
##    cfg_prefix = 'Test3'

    #------------------------------------------------
    # Use this when running on beach, but copy into
    # a separate sub directory when finished.
    #------------------------------------------------
    cfg_directory = '/home/beach/faculty/peckhams/CMT_Output/Erode_Global/'
    cfg_prefix    = 'Test'

##    cfg_directory   = '/home/beach/faculty/peckhams/ERODE3/Tests/Test_100x100/' # (9/1/11)
##    cfg_prefix = 'Test'

##    cfg_directory   = '/home/beach/faculty/peckhams/ERODE3/Tests/Test_200x200/' # (10/4/11)
##    cfg_prefix = 'Test'

##    cfg_directory   = '/home/beach/faculty/peckhams/ERODE3/Test_400x400/' # (10/4/11)
##    cfg_prefix = 'Test'

    #--------------------------------------------------------
    # Note: CSDMS_base.run_model() changes to cfg_directory
    #--------------------------------------------------------    
    c.run_model( cfg_directory=cfg_directory,
                 cfg_prefix=cfg_prefix,
                 n_steps=n_steps )
    
    #-------------------------------------------
    # Call initialize() and call update() once
    #-------------------------------------------
##    print 'STATUS =', c.get_status()
##    c.initialize(cfg_prefix=cfg_prefix, mode="driver")
##    print 'STATUS =', c.get_status()
##    time = float64(0)
##    c.update(time)
##    print 'STATUS =', c.get_status()

#   unit_test()
#-----------------------------------------------------------------------


