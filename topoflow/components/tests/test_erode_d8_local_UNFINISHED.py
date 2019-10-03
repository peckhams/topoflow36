#
# Copyright (c) 2001-2013, Scott D. Peckham
#
#-------------------------------------------------------------
#
#  get_test_data_info()             ## (3/1/12)
#  unit_test()
#
#-------------------------------------------------------------
def get_test_data_info( data_name='100' ):

    #------------------------------------------------------
    # Get name of platform (Darwin, Linux, Windows, Java)
    #------------------------------------------------------
    os_name = platform.system()
    beach_test_dir = '~/ERODE3/Tests/'
    cfg_prefix = 'Test'
    
    #------------------------------------------------
    # (100 x 100) test case
    #-----------------------------------------------
    if (data_name.upper() == '100'):
        # cfg_prefix    = 'New100'
        if (os_name == 'Darwin'):
            cfg_directory = '/Applications/Erode/Data/Test_100x100/'
        else:
            cfg_directory = beach_test_dir + 'Test_100x100/'

    #------------------------------------------------
    # (50 x 50) test case
    #-----------------------------------------------
    if (data_name.upper() == '50'):
        # cfg_prefix    = 'New50'
        if (os_name == 'Darwin'):
            cfg_directory = '/Applications/Erode/Data/Test_50x50/'
        else:
            cfg_directory = beach_test_dir + 'Test_50x50/'
            
    #------------------------------------------------
    # (200 x 200) test case
    #-----------------------------------------------
    if (data_name.upper() == '200'):
        # cfg_prefix    = 'New200'
        if (os_name == 'Darwin'):
            cfg_directory = '/Applications/Erode/Data/Test_200x200/'
        else:
            cfg_directory = beach_test_dir + 'Test_200x200/'

    #------------------------------------------------
    # (400 x 400) test case
    #-----------------------------------------------
    if (data_name.upper() == '400'):
        # cfg_prefix    = 'New400'
        if (os_name == 'Darwin'):
            cfg_directory = '/Applications/Erode/Data/Test_400x400/'
        else:
            cfg_directory = beach_test_dir + 'Test_400x400/'

    #------------------------------------------------
    # (800 x 800) test case
    #-----------------------------------------------
    if (data_name.upper() == '800'):
        # cfg_prefix    = 'New800'
        if (os_name == 'Darwin'):
            cfg_directory = '/Applications/Erode/Data/Test_800x800/'
        else:
            cfg_directory = beach_test_dir + 'Test_800x800/'
            
    #-----------------------------------
    return cfg_directory, cfg_prefix

#   get_test_data_info()
#-----------------------------------------------------------------------
def unit_test( data_name='100', n_steps=10 ):

    #------------------------------------------------------------------
    # Notes: Make sure that "Area grid units" in <prefix>_d8.cfg
    #        file are "m^2" vs. "km^2", except in the case of
    #        using the tests in d8_local.py.  Otherwise, the
    #        model will crash at the sys.exit() in update_dt_grid5()
    #        after only a few steps.
    #------------------------------------------------------------------
    c = erosion_component()
    c.CCA = False
    c.REPORT = False

    c.DEBUG  = False
    c.SILENT = True     # (now the default in CSDMS_base.py)

##    c.DEBUG  = True
##    c.SILENT = False

    #----------------------------------------------
    # Get cfg directory and prefix for unit tests
    #----------------------------------------------
    cfg_directory, cfg_prefix = get_test_data_info( data_name )
    
    #--------------------------
    # Change to CFG directory
    #--------------------------
    # os.chdir( cfg_directory )

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

