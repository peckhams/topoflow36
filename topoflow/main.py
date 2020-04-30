#! /usr/bin/env python
#
# Copyright (c) 2001-2020, Scott D. Peckham
#
# run_model()
# run_test()
#
#-----------------------------------------------------------------------

from topoflow.framework import emeli

#-----------------------------------------------------------------------
def run_model( driver_comp_name ='topoflow_driver',
               cfg_prefix='June_20_67',
               cfg_directory=None,
               time_interp_method='Linear'):

    #----------------------------------------------------------
    # Note: The "driver_comp_name " defaults to using a
    #       component of "topoflow_driver" type as the driver.
    #       The component of this type specified in the
    #       provider_file will be the driver component.
    #
    #       Any other component in the provider_file can
    #       also be used as the driver.  Examples are:
    #       meteorology, channels, snow, satzone, evap,
    #       infil, diversions, ice
    #-----------------------------------------------------
    f = emeli.framework()
    examples_dir = emeli.paths['examples']
    ## examples_dir = f.paths['examples']
    
    if (cfg_directory == None):
        cfg_directory = examples_dir + 'Treynor_Iowa_30m/'

    f.run_model( driver_comp_name =driver_comp_name ,
                 cfg_prefix=cfg_prefix,
                 cfg_directory=cfg_directory,
                 time_interp_method=time_interp_method )

#   run_model()
#-------------------------------------------------------------------
def run_test(test='treynor_kin_wave'):

    #---------------------------------------------------
    # Notes: These all use the Treynor, Iowa watershed
    #        and historic rain storm of June 20, 1967.
    #        Treynor is within the Nishnabotna River.
    #---------------------------------------------------
    cfg_prefix = 'June_20_67'
    test_dir = emeli.paths['examples']
    test_dir += 'Treynor_Iowa_30m/__Tests/'
 
    if (test == 'treynor_kin_wave'):
        cfg_dir = test_dir + 'Test1_kin_wave/'
    elif (test == 'treynor_diff_wave'):
        cfg_prefix = 'June_20_67'
        cfg_dir  = test_dir + 'Test2_diff_wave'
    elif (test == 'treynor_infil_smith_parlange'):
        cfg_prefix = 'June_20_67'
        cfg_dir  = test_dir + 'Test3_infil_smith_parlange'
        
    #----------------------       
    # Run the chosen test
    #----------------------  
    run_model(cfg_prefix=cfg_prefix, cfg_directory=cfg_dir)
   
#   run_test()
#-------------------------------------------------------------------
