
## Copyright (c) 2013, Scott D. Peckham
## August 2013
## Unit tests for "template_file.py" in "utils" folder.

import os
from topoflow.utils import template_files

#-------------------------------------------------------------------------
#
# get_replacement_dictionary()
# make_test_template()
# test1()
#
#-------------------------------------------------------------------------
def get_replacement_dictionary( var_names=None, values=None,
                                method='RE'):

    #-----------------------
    # Defaults for testing
    #-----------------------
    if (var_names == None):
        if (method == 'RE'):
            var_names=['${var1}','${var2}']
        else:
            var_names=['var1','var2']
    if (values == None):
        values = [3.141, 2.718]
        
    #-------------------------------------
    # Store replacements in a dictionary
    #-------------------------------------
    d = dict( list(zip( var_names, values )) )
    return d

#   get_replacement_dictionary()
#-------------------------------------------------------------------------
def make_test_template( template_file=None ):

    if (template_file == None):
        template_file='Model_config_file_template.cfg.in'
        
    template_unit = open( template_file, 'w' )
    lines = ['Configuration file template for a fake model\n',
             '\n',
             'Test case of standard key-value pairs:\n',
             'var1 = ${var1}\n',
             'var2 = ${var2}\n',
             'Test case of no replacement in dictionary:\n',
             'var3 = ${var3}\n',
             'Test case of two replacements on one line:\n',
             'var1, var2 = ${var1}, ${var2}\n',
             'More lines\n']
    template_unit.writelines( lines )
    template_unit.close()

#   make_test_template()
#-------------------------------------------------------------------------
def test1( method='RE' ):

    #-----------------------------
    # Change to a test directory
    #-----------------------------
    test_dir = '/Users/peckhams/Documents'
    os.chdir( test_dir )

    #---------------------------------------------------
    # Get dictionary of replacements (values for vars)
    #---------------------------------------------------
    dictionary = get_replacement_dictionary( method=method )

    #-----------------------------------------
    # Make a fake model config file template
    #-----------------------------------------
    cfg_template_file = 'Model_config_file_template.cfg.in'
    make_test_template( cfg_template_file )

    #-----------------------------------------------------
    # Make replacements to create fake model config file
    #-----------------------------------------------------
    new_cfg_file = 'Model_config_file.cfg'
    template_files.replace( cfg_template_file, new_cfg_file,
                            dictionary, method=method )
    
#   test1()
#-------------------------------------------------------------------------

