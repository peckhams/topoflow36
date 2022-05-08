# Copyright (c) 2022, Scott D. Peckham
# April 2022
# Tools for working with CFG "template files"
# See old "template_files.py" in "utils" folder.
# There is a related Python package called jinja.
#
#-------------------------------------------------------------------------
#
# read_cfg_as_string()
# write_string_as_cfg()
# make_new_cfg_file()
#
#-------------------------------------------------------------------------
def read_cfg_as_string( cfg_template ):

    #------------------------------------------- 
    # Note: CFG files are usually quite small.
    #-------------------------------------------
    cfg_unit = open(cfg_template, 'r')
    string = cfg_unit.read()
    cfg_unit.close()
    return string

#   read_cfg_as_string()
#-------------------------------------------------------------------------
def write_string_as_cfg( string, cfg_file ):

    cfg_unit = open(cfg_file, 'w')
    result = cfg_unit.write( string )
    cfg_unit.close()
  
#   write_string_as_cfg()
#-------------------------------------------------------------------------
def make_new_cfg_file( vars, vals, cfg_template, cfg_file ):

    text = read_cfg_as_string( cfg_template )
    k = 0
    for var in vars:
        old_str = '{{' + var + '}}'
        new_str = str( vals[k] )
        text = text.replace( old_str, new_str )
        k += 1
    write_string_as_cfg( text, cfg_file )
 
#   make_new_cfg_file()
#-------------------------------------------------------------------------

