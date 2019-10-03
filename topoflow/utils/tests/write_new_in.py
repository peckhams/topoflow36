#!/usr/bin/env python

# Copyright (c) 2011, Scott D. Peckham

#------------------------------------------------------
# S.D. Peckham
# May 9, 2011
#
# Tool to convert my old ".cfg" files to a new format
# which no longer has GUI info.
#
# Example of use at a Unix prompt:
#
#    % ./write_new_in.py my_file.cfg
#------------------------------------------------------
#
# Functions:
#
# write_new_in()  # (write a ".in" template file)
#
#------------------------------------------------------

import os.path
import sys

#------------------------------------------------------
def write_new_in( cfg_file, word_delim='|' ):
 
    #--------------------
    # Open the CFG file
    #--------------------
    try:
        cfg_unit = open( cfg_file, 'r' )
    except:
        print('SORRY: Could not open CFG file named:')
        print('       ' + cfg_file)

    #-----------------------------
    # Open new CFG template file
    #-----------------------------
    pos    = cfg_file.rfind('.')
    prefix = cfg_file[0:pos]         # (still used below)
    new_cfg_file = cfg_file + '.in'
    ## new_cfg_file = prefix + '_in.cfg'
    #----------------------------------------------
    NEW_CFG_EXISTS = os.path.exists( new_cfg_file )
    if (NEW_CFG_EXISTS):
        print('SORRY, A "new" CFG file with the name')
        print('       ' + new_cfg_file)
        print('       already exists.')
        return
    new_cfg_unit = open( new_cfg_file, 'w' )

    #-----------------
    # Write a header
    #-----------------
    new_cfg_unit.write('#'.ljust(80, '=') + '\n')
    header = '# TopoFlow Config File for: ' + prefix
    new_cfg_unit.write( header + '\n' )
    
    #-----------------------------
    # Read CFG file header lines
    #-----------------------------
    for k in range(4):
        line = cfg_unit.readline()

    ## FIRST_TAB = True
    DONE = False
    while (True):
        #-------------------------------
        # Read data line from CFG file
        #-------------------------------   
        line = cfg_unit.readline()
        if (line == ''):
            break
        words = line.split( word_delim )
        first_word = words[0].strip()
        
        if (first_word == 'PPF_GROUP_NAME'):
            #-------------------
            # Create a new tab
            #-------------------
            new_cfg_unit.write('#'.ljust(80, '=') + '\n')
            new_cfg_unit.write('# ' + words[1].strip() + '\n')
            ## FIRST_TAB = False
        elif (first_word == 'HTML_HELP_FILE'):
            pass
            ## help_url = words[3].strip()
        elif (first_word == ''):
            pass
        elif (first_word[0] == '#'):
            pass
        else:              
            #---------------------------
            # Prepare values to write
            #-------------------------------------------------
            # var_type will be string, long, int or float;
            # NOT scalar, time_series, grid or grid_sequence.
            # Latter types are stored in separate vars now
            # (before the var itself), with type "string".
            # CSDMS_base.read_config_file() will handle
            # things accordingly.
            #-------------------------------------------------
            var_name = first_word.ljust(20) + '| '  # (pad to length 20 with spaces)
            var_type = words[2].strip().ljust(10) + '| '
            value    = ('${' + first_word + '}').ljust(20) + '| '
            help_str = words[6].strip()

            #-------------------------------
            # Write a line in new CFG file
            #-------------------------------
            new_cfg_unit.write(var_name + value + var_type + help_str + '\n')

    #---------------------
    # Close the CFG file
    #---------------------
    cfg_unit.close()
          
    #-------------------------------
    # Close the new-style CFG file
    #-------------------------------
    new_cfg_unit.close()
    print('Finished writing new TopoFlow CFG template file.')
    print(' ')
    
#   write_new_in()
#-----------------------------------------------------------------------
if (__name__ == "__main__"):
    
    #-----------------------------------------------------
    # Note: First arg in sys.argv is the command itself.
    #-----------------------------------------------------
    n_args = len(sys.argv)
    if (n_args < 2):
        print('ERROR: This tool requires a CFG filename')
        print('       argument.')
        print('sys.argv =', sys.argv)
        print(' ')
    elif (n_args == 2):
        write_new_in( sys.argv[1] )
    else:
        print('ERROR: Invalid number of arguments.')
        
#-----------------------------------------------------------------------

        
