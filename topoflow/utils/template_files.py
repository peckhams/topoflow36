
## Copyright (c) 2013, Scott D. Peckham
## August 2013
## Tools for working with CFG "template files"
## See "test_template_file.py" in "utils/tests" folder.

import glob
import re      # (for re.compile, re.findall)
import string  # (for string.Template)

#-------------------------------------------------------------------------
#
#   get_replacements()
#   replace()
#
#-------------------------------------------------------------------------
def get_replacements( var_names=None, values=None, method='RE'):

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
        
    #---------------------------------------------
    # Read var_names & values from a file here ?
    #---------------------------------------------

    #-------------------------------------
    # Store replacements in a dictionary
    #-------------------------------------
    d = dict( list(zip( var_names, values )) )
    return d

#   get_replacements()
#-------------------------------------------------------------------------
def replace( cfg_template_file, new_cfg_file, dictionary=None,
             method='RE'):

    #---------------------------------------------------------
    # Note: Rename this to "make_cfg_file" ??
    #---------------------------------------------------------
    # Note: If (method == 'RE'), use the re package.
    #       If (method == 'STRING'), use the string package.
    #---------------------------------------------------------
    # If (method == 'STRING'), var_names used to build the
    # mapping dictionary should not have the ${.} wrapper.
    #---------------------------------------------------------
    RE_METHOD     = (method == 'RE')
    STRING_METHOD = not(RE_METHOD)
    
    if (dictionary == None):
        dictionary = get_replacements()
        
    #---------------------------------
    # Open the template_file to read
    #---------------------------------
    result = glob.glob( cfg_template_file )
    if (len(result) == 0):
        print('SORRY: The template_file:')
        print('   ' + cfg_template_file)
        print('does not exist.')
        return
    template_unit = open( cfg_template_file, 'r' )

    #---------------------------------
    # Open the new_cfg_file to write
    #---------------------------------
    result = glob.glob( new_cfg_file )
    if (len(result) > 0):
        print('WARNING: There is already a CFG file called:')
        print('     ' + new_cfg_file)
        # print '       Are you sure you want to overwrite it?'
        # return
    cfg_unit = open( new_cfg_file, 'w' )

    #------------------------------------------
    # Compile a regex pattern to match "${*}"
    #--------------------------------------------------
    # Use "+" vs. "*"; there must be > 0 chars in "*"
    #--------------------------------------------------
    if (RE_METHOD):
        pattern = re.compile('\$\{[^\}]+\}')
        print('Using RE method...')
    else:
        print('Using STRING method...')

    #---------------------------------------------------------
    # Scan template_file, make replacements & write CFG file
    #---------------------------------------------------------
    while (True):
        line = template_unit.readline()
        if (line == ''):
            break  # (end of file)

        #------------------------------
        # Method that uses re package
        #------------------------------
        if (RE_METHOD):
            matches = re.findall( pattern, line )
            for in_str in matches:
                if (in_str in dictionary):
                    out_str = str( dictionary[ in_str ] )
                    line = line.replace( in_str, out_str )
            cfg_unit.write( line )
            ## cfg_unit.write( line + '\n' )  # (not needed)

        #----------------------------------
        # Method that uses string package
        #---------------------------------------------
        # Use "safe_substitute" to skip over entries
        # without a replacement in the dictionary.
        #---------------------------------------------
        if (STRING_METHOD):
            s = string.Template( line )
            line = s.safe_substitute( dictionary )
            cfg_unit.write( line )
        
    #------------------
    # Close the files
    #------------------
    template_unit.close()
    cfg_unit.close()
    print('Finished creating new CFG file.')
    print(' ')
    
#   replace()
#-------------------------------------------------------------------------



