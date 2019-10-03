#!/usr/bin/env python

# Copyright (c) 2011, Scott D. Peckham

#------------------------------------------------------
# S.D. Peckham
# May 6, 2011
# Tool to convert my ".cfg" (gui info) files to XML.

# Example of use at a Unix prompt:
#
#    % ./cfg2xml.py my_file.cfg model_name
#------------------------------------------------------
#
# Functions:
#
# Test_EOF()
# cfg2xml()
#
#------------------------------------------------------

import os.path
import sys

#------------------------------------------------------
def Test_EOF( text_file ):

    #--------------------------------------------------
    # Recall that a "blank line", with just a (hidden)
    # newline character will not be null and will
    # have len(line) = 1.
    #--------------------------------------------------
    text_unit = open( text_file, 'r' )
    while (True):
        line = text_unit.readline()
        if (line == ''):
            break
        print('len(line) =', len(line))
        print('line = ', line + '_____')
        
#   Test_EOF()
#------------------------------------------------------
def cfg2xml( cfg_file, model_name, word_delim='|' ):

    vp = '/' + model_name + '/Input/Var/'     # (vp = variable prefix)
    project_name = 'TopoFlow'
    
    #--------------------
    # Open the CFG file
    #--------------------
    try:
        cfg_unit = open( cfg_file, 'r' )
    except:
        print('SORRY: Could not open CFG file named:')
        print('       ' + cfg_file)

    #--------------------
    # Open new XML file
    #--------------------
    pos    = cfg_file.rfind('.')
    prefix = cfg_file[0:pos]
    xml_file = prefix + '.xml'
    #----------------------------------------------
    XML_EXISTS = os.path.exists( xml_file )
    if (XML_EXISTS):
        print('SORRY, An XML file with the name')
        print('       ' + xml_file)
        print('       already exists.')
        return
    xml_unit = open( xml_file, 'w' )
    xml_unit.write('<dialog>\n')

    #-----------------------------
    # Read CFG file header lines
    #-----------------------------
    for k in range(4):
        line = cfg_unit.readline()

    FIRST_TAB = True
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
            #---------------------
            # Close previous tab
            #---------------------
            if not(FIRST_TAB):
                xml_unit.write('  </tab>\n\n')
            #-------------------
            # Create a new tab
            #-------------------
            xml_unit.write('<!-- ====================================================================== -->\n')
            xml_unit.write('  <tab name="' + words[1].strip() + '">\n')
            FIRST_TAB = False
        elif (first_word == 'HTML_HELP_FILE'):
            help_url = words[3].strip()
            ### xml_unit.write('  <html_help_file>' + help_url + '</html_help_file>\n')
        elif (first_word == ''):
            pass
        elif (first_word[0] == '#'):
            pass
        else:              
            #------------------
            # Create an entry
            #------------------
            var_name = vp + words[0].strip()
            xml_unit.write('    <entry name="' + var_name + '">\n')
            #---------------------------------------------------------------
            label    = words[1].strip()
            xml_unit.write('      <label>' + label + '</label>\n')
            #---------------------------------------------------------------
            value = words[3].strip()
            if (first_word == 'in_directory'):
                #----------------------------------
                # Change the input file directory     ##########################
                #----------------------------------
                value = value.replace('/data/sims/topoflow', '/home/csdms/models/topoflow/3.1/share/data')
            xml_unit.write('      <default>' + value + '</default>\n')
            #---------------------------------------------------------------  
            type_str = words[2].strip()  # (does first char need to be cap?)
            xml_unit.write('      <type>' + type_str + '</type>\n')
            if (type_str.lower() != 'string'):
                xml_unit.write('      <range>\n')
                min_str = words[4].strip()
                xml_unit.write('        <min>' + min_str + '</min>\n')
                max_str = words[5].strip()
                xml_unit.write('        <max>' + max_str + '</max>\n')
                xml_unit.write('      </range>\n')
            #---------------------------------------------------------------
            help_str = words[6].strip()
            xml_unit.write('      <help_brief>' + help_str + '</help_brief>\n')
            #---------------------------------------------------------------
            xml_unit.write('    </entry>\n\n')

    #---------------------
    # Close the CFG file
    #---------------------
    cfg_unit.close()
    
    #----------------------
    # Create an About tab
    #----------------------
    ABOUT_TAB = True
    if (ABOUT_TAB):
        xml_unit.write('<!-- ====================================================================== -->\n')
        xml_unit.write('  <tab name="About">\n')
        
        xml_unit.write('    <entry name="/' + model_name + '/ModelName">\n')
        xml_unit.write('      <label>Model name:</label>\n')
        xml_unit.write('      <help_brief>Name of the model</help_brief>\n')
        xml_unit.write('      <default>' + model_name + '</default>\n')
        xml_unit.write('      <type>String</type>\n')
        xml_unit.write('    </entry>\n\n')

        xml_unit.write('    <entry name="/' + model_name + '/ModelAuthor">\n')
        xml_unit.write('      <label>Author name:</label>\n')
        xml_unit.write('      <help_brief>Name of the model author</help_brief>\n')
        xml_unit.write('      <default>Scott D. Peckham </default>\n')
        xml_unit.write('      <type>String</type>\n')
        xml_unit.write('    </entry>\n\n')

        xml_unit.write('    <entry name="HTML_HELP_FILE">\n')
        xml_unit.write('      <label>HTML help file:</label>\n')
        xml_unit.write('      <help_brief>URL for HTML help file</help_brief>\n')
        xml_unit.write('      <default>' + help_url + '</default>\n')
        xml_unit.write('      <type>String</type>\n')
        xml_unit.write('    </entry>\n\n')
        
        xml_unit.write('  </tab>\n')
            
    #---------------------
    # Close the XML file
    #---------------------
    xml_unit.write('</dialog>\n\n')
    xml_unit.close()
    print('Finished writing TopoFlow CFG file to XML format.')
    print(' ')
    
#   cfg2xml()
#-----------------------------------------------------------------------
if (__name__ == "__main__"):
    
    #-----------------------------------------------------
    # Note: First arg in sys.argv is the command itself.
    #-----------------------------------------------------
    n_args = len(sys.argv)
    if (n_args < 3):
        print('ERROR: This tool requires CFG filename')
        print('       and model_name arguments.')
        print('sys.argv =', sys.argv)
        print(' ')
    elif (n_args == 3):
        cfg2xml( sys.argv[1], sys.argv[2] )
    else:
        print('ERROR: Invalid number of arguments.')
        
#-----------------------------------------------------------------------

        
