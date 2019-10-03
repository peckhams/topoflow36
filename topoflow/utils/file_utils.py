
# S.D. Peckham
# June 4, 2010
# July 26, 2010 (file_exists() & count_lines() from tf_utils.)

import glob
## import os     (for os.remove)
import os.path
import numpy
import time

#-------------------------------------------------------------------
#
#   add_extension()
#   replace_extension()
#   get_prefix_and_extension()
#   check_overwrite()
#
#   file_exists()
#   count_lines()
#
#-------------------------------
#   These are not used anymore
#-------------------------------
#   insert_time_stamp()
#   check_overwrite_last()
#
#-------------------------------------------------------------------
def add_extension( file_name, extension='.nc' ):

    #------------------------------------------------------
    # Note: This function is repeated in model_output.py.
    #------------------------------------------------------
    n_ext = len(extension)
    if (file_name[-n_ext:].lower() != extension.lower()):
        return (file_name + extension)
    else:
        return file_name

#   add_extension()
#-------------------------------------------------------------------
def replace_extension( file_name, extension='.nc' ):

    prefix, old_extension = get_prefix_and_extension( file_name )

    return (prefix + extension)

#   replace_extension()
#-------------------------------------------------------------------
def get_prefix_and_extension( file_name ):

    #---------------------------------
    # Find the last dot in file_name
    #---------------------------------
    pos       = file_name.rfind('.')
    prefix    = file_name[:pos]
    extension = file_name[pos:]

    return (prefix, extension)

#   get_prefix_and_extension()
#-------------------------------------------------------------------
def check_overwrite( file_name ):

    #------------------------------------------------------
    # Notes: This routine checks if file_name is already
    #        in use and if so, it returns a new filename
    #        with a "version number" inserted.
    #------------------------------------------------------
    if not(os.path.exists( file_name )):
        return file_name

    #--------------------------------------------------
    # File already exists, so insert a version number
    #--------------------------------------------------
    prefix, extension = get_prefix_and_extension( file_name )
    k = 1
    while (True):
        kstr = '_' + str(k)
        new_name = (prefix + kstr + extension)
        if (os.path.exists( new_name )):
            k += 1
        else:
            break
    return new_name

#   check_overwrite()
#-------------------------------------------------------------------
def file_exists( filename ):

    f = glob.glob( filename )
    if (len(f) == 0):
        print('SORRY, The file:')
        print('  ' + filename)
        print('was not found in the working directory.')
        print(' ')
        found = False
    else:
        found = True

    return found

#   file_exists()
#-------------------------------------------------------------------
def count_lines( filename, SILENT=False ):

    #----------------
    # Open the file
    #----------------
    file_unit = open(filename, 'r')
    
    #------------------
    # Count the lines
    #------------------
    n_lines = numpy.int32(0)
    n_total = numpy.int32(0)
    while (True):
        line = file_unit.readline()
        if (line == ''):
            break
        n_total += 1
        #---------------------------
        # Count the nonblank lines
        #---------------------------
        n_chars = len(line.strip())
        if (n_chars != 0):    
            n_lines += 1
                 
##    line = ''
##    while logical_not(idl_func.eof(file_unit)):
##        line = idl_func.readf(file_unit, line)
##        n_total += 1
##        #---------------------------
##        # Count the nonblank lines
##        #---------------------------
##        _len = len(line.strip())
##        if (_len != 0):    
##            n_lines += 1
    
    #-----------------
    # Close the file
    #-----------------
    file_unit.close()
    
    #--------------------
    # Print a message ?
    #--------------------
    if not(SILENT):    
        print('For the file: ' + filename)
        print('Total number of lines    = ' + str(n_total))
        print('Number of nonblank lines = ' + str(n_lines))
        print(' ')
    
    return n_lines

#   count_lines()
#-------------------------------------------------------------------
##def insert_time_stamp( file_name ):
##
##    prefix, extension = get_prefix_and_extension( file_name )
##    
##    jstr = '-'
##    ## jstr = '.'
##    
##    time_str = jstr.join("%s" % k for k in time.localtime()[:6])
##
##    return (prefix + '_' + time_str + extension)
##            
###   insert_time_stamp()
###-------------------------------------------------------------------
##def check_overwrite_last( file_name ):
##
##    #----------------------------
##    # Does file already exist ?
##    #----------------------------
##    if not(os.path.exists( file_name )):
##        return file_name
##    
##    #----------------------------------------
##    # Append timestamp to make new filename
##    #----------------------------------------
##    print 'File already exists: ' + file_name
##    print 'Inserting timestamp in filename.'
##    new_name = insert_time_stamp( file_name )
##    return new_name
##    
##    #--------------------------------------------
##    # Overwrite not allowed; must remove first.
##    #--------------------------------------------
####    print 'Deleting existing file: ' + file_name
####    os.remove( file_name )
##            
###   check_overwrite_last()
#-------------------------------------------------------------------

