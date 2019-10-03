
## Copyright (c) 2001-2019, Scott D. Peckham
## October 2019  (updated for Python 3)
## January 2009  (converted from IDL)
## November 2009 (collected into cfg_files.py
## May 2010 (added read_key_value_pair())
## July 2010 (added read_list()

import numpy as np

#---------------------------------------------------------------------
#
#   unit_test()
#
#   skip_header()
#   get_yes_words()
#   read_words()
#   read_list()             # (7/27/10)
#
#   read_key_value_pair()   # (5/7/10)
#   read_line_after_key()
#   read_words_after_key()
#   read_list_after_key()
#   read_value()
#
#   var_type_code()
#   read_input_option()
#   read_output_option()    # (boolean and string)
#
#---------------------------------------------------------------------
def unit_test():

    import d8_base
    comp = d8_base.d8_component()
    #------------------------------
    comp.CCA = False
    comp.directory   = '/Applications/Erode/Data/Test1/'
    comp.data_prefix = 'Test1'
    comp.case_prefix = 'Test1'

    comp.read_config_file()
    print('comp.method         =', comp.method)
    print('comp.method_name    =', comp.method_name)
    print('comp.dt             =', comp.dt)
    print('comp.dt.dtype       =', comp.dt.dtype)
    print('comp.LINK_FLATS     =', comp.LINK_FLATS)
    print('comp.LR_PERIODIC    =', comp.LR_PERIODIC)
    print('comp.TB_PERIODIC    =', comp.TB_PERIODIC)
    print('comp.SAVE_DW_PIXELS =', comp.SAVE_DW_PIXELS)
    
    print('Finished with cfg_files.unit_test().')
    print(' ')
    
#   unit_test()
#---------------------------------------------------------------------
def skip_header(file_unit, n_lines=4):

    #-------------------------
    # Skip over header lines
    #-------------------------
    for k in range(n_lines):
        line = file_unit.readline()
            
#   skip_header()
#---------------------------------------------------------------------
def get_yes_words():

    yes_words = ['1', 'true', 'on', 'yes', 'ok'] 
    return yes_words

#   get_yes_words()  
#---------------------------------------------------------------------
def read_words(file_unit, word_delim=None):

    #----------------------------------------------------
    # Note: If (word_delim == None), then the "split()"
    #       method for strings will use any whitespace
    #       string as a separator.
    #----------------------------------------------------
    line = file_unit.readline()

    words = line.split( word_delim )
    return words

#   read_words()
#---------------------------------------------------------------------
def read_list(file_unit, dtype_list=None, dtype='string',
              word_delim=None):
    
    #-------------------------------------------------------------
    # Notes:  Example (read boolean and string/filename):
    #         vlist = read_list_after_key(file_unit,
    #                                     ['boolean', 'string'])
    #-------------------------------------------------------------
    words = read_words(file_unit, word_delim=word_delim)
    
    #----------------------------------------------
    # If "dtype_list" is None, then assume that
    # every element in the list has type "dtype".
    #----------------------------------------------------
    # NB!  "dtype_list" cannot default to "[]", because
    #      then values set from a previous call to this
    #      function are remembered and used.     
    #----------------------------------------------------
##    if (dtype_list == []):
##    if (len(dtype_list) == 0):
    if (dtype_list == None):
        dtype_list = []
        for k in range(len(words)):
            dtype_list.append( dtype.lower() )
    elif (len(dtype_list) > len(words)):
        print('ERROR in cfg_files.read_list_after_key().')
        print('   Not enough values in the line.')
        return

    k = 0
    yes_words = get_yes_words()  
    var_list = []

    for type_str in dtype_list:
        vtype = type_str.lower()
        word  = words[k].strip()
        if   (vtype == 'string'):
            var = word
        elif (vtype == 'boolean'):
            var = (word in yes_words)
        else:
            value = eval( word )
            var = eval( "np." + vtype + '( value )' )  #(2019-10-03)
            ### exec('var = np.' + vtype + '( value )')
        var_list.append( var )
        k += 1
        
    return var_list

#   read_list()
#---------------------------------------------------------------------
def read_key_value_pair(file_unit, key_delim=':', SILENT=True):

    line = file_unit.readline()
  
    #--------------------------------------
    # Extract the variable name or label,
    # which may contain blank spaces
    #--------------------------------------
    p = line.find( key_delim )
    if (p == -1):
        if not(SILENT):
            print('ERROR in cfg_files.read_line_after_key():')
            print('   Key-value delimiter not found.')
        return ('', '')
    
    key   = line[:p]
    value = line[p + 1:]
    value = value.strip() # (strip leading & trailing whitespace)
    return (key, value)

#   read_key_value_pair()
#---------------------------------------------------------------------
def read_line_after_key(file_unit, key_delim=':'):

    line = file_unit.readline()
  
    #--------------------------------------
    # Extract the variable name or label,
    # which may contain blank spaces
    #--------------------------------------
    p = line.find( key_delim )
    if (p == -1):
        print('ERROR in cfg_files.read_line_after_key():')
        print('   Key-value delimiter not found.')
        return ''
    label = line[:p]
    line  = line[p + 1:]
    
    return line.strip()   # (strip leading and trailing whitespace)

#   read_line_after_key()
#---------------------------------------------------------------------
def read_words_after_key(file_unit, key_delim=':',
                         word_delim=None, n_words=None):

    #----------------------------------------------------
    # Note: If (word_delim == None), then the "split()"
    #       method for strings will use any whitespace
    #       string as a separator.
    #----------------------------------------------------
    line = read_line_after_key( file_unit, key_delim=key_delim)
    
    #-------------------------------
    # Extract variables as strings
    #-------------------------------
    words = line.split( word_delim )

    #-----------------------------------
    # Option to check for enough words
    #-----------------------------------
    if (n_words == None):
        return words
    if (len(words) < n_words):
        print('ERROR in read_words_after_key():')
        print('  Not enough words found.')
        return words

#   read_words_after_key()
#---------------------------------------------------------------------
def read_list_after_key(file_unit, dtype_list=None, dtype='string',
                        key_delim=':', word_delim=None):

##    print 'before: dtype =', dtype
##    print 'before: dtype_list =', dtype_list
    
    #-------------------------------------------------------------
    # Notes:  Example (read boolean and string/filename):
    #         vlist = read_list_after_key(file_unit,
    #                                     ['boolean', 'string'])
    #-------------------------------------------------------------
    words = read_words_after_key(file_unit, key_delim=key_delim,
                                 word_delim=word_delim)

    #----------------------------------------------
    # If "dtype_list" is None, then assume that
    # every element in the list has type "dtype".
    #----------------------------------------------------
    # NB!  "dtype_list" cannot default to "[]", because
    #      then values set from a previous call to this
    #      function are remembered and used.     
    #----------------------------------------------------
##    if (dtype_list == []):
##    if (len(dtype_list) == 0):
    if (dtype_list == None):
        dtype_list = []
        for k in range(len(words)):
            dtype_list.append( dtype.lower() )
    elif (len(dtype_list) > len(words)):
        print('ERROR in cfg_files.read_list_after_key().')
        print('   Not enough values in the line.')
        return

##    print 'after: dtype =', dtype
##    print 'after: dtype_list =', dtype_list
##    print '--------------------------------'
    k = 0
    yes_words = get_yes_words()  
    var_list = []

    for type_str in dtype_list:
        vtype = type_str.lower()
        word  = words[k].strip()
        if   (vtype == 'string'):
            var = word
        elif (vtype == 'boolean'):
            var = (word in yes_words)
        else:
            value = eval( word )
            var = eval( "np." + vtype + '( value )' )  # (2019-10-03)
            #### exec('var = np.' + vtype + '( value )')
        var_list.append( var )
        k += 1
        
    return var_list

#   read_list_after_key()
#---------------------------------------------------------------------
def read_value(file_unit, dtype='string', key_delim=':'):

    #--------------------------------------------------------------
    # Notes: Valid "var_types" are:
    #           'file', 'string', 'boolean' and any numpy dtype,
    #        such as:
    #           'uint8', 'int16', 'int32', 'float32', 'float64'

    #  If (var_type eq 'file'), then we want to read everything
    #  after the ":", which may contain space characters in
    #  the interior (e.g a directory), but with leading and
    #  trailing spaces removed.
    #--------------------------------------------------------------
    vtype = dtype.lower()
    if (vtype == 'file'):
        return read_line_after_key(file_unit, key_delim=key_delim)

    words = read_words_after_key( file_unit )

    #--------------------    
    # Return a string ?
    #--------------------
    if (vtype == 'string'):
        return words[0]

    #-------------------------------------
    # Return a boolean (True or False) ?
    #-------------------------------------
    if (vtype == 'boolean'):
        yes_words = get_yes_words()
        return (words[0].lower() in yes_words)
    
    #------------------------------------    
    # Try to convert string to a number
    #------------------------------------
    try:
        value = eval(words[0])
    except:
        print('ERROR in cfg_files.read_value():')
        print('   Unable to convert string to number.')
        return words[0]

    #----------------------------------
    # Return number of requested type
    #----------------------------------
    result = eval( "np." + vtype + '( value )' )  # (2019-10-03)
    ### exec('result = np.' + vtype + '( value )')
    return result
        
#   read_value()
#---------------------------------------------------------------------
def var_type_code(var_type):

    cmap = {'scalar':        0, \
            'time_series':   1, \
            'time series':   1, \
            'grid'       :   2, \
            'grid_stack' :   3, \
            'grid stack' :   3, \
            'grid_sequence': 3, \
            'grid sequence': 3 }

    code = cmap[ var_type.lower() ]
    return np.int16( code )
    
#   var_type_code()
#---------------------------------------------------------------------
def read_input_option(file_unit, key_delim=':', word_delim=None):

    words= read_words_after_key( file_unit, key_delim=key_delim,
                                 word_delim=word_delim )
    
    #-----------------------------------------------
    # TopoFlow "var types" are:
    #     Scalar, Time_Series, Grid, Grid_Sequence
    #-----------------------------------------------
    var_type = words[0].lower()
    if (var_type == 'scalar'):
        type_code = np.int16(0)
        scalar    = np.float64( eval( words[1] ) )
        filename  = ''   # (or use None ??)
    else:
        type_code = var_type_code( var_type )
        scalar    = None
        filename  = words[1]
        
    return (type_code, scalar, filename)
    
#   read_input_option()
#---------------------------------------------------------------------
def read_output_option(file_unit, key_delim=':', word_delim=None):

    #-------------------------------
    # Extract variables as strings
    #-------------------------------
    words = read_words_after_key( file_unit, key_delim=key_delim,
                                  word_delim=word_delim )
    count = len(words)
    if (count == 0):
        print('ERROR in cfg_files.read_output_option():')
        print('   No value found after delimiter.')
        return ''
    if (count == 1):
        print('ERROR in cfg_files.read_output_option():')
        print('   No filename provided after option.')
        return ''

    yes_words = ['1','true','yes','on']
    option    = (words[0].lower() in yes_words)
    filename  = words[1]

    return option, filename

#   read_output_option() 
#---------------------------------------------------------------------


