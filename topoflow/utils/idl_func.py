
#   idl_func.py

#   Author:   Dr. Scott D. Peckham, INSTAAR, Univ. of Colorado
#   Created:  August 14, 2008

from numpy import *
import numpy
import os
import os.path
import re
import platform
 
#------------------------------------------
#  Define the "bunch" class for structures
#----------------------------------------------
#  Note:  IDL allows fields in a structure to
#  be accessed by their "tag index", as in:
#  <struct>.(<index>).
#  We can get the same thing here as follows:
#     y = bunch(a=5.0, b=1.0)
#     print y.__dict__.values()[0]
#----------------------------------------------
class bunch:
    """Used by I2PY to simulate IDL structures"""
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
       
#------------------------
#  IDL's BYTE function
#------------------------
def byte(s):
    """Reproduces IDL's BYTE function"""
    #---------------------------------------------------------
    # Note:  IDL's BYTE function is overloaded.
    #        If arg is a string, then a byte array of ASCII
    #           ordinal values is returned.
    #        Otherwise, if arg is scalar or array, type
    #           is converted to (unsigned) byte.  
    #---------------------------------------------------------    
    #        If called with more than 1 arg, then it can
    #        be used to extract 1 byte of data from 1st arg.
    #        This part is not supported yet.
    #---------------------------------------------------------
    # is_string = not(str(s).isdigit)  # (not as general)
    is_string = (type(s) == type('abc'))
    if (is_string):
        result = array(list(map(ord, s)), copy=0)
    else:
        result = array(s, copy=0).astype('UInt8')
        # Is "mod" redundant in next line ??
        # result = array(mod(s,256), copy=0).astype('UInt8')
        #----------------------------
        # This crashes Python/NumPy
        #----------------------------
        # result = mod(array(s, copy=0).astype('UInt8'), 256)
        # result = array(s, copy=0).astype('UInt8') % 256
    return result


#--------------------------
#  IDL's BYTSCL function
#--------------------------
def bytscl(a, max=None, min=None, top=255):
    """Reproduces IDL's BYTSCL function"""
    amax = max
    amin = min
    if (amax is None): amax = a.max()
    if (amin is None): amin = a.min()
    top = (top % 256)  # (force to byte range)

    #--------------------------------------
    # Case of integer or float argument ?
    #--------------------------------------
    type_str = str(a.dtype)
    if (type_str[0] == 'i'): 
        b = ((top + 1) * (a - amin) - 1)/(amax - amin)
    else:
        b = (top + 0.9999) * (a - amin) / (amax - amin)
    return b

#---------------------------------
#  IDL's !d.name system variable
#---------------------------------
def device_name():
    """Reproduces IDL's !d.name system variable"""
    #---------------------------------------------------
    # Note: Checked that !d.name on new Mac is 'X'.
    #       Presumably will get 'X' on any Unix/Linus.
    #---------------------------------------------------
    sys_name = platform.system()
    if   (sys_name == 'Windows'): return 'WIN'
    elif (sys_name == 'Darwin'):  return 'X'
    else: return 'X'

#---------------------------------------
#  IDL's EOF function (no longer used)
#---------------------------------------
def eof(file_obj):
    """Reproduces IDL's EOF function"""
    #-------------------------------------------
    #  This recomputes file_size every time it
    #  is called (as in a while loop), so it
    #  has an extra cost.  It is also not the
    #  "Python way", but seems to work.
    #-------------------------------------------
    file_size = os.path.getsize(file_obj.name)
    file_pos  = file_obj.tell()
    return (file_pos == file_size)

#-------------------------------
#  IDL's FILE_DELETE procedure
#-------------------------------
def file_delete(*files):
    """Reproduces IDL's FILE_DELETE procedure"""
    for path in files:
        print(path)
        if (os.path.exists(path)):
            if (os.path.isdir(path)):
                os.rmdir(path)
            else:
                os.remove(path)
                
#------------------------
#  IDL's FSTAT function
#------------------------
def fstat(file_obj):
    """Reproduces IDL's FSTAT function"""
    #-----------------------------------------------
    # Note:  For IDL routines like FSTAT that take
    #        the file unit number as an argument,
    #        we instead pass a Python file object,
    #        with a name constructed from the unit
    #        variable name.
    #-----------------------------------------------
    # Note:  Unsupported fields are:
    #        uisagui, interactive, xdr, compress,
    #        atime, ctime, mtime, transfer_count,
    #        cur_ptr, rec_len.
    #-----------------------------------------------
    #        If file is open for update, mode may
    #        be 'U', and both "read" and "write"
    #        should probably be True. ???
    #-----------------------------------------------
    file_size = os.path.getsize(file_obj.name)
    stat = bunch(size=file_size, \
                 name=file_obj.name, \
                 open=not(file_obj.closed), \
                 read=('r' in file_obj.mode), \
                 write=('w' in file_obj.mode), \
                 unit=(file_obj.fileno()), \
                 isatty=file_obj.isatty() )
    return stat

#------------------------------
#  IDL's KEYWORD_SET function
#   (this one not used now)
#------------------------------
def keyword_set(arg):
    """Simulates IDL's KEYWORD_SET function"""
    return (arg is not None) and (arg != 0)

#-----------------------------
#  IDL's N_ELEMENTS function
#-----------------------------
def n_elements(arg):
    """Simulates IDL's N_ELEMENTS function"""
    
    zero = numpy.int32(0)      # returns a numpy object
    if ('arg' in locals()):
        #--------------------------------------------------------
        # Note: An unset IDL keyword will be assigned a default
        #       value of "None" by I2PY and will therefore be
        #       found in locals().  We want to return 0 in this
        #       case as well.  (numpy.size(None) equals 1.)
        #--------------------------------------------------------
        if (arg is None):
            return zero
        else:
            return numpy.size(arg)
    else:
        return zero
    
#------------------------
#  IDL's READF function
#------------------------
def readf(file_obj, *args, **keys):
    """Simulates IDL's READF function"""
    #-----------------------------------------
    # Note: Also check out the "csv" module.
    #-----------------------------------------
##    args_out = list(args)  # to allow assignment
 
##    args_out = []
##    for aa in args:
##        args_out.append(aa)

    str_type = type('abc')
    
    if not('format' in keys):    
        line = file_obj.readline()
        if (len(args) == 1) and (type(args[0]) == str_type):
            #-----------------------------------------
            # Read entire line into a string
            # This assumes FORMAT keyword is not set
            #----------------------------------------------
            # 12/2/08. Bug fix. Return line itself to
            # avoid changing "line" from 'str' to 'list'
            #---------------------------------------------
            return line
            ## args_out[0] = line
        else: 
            args_out = list(args)  # to allow assignment
            vals = line.split()
            
            for k in range(len(args)):
                #---------------------------------------------
                # This works for scalar numbers and strings,
                # but doesn't work for arrays yet.
                #---------------------------------------------
                ## print 'type(args[k]) =', type(args[k])  #########
                
                if (type(args[k]) == str_type):
                    args_out[k] = numpy.array(vals[k])
                else:
                    dtype = str(args[k].dtype)
                    # print 'dtype =', dtype
                    
                    if (dtype[:2] == '|S'):
                        args_out[k] = numpy.array(vals[k])
                    else:
                        args_out[k] = numpy.array(vals[k]).astype(dtype)
                    # print 'args_out[k] =', args_out[k]
    else:
        #-----------------------------------------------------
        # Assume values to read are separated by white space
        #-----------------------------------------------------
        # Use the string() function below somehow ??
        #----------------------------------------------
        args_out = list(args)  # to allow assignment
        format = keys['format']
        format = format.replace('(','')
        format = format.replace(')','')
        format = format.replace(' ','')
        parts  = format.split(',')
        for k in range(len(args)):
            #-----------------------------------
            # Extract n from formatting string
            #-----------------------------------------------
            # There shouldn't be any "x" codes for reading ?
            #-----------------------------------------------
            # if ('x' in parts[k]): ...
            r = re.split('[a-zA-Z]', parts[k])
            n = int(eval(r[1]))
            val_str = file_obj.read(n)
            dtype = str(args[k].dtype)
            if (dtype[:2] == '|S'):
                args_out[k] = numpy.array(val_str)
            else:
                args_out[k] = numpy.array(val_str).astype(dtype)

##    print 'args_out =', args_out
##    print ' '

    # Need this.               
    if (len(args) == 1):
        return args_out[0]
    else:
        return args_out

#------------------------
#  IDL's READS function
#------------------------
def reads(line, *args, **keys):
    """Simulates IDL's READS function"""
    #-----------------------------------------
    # Note: Also check out the "csv" module.
    #-----------------------------------------
    args = list(args)  # to allow assignment
    if not('format' in keys):    
        if (len(args) == 1) and (type(args[0]) == type('abc')):
            #-----------------------------------------
            # Read entire line into a string
            # This assumes FORMAT keyword is not set
            #-----------------------------------------
            args[0] = line
        else:
            vals = line.split()
            for k in range(len(args)):
                #---------------------------------------------
                # This works for scalar numbers and strings,
                # but doesn't work for arrays yet.
                #---------------------------------------------
                dtype = str(args[k].dtype)
                if (dtype[:2] == '|S'):
                    args[k] = numpy.array(vals[k])
                else:
                    args[k] = numpy.array(vals[k]).astype(dtype)
    else:           
        #-------------------------------------
        # Use the string() function below ??
        #-------------------------------------
        format = keys['format']
        format = format.replace('(','')
        format = format.replace(')','')
        format = format.replace(' ','')
        parts  = format.split(',')
        pos = 0
        for k in range(len(args)):
            #-----------------------------------
            # Extract n from formatting string
            #-----------------------------------------------
            # There shouldn't be any "x" codes for reading ?
            #-----------------------------------------------
            # if ('x' in parts[k]): ...
            r = re.split('[a-zA-Z]', parts[k])
            n = int(eval(r[1]))
            val_str = line[pos: pos+n-1]
            pos += n
            dtype = str(args[k].dtype)
            if (dtype[:2] == '|S'):
                args[k] = numpy.array(val_str)
            else:
                args[k] = numpy.array(val_str).astype(dtype)
            
    return args

#--------------------------------------------
#  IDL's SINDGEN function  (no longer used)
#  (only 1 argument supported now)
#--------------------------------------------
def sindgen(n):
    """Reproduces IDL's SINDGEN function"""
    a = arange(n, dtype='Int32')
    return list(map(str, a))

##    s = []
##    for k in a: s.append(str(a[k]).rjust(12))
##    return s

#---------------------------------------
#  IDL's SIZE function  (not used yet)
#---------------------------------------
def size(a, n_dimensions=False, n_elements=False, dimensions=False, \
         _type=False):
    """Reproduces IDL's SIZE function"""
    #---------------------------------------------------------
    #  By default, IDL's SIZE function returns a vector of
    #  values that describe an object's dimensions and type.
    #---------------------------------------------------------
    ndim_a = numpy.ndim(a)
    if (n_dimensions): return ndim_a
    if (n_elements):   return numpy.size(a)
    if (dimensions):
        s  = numpy.zeros(ndim_a, dtype='int32')
        #----------------------------------------
        # IDL returns dimension in reverse order
        #----------------------------------------
        for k in range(ndim_a):
            j    = (ndim_a - k - 1)
            s[j] = numpy.size(a, k)
        return s
    if (_type):
        type_code = {'uint8':1, 'int16':2, 'int32':3, 'float32':4, \
                     'float64':5, 'complex32':6, 'str':7}
        type_str = str(a.dtype)
        if (type_str[:2] == '|S'): type_str = 'str'  ###
        return type_code[type_str]
        
    #---------------------------------------------------
    # Otherwise, return an integer array like IDL does
    # NB!  IDL returns dimension in reverse order.
    #---------------------------------------------------
    s = numpy.zeros(ndim_a + 3, dtype='int32')
    s[0] = ndim_a
    for k in range(ndim_a):
        j    = (ndim_a - k)   # (now, don't subtract 1)
        s[j] = numpy.size(a, k)

    type_code = {'uint8':1, 'int16':2, 'int32':3, 'float32':4, \
                 'float64':5, 'complex32':6, 'str':7}
    type_str = str(a.dtype)
    if (type_str[:2] == '|S'): type_str = 'str'  ###
    s[ndim_a + 1] = type_code[type_str]  
    s[ndim_a + 2] = numpy.size(a)
    return s   

#-------------------------
#  IDL's STRING function
#-------------------------
def string(*args, **keys):
    """Reproduces IDL's STRING function"""
    #-------------------------------------------------------------
    # Note: Now supports case where there is a single,
    #       byte array argument and converts to string. i.e.
    
    #       a = numpy.array([72,101,108,108,111], dtype='uint8')
    #       string(a) returns "Hello"
    #-------------------------------------------------------------
    #       Only the FORMAT keyword is supported so far.
    #       IDL's STRING has the additional keywords:
    #       AM_PM, DAYS_OF_WEEK, MONTHS, STDIO_NON_FINITE
    #-------------------------------------------------------------
    #print 'keys =', keys  # equals "{}" if no keys
    FORMATTED = 'format' in keys
    if (FORMATTED):
        format = keys['format']
        if (format == None): FORMATTED=False
        
##    print 'format =', format
##    print 'type(format) =', type(format)   ########
    
    if not(FORMATTED):
        SINGLE  = (len(args) == 1)
        NDARRAY = (str(type(args[0])) == "<type 'numpy.ndarray'>")
        if (SINGLE and NDARRAY):
            BYTE_TYPE = (args[0].dtype == 'uint8')
            if (BYTE_TYPE):
                result_str = ''.join(map(chr, args[0]))
            else:
                result_str = '  '.join(map(str, args))
        else:
            result_str = '  '.join(map(str, args))
        return result_str

    #--------------------------------------------------
    # Convert IDL formatting string to Python version
    #--------------------------------------------------
    # First, remove all parentheses, not just ones at
    # the beginning and end.  Remove all spaces, too.
    #--------------------------------------------------
    # Convert "args" from tuple to list.  We need it
    # to be mutable so we can use insert() method to
    # handle case of "4x", etc.  But near the end we
    # need a tuple again.
    #--------------------------------------------------
    ###  format = str(format)  #######
    
    args  = list(args)
    f_str = format.replace(')','')   ###
    f_str = f_str.replace('(','')
    f_str = f_str.replace(' ','')
    f_str = f_str.lower()            ###

    parts = f_str.split(',')
    for m in range(len(parts)):
        s   = parts[m]
        pos = s.find('x')
        if (pos == -1):
            if (s[0].isalpha()):
                #-----------------------------------------
                # Move "format letter" from start to end
                #-----------------------------------------
                parts[m] = s[1:len(s)] + s[0]
            else:
                #--------------------------------------               
                # Process a "repetition count" before
                # the "format letter".
                #--------------------------------------
                rep_str = ''
                len_str = ''
                FOUND = False
                for c in s:
                    if not(FOUND):
                        if c.isalpha():
                            code = c
                            FOUND = True
                        else:
                            rep_str += c
                    else:
                        len_str += c
                nreps = eval(rep_str)
                parts[m] = '%'.join(numpy.repeat(len_str + code, nreps))
        else:      
            #-----------------------------------------
            # Insert "blank" strings into i to match
            # places with "<n>x" string formatting
            #-----------------------------------------
            nstr  = s.replace('x','')
            nstr  = nstr.replace(' ','')
            blank = ''.ljust(eval(nstr))
            args.insert( m+1, blank )          ### (m+1) ###

    # print 'args =', args
    # print 'len(args) =', len(args)
    
    f_str = "%" + "%".join(parts)
    f_str = f_str.replace('d', 'f')  # for floats
    f_str = f_str.replace('i', 'd')  # for integers
    f_str = f_str.replace('a', 's')  # for strings
    f_str = f_str.replace('x', 's')  # see above        

##    print 'f_str =', f_str
##    print 'args  =', args

    result_str = f_str % tuple(args)   
    return result_str

