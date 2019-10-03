
## Copyright (c) 2001-2010, Scott D. Peckham
## January 2009  # (converted from IDL)
## August 2009
## May 2010 (renamed *Data_Prefix to *Site_Prefix)

#-------------------------------------------------------------------

#  Functions:
#      TF_Use_CCA()
#      TF_Set_Test_Info()       # (6/23/10)
#      TF_Test_Directory()
#      TF_Test_Site_Prefix()
#      TF_Test_Case_Prefix()
#      TF_Build_Date()
#      TF_Version_Number()
#      TF_Version()
#      TF_Print()
#      TF_String()

#      file_exists() # (2/16/09, used in diversions_fraction_method.py)

#      Count_Lines()
#      Current_Directory()
#      Resize()
#      Make_Savefile()
#      Get_site_prefix()
#      Get_case_prefix()  # (gets next available case number)
#      Not_Same_Byte_Order()
#      Courant_Condition_OK()
#      Stable_Timestep()
#      TF_Tan()

#      Convert_Flow_Grid()  # (move this to d8_base.py)

#-------------------------------------------------------------------
from numpy import *
import numpy
import platform, sys, os, glob  # (all are used)

from . import idl_func

#-------------------------------------------------------------------
def TF_Use_CCA():

    return True    # (to use as a set of CCA components)
    ## return False   # (to use as stand-along Python program)

#   TF_Use_CCA()
#-------------------------------------------------------------------
def TF_Set_Test_Info(self):

    system = platform.system()

    self.site_prefix = 'Treynor'
    self.case_prefix = 'Case5'
    self.cfg_prefix  = 'Case5'
    
    if (system == 'Darwin'):
        #----------------------------
        # For my MacBook and MacPro
        #----------------------------
        self.in_directory  = '/Applications/TopoFlow/Data/Treynor_Iowa/'
        # self.in_directory = '/Applications/TopoFlow/Data/Animas/'
        self.out_directory = '~/CMT_Output/'
    elif (system == 'Linux'):
        self.in_directory   = '/data/sims/topoflow/treynor_iowa/'
        # self.in_directory = '/data/sims/topoflow/animas/'
        self.out_directory  = '~/CMT_Output/'
    else:
        #-------------------------------------------------
        # Can "system" be either "Windows" or "Vista" ??
        #-------------------------------------------------
        self.in_directory  = 'C:\Program Files\TopoFlow\Data\Treynor_Iowa\\'
        # self.in_directory = 'C:\Program Files\TopoFlow\Data\Animas\\'       
        self.out_directory = self.in_directory

    self.cfg_directory = self.in_directory
        
#   TF_Set_Test_Info()
#-------------------------------------------------------------------
def TF_Test_Directory():

    system = platform.system()
    
    if (system == 'Darwin'):
        #----------------------------
        # For my MacBook and MacPro
        #----------------------------  
        return '/Applications/TopoFlow/Data/Treynor_Iowa/'
        ## return '/Applications/TopoFlow/Data/Animas/'
    elif (system == 'Linux'):
        return '/data/progs/topoflow/3.0/data/treynor_iowa/'
        ## return '/data/progs/topoflow/3.0/data/animas/'
        ## return '/home/beach/faculty/hannonmt/TFtest/'
    else:
        #-------------------------------------------------
        # Can "system" be either "Windows" or "Vista" ??
        #-------------------------------------------------
        return 'C:/Program Files/TopoFlow/Data/Treynor_Iowa/'
        ## return 'C:/Program Files/TopoFlow/Data/Animas/'       

#   TF_Test_Directory()
#-------------------------------------------------------------------
def TF_Test_Site_Prefix():

    ## return 'sb_nger'
    return 'Treynor'
    ## return 'Animas_200'

#   TF_Test_Site_Prefix()
#-------------------------------------------------------------------
def TF_Test_Case_Prefix():

    return 'Case5'

#   TF_Test_Case_Prefix()
#-------------------------------------------------------------------   
def TF_Build_Date():

    #---------------------------------------------------------
    # Notes:  Update this whenever a new version is released
    #---------------------------------------------------------
    return '11/9/11'

#   TF_Build_Date()
#-------------------------------------------------------------------
def TF_Version_Number():

    return 3.1

#   TF_Version_Number()
#-------------------------------------------------------------------
def TF_Version():

    date_str   = '(' + TF_Build_Date() + ')'
    num_str    = str(TF_Version_Number()) + ' beta '
    ver_string = 'TopoFlow Version ' + num_str + date_str

    return ver_string

#   TF_Version()
#-------------------------------------------------------------------
def TF_Print(string):
    """Print a string to the TopoFlow output log window."""

    #------------------------------------------
    # For now, just print to terminal window.
    # Later, this will be method in tf_gui.py.
    #------------------------------------------
    print(string)

#   TF_Print()  
#-------------------------------------------------------------------
def TF_String(number, FORMAT=None):

    return idl_func.string(number, format=FORMAT).strip()

#   TF_String()
#-------------------------------------------------------------------
def file_exists(filename):

    f = glob.glob(filename)
    if (len(f) == 0):
        print('SORRY, The file:')
        print('  ' + filename)
        print('was not found in the working directory.')
        print(' ')
        found = False
    else:
        found = True

    return found

#   find_file()
#-------------------------------------------------------------------
def Count_Lines(filename, SILENT=False):

    #----------------
    # Open the file
    #----------------
    file_unit = open(filename, 'r')
    
    #------------------
    # Count the lines
    #------------------
    n_lines = int32(0)
    n_total = int32(0)
    line = ''
    while logical_not(idl_func.eof(file_unit)):
        line = idl_func.readf(file_unit, line)
        n_total = (n_total + int32(1))
        #---------------------------
        # Count the nonblank lines
        #---------------------------
        _len = len(line.strip())
        if (_len != 0):    
            n_lines = (n_lines + int32(1))
    
    #-----------------
    # Close the file
    #-----------------
    file_unit.close()
    
    #--------------------
    # Print a message ?
    #--------------------
    if not(SILENT):    
        print('For the file: ' + filename)
        print('Total number of lines    = ' + TF_String(n_total))
        print('Number of nonblank lines = ' + TF_String(n_lines))
        print(' ')
    
    return n_lines

#  Count_Lines
#-------------------------------------------------------------------
def Current_Directory():

    #-----------------------------------------------------------
    #NOTES:  IDL's CD routine behaves differently on different
    #        platforms.  This routine provides a wrapper around
    #        CD to avoid this problem.  On Macs, the directory
    #        separator is ":" and it is automatically appended.

    #        With the separator, filepaths can be constructed
    #        as follows:   filepath = (directory + filename)
    #-----------------------------------------------------------   
    directory = os.getcwd()
    
    #---------------------------------
    # Append a directory separator ?
    #---------------------------------
    I2PY_expr = idl_func.device_name()
    if I2PY_expr == 'WIN':    
        directory = (directory + '\\')
    elif I2PY_expr == 'X':    
        directory = (directory + '/')
    elif I2PY_expr == 'MAC':    
        dum = int16(0)   #(no action needed on Macs)
    else:    
        dum = int16(0)
      
    return directory
#    Current_Directory
#-------------------------------------------------------------------
def Resize(array, factor, REDUCE, new_ncols, new_nrows, SAMP_TYPE=None):

    #-------------------------------------------------------------
    #NOTES:  This routine allows any 2D array to be "rebinned"
    #        to n times larger or smaller than its original size,
    #        where n is a positive integer.

    #        If (REDUCE eq 1b) then the array is reduced in size
    #        by sampling at regular, equal-spaced intervals.

    #        If (REDUCE eq 0b) then the array is enlarged in size
    #        by pixel replication.

    #        This routine reduces to IDL's REBIN routine if:
    #           (1) (REDUCE eq 0), or
    #           (2) (REDUCE eq 1) AND the requested array
    #                dimensions are multiples of the original
    #                dimensions.
    #        When these conditions are not met, the original
    #        array is cropped by the smallest amount necessary
    #        to meet them.  This amount will never exceed n.

    #        Resizing by arbitrary scale factors is avoided

    #        because it distorts the data too much.
    #-------------------------------------------------------------

    if (factor <= 0):    
        GUI_Error_Message(array(['Scale factor must be greater than zero.']))
        return 0
    
    if (SAMP_TYPE is None):    
        SAMP_TYPE = int16(0)
    
    #--------------------------
    # Get dimensions of array
    #--------------------------
    s = idl_func.size(array, dimensions=True)
    ncols = s[0]
    nrows = s[1]
    
    if (factor != 1):    
        if (REDUCE):    
            new_ncols = maximum(int32(ncols / factor), int32(1))
            new_nrows = maximum(int32(nrows / factor), int32(1))
            
            x_rem = (ncols % new_ncols)
            y_rem = (nrows % new_nrows)
            #print,'x_rem = ' + string(x_rem)
            #print,'y_rem = ' + string(y_rem)
            #print,' '
            
            
            if (logical_and((x_rem == 0), (y_rem == 0))):    
                a = idl_func.rebin(array, new_ncols, new_nrows, SAMPLE=SAMP_TYPE)
            else:    
                a = idl_func.rebin(array[0:((nrows - int32(1) - y_rem))+1,0:((ncols - int32(1) - x_rem))+1], new_ncols, new_nrows, SAMPLE=SAMP_TYPE)
        else:    
            new_ncols = int32(factor) * ncols
            new_nrows = int32(factor) * nrows
            a = idl_func.rebin(array, new_ncols, new_nrows, SAMPLE=SAMP_TYPE)
    else:    
        new_ncols = ncols
        new_nrows = nrows
        return array
    
    return a
#    Resize
#-------------------------------------------------------------------
def Make_Savefile():
  
    #--------------------------------------
    # Get name and location for save file
    #--------------------------------------
    save_file = I2PY_filepath = []
    app = wx.PySimpleApp()
    I2PY_dialog = wx.FileDialog(parent=None, defaultDir=os.getcwd(), defaultFile='topoflow.sav', style=wx.SAVE, wildcard='(*.sav)|*.sav')
    if (I2PY_dialog.ShowModal() == wx.ID_OK):
        I2PY_filepath.append(I2PY_dialog.GetPath())
    I2PY_dialog.Destroy()
    if (save_file == ''):    
        return _ret()
    
    #----------------------
    # Create the SAV file
    #----------------------
    print(' ')
    print('Saving TopoFlow model code to: ')
    print(save_file + '...')
    
    save(filename=save_file, routines=True)
    print('Finished.')
    print(' ')
    
    #----------------------------------
    # Print the size of the save file
    #----------------------------------
    file_unit = open(save_file, 'r')
    info = idl_func.fstat(file_unit)
    file_unit.close()
    print('Size of new SAV file is ' + TF_String(info.size) + ' bytes.')
    print(' ')
    
#  Make_Savefile
#-------------------------------------------------------------------
def Get_Site_Prefix(filepath):
 
    #----------------------------------
    # Get file separator for platform
    #----------------------------------
    I2PY_expr = idl_func.device_name()
    if (I2PY_expr == 'WIN'):    
        sep = uint8(92)    #byte('\')
    elif (I2PY_expr == 'X'):    
        sep = uint8(47)    #byte('/')
    elif (I2PY_expr == 'MAC'):    
        sep = uint8(58)    #byte(':')
    else:    
        sep = uint8(92)
       
    #--------------------------------------
    # Extract filename from filepath as
    # everything after last dir separator
    #--------------------------------------
    b = idl_func.byte(filepath)
    w = where(b == sep)
    ns = size(w[0])
    if (ns != 0):    
        start = w[ns - int32(1)] + int32(1)
        filename = filepath[start:]
    else:    
        filename = filepath
    
    #-------------------------------------
    # Extract prefix from filename as
    # everything up to underscore or dot
    #-------------------------------------
    b = idl_func.byte(filename)
    w = where(b == uint8(95))
    nw = size(w[0])
    if (nw == 0):    
        w = where(b == uint8(46))
        nw = size(w[0])  #(look for dots)
    if (nw != 0):    
        p = w[nw - 1]
        prefix = filename[0:0+p]
    else:    
        prefix = filename
        
    return (prefix, filename)

#   Get_Site_Prefix
#-------------------------------------------------------------------
def Get_Case_Prefix():

    #--------------------------------------------------------------
    #Notes:   7/13/06.  Get a default case prefix as the lowest
    #         unused case number.  This currently handles numbers
    #         1 to 999.

    #         3/8/07.  Fixed bug by using a new method.
    #--------------------------------------------------------------

    #------------------------------------------
    # Doesn't matter now whether FILE_SEARCH
    # returns sorted names, but it seems to
    #------------------------------------------
    digits = array(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])
    files = glob.glob('Case*')
    
    #print,'files = ', files
    #;files  = files[sort(files)]
    #;print,'files = ', files
    
    #-------------------------------
    # Find the maximum case number
    #-------------------------------
    max_case = -1
    ## n_files = idl_func.n_elements(files)
    n_files = size(files)
    for k in arange(int32(0), n_files):
        #--------------------------
        # Get characters 5, 6 & 7
        #--------------------------
        char5 = files[k][4:4+1]
        char6 = files[k][5:5+1]
        char7 = files[k][6:6+1]
        #--------------------------------
        w1 = where(digits == char5)
        n1 = size(w1[0])
        if (n1 != 0):    
            _str = char5
            w2 = where(digits == char6)
            n2 = size(w2[0])
            if (n2 != 0):    
                _str = _str + char6
                w3 = where(digits == char7)
                n3 = size(w3[0])
                if (n3 != 0):    
                    _str = _str + char7
            max_case = (maximum(max_case, int16(_str)))
    
    #---------------------------------------------
    # Construct case prefix from max case number
    #---------------------------------------------
    if (max_case == -1):    
        case_prefix = 'Case1'
    else:    
        case_prefix = 'Case' + TF_String(max_case + 1)
    
    #print, 'case_prefix = ' + case_prefix
      
    return case_prefix

#   Get_Case_Prefix()
#-------------------------------------------------------------------
def Not_Same_Byte_Order(byte_order):
  
    if (idl_func.n_elements(byte_order) == 0):    
        print('ERROR: Byte order argument is required in')
        print('        the Not_Same_Byte_Order function.')
        print(' ')
        return -1
    
    #------------------------------------
    # Get byte order of user's computer
    #------------------------------------
    big_endian = (sys.byteorder == 'big')
    if (big_endian):    
        machine_byte_order = 'MSB'
    else:    
        machine_byte_order = 'LSB'
    
    #-----------------------------------------
    # Compare machine byte order to the byte
    # order of the current binary grid.
    #-----------------------------------------
    NOT_SAME = (machine_byte_order != byte_order)
    
    return NOT_SAME

#  Not_Same_Byte_Order
#-------------------------------------------------------------------
def Courant_Condition_OK(dx, dt, vmax):

    #-------------------------------------------------------------
    # Notes: This condition for numerical stability dictates that
    #        the maximum distance travelled by water anywhere in
    #        the basin in one timestep (vmax * dt) must be less
    #        than one pixel width (dx).

    #        vmax will often be an estimate.
    #        vmax has units of meters/second.
    #        dx has units of meters & dt has units of seconds.

    #        Typically, vmax will have a built-in factor
    #        of safety, so set factor = 1d here.

    # NB!    We must use LE vs. LT here and the same factor
    #        as in Stable_Timestep function.
    #-------------------------------------------------------------
    factor = float64(1)
    
    OK = ((vmax * dt * factor) <= dx)
    
    return OK
    
#   Courant_Condition_OK
#-------------------------------------------------------------------
def Stable_Timestep(dx, vmax):

    #------------------------------------------------------
    # Notes: See notes for Courant_Condition_OK function.
    #        Typically, vmax will have a built-in factor
    #        of safety, so set factor = 1d here.
    #-----------------------------------------------------
    factor = float64(1)
    dt = dx / (vmax * factor)
    
    return dt
    
#   Stable_Timestep
#-------------------------------------------------------------------
def TF_Tan(angle):

    return tan(angle)

##        #----------------------------------------------------------
##        # Notes: IDL for 64-bit Linux has problem with the TAN
##        #        function so use SIN()/COS() instead  (B. Bolton)
##
##        #        Seems to only be a problem for 64-bit Linux, so
##        #        this function tests for that platform and uses
##        #        TAN() for all other platforms.
##        #----------------------------------------------------------
##        LINUX_OS = (platform.system().upper() == 'LINUX')
##        
##        if (LINUX_OS):    
##            return sin(angle) / cos(angle)
##        else:    
##            return tan(angle)
    
#  TF_Tan
#-------------------------------------------------------------------    
def Convert_Flow_Grid(infile, outfile, itype, otype, byte_order, \
                      icodes=None, ocodes=None, REPORT=False):

    #--------------------------------------------------------
    # NOTES:  The input filename "infile" and output filename
    #         "outfile" must be in single quotes.

    #         Any flow codes other than incodes get mapped
    #         to zero.

    #         Conversion is applied in blocks for speed.

    #         RiverTools flow codes are:
    #         codes = [1,   2,  4,   8, 16, 32, 64, 128]

    #         ARC/INFO flow codes are:
    #         codes = [128,  1,  2,  4,   8, 16, 32, 64]

    #         TOPAZ flow codes are:
    #         codes = ???????????
    #--------------------------------------------------------

    n_params = 7 - [icodes, ocodes].count(None)
    
    #------------------------------
    # Is outfile same as infile ?
    #------------------------------
    if (infile == outfile):    
        msg = array(['SORRY, ', ' ', \
                     'You must use different names for the input ', \
                     'and output files. ', ' '])
        GUI_Error_Message(msg)
        return _ret()
    
    #--------------------------------------------
    # Default is to convert from ARC flow codes
    #--------------------------------------------
    if (idl_func.n_elements(icodes) == 0):    
        icodes = array([128, 1, 2, 4, 8, 16, 32, 64])
    
    #-----------------------------------------
    # Default is to convert to RT flow codes
    #-----------------------------------------
    if (idl_func.n_elements(ocodes) == 0):    
        ocodes = array([1, 2, 4, 8, 16, 32, 64, 128])
    
    #-------------------------------------
    # Convert type strings to upper case
    #-------------------------------------
    itype = itype.upper()
    otype = otype.upper()
    
    #-----------------
    # Create the map
    #-----------------
    _map = zeros([255], dtype='UInt8')
    _map[icodes] = ocodes
    
    #-----------------
    # Open the files
    #-----------------
    file_iunit = open(infile, 'rb')
    file_ounit = open(outfile, 'wb')
    SWAP_ENDIAN = Not_Same_Byte_Order(byte_order)
    print('Writing new flow grid to ' + outfile + '...')
    
    #-------------------------------------
    # Compute number of blocks in infile
    #-------------------------------------
    blocksize = int32(1024)
    temp = idl_func.fstat(file_iunit)
    filesize = temp.size
    n_blocks = (filesize / blocksize)
    rem_size = (filesize % blocksize)
    
    if (REPORT):    
        print(' > blocksize = ' + TF_String(blocksize))
        print(' > filesize  = ' + TF_String(filesize))
        print(' > n_blocks  = ' + TF_String(n_blocks))
        print(' > rem_size  = ' + TF_String(rem_size))
        print(' >')
    
    #------------------------
    # Initialize block array
    #------------------------
    if itype == 'BYTE':    
        block = zeros([blocksize], dtype='UInt8')
    elif itype == 'INTEGER':    
        block = zeros([blocksize], dtype='Int16')
    elif itype == 'LONG':    
        block = zeros([blocksize], dtype='Int32')
    else:    
        block = zeros([blocksize], dtype='UInt8')
    
    
    #-------------------------------
    # Initialize last block array ?
    #-------------------------------
    if (rem_size != 0):    
        if itype == 'BYTE':    
            last_block = zeros([rem_size], dtype='UInt8')
        elif itype == 'INTEGER':    
            last_block = zeros([rem_size], dtype='Int16')
        elif itype == 'LONG':    
            last_block = zeros([rem_size], dtype='Int32')
        else:    
            last_block = zeros([rem_size], dtype='UInt8')
        
    
    #-----------------------------------
    # Apply conversion, block-by-block
    #-----------------------------------
    for k in arange(1, (n_blocks)+(1)):
        block = fromfile(file_iunit, count=size(block), dtype=str(block.dtype))
        if (SWAP_ENDIAN):
            array(block, copy=0).byteswap(True)
        block = _map[block]    #(block must have integer type)
        if otype == 'BYTE':    
            block = idl_func.byte(block)
        elif otype == 'INTEGER':    
            block = int16(block)
        elif otype == 'LONG':    
            block = int32(block)
        else:    
            dum = int16(0)
        
        if (SWAP_ENDIAN):
            array(block, copy=0).byteswap(True)
        block.tofile(file_ounit)
    #------------------
    # Remainder block
    #------------------
    if (rem_size != 0):    
        last_block = fromfile(file_iunit, count=size(last_block), \
                              dtype=str(last_block.dtype))
        if (SWAP_ENDIAN):
            array(last_block, copy=0).byteswap(True)
        last_block = _map[last_block]
        if otype == 'BYTE':    
            last_block = idl_func.byte(last_block)
        elif otype == 'INTEGER':    
            last_block = int16(last_block)
        elif otype == 'LONG':    
            last_block = int32(last_block)
        else:    
            dum = int16(0)
        
        if (SWAP_ENDIAN):
            array(last_block, copy=0).byteswap(True)
        last_block.tofile(file_ounit)
    
    #------------------
    # Close the files
    #------------------
    print('Finished.')
    print(' ')
    file_iunit.close()
    file_ounit.close()
    
#  Convert_Flow_Grid
#-------------------------------------------------------------------




    
