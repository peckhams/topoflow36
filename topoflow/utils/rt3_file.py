
## Copyright (c) 2009, Scott D. Peckham
## January 7, 2009  (Converted from IDL)

#----------------------------------------------------------------------------

from numpy import *
import numpy

from . import idl_func
from . import tf_utils

#----------------------------------------------------------------------------
class rt3_file:
    def Open(RT3_file, nz=None, n_stacks=None, \
             READ=False, WRITE=False, VERBOSE=False, \
             nx=None, ny=None, RTI_file=None):

        #------------------------------------------------------------
        # NOTE:  nx and ny keywords are used only to return values.
        #------------------------------------------------------------
        n_params = 4 - [nz, n_stacks].count(None)
        unit = None
        _opt = (nz, n_stacks)
        def _ret():
            _optrv = list(zip(_opt, [nz, n_stacks]))
            _rv = [unit]
            _rv += [_o[1] for _o in _optrv if _o[0] is not None]
            return tuple(_rv)
        
        READ = (READ or not(WRITE))
        
        #-------------------------------------
        # Construct RTI filename, if missing
        #-------------------------------------
        if (RTI_file is None):    
            RTI_file = Get_RTI_Filename(RT3_file)
        
        if (VERBOSE):    
            TF_Print('RT3_file = ' + RT3_file)
            TF_Print('Reading RTI_file: ' + RTI_file)
        
        #--------------------------
        # Read info from RTI file
        #--------------------------
        info = Read_RTI_File(RTI_file)
        nx = info.ncols
        ny = info.nrows
        byte_order = info.byte_order
        self.SWAP_ENDIAN = Not_Same_Byte_Order(byte_order)  ########

        #---------------------------------
        # Open RT3 file to read or write
        #---------------------------------
        if (WRITE):    
            file_unit = open(RT3_file, 'wb')
        else:    
            file_unit = open(RT3_file, 'rb')
        
        #-------------------------------------
        # How many (nx x ny x nz) stacks in
        # file assuming data type is FLOAT ?
        #-------------------------------------
        if (READ):    
            temp = idl_func.fstat(file_unit)
            filesize = temp.size
            n_stacks = (filesize / (int32(4) * int32(nx) * int32(ny) * int32(nz)))
        else:    
            filesize = int32(0)
            n_stacks = int32(1)  
        
        return _ret()

    #  Open
    #-------------------------------------------------------------------

    

