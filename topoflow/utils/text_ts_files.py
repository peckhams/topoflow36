
## Copyright (c) 2010-2014, Scott D. Peckham
## June 2010  (collected here from model_output.py)

#-------------------------------------------------------------------
# This class is for I/O of time series data to multi-column
# text files.  e.g. hydrographs for different grid cells (given
# by IDs) within a DEM.
#-------------------------------------------------------------------

#   unit_test()
#
#   class ts_file():    # (time series as multi-column text)
#
#       open_file()
#       open_new_file()
#       write_header()
#       add_values()
#       add_values_at_IDs()
#       values_at_IDs()
#       close_file()
#       close()
#
#-------------------------------------------------------------------

import numpy as np

from . import file_utils

#-------------------------------------------------------------------
def unit_test(file_name='MC_TEXT_FILE_TEST.txt'):

    print(' ')
    print('Running unit_test()...')

    #---------------------------------------
    # Create a grid of values and some IDs
    #---------------------------------------
    var   = np.arange(25, dtype='Float32').reshape(5,5)
    rows  = np.array([1,2,4])
    cols  = np.array([1,2,4])
    IDs   = (rows, cols)
    n_IDs = np.size(rows)
    var_names = ['Q[1,1]', 'Q[2,2]', 'Q[4,4]']
    n_times = 7
    
    #-----------------------------------
    # Make instance of ts_file() class
    #-----------------------------------
    ts_unit = ts_file()

    OK = ts_unit.open_new_file( file_name, var_names )
    if not(OK):
        print('ERROR during open_new_file().')
        return
   
    for k in range(n_times):
        time = np.float32(k)
        ts_unit.add_values_at_IDs( time, var, IDs )
        var += 1.0
        
    ts_unit.close_file()
    print('Finished writing multi-column text file: ' + file_name)
    print(' ')
    
#   unit_test()
#-------------------------------------------------------------------
class ts_file():

    #-------------------------------------------------------------------   
    def open_file(self, file_name):

        try:
            self.ts_unit = open( file_name, 'r' )
            return True
        except:
            print('SORRY, Could not open multi-column text file:')
            print('   ' + file_name)
            self.mct_unit = False
            return False

    #   open_file()
    #-------------------------------------------------------------------   
    def open_new_file( self, file_name,
                       var_names=['X'],
                       dtype='float64',
                       time_units='minutes' ):

        #-----------------------------------------------------------
        # Note:  The "dtype" argument is included to match similar
        #        routines, but is not currently used. (9/18/14)
        #-----------------------------------------------------------
        
        #----------------------------
        # Does file already exist ?
        #----------------------------
        file_name = file_utils.check_overwrite( file_name )
        self.file_name = file_name
        
        self.var_names  = var_names
        self.time_units = time_units
        
        try:
            self.ts_unit = open( file_name, 'w' )
            self.write_header( var_names )
            return True
        except:
            print('SORRY, Could not open new multi-column text file:')
            print('   ' + file_name)
            self.ts_unit = False
            return False
        
    #   open_new_file()
    #-------------------------------------------------------------------
    def write_header(self, var_names, SILENT=True):

        #-------------------------------------------------------
        # Notes:  This is currently hardwired to have columns
        #         that are each 15 characters wide.
        #-------------------------------------------------------
        n_vars     = len(var_names)
        col_width = 15
        ustr_map  = {'minutes':'min', 'seconds':'sec', \
                     'hours':'hrs', 'days':'days', \
                     'months':'mon', 'years':'yrs'}
        units_str = ustr_map[ self.time_units.lower() ]
        time_str  = 'Time [' + units_str + ']'
        
        #---------------------------
        # Create the header labels
        #---------------------------
        self.ts_unit.write( time_str.rjust(col_width) )
        for k in range(n_vars):
            self.ts_unit.write( var_names[k].rjust(col_width) ) 
        self.ts_unit.write("\n")
        
        #-------------------------------------
        # Print a horizontal line of hyphens
        #-------------------------------------
        width = (n_vars + 1) * col_width
        hline = ''.ljust(width, '-')
        self.ts_unit.write(hline + "\n")

    #   write_header()
    #-------------------------------------------------------------------
    def add_values(self, time, values): 

        col_width = 15        
        tstr = ('%15.7f' % time)
        self.ts_unit.write( tstr.rjust(col_width)  )

        n_values = np.size(values)
        for k in range(n_values):
            vstr = ('%15.7f' % values[k])
            self.ts_unit.write( vstr.rjust(col_width) ) 
        self.ts_unit.write("\n")
        
    #   add_values()
    #-------------------------------------------------------------------
    def add_values_at_IDs(self, time, var, IDs):

        values = self.values_at_IDs( var, IDs )
        self.add_values( time, values )

    #   add_values_at_IDs()
    #-------------------------------------------------------------------
    def values_at_IDs(self, var, IDs):

        #----------------------------------------------------------
        # Notes:  If "var" is a grid, subscript with self.IDs to
        #         get a 1D array of values.  If "var" is scalar,
        #         return a vector with the scalar value repeated
        #         once for each ID in self.IDs.
        #----------------------------------------------------------
        
        #---------------------------------
        # Is variable a grid or scalar ?
        #---------------------------------
        if (np.ndim(var) > 0):
            values = np.float32( var[ IDs ] )
        else:
            #-----------------------------------------------------
            # (3/16/07) Bug fix.  This gets used in case of q0,
            # which is a scalar when INFIL_ALL_SCALARS is true.
            # Without this, don't get a value for every ID.
            #-----------------------------------------------------
            n_IDs  = np.size(IDs[0])
            vector = np.zeros( n_IDs, dtype='Float32')
            values = (vector + np.float32(var)) 

##        print 'VALUES ='
##        print values
        return values
    
    #   values_at_IDs()
    #-------------------------------------------------------------------
    def close_file(self):
        
        self.ts_unit.close()

    #   close_file()
    #-------------------------------------------------------------------
    def close(self):
        
        self.ts_unit.close()

    #   close()
    #-------------------------------------------------------------------
    
