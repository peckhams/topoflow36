
# Copyright (c) 2022, Scott D. Peckham

# Dec 2022.  For working with NextGen (and other) CSV files

#---------------------------------------------------------------------
# This class is for reading data from multi-column text (mct) files.
# This includes CSV files, if delimeter is ",".
#---------------------------------------------------------------------
#
#   unit_test()
#
#   class mct_file():
#
#       open_file()
#       count_lines()
#       read_column()
#       find_wettest_month()
#       close_file()
#
#-------------------------------------------------------------------

import numpy as np

#-------------------------------------------------------------------
def unit_test(file_name='MCT_FILE_TEST.txt'):

    print('Running unit_test() for mct_files.py...')

    f = mct_file()
    f.count_lines( file_name )
    f.open_file( file_name )
    vals = f.read_column(col=1, skip=1, delim=',')    
    f.close_file()
    
    print('Finished reading multi-column text file: ' + file_name)
    print(' ')
    
#   unit_test()
#-------------------------------------------------------------------
class mct_file():

    #-------------------------------------------------------------------   
    def open_file(self, mct_file):

        try:
            self.mct_unit = open( mct_file, 'r' )
        except:
            print('SORRY, Could not open multi-column text file:')
            print('   ' + file_name)
            self.mct_unit = False

    #   open_file()
    #-------------------------------------------------------------------
    def count_lines(self, mct_file, SILENT=False ):

        #----------------
        # Open the file
        #----------------
        mct_unit = open(mct_file, 'r')
    
        #------------------
        # Count the lines
        #------------------
        n_lines = np.int32(0)
        n_total = np.int32(0)
        while (True):
            line = mct_unit.readline()
            if (line == ''):
                break
            n_total += 1
            #---------------------------
            # Count the nonblank lines
            #---------------------------
            n_chars = len(line.strip())
            if (n_chars != 0):    
                n_lines += 1
    
        #-----------------
        # Close the file
        #-----------------
        mct_unit.close()
    
        #--------------------
        # Print a message ?
        #--------------------
        if not(SILENT):    
            print('For the file: ' + mct_file)
            print('Total number of lines    = ' + str(n_total))
            print('Number of nonblank lines = ' + str(n_lines))
            print(' ')
        
        self.n_lines = n_total
        return n_total

    #   count_lines()
    #-------------------------------------------------------------------   
    def read_column(self, mct_file, col=0, skip_lines=0, delim=' ',
                    dtype=None, DATETIME=False, REPORT=True ):

        self.open_file( mct_file )

        #---------------------
        # Skip header rows ?
        #---------------------
        for k in range(skip_lines):
            line = self.mct_unit.readline()        
    
        k = 0
        if (DATETIME):
            dtype = '<U19'   # (string array <= 19 chars each)
        if (dtype is None):
            dtype = 'float64'

        values = np.zeros( self.n_lines, dtype=dtype )            
            
        while True:
            line  = self.mct_unit.readline()
            if (line == ''):
                break
            words = line.split( delim )
            value = np.array( words[col], dtype=dtype )
            values[k] = value
            k += 1

        values = values[:k]
        if (REPORT):
            vmin = values.min()
            vmax = values.max()
            print('vmin, vmax =', vmin, ',', vmax)
            print()
        self.close_file()

        return values
        
    #   read_column()
    #-------------------------------------------------------------------
    def find_wettest_month(self, mct_file, max_years=20, REPORT=True):    
    
        self.count_lines( mct_file )
        datetimes = self.read_column( mct_file, col=0, skip_lines=1,
                         delim=',', DATETIME=True, REPORT=False )
        rainrates = self.read_column( mct_file, col=1, skip_lines=1,
                         delim=',', REPORT=False )
        n_lines = len( datetimes )
        last_month  = None
        rain_sums   = np.zeros((20, 12), dtype='float64')  # (year, month)
        rain_sum    = 0
        year_index  = 0
        month_index = 0
        
        for k in range(n_lines):
            datetime = datetimes[k]
            rainrate = rainrates[k]
            date = datetime[0:10]
            (year, month, day) = date.split('-')
            if (last_month is None):
                last_month = month
                last_year  = year
                min_year   = np.int16(year)
            DONE      = (k == n_lines-1)
            NEW_MONTH = (month != last_month) or DONE
            NEW_YEAR  = (year != last_year) or DONE
            if (NEW_MONTH):
                rain_sums[year_index, month_index ] = rain_sum
                if (REPORT):
                    print('Done with month =', last_month)
                    print('      rain sum  =', rain_sum)
                    print('          index =', month_index)
                    print('----------------------------------')
                month_index = (month_index + 1) % 12
                last_month = month  ####
                rain_sum = 0
            if (NEW_YEAR):
                if (REPORT):
                    print('## Done with year  =', last_year)
                    ## print('          index =', year_index)
                year_index += 1
                last_year = year   ####
            rain_sum += rainrate

        #------------------------------
        # Find wettest year and month
        #------------------------------
        rmax = rain_sums.max()
        w1 = np.where( rain_sums == rmax )
        year_index  = w1[0][0]
        month_index = w1[1][0]
        print()
        print('Wettest year, month =', year_index + min_year, ',', month_index + 1 )
        print('  rain sum =', rmax)
        print()

        #-----------------------------
        # Find driest year and month
        #-----------------------------
#         rmin = rain_sum.min()
#         w2 = np.where( rain_sum == rmin )  ## HANDLE ZEROS
                                
    #   find_wettest_month()
    #-------------------------------------------------------------------
    def close_file(self):
        
        self.mct_unit.close()

    #   close_file()
    #-------------------------------------------------------------------

    
