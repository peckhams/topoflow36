
# S.D. Peckham
# October 14, 2009

import os
import os.path
import sys

import numpy as np

from . import bov_files
from . import file_utils
from . import rti_files

#-------------------------------------------------------------------
#
#   unit_test()
#
#   class rts_file():
#
#       open_file()
#       get_bounds()    # (2020-05-26)
#       open_new_file()
#       add_grid()
#       get_grid()
#       read_grid()  # alias to get_grid()
#       close_file()
#       close()
#       --------------------
#       byte_swap_needed()
#       number_of_grids()
#
#-------------------------------------------------------------------
def unit_test(nx=4, ny=5, n_grids=6, VERBOSE=False,
              file_name="TEST_FILE.rts"):

    print('Running unit_test()...')

    #------------------------------------
    # Make instance of rts_file() class
    #------------------------------------
    rts = rts_file()  
    dx = 100
    dy = 100

    #---------------------------------
    # These are unused for RTS files
    #---------------------------------
##    grid_name  = "depth"
##    long_name  = "depth of water"
##    units_name = "meters"

    info = rti_files.make_info( file_name, nx, ny, dx, dy )
    OK = rts.open_new_file( file_name, info )

    if not(OK):
        print('ERROR during open_new_file().')
        return
    
    grid = np.arange(nx * ny, dtype='Float32')
    grid = grid.reshape( (ny, nx) )
    
    #----------------------------------
    # Add some test grids to the file
    #----------------------------------
    for time_index in range(n_grids):
        rts.add_grid( grid )                         
        grid = (grid + 1)
        
    rts.close_file()
    print('Finished writing file: ' + file_name)
    print(' ')

    #---------------------------------------------
    # Re-open the file and read grids one-by-one 
    #---------------------------------------------
    OK = rts.open_file( file_name )
    if not(OK): return
    n_grids = rts.number_of_grids()
    print('Reading grids from RTS file: ')
    print('rts.number_of_grids()  =', n_grids)
    print('rts.byte_swap_needed() =', rts.byte_swap_needed())
    print(' ')
    for time_index in range(n_grids):
        grid = rts.get_grid( time_index )
        print('grid[' + str(time_index) + '] = ')
        print(grid)
        print('-----------------------------------------------')

    #----------------------------
    # Go back and read 2nd grid
    #----------------------------
    grid = rts.get_grid( 1 )
    print(' ')
    print('Reading second grid again...')
    print('Second grid =')
    print(grid)
    print('-----------------------------------------------')
    rts.close_file()
    print('Finished reading file: ' + file_name)
    print(' ')

    #---------------------------------------
    # Re-open the file and change one grid
    #---------------------------------------
    print('Updating RTS file:', file_name)
    grid = np.ones( (ny, nx), dtype='Float32' )
    OK = rts.open_file( file_name, UPDATE=True )
    if not(OK): return
    rts.add_grid( grid, time_index=0 )
    rts.close_file()
    print('Finished updating RTS file.')
    print(' ')
    
    #---------------------------------------------
    # Re-open the file and read grids one-by-one 
    #---------------------------------------------
    OK = rts.open_file( file_name )
    if not(OK): return
    n_grids = rts.number_of_grids()
    print('Reading grids from RTS file: ')
    print('rts.number_of_grids()  =', n_grids)
    print('rts.byte_swap_needed() =', rts.byte_swap_needed())
    print(' ')
    for time_index in range(n_grids):
        grid = rts.get_grid( time_index )
        print('grid[' + str(time_index) + '] = ')
        print(grid)
        print('-----------------------------------------------')
    rts.close_file()
    print('Finished reading file: ' + file_name)
    print(' ')    
    
#   unit_test()
#-------------------------------------------------------------------
class rts_file():

    #----------------------------------------------------------
    def open_file(self, file_name, UPDATE=False):

        info = rti_files.read_info( file_name )
        if (info == -1): return

        #----------------------
        # Store info in state
        #----------------------
        self.info = info
        self.nx   = info.ncols
        self.ny   = info.nrows
        self.dx   = info.xres
        self.dy   = info.yres

        BPE = rti_files.get_bpe( info.data_type )
        self.grid_size   = (self.nx * self.ny * BPE)
        self.SWAP_ENDIAN = self.byte_swap_needed()
        self.file_name   = file_name
        self.time_index  = 0

        #------------------------------------------
        # Compute map bounds in different "styles"
        #------------------------------------------
        minlon = info.x_west_edge
        maxlon = info.x_east_edge
        minlat = info.y_south_edge
        maxlat = info.y_north_edge
        #----------------------------------
        # For matplotlib.pyplot.imshow():
        self.bounds = [minlon, maxlon, minlat, maxlat]                
        #-----------------------------------------------------
        # For "ipyleaflet":
        # self.bounds = [[minlat, maxlat], [minlon, maxlon]]
        #-----------------------------------------------------
        # For "sw_and_ne_corner":    
        # self.bounds = [minlon, minlat, maxlon, maxlat]
        #-------------------------------------------------
                
        #-----------------------------------
        # Open file to read only or update
        #-----------------------------------        
        try:
            if (UPDATE):
                rts_unit = open(file_name, 'rb+')
                self.rts_unit = rts_unit
            else:
                rts_unit = open(file_name, 'rb')
                self.rts_unit = rts_unit
            ### return rts_unit
            return True
        except:
            print('ERROR during rts.open_file().')
            return False
    
    #   open_file()
    #----------------------------------------------------------
    def get_bounds(self):

        # Note:  Saved into self by open_file().
        if (hasattr(self, 'bounds')):
            return self.bounds
        else:
            return None
            
    #   get_bounds()    
    #----------------------------------------------------------
    def check_and_store_info(self, file_name, info=None,
                             var_name='UNKNOWN',
                             dtype='float32',
                             MAKE_RTI=True, MAKE_BOV=False):

        #-----------------------------------------------------
        # Note: This object (self) may be new or it may have
        #       been used previously.  In the latter case,
        #       "info" should still be available in "self".
        #       We only need info if MAKE_RTI or MAKE_BOV.
        #-----------------------------------------------------
        self.format     = 'RTS'
        self.file_name  = file_name
        self.dtype      = dtype
        self.time_index = 0  # (need here for RTS files)
        if not(MAKE_RTI or MAKE_BOV): return

        #---------------------------------
        # Was "info" argument provided ?
        #---------------------------------
        NEW_INFO = True
        if (info == None):
            try:
                info = self.info
                NEW_INFO = False
                ## print 'Found info in state.'
            except:
                #------------------------------------------
                # Try to find RTI file to copy info from.
                # Don't create a new RTI file.
                #------------------------------------------
                RTI_file = rti_files.try_to_find_rti_file( file_name )
                if (RTI_file != 'none'):
                    info = rti_files.read_info( RTI_file )
                    ## print 'Reading info from: ' + RTI_file
                else:
                    print('ERROR during open_new_file():')
                    print('   Could not find RTI file and "info"')
                    print('   argument was not provided.')
                    print(' ')
                    return

        #-----------------------------
        # Update "info" as necessary
        #-----------------------------
        info.grid_file   = file_name
        info.data_type   = rti_files.get_rti_data_type( dtype )
        info.data_source = 'TopoFlow 3.0'
        info.gmin        = -9999.0
        info.gmax        = -9999.0
        
        #---------------------------------------
        # If new "info" was provided, store it
        #---------------------------------------
        if (NEW_INFO):
            self.info = info
            self.nx   = info.ncols
            self.ny   = info.nrows           
            ## print 'Stored new info in state.'
            
##        #---------------------------------
##        # Was "info" argument provided ?
##        #---------------------------------
##        if (info is not None):
##            #------------------------------
##            # Save info to a new RTI file
##            #------------------------------
##            prefix   = rti_files.get_file_prefix( file_name )
##            RTI_file = (prefix + '.rti')
##            rti_files.write_info( RTI_file, info )          
##
##        else:
##            #------------------------------------------
##            # Try to find RTI file to copy info from.
##            # Don't create a new RTI file.
##            #------------------------------------------
##            RTI_file = rti_files.try_to_find_rti_file( file_name )
##            if (RTI_file != 'none'):
##                info = rti_files.read_info( RTI_file )
##                info.file_name = file_name
##                info.data_type = rti_files.get_rti_data_type( dtype )
##            else:
##                print 'ERROR during open_new_file():'
##                print '   Could not find RTI file and "info"'
##                print '   argument was not provided.'
##                print ' '
##                return

        #-------------------
        # Write RTI file ?
        #-------------------
        if (MAKE_RTI):
            prefix   = rti_files.get_file_prefix( file_name )
            RTI_file = (prefix + '.rti')
            rti_files.write_info( RTI_file, info )
            # print 'Wrote grid info to: ' + RTI_file   ######
            
        #-------------------
        # Write BOV file ?
        #-------------------
        if (MAKE_BOV):
            bov_files.write_info_as_bov( file_name, info, var_name)
                                         ###  time )
        
    #   check_and_store_info()
    #----------------------------------------------------------
    def open_new_file(self, file_name, info=None,
                      var_name='UNKNOWN',
                      dtype='float32',
                      VERBOSE=False,
                      MAKE_RTI=True, MAKE_BOV=False):

        #----------------------------
        # Does file already exist ?
        #----------------------------
        file_name = file_utils.check_overwrite( file_name )
        self.file_name = file_name
        
        #---------------------------------------
        # Check and store the grid information
        #---------------------------------------
        self.check_and_store_info( file_name, info, var_name,
                                   dtype, MAKE_RTI, MAKE_BOV )
        
        #------------------------------------
        # Try to open new RTS file to write
        #------------------------------------
        try:            
            if (VERBOSE):
                print('Preparing to write new RTS file:')
                print('   ' + file_name)  
            self.rts_unit = open(file_name, 'wb')
            return True
        except:
            return False
        
    #   open_new_file()
    #----------------------------------------------------------
    def add_grid(self, grid, time_index=-1):

        #------------------------------------------------------
        # Notes:  If the "grid" argument is actually a scalar
        #         then we still make sure to write a grid.
        #
        #         In order to maintain a consistent byte
        #         order within a set of RTG and RTS files,
        #         the grid may be byte-swapped before it is
        #         written to the file.  However, this must
        #         be done carefully, without setting the
        #         "inplace" argument to the byteswap method
        #         to True on the original grid.  Otherwise,
        #         the byte order of the original grid will
        #         "flip-flop" every time this function is
        #         called.
        #------------------------------------------------------
        dtype = self.dtype   # (set in open_new_file())
        
        #---------------------------------------------
        # Can use time_index to move file pointer
        # and overwrite an existing grid vs. append.
        #---------------------------------------------
        if (time_index >= 0):
            offset = (time_index * self.grid_size)
            self.rts_unit.seek( offset )

        #-------------------------------------------------
        # Convert grid to Float32 (usually from Float64)
        #-------------------------------------------------
        out_grid = grid.copy().astype(dtype)

        #-------------------------------------------------
        # Convert byteorder, if needed to match the byte
        # order of other RTG and RTS files in data set.
        #-------------------------------------------------
        if (self.info.SWAP_ENDIAN):
            inplace = True
            out_grid.byteswap(inplace)

        #---------------------------------------            
        # Convert "grid" from scalar to grid ?
        #---------------------------------------
        if (np.ndim(out_grid) == 0):
            out_grid += np.zeros([self.ny, self.nx], dtype=dtype)

        #--------------------------------------------
        # Write grid as binary to existing RTS file
        #--------------------------------------------            
        out_grid.tofile( self.rts_unit )

        #-------------------------
        # Advance the time_index
        #-------------------------
        self.time_index += 1
        
##        #--------------------------------------------
##        # Write grid as binary to existing RTS file
##        #--------------------------------------------
##        if (np.ndim(grid) == 0):
##            #-----------------------------------------------
##            # "grid" is actually a scalar (dynamic typing)
##            # so convert it to a grid before saving
##            #-----------------------------------------------
##            grid2 = grid + np.zeros([self.ny, self.nx], dtype='Float32')
##            if (self.info.SWAP_ENDIAN):
##                grid2.byteswap().tofile( self.rts_unit )
##            else:
##                grid2.tofile( self.rts_unit )
##        else:
##            if (self.info.SWAP_ENDIAN):
##                grid.byteswap().tofile( self.rts_unit )
##            else:
##                grid.tofile( self.rts_unit )
##        self.time_index += 1
        
        #-------------------------------
        # Write grid as binary to file
        #-------------------------------
##        grid.tofile( self.rts_unit )
##        self.time_index += 1
        
    #   add_grid()
    #----------------------------------------------------------
    def get_grid(self, time_index, dtype='float32'):

        #-----------------------------------------------
        # Compute offset from time_index and grid_size
        #-----------------------------------------------
        n_values = self.nx * self.ny
        offset   = (time_index * self.grid_size)
        self.rts_unit.seek( offset )

        grid = np.fromfile( self.rts_unit, count=n_values,
                               dtype=dtype )
        grid = grid.reshape( self.ny, self.nx )

        #--------------------------------
        # Swap byte order, if necessary
        #--------------------------------
        if (self.info.SWAP_ENDIAN):
            inplace = True
            grid.byteswap(inplace)
            
        return grid
    
    #   get_grid()
    #-------------------------------------------------------------------
    def read_grid(self, time_index, dtype='float32'):
    
        # Note:  This is an alias to get_grid().
        grid = self.get_grid( time_index, dtype=dtype )
        return grid
        
    #   read_grid()
    #-------------------------------------------------------------------
    def close_file(self):

        self.rts_unit.close()

    #   close_file()
    #-------------------------------------------------------------------
    def close(self):

        self.rts_unit.close()

    #   close()
    #-------------------------------------------------------------------
    def byte_swap_needed(self):

        machine_byte_order = rti_files.get_rti_byte_order()
        SWAP =  (machine_byte_order != self.info.byte_order)
        
        return SWAP
    
    #   byte_swap_needed()
    #-------------------------------------------------------------------
    def number_of_grids(self):

        file_size = os.path.getsize( self.file_name )
        n_grids   = (file_size / self.grid_size)

        # self.file_size = file_size
        # self.n_grids   = n_grids
        
        return n_grids
        
    #   number_of_grids()
    #-------------------------------------------------------------------

    
