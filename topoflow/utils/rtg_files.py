
# S.D. Peckham
# October 16, 2009

#-------------------------------------------------------------------
# Notes: This could be implemented as a set of functions or
#        a class, as done here.  With a class, variables
#        that describe the grid can be read or computed
#        once and then stored across multiple calls.
#        Note that stand-alone "read_grid()" and "write_grid()"
#        functions are provided outside of the class.
#
#        The class methods (functions) used here have the same
#        names (sometimes through an alias) as the classes
#        that are used for RTS and NCGS files (grid stacks).
#-------------------------------------------------------------------

import numpy as np
from . import bov_files
from . import rti_files

#-------------------------------------------------------------------
#
#   unit_test()   # (function)
#
#   class rtg_file():
#
#       open_file()
#       get_bounds()    # (2020-05-26)
#       check_and_store_info()
#       open_new_file()
#       write_grid()
#           add_grid()  # (alias for write_grid())
#       read_grid()
#           get_grid()  # (alias for read_grid())
#       close_file()
#       close()
#       --------------------
#       byte_swap_needed()
#
#   ------------------------------------
#   These functions are outside of class
#   ------------------------------------
#   read_grid()
#   write_grid()
#   read_xyz_as_grid()
#
#-------------------------------------------------------------------
def unit_test(nx=4, ny=5, VERBOSE=False,
              file_name="TEST_FILE.rtg"):

    print('Running unit_test()...')

    #------------------------------------
    # Make instance of rtg_file() class
    #------------------------------------
    rtg = rtg_file()  
    dx = 100
    dy = 100

    #---------------------------------
    # These are unused for RTG files
    #---------------------------------
##    grid_name  = "depth"
##    long_name  = "depth of water"
##    units_name = "meters"

    info = rti_files.make_info( file_name, nx, ny, dx, dy )
    OK   = rtg.open_new_file( file_name, info, ADD_INDEX=True )

    #----------------------------------------------
    # This also works if RTI file already exists.
    #----------------------------------------------
    # OK = rtg.open_new_file( file_name )
    
    if not(OK):
        print('ERROR during open_new_file().')
        return
    
    grid = np.arange(nx * ny, dtype='Float32')
    grid = grid.reshape( (ny, nx) )
    ### print 'AT START: (nx, ny) =', nx, ny
    
    #---------------------------
    # Write a grid to the file
    #---------------------------
    rtg.write_grid( grid, VERBOSE=True )                             

    #-------------------------------------
    # Re-open the file and read the grid
    #-----------------------------------------------------------
    # The file_name of the file we created is not "file_name",
    # but now includes a time index.  It was saved in rtg, so
    # call rtg.open_file() with no argument to use that name.
    #-----------------------------------------------------------    
    OK = rtg.open_file()
    ## print 'NOW: (nx, ny) =', rtg.nx, rtg.ny
    
    if not(OK): return
    print('rtg.byte_swap_needed() =', rtg.byte_swap_needed())
    print(' ')
    grid = rtg.read_grid( VERBOSE=True )
    print('grid = ')
    print(grid)  

    #-------------------------------------------
    # Write another grid, with next time index
    #-------------------------------------------
    OK = rtg.open_new_file( file_name, info, ADD_INDEX=True)
    rtg.write_grid( grid + 1, VERBOSE=True)
    
#   unit_test()
#-------------------------------------------------------------------
class rtg_file():

    #----------------------------------------------------------
    def __init__(self):
        
        self.time_index  = 1  # (Need here for RTG files.)
        self.index_width = 5
        
    #   __init__()   
    #----------------------------------------------------------
    def open_file(self, file_name=None, UPDATE=False):

        if (file_name == None):
            #-----------------------------------
            # Re-open a previously opened file
            #-----------------------------------
            file_name = self.file_name

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
        self.file_name = file_name

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
                rtg_unit = open(file_name, 'rb+')
                self.rtg_unit = rtg_unit
            else:
                rtg_unit = open(file_name, 'rb')
                self.rtg_unit = rtg_unit
            ### return rtg_unit
            return True
        except:
            print('ERROR during rtg.open_file().')
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
        self.format = 'RTG'
        self.file_name = file_name
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
    def add_index_to_file_name(self, file_name): 

        #---------------------------------------------------
        # Insert a 5-digit time index into new file name ?
        #---------------------------------------------------
        n_digits = self.index_width
        suffix = '_' + str(self.time_index).zfill(n_digits)
        pos = file_name.rfind('.')
        if (pos != -1):
            new_file_name  = file_name[:pos] + suffix 
            new_file_name += file_name[pos:]  # (add extension)
        else:
            new_file_name = file_name + suffix

        return new_file_name

    #   add_index_to_file_name() 
    #----------------------------------------------------------
    def open_new_file(self, file_name, info=None,
                      var_name='UNKNOWN', dtype='float32',
                      VERBOSE=False, ADD_INDEX=False,
                      MAKE_RTI=True, MAKE_BOV=False):

        #---------------------------------------------------
        # Insert a 5-digit time index into new file name ?
        #---------------------------------------------------
        if (ADD_INDEX):
            grid_file = self.add_index_to_file_name( file_name )
        else:
            grid_file = file_name

        #---------------------------------------
        # Check and store the grid information
        #---------------------------------------
        self.check_and_store_info( grid_file, info, var_name,
                                   dtype, MAKE_RTI, MAKE_BOV )
              
        #------------------------------------
        # Try to open new RTG file to write
        #------------------------------------
        try:  
            if (VERBOSE):
                print('Preparing to write new RTG file:')
                print('   ' + grid_file)  
            rtg_unit = open(grid_file, 'wb')
            self.rtg_unit = rtg_unit
            return True
        except:
            return False
        
    #   open_new_file()
    #----------------------------------------------------------
    def write_grid(self, grid, VERBOSE=False):

        if (VERBOSE):    
            print('Writing grid values...')
            
        #--------------------------------
        # Swap byte order, if necessary
        #--------------------------------------------
        # The goal is to use the same byte order as
        # the other binary files in the data set.
        # Byte order is recorded in the RTI file.
        #--------------------------------------------
        if (self.info.SWAP_ENDIAN): grid.byteswap(True)

        #-------------------------------
        # Write grid as binary to file
        #-------------------------------
        grid.tofile( self.rtg_unit )

        #----------------
        # Close the file
        #----------------
        self.rtg_unit.close()
        
        self.time_index += 1   ##############
        
        if (VERBOSE):     
            print('Finished writing grid to:')
            print('  ' + self.file_name)
        
    #   write_grid()
    #----------------------------------------------------------
    def add_grid(self, grid):

        self.write_grid( grid )

    #   add_grid()
    #----------------------------------------------------------
    def read_grid(self, dtype='float32', rtg_type=None,
                  REPORT=False, VERBOSE=False):

        if (VERBOSE):    
            print('Reading grid values....')

        #--------------------------------------
        # Different ways to specify data type
        #--------------------------------------
        if (rtg_type is not None):
            dtype = rti_files.get_numpy_data_type( rtg_type )

        #--------------------------------
        # Read grid as binary from file
        #--------------------------------
        n_values = self.nx * self.ny
        grid = np.fromfile( self.rtg_unit, count=n_values,
                               dtype=dtype )
        grid = grid.reshape( self.ny, self.nx )

        #----------------
        # Close the file
        #----------------
        self.rtg_unit.close()
        
        #--------------------------------
        # Swap byte order, if necessary
        #--------------------------------
        if (self.info.SWAP_ENDIAN): grid.byteswap( True )

        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            gmin = grid.min()
            gmax = grid.max()
            print('    min(grid) =', gmin)
            print('    max(grid) =', gmax)
        if (VERBOSE):    
            print('Finished reading grid from:')
            print('  ' + self.file_name)
        
        return grid
    
    #   read_grid()
    #----------------------------------------------------------
    def get_grid(self, grid):

        self.read_grid( grid )

    #   get_grid()
    #----------------------------------------------------------
    def close_file(self):

        self.rtg_unit.close()

    #   close_file()
    #----------------------------------------------------------
    def close(self):

        self.rtg_unit.close()

    #   close()
    #----------------------------------------------------------
    def byte_swap_needed(self):

        machine_byte_order = rti_files.get_rti_byte_order()
        SWAP =  (machine_byte_order != self.info.byte_order)
        
        return SWAP
    
    #   byte_swap_needed()
    #----------------------------------------------------------
#--------------------------------------------------------------------
def read_grid( RTG_file, rti, RTG_type='FLOAT',
               REPORT=False, SILENT=True ):

    #----------------------------------------------------------
    # Notes: This function is outside of the "rtg_file" class
    #        and will be more convenient in many cases.
    #----------------------------------------------------------
    if not(SILENT):    
        print('Reading grid values...')
        
    #-------------------
    # Read in the grid
    #-------------------
    file_unit = open( RTG_file, 'rb' )
    RTG_type = rti.data_type   ######### (07/12/18) ########
    dtype = rti_files.get_numpy_data_type( RTG_type )
    grid  = np.fromfile( file_unit, count=rti.n_pixels,
                            dtype=dtype )
   
#     print 'n_pixels = ', rti.n_pixels                         
#     print 'nrows = ', rti.nrows
#     print 'ncols = ', rti.ncols
#     print 'grid.shape =', grid.shape
    
    grid  = grid.reshape( rti.nrows, rti.ncols )
    
    if (rti.SWAP_ENDIAN):
        grid.byteswap( True )
    file_unit.close()
    
    if not(SILENT):    
        print('Finished reading grid from:')
        print('  ' + RTG_file)
    
    #------------------
    # Optional report
    #------------------
    if (REPORT):    
        print('    min(grid) =', grid.min())
        print('    max(grid) =', grid.max())
        print(' ')

    return grid

#   read_grid()
#--------------------------------------------------------------------
def write_grid( grid, RTG_file, rti, RTG_type='FLOAT',
                REPORT=False, SILENT=True ):

    if not(SILENT):    
        print('Writing grid values...')
        
    #---------------------------------
    # Convert to specified data type
    #---------------------------------
    if   (RTG_type == 'BYTE'):    
        grid = np.uint8(np.int16(grid))
    elif (RTG_type == 'INTEGER'):    
        grid = np.int16(grid)
    elif (RTG_type == 'LONG'):    
        grid = np.int32(grid)
    elif (RTG_type == 'LONG64'):    
        grid = np.int64(grid)
    elif (RTG_type == 'FLOAT'):    
        grid = np.float32(grid)
    elif (RTG_type == 'DOUBLE'):    
        grid = np.float64(grid)
    else:
        raise RuntimeError('No match found for expression')
    
    #--------------------------------------
    # Note that we use same byte order as
    # the other binary files in data set
    #--------------------------------------
    file_unit = open( RTG_file, 'wb' )
    if (rti.SWAP_ENDIAN):
        grid.byteswap(True)
    grid.tofile( file_unit )
    file_unit.close()
    
    if not(SILENT):     
        print('Finished writing grid to:')
        print('    ' + RTG_file)
    
    #------------------
    # Optional report
    #------------------
    if (REPORT):
        print('    min(grid) =', grid.min())
        print('    max(grid) =', grid.max())
        print(' ')
    
#   write_grid()
#--------------------------------------------------------------------
def read_xyz_as_grid(XYZ_file, RTG_file):

    import sys
    from . import tf_utils
    
    #--------------------------------------------
    # Notes:  It is assumed that X and Y values
    # are for the center of the pixel.
    #--------------------------------------------
    n_lines = tf_utils.Count_Lines(XYZ_file)
    
    #--------------------------------
    # Read XYZ data into a 2D array
    #--------------------------------
    data = np.zeros([n_lines, 3], dtype='Float32')
    xyz  = np.zeros([3], dtype='Float32')
    k    = np.int32(0)
    file_unit = open( XYZ_file, 'r' )
    while not(idl_func.eof(file_unit)):
        xyz = idl_func.readf(file_unit, xyz)
        data[k,:] = xyz
        k = (k + int32(1))
    file_unit.close()
    data = data[:k, :]   #(discard extra slots)
    
    #--------------------------------
    # Extract x, y & z as 1D arrays
    #-----------------------------------
    # xres,xmin,xmax will have DOUBLE
    # type, and so will yres,ymin,ymax
    #-----------------------------------
    x = np.float64(data[:,0])
    y = np.float64(data[:,1])
    z = np.float32(data[:,2])
    
    #-----------------------------
    # Compute the mins and maxes
    #-----------------------------
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()
    zmin = z.min()
    zmax = z.max()
    
    #-----------------------------
    # How many cols are there ?
    # Assume row major ordering.
    # Compute number of rows.
    #-----------------------------
    w = np.where(y == y[0])
    ncols = np.size(w[0])
    n     = np.size(x)
    nrows = (n / ncols)
    
    #-------------------------
    # Find xres, yres & zres
    #-------------------------
    yres = np.absolute(y[ncols] - y[0])
    xres = np.absolute(x[1] - x[0])
    zres = np.float32(1.0)
    
    #--------------------
    # This isn't needed
    #--------------------
    #x = np.reform(x, ncols, nrows)
    #y = np.reform(y, ncols, nrows)
    #z = np.reform(z, ncols, nrows)

    #----------------------------------
    # Get byte order of user's computer
    #----------------------------------
    big_endian = (sys.byteorder == 'big')
    if (big_endian):    
        byte_order = 'MSB'
    else:    
        byte_order = 'LSB'
        
    #----------------------------------
    # Write data as a binary RTG file
    #----------------------------------
    file_unit = open( RTG_file, 'wb' )
    SWAP_ENDIAN = tf_utils.Not_Same_Byte_Order(byte_order)
    if (SWAP_ENDIAN):
        np.array(z, copy=0).byteswap(True)
    z.tofile(file_unit)
    file_unit.close()
    
    #---------------------------
    # Get name of new RTI file
    #---------------------------
    RTI_file = Get_RTI_Filename(RTG_file)
    
    #---------------------------------
    # Create an RTI file for new DEM
    #---------------------------------
    info = idl_func.bunch(DEM_file=RTG_file, RTI_file=RTI_file, \
                          data_source='Created by XYZ to Grid routine.', \
                          ncols=ncols, nrows=nrows, \
                          data_type='FLOAT', byte_order=byte_order, \
                          pixel_geom=uint8(1), \
                          xres=xres, yres=yres, zres=zres, \
                          z_units='METERS', \
                          y_north_edge=ymax + (yres / float64(2)), \
                          y_south_edge=ymin - (yres / float64(2)), \
                          x_west_edge=xmin - (xres / float64(2)), \
                          x_east_edge=xmax + (xres / float64(2)), \
                          box_units='METERS', emin=zmin, emax=zmax, \
                          UTM_zone='UNKNOWN')
    
    #-------------------------
    # Write the new RTI file
    #-------------------------
    tf_utils.Write_RTI_File(RTI_file, info, SILENT=True)
    
    #----------------------
    # Print final message
    #----------------------
    print('Finished creating RTG file.')
    print(' ')
    
#   read_xyz_as_grid()
#--------------------------------------------------------------------



    
