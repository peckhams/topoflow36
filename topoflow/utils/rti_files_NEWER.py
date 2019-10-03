
## Copyright (c) 2001-2010, Scott D. Peckham
## January 2009  (converted from IDL)

import glob   # (for exists())
import numpy  # (for things like uint8(), int16(), float64())
import sys    # (for sys.byteorder)

#-------------------------------------------------------------------------
#
#   -------------
#   Functions:
#   -------------
#   unit_test()
#   get_file_prefix()
#   try_to_find_rti_file()
#   get_rti_file_name()
#   get_rti_data_type()
#   get_numpy_data_type()
#   get_rti_byte_order()
#   get_bpe()
#   get_grid_size()
#   byte_swap_needed()

#   exists()
#   read_value()
#   read_info()
#   make_info()
#   write_info()
#
#   make_new_if_needed()  ## future ????????
#   make_new_from_old()   ## future ????????
#
#-------------------------------------------------------------------------  
def unit_test():

    #------------------------
    # Create a new RTI file
    #------------------------
    grid_file = 'TEST_FILE.rts'
    RTI_file  = 'TEST_FILE.rti'
    info = make_info(grid_file=grid_file,
                     ncols=10, nrows=10,
                     xres=101, yres=101, SILENT=False)
    write_info(RTI_file, info, SILENT=False)

    #-----------------------------
    # Read info from an RTI file
    #-----------------------------
    info = read_info(RTI_file, REPORT=True)
    # info = read_info(grid_file, REPORT=True)  # (also works)
    
    #---------------------------
    # Test the "get" functions
    #---------------------------
    print('RTI byte order =', get_rti_byte_order())
    print('RTI file BPE   =', get_bpe( info.data_type ))
    print('RTI grid size  =', get_grid_size( RTI_file ))
    print(' ')
    print('Finished with unit_test().')
    print(' ')
    
#   unit_test()
#-------------------------------------------------------------------
def get_file_prefix( file_name ):

    pos    = file_name.rfind('.')   # (find the last dot)
    prefix = file_name[:pos]

    return prefix

#   get_file_prefix()
#---------------------------------------------------------------------
def get_rti_file_name( file_name ):

    if (file_name[-4:].lower == '.rti'):
        return file_name
    else:
        prefix = get_file_prefix( file_name )
        return (prefix + '.rti')

#   get_rti_file_name()
#---------------------------------------------------------------------
def try_to_find_rti_file( file_name, SILENT=False ):

    #--------------------------------------------------------
    # Notes: If file_name is not an RTI file, try to guess.
    #        Start by replacing the file extension (i.e.
    #        after the dot) with ".rti".  If resulting
    #        RTI file is not found in CWD, try removing
    #        "compound extensions" (everything after "_")
    #        until found or there are no underscores left.
    #--------------------------------------------------------
    RTI_file = get_rti_file_name( file_name )

    if (exists(RTI_file, SILENT=True)):
        return RTI_file
    else:
        prefix = get_file_prefix( file_name )  # (strip extension)
        while (True):
            pos = prefix.rfind('_')
            if (pos == -1): break
            prefix = prefix[:pos]
            RTI_file = (prefix + '.rti')
            if (exists(RTI_file, SILENT=True)):
                return RTI_file

        #-----------------------------------
        # If we didn't return before now,
        # then the RTI file was not found.
        #-----------------------------------
        if not(SILENT):
            print('ERROR: Could not find RTI file for:')
            print('      ' + file_name)
            print(' ')
        #---------------
        return None
    
#   try_to_find_rti_file()
#---------------------------------------------------------------------
def get_rti_data_type( dtype ):

    type_map = {'uint8':'BYTE', 'int16':'INTEGER',
                'int32':'LONG', 'int64':'LONG64',
                'float32':'FLOAT', 'float64':'DOUBLE'}

    return type_map[ dtype.lower() ]

#   get_rti_data_type()
#---------------------------------------------------------------------
def get_numpy_data_type( data_type ):

    type_map = {'BYTE':'uint8', 'INTEGER':'int16',
                'LONG':'int32', 'LONG64':'int64',
                'FLOAT':'float32', 'DOUBLE':'float64'}

    return type_map[ data_type.upper() ]

#   get_numpy_data_type()
#---------------------------------------------------------------------
def get_rti_byte_order( python_byte_order=sys.byteorder ):
    
    order_map = {'big':'MSB', 'little':'LSB'}

    return order_map[ python_byte_order ]

#   get_rti_byte_order()
#---------------------------------------------------------------------
def get_bpe( data_type ):

    #--------------------------------------------
    # Get the number of Bytes Per Element (BPE)
    #--------------------------------------------
    BPE_map = {'BYTE':1, 'INTEGER':2, 'LONG':4,
               'LONG64':8, 'FLOAT':4, 'DOUBLE':8}
    
    return BPE_map[ data_type.upper() ]

#   get_bpe()
#---------------------------------------------------------------------
def get_grid_size( file_name ):

    RTI_file = get_rti_file_name( file_name )

    info = read_info( RTI_file )
    nx   = info.ncols
    ny   = info.nrows
    BPE  = get_bpe( info.data_type )

    return (nx * ny * BPE)

#   get_grid_size()
#---------------------------------------------------------------------
def byte_swap_needed( file_name, info=None ):

    if (info == None):
        RTI_file = get_rti_file_name( file_name )
        info = read_info( RTI_file )
        
    machine_byte_order = get_rti_byte_order()
    
    return (machine_byte_order != info.byte_order)

#   byte_swap_needed()
#-------------------------------------------------------------------------  
def exists(RTI_file, SILENT=False):
 
    #------------------------
    # Does RTI file exist ?
    #------------------------
    result = glob.glob( RTI_file )
    count  = len(result)
    if (count > 0):
        return True
    else:
        if not(SILENT):
            msg = ['ERROR:  RTI file not found. ', ' ', \
                   'The file: ', '  ' + RTI_file, \
                   'was not found in the working directory. ', \
                   ' ']
            for line in msg: print(line)    
        return False

#   exists()
#-------------------------------------------------------------------------  
def read_value(RTI_unit, data_type, UPCASE=False):

    #-----------------------------------------------------------------
    # Notes: A line can just be the newline character, "\n".
    #        "\n".strip() = '' (a null string).
    #        At EOF, readline() returns a null string.
    #        Using "while (True):" with "breaks" is the
    #          "Pythonic way".
    #-----------------------------------------------------------------
    # NB!    Needed to add "strip()" when getting "s" below
    #        to avoid strange bug that only occurred on my
    #        MacBook, running NumPy version 1.0.4.dev4155.
    #        Otherwise kept getting the following error:
    #        "ValueError: setting an array element with a sequence."
    #-----------------------------------------------------------------
    while (True):
        line = RTI_unit.readline()
        EOF  = (line == '')
        if (EOF): break
        char1 = line[0]  # (could be '\n')
        pos   = line.find(':')
        if (char1 != ';') and (pos != -1):
            s = line[pos + 1:].strip()   #### Need strip() !! ####
            break

    #############################
    if (data_type != 'STRING'):
        val = eval(s)
    else:
        val = s
    #############################
    if   (data_type == 'BYTE'):
        ## value = numpy.int16( val )
        value = numpy.uint8(val)
    elif (data_type == 'INTEGER'):    
        value = numpy.int16(val)
    elif (data_type == 'FLOAT'):    
        value = numpy.float32(val)
    elif (data_type == 'DOUBLE'):    
        value = numpy.float64(val)
    elif (data_type == 'LONG'):    
        value = numpy.int32(val)
    elif (data_type == 'LONG64'):    
        value = numpy.int64(val)
    elif (data_type == 'STRING'):    
        value = s.strip()
        if (UPCASE):    
            value = value.upper()
    else:
        raise RuntimeError('no match found for expression')

    ### print '### RTI value =', value
    
    return value

#   read_value()
#---------------------------------------------------------------------
def read_info(file_name, SILENT=False, REPORT=False):

    #----------------------------------------------------------
    #NOTES:  This routine reads information from a RiverTools
    #        Information (RTI) file that describes the grids
    #        in the current directory.

    #        If the value -9999 is used for unknown values in
    #        the RTI file, then it will cause this routine to
    #        generate a warning message in cases where this
    #        value is not reasonable.  However, gmin and gmax
    #        are not checked for validity.

    #        "type" must be INTEGER, FLOAT, LONG, or DOUBLE.
    #        "byte_order" must be MSB or LSB.
    #        Type and byte order are converted to upper case.
    #----------------------------------------------------------
    RTI_file = try_to_find_rti_file( file_name )
    if (RTI_file == None):
        # Message already prints before this.
        # if not(SILENT): print 'ERROR: Unable to find RTI file.'
        return None
 
    #----------------------------
    # Open the RTI file to read
    #----------------------------
    RTI_unit = open(RTI_file, 'r')
    
    #-------------------------------
    # Does first line look right ?
    #-------------------------------
    line1 = RTI_unit.readline()
    line1 = line1.upper().strip()
    if (line1 != 'RIVERTOOLS INFO FILE'):
        msg = ['ERROR:  Invalid RTI file. ', ' ', \
               'The file: ', '  ' + self.RTI_file, \
               'is not a valid RTI file. ', ' ']
        if not(SILENT):
            for line in msg: print(line)
        return None
    
    #--------------------------------
    # Read values from the RTI file
    #--------------------------------
    class info:
        grid_file    = read_value(RTI_unit, 'STRING')
        data_source  = read_value(RTI_unit, 'STRING')
        ncols        = read_value(RTI_unit, 'LONG')
        nrows        = read_value(RTI_unit, 'LONG')
        data_type    = read_value(RTI_unit, 'STRING', UPCASE=True)
        byte_order   = read_value(RTI_unit, 'STRING', UPCASE=True)
        pixel_geom   = read_value(RTI_unit, 'BYTE')
        xres         = read_value(RTI_unit, 'DOUBLE')
        yres         = read_value(RTI_unit, 'DOUBLE')
        zres         = read_value(RTI_unit, 'FLOAT')
        z_units      = read_value(RTI_unit, 'STRING', UPCASE=True)
        y_south_edge = read_value(RTI_unit, 'DOUBLE')
        y_north_edge = read_value(RTI_unit, 'DOUBLE')
        x_east_edge  = read_value(RTI_unit, 'DOUBLE')
        x_west_edge  = read_value(RTI_unit, 'DOUBLE')
        box_units    = read_value(RTI_unit, 'STRING', UPCASE=True)
        gmin         = read_value(RTI_unit, 'FLOAT')
        gmax         = read_value(RTI_unit, 'FLOAT')
        UTM_zone     = read_value(RTI_unit, 'STRING', UPCASE=True)
        
    if (info.UTM_zone == '-1'):    
        info.UTM_zone = 'unknown'

    #---------------------
    # Close the RTI file
    #---------------------
    RTI_unit.close()

    #-----------------------------------------
    # Add some other useful things to "info"
    #-----------------------------------------
    info.n_pixels    = info.ncols * info.nrows
    info.bpe         = get_bpe( info.data_type )
    info.grid_size   = info.bpe * info.n_pixels
    info.SWAP_ENDIAN = byte_swap_needed( RTI_file, info )

    #------------------------------------------------
    # Get grid cell areas, "da", which is either a
    # scalar (if same for all grid cells) or a grid
    #------------------------------------------------
    ## info.da = pixels2.get_da( info )
        
    #------------------
    # Optional report
    #------------------
    if (REPORT):    
        print('-----------------')
        print('Grid Information')
        print('-----------------')
        print('grid_file    =', info.grid_file)
        print('data_source  =', info.data_source)
        print('ncols        =', info.ncols)
        print('nrows        =', info.nrows)
        print('data_type    =', info.data_type)
        print('byte_order   =', info.byte_order)
        print('pixel_geom   =', info.pixel_geom)
        print('xres         =', info.xres)
        print('yres         =', info.yres)
        print('zres         =', info.zres)
        print('z_units      =', info.z_units)
        print('y_south_edge =', info.y_south_edge)
        print('y_north_edge =', info.y_north_edge)
        print('x_east_edge  =', info.x_east_edge)
        print('x_west_edge  =', info.x_west_edge)
        print('box_units    =', info.box_units)
        print('gmin         =', info.gmin)
        print('gmax         =', info.gmax)
        print('UTM_zone     =', info.UTM_zone)
        print(' ')

    return info

#   read_info()
#---------------------------------------------------------------------
def make_info(grid_file=None,
              ncols=None, nrows=None,
              xres=None, yres=None,
              #---------------------------------
              data_source='Unknown',
              data_type='FLOAT',
              byte_order=get_rti_byte_order(),
              pixel_geom=1,
              zres=1.0, z_units='METERS',
              y_south_edge=0.0,
              y_north_edge=None,
              x_west_edge=0.0,
              x_east_edge=None,
              box_units='METERS',
              gmin=-9999.0, gmax=-9999.0,
              UTM_zone='UNKNOWN',
              #----------------------------
              dtype=None, SILENT=True):

    if (grid_file == None):
        print('---------------------------------------------------')
        print('ERROR: make_info() requires "grid_file" argument.')
        print('---------------------------------------------------')
        print(' ')
        return None

    if (ncols == None) or (nrows == None) or \
       (xres == None)  or (yres == None):
        RTI_file = try_to_find_rti_file( grid_file )
        if (RTI_file is not None):
            print('Copying info from RTI_file: ')
            print('    ' + RTI_File)
            info = read_info( grid_file )
            if (ncols == None): ncols = info.ncols
            if (nrows == None): nrows = info.nrows
            if (xres  == None): xres  = info.xres
            if (yres  == None): yres  = info.yres
        else:
            print('----------------------------------------------')
            print('ERROR: "make_info()" requires the arguments:')
            print('        ncols, nrows, xres and yres.')
            print('        Could not find any RTI_file to copy.')
            print('----------------------------------------------')
            print(' ')
            return None

    #-------------------------------------------------
    # Was data type specified with "dtype" keyword ?
    #-------------------------------------------------
    if (dtype is not None):
        data_type = get_rti_data_type( dtype )
        
    #---------------------------------------------------------
    # Compute x_east_edge and y_north_edge from other info ?
    #---------------------------------------------------------
    if (x_east_edge == None):
        if (pixel_geom == 1):
            xsize = (ncols * xres)
        else:
            xres_deg = (xres / numpy.float64(3600))  # (arcseconds to degrees)
            xsize    = (ncols * xres_deg)  
        x_east_edge  = (x_west_edge + xsize)
    #--------------------------------------------------------------------------
    if (y_north_edge == None):
        if (pixel_geom == 1):
            ysize = (nrows * yres)
        else:
            yres_deg = (yres / numpy.float64(3600))  # (arcseconds to degrees)
            ysize    = (nrows * yres_deg)
        y_north_edge = y_south_edge + ysize

    #-----------------------------------------
    # Add some other useful things to "info"
    #-----------------------------------------             
    machine_byte_order = get_rti_byte_order()
    SWAP_ENDIAN        = (machine_byte_order != byte_order)
    n_pixels           = (ncols * nrows)
    bpe                = get_bpe( data_type )
    
    #------------------------------------------
    # Class for putting data into a structure
    #------------------------------------------
    class bunch:
        def __init__(self, **kwds):
            self.__dict__.update(kwds)

    info = bunch( grid_file    = grid_file,
                  data_source  = data_source,
                  ncols        = ncols,
                  nrows        = nrows,
                  data_type    = data_type,
                  byte_order   = byte_order,
                  pixel_geom   = pixel_geom,
                  xres         = xres,
                  yres         = yres,
                  zres         = zres,
                  z_units      = z_units,
                  y_south_edge = y_south_edge,
                  y_north_edge = y_north_edge,
                  x_west_edge  = x_west_edge,
                  x_east_edge  = x_east_edge,
                  box_units    = box_units,
                  gmin         = gmin,
                  gmax         = gmax,
                  UTM_zone     = UTM_zone,
                  #-------------------------------------
                  n_pixels     = n_pixels,
                  bpe          = bpe,
                  grid_size    = bpe * n_pixels,
                  SWAP_ENDIAN  = SWAP_ENDIAN )
    
    if not(SILENT):
        print('Finished creating info structure.')
        print(' ')
        
    return info

#   make_info()
#---------------------------------------------------------------------
def write_info(file_name, info, SILENT=True):

    RTI_file = get_rti_file_name( file_name )
    
    #--------------------------------------------------
    # Prepare floating-point values as output strings
    #--------------------------------------------------
    format = "%22.12f"
    y_south_string = (format % info.y_south_edge).strip()
    y_north_string = (format % info.y_north_edge).strip()
    x_east_string  = (format % info.x_east_edge).strip()
    x_west_string  = (format % info.x_west_edge).strip()
    #------------------------------------------------------
    xres_string    = (format % info.xres).strip()
    yres_string    = (format % info.yres).strip()
    zres_string    = (format % info.zres).strip()
    #------------------------------------------------------    
    pixel_geom_string = str(numpy.int16(info.pixel_geom))
    #------------------------------------------------------ 
    gmin_string    = (format % info.gmin).strip()
    gmax_string    = (format % info.gmax).strip()
    
    #-------------------------
    # Open RTI file to write
    #-------------------------
    RTU = open(RTI_file, 'w')
    
    #-------------------------
    # Write info to RTI file
    #-------------------------
    RTU.write('RiverTools Info File\n')
    RTU.write('\n')
    RTU.write(';--------------------\n')
    RTU.write(';Description of Grid\n')
    RTU.write(';--------------------\n')
    RTU.write('Grid filename:  ' + info.grid_file + '\n')
    RTU.write('Grid source:    ' + info.data_source + '\n')
    RTU.write(' \n')
    RTU.write(';----------------------\n')
    RTU.write(';Grid dimensions, etc.\n')
    RTU.write(';----------------------\n')
    RTU.write('Number of columns:  ' + str(info.ncols) + '\n')
    RTU.write('Number of rows:     ' + str(info.nrows) + '\n')
    RTU.write('Data type:          ' + info.data_type  + '\n')
    RTU.write('Byte order:         ' + info.byte_order + '\n')
    RTU.write('\n')
    RTU.write(';--------------------\n')
    RTU.write(';Pixel geometry info\n')
    RTU.write(';--------------------\n')
    RTU.write('Pixel geometry code:  ' + pixel_geom_string + '\n')
    RTU.write('Pixel x-resolution:   ' + xres_string + '\n')
    RTU.write('Pixel y-resolution:   ' + yres_string + '\n')
    RTU.write('Pixel z-resolution:   ' + zres_string + '\n')
    RTU.write('Z-resolution units:   ' + info.z_units + '\n')
    RTU.write('\n')
    RTU.write(';--------------------------------------------------\n')
    RTU.write(';Bounding box coordinates (degrees lat/lon or UTM)\n')
    RTU.write(';--------------------------------------------------\n')
    RTU.write('South edge y-coordinate:  ' + y_south_string + '\n')
    RTU.write('North edge y-coordinate:  ' + y_north_string + '\n')
    RTU.write('East  edge x-coordinate:  ' + x_east_string  + '\n')
    RTU.write('West  edge x-coordinate:  ' + x_west_string  + '\n')
    RTU.write('Measurement units:        ' + info.box_units + '\n')
    RTU.write('\n')
    RTU.write(';------------------\n')
    RTU.write(';Min and max value\n')
    RTU.write(';------------------\n')
    RTU.write('Minimum value:  ' + gmin_string + '\n')
    RTU.write('Maximum value:  ' + gmax_string + '\n')
    RTU.write('\n')
    RTU.write(';---------\n')
    RTU.write(';UTM zone\n')
    RTU.write(';---------\n')
    RTU.write('UTM zone: ' + info.UTM_zone + '\n')
    RTU.write('\n')
    RTU.write('\n')
    
    #-------------------------
    # Write RTI file trailer
    #-------------------------
    RTU.write(';----------------------------------\n')
    RTU.write(';Notes About RiverTools Info Files\n')
    RTU.write(';----------------------------------\n')
    RTU.write(";(1)  RiverTools grid info filenames end with '.rti'.\n")
    RTU.write(';(2)  The first line should be:  RiverTools Info File\n')
    RTU.write(';(3)  Lines starting with a semi-colon are ignored.\n')
    RTU.write(';(4)  Colons are used to delimit labels from values.\n')
    RTU.write(';(5)  The order of the numbers above is important.\n')
    RTU.write(';(6)  Number of rows and columns are required.\n')
    RTU.write(';(7)  Pixel geometry codes are: 0=Fixed-Angle, ' + '1=Fixed-Length.\n')
    RTU.write(';(8)  Pixel x-resolution and y-resolution are required.\n')
    RTU.write(';(9)  Measurement units must be METERS or DEGREES.\n')
    RTU.write(';(10) Elevation data type is required.\n')
    RTU.write(';     Allowed types are BYTE, INTEGER, LONG, LONG64, FLOAT, DOUBLE.\n')
    RTU.write(';(11) Byte order is required; either LSB or MSB.\n')
    RTU.write(';     (LSB = Least Significant Byte = little-endian)\n')
    RTU.write(';     (MSB = Most Significant Byte  = big-endian)\n')
    RTU.write(";(12) For 'fixed-angle' pixels, bounding lats and lons\n")
    RTU.write(';     are required to compute lengths and areas correctly.\n')
    RTU.write(';(13) Latitudes south of equator are negative and\n')
    RTU.write(';     longitudes west of prime meridian are negative.\n')
    RTU.write(';(14) This file is best modified with the View DEM Info\n')
    RTU.write(';     dialog in the File menu.\n')
    RTU.write('\n')
    RTU.write('\n')
    
    #---------------------
    # Close the RTI file
    #---------------------
    RTU.close()
    
    #-----------------------------
    # Print or display a message
    #-----------------------------
    if not(SILENT):    
        print('Finished writing new RTI file: ')
        print('   ' + RTI_file)
        print(' ')
    
#   write_info()
#---------------------------------------------------------------------

   
