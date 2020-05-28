
# Copyright (c) 2014, Scott D. Peckham
# September 2014

#------------------------------------------------------------------------
# Note: Wrote this on 9/19/14 to eliminate a cyclic dependency between
#       BMI_base.py and basins.py.  basins.py inherits from BMI_base.py
#       and BMI_base.py has initialize_basin_vars.
#------------------------------------------------------------------------
# 
# write_outlet_file()
# read_outlet_file()
# check_outlet_IDs()
# read_main_basin_IDs()
# 
#-----------------------------------------------------------------------

import numpy as np
import os
import os.path

#-----------------------------------------------------------------------
def write_outlet_file( self ):

    #---------------------------------------------------
    # Written to provide default if missing. (11/7/11)
    #---------------------------------------------------
    dash_line = ''.rjust(60, "-")  #########
    unit = open( self.outlet_file, 'w' )
    unit.write('\n')
    unit.write( dash_line + '\n')
    unit.write(' Monitored Grid Cell (Outlet) Information\n')
    unit.write( dash_line + '\n')
    format1 = '%10s%10s%16s%16s'
    format2 = '%10i%10i%16.4f%16.4f'
    header = format1 % ('Column', 'Row', 'Area [km^2]', 'Relief [m]')
    unit.write( header + '\n' )
    unit.write( dash_line + '\n')
    info_str = format2 % (self.nx/2, self.ny/2, 0.0, 0.0)
    unit.write( info_str + '\n')
    unit.close()

#   write_outlet_file()
#-------------------------------------------------------------------
def read_outlet_file( self ):

    #----------------------------------------------------------
    # Note: The self argument is a TopoFlow BMI object.
    #----------------------------------------------------------
    # Note: outlet_IDs and basin_IDs are 1D arrays or vectors
    #       but have very different sizes.  Recall that they
    #       are (row,col) tuples and not just long ints.
    #----------------------------------------------------------
    # Note:  "basin_IDs" gives IDs of all grid cells that lie
    #        within the basin that drains to outlet_ID.
    #----------------------------------------------------------
    # Notes: outlet data is stored in a multi-column
    #        textfile as:
    #            Column, Row, Area [km^2], Relief [m]
    #------------------------------------------------------
    # Make sure outlet_file is a full path name (9/21/14)
    #------------------------------------------------------
    # print('### self.case_prefix   =', self.case_prefix )
    # print('### self.cfg_directory =', self.cfg_directory )
    # print('### self.in_directory  =', self.in_directory )
    
    if (hasattr(self, 'pixel_file')):
        outlet_file = self.pixel_file    # (10/25/11)
    else:
        outlet_file = (self.case_prefix + '_outlets.txt')
    # 2020-05-03. Now using cfg_dir vs. in_dir.
    self.outlet_file = (self.cfg_directory + outlet_file)
    ## self.outlet_file = (self.in_directory + outlet_file)
    ## print '### self.outlet_file =', self.outlet_file

    #--------------------------------------
    # Does outlet_file exist ? (10/25/11)
    #--------------------------------------
    if not(os.path.exists( self.outlet_file )):
        hash_line = ''.rjust(70, "#")
        print(hash_line) 
        print(' ERROR: Could not find monitored pixel file:')
        print(' ' + self.outlet_file)
        #----------------------------------------------
        print(hash_line) 
        print(' ')
        raise IOError('Could not find outlet_file.')

        #------------------------------------------------------
        # This led to more problems than it fixed. (11/11/16)
        #------------------------------------------------------
        # print ' Creating default version of file. '
        # print hash_line 
        # print ' '
        # write_outlet_file( self ) 

##        file_unit = open(self.outlet_file, 'r')
##        cfg.skip_header(file_unit, n_lines=6)
    
    file_unit = open(self.outlet_file, 'r')    
    lines     = file_unit.readlines()
    file_unit.close()
    lines = lines[6:]   # (skip over 6 header lines)
    n_lines = len(lines)
    
    outlet_cols   = np.zeros(n_lines, dtype='Int64')  ###
    outlet_rows   = np.zeros(n_lines, dtype='Int64')
    basin_areas   = np.zeros(n_lines, dtype='Float64')
    basin_reliefs = np.zeros(n_lines, dtype='Float64')
    
    n = 0
    for line in lines:
        words = line.split()    
        if (len(words) >= 4):
            outlet_cols[n]   = np.int64(np.float64(words[0]))
            outlet_rows[n]   = np.int64(np.float64(words[1]))
            basin_areas[n]   = np.float64(words[2])
            basin_reliefs[n] = np.float64(words[3])
            n += 1
    n_outlets = n

    #------------------------------------------------    
    # Save area and relief of first basin into self
    #------------------------------------------------
    self.basin_area   = basin_areas[0]
    self.basin_relief = basin_reliefs[0]

    #----------------------------------------------------
    # First, compute IDs as a 1D array of long-integer,
    # calendar-style indices
    #----------------------------------------------------
    # NB! Must return as Int32 vs. Int64 currently !!!!     SEEMS OBSOLETE.
    #----------------------------------------------------
    # (11/20/11) If outlet_file only has one outlet,
    # then using int32() here converts vector with one
    # element to a scalar and produces an error.
    #-----------------------------------------------------
    ## self.outlet_IDs = int32(self.outlet_rows * self.nx) + int32(self.outlet_cols)
    outlet_IDs = (outlet_rows * self.nx) + outlet_cols
    outlet_ID  = outlet_IDs[0]
    ## print 'outlet_ID =', outlet_ID

    #------------------------------------------
    # Are all the outlet IDs inside the DEM ?
    #------------------------------------------
    OK = check_outlet_IDs( outlet_IDs, self.rti.n_pixels )
    if not(OK):
        print('ERROR: Some outlet_IDs lie outside of DEM.')
        print(' ')
        return
        
    #-------------------------------------------
    # Save IDs as a tuple of row indices and
    # calendar indices, "np.where" style
    #-------------------------------------------
    self.n_outlets  = n_outlets   ## (new; 9/19/14)
    self.outlet_IDs = (outlet_rows,    outlet_cols)
    self.outlet_ID  = (outlet_rows[0], outlet_cols[0])
    ## self.outlet_IDs = (outlet_IDs / self.nx, outlet_IDs % self.nx)
    ## self.outlet_ID  = (outlet_ID  / self.nx, outlet_ID  % self.nx)

#   read_outlet_file()
#-------------------------------------------------------------------
def check_outlet_IDs( outlet_IDs, n_pixels ):

    OK = True
    wbad  = np.where(np.logical_or((outlet_IDs < 0), \
                                    (outlet_IDs > (n_pixels - 1))))
    n_bad = np.size(wbad[0])
    
    if (n_bad != 0):    
        msg = np.array(['SORRY, ', ' ', \
                     'One or more of the monitored pixel IDs ', \
                     'are not in the range of valid values. ', \
                     ' ', \
                     'You can use hydrologic GIS software to get ', \
                     'the outlet ID for the main basin. ', ' '])
        result = GUI_Message(msg, INFO=True, TITLE='Missing Input')
        OK = False
        
    return OK

#   check_outlet_IDs()
#-------------------------------------------------------------------
def read_main_basin_IDs( self ):

    #----------------------------------------------------------
    # Note: outlet_IDs and basin_IDs are 1D arrays or vectors
    #       but have very different sizes.  Recall that they
    #       are (row,col) tuples and not just long ints.
    #----------------------------------------------------------
    # Note:  "basin_IDs" gives IDs of all grid cells that lie
    #        within the basin that drains to outlet_ID.
    #----------------------------------------------------------
    self.basin_RTM_file = (self.in_directory +
                           self.site_prefix + '_basin.rtm')
    
    #----------------------------------
    # Read basin pixels from RTM_file
    #----------------------------------
    RTM_unit = open(self.basin_RTM_file, 'rb')
    RTM_filesize = os.path.getsize(RTM_unit.name)
    n_IDs = RTM_filesize / int32(4)

    basin_IDs = np.fromfile(RTM_unit, count=n_IDs, dtype='Int32')
    if (self.rti.SWAP_ENDIAN):
        basin_IDs.byteswap(True)
    RTM_unit.close()

    #-----------------------------------------
    # NB! basin_IDs is now 1D array of longs
    #-----------------------------------------
    wb  = np.where(basin_IDs >= 0)
    nwb = np.size(wb)   # (see note)
    ### nwb = size(wb[0])
    if (nwb != 0):    
        basin_IDs = basin_IDs[wb]
    else:    
        basin_IDs = self.outlet_ID  ###

    #------------------------------------------
    # Save IDs as a 1D array of long-integer,
    # calendar-style indices
    #------------------------------------------
    self.basin_IDs = basin_IDs    # (use Q.flat[basin_IDs])

    #-------------------------------------------
    # Save IDs as a tuple of row indices and
    # calendar indices, "np.where" style
    #-------------------------------------------        
    ### nx = self.nx
    ### self.basin_IDs = (basin_IDs / nx, basin_IDs % nx)

#   read_main_basin_IDs()
#-------------------------------------------------------------------