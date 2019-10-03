
# (2/3/13) read_source_data() and read_sink_data() get "dt" from
#          source_file or sink_file vs. channels comp, but what
#          about canals ? Can't use "channel_dt" as before.
# (2/3/13) "vol@channel" is now obtained as a reference.

## THIS STILL USES FUNCTIONS IN idl_func (11/14/11)
##
## NB! CFG file does not have "n_steps" so it is set
##     to 1 in BMI_base.py.  Maybe it should be
##     set to max of "nd_max" values that appear in
##     read_source_data(), read_sink_data(), etc. (11/14/11)
##============================================================
#
#  Copyright (c) 2001-2014, Scott D. Peckham
#
#  Sep 2014.  New standard names and BMI updates and testing.
#  Nov 2013.  Converted TopoFlow to a Python package.
#  Feb 2013.  Adapted to use EMELI framework.
#  Oct 2012.  CSDMS Standard Names and BMI.
#  May 2010.  Changes to initialize() and read_cfg_file().
#  Jul 2009.  Updates.
#  May 2009.  Updates.
#  Jan 2009.  Converted from IDL to Python with I2PY.
#
#---------------------------------------------------------------------
# Notes:  Maybe replace "dur_sums" approach with the same approach
#         now used in precip.py ??

#         Make sure volume in caller gets updated correctly.
#---------------------------------------------------------------------
#
#  class diversions_component:  (inherits from BMI_base.py)
# 
#      initialize()
#      update()
#      finalize()
#      set_computed_input_vars()
#---------------------------------
#      read_input_files()
#      read_source_data()
#      read_sink_data()
#      read_canal_data()
#----------------------------
#      update_sources()
#      update_sinks()
#      update_canals()
#
#--------------------------------------------------------------------

import numpy as np
import glob
import os

from topoflow.utils import BMI_base
from topoflow.utils import cfg_files as cfg
from topoflow.utils import idl_func  # (still needed for idl.readf...)

#---------------------------------------------------------------------
class diversions_component( BMI_base.BMI_component ):
 
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        if not(SILENT):
            print(' ')
            print('Diversions component: Initializing...')
        
        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file 
        #-----------------------------------------------
        # self.set_constants()
        self.initialize_config_vars() 
        self.read_grid_info()  # (need this, 5/19/10)
        self.initialize_basin_vars()  # (5/14/10)
        
        ############################################################   
        # With new framework approach, we can't request the time
        # step for a specific component as "dt" (due to conflicts).
        ############################################################   
        # The Diversions component "dt" should probably match that
        # of the Channels component.  (Must be named "dt" vs.
        # "canal_dt".)
        ############################################################
        self.initialize_time_vars()
        
        #-----------------------------------------------------
        # These are used by the Channels component (9/22/14)
        # to tell if diversions are available and "on".
        #-----------------------------------------------------
        if not(self.use_canals):
            self.n_canals = self.initialize_scalar( 0, dtype='int32')
        if not(self.use_sinks):
            self.n_sinks = self.initialize_scalar( 0, dtype='int32')
        if not(self.use_sources):
            self.n_sources = self.initialize_scalar( 0, dtype='int32')

        #-------------------------------
        # Return if no method selected
        #-------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print('Diversions component: Disabled in CFG file.')
            self.n_canals          = self.initialize_scalar( 0, dtype='int32')  # int
            self.n_sinks           = self.initialize_scalar( 0, dtype='int32')
            self.n_sources         = self.initialize_scalar( 0, dtype='int32')
            self.Q_canals_fraction = self.initialize_scalar( 0, dtype='float64')
            self.Q_canals_out      = self.initialize_scalar( 0, dtype='float64')
            self.Q_sinks           = self.initialize_scalar( 0, dtype='float64')
            self.Q_sources         = self.initialize_scalar( 0, dtype='float64')
            self.canals_in_x       = self.initialize_scalar( 0, dtype='float64')
            self.canals_in_y       = self.initialize_scalar( 0, dtype='float64')
            self.canals_out_x      = self.initialize_scalar( 0, dtype='float64')
            self.canals_out_y      = self.initialize_scalar( 0, dtype='float64')
            self.sinks_x           = self.initialize_scalar( 0, dtype='float64')
            self.sinks_y           = self.initialize_scalar( 0, dtype='float64')
            self.sources_x         = self.initialize_scalar( 0, dtype='float64')
            self.sources_y         = self.initialize_scalar( 0, dtype='float64')
            self.DONE   = True
            self.status = 'initialized'  # (OpenMI 2.0 convention)
            return

        #----------------------------------------
        # Initialize all time-related variables
        #----------------------------------------
        self.initialize_time_vars()
        
        #-----------------------------------------------
        # Read from files as needed to initialize vars 
        #--------------------------------------------------
        # source_files and sink_files have their own "dt"
        # which will override the default above. (2/3/13)
        #--------------------------------------------------
        self.read_input_files()
        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
    #   initialize()
    #--------------------------------------------------------------------------
    ## def update(self, dt=-1.0, time_seconds=None):
    def update(self, dt=-1.0):

        if (self.comp_status == 'Disabled'):
            return
        self.status = 'updating'  # (OpenMI 2.0 convention)

        #-----------------------------------------------------
        # Update info from all sources, sinks and diversions
        #-----------------------------------------------------
        # print '### Calling update_sources()...'
        self.update_sources()
        # print '### Calling update_sinks()...'
        self.update_sinks()
        # print '### Calling update_canals()...'
        self.update_canals()

        #------------------------
        # Update internal clock
        #------------------------
        self.update_time()
        self.status = 'updated'  # (OpenMI 2.0 convention)
        
    #   update()
    #--------------------------------------------------------------------------
    def finalize(self):
        
        self.status = 'finalizing'  # (OpenMI)
##        if (self.comp_status == 'Enabled'):
##            self.close_input_files()
##            self.close_output_files()
        self.status = 'finalized'  # (OpenMI)

        self.print_final_report(comp_name='Diversions component')
        ## print 'Diversions component: Finished.'
        
    #   finalize()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        #---------------------------------------------------------------    
        # Note: The initialize() method calls initialize_config_vars()
        #       (in BMI_base.py), which calls this method at the end.
        #       But read_input_files() has not been called yet.
        #       See initialize_computed_vars().
        #--------------------------------------------------------------
        self.use_sources = (self.use_sources == 'Yes')
        self.use_sinks   = (self.use_sinks   == 'Yes')
        self.use_canals  = (self.use_canals  == 'Yes')

    #   set_computed_input_vars()
    #--------------------------------------------------------------------------
    def read_input_files(self):
    
        self.source_file = self.in_directory + self.source_file
        self.sink_file   = self.in_directory + self.sink_file
        self.canal_file  = self.in_directory + self.canal_file
 
        self.read_source_data()
        self.read_sink_data()
        self.read_canal_data()

    #   read_input_files()
    #--------------------------------------------------------------------------                
    def read_source_data(self):

        #-----------------------------------------------------------------
        # Notes:  This assumes that source_file is organized as follows:
        #            ID         (source pixel ID as long integer)
        #            nd         (number of durations and Q values)
        #            durations  (vector of durations in minutes)
        #            Q_sources  (vector of discharges in m^3/sec)
        #-----------------------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        if not(self.use_sources): return
        
        #-----------------------------
        # Can source_file be found ?
        #-----------------------------
        f = glob.glob(self.source_file)
        if (len(f) == 0):    
            msg = array(['ERROR: source file not found. ', ' ', \
                         'The file: ', '  ' + self.source_file, \
                         'was not found in the working directory. ', \
                         ' '])
            self.use_sources = False
            return
        
        #-----------------------------
        # Count lines in source_file
        #-----------------------------
        #Count_Lines, n_lines, self.source_file, /SILENT
        #n_IDs = (n_lines / 4L)
        
        #------------------------------------
        # Read the lines with the nd values
        # and find the max, called nd_max
        #------------------------------------
        nd     = np.int32(0)
        nd_max = np.int32(0)
        n_IDs  = np.int32(0)
        #------------------------------------------
        file_unit = open(self.source_file, 'r')
##        while (True):
##            line1 = file_unit.readline()
##            if (line1 == ''): break
##            nd    =
##            line3 = file_unit.readline()
##            line4 = file_unit.readline()
##            #-----------------------------
##            nd_max = np.maximum(nd_max, nd)
##            n_IDs += 1
##        file_unit.close()

        line1 = ''
        line3 = ''
        line4 = ''        
        while not(idl_func.eof(file_unit)):
            line1 = idl_func.readf(file_unit, line1)
            nd    = idl_func.readf(file_unit, nd)
            line3 = idl_func.readf(file_unit, line3)
            line4 = idl_func.readf(file_unit, line4)
            #-----------------------------------------
            nd_max = np.maximum(nd_max, nd)
            n_IDs = (n_IDs + np.int32(1))   #***
        file_unit.close()
  
        #--------------------
        # Initialize arrays
        #--------------------
        self.source_IDs  = np.zeros([n_IDs], dtype='Int32')
        self.nd_vals     = np.zeros([n_IDs], dtype='Int32')   #*****
        self.dur_sources = np.zeros([nd_max, n_IDs], dtype='Float32')
        self.Q_sources   = np.zeros([nd_max, n_IDs], dtype='Float32')
        
        #--------------------
        # Open file to read
        #--------------------
        file_unit = open(self.source_file, 'r')
        
        #----------------------
        # Read data from file
        #----------------------
        ID = np.int32(0)
        nd = np.int32(0)
        k  = np.int32(0)
        while not(idl_func.eof(file_unit)):
            ID = idl_func.readf(file_unit, ID)
            self.source_IDs[k] = ID
            #----------------------------------
            nd = idl_func.readf(file_unit, nd)
            self.nd_vals[k] = nd
            #----------------------------------
            durs = np.zeros([nd], dtype='Float32')
            durs = idl_func.readf(file_unit, durs)
            self.dur_sources[np.int32(0):((nd - np.int32(1)))+1,k] = durs
            #----------------------------------
            Q = np.zeros([nd], dtype='Float32')
            Q = idl_func.readf(file_unit, Q)
            self.Q_sources[np.int32(0):((nd - np.int32(1)))+1,k] = Q
            #----------------------------------
            k = (k + np.int32(1))
        
        #-----------------------
        # Close the input file
        #-----------------------
        file_unit.close()
        
        #------------------------------------
        # Compute partial sums of durations
        #------------------------------------
        dur_sum_sources = np.zeros([nd_max + np.int32(1), n_IDs], dtype='Float32')
        for k in range(n_IDs):
            for i in range(nd_max):
                dur_sum_sources[i + 1,k] = dur_sum_sources[i,k] + dur_sources[i,k]
        dur_sum_sources = dur_sum_sources[1:(nd_max)+1,:]  #(remove leading zero)
        self.dur_sum_sources  = dur_sum_sources       
   
    #   read_source_data()
    #--------------------------------------------------------------------------
    def read_sink_data(self):
        
        #----------------------------------------------------------------
        # Notes:  This assumes that sink_file is organized as follows:
        #            ID         (sink pixel ID as long integer)
        #            nd         (number of durations and Q values)
        #            durations  (vector of durations in minutes)
        #            Q_sinks    (vector of discharges in m^3/sec)
        #----------------------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        if not(self.use_sinks): return
    
        #---------------------------
        # Can sink_file be found ?
        #---------------------------
        f = glob.glob(self.sink_file)
        c = len(f)
        if (c == 0):    
            msg = array(['ERROR: sink file not found. ', ' ', \
                         'The file: ', '  ' + self.sink_file, \
                         'was not found in the working directory. ', \
                         ' '])
            result = GUI_Message(msg, INFO=True)
            self.use_sinks = False
            return
        
        #-----------------------------
        # Count lines in source_file
        #-----------------------------
        #Count_Lines, n_lines, self.sink_file, /SILENT
        #n_IDs = (n_lines / 4L)
        
        #------------------------------------
        # Read the lines with the nd values
        # and find the max, called nd_max
        #------------------------------------
        nd = np.int32(0)
        nd_max = np.int32(0)
        n_IDs = np.int32(0)  #******
        line1 = ''
        line3 = ''
        line4 = ''
        #---------------------------
        file_unit = open(self.sink_file, 'r')
        while not(idl_func.eof(file_unit)):
            line1 = idl_func.readf(file_unit, line1)
            nd    = idl_func.readf(file_unit, nd)
            line3 = idl_func.readf(file_unit, line3)
            line4 = idl_func.readf(file_unit, line4)
            #-----------------------------------------
            nd_max = (np.maximum(nd_max, nd))
            n_IDs = (n_IDs + np.int32(1))   #***
        file_unit.close()
        
        #--------------------
        # Initialize arrays
        #--------------------
        self.sink_IDs  = np.zeros([n_IDs], dtype='Int32')
        self.nd_vals   = np.zeros([n_IDs], dtype='Int32')   #*****
        self.dur_sinks = np.zeros([nd_max, n_IDs], dtype='Float32')
        self.Q_sinks   = np.zeros([nd_max, n_IDs], dtype='Float32')
        
        #--------------------
        # Open file to read
        #--------------------
        file_unit = open(self.sink_file, 'r')
        
        #----------------------
        # Read data from file
        #----------------------
        ID = np.int32(0)
        nd = np.int32(0)
        k  = np.int32(0)
        while not(idl_func.eof(file_unit)):
            ID = idl_func.readf(file_unit, ID)
            self.sink_IDs[k] = ID
            #----------------------------------
            nd = idl_func.readf(file_unit, nd)
            self.nd_vals[k] = nd
            #----------------------------------
            durs = np.zeros([nd], dtype='Float32')
            durs = idl_func.readf(file_unit, durs)
            self.dur_sinks[np.int32(0):((nd - np.int32(1)))+1,k] = durs
            #----------------------------------
            Q = zeros([nd], dtype='Float32')
            Q = idl_func.readf(file_unit, Q)
            self.Q_sinks[np.int32(0):((nd - np.int32(1)))+1,k] = Q
            #----------------------------------
            k = (k + np.int32(1))
        
        #-----------------------
        # Close the input file
        #-----------------------
        file_unit.close()
        
        #------------------------------------
        # Compute partial sums of durations
        #------------------------------------
        dur_sum_sinks = zeros([nd_max + np.int32(1), n_IDs], dtype='Float32')
        for k in range(n_IDs):
            for i in range(nd_max):
                dur_sum_sinks[i + 1,k] = dur_sum_sinks[i,k] + dur_sinks[i,k]
        dur_sum_sinks = dur_sum_sinks[1:(nd_max)+1,:]  #(remove leading zero)
        self.dur_sum_sinks = dur_sum_sinks
        
    #   read_sink_data()
    #--------------------------------------------------------------------------
    def read_canal_data(self):

        #-------------------------------------------------------------------
        # Note:  Q_canals is same at upstream and downstream ends, but the
        #        downstream end lags the upstream end by the travel time
        #        from in_ID to out_ID.  As a result, the duration and Q
        #        vector for the downstream end are computed from those of
        #        the upstream end, and the travel time, td, as:
        #            Q_out   = [0,  Q_in]
        #            dur_out = [td, dur_in]
        #            dur_sum_out = [0, dur_sum_in] + td

        #        Rather than create the dur_sum_canals_out and
        #        Q_canals_out vectors, can construct them in Update_Canals.
        #-------------------------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        if not(self.use_canals): return
    
        #---------------------------
        # Can canal_file be found ?
        #---------------------------
        f = glob.glob(self.canal_file)
        c = len(f)
        if (c == 0):    
            msg = array(['ERROR: canal file not found. ', ' ', \
                         'The file: ', '  ' + self.canal_file, \
                         'was not found in the working directory. ', \
                         ' '])
            result = GUI_Message(msg, INFO=True)
            self.use_canals = False
            return
        
        #-----------------------------
        # Count lines in source_file
        #-----------------------------
        #Count_Lines, n_lines, canal_file, /SILENT
        #n_IDs = (n_lines / 6L)
        
        #------------------------------------
        # Read the lines with the nd values
        # and find the max, called nd_max
        #------------------------------------
        nd = np.int32(0)
        nd_max = np.int32(0)
        n_IDs = np.int32(0)  #******
        line1 = ''
        line2 = ''
        line3 = ''
        line5 = ''
        line6 = ''

        file_unit = open(self.canal_file, 'r')
        while not(idl_func.eof(file_unit)):
            line1 = idl_func.readf(file_unit, line1)
            line2 = idl_func.readf(file_unit, line2)
            line3 = idl_func.readf(file_unit, line3)
            nd    = idl_func.readf(file_unit, nd)
            line5 = idl_func.readf(file_unit, line5)
            line6 = idl_func.readf(file_unit, line6)
            #-----------------------------------------
            nd_max = np.maximum( nd_max, nd )
            n_IDs = (n_IDs + np.int32(1))   #***
        file_unit.close()
        
        #--------------------
        # Initialize arrays
        #--------------------
        self.canal_in_IDs  = np.zeros([n_IDs], dtype='Int32')
        self.canal_out_IDs = np.zeros([n_IDs], dtype='Int32')
        self.canal_t_vals  = np.zeros([n_IDs], dtype='Int32')
        self.nd_vals       = np.zeros([n_IDs], dtype='Int32')   #*****
        self.dur_canals    = np.zeros([nd_max, n_IDs], dtype='Float32')
        self.Q_canals_in   = np.zeros([nd_max, n_IDs], dtype='Float32')
        
        #--------------------
        # Open file to read
        #--------------------
        file_unit = open(self.canal_file, 'r')
        
        #----------------------
        # Read data from file
        #----------------------
        ID   = np.int32(0)
        nd   = np.int32(0)
        k    = np.int32(0)
        tval = np.float32(0.0)
        while not(idl_func.eof(file_unit)):
            ID = idl_func.readf(file_unit, ID)
            self.canal_in_IDs[k] = ID
            #-----------------------------------
            ID = idl_func.readf(file_unit, ID)
            self.canal_out_IDs[k] = ID
            #--------------------------------------
            tval = idl_func.readf(file_unit, tval)
            self.canal_t_vals[k] = tval
            #-----------------------------------
            nd = idl_func.readf(file_unit, nd)
            self.nd_vals[k] = nd
            #-----------------------------------
            durs = np.zeros([nd], dtype='Float32')
            durs = idl_func.readf(file_unit, durs)
            self.dur_canals[np.int32(0):((nd - np.int32(1)))+1,k] = durs
            #-----------------------------------
            Q = np.zeros([nd], dtype='Float32')
            Q = idl_func.readf(file_unit, Q)
            self.Q_canals_in[np.int32(0):((nd - np.int32(1)))+1,k] = Q
            #-----------------------------------
            k = (k + np.int32(1))
        
        #-----------------------
        # Close the input file
        #-----------------------
        file_unit.close()
        
        #------------------------------------
        # Compute partial sums of durations
        #------------------------------------
        dur_sum_canals_in = np.zeros([nd_max + np.int32(1), n_IDs], dtype='Float32')
        for k in range(n_IDs):
            for i in range(nd_max):
                dur_sum_canals_in[i + 1,k] = dur_sum_canals_in[i,k] + dur_canals[i,k]
        dur_sum_canals_in = dur_sum_canals_in[1:(nd_max)+1,:]  #(remove leading 0)
        self.dur_sum_canals_in = dur_sum_canals_in
        
        #------------------------------------------------
        # Compute Q_canals_out (Q[t] at downstream end)
        # See Notes;  now computed in Update_Canals
        #------------------------------------------------
        #*** Q_canals_out =
        #*** dum_sum_canals_out =
        
    #   read_canal_data()
    #--------------------------------------------------------------------------
    def update_sources(self):

        if not(self.use_sources): return     
        # n = size(self.source_IDs)
        n = self.source_IDs.size
        
        for k in range(n):
            dvol = ( self.dt * self.Q_sources[self.time_index, k] )
            self.vol.flat[self.source_IDs[k]] += dvol
        
##        for k in xrange(n):           
##            w  = np.where( time_min < self.dur_sum_sources[:,k] )
##            nw = np.size( w[0] )
##            if (nw != 0):    
##                dvol = (dt * self.Q_sources[w[0], k])
##                vol[self.source_IDs[k]] += dvol
        
    #   update_sources()
    #--------------------------------------------------------------------------
    def update_sinks(self):

        if not(self.use_sinks): return     
        # n = size(self.sink_IDs)
        n = self.sink_IDs.size

        for k in range(n):
            dvol = (self.dt * self.Q_sinks[self.time_index, k])
            self.vol.flat[self.sink_IDs[k]] -= dvol
            self.vol = np.maximum(self.vol, 0.0)
                
##        for k in xrange(n):
##            w  = np.where( time_min < self.dur_sum_sinks[:,k] )
##            nw = np.size( w[0] )
##            if (nw != 0):    
##                dvol = (dt * self.Q_sinks[w[0],k])
##                vol[self.sink_IDs[k]] = np.maximum((vol[self.sink_IDs[k]] - dvol), 0.0)
        
    #   update_sinks()
    #--------------------------------------------------------------------------
    def update_canals(self):

        #------------------------------------------------------------------
        # Note: Q_canals is same at upstream and downstream ends, but the
        #       downstream end lags the upstream end by the travel time
        #       from in_ID to out_ID.  As a result, the duration and Q
        #       vector for the downstream end are computed from those of
        #       the upstream end, and the travel time, td, as:
        #           Q_out   = [0,  Q_in]
        #           dur_out = [td, dur_in]
        #           dur_sum_out = [0, dur_sum_in] + td
        #------------------------------------------------------------------
        if not(self.use_canals): return
        n = np.size( self.canal_in_IDs )

        #--------------------------------
        # Process upstream IDs as sinks
        #--------------------------------
        for k in range(n):   
            dvol = (self.dt * self.Q_canals_in[self.time_index, k])
            self.vol.flat[self.canal_in_IDs[k]] -= dvol
            self.vol = np.maximum( self.vol, 0.0 )
        
        #------------------------------------
        # Process downstream IDs as sources
        # Must account for time lag in Q.
        #------------------------------------
        for k in range(n):
            #--------------------------------------------
            # Compute Q_canals_out & dur_sum_canals_out
            #--------------------------------------------
            dur_sum_canals_out = np.array([0.0, self.dur_sum_canals_in[:,k]]) + self.canal_t_vals[k]
            Q_canals_out       = np.array([0.0, self.Q_canals_in[:,k]])
  
            dvol = (self.dt * Q_canals_out[self.time_index])
            self.vol[self.canal_out_IDs[k]] += dvol

                
##        #--------------------------------
##        # Process upstream IDs as sinks
##        #--------------------------------
##        for k in xrange(n):
##            w  = np.where( time_min < self.dur_sum_canals_in[:,k] )
##            nw = np.size( w[0] )
##            if (nw != 0):    
##                dvol = (dt * self.Q_canals_in[w[0],k])
##                vol[self.canal_in_IDs[k]] = np.maximum((vol[self.canal_in_IDs[k]] - dvol), 0.0)
##        
##        #------------------------------------
##        # Process downstream IDs as sources
##        # Must account for time lag in Q.
##        #------------------------------------
##        for k in xrange(n):
##            #--------------------------------------------
##            # Compute Q_canals_out & dur_sum_canals_out
##            #--------------------------------------------
##            dur_sum_canals_out = np.array([0.0, self.dur_sum_canals_in[:,k]]) + self.canal_t_vals[k]
##            Q_canals_out       = np.array([0.0, self.Q_canals_in[:,k]])
##
##            w  = np.where( time_min < dur_sum_canals_out )
##            nw = np.size( w[0] )
##            if (nw != 0):    
##                dvol = (dt * Q_canals_out[w[0]])
##                vol[self.canal_out_IDs[k]] += dvol

    #   update_canals()
    #--------------------------------------------------------------------------


    

