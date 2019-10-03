
# (2/3/13) Get "dt" from source_file or sink_file vs.
#          channels comp, but what about canals ?
#
# (2/3/13) Still need to remove calls to "idl_funcs".

########################################################
#        
# Copyright (c) 2010-2014, Scott D. Peckham
#
# Sept 2014.  Big changes so Channels component now requests
#             what is needed from Diversions component.
#
# January 2013   (Revised handling of input/output names).
# October 2012   (CSDMS Standard Names with BMI)
# Jan-Feb 2010   (started from diversions_base.py)
# May 2010 (changes to unit_test())

#---------------------------------------------------------------------
# Notes:  This component is written so that only a small amount
#         of data is retrieved from, altered and then passed
#         back to the Channels component.  This is accomplished
#         by using new interface functions in CSDMS_base.py,
#         namely:
#             get_values_in_grid_double()
#             set_values_in_grid_double()
#             get_values_in_grid_int()
#             set_values_in_grid_int()
#         Note that these also had to be added to the IRFPort
#         for the TopoFlow CCA project, as defined in the file
#         topoflow3.IRFPort.sidl.
#
#         The old method required the Diversion component to
#         retrieve the entire Q and vol grid from Channel component
#         at each timestep, alter it, and then pass it back.
#         (02/18/10)
#
#         Part of the philosophy of this version is that only
#         tiny changes to the code of the Channels component should
#         be necessary.  An intermediate approach (DIV_METHOD1)
#         (before the above functions were added), required a new
#         function "update_diversions()" to be added to the
#         Channels component.  However, that version seems to be
#         somewhat faster.  For the "test_plane_canal" test, the
#         run time was 0.38 secs vs. 0.44 secs (on beach) for this
#         new method.  This extra cost should only be due to the
#         extra interface function calls, and should therefore be
#         a fixed cost that doesn't increase with grid size.  This
#         remains to be tested, however.  To test, simply swap
#           channels_base_DIV_METHOD1.py and
#           diversions_fraction_method_DIV_METHOD1.py for
#         channels_base.py and this file.
#
#         cp.update_discharge() now calls dp.update(). (2/1/10)
#---------------------------------------------------------------------
#
#  class diversions_component:  (inherits from diversions_base.py)
#
#      get_component_name()
#      get_attribute()          # (10/26/11)
#      get_input_var_names()    # (5/16/12, Bolton)
#      get_output_var_names()   # (5/16/12, Bolton)
#      get_var_name()           # (5/16/12, Bolton)
#      get_var_units()          # (5/16/12, Bolton)
#      update()
#----------------------------
#      read_input_files()
#      read_source_data()
#      read_sink_data()
#      read_canal_data()
#----------------------------
#      update_sources()
#      update_sinks()
#      update_canals()
#
#---------------------------------------------------------------------

import numpy as np

import glob
import os

from topoflow.components import diversions_base

from topoflow.utils import cfg_files as cfg
from topoflow.utils import tf_utils

#---------------------------------------------------------------------
class diversions_component( diversions_base.diversions_component ):

    #-----------------------------------------------------------------
    _att_map = {
        'model_name':         'Diversions_Fraction_Method',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'uniform',   ## (or "none" ?)
        'time_step_type':     'fixed',
        'step_method':        'explicit',  ## (or "none" ?)
        #------------------------------------------------------
        'comp_name':          'Diversions',     # CHANGE LATER ?
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Diversions_Fraction_Method.cfg.in',
        'cfg_extension':      '_diversions_fraction_method.cfg',
        # 'cmt_var_prefix':     '/Diversions/Input/Var/',     # (see Diversions_Fraction_Method.xml)
        'cmt_var_prefix':     '/DiversionsFraction/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Diversions_Fraction_Method.xml',
        'dialog_title':       'Diversions: Fraction Method Parameters',
        'time_units':         'seconds' }

    _input_var_names = []
#         'canals_entrance_water__volume_flow_rate' ]        # Q_canals_in

    _output_var_names = [
        'canals__count',                                   # n_canals
        'canals_entrance_water__volume_fraction',          # Q_canals_fraction
        'canals_entrance__x_coordinate',                   # canals_in_x
        'canals_entrance__y_coordinate',                   # canals_in_y
        'canals_exit_water__volume_flow_rate',             # Q_canals_out
        'canals_exit__x_coordinate',                       # canals_out_x
        'canals_exit__y_coordinate',                       # canals_out_y
        'model__time_step',                                # dt
        'sinks_water__volume_flow_rate',                   # Q_sinks
        'sinks__count',                                    # n_sinks
        'sinks__x_coordinate',                             # Q_sinks_x
        'sinks__y_coordinate',                             # Q_sinks_y
        'sources_water__volume_flow_rate',                 # Q_sources
        'sources__count',                                  # n_sources
        'sources__x_coordinate',                           # Q_sources_x
        'sources__y_coordinate' ]                          # Q_sources_y

    _var_name_map = {
        'canals_entrance_water__volume_flow_rate': 'Q_canals_in',
        #-----------------------------------------------------------
        'canals__count':                          'n_canals',
        'canals_entrance__x_coordinate':          'canals_in_x',
        'canals_entrance__y_coordinate':          'canals_in_y',
        'canals_entrance_water__volume_fraction': 'Q_canals_fraction',
        'canals_exit__x_coordinate':              'canals_out_x',
        'canals_exit__y_coordinate':              'canals_out_y',
        'canals_exit_water__volume_flow_rate':    'Q_canals_out',
        'model__time_step':                       'dt',
        'sinks__count':                           'n_sinks',
        'sinks__x_coordinate':                    'sinks_x',
        'sinks__y_coordinate':                    'sinks_y',
        'sinks_water__volume_flow_rate':          'Q_sinks',
        'sources__count':                         'n_sources',
        'sources__x_coordinate':                  'sources_x',
        'sources__y_coordinate':                  'sources_y',
        'sources_water__volume_flow_rate':        'Q_sources' }

    _var_units_map = {
        'canals_entrance_water__volume_flow_rate': 'm3 s-1',
        #------------------------------------------------------
        'canals__count':                          '1',
        'canals_entrance__x_coordinate':          'm',
        'canals_entrance__y_coordinate':          'm',
        'canals_entrance_water__volume_fraction': '1',
        'canals_exit__x_coordinate':              'm',
        'canals_exit__y_coordinate':              'm',
        'canals_exit_water__volume_flow_rate':    'm3 s-1',
        'model__time_step':                       's',
        'sinks__count':                           '1',
        'sinks__x_coordinate':                    'm',
        'sinks__y_coordinate':                    'm',
        'sinks_water__volume_flow_rate':          'm3 s-1',
        'sources__count':                         '1',
        'sources__x_coordinate':                  'm',
        'sources__y_coordinate':                  'm',
        'sources_water__volume_flow_rate':        'm3 s-1' }
    
    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Diversions_Fraction_Method'

    #   get_component_name()     
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print('###################################################')
            print(' ERROR: Could not find attribute: ' + att_name)
            print('###################################################')
            print(' ')

    #   get_attribute() 
    #-------------------------------------------------------------------
    def get_input_var_names(self):
        
        #--------------------------------------------------------
        # Note: These are currently variables needed from other
        #       components vs. those read from files or GUI.
        #--------------------------------------------------------
        
        return self._input_var_names
    
    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):
  
        return self._output_var_names

    #   get_output_var_names()
    #-------------------------------------------------------------------
    def get_var_name(self, long_var_name):

        return self._var_name_map[ long_var_name ]

    #   get_var_name()
    #-------------------------------------------------------------------
    def get_var_units(self, long_var_name):
 
        return self._var_units_map[ long_var_name ]
    
    #   get_var_units()
    #-------------------------------------------------------------------
    def update(self, time_seconds=None):
  
        if (self.comp_status == 'Disabled'):
            return
        self.status = 'updating'  # (OpenMI 2.0 convention)
        
        #-----------------------------------------------
        # Update self.vol with inputs/outputs from all
        # sources, sinks and diversions
        #-----------------------------------------------
        # print '#### Calling update_sources()...'
        self.update_sources()
        # print '#### Calling update_sinks()...'
        self.update_sinks()
        # print '#### Calling update_canals()...'
        self.update_canals()

        #------------------------
        # Update internal clock
        #------------------------
        # print '#### Calling update_time()...'
        self.update_time()
        self.status = 'updated'  # (OpenMI 2.0 convention)
        
    #   update()
    #--------------------------------------------------------------------------
    def read_source_data(self):

        #------------------------------------------------------------
        # Notes:  Assume that source_file contains key-value pairs,
        #         starting with "n_sources:", "nt_max" and "dt:",
        #         followed by "n_sources" blocks of the form:
        #
        #         source_ID: (source pixel ID as long integer)
        #         nt:        (number of discharge (Q) values)
        #         Q:         (vector of discharges in m^3/s)
        #------------------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        if not(self.use_sources): return
        
        #-----------------------------
        # Can source_file be found ?
        #-----------------------------
        FOUND = tf_utils.file_exists( self.source_file )
        if not(FOUND):
            self.use_sources = False
            return

        #-------------------------
        # Open the "source_file"
        #-------------------------
        file_unit = open(self.source_file, 'r')

        #----------------------------------------------------
        # Read number of sources, max number of timesteps
        # for any source and the common timestep, source_dt
        #----------------------------------------------------
        n_sources = cfg.read_value(file_unit, dtype='Int32')
        nt_max    = cfg.read_value(file_unit, dtype='Int32')
        source_dt = cfg.read_value(file_unit, dtype='Float64')
        self.source_dt =source_dt
        
        #--------------------
        # Initialize arrays
        #--------------------
        self.source_IDs     = np.zeros([n_sources], dtype='Int32')
        self.nt_sources     = np.zeros([n_sources], dtype='Int32')
        self.Q_sources_all  = np.zeros([n_sources, nt_max], dtype='Float64')
        self.n_sources      = n_sources
        self.nt_max_sources = nt_max
        
        #-----------------------------------
        # Read information for each source
        #-----------------------------------
        for k in range(n_sources):
            source_ID = cfg.read_value(file_unit, dtype='Int32')
            nt        = cfg.read_value(file_unit, dtype='Int32')
            Q_values  = cfg.read_list_after_key(file_unit, dtype='Float64')
            #---------------------------------------------------------------
            nQ        = np.size(Q_values)
            print('Diversions component: Read', nQ, 'Q_values for source.')
            #---------------------------------------------------------------            
            self.source_IDs[k]     = source_ID
            self.nt_sources[k]     = nt
            self.Q_sources_all[k,0:nt] = Q_values
 
        #-----------------------
        # Close the input file
        #-----------------------
        file_unit.close()

        #-------------------------------------
        # Compute xy coordinates for sources
        #-------------------------------------
        source_rows    = (self.source_IDs / self.nx)
        source_cols    = (self.source_IDs % self.nx)
        self.sources_x = (source_cols * self.dx)
        self.sources_y = (source_rows * self.dy)  
                    
    #   read_source_data()
    #--------------------------------------------------------------------------
    def read_sink_data(self):
        
        #------------------------------------------------------------
        # Notes:  Assume that source_file contains key-value pairs,
        #         starting with "n_sinks:", "nt_max" and "dt:",
        #         followed by "n_sinks" blocks of the form:
        #
        #         sink_ID: (sink pixel ID as long integer)
        #         nt:      (number of discharge (Q) values)
        #         Q:       (vector of discharges in m^3/s)
        #------------------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        if not(self.use_sinks): return
    
        #---------------------------
        # Can sink_file be found ?
        #---------------------------
        FOUND = tf_utils.file_exists( self.sink_file )
        if not(FOUND):
            self.use_sinks = False
            return

        #-----------------------
        # Open the "sink_file"
        #-----------------------
        file_unit = open(self.sink_file, 'r')

        #------------------------------------------------
        # Read number of sinks, max number of timesteps
        # for any sink and the common timestep, dt
        #------------------------------------------------
        n_sinks = cfg.read_value(file_unit, dtype='Int32')
        nt_max  = cfg.read_value(file_unit, dtype='Int32')
        sink_dt = cfg.read_value(file_unit, dtype='Float64')
        self.sink_dt = sink_dt
        
        #--------------------
        # Initialize arrays
        #--------------------
        self.sink_IDs     = np.zeros([n_sinks], dtype='Int32')
        self.nt_sinks     = np.zeros([n_sinks], dtype='Int32')
        self.Q_sinks_all  = np.zeros([n_sinks, nt_max], dtype='Float64')
        self.n_sinks      = n_sinks
        self.nt_max_sinks = nt_max
        
        #---------------------------------
        # Read information for each sink
        #---------------------------------
        for k in range(n_sinks):
            sink_ID  = cfg.read_value(file_unit, dtype='Int32')
            nt       = cfg.read_value(file_unit, dtype='Int32')
            Q_values = cfg.read_list_after_key(file_unit, dtype='Float64')
            #---------------------------------------------------------------
            nQ        = size(Q_values)
            print('Diversions component: Read', nQ, 'Q_values for sink.')
            #--------------------------------------------------------------- 
            self.sink_IDs[k]     = sink_ID
            self.nt_sinks[k]     = nt
            self.Q_sinks_all[k,0:nt] = Q_values

        #-----------------------
        # Close the input file
        #-----------------------
        file_unit.close()              

        #-----------------------------------
        # Compute xy coordinates for sinks
        #-----------------------------------
        sink_rows    = (self.sink_IDs / self.nx)
        sink_cols    = (self.sink_IDs % self.nx)
        self.sinks_x = (sink_cols * self.dx)
        self.sinks_y = (sink_rows * self.dy)  
        
    #   read_sink_data()
    #--------------------------------------------------------------------------
    def read_canal_data(self):

        #-------------------------------------------------------------------
        # Notes:  Assume that canal_file contains key-value pairs,
        #         starting with "n_canals:" and followed by "n_canals"
        #         blocks of the form:
        #            canal_in_ID:  (pixel ID as long integer)
        #            canal_out_ID: (pixel ID as long integer)
        #            Q_fraction:   (fraction to take from in_ID in [0,1])
        #            travel_time:  (canal travel time, in minutes)
        #
        #         nt_canals is computed as ceil(travel_time / cp.dt)
        #-------------------------------------------------------------------
        # Note:  Q_canals is same at upstream and downstream ends, but the
        #        downstream end lags the upstream end by the travel time
        #        from in_ID to out_ID.  As a result, the duration and Q
        #        vector for the downstream end are computed from those of
        #        the upstream end, and the travel time, td, as:
        #            Q_out   = [0,  Q_in]
        #            dur_out = [td, dur_in]
        #            dur_sum_out = [0, dur_sum_in] + td
        #
        #        Rather than create the dur_sum_canals_out and
        #        Q_canals_out vectors, can construct them in Update_Canals.
        #-------------------------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        if not(self.use_canals): return
    
        #---------------------------
        # Can canal_file be found ?
        #---------------------------
        FOUND = tf_utils.file_exists( self.canal_file )
        if not(FOUND):
            self.use_canals = False
            return

        #------------------------
        # Open the "canal_file"
        #------------------------
        file_unit = open(self.canal_file, 'r')

        #------------------------
        # Read number of canals
        #------------------------
        n_canals      = cfg.read_value(file_unit, dtype='Int32')
        self.n_canals = n_canals

        #--------------------
        # Initialize arrays
        #--------------------
        self.canal_in_IDs      = np.zeros([n_canals], dtype='Int32')
        self.canal_out_IDs     = np.zeros([n_canals], dtype='Int32')
        self.canal_Q_fractions = np.zeros([n_canals], dtype='Float64')
        self.canal_times       = np.zeros([n_canals], dtype='Float64')
        
        #----------------------------------
        # Read information for each canal
        #----------------------------------
        for k in range(n_canals):
            canal_in_ID  = cfg.read_value(file_unit, dtype='Int32')
            canal_out_ID = cfg.read_value(file_unit, dtype='Int32')
            Q_fraction   = cfg.read_value(file_unit, dtype='Float64')
            travel_time  = cfg.read_value(file_unit, dtype='Float64')
            #----------------------------------------------------------
            self.canal_in_IDs[k]      = canal_in_ID
            self.canal_out_IDs[k]     = canal_out_ID
            self.canal_Q_fractions[k] = Q_fraction
            self.canal_times[k]       = travel_time

        #--------------------------------------------------------
        # Compute "nt_canals", which is the number of timesteps
        # it takes for flow to travel from end to end.
        #--------------------------------------------------------
        # This depends on "self.dt", which is now read from the
        # Diversion component CFG file.  ## (9/22/14)
        #--------------------------------------------------------
        self.nt_canals = np.ceil(self.canal_times / self.dt)
        
        #-----------------------
        # Close the input file
        #-----------------------
        file_unit.close()

        #-----------------------------------------------------
        # Compute xy coordinates for canal entrance and exit
        #-----------------------------------------------------
        canal_in_rows    = (self.canal_in_IDs / self.nx)
        canal_in_cols    = (self.canal_in_IDs % self.nx)
        self.canals_in_x = (canal_in_cols * self.dx)
        self.canals_in_y = (canal_in_rows * self.dy)       
        #-----------------------------------------------------
        canal_out_rows    = (self.canal_out_IDs / self.nx)
        canal_out_cols    = (self.canal_out_IDs % self.nx)
        self.canals_out_x = (canal_out_cols * self.dx)
        self.canals_out_y = (canal_out_rows * self.dy)         
                        
        #-----------------------------------------------------
        # Create a 2D array to store the discharge values as
        # they are moving toward downstream end of canal.
        #-----------------------------------------------------
        # update_canals() will "roll" this array downstream
        # by one array element each time step
        #-----------------------------------------------------
        nt_max       = np.int(self.nt_canals.max())
        nt_min       = np.int(self.nt_canals.min())
        self.canal_Q = np.zeros([n_canals, nt_max], dtype='Float64')
        self.nt_max  = nt_max
        print('Diversions component: Min steps per canal =', nt_min)
        print('Diversions component: Max steps per canal =', nt_max)
        
    #   read_canal_data()
    #--------------------------------------------------------------------------
    def update_sources(self):

        #---------------------------------------------------------
        # Notes: This function avoids loops in favor of array
        #        operations to increase speed.
        #---------------------------------------------------------
        #        The number of Q-values for each source ID are
        #        stored as "self.nt_sources".  However, for any
        #        given source ID, the Q-values beyond that index
        #        are set to zero.
        #---------------------------------------------------------        
        if not(self.use_sources): return
        
        #-------------------------------------------  
        # Update discharge, Q, for every source ID
        #-------------------------------------------
        if (self.time_index < self.nt_max_sources):
            self.Q_sources[:] = self.Q_sources_all[ :, self.time_index ]
        else:
            self.Q_sources[:] = np.zeros(self.n_sources)

        #------------------------------------------------------------
        # Update flow volumes, vol, in CHANNELS component (2/17/10)
        #------------------------------------------------------------
        # If a grid cell contains a "source", then an additional Q
        # will flow *into* that grid cell and increase flow volume.
        #------------------------------------------------------------ 
        ## vols = self.cp.get_values_in_grid_double( 'vol', self.source_IDs )
        ## vols += (Q_sources * self.dt)
        ## self.cp.set_values_in_grid_double( 'vol', self.source_IDs, vols )
        
        #--------------
        # For testing
        #--------------
        ## print 'Finished with update_sources() in Diversions.'
        ## print 'Q_sources ='
        ## print Q_sources
        
    #   update_sources()
    #--------------------------------------------------------------------------
    def update_sinks(self):

        #-------------------------------------------------------
        # Notes: This function avoids loops in favor of array
        #        operations to increase speed.
        #-------------------------------------------------------
        #        The number of Q-values for each sink ID are
        #        stored as "self.nt_sinks".  However, for any
        #        given sink ID, the Q-values beyond that index
        #        are set to zero.
        #-------------------------------------------------------
        # NB!    This changes Q grid, so must be called before
        #        cp.update_flow_volume() uses the Q grid.
        #-------------------------------------------------------        
        if not(self.use_sinks): return

        #-----------------------------------------  
        # Update discharge, Q, for every sink ID
        #-----------------------------------------
        # Make sure sink cannot produce negative
        # discharge values.
        #-----------------------------------------
        if (self.time_index < self.nt_max_sinks):
            Q_sinks[:] = self.Q_sinks_all[ :, self.time_index ]
        else:
            Q_sinks[:] = np.zeros(self.n_sinks)

        #--------------------------------------------------------
        # Update discharges, Q, in CHANNELS component (2/17/10)
        #--------------------------------------------------------
        # NB! If a grid cell contains a "sink" or the upstream
        #     end of a "canal", then the discharge *leaving*
        #     that grid cell must be reduced accordingly.
        #     This changes the Q grid and must be done before
        #     using the Q grid for further calculations below.
        #--------------------------------------------------------
        ## Q_vals = self.cp.get_values_in_grid_double( 'Q', self.sink_IDs )
        ## Q_vals -= Q_sinks
        ## self.cp.set_values_in_grid_double( 'Q', self.sink_IDs, Q_vals )

        #------------------------------------------------------------
        # Update flow volumes, vol, in CHANNELS component (2/17/10)
        # NB!  We MUST update "vol" also and not just "Q".
        #------------------------------------------------------------
        ## vols = self.cp.get_values_in_grid_double( 'vol', self.sink_IDs )
        ## vols -= (Q_sinks * self.dt)
        ## self.cp.set_values_in_grid_double( 'vol', self.sink_IDs, vols )
        
    #   update_sinks()
    #--------------------------------------------------------------------------
    def update_canals(self):

        #----------------------------------------------------------
        # Notes: Before 2/1/10, TopoFlow would update the channel
        #        component before the diversion component.  Now
        #        cp.update_discharge() calls dp.update() itself.
        #
        #        (2/16/10) Tested for single canal and seems to
        #        work as intended.
        #
        # NB!    Flow volumes are incremented (by the function
        #        update_flow_volume(), but discharges are
        #        always recomputed from d, v and channel geom.
        #        So changes to cp.Q by calling dp.update()
        #        after cp.update() would be lost.
        #
        #        cp.update() computes variables in this order:
        #           update_R()
        #           update_discharge()    (using d and v)
        #           update_flow_volume()  (using Q and R)
        #           update_flow_depth()   (using vol)
        #           update_velocity()
        #----------------------------------------------------------
        # Notes: This function avoids loops in favor of array
        #        operations to increase speed.
        #----------------------------------------------------------
        # Notes: In this version, the Channels component uses
        #        canal_Q_fractions, canal_in_IDs and its own
        #        Q grid to compute Q_canals_in.  It then sets
        #        this into the Diversion component.
        #----------------------------------------------------------
        # print '#### Starting update_canals()...'
        if not(self.use_canals): return
       
        #---------------------------------------------------
        # Update discharges, Q, at upstream ends of canals
        # in CHANNELS component (2/17/10)
        #---------------------------------------------------
        # print '#### update_canals(), Q_vals block...'
#         Q_vals = self.cp.get_values_in_grid_double( 'Q', self.canal_in_IDs )
#         Q_canals_in = (self.canal_Q_fractions * Q_vals)
#         Q_vals -= Q_canals_in
#         self.cp.set_values_in_grid_double( 'Q', self.canal_in_IDs, Q_vals )

        #-------------------------------------------------------
        # Update flow volumes, vol, at upstream ends of canals
        # in CHANNELS component (2/17/10)
        # NB!  We MUST update "vol" also and not just "Q".
        #-------------------------------------------------------
        # print '#### update_canals(), vols block...'
#         vols = self.cp.get_values_in_grid_double( 'vol', self.canal_in_IDs )
#         vols -= (Q_canals_in * self.dt)
#         self.cp.set_values_in_grid_double( 'vol', self.canal_in_IDs, vols )

        #----------------------------------------------
        # Add specified discharge (Q_in) to upstream
        # end of each canal (at array index 0)
        ######################################################
        # Diversions component now gets Q_canals_in from the
        # Channels component as a requested input_var.
        ######################################################
        # print '#### update_canals(), canal_Q block...'
        ## self.canal_Q[:, 0] = Q_canals_in        
        self.canal_Q[:, 0] = self.Q_canals_in   # (from Channels)
        
        #------------------------------------------------
        # Get "Q" at downstream end of each canal.
        # It will be zero until flow has had time to
        # travel the distance.
        #------------------------------------------------
        # NB!  Each canal can have a different travel
        # time and therefore its own "nt_canal" value.
        #------------------------------------------------
        # NB! Q_canals_out will be retrieved by the
        #     Channels component.
        #-------------------------------------------------------
        # Note that canal_Q is defined as:
        # canal_Q = zeros([n_canals, nt_max], dtype='Float64')
        #-------------------------------------------------------
        # print '#### update_canals(), for loop...'
        nc  = self.n_canals
        ## Q_canals_out = np.empty(nc, dtype='Float32')
        self.Q_canals_out[:] = np.empty(nc, dtype='Float32')
        for k in range(nc):
            nt_k  = self.nt_canals[k]
            self.Q_canals_out[k] = self.canal_Q[k, nt_k-1]
            ## Q_canals_out[k] = self.canal_Q[k, nt_k-1]
            ## Q_canals_out[k] = self.canal_Q[:, nt_k-1]
        ## self.canal_Q[:, nt_k:] = 0.0  # (not necessary)

        #---------------------------------------------------------
        # Update flow volumes, vol, at downstream ends of canals
        # in CHANNELS component (2/17/10)
        # NB!  We MUST update "vol" also and not just "Q".
        #---------------------------------------------------------
#         vols = self.cp.get_values_in_grid_double( 'vol', self.canal_out_IDs )
#         vols += (Q_canals_out * self.dt)
#         self.cp.set_values_in_grid_double( 'vol', self.canal_out_IDs, vols )

        #--------------
        # For testing
        #--------------
##        print 'self.canal_Q ='
##        print self.canal_Q
##        print 'self.Q_canals_out ='
##        print self.Q_canals_out
##        print ' '
        
        #----------------------------------------------------
        # "Roll" the canal_Q array downstream (along its
        # 2nd index, with index = 1) by one array element
        # in each time step. "canal_Q" starts with zeros;
        # i.e.  canal_Q = zeros([n_canals, n_steps])
        #----------------------------------------------------
        # In next call to update_canals(), we'll replace
        # the first Q-value in each canal.
        #----------------------------------------------------
        # print '#### update_canals(), roll...'
        self.canal_Q = np.roll( self.canal_Q, 1, axis=1 )

        # print '#### Exiting update_canals()...'
    
    #   update_canals()
    #--------------------------------------------------------------------------


    

