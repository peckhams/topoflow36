
## Copyright (c) 2010-2013, Scott D. Peckham
## July 2010 (created)

#-----------------------------------------------------------------------
#  Notes:  This file defines a "base class" for HIS Data components
#          that retrieve HIS data from a web service using suds.
#-----------------------------------------------------------------------
#
#  class HIS_component
#
#      get_component_name()
#      get_attribute()
#      get_input_var_names() 
#      get_output_var_names()
#      get_var_name()
#      get_var_units()
#      ------------------------
#      set_constants()
#      import_suds()             ##########
#      initialize()
#      update()
#      finalize()
#      set_computed_input_vars()
#      ---------------------------
#      check_input_types()
#      initialize_computed_vars()
#      -------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      ------------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()

#-----------------------------------------------------------------------

import numpy as np

import os
import socket          # (to get hostname)
import suds.client
import sys
import urllib.request, urllib.parse, urllib.error
# import urllib2

from topoflow.utils import BMI_base
from topoflow.utils import model_input
from topoflow.utils import model_output

#-----------------------------------------------------------------------
class HIS_component(BMI_base.BMI_component):

    _att_map = {
        'model_name':         'Data_HIS',
        'version':            '3.1',
        'author_name':        'Scott D. Peckham',
        'grid_type':          'none',
        'time_step_type':     'fixed',
        'step_method':        'none',  # (or 'explicit' ??)
        #---------------------------------------------------
        'comp_name':          'DataHIS',       
        'model_family':       'TopoFlow',
        'cfg_template_file':  'Data_HIS.cfg.in',
        'cfg_extension':      '_data_his.cfg',
        'cmt_var_prefix':     '/DataHIS/Input/Var/',
        'gui_xml_file':       '/home/csdms/cca/topoflow/3.1/src/share/cmt/gui/Data_HIS.xml',
        'dialog_title':       'Data: HIS Query Parameters',
        'time_units':         'seconds' }

    #--------------------------------------------------------
    # These must be determined dynamically, based on query.
    # Start with empty lists and dictionaries.
    #--------------------------------------------------------
    _input_var_names = []

    ### Possible examples:  ['Q', 'T_air', 'P']
    _output_var_names = []

    _var_name_map = {}

    _var_units_map = {}

    #------------------------------------------------    
    # Return NumPy string arrays vs. Python lists ?
    #------------------------------------------------
    ## _input_var_names  = np.array( _input_var_names )
    ## _output_var_names = np.array( _output_var_names )

    #-------------------------------------------------------------------
    def get_component_name(self):
  
        return 'TopoFlow_Data_HIS'

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
##    def get_var_type(self, long_var_name):
##
##        #---------------------------------------
##        # So far, all vars have type "double",
##        # but use the one in BMI_base instead.
##        #---------------------------------------
##        return 'float64'
##    
##    #   get_var_type()    
    #-------------------------------------------------------------------
    def set_constants(self):       

        #------------------------
        # Define some constants
        #--------------------------------------------------------------
        # Note: When running on one of our cluster's compute nodes,
        #       vs. on the head node (beach.colorado.edu), we get
        #       the following error message:
        #
        #       "Unexpected Python exception:
        #         <urlopen error [Errno -3]
        #         Temporary failure in name resolution>
        #--------------------------------------------------------------       
        URL = "http://water.sdsc.edu/hiscentral/webservices/hiscentral.asmx?WSDL"  # (same msg)
        self.HIS_Central_URL = URL

        #----------------------------------------------
        # URL for National Elevation Dataset (NED) ??
        #----------------------------------------------
        # URL2 = "http://extract.cr.usgs.gov/axis2/services/DownloadService?wsdl"
        # URL3 = "http://gisdata.usgs.gov/XMLWebServices2/Elevation_Service.asmx"
        # info ="http://imselev.cr.usgs.gov/WMS_Capabilities/USGS_EDC_Elev_NED/capabilities_1_3_0.xml"

        #----------------------------------------------------------
        # Create a Python dictionary with proxies.  The compute
        # nodes on our cluster don't have direct access so need
        # to specify the head node as proxy (beach.colorado.edu).
        #----------------------------------------------------------
        # Well-known TCP/IP port numbers:
        #     HTTP:        80   (or sometimes 8080, 7080, 9080)
        #     HTTPS:       443
        #     FTP-Control: 21
        #     FTP-Data:    20
        #     SSH:         22
        #     Telnet:      23
        # We can get this number on localhost with:
        #     socket.getservbyname('http', 'tcp')  # (service by name)
        # Or the other way around:
        #     socket.getservbyport(8080)  # (-> 'webcache' on beach)
        #--------------------------------------------------------------
        # Not used now.  Proxy is now set as environment variable
        # "http_proxy="http://128.138.77.51", but still needs to
        # be set for everyone.
        #--------------------------------------------------------------
##        self.proxies = {'http':'http://128.138.77.51',
##                        'https':'https://128.138.77.51',
##                        'ftp':'ftp://128.138.77.51'}
        self.proxies = {'http':'http://128.138.77.51:80',
                        'https':'https://128.138.77.51:443',
                        'ftp':'ftp://128.138.77.51:21'}

        #------------------------------------------
        # Get the hostname (head or compute node)
        #------------------------------------------
        self.hostname = socket.gethostname()
        
        #-----------------------------
        # Option to run a diagnostic
        #-----------------------------
##        if (self.DEBUG):
##            print 'Trying to open HIS Central URL directly...'
##            try:
##                HIS_URL = urllib.urlopen( self.HIS_Central_URL,
##                                          proxies=self.proxies )
##                # HIS_URL = urllib2.urlopen( self.HIS_Central_URL,
##                #                            proxies=self.proxies )
##                print HIS_URL
##                print 'SUCCESS: Opened URL with urllib.urlopen().'
##            except:
##                print 'FAILURE: Could not open URL with urllib.urlopen().'
##            print ' '
        
    #   set_constants()
    #-------------------------------------------------------------------
    def import_suds(self):

        python_version = sys.version[:5]    # (e.g. 2.6.2)
        # python_version = sys.version[:3]  # (e.g. 2.6)
         
        try:
            import suds
            print('Imported suds version: ' + suds.__version__)
            print('Python version:        ' + python_version)
            print('Hostname or node:      ' + self.hostname)
            print(' ')
            return True
        except:

            print(' ')
            print('SORRY, Cannot access HIS web service because')
            print('the "suds.client" package cannot be imported.')
            print(' ')
##            if (python_version != '2.6'):
##                print 'Note that "suds" is only installed for'
##                print 'Python version 2.6 on "beach".'
##                print 'The current Python version is:', python_version
##                print ' '
            return False

    #   import_suds()   
    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        if not(SILENT):
            print(' ')
            print('Data-HIS component: Initializing...')
            
        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()
        self.initialize_config_vars()
        ## self.read_grid_info()    # NOW IN initialize_config_vars()
        # self.initialize_basin_vars()
        self.initialize_time_vars()
        self.dt = 1.0      # (Needed for final report.)
        self.n_steps = 1   # (Need for BMI_base.update_time().)
    
        #----------------------------------
        # Has component been turned off ?
        #----------------------------------
        if (self.comp_status == 'Disabled'):
            if not(SILENT):
                print('HIS Data component: Disabled in CFG file.')
            self.T_air  = np.float64(0)   ####################################
            self.P      = np.float64(0)
            self.DONE   = True
            self.status = 'initialized'  # (OpenMI 2.0 convention)
            return
        else:
            self.DONE = False  ########
            
        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        # self.open_input_files()
        # self.read_input_files()
        # self.check_input_types()
        # self.initialize_computed_vars()

        #--------------------------------------
        # Import suds and print version, etc.
        #--------------------------------------
        result = self.import_suds()
        #-------------------------------        
        # These are not used currently
        #-------------------------------
        # import xml.etree.ElementTree as ET
        # import os, fnmatch
        # import datetime as DT
        
        self.open_output_files()
        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
    #   initialize()
    #-------------------------------------------------------------------
    ## def update(self, dt=-1.0, time_seconds=None):
    def update(self, dt=-1.0):

        if (self.comp_status == 'Disabled'): return
        self.status = 'updating'  # (OpenMI 2.0 convention)

        #----------------------------------------
        # Retrieve data values from web service
        #----------------------------------------
        xmin   = str(self.west_edge_lon)
        xmax   = str(self.east_edge_lon)
        ymax   = str(self.north_edge_lat)
        ymin   = str(self.south_edge_lat)
        networkIDs = ''
        
        print('Attempting to retrieve data from HIS Central...')
        print('     HIS_Central_URL =', self.HIS_Central_URL)
        print('     Data keyword    =', self.data_keyword)
##        print '     xmin, xmax      =', xmin, ',', xmax
##        print '     ymin, ymax      =', ymin, ',', ymax
##        print '     start_date      =', self.start_date
##        print '     stop_date       =', self.stop_date

        #----------------------------------------------------------------------
        # (7/19/10)  This works on both head and compute nodes now, after
        # reconfiguring the system (UnixOps) and adding the environment
        # variable "http_proxy" to ".bash_profile".
        #----------------------------------------------------------------------      
        client = suds.client.Client( self.HIS_Central_URL )

        # proxy = self.proxies    ######        
        #----------------------------------------------------------------------        
        # client = suds.client.Client( self.HIS_Central_URL )
        # ERROR: Unexpected Python exception:  <urlopen error [Errno -3]
        #        Temporary failure in name resolution>
        #----------------------------------------------------------------------
        # client = suds.client.Client( self.HIS_Central_URL, proxy=proxy )
        # ERROR: Unexpected Python exception: HTTP Error 404: Not Found
        #----------------------------------------------------------------------
        # client = suds.client.Client( self.HIS_Central_URL, service=proxy )
        # ERROR: Unexpected Python exception:
        #        "service" must be: (<type 'int'>, <type 'basestring'>
        #----------------------------------------------------------------------
        # from suds.transport.http import HttpAuthenticated
        # tran = HttpAuthenticated(proxy=proxy)
        # client = suds.client.Client( self.HIS_Central_URL, transport=tran )
        # ERROR: Unexpected Python exception: HTTP Error 404: Not Found
        #----------------------------------------------------------------------
        # from suds.transport.http import HttpTransport
        # tran   = HttpTransport( proxy=proxy )
        # client = suds.client.Client( self.HIS_Central_URL, transport=tran )
        # ERROR: Unexpected Python exception: HTTP Error 404: Not Found
        #----------------------------------------------------------------------
        # from suds.transport import HttpTransport
        # tran   = HttpTransport( proxy=proxy )
        # client = suds.client.Client( self.HIS_Central_URL, transport=tran )
        # ERROR: Unexpected Python exception: cannot import name HttpTransport
        #----------------------------------------------------------------------
        # client = suds.client.Client( url='' )
        # client.set_options( url=self.HIS_Central_URL, proxy=proxy )
        # ERROR: Unexpected Python exception: unknown url type: 
        #----------------------------------------------------------------------
        # client = suds.client.Client( url='http' )
        # client.set_options( url=self.HIS_Central_URL, proxy=proxy )
        # ERROR: Unexpected Python exception: unknown url type: http
        #----------------------------------------------------------------------
        # client = suds.client.Client( url='http://' )
        # client.set_options( url=self.HIS_Central_URL, proxy=proxy )
        # ERROR: Unexpected Python exception: <urlopen error no host given> 
        #----------------------------------------------------------------------
        # client = suds.client.Client( url=self.HIS_Central_URL )
        # client.set_options( url=self.HIS_Central_URL, proxy=proxy )
        # ERROR: Unexpected Python exception: <urlopen error [Errno -3]
        #        Temporary failure in name resolution>
        #----------------------------------------------------------------------
        # client = suds.client.Client( url=self.HIS_Central_URL )
        # client.set_options( url=self.HIS_Central_URL, proxy=proxy )
        # ERROR: Unexpected Python exception: <urlopen error [Errno -3]
        #        Temporary failure in name resolution>
        #----------------------------------------------------------------------
        # client = suds.client.Client( url="http://localhost" )
        # client.set_options( url=self.HIS_Central_URL, proxy=proxy )
        # ERROR: Unexpected Python exception: <urlopen error [Errno 111]
        #        Connection refused>
        #----------------------------------------------------------------------
        # client = suds.client.Client( url="http://127.0.0.1" )
        # client.set_options( url=self.HIS_Central_URL, proxy=proxy )
        # ERROR: Unexpected Python exception: <urlopen error [Errno 111]
        #        Connection refused>
        #----------------------------------------------------------------------
        # from suds.transport.http import HttpTransport
        # tran   = HttpTransport()
        # import urllib2
        # opener     = urllib2.build_opener( proxy )
        # tran.urlopener = opener
        # client = suds.client.Client( self.HIS_Central_URL, transport=tran )
        # ERROR: Unexpected Python exception:
        #        expected BaseHandler instance, got <type 'dict'>
        #----------------------------------------------------------------------
        # from suds.transport.http import HttpTransport
        # tran   = HttpTransport()
        # import urllib2
        # proxy_han  = urllib2.ProxyHandler( proxies=proxy )
        # opener     = urllib2.build_opener( proxy_han )
        # tran.urlopener = opener
        # client = suds.client.Client( self.HIS_Central_URL, transport=tran )
        # ERROR: Unexpected Python exception:
        #        expected BaseHandler instance, got <type 'dict'>
        #----------------------------------------------------------------------
        # from suds.transport.http import HttpTransport
        # tran   = HttpTransport()
        # import urllib2
        # proxy_han  = urllib2.ProxyHandler( self.proxies )
        # opener     = urllib2.build_opener( proxy_han )
        # tran.urlopener = opener
        # client = suds.client.Client( self.HIS_Central_URL, transport=tran )
        # ERROR: Unexpected Python exception:
        #        expected BaseHandler instance, got <type 'dict'>
        #----------------------------------------------------------------------
        # print client  #########
        #######################
        ## print 'Calling client.add_prefix() method...'
        ## client.add_prefix( 'hiscentral.asmx' )  #####
        ## print 'Created client object.  ########'
        response = client.service.GetSeriesCatalogForBox2(xmin, xmax, ymin, ymax, \
                                                          self.data_keyword, networkIDs, \
                                                          self.start_date, self.stop_date)
        rep_len = len(response)
        if (rep_len == 0):
            print('-------------------------------------------')
            print('SORRY, Your query returned no matches.')
            print('       You may want to try:')
            print('       (1) Changing the keyword.')
            print('       (2) Expanding the bounding box.')
            print('       (3) Expanding the range of dates.')
            print('-------------------------------------------')
            print(' ')
            self.DONE = True
            return
        # print 'len(response) =', rep_len

        series_array = response[0]
        n_series     = len( series_array )
        print('Number of time series retrieved =', n_series)
        print('Max allowed number of series    =', self.max_n_series) 
        print(' ')
   
        #------------------------------------------------ 
        # Return if n_series > max_n_series. (10/25/11)
        #------------------------------------------------ 
        if (n_series > self.max_n_series):
            print('###############################################')
            print(' ABORTING: n_series > max_n_series.')
            print('           n_series =', n_series)
            print('       max_n_series =', self.max_n_series)
            print('###############################################')
            print(' ')
            return
    
        #-----------------------------------
        # Check the data for nodata/NaN ??
        #-----------------------------------   

        #------------------------------------------
        # Read next infil vars from input files ?
        #------------------------------------------
        # if (self.time_index > 0):
        #     self.read_input_files()

        #-------------------------------------------
        # Extract data values for each time series
        #-------------------------------------------
        ## print 'Writing time series to local files...'
        series_num = 0  
        # print 'Processing array of time series...'
        for series in series_array:

            series_num += 1
            print('=========================================================')
            print('Working on series', series_num, 'of', n_series)
            try:
                client     = suds.client.Client( series.ServURL )
                values_obj = client.service.GetValuesObject(series.location, series.VarCode, \
                                                            self.start_date, self.stop_date)
                values   = values_obj.timeSeries.values
                n_values = values._count

                SKIP = False
                if (n_values < 1):
                    print('Skipping series: n_values =', n_values, ' ( < 1 )')
                    print(' ')
                    SKIP = True
                if (n_values > self.max_n_values):
                    print('Skipping series: n_values =', n_values)
                    print('         but max_n_values =', self.max_n_values)
                    print(' ')
                    SKIP = True
                #------------------------------------------------------
                # Get actual time series values and date-time strings
                #------------------------------------------------------
                if not(SKIP):
                    #-------------------------------------
                    # Is there a faster way than this ??
                    #-------------------------------------
                    vals = np.asarray([np.float(val.value) for val in values.value])
                    dts  = np.asarray([val._dateTime for val in values.value])
                    #-----------------------------------------------
                    # Write the data values to a text file here ??
                    #-----------------------------------------------
                    ## if (self.DEBUG):
                    if (True):
                        print('series ID         =', series_num)
                        print('size(series)      =', np.size(vals))
                        print('min(series)       =', vals.min())
                        print('max(series)       =', vals.max()) 
                        print('series.Sitename   =', series.Sitename)
                        print('series.location   =', series.location)
                        print('series.VarName    =', series.VarName)
                        print('series.VarCode    =', series.VarCode)
                        print('series.ValueCount =', series.ValueCount)
                        print('series.datatype   =', series.datatype)
                        print('series.valuetype  =', series.valuetype)
                        print('series.timeunits  =', series.timeunits)
                        # print 'len(series)       =', len(series)   # (always 18)
                        # print 'type(values_obj)  =', type(values_obj)   # ("instance")
                        # print 'type(values)      =', type(values)    # ("instance")
                        print(' ')
            except:
                print('ERROR: Skipping time series.')
                print(' ')

        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # Components use own self.time_sec by default.
        #-----------------------------------------------
        #self.write_output_files()
        ## self.write_output_files( time_seconds )

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        self.update_time( dt )
        self.status = 'updated'  # (OpenMI 2.0 convention)

        #-----------------------------------
        # Only allow one call to update ??
        #-----------------------------------
        self.DONE = True
        
    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalizing'  # (OpenMI 2.0 convention)
        ## self.close_input_files()   ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'  # (OpenMI 2.0 convention)

        # self.print_final_report(comp_name='HIS Data component')
     
    #   finalize()
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        #-----------------------------------------------
        # Convert start_month & stop_month from string
        # to integer.  January should be 1.
        #-----------------------------------------------
        month_list = ['January', 'February', 'March', 'April',
                      'May', 'June', 'July', 'August', 'September',
                      'October', 'November', 'December']
        self.start_month = month_list.index( self.start_month ) + 1
        self.stop_month  = month_list.index( self.stop_month )  + 1

        #-----------------------------------------------
        # Construct the query start date and stop date
        #-----------------------------------------------
        self.start_date = str(self.start_year) + '-' + \
                          str(self.start_month).zfill(2) + '-' + \
                          str(self.start_day).zfill(2)
        self.stop_date  = str(self.stop_year) + '-' + \
                          str(self.stop_month).zfill(2) + '-' + \
                          str(self.stop_day).zfill(2)

        #--------------------------------------------------------
        # Define these here, so all components can use the same
        # output file functions, like "open_output_files()".
        #--------------------------------------------------------
##        self.SAVE_Q0_GRIDS  = False
##        self.SAVE_ZW_GRIDS  = False

        #--------------------------------------------
        # Can maybe remove this when GUI info file
        # has the "save cubes" part put back in.
        #--------------------------------------------
##        self.q_cs_file = ''
##        self.p_cs_file = ''
        
        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
##        self.save_grid_dt    = np.maximum(self.save_grid_dt,    self.dt)
##        self.save_pixels_dt  = np.maximum(self.save_pixels_dt,  self.dt)
##        self.save_profile_dt = np.maximum(self.save_profile_dt, self.dt)
##        self.save_cube_dt    = np.maximum(self.save_cube_dt,    self.dt)
        
    #   set_computed_input_vars()      
    #-------------------------------------------------------------------
    def check_input_types(self):

        #------------------------------------------------------
        # Notes: Usually this will be overridden by a given
        #        method of computing ET.  But this one should
        #        work for Green-Ampt and Smith-Parlange.
        #------------------------------------------------------
##        are_scalars = np.array([
##                         self.mp.is_scalar('P'),
##                         self.sp.is_scalar('SM'),
##                         self.gp.is_scalar('h_table'),
##                         #------------------------------
##                         self.is_scalar('Ks'),
##                         self.is_scalar('Ki'),
##                         self.is_scalar('qs'),
##                         self.is_scalar('qi'),
##                         self.is_scalar('G'),
##                         self.is_scalar('gam')  ])
##
##        self.ALL_SCALARS = np.all(are_scalars)

        self.ALL_SCALARS = True   # (HIS only works with time series)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        pass
         
    #   initialize_computed_vars()
    #-------------------------------------------------------------------  
    def open_input_files(self):

        pass
    
##        self.Ks_unit  = []  # (empty lists to hold file objects)
##
##        for k in xrange(self.n_layers):
##            self.Ks_unit.append(  model_input.open_file(self.Ks_type[k],  self.Ks_file[k]) )
        
    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        if (self.DEBUG):
            print('Calling read_input_files()...')

        pass
    
##        rti = self.rti
##
##        #-------------------------------------------------------
##        # All grids are assumed to have data type of Float32.
##        #-------------------------------------------------------
##        # This method works for Green-Ampt and Smith-Parlange
##        # but must be overridden for Richards 1D.
##        #-------------------------------------------------------
##        # NB! Green-Ampt and Smith-Parlange currently only
##        #     support ONE layer (n_layers == 1).
##        #------------------------------------------------------- 
##        for k in xrange(self.n_layers):
##            Ks = model_input.read_next(self.Ks_unit[k], self.Ks_type[k], rti)
##            if (Ks is not None): self.Ks[k] = Ks
          
    #   read_input_files()       
    #-------------------------------------------------------------------  
    def close_input_files(self):

        pass
    
##        for k in xrange(self.n_layers):
##            if (self.Ks_type[k]  != 'Scalar'): self.Ks_unit[k].close()
          
    #   close_input_files()
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        pass
    
        #-------------------------------------------------
        # Notes:  Append out_directory to outfile names.
        #-------------------------------------------------
##        self.v0_gs_file = (self.out_directory + self.v0_gs_file)
##        #-------------------------------------------------------------
##        self.v0_ts_file = (self.out_directory + self.v0_ts_file)

    #   update_outfile_names()   
    #-------------------------------------------------------------------  
    def open_output_files(self):

        #-------------------------------------------------
        # Notes:  v0 = infiltration rate at surface
        #-------------------------------------------------
        model_output.check_netcdf()
        self.update_outfile_names()
        
##        #--------------------------------------
##        # Open new files to write grid stacks
##        #--------------------------------------
##        if (self.SAVE_V0_GRIDS):
##            model_output.open_new_gs_file( self, self.v0_gs_file, self.rti,
##                                           var_name='v0',
##                                           long_name='infiltration_rate_at_surface',
##                                           units_name='m/s')
##                                
##        #--------------------------------------
##        # Open new files to write time series
##        #--------------------------------------
##        IDs = self.outlet_IDs
##        if (self.SAVE_V0_PIXELS):
##            model_output.open_new_ts_file( self, self.v0_ts_file, IDs,
##                                           var_name='v0',
##                                           long_name='infiltration_rate_at_surface',
##                                           units_name='m/s')
            
    #   open_output_files()
    #-------------------------------------------------------------------
    def write_output_files(self, time_seconds=None):

        #---------------------------------------------------------
        # Notes:  This function was written to use only model
        #         time (maybe from a caller) in seconds, and
        #         the save_grid_dt and save_pixels_dt parameters
        #         read by read_cfg_file().
        #
        #         read_cfg_file() makes sure that all of
        #         the "save_dts" are larger than or equal to the
        #         process dt.
        #---------------------------------------------------------
        if (self.DEBUG):
            print('Calling write_output_files()...')

        #-----------------------------------------
        # Allows time to be passed from a caller
        #-----------------------------------------
        if (time_seconds is None):
            time_seconds = self.time_sec
        model_time = int(time_seconds)
        
        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
        if (model_time % int(self.save_grid_dt) == 0):
            self.save_grids()
            
        if (model_time % int(self.save_pixels_dt) == 0):
            self.save_pixel_values()
      
    #   write_output_files()           
    #-------------------------------------------------------------------
    def close_output_files(self):

        pass
    
##        if (self.SAVE_V0_GRIDS): model_output.close_gs_file( self, 'v0')   
##
##        if (self.SAVE_V0_PIXELS): model_output.close_ts_file( self, 'v0')
        
    #   close_output_files()
    #-------------------------------------------------------------------  
    def save_grids(self):
        
        #-----------------------------------
        # Save grid stack to a netCDF file
        #---------------------------------------------
        # Note that add_grid() methods will convert
        # var from scalar to grid now, if necessary.
        #---------------------------------------------
        if (self.DEBUG):
            print('Calling save_grids()...')
            
##        if (self.SAVE_V0_GRIDS):
##            model_output.add_grid( self, self.IN, 'v0', self.time_min )

    #   save_grids()                             
    #-------------------------------------------------------------------  
    def save_pixel_values(self):

        if (self.DEBUG):
            print('Calling save_pixel_values()...')
            
##        IDs  = self.outlet_IDs
##        time = self.time_min   ########
##        
##        #--------------------------------------------
##        # Save a subsequence of IN var pixel values
##        #--------------------------------------------  
##        if (self.SAVE_V0_PIXELS):
##            model_output.add_values_at_IDs( self, time, self.IN, 'v0', IDs )

    #   save_pixel_values()
    #-------------------------------------------------------------------

