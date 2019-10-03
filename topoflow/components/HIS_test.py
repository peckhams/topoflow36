
##import xml.etree.ElementTree as ET
##import os, fnmatch
##import pylab
##import numpy as NY
##import datetime as DT

import numpy
import socket
import suds.client
import suds.transport
import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse

#-------------------------------------------------------------------------------
def HIS_test():

    #--------------------------------------------------------------
    # Notes: The year 2000 was a leap year, and we're asking for
    #        one year's worth of data.  The timeunits are "days",
    #        so size(vals) = 366 when len(values) = 3.
    #        What does it mean when len(values) = 1 ??
    #
    #        Note that len(series) = 18, which seems to be just
    #        the number of fields or elements in a series object.
    #--------------------------------------------------------------
    # URL to the HIS Central API
    #-----------------------------
    HIS_Central_URL   = 'http://water.sdsc.edu/hiscentral/webservices/hiscentral.asmx?WSDL'
    HIS_location      = 'http://water.sdsc.edu/'
    # HIS_Central_URL   = 'http://water.sdsc.edu/hiscentral/webservices/hiscentral.asmx?wsdl'
        # (Changing case of WSDL as shown is OK for urllib2 test.)
    # HIS_Central_URL   = 'http://water.sdsc.edu/hiscentral/webservices/hiscentral.asmx?wsdl@128.138.77.51'
        # (Adding proxy IP address as shown is OK for urllib2 test.)
    USGS_Download_URL = 'http://extract.cr.usgs.gov/axis2/services/DownloadService?wsdl'
    USGS_location     = 'http://extract.cr.usgs.gov/'

    # web_service_URL = HIS_location
    web_service_URL = HIS_Central_URL
    # web_service_URL = USGS_Download_URL
    print('web_service_URL =', web_service_URL)
    
    #----------------------------------------------
    # NB!  Without the WSDL part, it doesn't work
    #----------------------------------------------
    # HIS_Central_URL = "http://water.sdsc.edu/hiscentral/webservices/hiscentral.asmx"

    #----------------------------------------------------------
    # Create a Python dictionary with proxies.  The compute
    # nodes on our cluster don't have direct access so need
    # to specify the head node as proxy (beach.colorado.edu).
    #----------------------------------------------------------
    # Well-known TCP/IP port numbers:
    #     HTTP:        80    (or sometimes 8080, 7080, 9080)
    #     HTTPS:       443
    #     FTP-Control: 21
    #     FTP-Data:    20
    #     SSH:         22
    #     Telnet:      23
    # When can get this number on localhost by:
    #     socket.getservbyname('http', 'tcp')
    # Or the other way around:
    #     socket.getservbyport(8080)  # (-> 'webcache' on beach)
    #--------------------------------------------------------------
    proxies = {'http':'http://128.138.77.51:80',
               'https':'https://128.138.77.51:443',
               'file':'file://128.138.77.51'}    #####
    # proxy = {'http':'http://128.138.77.51:80'}
    # proxy = {'http':'http://128.138.77.51:8080/'}
    proxy = {'http':'http://128.138.77.51'}
    
##    proxies = {'http':'http://128.138.77.51',
##               'https':'https://128.138.77.51',
##               'file':'file://128.138.77.51'}    #####
##    proxy = {'http':'http://128.138.77.51'}
    # proxy = {'http':'http://128.138.77.51/'}
    
    #------------------------------------------
    # Get the hostname (head or compute node)
    #------------------------------------------
    hostname = socket.gethostname()
    print('local hostname =', hostname)   # (e.g. compute node: cl1n020)
    print(' ')

    #----------------------------------------------------
    # Try to open web service URL directly with urllib
    #----------------------------------------------------
    # Using proxies when we're running on the head node
    # seems to be fine.
    #----------------------------------------------------
    print('Trying to open URL with urllib...')
    try:
        # HIS_URL = urllib.urlopen( HIS_location,
        HIS_URL = urllib.request.urlopen( web_service_URL,
                                  proxies=proxies )
        # print HIS_URL
        print('SUCCESS: Opened URL with urllib.urlopen().')
        ## print HIS_URL.read()
    except:
        print('FAILURE: Could not open URL with urllib.urlopen().')
    print(' ')

    ### return  #########################

    #--------------------------------------------------------
    # Try to open web service URL directly with urllib2
    #----------------------------------------------------------------
    # NB! urllib2.urlopen() does not support the keywords
    #     "proxies" or "proxy".  The procedure below, which uses
    #     ProxyHandler, is supposed to be used instead.
    #----------------------------------------------------------------
    # NB! urllib2 does NOT support fetching HTTPS locations
    #     through a proxy. See this link:
    #     http://www.voidspace.org.uk/python/articles/urllib2.shtml
    #----------------------------------------------------------------
    # NB! urllib2 is supposed to look for an environment variable
    #     called "http_proxy" (lower case).  This actually seems
    #     to work, otherwise get the message:
    #
    #     urllib2.URLError: <urlopen error [Errno -3]
    #     Temporary failure in name resolution>
    #----------------------------------------------------------------
    # NB! urllib has been deprecated in favor of urllib2 (at least
    #     in Python 3.0).  See the 2nd website just below.
    #----------------------------------------------------------------
    # NB! For more info on urllib and urllib2, see:
    #         http://docs.python.org/library/urllib2.html
    #         http://docs.python.org/library/urllib.html
    #----------------------------------------------------------------
    print('Trying to open URL with urllib2...')
    # HIS_URL2 = urllib2.urlopen( HIS_location )    # (using http_proxy env var)
    #--------------------------------------------------------
    # HIS_URL2 = urllib2.urlopen( web_service_URL )  # (using http_proxy env var)    
    #--------------------------------------------------------
    # proxy_support = urllib2.ProxyHandler({})   # (override any detected proxies)
    # proxy_support = urllib2.ProxyHandler( proxy )
    # opener = urllib2.build_opener( proxy_support )
    # urllib2.install_opener( opener )
    # HIS_URL2 = urllib2.urlopen( web_service_URL )
    # ERROR: urllib2.HTTPError: HTTP Error 404: Not Found
    # ERROR: urllib2.HTTPError: HTTP Error 404:
    #        /hiscentral/webservices/hiscentral.asmx
    #--------------------------------------------------------
##    proxy_support = urllib2.ProxyHandler( proxy )
##    opener = urllib2.build_opener( proxy_support )
##    opener.open( web_service_URL )
    # ERROR: urllib2.HTTPError: HTTP Error 404: Not Found
    #--------------------------------------------------------
    #print 'SUCCESS: Opened URL with urllib2.urlopen().'
    
    try:
        #--------------------------------------------
        # Use the "http_proxy" environment variable
        #--------------------------------------------
        HIS_URL2 = urllib.request.urlopen( web_service_URL )
        #--------------------------
        # Try to use ProxyHandler
        #--------------------------
##        proxy_support = urllib2.ProxyHandler( proxy )
##        opener = urllib2.build_opener( proxy_support )
##        urllib2.install_opener( opener )
##        HIS_URL2 = urllib2.urlopen( web_service_URL )
        print(HIS_URL2)
        print('SUCCESS: Opened URL with urllib2.urlopen().')
    except:
        print('FAILURE: Could not open URL with urllib2.urlopen().')
    print(' ')
    
    #----------------------------------------   
    # Search parameters (in South Carolina)
    #----------------------------------------
    xmin = '-81.25'
    xmax = '-80.84'
    ymin = '33.84'
    ymax = '34.24'
    keyword    = 'Streamflow'    # (conceptKeyword)
    networkIDs = ''
    start_date = '2000-01-01'    # (2000 is a leap year)
    end_date   = '2000-12-31'

    #------------------------------------     
    # Query HIS Central for time series
    #----------------------------------------------------------------------
    # Note: If calling from the head node (hostname = beach), then
    #       this call fails if we set the proxy=proxies.
    #---------------------------------------------------------------------- 
    #       If we set proxy=proxy (just one entry), then this call
    #       succeeds but get error when accessing the service:
    #       ERROR: No handlers could be found for logger "suds.client"
    #              Exception: (404, u'Not Found')
    #----------------------------------------------------------------------
    # (7/19/10)  This works on both head and compute nodes now, after
    # reconfiguring the system (UnixOps) and adding the environment
    # variable "http_proxy" to ".bash_profile".
    #----------------------------------------------------------------------
    print('Attempting to retrieve data from HIS Central...')
    client = suds.client.Client( web_service_URL )
    print('SUCCESS Opened URL with suds.')
     
    #----------------------------------------------------------------------    
##    proxy2 = {'http':'http://128.138.77.51:80'}
##    tran   = suds.transport.http.HttpTransport( proxy=proxy2, timeout=90 )
##    client = suds.client.Client( web_service_URL, transport=tran )
##    # ERROR:   
    #----------------------------------------------------------------------
##    proxy2 = {'http':'http://128.138.77.51'}
##    tran   = suds.transport.http.HttpTransport( proxy=proxy2, timeout=90 )
##    client = suds.client.Client( web_service_URL, transport=tran )
##    # ERROR: suds.transport.TransportError:
##    #        HTTP Error 404: Not Found
##    #----------------------------------------------------------------------
##    proxy2 = {'http':'http://128.138.77.51:8080'}
##    tran   = suds.transport.http.HttpTransport( proxy=proxy2, timeout=90 )
##    client = suds.client.Client( web_service_URL, transport=tran )
##    # ERROR: suds.transport.TransportError:
##    #        HTTP Error 404: /hiscentral/webservices/hiscentral.asmx
    #----------------------------------------------------------------------
##    print 'Getting transport object...'
##    tran       = suds.transport.http.HttpTransport()
##    print 'Got transport'
##    proxy_han  = urllib2.ProxyHandler( proxy )
##    ## proxy_han  = urllib2.ProxyHandler( proxies )
##    print 'Got proxy handler.'
##    opener     = urllib2.build_opener( proxy_han )
##    print 'Got opener.'
##    tran.urlopener  = opener
##    # tran.urlhandler = proxy_han  # (guessing about urlhandler)
##    client = suds.client.Client( web_service_URL, transport=tran )
    #----------------------------------------------------------------------
    # client = suds.client.Client( HIS_Central_URL, proxy=proxies )
    # ERROR: suds.transport.TransportError: HTTP Error 404: Not Found
    #----------------------------------------------------------------------
    # proxy_URL = HIS_Central_URL + '@128.138.77.51'
    # client = suds.client.Client( proxy_URL, proxy=proxies )
    # ERROR: suds.transport.TransportError: HTTP Error 404: Not Found
    #----------------------------------------------------------------------
    # client = suds.client.Client( HIS_Central_URL )
    # ERROR: Unexpected Python exception:  <urlopen error [Errno -3]
    #        Temporary failure in name resolution>    
##    #----------------------------------------------------------------------
##    # client = suds.client.Client( HIS_Central_URL, transport=tran )
##    # ERROR: suds.transport.TransportError: HTTP Error 404: Not Found
##    #----------------------------------------------------------------------
##    # test_URL = 'http://xxxxxxx/web.asmx?WSDL'
##    # test_URL = 'file:///home/beach/faculty/peckhams/SUDS_TEST.txt'
##    # test_URL = 'file://128.138.77.51/home/beach/faculty/peckhams/SUDS_TEST.txt'
##    # test_URL = 'file://127.0.0.1/home/beach/faculty/peckhams/SUDS_TEST.txt'
##    # test_URL = 'https://water.sdsc.edu/hiscentral/webservices/hiscentral.asmx?WSDL'
##    #   ERROR: urllib2.URLError: <urlopen error [Errno 111] Connection refused>
##    # client = suds.client.Client( test_URL, transport=tran )
##    client = suds.client.Client( HIS_Central_URL, transport=tran,
##                                 location=HIS_Central_location )
##    print 'Got client.'
    # client.set_options( url=HIS_Central_URL, proxy=proxies )
    
    # ERROR: suds.transport.TransportError: HTTP Error 404: Not Found 
    #----------------------------------------------------------------------
    # ERROR: urllib2.URLError: <urlopen error [Errno 2]
    # No such file or directory:
    # 'file://128.138.77.51/home/beach/faculty/peckhams/SUDS_TEST.txt'>
    #----------------------------------------------------------------------
    # ERROR: urllib2.URLError: <urlopen error [Errno 2]
    # No such file or directory:
    # 'file://127.0.0.1/home/beach/faculty/peckhams/SUDS_TEST.txt'>
    #----------------------------------------------------------------------
    # print 'SUCCESS Opened URL with suds.'
    
    #------------------------------------       
    # Get and process response to query
    #------------------------------------
    response = client.service.GetSeriesCatalogForBox2(xmin, xmax, ymin, ymax, \
                               keyword, networkIDs, start_date, end_date)
    series_array = response[0]
    print('Number of time series retrieved =', len(series_array))
    print(' ')

    #-----------------------------------------------------    
    # For each series, call GetValues on the data server
    #-----------------------------------------------------
    for series in series_array:
        client     = suds.client.Client(series.ServURL)
        values_obj = client.service.GetValuesObject(series.location, \
                                                    series.VarCode, \
                                                    start_date, end_date)
        values = values_obj.timeSeries.values
        vlen   = len(values)
        print('len(values)  =', vlen)
        print('type(values) =', type(values))
##        print 'values[0]    =', values[0]
##        if (vlen > 1):
##            print 'values[1]   =', values[1]
##        if (vlen > 2):
##            print 'values[2]   =', values[2]
##        # print 'dir(values)  =', dir(values)

        # vals  = numpy.asarray([float(val.value) for val in values_obj.timeSeries.values])
        # times = numpy.asarray([val._dateTime for val in values_obj.timeSeries.values])
        if (values._count > 5):
            vals  = numpy.asarray([float(val.value) for val in values.value])
            times = numpy.asarray([val._dateTime for val in values.value])

            #---------------------
            # None of these work
            #---------------------
            # vals  = numpy.array( float(values.value) )
            # vals  = numpy.array( numpy.float64(values.value) )
            # vals  = numpy.float64(values.value)
            # vals  = numpy.array( values.value[:], dtype='float64' )
            # times = numpy.array( values._dateTime[:] )
            # vals  = numpy.asarray( values.value, dtype='float64' )
            # times = numpy.asarray( values._dateTime )
            # times = numpy.array( values._dateTime )
            
            print('type(values.value) =', type(values.value))   # (a Python "list")
            print('len(values.value)  =', len(values.value))
            # print 'values.value =', values.value
            # print ' '
            print('size(vals) =', numpy.size(vals))
            print('min(vals), max(vals) =', vals.min(), vals.max())
            
##        values = values_obj.timeSeries.values
##        print 'len(values)      =', len(values)
##        print 'values           =', values
##        times  = values_obj.timeSeries.values._dateTime

        # values = numpy.array(val_obj.timeSeries.values.value)
        

        print('series.Sitename   =', series.Sitename)
        print('series.location   =', series.location)
        print('series.VarName    =', series.VarName)
        print('series.VarCode    =', series.VarCode)
        print('series.ValueCount =', series.ValueCount)
        print('series.datatype   =', series.datatype)
        print('series.valuetype  =', series.valuetype)
        print('series.timeunits  =', series.timeunits)
        print('len(series)       =', len(series))
        print('type(values_obj)  =', type(values_obj))   # ("instance")
        ## print 'type(values)      =', type(values)
##        print 'min(values)      =', min(values)
##        print 'max(values)      =', max(values)    
        print(' ')

#   HIS_test()
#-------------------------------------------------------------------------------
