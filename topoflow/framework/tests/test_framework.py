#! /usr/bin/env python
#
# Copyright (c) 2001-2017, Scott D. Peckham
#
#-----------------------------------------------------------------------

import numpy as np
import os
import tempfile
# See:  http://docs.python.org/2/library/tempfile.html

from topoflow.framework import emeli  ###########
## from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
#
#  topoflow_test()    # Use framework to run TopoFlow.
#  erode_test()
#
#  ref_test()         # For passing references between components.
#
#  framework_test1()  # Test some basic framework functions.
#  framework_test2()
#
#  bobs_erode_test()  # Use framework to run Erode.
#
#-----------------------------------------------------------------------
def topoflow_test( driver_comp_name ='topoflow_driver',
                   cfg_prefix=None, cfg_directory=None,
                   time_interp_method='Linear'):

    #----------------------------------------------------------
    # Note: The "driver_comp_name " defaults to using a
    #       component of "topoflow_driver" type as the driver.
    #       The component of this type specified in the
    #       provider_file will be the driver component.
    #
    #       Any other component in the provider_file can
    #       also be used as the driver.  Examples are:
    #          meteorology
    #          channels
    #          snow
    #          satzone
    #          evap
    #          infil
    #          diversions
    #          ice
    #-----------------------------------------------------

    #-------------------------------------------------------
    # (2/6/13) Since the framework runs the clock now, do
    # we still need to specify a "driver_comp" ??
    # Might still be necessary for use in CSDMS framework.
    #-------------------------------------------------------
    f = emeli.framework()
    examples_dir = emeli.paths['examples']
    ## examples_dir = f.paths['examples']
    
    #--------------------
    # Default arguments
    #--------------------
    if (cfg_prefix == None):
        cfg_prefix = 'June_20_67'
    if (cfg_directory == None):
        cfg_directory = examples_dir + 'Treynor_Iowa_30m/'
    if (driver_comp_name == None):
        driver_comp_name = 'topoflow_driver'   ### (1/30/18)
            
    #------------------------------
    # Run the full TopoFlow model
    #------------------------------
    f.run_model( driver_comp_name =driver_comp_name ,
                 cfg_prefix=cfg_prefix,
                 cfg_directory=cfg_directory,
                 time_interp_method=time_interp_method )

#   topoflow_test()
#-----------------------------------------------------------------------
def topoflow_test2( driver_comp_name ='topoflow_driver',
                   cfg_prefix=None, cfg_directory=None,
                   time_interp_method='Linear'):

    #-----------------------------------------------------------
    # Note: The "driver_comp_name " defaults to using a
    #       component of "topoflow_driver" type as the driver.
    #       The component of this type specified in the
    #       provider_file will be the driver component.
    #
    #       Any other component in the provider_file can
    #       also be used as the driver.  Examples are:
    #          meteorology
    #          channels
    #          snow
    #          satzone
    #          evap
    #          infil
    #          diversions
    #          ice
    #-----------------------------------------------------

    #-------------------------------------------------------
    # (2/6/13) Since the framework runs the clock now, do
    # we still need to specify a "driver_comp" ??
    # Might still be necessary for use in CSDMS framework.
    #-------------------------------------------------------
    f = emeli.framework()
    examples_dir = emeli.paths['examples']
    ## examples_dir = f.paths['examples']
    
    #--------------------
    # Default arguments
    #--------------------
    if (cfg_prefix == None):
        cfg_prefix = 'Test1'
    if (cfg_directory == None):
        cfg_directory = examples_dir + 'C2_basin/'
    
    #------------------------------
    # Run the full TopoFlow model
    #------------------------------
    f.run_model( driver_comp_name =driver_comp_name ,
                 cfg_prefix=cfg_prefix,
                 cfg_directory=cfg_directory,
                 time_interp_method=time_interp_method )

#   topoflow_test2()
#-----------------------------------------------------------------------
def erode_test( cfg_prefix=None, cfg_directory=None,
                time_interp_method='Linear'):
         
    f = emeli.framework()
    driver_comp_name  = 'LEM'
    examples_dir = emeli.paths['examples']
    ## examples_dir = f.paths['examples']

    #--------------------
    # Default arguments
    #--------------------
    if (cfg_prefix == None):
        cfg_prefix = 'Test'
    if (cfg_directory == None):
        cfg_directory = examples_dir + 'Erode_Test/'

    #----------------------
    # Run the Erode model
    #----------------------   
    f.run_model( driver_comp_name =driver_comp_name ,
                 cfg_prefix=cfg_prefix,
                 cfg_directory=cfg_directory,
                 time_interp_method=time_interp_method )

#   erode_test()
#-----------------------------------------------------------------------
def ref_test():

    #---------------------------------------------------------
    # Notes: See "initialize_scalar()" method in BMI.base.py
    #        which is based on this test.
    #---------------------------------------------------------
    import random
    
    class comp1():
        def initialize(self):
            #---------------------------------
            # Case 1. 0-d numpy array.
            # Mutable, if done carefully.
            # Prints without array brackets.
            #---------------------------------
            self.x = np.array(1.0, dtype='float64')

            #--------------------------------------
            # Case 2. 1-d numpy array, 1 element.
            # Mutable, if done carefully.
            # Prints with array brackets.
            #--------------------------------------
            ## self.x = np.array([1.0], dtype='float64')

            #------------------------------------
            # Case 3. numpy scalar;  immutable.
            #------------------------------------           
            ## self.x = np.float64(1)
            
        #-------------------------------------------------
        def update(self):
            #--------------------------
            # These work with Case 1.
            #--------------------------
            # self.x += 1
            # self.x += random.random()
            # self.x *= random.random()
            self.x.fill( random.random() )
            
            #----------------------------
            # Doesn't work with Case 1.
            #----------------------------            
            # self.x = self.x + 1
            # self.x = random.random()
            # self.x[:] = random.random()   #(can't slice 0D)
            # setattr( self, 'x', random.random() )
            #----------------------------------------------------------
            # new_x = np.array( random.random(), dtype='float64' )
            # setattr( self, 'x', new_x )
            #----------------------------------------------------------
            # new_x = np.array( random.random(), dtype='float64' )
            # setattr( self, 'x[:]', new_x )
            
            #----------------------------
            # This works with Case 2.
            #----------------------------            
            # self.x[0] = random.random()
            # self.x[0] += 1  # (works)
            # self.x[:] = random.random()  # (works for 1D, 2D, etc.)
            
            ## self.x = self.x + 1  # (ref is lost; not in place)
        #-------------------------------------------------
        def get_values(self, var_name):
            return getattr( self, var_name )
            #-----------------------------------
            # Using exec like this also works.
            #-----------------------------------
            # exec("result = self." + var_name) in globals(), locals()
            # return result
            
            #------------------------------------------
            # These next three "break" the reference.
            #------------------------------------------
            # return np.float64( result )
            # return result.astype('float64')
            # return np.array( result, dtype='float64' )

    #-----------------------------------------------------
    class comp2():
        def print_x(self):
            print('After c1.update(), c2.x =', self.x)
        #-------------------------------------------------
        def set_values(self, var_name, scalar):
            setattr( self, var_name, scalar )
            #-----------------------------------
            # Using exec like this also works.
            #-----------------------------------
            # exec("self." + var_name + " = scalar") in globals(), locals()
    #-----------------------------------------------------
        
    #---------------------------
    # Instantiate 2 components
    #---------------------------
    c1 = comp1()
    c2 = comp2()
    c1.initialize()
    
    #-------------------------------------
    # Copy reference to x from c1 to c2.
    #-------------------------------------
    vector = c1.get_values('x')
    c2.set_values('x', vector)
    #------------------------------------------    
    # c2.x = c1.x   # (works)
    # exec("c2.x = c1.x") in globals(), locals()  # (works)
    # exec("c2.x = c1.x")       # (problem?)
    # exec("c2." + "x = c1.x")  # (problem?)
    
    for k in range(3):
        c1.update()
        c2.print_x()
    
#   ref_test()
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def framework_test1():

    f = emeli.framework()
    examples_dir = emeli.paths['examples']
    ## examples_dir = f.paths['examples']
    ## driver_comp_name  = 'topoflow_driver'  # (TopoFlow Driver test)

    #-----------------------------------------
    # Set the working directory for test run
    #-----------------------------------------
    f.cfg_directory = examples_dir + 'Treynor_Iowa/'
    f.cfg_prefix = 'June_20_67'
    
    f.read_repository( SILENT=False )  
    print('Components in repository:')
    for comp_name in f.repo_list:
        print('   ' + comp_name)
    print(' ')

    tf_comp_info = f.comp_info[ 'topoflow_driver' ]  
    print('Checking some info for "topoflow_driver":')
    print('    comp_name  =', tf_comp_info.comp_name)
    print('    model_name =', tf_comp_info.model_name)
    print('    uses_ports =', tf_comp_info.uses_ports)
    print(' ')

    #--------------------------------------------
    # Instantiate a complete set of components.
    #--------------------------------------------
    # Now the instantiate() method only allows
    # one component of each "type" (port_name).
    #--------------------------------------------
    comp_list = ['tf_meteorology', 
                 'tf_channels_kin_wave',
                 'tf_infil_green_ampt',
                 'tf_evap_priestley_taylor',
                 'tf_snow_degree_day',
                 'tf_ice_gc2d',
                 'tf_satzone_darcy_layers',
                 'tf_diversions_fraction_method',
                 'topoflow_driver' ]
    for comp_name in comp_list:
        f.instantiate( comp_name, SILENT=False )
    print(' ')
    print('ALL_PYTHON =', f.ALL_PYTHON)
    print(' ')

    #------------------------------------------
    # Check if all uses ports have a provider
    #------------------------------------------
    COMPLETE = f.comp_set_complete( SILENT=False, REPORT=True )
    
    #------------------------------------------------------
    # Try to instantiate another with port_name == snow.
    # It should issue an error message.
    # Also, Treynor_Iowa doesn't have this CFG file yet.
    #------------------------------------------------------
    f.instantiate( 'tf_snow_energy_balance', SILENT=False )

    #------------------------------------------------
    # Remove existing snow component and try again.
    #------------------------------------------------
    f.remove( 'tf_snow_degree_day', SILENT=False )
    f.instantiate( 'tf_snow_energy_balance', SILENT=False )
    
    #-----------------------------------
    # Instantiate all components and
    # store instances in self.comp_set.
    #-----------------------------------
##    skip_list = ['tf_data_his', 'erode_d8_global',
##                 'erode_d8_local', 'd8_global', 'dem_smoother']
##    for comp_name in f.repo_list:
##        if (comp_name not in skip_list):
##            f.instantiate( comp_name )
##            print('Instantiated component named: ' + comp_name)
##    print ' '
    #-------------------------------------------------------------
##    try:
##        for comp_name in f.repo_list:
##            f.instantiate( comp_name )
##            print('Instantiated component named: ' + comp_name)
##    except:
##        print 'ERROR: Could not instantiate component.'
##    print ' '
 
    #-----------------------------
    # Initialize some components
    #-------------------------------------------------------------
    # (3/15/12) The TopoFlow components currently call a method:
    # "initialize_required_components()" that in turn calls:
    # "initialize_ports()", inherited from CSDMS_base.py.
    # But this won't be done with the new approach.
    #-------------------------------------------------------------
    ## cfg_file = None
    ## f.initialize( 'meteorology', cfg_file )
    f.initialize( 'meteorology' )
    print('Initialized component of type: meteorology.')
    ## f.initialize( 'snow', cfg_file)
    f.initialize( 'snow' )
    print('Initialized component of type: snow.')   
    print(' ')

    #--------------------------------------------
    # Copy references from Meteorology to Snow.
    #--------------------------------------------
    print('Copying references from Meteorology to Snow component...')
    provider_name = 'meteorology'
    user_name     = 'snow'
    var_name_list = ['atmosphere_bottom_air__temperature',
                     'land_surface__temperature',
                     'land_surface_net-total-energy__energy_flux',
                     'water-liquid__mass-per-volume_density' ]

    for long_var_name in var_name_list:
        f.connect( provider_name, user_name,
                   long_var_name )
    
    #--------------------------------------------
    # Copy references from Snow to Meteorology.
    #--------------------------------------------
    print('Copying references from Snow to Meteorology component...')
    provider_name = 'snow'
    user_name     = 'meteorology'
    var_name_list = [ 'snowpack__depth',
                      'snowpack__liquid-equivalent_depth',
                      'snowpack__melt_volume_flux',
                      'snowpack__z_mean_of_mass-per-volume_density']
    for long_var_name in var_name_list:
        f.connect( provider_name, user_name,
                   long_var_name )

    #---------------------------------------
    # Check whether the references worked.
    #---------------------------------------
    snow_comp = f.comp_set['snow']
    met_comp  = f.comp_set['meteorology']
    print(' ')
    print('From meteorology, snow got: density of water =', \
          snow_comp.rho_H2O)
          ### snow_comp.met_vars.rho_H2O
    print('From snow, meteorology got: density of snow  =', \
          met_comp.rho_snow)
          ### met_comp.snow_vars.rho_snow
    
    #----------------
    # Final message
    #----------------
    print('Finished with framework_test1.')
    print(' ')
    
#   framework_test1()
#-----------------------------------------------------------------------
def framework_test2():

    f = emeli.framework()
    examples_dir = emeli.paths['examples']
    ## examples_dir = f.paths['examples']    
    ## driver_comp_name  = 'topoflow_driver'  # (TopoFlow Driver test)

    f.read_repository( SILENT=False )

    #-----------------------------------------
    # Set the working directory for test run
    #-----------------------------------------
    f.cfg_directory = examples_dir + 'Treynor_Iowa/'
    f.cfg_prefix = 'June_20_67'  # (needed by initialize())
        
    print('Components in repository:')
    for comp_name in f.repo_list:
        print('   ' + comp_name)
    print(' ')
 
    #-----------------------------------------------------
    # Set self.comp_set_list and self.provider_list
    # from info in the provider file.
    #-----------------------------------------------------
    filename = (f.cfg_prefix + '_providers.txt')
    f.provider_file = (f.cfg_directory + filename)
    f.read_provider_file()
    
    #--------------------------------------------
    # Instantiate a complete set of components.
    #--------------------------------------------
    # Now the instantiate() method only allows
    # one component of each "type" (port_name).
    #--------------------------------------------    
    for comp_name in f.comp_set_list:
        f.instantiate( comp_name, SILENT=False )

    #---------------------------------------------
    # Try to automatically connect every user to
    # a corresponding provider in the comp_set.
    #---------------------------------------------------
    # Provider components are initialized in the order
    # of provider_list and then set references in each
    # component that uses one or more of their vars.
    #---------------------------------------------------
    OK = f.initialize_and_connect_comp_set( REPORT=True )
    if not(OK):
        return
    
    #------------------------------------------------
    # Call update() for each provider, in order.
    # Don't worry about reconciling time steps yet.
    #------------------------------------------------
    print(' ')
    # f.update_all()  # (not ordered)

##    for k in xrange(10):
##        for comp_name in f.provider_list:
##            print 'Calling update() for comp_name =', comp_name
##            if (comp_name != 'topoflow_driver'):
##                bmi = f.comp_set[ comp_name ]
##                bmi.update( -1.0 )
            
    #----------------------------------------------
    # Call the driver component's "update" method
    #----------------------------------------------
##    print ' '
##    print "Calling driver's update() method..."
##    f.update( 'topoflow_driver' )
    
    #----------------
    # Final message
    #----------------
    print('Finished with framework_test2.')
    print(' ')

#   framework_test2()
#-----------------------------------------------------------------------
def bobs_erode_test( driver_comp_name ='LEM',
                     cfg_prefix=None, cfg_directory=None,
                     time_interp_method='Linear'):

    #-------------------------------------------------------
    # (2/6/13) Since the framework runs the clock now, do
    # we still need to specify a "driver_comp" ??
    # Might still be necessary for use in CSDMS framework.
    #-------------------------------------------------------
    f = emeli.framework()
    examples_dir = emeli.paths['examples']
    ## examples_dir = f.paths['examples']   
    
    #--------------------
    # Default arguments
    #--------------------
    if (cfg_prefix == None):
        cfg_prefix = 'Erode_Test_LCP'
    if (cfg_directory == None):
        cfg_directory = os.getenv("HOME") + '/Erode_Tests/Bob_LCP/'
        ## cfg_directory = '~/Erode_Tests/Bob_LCP/'  # (doesn't work)
           
    #------------------------------
    # Run the full TopoFlow model
    #------------------------------
    f.run_model( driver_comp_name =driver_comp_name ,
                 cfg_prefix=cfg_prefix,
                 cfg_directory=cfg_directory,
                 time_interp_method=time_interp_method )

#   bobs_erode_test()
#-----------------------------------------------------------------------


