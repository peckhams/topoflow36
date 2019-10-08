#! /usr/bin/env python
#
#  EMELI = Experimental Modeling Environment for Linking and Interoperability
#
#-------------------------------------------------------
#  Before passing "vars_provided" to time interpolator,
#  we should also restrict list to those that actually
#  change over time, i.e. remove static vars.
#------------------------------------------------------- 
#  Profile the code as follows:
#
#  >>> import cProfile
#  >>> cProfile.run('topoflow.framework.tests.test_framework.topoflow_test()')  
#      
#-----------------------------------------------------------------------      
## Copyright (c) 2012-2016, Scott D. Peckham
##
## Nov 2016.  Updated to use comp_names vs. port_names throughout.
##            EMELI doesn't link components by type, but couples
##            components to others that provide required variables.
##            Changed "port_name" as "comp_type", but not used. 
##            Also removed old XML tags from component repository.
##            Also added unit conversion that uses cfunits.Units.
## 
## Apr 2013.  Added automatic time interpolation, using new
##            time_interpolator class in time_interpolation.py.
##
## Feb 2013.  New approach, starting from framework.py.
## Jan 2013.  Testing with full standard names and removal of
##               all calls to get_port_data(), get_grid_*, etc.
## Mar 2012.  Started from port_queue.py.
## Apr 2012.
## May 2012.  Can remove embed_name tag from XML file now.
##
#-----------------------------------------------------------------------
# Notes: In this version, the "run_model()" method calls a function
#        called "get_required_vars()" within its time loop that gets
#        vars from providers and sets them into all components that
#        need them.  This allows unit conversion and/or regridding,
#        etc. to be applied between the get_values() and set_values()
#        calls.
#
#        This version also supports time interpolation methods of
#        "None" (step-like) and "Linear".  The "run_model()" method
#        calls "initialize_time_interpolation()" to set things up
#        and then calls "update_time_interpolation()" within the
#        time loop.  It doesn't use embedded references.
#
#        See framework0.py for a version that uses only embedded
#        references. (2/18/13)
#
#        See framework2.py for a version without time interpolation.
#
#-----------------------------------------------------------------------
# Notes: The "cfg_directory" is the directory which contains the
#        configuration files for a given model run.  Similarly, the
#        "cfg_prefix" is the filename prefix for the CFG files.
#        The "cfg_directory" should contain a "provider_file" that
#        provides the repository_path ("repo_path") followed by a
#        list of port names and component names, where the component
#        name tells the framework which component to use as the
#        provider of the corresponding portname.
#
#        The "unit_test()" function is used to run the model that
#        is associated with a given provider_file.  It creates an
#        instance of the framework and then calls its "run_model()"
#        method.  See the run_model() method code for details.
#  
#-----------------------------------------------------------------------
#
#  See the "tests" subfolder for tests.
#
#  class comp_data()
#      __init__
#
#  class framework()
#
#      read_repository()
#      read_provider_file()
#      comp_name_valid()
#      comp_set_complete()          # NOT USED. Still uses ports.
#      -------------------------
#      instantiate()
#      remove()
#      connect()
#      disconnect()
#      ----------------------------
#      These are from CMI_base.py
#      ----------------------------
#      get_values()                  # NOT USED. See get_required_vars()
#      set_values()     
#      get_values_at_indices()
#      set_values_at_indices()
#      -------------------------
#      initialize()
#      update()
#      finalize()
#      -------------------------
#      get_cfg_filename()            # (9/18/14)
#      initialize_all()
#      update_all()
#      finalize_all()
#      -------------------------
#      go()
#      run_model()                   # (4/18/13. New way to set refs.)
#      run_rc_script()               # Not ready yet.
#      -------------------------
#      initialize_time_vars()
#      convert_time_units()
#      initialize_framework_dt()
#      update_time()
#
#      ---------------------------------
#      These provide "autoconnection"
#      ---------------------------------
#      find_var_users_and_providers()
#      check_var_users_and_providers()
#      initialize_comp_set()              ## (2/18/13)
#      get_required_vars()                ## (4/18/13)

#      ----------------------------------
#      Alternate approach, not used now
#      ----------------------------------
#      initialize_and_connect_comp_set()  ## (used before 2/18/13)
#      run_model_old()
#      set_provided_vars()                ## (2/18/13, not used now)
#
#-----------------------------------------------------------------------

import numpy as np
import os
# import sys
import time
# import traceback
# import wx
import xml.dom.minidom

from topoflow.framework import time_interpolation    # (time_interpolator class)
# from topoflow.framework import unit_conversion  # (unit_convertor class)
# from topoflow.framework import grid_remapping

# import OrderedDict_backport  # (for Python 2.4 to 2.7)

#--------------------------------------------------------------
# TOGGLED THIS OFF SINCE NOT NEEDED AND INSTALLING
# UDUNITS-2 IS A PAIN, ESPECIALLY ON WINDOWS.  NEED
# TO FIND A PURE PYTHON SOLUTION.
#--------------------------------------------------------------
# Load the unit conversion package.
# This requires installing UDUnits2.2 on your Mac first, e.g.
#     brew install udunits
# and then installing the Python API package: cfunits
#--------------------------------------------------------------
# from cfunits import Units

import sys    #### for testing

#-----------------------------------------------------------------------
# Embed some path info directly into EMELI.
#--------------------------------------------
##################################################################
# NOTE:  "get_package_paths" will not work as intended on Python
# versions less than 3.4 if os.chddir() is called between two
# calls to get_package_paths().  This is a known Python issue
# that has been discussed frequently online.  (9/17/14)
##################################################################

#------------------------------------------------------
# Get path to the current file (emeli.py).  (7/29/13)
# At top need: "#! /usr/bin/env python" ??
# See: https://docs.python.org/2/library/os.path.html
#------------------------------------------------------
framework_dir = os.path.dirname( __file__ )
parent_dir    = os.path.join( framework_dir, '..' )
# parent_dir    = os.path.join( framework_dir, os.path.pardir )
examples_dir  = os.path.join( parent_dir, 'examples' )
#-------------------------------------------------------
framework_dir = os.path.abspath( framework_dir )
parent_dir    = os.path.abspath( parent_dir )
examples_dir  = os.path.abspath( examples_dir )
#-------------------------------------------------------
#     framework_dir = os.path.realpath( framework_dir )
#     parent_dir    = os.path.realpath( parent_dir )
#     examples_dir  = os.path.realpath( examples_dir )
#-------------------------------------------------------
framework_dir = framework_dir + os.sep
parent_dir    = parent_dir    + os.sep
examples_dir  = examples_dir  + os.sep

SILENT = False
if not(SILENT):
    # print ' '
    print('Paths for this package:')
    print('framework_dir = ' + framework_dir)
    print('parent_dir    = ' + parent_dir)
    print('examples_dir  = ' + examples_dir)
    print('__file__      = ' + __file__)
    print('__name__      = ' + __name__)
    print(' ')

#--------------------------------------
# Save the full paths in a dictionary
#--------------------------------------
paths = dict()
paths['framework']        = framework_dir
paths['examples']         = examples_dir
paths['framework_parent'] = parent_dir

## This works fine.
# from topoflow.framework.tests import rti_files   ## (works fine)

## This doesn't work because test_framework.py already imported emeli.py.
# from topoflow.framework.tests import test_framework

  
#-----------------------------------------------------------------------
class comp_data():
    
    def __init__(self, comp_name=None, module_path=None,
                 module_name=None, class_name=None,
                 model_name=None, version=None,
                 language=None, author=None, 
                 help_url=None, cfg_template=None, 
                 time_step_type=None, time_units=None,
                 grid_type=None, description=None,
                 comp_type=None, uses_types=None):
     
        #-------------------------------------        
        # Store static data for a component.
        #-------------------------------------       
        self.comp_name      = comp_name
        self.module_path    = module_path
        self.module_name    = module_name
        self.class_name     = class_name    
        self.model_name     = model_name
        self.version        = version
        self.language       = language
        self.author         = author
        #---------------------------------------       
        self.help_url       = help_url
        self.cfg_template   = cfg_template
        #---------------------------------------
        self.time_step_type = time_step_type
        self.time_units     = time_units
        self.grid_type      = grid_type
        self.description    = description
        self.comp_type      = comp_type
        #-------------------------------
        # This should be a Python list
        #-------------------------------
        self.uses_types     = uses_types

        
#   comp_data (class)
#-----------------------------------------------------------------------
class framework():

    #----------------------------------------
    # Define some unit-conversion constants
    #----------------------------------------
    secs_per_min   = 60
    secs_per_hour  = 60  * secs_per_min
    secs_per_day   = 24  * secs_per_hour 
    secs_per_year  = 365 * secs_per_day
    secs_per_month = secs_per_year / 12    #########

#     ##################################################################
#     # NOTE:  "get_package_paths" will not work as intended on Python
#     # versions less than 3.4 if os.chddir() is called between two
#     # calls to get_package_paths().  This is a known Python issue
#     # that has been discussed frequently online.  (9/17/14)
#     ##################################################################
# 
#     #------------------------------------------------------
#     # Get path to the current file (emeli.py).  (7/29/13)
#     # At top need: "#! /usr/bin/env python" ??
#     # See: https://docs.python.org/2/library/os.path.html
#     #------------------------------------------------------
#     framework_dir = os.path.dirname( __file__ )
#     parent_dir    = os.path.join( framework_dir, '..' )
#     # parent_dir    = os.path.join( framework_dir, os.path.pardir )
#     examples_dir  = os.path.join( parent_dir, 'examples' )
#     #-------------------------------------------------------
#     framework_dir = os.path.abspath( framework_dir )
#     parent_dir    = os.path.abspath( parent_dir )
#     examples_dir  = os.path.abspath( examples_dir )
#     #-------------------------------------------------------
# #     framework_dir = os.path.realpath( framework_dir )
# #     parent_dir    = os.path.realpath( parent_dir )
# #     examples_dir  = os.path.realpath( examples_dir )
#     #-------------------------------------------------------
#     framework_dir = framework_dir + os.sep
#     parent_dir    = parent_dir    + os.sep
#     examples_dir  = examples_dir  + os.sep
# 
#     SILENT = False
#     if not(SILENT):
#         print ' '
#         print 'Paths for this package:'
#         print 'framework_dir =', framework_dir
#         print 'parent_dir    =', parent_dir
#         print 'examples_dir  =', examples_dir
#         print '__file__      =', __file__
#         print '__name__      =', __name__
#         print ' '
#     
#     #--------------------------------------
#     # Save the full paths in a dictionary
#     #--------------------------------------
#     paths = dict()
#     paths['framework']        = framework_dir
#     paths['examples']         = examples_dir
#     paths['framework_parent'] = parent_dir
        
    #-------------------------------------------------------------------
    def read_repository( self, SILENT=True ):

        #---------------------------------------------------
        # Read a file that contains static information for
        # all components in the repository.
        #---------------------------------------------------------
        # Notes:  It is helpful to have a look at the DOM specs,
        #         which can be found online at:
        #         http://www.w3.org/TR/1998/REC-DOM-Level-1-
        #                19981001/introduction.html
        #                (see the tree diagram)
        #         http://www.w3.org/TR/1998/REC-DOM-Level-1-
        #                19981001/level-one-core.html
        #                (search for "firstChild")
        #---------------------------------------------------------

        #---------------------------------------------------------
        # (11/4/13) Now the component repository file is always
        # stored in the same directory as the framework.py file,
        # so we don't need to pass it or store it in the
        # provider_file.  At top of "framework.py" file we need
        # "#! /usr/bin/env python" for this to work.
        #---------------------------------------------------------
        # framework_dir = self.paths['framework']
        framework_dir = paths['framework']   # (paths was saved at top of file)
        repo_dir  = framework_dir
        repo_file = 'component_repository.xml'
        comp_repo_file = repo_dir + repo_file

        #----------------------------------------
        # Make sure comp_repo_file was found
        #----------------------------------------
        ########
    
        if not(SILENT):
            print('Reading info from comp_repo_file:')
            print('    ' + comp_repo_file)
            print(' ')
        
        #-------------------------------------------
        # Read all component info from an XML file
        # into a big string called "doc_string"
        #-------------------------------------------
        repo_unit = open( comp_repo_file, 'r' )
        doc_string = repo_unit.read()
        dom = xml.dom.minidom.parseString( doc_string )

        #----------------------------------------------
        # Count all tags in XML file of various types
        #----------------------------------------------
        C_elements = dom.firstChild.getElementsByTagName("component") 
        n_comps    = len(C_elements)
        if (n_comps == 0):
            print('########################################')
            print(' ERROR: Component repository XML file')
            print('        has no "component" tags.')
            print('########################################')
            print(' ')
            return
            
        #------------------------------------------------------
        # We'll store info for all components in a dictionary
        # where the key is "comp_name" and the return value
        # is a struct-like class.
        #------------------------------------------------------
        self.comp_info = dict()
        self.repo_list = []      # (for all component names)
        
        #------------------------------------------------
        # For each component, get all of its attributes
        # and store them in a "comp_data" object.
        #------------------------------------------------
        for comp in C_elements:
            nodes = comp.getElementsByTagName("comp_name")
            comp_name = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("module_path")
            module_path = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("module_name")
            module_name = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("class_name")
            class_name = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("model_name")
            model_name = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("version")
            version = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("language")
            language = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("author")
            author = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("help_url")
            help_url = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("cfg_template")
            cfg_template = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("time_step_type")
            time_step_type = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("time_units")
            time_units = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("grid_type")
            grid_type = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("description")
            description = nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            nodes = comp.getElementsByTagName("comp_type")
            comp_type = nodes[0].firstChild.data.strip()
            #-------------------------------------------------------------
            nodes = comp.getElementsByTagName("uses_types")
            uses_types_list = nodes[0].firstChild.data.strip().split(",")

            #-----------------------------------------------------
            # Without the "str()", get extra "u'" when printing.
            #-----------------------------------------------------           
#             for k in xrange( len(uses_types_list) ):
#                 uses_types_list[k] = str(uses_types_list[k].strip())
                
            #----------------------------------------
            # Store comp data in "comp_data" object
            # Note: "uses_ports" is a Python list.
            #----------------------------------------
            new_comp_data = comp_data(comp_name=comp_name,
                                      module_path=module_path,
                                      module_name=module_name,
                                      class_name=class_name,
                                      model_name=model_name,
                                      version=version,
                                      language=language,
                                      author=author,
                                      help_url=help_url,
                                      cfg_template=cfg_template,
                                      time_step_type=time_step_type,
                                      time_units=time_units,
                                      grid_type=grid_type,
                                      description=description,
                                      comp_type=comp_type,
                                      uses_types=uses_types_list )
           
            #-----------------------------------------------        
            # Put comp_data in dictionary; key = comp_name
            #-----------------------------------------------
            self.comp_info[ comp_name ] = new_comp_data
            self.repo_list.append( comp_name )

    #   read_repository()
    #-------------------------------------------------------------------
    def read_provider_file(self, SILENT=False):

        #-----------------------------------------------------
        # Notes:  The "provider_file" specifies the specific
        #         components that are to be used to build a
        #         new composite model.  Attributes like
        #         "in_directory" and "out_directory" should
        #         arguably be thought of as config vars for
        #         this new composite model.
        #-----------------------------------------------------

        if not(SILENT):
            print('Reading info from provider_file:')
            print('    ' + self.provider_file)

        self.provider_list = []
        self.comp_set_list = []
        
        file_unit = open( self.provider_file, 'r' )

        #--------------------------------
        # Read provider info into lists
        #--------------------------------
        while (True):
            line = file_unit.readline()

            ## print '## len(line) =', len(line)  #########
            
            if (line == ''):
                break
            if (line[0] != '#'):       # (skip comment/header lines)
                words = line.split()   # (split on white space)
                if (len(words) >= 2):
                    ##########  self.provider_list.append( words[0] ) ##########
                    self.provider_list.append( words[1] )
                    self.comp_set_list.append( words[1] )
        
    #   read_provider_file()
    #-------------------------------------------------------------------
    def comp_name_valid( self, comp_name ):

        #------------------------------
        # Is comp_name in repo_list ?
        #------------------------------
        if (comp_name not in self.repo_list):
            print('#########################################')
            print(' ERROR: There is no component named:')
            print('    ' + comp_name)
            print(' in the component repository.')
            print('#########################################')
            print(' ')
            return False
        else:
            return True
        
    #   comp_name_valid()
    #-------------------------------------------------------------------
#     def comp_set_complete( self, SILENT=True, REPORT=False ):
# 
#         #-----------------------------------------------------
#         # Note:  This checks if all of the "uses" ports for
#         #        components in comp_set have a corresponding
#         #        "provider" port in the comp_set.
#         #-----------------------------------------------------
#         # Note:  This currently ignores the "ppf" port.
#         #-----------------------------------------------------
#         
#         #---------------------------------------------------
#         # Build up a list of all the uses ports and a list
#         # of all the provider ports for the components in
#         # the comp_set (configuration).
#         #---------------------------------------------------
#         uses_ports_found     = []
#         provides_ports_found = []
#         for port_name in self.comp_set:
#             provides_ports_found.append( port_name )
#             #-----------------------------------------
#             info = self.port_info[ port_name ]
#             for name in info.uses_ports:
#                 if (name not in uses_ports_found):
#                     uses_ports_found.append( name )
# 
#         #-------------------------------------
#         # Any uses and provides ports found?
#         #-------------------------------------
#         n_uses     = len( uses_ports_found )
#         n_provides = len(provides_ports_found )
#         if (n_uses == 0):
#             print 'ERROR: No uses ports found.'
#             return
#         if (n_provides == 0):
#             print 'ERROR: No provides ports found.'
#             return
# 
#         #---------------------
#         # Sort the two lists
#         #---------------------
#         uses_ports_found.sort()
#         provides_ports_found.sort()
#         
#         #----------------------------------------------
#         # Option to print all uses and provides ports
#         #----------------------------------------------
#         if (REPORT):
#             print 'Uses ports in comp_set are:'
#             for name in uses_ports_found:
#                 print '    ' + name
#             print 'Provides ports in comp_set are:'
#             for name in provides_ports_found:
#                 print '    ' + name
#             print ' '
#             
#         #------------------------------------------------
#         # Check if every uses_port is in providers
#         #------------------------------------------------
#         missing_providers = []
#         for uses_port_name in uses_ports_found:
#             FOUND = (uses_port_name in provides_ports_found)
#             ######## if not(FOUND):
#             if not(FOUND) and (uses_port_name != 'ppf'):  #########
#                 missing_providers.append( uses_port_name )
#         if (len(missing_providers) > 0):
#             if not(SILENT):
#                 print 'These uses ports have no matching provider:'
#             for name in missing_providers:
#                 print '   ' + name
#             print ' '
#             return False
#         else:
#             if not(SILENT):
#                 print 'All uses ports have matching provider port.\n'
#             return True
#         
#     #   comp_set_complete()  
    #-------------------------------------------------------------------
    def instantiate( self, comp_name, SILENT=True ):

        #-------------------------------------------------------
        # Note: This version only allows one component of each
        #       "type" (given by comp_type) to be instantiated
        #       and included in the comp_set.  However, this
        #       isn't necessary, because EMELI also checks
        #       if a var has multiple providers.
        #-------------------------------------------------------
    
        #------------------------------
        # Is comp_name in repo_list ?
        #------------------------------
        if not(self.comp_name_valid( comp_name )):
            return
        
        #-------------------------------------------------
        # Get info for this component that was read from
        # the component_info_file (repository).
        #-------------------------------------------------
        module_name = self.comp_info[ comp_name ].module_name
        class_name  = self.comp_info[ comp_name ].class_name
        comp_type   = self.comp_info[ comp_name ].comp_type

        #----------------------------------------------
        # Create empty dictionary if not created yet.
        #----------------------------------------------
        if not(hasattr(self, 'comp_set')):
            self.comp_set   = dict()
            # self.port_info  = dict()
            self.ALL_PYTHON = True  # (see below)
            
        #------------------------------------------------
        # Do we already have a component of the type
        # "comp_type" in this comp_set (configuration)?
        #------------------------------------------------
#         if (comp_type in self.comp_set):
#             print '##############################################'
#             print ' ERROR: Cannot instantiate component named:'
#             print '        ' + comp_name
#             print '    because this configuration already has'
#             print '    a component of the type:', comp_type + '.'
#             print '##############################################'
#             print ' '
#             return

        #----------------------------------------------
        # (6/20/13) Add package prefix to module_name
        #----------------------------------------------
        module_prefix    = 'topoflow.components.'          ## COULD GET FROM module_path
        full_module_name = module_prefix + module_name
        
        #--------------------------------------------
        # Import the module (no .py extension) and
        # then create an instance called "comp" and
        # place it in the framework.
        #--------------------------------------------
        ## print '### full_module_name = ', full_module_name
        cmd = 'from topoflow.components import ' + module_name
        exec( cmd, globals(), locals())
        comp = eval( module_name + '.' + class_name + '()' )   ### (2019-10-03)
        ### exec( 'comp = ' + module_name + '.' + class_name + '()', globals())
                        
        #--------------------------------------------------
        # Add new component to the "comp_set" dictionary.
        #--------------------------------------------------
        self.comp_set[ comp_name ] = comp
        #-------------------------------------------
        # Use "comp_type" vs. "comp_name" for key.
        #-------------------------------------------
        ## self.comp_set[ comp_type ] = comp

        #----------------------------------------
        # Copy info from comp_info to port_info
        #----------------------------------------
#         info = self.comp_info[ comp_name ]             #### NO LONGER NEEDED
#         self.port_info[ comp_type ] = info

        #-----------------------------------------
        # Are all components written in Python ?
        #-----------------------------------------
        info = self.comp_info[ comp_name ]
        self.ALL_PYTHON = self.ALL_PYTHON and \
                           (info.language.lower() == 'python')

        #----------------
        # Final message
        #----------------
        if not(SILENT):
            print('Instantiated component: ' + comp_name)
            ## print('        of comp_type: ' + comp_type)

    #   instantiate()
    #-------------------------------------------------------------------
    def remove( self, comp_name, SILENT=True ):

        #-----------------------------------------------------
        # Note: Remove a component from comp_set, perhaps to
        #       replace with another with same comp_type.
        #
        # NB!   This previously assumed that there can only
        #       be one component in comp_set that has a
        #       given comp_type.
        #-----------------------------------------------------

        #----------------------------------------------
        # Return if there is nothing in comp_set yet.
        #----------------------------------------------
        if not(hasattr(self, 'comp_set')):
            return

        #---------------------------------------------   
        # Delete component from comp_set dictionary.
        #---------------------------------------------   
        del self.comp_set[ comp_name ]
    
        #------------------------------------------
        # Get comp_type for this component and
        # delete it from the comp_set dictionary.
        #------------------------------------------
#         comp_type = self.comp_info[ comp_name ].comp_type
#         del self.comp_set[ comp_type ]

        #----------------
        # Final message
        #----------------
        if not(SILENT):
            print('Removed component:', comp_name)
            
    #   remove()
    #-------------------------------------------------------------------
    def connect( self, provider_name, user_name,
                 long_var_name, REPORT=False ):

        if (REPORT):
            print('Connecting user: ' + user_name)
            print('    to provider: ' + provider_name)
            print('    for the variable: ' + long_var_name)

        #---------------------------------------------
        # Get a reference to long_var_name from the
        # component with provider_name (a comp_name)
        #---------------------------------------------
        values = self.get_values( long_var_name, provider_name )
      
        #-------------------------------------------------------
        # This is where service components could be called to
        # do regridding, semantic mediation, unit conversion,
        # and so on.  Right now, we assume that all of these
        # already match.  Note that self.get_values() returns
        # vars on provider's grid, with provider's units, etc.
        #-------------------------------------------------------
        # Call service components here.

        #--------------
        # For testing
        #--------------
        # print 'provider_name =', provider_name
        # print 'user_name =', user_name
        # print 'long_var_name =', long_var_name
        # print '============================================'

        #------------------------------
        # Convert units, if necessary
        #------------------------------
        p_bmi   = self.comp_set[ provider_name ]
        u_bmi   = self.comp_set[ user_name ]
        p_units = p_bmi.get_var_units( long_var_name )
        u_units = u_bmi.get_var_units( long_var_name )
        if (u_units != p_units):
            values = Units.conform( values, Units(p_units), Units(u_units) )

        #---------------------------------------------------
        # Embed a reference to long_var_name from the
        # provider into the (BMI level of) user component.
        #---------------------------------------------------
        self.set_values( long_var_name, values, user_name )

        #-------------------
        # Test on Q_outlet
        #-------------------
##        if (long_var_name == 'basin_outlet_water_x-section__volume_flow_rate'):
##            comp = self.comp_set[ 'hydro_model' ]
##            vals = comp.Q_outlet
##            print '########################################################'
##            print ' After get_values, rank( Q_outlet ) =', np.ndim(values)
##            print ' After set_values, rank( Q_outlet ) =', np.ndim(vals)
##            print '########################################################'

        #---------------------------------------------------
        # This needs to be done using a component's own
        # get_var_name() function, not here.  (5/17/12)
        #---------------------------------------------------
        # embed_name could be null string or the name of a
        # named container (struct) such as "met_vars."
        #---------------------------------------------------
##        embed_name     = self.port_info[ user_name ].embed_name
##        long_var_name2 = (embed_name + long_var_name)
##        self.set_values( long_var_name2, values, user_name )

    #   connect()
    #-------------------------------------------------------------------
    def disconnect( self, provider_name, user_name,
                    long_var_name ):

        #--------------------------------------------
        # Note: We may not need this.  Use "del" ??
        #--------------------------------------------
        pass
        
    #   disconnect()
    #-------------------------------------------------------------------
    def get_values( self, long_var_name, comp_name ):

        #---------------------------------------------------
        # Note: We will later have a PYTHON_ONLY flag that
        #       will allow us to bypass the static-type
        #       BMI methods.
        #---------------------------------------------------
        
        #---------------------------------------------------------
        # Note: This method returns a NumPy "ndarray" object
        #       that Babel is able to pass to other components
        #       as a SIDL generic array.
        #---------------------------------------------------------
        #       The functions called by this one below are "port
        #       functions" (or methods) and therefore have
        #       static return types.
        #---------------------------------------------------------
        #       This component (self) often won't have a port
        #       connection to the caller and so can't call the
        #       methods of its caller.  Regridding and unit
        #       conversion (if supported) would be done just
        #       prior to calling set_values() in connect().
        #---------------------------------------------------------
        bmi = self.comp_set[ comp_name ]

        ###############################################
        # THIS IS A LOT FASTER.  (2/20/13)
        ###############################################
        return bmi.get_values( long_var_name )
    
#         #----------------------------------------------
#         # Get data type and rank for long_var_name.
#         # Assume that NumPy dtype string is returned.
#         #----------------------------------------------
#         dtype = bmi.get_var_type( long_var_name )
#         rank  = bmi.get_var_rank( long_var_name )
# 
#         #-----------------------------------------------------------
#         # Call the appropriate BMI method.  We could do this other
#         # ways in Python, but this is the general approach.
#         #-----------------------------------------------------------
#         if (dtype == 'float64'):
#             if (rank == 0):
#                 values = bmi.get_0d_double( long_var_name )
#             elif (rank == 1):
#                 values = bmi.get_1d_double( long_var_name )
#             elif (rank == 2):
#                 values = bmi.get_2d_double( long_var_name )
#             elif (rank == 3):
#                 values = bmi.get_3d_double( long_var_name )
# ##            elif (rank == 4):
# ##                values = bmi.get_4d_double( long_var_name )
#                 #-------------------------------------------------------------
#         elif (dtype == 'int32'):
#             if (rank == 0):
#                 values = bmi.get_0d_int( long_var_name )
#             elif (rank == 1):
#                 values = bmi.get_1d_int( long_var_name )
#             elif (rank == 2):
#                 values = bmi.get_2d_int( long_var_name )
#             elif (rank == 3):
#                 values = bmi.get_3d_int( long_var_name )
# ##            elif (rank == 4):
# ##                values = bmi.get_4d_int( long_var_name )
#                 #-------------------------------------------------------------
#         else:
#             print '############################################'
#             print ' ERROR: In framework.get_values():'
#             print '        Type_name of "' + dtype + '"'
#             print '        is not yet supported.'
#             print '############################################'
#             print ' '
#             
#         #-----------------------------------------------
#         # Return values on this component's own grid
#         # and with its own units for variables.
#         # Set values will do conversions when needed ?
#         #-----------------------------------------------
#         return values

    #   get_values()   
    #-------------------------------------------------------------------
    def set_values( self, long_var_name, values, comp_name):

        #---------------------------------------------------
        # Note: We will later have a PYTHON_ONLY flag that
        #       will allow us to bypass the static-type
        #       BMI methods.
        #-----------------------------------------------------
        
        #--------------------------------------------
        # Get data type and rank for long_var_name.
        #--------------------------------------------
        bmi = self.comp_set[ comp_name ]

        ###############################################
        # THIS IS A LOT FASTER.  (2/20/13)
        ###############################################
        bmi.set_values( long_var_name, values )
        return
    
#         dtype = str( values.dtype )
#         rank  = np.ndim( values )
# 
#         #------------------------------------------
#         # Use dtype and rank to call appropriate,
#         # static-type BMI functions.
#         #------------------------------------------
#         if (dtype == 'float64'):
#             if (rank == 0):
#                 bmi.set_0d_double( long_var_name, values )
#             elif (rank == 1):
#                 bmi.set_1d_double( long_var_name, values )
#             elif (rank == 2):
#                 bmi.set_2d_double( long_var_name, values )
#             elif (rank == 3):
#                 bmi.set_3d_double( long_var_name, values )
#         #------------------------------------------------------
#         elif (dtype == 'int32'):
#             if (rank == 0):
#                 bmi.set_0d_int( long_var_name, values )
#             elif (rank == 1):
#                 bmi.set_1d_int( long_var_name, values )
#             elif (rank == 2):
#                 bmi.set_2d_int( long_var_name, values )
#             elif (rank == 3):
#                 bmi.set_3d_int( long_var_name, values )
#         #------------------------------------------------------
#         else:
#             print '############################################'
#             print ' ERROR: In framework.set_values():'
#             print '        Type_name of "' + dtype + '"'
#             print '        is not yet supported.'
#             print '############################################'
#             print ' '
                
    #   set_values()
    #-------------------------------------------------------------------
    def get_values_at_indices( self, long_var_name, indices,
                               comp_name ):

        #-------------------------------------------------------------
        # Notes:  For the specified variable, get the values at the
        #         specified indices and return them.  If the wrapped
        #         model is raster, then each raster cell has a long
        #         integer, calendar-style index and these are used
        #         for indices.  If the wrapped model is ugrid, then
        #         each cell in the grid has a unique, long-integer
        #         ID and these are used for indices.
        #-------------------------------------------------------------

        #----------------------------------------------
        # Get data type and rank for long_var_name.
        # Assume that NumPy dtype string is returned.
        #----------------------------------------------
        bmi   = self.comp_set[ comp_name ]
        return bmi.get_values_at_indices( long_var_name, indices )
    
#         dtype = bmi.get_var_type( long_var_name )
#         rank  = bmi.get_var_rank( long_var_name )
# 
#         if (dtype == 'float64'):
#             if (rank == 2):
#                 return bmi.get_2d_double_at_indices( long_var_name, indices )
#         elif (dtype == 'int32'):
#             if (rank == 2):
#                 return bmi.get_2d_int_at_indices( long_var_name, indices )
            
    #   get_values_at_indices()
    #-------------------------------------------------------------------
    def set_values_at_indices( self, long_var_name, indices, values,
                               comp_name ):

        #----------------------------------------------
        # Get data type and rank for long_var_name.
        # Assume that NumPy dtype string is returned.
        #----------------------------------------------
        bmi   = self.comp_set[ comp_name ]
        bmi.set_values_at_indices( long_var_name, indices, values )
        
#         dtype = bmi.get_var_type( long_var_name )
#         rank  = bmi.get_var_rank( long_var_name )
# 
#         if (dtype == 'float64'):
#             if (rank == 2):
#                 bmi.set_2d_double_at_indices( long_var_name, indices, values )
#                 return
#         elif (dtype == 'int32'):
#             if (rank == 2):
#                 bmi.set_2d_int_at_indices( long_var_name, indices, values )
#                 return
# 
#         print '######################################################'
#         print ' ERROR: framework.set_values_at_indices() currently'
#         print '        requires target variable to be a 2D array.'
#         print '######################################################'
#         print ' '
                
    #   set_values_at_indices()
    #-------------------------------------------------------------------
    def initialize( self, comp_name, cfg_file=None,
                    mode='nondriver'):

        #------------------------------
        # Is comp_name in comp_list ?
        #------------------------------
##        if not(self.comp_name_valid( comp_name )):
##            return
           
        bmi = self.comp_set[ comp_name ]
        if (cfg_file == None):
            cfg_file = self.get_cfg_filename( bmi )
        bmi.initialize( cfg_file=cfg_file, mode=mode )
            
    #   initialize()
    #-------------------------------------------------------------------
    def update( self, comp_name ):

        #------------------------------
        # Is comp_name in comp_list ?
        #------------------------------
##        if not(self.comp_name_valid( comp_name )):
##            return

        #----------------------------------------------------
        # Note: Passing dt=-1.0 to update() method tells
        #       component to use its own timestep.
        #       And if we don't pass "time_sec" to update()
        #       and then on to write_output_files() it will
        #       use its internal self.time_sec for output.      
        #----------------------------------------------------      
        bmi = self.comp_set[ comp_name ]
        bmi.update( -1.0 )
        ## bmi.update( -1.0, time_seconds=self.time_sec )

        #---------------------------------------------
        # Was this update was successful ?
        # Status could be "updated", "failed" or
        # "initialized" (if comp_status == Disabled)
        #---------------------------------------------
        status = bmi.get_status()
        if (status == 'failed'):
            print('================================================')
            print('ERROR: Model run aborted.')
            print('  Update failed on component: ' + comp_name + '.')
            print('================================================')
            print(' ')
            self.DONE = True
    
    #   update()
    #-------------------------------------------------------------------
    def finalize( self, comp_name ):

        #------------------------------
        # Is comp_name in comp_list ?
        #------------------------------
##        if not(self.comp_name_valid( comp_name )):
##            return
        
        bmi = self.comp_set[ comp_name ]
        bmi.finalize()
            
    #   finalize()
    #-------------------------------------------------------------------
##    def instantiate_all( self ):
##
##        #--------------------------------------------------
##        # Note: The "comp_set_list" has the comps in
##        #       a particular order, while self.comp_set
##        #       is an unordered dictionary.
##        #--------------------------------------------------
##        if not(hasattr(self, 'comp_set_list')):
##            print 'Providers not yet read from provider_file.'
##            return
##        
##        for comp_name in self.comp_set_list:
##            self.instantiate( comp_name, SILENT=False )
##            
##    #   instantiate_all()
    #-------------------------------------------------------------------
    def get_cfg_filename( self, bmi ):  ###################
    
        #---------------------------
        # Get name of the cfg_file
        #---------------------------
        cfg_extension = bmi.get_attribute( 'cfg_extension' )   ###### NO LONGER IN BMI
        file_name     = (self.cfg_prefix + cfg_extension)
        #-------------------------------------------------
        ## cfg_directory = (os.getcwd() + os.sep)
        cfg_directory = self.cfg_directory + os.sep
        cfg_file      = (cfg_directory + file_name) 

        #--------------
        # For testing
        #--------------
#         print 'In EMELI, cfg_file =', cfg_file

        return cfg_file
          
    #   get_cfg_filename()
    #-------------------------------------------------------------------
    def initialize_all( self, cfg_prefix=None, mode='nondriver'):

        ####### THIS IS NOT USED YET #########

        #--------------------------------------------------
        # Note: The "provider_list" has the ports in
        #       a particular order, while self.comp_set
        #       is an unordered dictionary.
        #--------------------------------------------------
        if not(hasattr(self, 'provider_list')):
            print('Providers not yet read from provider_file.')
            return
        
        for comp_name in self.provider_list:
            bmi = self.comp_set[ comp_name ]
            cfg_file = self.get_cfg_filename( bmi )
            bmi.initialize( cfg_file=cfg_file, mode=mode )
            
    #   initialize_all()
    #-------------------------------------------------------------------
    def update_all( self ):

        #--------------------------------------------------
        # Note: The "provider_list" has the ports in
        #       a particular order, while self.comp_set
        #       is an unordered dictionary.
        #--------------------------------------------------
        if not(hasattr(self, 'provider_list')):
            print('Providers not yet read from provider_file.')
            return
        
        for comp_name in self.provider_list:
            bmi = self.comp_set[ comp_name ]
            bmi.update( -1.0 )
            
    #   update_all()
    #-------------------------------------------------------------------
    def finalize_all0( self ):

        #--------------------------------------------------
        # Note: The "provider_list" has the ports in
        #       a particular order, while self.comp_set
        #       is an unordered dictionary.
        #--------------------------------------------------
        if not(hasattr(self, 'provider_list')):
            print('Providers not yet read from provider_file.')
            return
        
        for comp_name in self.provider_list:
            bmi = self.comp_set[ comp_name ]
            bmi.finalize()
          
    #   finalize_all0()
    #-------------------------------------------------------------------
    def finalize_all( self ):

        #------------------------------------------------------
        # Note: This version calls "get_required_vars()"
        #       one more time for each component.  This
        #       allows a component later in provider_list,
        #       like "topoflow_driver" to get and use
        #       final values from other components. (8/20/13)
        #------------------------------------------------------
        # Note: The "provider_list" has the ports in
        #       a particular order, while self.comp_set
        #       is an unordered dictionary.
        #--------------------------------------------------
        if not(hasattr(self, 'provider_list')):
            print('Providers not yet read from provider_file.')
            return
        
        for comp_name in self.provider_list:
            bmi = self.comp_set[ comp_name ]
            
            #-----------------------------------------------------
            # Get current time of component with this comp_name.
            # Convert units to framework time units, if needed.
            #-----------------------------------------------------
            bmi_time_units = bmi.get_time_units()
            bmi_time       = bmi.get_current_time()
            bmi_time = self.convert_time_units( bmi_time, bmi_time_units )

            #---------------------------------------------
            # Use get_values()/set_values() calls to get
            # latest vars that this component needs from
            # other components.
            #---------------------------------------------
            self.get_required_vars( comp_name, bmi_time )

            #-----------------------------------------------
            # This finalize() call will now have access to
            # final variables from other components, as
            # long as their finalize() was called first.
            #-----------------------------------------------
            bmi.finalize()
        
    #   finalize_all()
    #-------------------------------------------------------------------
    def run_model( self, driver_comp_name='topoflow_driver',
                   ## driver_comp_name='hydro_model',
                   cfg_directory=None, cfg_prefix=None,
                   time_interp_method='Linear'):
        ## (rename to run_comp_set ????)
        
        #-------------------
        # Default settings
        #-------------------
        DEBUG = True
        ## DEBUG = False
        if (cfg_prefix == None):
            print('ERROR: The "cfg_prefix" argument is required.')
            return
        if (cfg_directory == None):
            print('ERROR: The "cfg_directory" argument is required.')
            return
        
        #--------------------------------------------------
        # (11/4/13) Expand things like ".." and "~", then
        # check if need to add path separator at the end.
        #--------------------------------------------------
        cfg_directory = os.path.realpath( cfg_directory )
        # cfg_directory = cfg_directory + os.sep     #########
        self.cfg_prefix    = cfg_prefix
        self.cfg_directory = cfg_directory

        #-----------------------------------------------------
        # Set self.comp_set_list and self.provider_list
        # from info in the provider file, including the
        # repository path "repo_path".
        #-----------------------------------------------------
        provider_file = (cfg_directory + os.sep + cfg_prefix + '_providers.txt')
        # self.provider_file = (cfg_prefix + '_providers.txt')
        self.provider_file = provider_file
        self.read_provider_file()

        #------------------------------------
        # Get the component repository info
        #---------------------------------------------------
        # repo_path is no longer stored in provider_file
        # and doesn't have to be passed in. It is now the
        # directory that contains "framework.py" (11/4/13)
        #---------------------------------------------------
        self.read_repository( SILENT=False )
        
        #------------------------------------
        # Get the component repository info
        #------------------------------------
        ## repo_file = self.repo_path + 'component_repositor'
        ## self.read_repository( repo_file, SILENT=False )
##        print 'Components in repository:'
##        for comp_name in f.repo_list:
##            print '   ' + comp_name
##        print ' '
            
        #--------------------------------------------
        # Instantiate a complete set of components.
        #--------------------------------------------
        # Now the instantiate() method only allows
        # one component of each "type" (comp_type).
        #--------------------------------------------    
        for comp_name in self.comp_set_list:
            self.instantiate( comp_name, SILENT=False )
        ### self.instantiate_all()   ### Later; change it first.
       
        #---------------------------------------------
        # Try to automatically connect every user to
        # a corresponding provider in the comp_set.
        #---------------------------------------------------
        # Provider components are initialized in the order
        # of provider_list and then set references in each
        # component that uses one or more of their vars.
        #---------------------------------------------------
##        OK = self.initialize_and_connect_comp_set( REPORT=True )
##        if not(OK):
##            return

        #------------------------------------------------------------
        # (2/18/13) Previous version of framework class initialized
        # and then connected the components in the comp_set using
        # embedded references.  In this version we will call a
        # "get_required_vars()" method within the time loop, which
        # will in turn call the get_values() and set_values()
        # methods.  But we still need to initialize the comp set.
        #------------------------------------------------------------
        ## OK = self.initialize_comp_set( REPORT=True )
        OK = self.initialize_comp_set( REPORT=False )
        if not(OK):
            return
        
        #---------------------------------------
        # Set mode of the driver component.
        # Note: Must happen before next block.
        #---------------------------------------
        driver = self.comp_set[ driver_comp_name ]
        driver.mode = 'driver'
        print('Driver component name = ' + driver_comp_name)
        print(' ')
        
        #-----------------------------------
        # Initialize all time-related vars
        #-----------------------------------
        self.initialize_time_vars()
        self.initialize_framework_dt()

        #------------------------------------
        # Instantiate a "time_interpolator"
        #------------------------------------
        time_interpolator = time_interpolation.time_interpolator(
                                          self.comp_set,
                                          self.provider_list,
                                          self.vars_provided,
                                          time_interp_method )
        #------------------------------------------------------
        # This calls bmi.update(-1) on every comp in comp_set
        #------------------------------------------------------
        time_interpolator.initialize()
        #---------------------------------------------------------
        # This is used by set_provided_vars() in run_model_old()
        # and by get_required_vars() in run_model().
        #---------------------------------------------------------
        self.time_interpolator = time_interpolator
        
        while not(self.DONE):

            # try:
            #-------------------------------------------------
            # Update components that are ready to be updated
            #----------------------------------------------------
            # (2/18/13) The "provider_list" we loop over
            # here refers to entries in the provider_file; they
            # don't necessarily provide anything to another
            # component in the comp_set.
            #----------------------------------------------------
            # It might be more clear to change these names:
            #     comp_name     -> provider_name
            #----------------------------------------------------
            ## for bmi in self.comp_set:
            for comp_name in self.provider_list:
                bmi = self.comp_set[ comp_name ]

                #-----------------------------------------------------
                # Get current time of component with this comp_name.
                # Convert units to framework time units, if needed.
                #-----------------------------------------------------
                bmi_time_units = bmi.get_time_units()
                bmi_time       = bmi.get_current_time()
                bmi_time = self.convert_time_units( bmi_time, bmi_time_units )

                #------------------------------------
                # Is it time to call bmi.update() ?
                #------------------------------------
                if (self.time > bmi_time):
                    #---------------------------------------------
                    # Use get_values()/set_values() calls to get
                    # latest vars that this component needs from
                    # other components.
                    #---------------------------------------------
                    self.get_required_vars( comp_name, bmi_time )
                    
                    bmi.update( -1.0 )
                    
                    #--------------------------------------------------
                    # Update time interpolation vars for every
                    # long_var_name that is provided by this provider.
                    # Interpolation methods = 'None', 'Linear', etc.
                    #--------------------------------------------------
                    time_interpolator.update2( comp_name )
        
                #------------------------------------------------
                # (2/18/13) Use get_values()/set_values() calls
                # here to set latest vars from this component
                # into all user components that need it.
                #------------------------------------------------
                # This also calls service components as needed.
                #------------------------------------------------
                # self.set_provided_vars( comp_name )
     
            #--------------------
            # Are we done yet ?
            #--------------------
            self.DONE = (driver.DONE or self.DONE)    ####
            self.update_time()
            ## print 'time =', self.time
                
##            except:
##                print 'ERROR in run_model() method at:'
##                print '   time_index =', self.time_index
##                self.status = 'failed'
##                self.DONE = True
                
        #-------------------------
        # Finalize the model run
        #-------------------------
        self.finalize_all()
            
    #   run_model()
    #-------------------------------------------------------------------
    def run_rc_script( self ):

        #----------------------------------------------------------
        # Note: This method will eventually be used to run all
        #       of the framework commands in a Ccaffeine RC file
        #       such as set_path(), instantiate(), connect(),
        #       remove(), etc.
        #----------------------------------------------------------
        pass
    
    #   run_rc_script()
    #-------------------------------------------------------------------
    def initialize_time_vars(self, units='seconds'):

        #-------------------------------------------------
        # Note: Copied here from BMI_base.py on 5/18/12.
        #-------------------------------------------------
        
        #------------------
        # Start the clock
        #------------------
        self.start_time = time.time()
        
        #--------------------------------
        # Initialize the time variables
        #--------------------------------
        self.time_units = units.lower()
        self.time_index = np.int32(0)
        self.time       = np.float64(0)
        self.DONE       = False

        #-------------------------------------------
        # (2/5/12) Set default stopping method:
        # (0=n_steps, 1=stop_time, 2=steady-state)
        #-------------------------------------------
        # Don't do this in set_constants().
        #-------------------------------------------
        if not(hasattr(self, 'stop_code')):
            self.stop_code = 0

        #---------------------------------------
        # Make sure n_steps is defined for new
        # use in update_time().  Was not set
        # for diversions_base. (11/14/11)
        #---------------------------------------
        if not(hasattr(self, 'n_steps')):
            self.n_steps = 1
 
        #--------------------------
        # Time conversion factors
        #--------------------------
        self.sec_per_year = np.float64(365) * 24 * 3600
        self.min_per_year = np.float64(365) * 24 * 60
        
        #-------------------------------------------
        # For backward compatibility with TopoFlow
        #-------------------------------------------
        self.time_sec = np.float64(0)
        self.time_min = np.float64(0)
            
        #--------------------------------------------
        # For print_time_and_value() function below
        #--------------------------------------------
        # Substract 100 seconds so we'll always
        # print values at time zero. (6/29/10)
        #--------------------------------------------
        self.last_print_time = time.time() - 100.0
        
##        self.last_check_time  = time.time()  # (for user interrupt)
##        self.last_plot_time   = float32(0.0)   ### CHECK ###
        
    #   initialize_time_vars()
    #-------------------------------------------------------------------
    def convert_time_units( self, in_time, in_units ):

        #-----------------------------------------------
        # Note:  Conversion constants are defined just
        #        inside (at top of) class declaration.
        #-----------------------------------------------

        #----------------------------------
        # Convert "in_units" to "seconds"
        #----------------------------------
        if (in_units in ['years', 'y']):
            time = in_time * self.secs_per_year
        elif (in_units == 'months'):            ### Use 'm' ????
            time = in_time * self.secs_per_month
        elif (in_units in ['days', 'd']):
            time = in_time * self.secs_per_day
        elif (in_units in ['hours', 'h']):
            time = in_time * self.secs_per_hour
        elif (in_units in ['minutes','m']):     ### month?
            time = in_time * self.secs_per_min
        else:
            time = in_time.copy()
            ## time = in_time

        return time
    
    #   convert_time_units()
    #-------------------------------------------------------------------
    def initialize_framework_dt( self ):
        
        #-----------------------------------------------
        # Get the time steps of each comp in comp_set.
        #-----------------------------------------------
        n_comps = len( self.provider_list )
        dt_array = np.zeros( n_comps )
        #-----------------------------------------------------------
        # Next line worked in Python 2.7, but does not work in 3.7
        #-----------------------------------------------------------
        ### dt_units = np.zeros( n_comps, dtype='|S30')   ###
        dt_units = np.zeros( n_comps, dtype='<U30')   # (2019-10-03)
        k = 0
        print('Original component time step sizes =')
        for comp_name in self.provider_list:    
            bmi      = self.comp_set[ comp_name ]
            dt       = bmi.get_time_step()
            units    = bmi.get_time_units()
            unit_str = ' [' + units + ']'
            #---------------------------------------------------            
            DISABLED = (bmi.comp_status.lower() == 'disabled')
            if not(DISABLED):
                print('    ' + comp_name + ' = ' + str(dt) + unit_str) 
                dt_array[ k ] = dt
                dt_units[ k ] = units
            else:
                print('    ' + comp_name + ' is Disabled.' )
                #--------------------------------------------
                # Set to large value so won't affect min dt
                #--------------------------------------------
                dt_array[ k ] = 1000000.0
                dt_units[ k ] = 'seconds'
            k += 1 

        #-----------------------------------------
        # Convert all time step units to seconds
        #-----------------------------------------
        print('Converting all time step units to seconds...')
        for k in range( n_comps ):
            dt    = dt_array[ k ]
            units = dt_units[ k ]
            dt2 = self.convert_time_units( dt, units )
#             print('  dt       = ' + str(dt) )
#             print('  in_units = ' + units )
#             print('  dt2      = ' + str(dt2) )
#             print('--------------------------------')
            dt_array[ k ] = dt2

        #--------------------------------
        # Set the "framework timestep".
        #------------------------------------------------------
        # Which component has the smallest timestep ?
        # There may be more than one, e.g. topoflow, channels
        # dt_min is used as the "framework timestep".
        #------------------------------------------------------
        dt_min = dt_array.min()
        self.dt = dt_min        # (framework timestep)
        dt_min_name = self.provider_list[ dt_array.argmin() ]
        print('Component with smallest time step is: ' + dt_min_name)
        print(' ')

        #----------------------------------------------
        # This is not used with new method. (4/13/13)
        #----------------------------------------------
        # Compute an "update_step" for each component
        #----------------------------------------------
##        self.comp_update_steps = np.zeros( n_comps )
##        for k in xrange( n_comps ):
##            dt = dt_array[ k ]
##            self.comp_update_steps[ k ] = np.ceil( dt / dt_min ).astype('Int32')  
##        self.n_comps = n_comps  # (used by run_model().)
        
    #   initialize_framework_dt()
    #-------------------------------------------------------------------
    def update_time(self, dt=-1):

        #-------------------------------------------------
        # Note: Copied here from BMI_base.py on 5/18/12.
        #-------------------------------------------------
        
        #---------------------
        # Increment the time
        #---------------------
        self.time_index += 1
        if (dt == -1):
            self.time += self.dt  # (use same units as dt)
        else:
            self.time += dt

##        print '---------------------------------'
##        print 'FRAMEWORK DT =', self.dt
##        print 'FRAMEWORK TIME =', self.time
        
        #--------------------
        # Are we DONE yet ?
        #-------------------------------------------------
        # Note that the components inherit their own
        # "update_time()" method from BMI_base.py which
        # calls "check_finished()".  Don't do this here.
        #-------------------------------------------------

        #------------------------------------------
        # Compute and store time_sec and time_min
        #------------------------------------------
        if (self.time_units == 'seconds'):
            self.time_sec = self.time                          # [seconds]
            self.time_min = self.time_sec / np.float64(60)  # [minutes]
        elif (self.time_units == 'years'):
            #-----------------------------------
            # Used by GC2D and Erode (12/4/09)
            #-----------------------------------
            self.time_sec = self.time * self.sec_per_year  ####
            self.time_min = self.time_sec / np.float64(60)  # [minutes]
            
    #   update_time()
    #-------------------------------------------------------------------
    def find_var_users_and_providers( self, REPORT=False ):

        self.var_providers   = dict()
        self.var_users       = dict()
        #--------------------------------
        self.all_input_var_names  = []
        self.all_output_var_names = []
        
        for comp_name in self.comp_set:
            bmi = self.comp_set[ comp_name ]
            #--------------------------------------------------------
            # Create a dictionary that maps every long_var_name
            # that can be provided by this comp_set to a comp_name.
            #--------------------------------------------------------
            output_var_names = bmi.get_output_var_names()
            for long_var_name in output_var_names:
                if (long_var_name != ''):
                    if (long_var_name not in self.all_output_var_names):
                        self.all_output_var_names.append( long_var_name )
                        self.var_providers[ long_var_name ] = []
                    self.var_providers[ long_var_name ].append( comp_name )
                 
            #--------------------------------------------------------
            # Create a dictionary that maps every long_var_name
            # that is required by this comp_set to a comp_name.
            #--------------------------------------------------------                 
            input_var_names = bmi.get_input_var_names()
            for long_var_name in input_var_names:
                ## print '### long_var_name =' + long_var_name + '___'
                if (long_var_name != ''):
                    if (long_var_name not in self.all_input_var_names):
                        self.all_input_var_names.append( long_var_name )
                        self.var_users[ long_var_name ] = []
                    self.var_users[ long_var_name ].append( comp_name )

        #---------------------------------------------------
        # Sort the lists of all input and output var names
        # The "sort()" method doesn't remove duplicates.
        #---------------------------------------------------
        self.all_input_var_names.sort()
        self.all_output_var_names.sort()
        if (REPORT):
            print('Input variables required by this comp_set:')
            for name in self.all_input_var_names:
                print('    ' + name)
            print(' ')
            print('Output variables provided by this comp_set:')
            for name in self.all_output_var_names:
                print('    ' + name)
            print(' ')

        #--------------------------------------------------- 
        # Make list of output_vars actually used by others
        # (2/18/13) This may improve runtime performance.
        #-------------------------------------------------------
        # Before passing "vars_provided" to time interpolator,
        # we should also restrict list to those that actually
        # change over time, i.e. remove static vars.
        #-------------------------------------------------------   
        self.vars_provided = dict()
        for comp_name in self.comp_set:
            bmi = self.comp_set[ comp_name ]
            output_var_names = bmi.get_output_var_names()
            self.vars_provided[ comp_name ] = []   #####
            for long_var_name in output_var_names:
                if (long_var_name in self.all_input_var_names):
                    self.vars_provided[ comp_name ].append( long_var_name )
            
    #   find_var_users_and_providers()
    #-------------------------------------------------------------------    
    def check_var_users_and_providers( self ):

        #------------------------------------------------------
        # Note: Check that there is at least one provider for
        #       each long_var_name in all_long_var_names.
        #       If there isn't, or if there is more than one
        #       provider, print a warning message.
        #------------------------------------------------------
        OK = True
        for long_var_name in self.all_input_var_names:
            if (long_var_name in self.var_providers):
                provider_list = self.var_providers[ long_var_name ]
                n_providers   = len( provider_list )
                if (n_providers > 1):
                    print('========================================')
                    print(' WARNING: Found multiple providers in')
                    print('          comp_set for long_var_name:')
                    print('            ' + long_var_name)
                    print(' Providers are:', provider_list)
                    OK = False
            else:
                print('========================================')
                print(' WARNING: Found no providers in this')
                print('          comp_set for long_var_name:')
                print('            ' + long_var_name)
                ### print,'          requested by: ' + ??????
                OK = False

        return OK
    
    #   check_var_users_and_providers()   
    #-------------------------------------------------------------------
    def initialize_comp_set( self, REPORT=False ):

        #-------------------------------------------------------------
        # Note: This first calls find_var_users_and_providers()
        #       to find and save self.all_input_var_names and
        #       self.all_output_var_names.  It then calls
        #       check_var_users_and_providers() to determine if
        #       self.comp_set has a provider (in self.var_providers)
        #       for every long_name in self.all_input_var_names.
        #-------------------------------------------------------------               
        # (5/17/12) The TopoFlow components used to call a method:
        # "initialize_required_components()" that in turn called:
        # "initialize_ports()", inherited from CSDMS_base.py.
        # But this is not done with the new approach.
        #-------------------------------------------------------------
        
        #----------------------------------------------
        # Find and check variable users and providers
        #----------------------------------------------
        self.find_var_users_and_providers()
        OK = self.check_var_users_and_providers()
        if not(OK):
            return OK    # (a report was just printed)

        #------------------------
        # For testing (5/17/12)
        #------------------------
##        print '##########################################'
##        print "var_users[ 'channel_time_step_size' ] ="
##        for each in self.var_users[ 'channel_time_step_size' ]:
##            print '    ' + each
##        print '##########################################'              
        
        #################################################################
        # We used to call "os.chdir()" here to change to the driver's
        # CFG directory.  This made the CFG directory the default path
        # for finding files, consistent with setting "in_directory" in
        # a CFG file to ".".  However, this created several unwanted
        # side effects.  For example, the method used to get the full
        # pathnames to the topoflow package and its subdirectories
        # (like "examples"), which uses "__file__" does not work as
        # intended for Python versions less than 3.4 when os.chdir()
        # is called.  For example, if we run topoflow_test() and then
        # erode_test(), or try to run topoflow_test() twice, it will
        # fail with a path problem.  To avoid these problems, we now
        # (9/21/14):
        # (1) Set CFG directory from CFG file in check_directories()
        #     in BMI_base.py.
        # (2) Change an "in_directory" read from a CFG file as "." to
        #     the CFG directory from (1).
        # (3) Avoid all use of os.chdir() in EMELI.
        # (4) Use full pathnames to files everywhere.
        ################################################################
#         if (self.cfg_directory is not None):
#             os.chdir( self.cfg_directory )

        #------------------------------------------------
        # Loop over providers in order of provider_list.
        # Note that the dictionary, self.comp_set, is
        # not ordered.  Order of initialize() matters.
        #------------------------------------------------         
        for provider_name in self.provider_list:
            bmi = self.comp_set[ provider_name ]
            #-------------------------------------------
            # Initialize the provider component to set
            # all of its variables, etc.
            #-------------------------------------------
            cfg_file = self.get_cfg_filename( bmi )
            self.initialize( provider_name, cfg_file )
            print('Initialized component: ' + provider_name + '.')
            ## print 'Initialized component of type: ' + provider_name + '.'
            ## print 'Initialized: ' + comp_name + '.'

            ####################################################
            # (2/18/13) This connection step which creates the
            # embedded references is excluded in this version.
            #
            # Put it back for initialization only.
            ####################################################
            
            #--------------------------------------------------------
            # Create a dictionary that maps every long_var_name
            # that can be provided by this comp_set to a comp_name.
            #--------------------------------------------------------
            output_var_names = bmi.get_output_var_names()
            for long_var_name in output_var_names:
                #-----------------------------------------------------
                # Is this long_var_name used by any other component?
                #-----------------------------------------------------
                try:
                    user_list = self.var_users[ long_var_name ]
                except:
                    user_list = []
                for user_name in user_list:
                    #---------------------------------------------
                    # Embed a reference in the user component to
                    # a variable in the provider component.
                    #---------------------------------------------
                    self.connect( provider_name, user_name,
                                  long_var_name, REPORT=REPORT )

        return OK
    
    #   initialize_comp_set()
    #-------------------------------------------------------------------
    def get_required_vars( self, user_name, bmi_time ):    

        #----------------------------------------------------------
        # Note:  This routine loops through all of the components
        #        that the component given by "user_name"
        #        neeeds and gets/sets the required variables.
        #        It is called just *before* a component update().
        #----------------------------------------------------------
        u_bmi = self.comp_set[ user_name ]  # (or pass bmi)
        
        input_var_names = u_bmi.get_input_var_names()
        
        for long_var_name in input_var_names:
            #-----------------------------------------------------
            # Note: "check_var_users_and_providers()" made sure
            # there is only one provider for each long_var_name.
            #-----------------------------------------------------
            provider_list = self.var_providers[ long_var_name ]
            provider_name = provider_list[0]
            p_bmi  = self.comp_set[ provider_name ]
            
            #---------------------------------------------
            # Get a reference to long_var_name from the
            # component with provider_name (a comp_name)
            #---------------------------------------------
            # values = self.get_values( long_var_name, provider_name )

            #------------------------------------
            # Break the reference for testing ?
            #------------------------------------
            # values = np.float64( values )
            
            #------------------------------------------------
            # Call Time Interpolator to get values that are
            # time interpolated to user's current time.
            #--------------------------------------------------
            # Do we need to process providers in run_model2()
            # in order of decreasing time_step size to ensure
            # that we don't request values at a time that is
            # beyond the provider's time ?  See the tests in
            # time_interpolator.get_values(). (4/18/13)
            # Do providers need to be "out in front" of users?
            #--------------------------------------------------
            values = self.time_interpolator.get_values( long_var_name,
                                                        provider_name,
                                                        bmi_time )
            #-----------------------------------------
            # Call Unit Converter to convert from
            # provider's units to this user's units.
            #---------------------------------------------------
            # The unit_converter should have an initialize()
            # method that gets and stores (in a dictionary),
            # all of the factors and offsets that will be
            # needed during the (composite) model run.
            # Its "convert" method will then use these, but
            # only when values change.
            #---------------------------------------------------
            # Better to pass u_bmi and p_bmi and make the
            #   calls to their get_var_units() in convert()?
            #---------------------------------------------------
            # new_values = a(units, new_units)*values +
            #              b(units, new_units)
            #---------------------------------------------------           
            # u_units = u_bmi.get_var_units( long_var_name )
            # p_units = p_bmi.get_var_units( long_var_name )
            # values  = self.unit_converter.convert( values, p_units, u_units )

            #---------------------------------------------------
            # Convert units, if necessary, with cfunits.Units
            #---------------------------------------------------
            # Storing scale factors and offsets may be faster.
            #---------------------------------------------------
#             p_units = p_bmi.get_var_units( long_var_name )
#             u_units = u_bmi.get_var_units( long_var_name )
#             if (u_units != p_units):
#                 #---------------------------------------------
#                 # Note: conform fails if units are the same.
#                 #---------------------------------------------
#                 values = Units.conform( values, Units(p_units), Units(u_units) )

            #-------------------------------------------
            # Call Regridder to regrid values from the
            # provider's grid to this user's grid.
            #    regridder -> grid_mapper ??
            #-------------------------------------------
            # values = self.regridder.regrid( values, user_bmi, provider_bmi)
            
            #---------------------------------------------------
            # Embed a reference to long_var_name from the
            # provider into the (BMI level of) user component.
            #---------------------------------------------------
            self.set_values( long_var_name, values, user_name )

            #------------------        
            # Optional report
            #------------------
            REPORT = False
            if (REPORT):
                print('user: ' + user_name)
                print('    just obtained var:  ' + long_var_name)
                print('    from provider: ' + provider_name)
           
    #   get_required_vars()
    #-------------------------------------------------------------------
    #  Not used currently.  Alternative to get_required_vars.
    #-------------------------------------------------------------------
#     def initialize_and_connect_comp_set( self, REPORT=False ):
# 
#         #-------------------------------------------------------------
#         # Note: This first calls find_var_users_and_providers()
#         #       to find and save self.all_input_var_names and
#         #       self.all_output_var_names.  It then calls
#         #       check_var_users_and_providers() to determine if
#         #       self.comp_set has a provider (in self.var_providers)
#         #       for every long_name in self.all_input_var_names.
#         #       If it does, then self.connect is used to create
#         #       all necessary connections (as embedded references)
#         #       between users and providers.
#         #-------------------------------------------------------------               
#         # (5/17/12) The TopoFlow components used to call a method:
#         # "initialize_required_components()" that in turn called:
#         # "initialize_ports()", inherited from CSDMS_base.py.
#         # But this is not done with the new approach.
#         #-------------------------------------------------------------
#         
#         #----------------------------------------------
#         # Find and check variable users and providers
#         #----------------------------------------------
#         self.find_var_users_and_providers()
#         OK = self.check_var_users_and_providers()
#         if not(OK):
#             return OK    # (a report was just printed)
# 
#         #------------------------
#         # For testing (5/17/12)
#         #------------------------
# ##        print '##########################################'
# ##        print "var_users[ 'channel_time_step_size' ] ="
# ##        for each in self.var_users[ 'channel_time_step_size' ]:
# ##            print '    ' + each
# ##        print '##########################################'              
#         
#         #------------------------------------------------
#         # Loop over providers in order of provider_list.
#         # Note that the dictionary, self.comp_set, is
#         # not ordered.  Order of initialize() matters.
#         #------------------------------------------------
#         for provider_name in self.provider_list:
#             bmi = self.comp_set[ provider_name ]
#             #-------------------------------------------
#             # Initialize the provider component to set
#             # all of its variables, etc.
#             #-------------------------------------------
#             ## print 'cfg_prefix =', cfg_prefix
#             cfg_file = None
#             self.initialize( provider_name, cfg_file )
#             print 'Initialized component of type: ' + provider_name + '.'
#                 
#             #--------------------------------------------------------
#             # Create a dictionary that maps every long_var_name
#             # that can be provided by this comp_set to a comp_name.
#             #--------------------------------------------------------
#             output_var_names = bmi.get_output_var_names()
#             for long_var_name in output_var_names:
#                 #-----------------------------------------------------
#                 # Is this long_var_name used by any other component?
#                 #-----------------------------------------------------
#                 try:
#                     user_list = self.var_users[ long_var_name ]
#                 except:
#                     user_list = []
#                 for user_name in user_list:
#                     #---------------------------------------------
#                     # Embed a reference in the user component to
#                     # a variable in the provider component.
#                     #---------------------------------------------
#                     self.connect( provider_name, user_name,
#                                   long_var_name, REPORT=REPORT )
# 
#         return OK
#     
#     #   initialize_and_connect_comp_set()
    #-------------------------------------------------------------------
#     def run_model_old( self, driver_comp_name='hydro_model',
#                        cfg_directory=None, cfg_prefix=None,
#                        time_interp_method='Linear'):
#         ## (rename to run_comp_set ????)
#       
#         #-------------------
#         # Default settings
#         #-------------------
#         DEBUG = True
#         ## DEBUG = False
#         if (cfg_prefix == None):
#             print 'ERROR: The "cfg_prefix" argument is required.'
#             return
#         if (cfg_directory == None):
#             print 'ERROR: The "cfg_directory" argument is required.'
#             return
#         
#         #--------------------------------------------------
#         # All components, including this one (the driver)
#         # will look in the CWD for their CFG file.
#         #--------------------------------------------------
#         # This must come after "repo" stuff, which also
#         # changes the directory.
#         #--------------------------------------------------        
#         if (cfg_directory is not None):
#             os.chdir( cfg_directory )
#         self.cfg_prefix    = cfg_prefix
#         self.cfg_directory = cfg_directory
# 
#         #-----------------------------------------------------
#         # Set self.comp_set_list and self.provider_list
#         # from info in the provider file, including the
#         # repository path "repo_path".
#         #-----------------------------------------------------
#         self.provider_file = (cfg_prefix + '_providers.txt')
#         self.read_provider_file()
#         
#         #------------------------------------
#         # Get the component repository info
#         #------------------------------------
#         ## self.repo_path = '/Users/peckhams/Dropbox/00_New_Framework/'
#         ## repo_file = self.repo_path + 'component_repository.xml'
#         ## repo_file = self.repo_dir + 'component_repository.xml'  # (6/19/13)
#         
#         self.read_repository( SILENT=False )
# ##        print 'Components in repository:'
# ##        for comp_name in f.repo_list:
# ##            print '   ' + comp_name
# ##        print ' '
#         
#         #--------------------------------------------
#         # Instantiate a complete set of components.
#         #--------------------------------------------
#         # Now the instantiate() method only allows
#         # one component of each "type" (comp_type).
#         #--------------------------------------------    
#         for comp_name in self.comp_set_list:
#             self.instantiate( comp_name, SILENT=False )
#         ### self.instantiate_all()   ### Later; change it first.
#        
#         #---------------------------------------------
#         # Try to automatically connect every user to
#         # a corresponding provider in the comp_set.
#         #---------------------------------------------------
#         # Provider components are initialized in the order
#         # of provider_list and then set references in each
#         # component that uses one or more of their vars.
#         #---------------------------------------------------
# ##        OK = self.initialize_and_connect_comp_set( REPORT=True )
# ##        if not(OK):
# ##            return
# 
#         #------------------------------------------------------------
#         # (2/18/13) Previous version of framework class initialized
#         # and then connected the components in the comp_set using
#         # embedded references.  In this version we will call a
#         # "get_required_vars()" method within the time loop, which
#         # will in turn call the get_values() and set_values()
#         # methods.  But we still need to initialize the comp set.
#         #------------------------------------------------------------
#         ## OK = self.initialize_comp_set( REPORT=True )
#         OK = self.initialize_comp_set( REPORT=False )
#         if not(OK):
#             return
#         
#         #---------------------------------------
#         # Set mode of the driver component.
#         # Note: Must happen before next block.
#         #---------------------------------------
#         driver = self.comp_set[ driver_comp_name ]
#         driver.mode = 'driver'
#         print('Driver component name = ' + driver_comp_name)
#         print()
#         
#         #-----------------------------------
#         # Initialize all time-related vars
#         #-----------------------------------
#         self.initialize_time_vars()
#         self.initialize_framework_dt()
# 
#         #------------------------------------
#         # Instantiate a "time_interpolator"
#         #------------------------------------
#         time_interpolator = time_interpolation.time_interpolator(
#                                           self.comp_set,
#                                           self.provider_list,
#                                           self.vars_provided,
#                                           time_interp_method )
#         time_interpolator.initialize()
#         #--------------------------------------------
#         # This will be used by set_provided_vars().
#         #--------------------------------------------
#         self.time_interpolator = time_interpolator
#         
#         while not(self.DONE):
# 
#             # try:
#             #-------------------------------------------------
#             # Update components that are ready to be updated
#             #----------------------------------------------------
#             # (2/18/13) The "provider_list" we loop over
#             # here refers to entries in the provider_file; they
#             # don't necessarily provide anything to another
#             # component in the comp_set.
#             #----------------------------------------------------
#             # It might be more clear to change these names:
#             #     comp_name     ->    provider_name
#             #----------------------------------------------------  
#             for comp_name in self.provider_list:
# 
#                 #--------------------------------------------------
#                 # Update time interpolation vars for every
#                 # long_var_name that is provided by this provider.
#                 # Interpolation methods = 'None', 'Linear', etc.
#                 #--------------------------------------------------
#                 # This calls bmi.update() whenever necessary.
#                 # It does so for the driver as well. (4/14/13)
#                 #--------------------------------------------------
#                 time_interpolator.update( comp_name, self.time )
#         
#                 #------------------------------------------------
#                 # (2/18/13) Use get_values()/set_values() calls
#                 # here to set latest vars from this component
#                 # into all user components that need it.
#                 #------------------------------------------------
#                 # This also calls service components as needed.
#                 #------------------------------------------------
#                 self.set_provided_vars( comp_name )
#      
#             #--------------------
#             # Are we done yet ?
#             #--------------------
#             self.DONE = (driver.DONE or self.DONE)    ####
#             self.update_time()
#             ## print 'time =', self.time
#                 
# ##            except:
# ##                print 'ERROR in run_model() method at:'
# ##                print '   time_index =', self.time_index
# ##                self.status = 'failed'
# ##                self.DONE = True
# 
#         #-------------------------
#         # Finalize the model run
#         #-------------------------
#         self.finalize_all()
#         
#     #   run_model_old()
    #-------------------------------------------------------------------
#     def set_provided_vars( self, provider_name ):    
# 
#         #----------------------------------------------------------
#         # Note:  This routine loops through all of the components
#         #        that the component given by "provider_name"
#         #        provides vars to and gets/sets the required
#         #        variables.
#         #----------------------------------------------------------
#         # Example: The variable, P, (precip rate) is computed by
#         #          the Meteorology component with a timestep of
#         #          60 seconds.  The framework uses the same
#         #          timestep as the Channels component, namely 6
#         #          seconds.  Since the Channels component uses P,
#         #          the framework sets a value of P into Channels
#         #          that has been time-interpolated to framework
#         #          time.  So every time the Channels component
#         #          gets updated, it is using a new (interpolated)
#         #          value of P.
#         #----------------------------------------------------------
#         
#         #------------------------------------------------------
#         # Restrict attention to output_var_names that are
#         # actually needed by other components.  (2/18/13)
#         # This uses "vars_provided" that was computed earlier
#         # and just once by "find_var_users_and_providers()".
#         #------------------------------------------------------
#         output_var_names = self.vars_provided[ provider_name ]
#         
#         for long_var_name in output_var_names:
# 
#             #---------------------------------------------
#             # Get a reference to long_var_name from the
#             # component with provider_name (a comp_name)
#             #---------------------------------------------           
#             # values = self.get_values( long_var_name, provider_name )
#             values = self.time_interpolator.get_values( long_var_name,
#                                                         provider_name,
#                                                         self.time )
#                 
#             #------------------------------------
#             # Break the reference for testing ?
#             #------------------------------------
#             # values = np.float64( values )
# 
#             user_list = self.var_users[ long_var_name ]
#             for user_name in user_list:
#             
#                 #------------------------------------------------------
#                 # Note that values from above are on provider's grid,
#                 # with provider's units, etc.  So need to call some
#                 # service components here, before self.set_values().
#                 #------------------------------------------------------
#                 
#                 #-------------------------------------
#                 # Call Unit Converter component here
#                 #-------------------------------------
#                 # convert to units used by this user
#                 
#                 #----------------------
#                 # Call Regridder here
#                 #----------------------
#                 # convert to grid used by this user
# 
#                 #---------------------------------------------------
#                 # Embed a reference to long_var_name from the
#                 # provider into the (BMI level of) user component.
#                 #---------------------------------------------------
#                 # Note that all values have been interpolated to
#                 # the current framework time.
#                 #---------------------------------------------------
#                 self.set_values( long_var_name, values, user_name )
# 
#                 #------------------        
#                 # Optional report
#                 #------------------
#                 REPORT = False
#                 if (REPORT):
#                     print 'provider: ' + provider_name
#                     print '    just shared var: ' + long_var_name
#                     print '    with user:       ' + user_name
# 
#         #---------------
#         # For Testing.
#         #---------------
# ##        bmi = self.comp_set[ provider_name ] 
# ##        comp_time = bmi.get_current_time()
# ##        print (provider_name + ' time ='), comp_time
#             
#     #   set_provided_vars()
    #-------------------------------------------------------------------              

    
