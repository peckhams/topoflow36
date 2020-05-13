
#-------------------------------------------------------------------     
# Copyright (c) 2013-2016, Scott D. Peckham
#
# Apr 2013. New time interpolator class from/for emeli.py.
#
#-------------------------------------------------------------------
#
#  class time_interp_data()
#      __init__()
#      update()
#
#  class time_interpolator()
#      __init__()
#      initialize()
#      update()
#      update_all()
#      get_values()
#      convert_time_units()
#
#-------------------------------------------------------------------
import numpy as np

#-------------------------------------------------------------------
class time_interp_data():

    #--------------------------------------------------------
    # Note: This is a small "utility class".  We create an
    #       instance of this class for every long_var_name
    #       that is shared between components.
    #--------------------------------------------------------
    # Note: Additional arguments will need to be added in
    #       order to perform time interpolation by a method
    #       other than "Linear" (or a new class?).
    #--------------------------------------------------------
    def __init__( self, v1=None, t1=None, long_var_name=None ):
                 ## interp_method='Linear'):
            
        #-------------------------------------------
        # Save (v1,t1) to (v2,t2) because update()
        # first sets (v1,t1) from old (v2,t2).
        #-------------------------------------------
        self.v2 = v1
        self.t2 = t1
        self.long_var_name = long_var_name
        ## self.interp_method = interp_method

        #--------------------------------------
        # Need this, too, for in-place updates
        #--------------------------------------
        self.v1 = v1.copy()
        self.t1 = t1.copy()

        #--------------
        # For testing
        #--------------
##        if (self.long_var_name == 'atmosphere_water__precipitation_leq-volume_flux'):
##            print 'In __init__():'
##            print '(P1,P2, t1,t2) =', self.v1, self.v2, self.t1, self.t2
            
    #   __init__()
    #----------------------------------------------------------
    def update( self, v2=None, t2=None ):
                
        #---------------------------------------------------- 
        # Note: v1 and v2 could be 0D, 1D, 2D or 3D arrays.
        #       However, since they are NumPy ndarrays, the
        #       equations used below will work regardless
        #       of the array's rank.
        #----------------------------------------------------
        # Note: We can use "in-place" assignments for v1
        #       and v2 as long as their rank is > 0.
        #----------------------------------------------------
        
        #---------------------------------------------
        # Update the "start values" (old end values)
        # (in-place, if possible)
        # Note: try/except is slightly faster.
        # Note: Need to use "copy()" as shown.
        #---------------------------------------------
        self.t1 = self.t2.copy()
        try:
            self.v1[:] = self.v2.copy()
        except:
            self.v1 = self.v2.copy()
        #-----------------------------------           
##        if (np.ndim( self.v1 ) > 0):
##            self.v1[:] = self.v2.copy()
##        else:
##            self.v1 = self.v2.copy()
            
        #--------------------------
        # Update the "end values"
        # (in-place, if possible)
        # Note: Need to use "copy()" as shown.
        #---------------------------------------------
        self.t2 = t2
        try:
            self.v2[:] = v2.copy()  ## NEED THIS!
        except:
            self.v2 = v2.copy()     ## NEED THIS!
        #-----------------------------------  
##        if (np.ndim( self.v2 ) > 0):
##            self.v2[:] = v2
##        else:
##            self.v2 = v2
        
        #---------------------------------------------
        # Update the interpolation parameters, a & b
        # They are used in get_values2().
        #---------------------------------------------
        # This would also work:
        #    v1_ne_v2 = (v2 - self.v1) != 0
        #    if np.any( v1_ne_v2 ) and (t2 != self.t1):
        #----------------------------------------------------       
        dv = np.abs(v2 - self.v1)
        dv_min = dv.min()
        if (dv_min != 0) and (t2 != self.t1):
            self.a = (v2 - self.v1) / (t2 - self.t1)
            self.b = v2 - (self.a * t2)
        else:
            #------------------------------------------
            # Variables that don't vary in time will
            # have v1 = v2, but t2 > t1.
            # Disabled TopoFlow components will have
            # v2 = v1 and t2 = t1, but they may still
            # provide default values (e.g. precip=0).
            #------------------------------------------            
            # This a and b gives "no interpolation",
            # that is, v[t] = v1 = v2.
            #------------------------------------------
            self.a  = np.float64(0)
            self.b  = v2

        #--------------
        # For testing
        #--------------
##        if (self.long_var_name == 'atmosphere_water__precipitation_leq-volume_flux'):
##            print '(P1,P2, t1,t2) =', self.v1, self.v2, self.t1, self.t2
            
    #   update()
    #----------------------------------------------------------

#     time_interp_data() (class)
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
class time_interpolator():

    #----------------------------------------
    # Define some unit-conversion constants
    #----------------------------------------
    secs_per_min   = 60
    secs_per_hour  = 60  * secs_per_min
    secs_per_day   = 24  * secs_per_hour 
    secs_per_year  = 365 * secs_per_day
    secs_per_month = secs_per_year / 12    #########
    
    #---------------------------------------------------------- 
    def __init__( self, comp_set, comp_names, vars_provided,
                  method='Linear' ):

        #-------------------------------------------------------
        # Note: These are currently passed in from framework.
        #
        #       comp_set = a dictionary that takes a comp_name
        #                  key and returns a reference to a
        #                  BMI model instance.
        #
        #       comp_names = a list of all comp_names that
        #                    provide vars to other components
        #                    (and all of the keys in comp_set)
        #
        #       vars_provided = a dictionary that takes a
        #                       comp_name key and returns a
        #                       list of all the long_var_names
        #                       that the comp actually provides
        #                       to another component in the set
        #
        #       method = 'None' or 'Linear' (so far)
        #-------------------------------------------------------
        self.comp_set             = comp_set
        self.provider_comp_list   = comp_names
        self.vars_provided        = vars_provided
        self.interpolation_method = method
        print('Time interpolation method = ' + method)
        print(' ')
        
    #   __init__()
    #----------------------------------------------------------   
    def initialize( self ):

        #------------------------------------------------------------
        # Note: This function initializes a dictionary called:
        #       self.time_interp_vars and should be called from the
        #       framework's initialize() function.
        #
        #       Given "long_var_name" as a key, the dictionary
        #       returns a bundle of variables that are needed to
        #       perform time interpolation for that variable when
        #       it is requested from other components.
        #
        #       Note that self.vars_provided[ provider_comp_name ]
        #       contains a list of vars that are actually provided
        #       by that provider port to some other component in
        #       the current set of components (i.e. "comp_set").
        #-------------------------------------------------------------
        # Note: If we could somehow distinguish between provided
        #       vars that vary in time and those that don't, then
        #       we could avoid some extra work.  But this works.
        #-------------------------------------------------------------
        # Note: provider_comp_list always includes the Driver.
        #-------------------------------------------------------------        
        method = self.interpolation_method
        
        #----------------------------     
        # Case of no interpolation
        # (i.e. "steps" or "jumps")
        #----------------------------
        if (method == 'None'):
            self.time_interp_vars = None   ####
            #--------------------------------------------
            # For new method, we must call bmi.update()
            # for every provider here. (4/13/13)
            #--------------------------------------------
            for comp_name in self.provider_comp_list:
                bmi = self.comp_set[ comp_name ]                
                bmi.update( -1.0 )
##                print 'Updated component: ' + comp_name
##            print ' '
            return
        
        #-------------------------------      
        # Case of Linear interpolation
        #-------------------------------
        if (method == 'Linear'):
            self.time_interp_vars = dict()
            
            for comp_name in self.provider_comp_list:
                bmi = self.comp_set[ comp_name ]
                #---------------------------------------
                # Get t1 and convert units, if needed.
                #---------------------------------------
                comp_time_units = bmi.get_time_units()
                t1 = bmi.get_current_time()
                t1 = self.convert_time_units( t1, comp_time_units )
                    
                #---------------------------------------------------
                # Get vars at start of interpolation time interval
                #---------------------------------------------------                
                for long_var_name in self.vars_provided[ comp_name ]:
                    v1 = bmi.get_values( long_var_name )
                    data = time_interp_data( v1=v1, t1=t1, \
                                long_var_name=long_var_name )
                    self.time_interp_vars[ long_var_name ] = data
                    
                #--------------------------------------------
                # Call this component's update() just once.
                #---------------------------------------------
                # Note: Driver is updated here, too, even if
                #       it doesn't provide vars to others.
                #---------------------------------------------           
                bmi.update( -1.0 )

                #---------------------------------------
                # Get t2 and convert units, if needed.
                #---------------------------------------
                t2 = bmi.get_current_time()
                t2 = self.convert_time_units( t2, comp_time_units )

                #--------------
                # For testing
                #--------------
##                print 'Updated component: ' + comp_name
##                print '   (t1, t2) =', t1, t2
                
                #-------------------------------------------------
                # Get vars at end of interpolation time interval
                #-------------------------------------------------                
                for long_var_name in self.vars_provided[ comp_name ]:       
                    v2 = bmi.get_values( long_var_name )
                    #-------------------------------------
                    # Save (v2,t2) and update the time
                    # interpolation parameters a and b.
                    #------------------------------------- 
                    self.time_interp_vars[ long_var_name ].update(v2, t2)
            return

        #-------------------------------------      
        # Case of Cubic Spline interpolation
        #-----------------------------------------------------
        # Note: Cubic spline interpolation with a natural or
        #       clamped boundary condition requires that all
        #       time interval endpoint values are available
        #       (i.e. for an entire model run).
        #       However, during a run with Linear or None
        #       interpolation we could compute the a0 and b0
        #       that are needed to compute a[n], b[n], c[n]
        #       and d[n] for a subsequent run that uses
        #       cubic spline interpolation.
        #-----------------------------------------------------
        # Note: We need to call bmi.update() 3 times here,
        #       and then just once below in update().
        #-----------------------------------------------------        
##        if (method == 'Cubic'):
##            self.time_interp_vars = dict()
##            return
        
    #   initialize()
    #-------------------------------------------------------------------        
    def update( self, comp_name, time ):

        #------------------------------------------------------------
        # Note: This function provides automatic time-interpolation
        #       for components that have different time steps.
        #------------------------------------------------------------
        # Note: The "framework time step" will be the same as the
        #       component with the smallest time step.
        #
        #       If follows that if the time argument is "framework
        #       time", then we only need to call bmi.update() once
        #       for any component to make its internal time greater
        #       than the framework time.
        #
        #       We must make sure that this method works when called
        #       for the component (comp_name) that has the smallest
        #       time step.  In that case, we don't need to do any
        #       time interpolation and should take the "None" branch
        #       below.  #### CHECK THAT THIS HAPPENS. ####
        #------------------------------------------------------------
        # Note: A component's current time is incremented every
        #       time its bmi.update() method is called, as when
        #       done by the initialize() method.
        #------------------------------------------------------------ 
        DEBUG = False
        ## DEBUG = True
        bmi = self.comp_set[ comp_name ]

        #-----------------------------------------------------
        # Get current time of component with this comp_name.
        # Convert units to framework time units, if needed.
        #-----------------------------------------------------
        comp_time_units = bmi.get_time_units()
        comp_time       = bmi.get_current_time()
        # comp_time0      = comp_time.copy()
        comp_time = self.convert_time_units( comp_time, comp_time_units )
        #--------------
        # For testing
        #--------------
##        print 'comp_name =', comp_name
##        print 'comp_time before =', comp_time0
##        print 'comp_time after  =', comp_time
##        print ' '
        
        if (DEBUG):
            print('=============================================')
            print('In update_time_interpolation():')
            print('  time (fmwk) =', time)
            print('  comp_name   =', comp_name)
            print('  comp_time   =', comp_time)

        #--------------------------------------------
        # Do we need to update interpolation vars ?
        #------------------------------------------------
        # Note: DISABLED components have comp_time = 0.
        #------------------------------------------------
        ### if (time < comp_time):   # (This works, too.)
        if (time <= comp_time):
            if (DEBUG):
                print('  NO update for: ' + comp_name + ' interp. vars')
            return

        #------------------------------------------------
        # The current "framework time" has passed this
        # model component's internal time so we need to
        # call the model's bmi.update() method and then
        # update the time interpolation vars.
        #------------------------------------------------
        if (DEBUG):
            print('  Framework updated: ' + comp_name + ' interp. vars')

        #------------------------------------------------
        # Note: We need to check the component status
        #       here because otherwise bmi.update() is
        #       called every time below for DISABLED
        #       components.
        #------------------------------------------------
        # Using (status = 'initialized') works because
        # the initialize() method caused all other
        # components to reach "updated" status.
        #------------------------------------------------       
        comp_status = bmi.get_status()
        # if (comp_status == 'disabled'):  # (not used/ready yet)
        if (comp_status == 'initialized'):  # (this works)
            return
        
        #---------------------------     
        # Case of no interpolation
        #---------------------------
        if (self.interpolation_method == 'None'):
            #-------------------------------------------------
            # Since the framework has the smallest timestep,
            # we should only need to call bmi.update() once
            # in order to get comp_time > time.
            #-------------------------------------------------
            bmi.update( -1.0 )
            ## self.update( comp_name )  # (Checks for failure.)
            return
        
        #-------------------------------      
        # Case of Linear interpolation
        #-------------------------------
        if (self.interpolation_method == 'Linear'):
                
            #--------------------------------------------
            # Call this component's update() just once.
            #--------------------------------------------
            bmi.update( -1.0 )
            ## self.update( comp_name )  # (has error messages)

            #---------------------------------------
            # Get t2 and convert units, if needed.
            #---------------------------------------
            comp_time_units = bmi.get_time_units()
            t2 = bmi.get_current_time()
            t2 = self.convert_time_units( t2, comp_time_units )

            #---------------------------------------------------
            # Get values at end of interpolation time interval
            #--------------------------------------------------- 
            for long_var_name in self.vars_provided[ comp_name ]:
                #------------------------------------------------
                # Note: bmi.get_values() works for any rank.
                #------------------------------------------------
                # self.time_interp_vars is a dictionary that is
                # initialized in the framework's initialize().
                #------------------------------------------------ 
                v2 = bmi.get_values( long_var_name )
                
                #--------------
                # For testing
                #--------------
##                if (long_var_name == 'atmosphere_water__precipitation_leq-volume_flux'):
##                    print '(time, P) =', t2, v2
                    
                i_vars = self.time_interp_vars[ long_var_name ]
                #-------------------------------------
                # This also updates v1 and t1 first.
                #-------------------------------------
                i_vars.update(v2, t2)

                #--------------
                # For testing
                #--------------
##                print 'Updated component: ' + comp_name
##                print '   (t1, t2) =', i_vars.t1, i_vars.t2
               
            return

        #-------------------------------------      
        # Case of Cubic Spline interpolation
        #-------------------------------------  
##         if (self.interpolation_method == 'Cubic'):
        

    #   update()
    #-------------------------------------------------------------------        
    def update2( self, comp_name ):

        #------------------------------------------------------------
        # Note: This function provides automatic time-interpolation
        #       for components that have different time steps.
        #------------------------------------------------------------
        # Note: The "framework time step" will be the same as the
        #       component with the smallest time step.
        #
        #       If follows that if the time argument is "framework
        #       time", then we only need to call bmi.update() once
        #       for any component to make its internal time greater
        #       than the framework time.
        #
        #       We must make sure that this method works when called
        #       for the component (comp_name) that has the smallest
        #       time step.  In that case, we don't need to do any
        #       time interpolation and should take the "None" branch
        #       below.  #### CHECK THAT THIS HAPPENS. ####
        #------------------------------------------------------------
        # Note: A component's current time is incremented every
        #       time its bmi.update() method is called, as when
        #       done by the initialize() method.
        #------------------------------------------------------------
        # Note: In this version, we assume that bmi.update() was
        #       already called by caller of this method. (4/18/13)
        #------------------------------------------------------------
        DEBUG = False
        ## DEBUG = True
        
        #---------------------------     
        # Case of no interpolation
        #---------------------------
        if (self.interpolation_method == 'None'):
            return

        #------------------------------------------------
        # Note: We need to check the component status
        #       here because otherwise bmi.update() is
        #       called every time below for DISABLED
        #       components.
        #------------------------------------------------
        # Using (status = 'initialized') works because
        # the initialize() method caused all other
        # components to reach "updated" status.
        #------------------------------------------------
        bmi = self.comp_set[ comp_name ]  # (or pass in bmi)
        comp_status = bmi.get_status()
        # if (comp_status == 'disabled'):  # (not used/ready yet)
        if (comp_status == 'initialized'):  # (this works)
            return
 
        #-------------------------------      
        # Case of Linear interpolation
        #-------------------------------
        if (self.interpolation_method == 'Linear'):

            #---------------------------------------
            # Get t2 and convert units, if needed.
            #---------------------------------------
            comp_time_units = bmi.get_time_units()
            t2 = bmi.get_current_time()
            t2 = self.convert_time_units( t2, comp_time_units )

            #---------------------------------------------------
            # Get values at end of interpolation time interval
            #--------------------------------------------------- 
            for long_var_name in self.vars_provided[ comp_name ]:
                #------------------------------------------------
                # Note: bmi.get_values() works for any rank.
                #------------------------------------------------
                # self.time_interp_vars is a dictionary that is
                # initialized in the framework's initialize().
                #------------------------------------------------ 
                v2 = bmi.get_values( long_var_name )
                
                #--------------
                # For testing
                #--------------
##                if (long_var_name == 'atmosphere_water__precipitation_leq-volume_flux'):
##                    print '(time, P) =', t2, v2
                    
                i_vars = self.time_interp_vars[ long_var_name ]
                #-------------------------------------
                # This also updates v1 and t1 first.
                #-------------------------------------
                i_vars.update(v2, t2)

                #--------------
                # For testing
                #--------------
##                print 'Updated component: ' + comp_name
##                print '   (t1, t2) =', i_vars.t1, i_vars.t2
               
            return

        #-------------------------------------      
        # Case of Cubic Spline interpolation
        #-------------------------------------  
##         if (self.interpolation_method == 'Cubic'):
        

    #   update2()
    #-------------------------------------------------------------------        
    def update_all( self, time ):

        for comp_name in self.provider_comp_list:
            self.update( comp_name, time )

    #   update_all()
    #-------------------------------------------------------------------        
    def get_values( self, long_var_name, comp_name, time ):

        #-------------------------------------------------------
        # Note: This method returns a NumPy "ndarray" object
        #       that Babel is able to pass to other components
        #       as a SIDL generic array.
        #-------------------------------------------------------
        # Note: The update() method is called for comp_name
        #       before this is called.
        #-------------------------------------------------------
        bmi = self.comp_set[ comp_name ]  # (pass in bmi ?)

        #-------------------------------------------------------
        # Has this component been disabled?  If so, it doesn't
        # advance time or update its initial values so time
        # interpolation is not needed.
        #------------------------------------------------------------
        # TopoFlow components currently have a "comp_status"
        # attribute that is either "Enabled" or "Disabled", set in
        # their CFG file and read by BMI_base.read_config_file().
        # They also have a "status" attribute that is from the
        # OpenMI status types (e.g. "initialized", "initializing").
        # Should we add "disabled" to the OpenMI status types?
        #------------------------------------------------------------
        # If a component has already be finalized, then just get
        # its current (final) values; do not interpolate.  This
        # is needed for framework.finalize_all() to work. (8/20/13)
        #------------------------------------------------------------
        comp_status = bmi.get_status()
        if (comp_status == 'initialized') or \
           (comp_status == 'finalized'):
            return bmi.get_values( long_var_name )

        #---------------------------     
        # Case of no interpolation
        #-------------------------------------------------------
        # Note that if (time < comp_time) then we are
        # just returning the same value that all users already
        # have.  Maybe we can avoid doing this somehow.
        #-------------------------------------------------------
        if (self.interpolation_method == 'None'):

            #------------------------------------
            # For testing. Is time in interval?
            #------------------------------------
            bmi_time = bmi.get_current_time()
            if (time > bmi_time):
                print('--------------------------------------------')
                print(' WARNING: (in time_interpolation.py)')
                print('     time > bmi_time in bmi.get_values().')
                print('     time, bmi_time =', time, bmi_time)
                print('     comp_name =', comp_name)
                print('--------------------------------------------')
                print(' ')
                
            return bmi.get_values( long_var_name )

        #-------------------------------      
        # Case of Linear interpolation
        #-------------------------------
        if (self.interpolation_method == 'Linear'):
            #------------------------------------------------
            # Compute and return a time-interpolated value.
            #------------------------------------------------
            i_vars = self.time_interp_vars[ long_var_name ]

            #------------------------------------
            # For testing. Is time in interval?
            #------------------------------------           
            if (time > i_vars.t2):
                print('#######################################')
                print(' ERROR: time > t2 in get_values().')
                print('        time, t2 =', time, i_vars.t2)
                print('        comp_name =', comp_name)
                print('#######################################')
                print(' ')

            value = (i_vars.a * time) + i_vars.b

            #--------------
            # For testing
            #--------------
##            if (long_var_name == 'atmosphere_water__precipitation_leq-volume_flux'):
##                print '(time, P, a, b) =', time, value, i_vars.a, i_vars.b

            return value
            
        #-------------------------------------      
        # Case of Cubic Spline interpolation
        #-------------------------------------  
##         if (self.interpolation_method == 'Cubic'):
##            #------------------------------------------------
##            # Compute and return a time-interpolated value.
##            #------------------------------------------------
##            value =  ??????
##            return value
                
                
    #   get_values()
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
            time = in_time * secs_per_hour
        elif (in_units in ['minutes','m']):     ### month?
            time = in_time * self.secs_per_min
        else:
            time = in_time.copy()
            ## time = in_time

        return time
    
    #   convert_time_units()
    #-------------------------------------------------------------------
