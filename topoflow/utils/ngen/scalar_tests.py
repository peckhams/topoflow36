
#  Copyright (c) 2023, Scott D. Peckham
#     
#--------------------------------------------------------------
# Notes:  This file contains code to test various methods
#         of creating and updating scalar variables in NumPy
#         and to see what happens when a component attempts
#         to access the current value from another component.
#         It also shows methods to update values "in-place".
#
#         To run the tests:
#         >>> import scalar_tests as st
#         >>> st.ref_test('case1')
#         >>> st.ref_test('case2')
#         >>> st.ref_test('case3')
#
#         Also see: np.shares_memory().
#         Also see: utils/BMI_base.py/initialize_scalar().
#
#--------------------------------------------------------------
# Case 1 Notes.  self.x = np.array(value, dtype=dtype)
#
# If we create scalars as 0D numpy arrays, then:
#    - Return var is MUTABLE and other components with
#      a reference to it WILL see it change.
#    - Values will print without square brackets.
#    - Trying to subscript it will generate the error:
#      "IndexError: 0-d arrays can't be indexed"
#    - In-place updates will preserve the reference, e.g.:
#         x.fill(3), x += 1, x *= 2
#         Use optional "out" argument to numpy ufuncs:
#            e.g. np.sin(x, out=x)
#    - Non in-place updates will BREAK the reference, e.g.:
#         x = x+1, x = x*2, x = np.sin(x)
#    - Data type changes will BREAK the reference, and may
#      require more or less memory, e.g.:
#         x = np.float32(x)
#    - We CAN compare this type of scalar to others with
#      "<", ">", "==", etc.
#    - We can get the data type with ".dtype".
#    - The data type will be "numpy.ndarray".
#    - numpy.rank(x) will return 0.
#
#--------------------------------------------------------------
# Case 2 Notes.  self.x = np.array([value], dtype=dtype)
#
# If we create scalars as 1D numpy arrays, then:
#    - They are really 1D arrays with one element.
#    - Return var is MUTABLE and other components with
#      a reference to it WILL see it change.
#    - Values will print WITH square brackets unless we
#      ask for element 0 "[0]".
#    - In-place updates will preserve the reference, e.g.:
#         x[0]=3, x[:]=3, x[0] += 1, x[0] *= 2
#         Use optional "out" argument to numpy ufuncs:
#            e.g. np.sin(x, out=x)
#    - Non in-place updates will BREAK the reference, e.g.:
#         x = x+1, x = x*2, x = np.sin(x)
#    - Data type changes will BREAK the reference, and may
#      require more or less memory, e.g.:
#         x = np.float32(x)
#    - We CAN compare this type of scalar to others with
#      "<", ">", "==", etc.
#    - We can get the data type with ".dtype".
#    - The data type will be "numpy.ndarray".
#    - numpy.rank(x) will return 1.
#
#--------------------------------------------------------------
# Case 3 Notes.  self.x = np.float64( value )
#
# If we create scalars as typed numpy scalars, then:
#    -  Return var is now IMMUTABLE and other components
#       with a reference to it WILL NOT see it change.
#    - We CAN compare this type of scalar to others with
#      "<", ">", "==", etc.
#    - We can get the data type with ".dtype".
#    - The data type will be "numpy.float64".
#    - numpy.rank(x) will return 0.
#
#-----------------------------------------------------------------------
#
# class comp1()
#     __init__()
#     initialize()
#     update()
#     get_value()
#
# class comp2()
#     __init__()
#     set_value()
#
# get_and_set_x()
# print_vals_and_types()
# ref_broken()
# repair_ref()
# ref_test()
#
#-----------------------------------------------------------------------

import numpy as np

#-----------------------------------------------------------------------  
class comp1():
    def __init__(self):
        option = 'case1'

    #-------------------------------------------------
    def initialize(self):
    
        if (self.option == 'case1'):
            #---------------------------------
            # Case 1. 0-d numpy array.
            # Mutable, if updated carefully.
            # Prints without array brackets.
            #---------------------------------
            self.x = np.array(1.0, dtype='float64')
        elif (self.option == 'case2'):
            #--------------------------------------
            # Case 2. 1-d numpy array, 1 element.
            # Mutable, if updated carefully.
            # Prints with array brackets.
            #--------------------------------------
            self.x = np.array([1.0], dtype='float64')
        elif (self.option == 'case3'):
            #------------------------------------
            # Case 3. numpy scalar;  immutable.
            #------------------------------------           
            self.x = np.float64(1)
        
    #-------------------------------------------------
    def update(self, option='opt1'):
    
        if (option == 'opt1'):
            print('Update method: x += 1')
            self.x += 1
        elif (option == 'opt2'):
            print('Update method: x *= 2')
            self.x *= 2
        elif (option == 'opt3'):
            print('Update method: x.fill(3)')
            self.x.fill(3)
        elif (option == 'opt4'):
            print('Update method: np.sin(x, out=x)')  #########
            np.sin(self.x, out=self.x)    # in-place 
        #------------------------------------------------
        # NOTE:  The remaining options change the type
        #        of c1.x from ndarray to scalar, and
        #        so mutable to immutable variable;
        #        therefore need to call repair_ref().
        #-------------------------------------------------
        elif (option == 'opt5'):
            print("Update method: setattr(self, 'x', np.float64(3))")
            # val = 3
            # val = np.float64(3)
            val = np.array(3, dtype='float64')
            setattr(self, 'x', val) 
        elif (option == 'opt6'):
            print('Update method: x = x + 1')
            self.x = self.x + 1
        elif (option == 'opt7'):
            print('Update method: x = x * 2')
            self.x = self.x * 2
        elif (option == 'opt8'):
            print('Update method: x = np.float64(5)')
            self.x = np.float64(5)  
        elif (option == 'opt9'):
            print('Update method: x = np.float32(x)')
            self.x = np.float32( self.x )

    #-------------------------------------------------
    def get_value(self, var_name):
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

#-----------------------------------------------------------------------
class comp2():

    def __init__(self):
        # Allocate memory for desired value "x"
        self.x = np.array(0.0, dtype='float64')

    #------------------------------------------
    def set_value(self, var_name, scalar):
        setattr( self, var_name, scalar )
        
        #-----------------------------------
        # Using exec like this also works.
        #-----------------------------------
        # exec("self." + var_name + " = scalar") in globals(), locals()

#-----------------------------------------------------------------------
def get_and_set_x(c1, c2):

    #-------------------------------------
    # Copy reference to x from c1 to c2.
    #-------------------------------------
    vector = c1.get_value('x')
    #----------------------------
#     if (hasattr(c2, 'x')):   # May need this if type(c2.x) != type(c1.x)
#         delattr(c2, 'x')
    #----------------------------
    c2.set_value('x', vector)

    #--------------------------
    # These methods also work
    #--------------------------    
    # c2.x = c1.x
    # exec("c2.x = c1.x") in globals(), locals()

#-----------------------------------------------------------------------
def print_vals_and_types(c1, c2):

    print('   c1.x =', c1.x, ',', 'type(c1.x) =', type(c1.x))
    print('   c2.x =', c2.x, ',', 'type(c2.x) =', type(c2.x))
        
#-----------------------------------------------------------------------
def ref_broken(c1, c2):
    
    return (c1.x is not c2.x)

#     Note: A different object can have same value & type. 
#     wrong_value = (c1.x != c2.x)
#     wrong_type  = type(c1.x) != type(c2.x)
#     if (wrong_type):
#         print('   type(c1.x) =', type(c1.x))
#         print('   type(c2.x) =', type(c2.x))
#     broken = (wrong_value or wrong_type) 
#     return broken 

#-----------------------------------------------------------------------
def repair_ref(c1, c2):

    """Try to fix broken reference to c1.x"""
    print('   Repairing reference:')
    c1.initialize()   ######## Get back to 0-d ndarray
    get_and_set_x(c1, c2)
    
    # Another idea     
    # delattr(c2, 'x')
    # c2.x = c1.x
    
    if not(c1.x is c2.x):
        print('      ERROR: Failed to fix broken reference.')
    if (str(type(c1.x)) != "<class 'numpy.ndarray'>"):
        print('      ERROR: Failed to fix broken reference.')
    print('After repair_ref:')
    print_vals_and_types(c1, c2)
                
#-----------------------------------------------------------------------
def ref_test( option='case1' ):

    #---------------------------
    # Instantiate 2 components
    #---------------------------
    c1 = comp1()
    c1.option = option
    c2 = comp2()
    c2.option = option
    c1.initialize()
    print('Initial value of c1.x =', c1.x)
    #-------------------------------------
    # Copy reference to x from c1 to c2.
    # Notice we only call this once,
    #   outside of the 2 for loops.
    #-------------------------------------
    get_and_set_x(c1, c2)
    
    for k in range(2):
        print('########################')
        print('# TIME STEP =', k)
        print('########################')
        for j in range(9):
            option = 'opt' + str(j+1)
            c1.update( option )
            print('After c1.update() with option=' + option + ':')
            print_vals_and_types(c1, c2)
            if (ref_broken(c1,c2)):
                print('   Reference broken: c1.x is not c2.x')
                repair_ref(c1, c2)
            print()
    
#   ref_test()
#-----------------------------------------------------------------------
