
#------------------------------------------------------
# S.D. Peckham
# March 25, 2009
# Speed tests for different uses of NumPy's "where"
#------------------------------------------------------

import time
import numpy
from numpy import *

def Speed_Test(n):
    a  = reshape(arange(n*n, dtype='Float64'), (n,n))
    nw = 0
    # SET_ZERO = False
    SET_ZERO = True
    
    #--------------------------
    # Standard where use test
    #--------------------------
    start_time = time.time()
    w  = where((a % 2) == 1)
    nw = w[0].size  # (fastest ?)
    ## nw = size(w)  # (This takes lots of time !! Reconfirmed: 7/26/10)
    # nw = w.size(0)  # (Doesn't work for tuples.)
    # nw = w.__len__()    # (gives tuple length)
    
    if (SET_ZERO and (nw != 0)):
        a[w] = 0.0
    run_time = (time.time() - start_time)
    print('size(w) =', nw)
    print('Standard Method run time =', run_time)
    print(' ')

    a = reshape(arange(n*n, dtype='Float64'), (n,n))
     
    #--------------------
    # Ravel method test
    #--------------------
    start_time = time.time()
    w  = where(ravel( (a % 2) == 1 ))[0]
    nw = w.size  # (These 2 methods take similar time.)
    # nw = size(w)
    
    if (SET_ZERO and (nw != 0)):
        a = a.flatten()   # (makes a copy)
        a[w] = 0.0
        a = a.reshape((n,n))
    run_time = (time.time() - start_time)
    print('size(w) =', nw)
    print('Ravel Method run time =', run_time)
    print(' ')

    a = reshape(arange(n*n, dtype='Float64'), (n,n))

    #-------------------------------------------------
    # Ravel method 2 test (50% faster than Standard)
    #-------------------------------------------------
    start_time = time.time()
    w  = where(ravel( (a % 2) == 1 ))[0]
    nw = w.size  # (These 2 methods take similar time.)
    # nw = size(w)
    
    if (SET_ZERO and (nw != 0)):
        a = ravel(a)  # (just changes "view", no copy made ?)
        a[w] = 0.0
        a = reshape(a, (n,n))
        # a = a.reshape((n,n))
    run_time = (time.time() - start_time)
    print('size(w) =', nw)
    print('Ravel Method 2 run time =', run_time)
    print(' ')

    a = reshape(arange(n*n, dtype='Float64'), (n,n))

    #---------------------
    # Ravel method 3 test
    #---------------------
    start_time = time.time()
    a  = ravel(a)  # (just changes "view", no copy made ?)
    # w  = where((a % 2) == 1 )
    w  = where(fmod(a,2) == 1 )
    nw = w[0].size
    # nw = size(w)
    
    if (SET_ZERO and (nw != 0)):
        a[w] = 0.0
        a = reshape(a, (n,n))
        # a = a.reshape((n,n))
    run_time = (time.time() - start_time)
    print('size(w) =', nw)
    print('Ravel Method 3 run time =', run_time)
    print(' ')

    a = reshape(arange(n*n, dtype='Float64'), (n,n))
    
    #----------------------
    # Flatten method test
    #----------------------
    start_time = time.time()
    w = where((a.flatten() % 2) == 1)[0]
    nw = w.size
    # nw = size(w)
    
    if (SET_ZERO and (nw != 0)):
        a = a.flatten()
        a[w] = 0.0
        a = a.reshape((n,n))
    run_time = (time.time() - start_time)
    run_time = (time.time() - start_time)
    print('size(w) =', nw)
    print('Flatten Method run time =', run_time)
    print(' ')

    
