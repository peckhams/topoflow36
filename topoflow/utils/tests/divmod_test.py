
import numpy
import time

def divmod_test(n):

    IDs = numpy.arange(n*n).reshape(n,n)

    start1    = time.time()
    rc_tuple1 = divmod(IDs, n)
    stop1     = time.time()

    start2    = time.time()
    rc_tuple2 = (IDs / n, IDs % n)
    stop2     = time.time()

    start3    = time.time()
    rc_tuple3 = (numpy.divide(IDs,n), numpy.mod(IDs,n))
    stop3     = time.time()
    
    print('Run time for divmod (builtin)    =', (stop1 - start1))
    print('Run time for manual method       =', (stop2 - start2))
    print('Run time for (np.divide, np.mod) =', (stop3 - start3))    
    
