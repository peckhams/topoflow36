
## Copyright (c) 2001-2009, Scott D. Peckham
## January 2009  (converted from IDL)
## October 2009

#------------------------------------------------
#  Implement a binary heap-based priority queue
#  in pure Python.  Later, compare to using the
#  built-in package "heapq".
#------------------------------------------------

from numpy import *
import numpy

import sys
import time

#-----------------------------------------------------------------------
#
#   unit_test()
#
#   class heap_base:
#       initialize()
#       insert()
#       get_min_ID()
#       get_min()
#
#-----------------------------------------------------------------------
def unit_test(nvals=100, REPORT=False):

    import heapq
    
    #----------------------
    # Initialize the heap
    #----------------------
    nx = int16(100)
    ny = int16(100)
    
    heap = heap_base()
    heap.initialize(nx, ny)
    
    if (nvals >= heap.size):    
        print('SORRY: Heap array is not big enough.')
        print('       Heap size = ' + str(heap.size))
        print('       but nvals = ' + str(nvals))
        print(' ')
        return -1
    
    #---------------------------------------
    # Generate some uniform random numbers
    #---------------------------------------
    seed = 36421
    numpy.random.seed( seed )
    ran    = -5 + (random.uniform(low=0.0, high=1.0, size=(nvals)) * 10)
    ranIDs = numpy.ceil(ran).astype('Int32')
    
    #----------------------------
    # Put random values on heap
    #----------------------------
    for k in range(nvals):
        x  = ran[k]
        ID = ranIDs[k]
        heap.insert(x, ID)

    #-----------------------------------
    # Put random values on "fast_heap"
    #--------------------------------------------------------
    # Note that "heapify()" requires argument to be a list.
    # Python's "zip()" is used to bundle values and IDs.
    #--------------------------------------------------------
    ## fast_heap = ran
    ## fast_heap = list(ran)
    fast_heap = list(zip(ran, ranIDs))
    heapq.heapify( fast_heap )
  
    #-----------------------------------
    # Retrieve values from heap, which
    # should result in a sorted array
    #-----------------------------------
    start1 = time.time()
    a = numpy.zeros([nvals], dtype='Float32')
    b = numpy.zeros([nvals], dtype='Float32')
    for k in range(nvals): 
        a[k] = heap.get_min()
        b[k] = heapq.heappop( fast_heap )[0]
    stop1 = time.time()
    #----------------------------------------------------
    # NB!  This timing isn't fair, for several reasons,
    #      such as the "for" loop required above.
    #----------------------------------------------------
    print('Time for heap sort  = ' + str(stop1 - start1))
    
    #---------------------------------
    # Now sort "ran" with IDL's SORT
    #---------------------------------
    start2 = time.time()
    ## a2 = ran[ numpy.argsort( ran ) ]  # (also works)
    a2 = numpy.sort(ran)
    a2 = numpy.float32(a2)
    stop2 = time.time()
    print('Time for quick sort = ' + str(stop2 - start2))

    #-----------------------------------------
    # Now compare the two heap-sorted arrays
    #-----------------------------------------
    w  = where(a != b)
    nw = size(w[0])
    if (nw == 0):    
        print('SUCCESS: The heap-sorted arrays are identical !!')
    else:    
        print('FAILURE: The heap-sorted arrays are not the same.')
        print('nvals  = ' + str(nvals))
        print('n_heap = ' + str(heap.size))
        print('nw     = ' + str(nw))
    print(' ')
    
    #-----------------------------
    # Now compare the two arrays
    #-----------------------------
    w = where(a != a2)
    nw = size(w[0])
    if (nw == 0):    
        print('SUCCESS: Heap-sort and quick-sort arrays are identical !!')
    else:    
        print('FAILURE: Heap-sort and quick-sort arrays are not the same.')
        print('nvals  = ' + str(nvals))
        print('n_heap = ' + str(heap.size))
        print('nw     = ' + str(nw))
    print(' ')

    #-----------------------------
    # Option to print the values
    #-----------------------------
    if (REPORT):
        print('Heap-sorted values  =')
        print(a)
        print(' ')
        print('"Fast heap"-sorted values  =')
        print(b)
        print(' ')
        print('Quick-sorted values =')
        print(a2)
        print(' ')
        
#   unit_test()
#-----------------------------------------------------------------------
class heap_base():

    #-------------------------------------------------------------------    
    def initialize(self, nx, ny, dtype='Float32'):

        self.USE_64_BITS = False

        if (self.USE_64_BITS):
            self.ID_type='Int64'
        else:
            self.ID_type='Int32'

        #-------------------------------------------
        # nheap is based on the number of boundary
        # pixels and a factor of safety
        #-------------------------------------------
        self.size = int64(500) * (nx + ny) + 1
        self.heap = numpy.zeros(self.size, dtype=dtype)
        self.IDs  = numpy.zeros(self.size, dtype=self.ID_type)
        
        #---------------------------
        # Type is important here !
        # Shouldn't be ULONG64.
        #---------------------------
        self.n      = numpy.int64(0)
        self.nmax   = numpy.int64(0)
        
    #   initialize()
    #------------------------------------------------------------------- 
    def insert(self, x, ID):
        
        #------------------------------------------------------
        # Notes:  Heap array is large and n gives the size of
        #         the part that is currently in use.  First
        #         value with index of 0 is not used, so the
        #         index (n + 1) points to first unused slot,
        #         which represents the bottom of the heap.
        #------------------------------------------------------
        #         i inherits data type from n (ULONG64)
        #------------------------------------------------------
        # Error trap here with CATCH had a time cost.
        # For n = 100,000 values, time dropped from
        # 0.531 to 0.358 seconds for each call.
        #------------------------------------------------------
        
        #----------------------------------------------
        # If we try to add a NaN value to the heap
        # then we'll get an infinite loop, because
        # (NaN) / 2 = NaN and  (NaN gt x) eq 0b.
        #----------------------------------------------
        # If heap[0]=0 and (x lt 0), then we get
        # an infinite while loop if we're not careful
        #----------------------------------------------
        if (numpy.isfinite(x) == 0):    
            return -1
            ## sys.exit()
        
        #-------------------------------------
        # This test actually costs some time
        # For inserting 100,000 values, time
        # decreases from 0.609 to 0.531.
        #-------------------------------------
        # Is heap array big enough ?
        #-----------------------------
        #if ((self.n + 1) > (self.size - 1)):
        #    print 'SORRY:  Heap array is not big enough.'
        #    print '        Heap size = ' + str(self.size)
        #    print '        but n = ' + str(n)
        #    sys.exit()
        
        #------------------------------------
        # Expand active range of heap array
        # to make room for the new value
        #------------------------------------
        # Store largest n so far in nmax
        #  (even this has a time cost)
        #------------------------------------
        self.n += 1
        self.nmax = numpy.maximum(self.nmax, self.n)

        #---------------------------
        # Insert value on the heap
        #--------------------------------------------------
        # Recall that heap[0] holds an initial value of 0
        # but is not supposed to be used, so we need:
        # "AND (i gt 1)".  Without this, we can get an
        # infinite loop, as in this example:
        #    Insert  1: heap = [0, 1, 0, 0, ...]
        #    Insert -1: Inf. loop, since heap[0] gt -1
        # We won't see this problem if all values gt 0.
        # If we use MAX_REPS to escape loop, then our
        # heap ends up looking like [-1, 0, 1, 0, 0,...].
        #--------------------------------------------------
        i = self.n
        #** n_reps = numpy.int64(0)
        #** MAX_REPS = numpy.int64(1000)
        #---------------------------------------------------
        #** while (heap[i/2] > x) AND (n_reps < MAX_REPS):
        #---------------------------------------------------
        #** while (heap[i/2] g> x):
        #---------------------------------------------------
        while ((self.heap[i / 2] > x) and (i > 1)):
            self.heap[i] = self.heap[i / 2]
            self.IDs[i]  = self.IDs[i / 2]
            i = i / 2
            ### n_reps += 1
        
        #--------------
        # For testing
        #--------------
        #if (n_reps >= MAX_REPS):
        #    nstr = str(n_reps)
        #    istr = '(' + str(i) + ', '
        #    hstr = str(heap[i]) + ')'
        #--------------------------------
        #    print 'ERROR:  Processing aborted'
        #    print ' '
        #    print 'Infinite loop encountered in heap.insert().'
        #    print 'Attempt to add value of:  ' + str(x)
        #    print ' '
        #    print '(i, heap[i]) = ' + istr + hstr
        #    print 'Loop counter = ' + nstr
        #    print ' '
        #    sys.exit()
        
        self.heap[i] = x
        self.IDs[i]  = ID
        
    #   insert()
    #-------------------------------------------------------------------
    def get_min_ID(self):

        #------------------------------------------------
        # Notes:  Separate function on 7/2/09 for I2PY.
        #         This must be called before get_min().
        #------------------------------------------------
        return self.IDs[1]
        # self.min_ID = self.IDs[1]
        
    #   get_min_ID()
    #-------------------------------------------------------------------
    def get_min(self):

        #--------------------------------------------------------
        #   Note:  Repeated calls to Heap_Get_Min pull values
        #          off the heap from smallest to largest.
        #          Note, however, that values in heap array
        #          itself are not in sorted order.  Also, any
        #          values past the nth slot in heap array (out
        #          of the active range) are ignored.

        #          **********************************
        #          Assume HAS_ID is always true now
        #          **********************************
        #          IDs is an optional array of DEM locations
        #          that are paired with the values in heap.
        #          IDs are also rearranged so that they remain
        #          paired to their values.
        #-------------------------------------------------------

        #----------------------------------
        # Heap array (queue) is now empty
        #----------------------------------
        #if (n == 0):
        #    print 'Priority queue is now empty.'
        #    print ' '
        #    return -9999

        #-------------------------------------------
        # Retrieve minimum value and decrease size
        # of active part of heap array by one
        #-------------------------------------------
        # NB! Last value in heap array is not max.
        #-------------------------------------------
        heap_min  = self.heap[1]
        heap_last = self.heap[self.n]
        last_ID   = self.IDs[self.n]
        self.n    = (self.n - 1)
        #*** min_ID    = self.IDs[1]   # (see get_min_ID())
        
        #---------------------------------------
        # Rearrange elements in the heap array
        # with new minimum at the top
        #---------------------------------------
        # NB!  Type of i is very important !!
        #---------------------------------------
        i = int64(1)
        while ((i * 2) <= self.n):
            #------------------------------------
            # Find a smaller child, at position
            # (2*i) or (2*i)+1
            #------------------------------------
            child = i * 2
            if (child != self.n):    
                if (self.heap[child + 1] < self.heap[child]):    
                    child += 1
            
            #----------------------
            # Percolate one level
            #----------------------
            if (heap_last > self.heap[child]):    
                #-----------------------------------
                # Move child up to parent position
                #-----------------------------------
                self.heap[i] = self.heap[child]
                self.IDs[i]  = self.IDs[child]
            else:    
                break
            
            i = child
        
        self.heap[i] = heap_last
        self.IDs[i]  = last_ID
        
        return heap_min
        ## self.heap_min = heap_min
    
    #   get_min()
    #-------------------------------------------------------------------


