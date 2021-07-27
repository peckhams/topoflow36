
#  Copyright (c) 2021, Scott D. Peckham
#  May 2021.  Started on May 17 for Localized Conflict Modeling.
#             First, read GPW population count data as GeoTIFF.
#             Write-up and theoretical results in late May.
#-------------------------------------------------------------------
#
#  conda activate tf36
#  python
#  from topoflow.utils import conflict
#  pop_grid = conflict.read_geotiff()  
#
#-------------------------------------------------------------------
#
#  test1()
#  test2()
#
#  class conflict()
#      initialize()
#      initialize_U()
#      initialize_C1()
#      initialize_C2()
#-----------------------
#      update()
#      update_U()
#      update_C1()
#      update_C2()
#-----------------------
#      update_p()
#      update_S()
#      update_S1()
#      update_S2()       # (obsolete soon?)
#      update_S_old()    # (obsolete soon?)
#      update_time()
#--------------------------------------
#      get_neighbor_cols_and_rows()
#      get_neighbor_values()
#      spread_conflicts1()
#      spread_conflicts2()
#      finalize()
#      run_model()
#
#  get_raster_cellsize()
#  get_raster_bounds()
#  bounds_disjoint()
#  read_geotiff()     # can also create RTG and RTI files
#  regrid_geotiff()
#
#  read_acled_data()
#
#-------------------------------------------------------------------
import numpy as np
import numpy.random as rn
#### import random as rn
import pandas as pd
import time

# For read_geotiff(), etc.
import gdal, osr  ## ogr
import glob, sys
import os, os.path
from . import rti_files
from . import rtg_files

#-------------------------------------------------------------------
def test1():

    cfg_file = 'conflict.cfg'
    c = conflict()
    c.run_model( cfg_file )

#   test1()
#-------------------------------------------------------------------
def test2():

    pop_grid = read_geotiff()

#   test2()
#-------------------------------------------------------------------
def test3( SUBSAMPLE=False ):

    #--------------------------------
    # Use Horn of Africa as a test.
    #--------------------------------
    in_dir   = '/Users/peckhams/Conflict/Data/GPW-v4/'
    in_file  = 'gpw_v4_population_count_rev11_2020_30_sec.tif'
    in_file  = in_dir + in_file
    # Bounds = [ minlon, minlat, maxlon, maxlat ]
    out_bounds = [ 25.0, -5.0, 55.0, 25.0]
            
    if not(SUBSAMPLE):
        #--------------------------------------
        # 3600 cols x 3600 rows, 30 arcseconds
        #--------------------------------------
        out_file = 'Horn_of_Africa_GPW-v4_pop_count_2020_30sec.tif'
        out_file = in_dir + out_file
        out_xres_sec = None    # will default to in_xres_sec
        out_yres_sec = None    # will default to in_yres_sec
        print('Reading & clipping GeoTIFF file...')
    else:
        #--------------------------------------
        # 360 cols x 360 rows, 300 arcseconds
        #--------------------------------------
        out_file = 'Horn_of_Africa_GPW-v4_pop_count_2020_450sec.tif'
        out_file = in_dir + out_file
        out_xres_sec = 450.0  # (15 times lower resolution)
        out_yres_sec = 450.0  # (15 times lower resolution)
        print('Reading, clipping & subsampling GeoTIFF file...') 
        
        #--------------------------------------
        # 360 cols x 360 rows, 300 arcseconds
        #--------------------------------------
#         out_file = 'Horn_of_Africa_GPW-v4_pop_count_2020_300sec.tif'
#         out_file = in_dir + out_file
#         out_xres_sec = 300.0  # (10 times lower resolution)
#         out_yres_sec = 300.0  # (10 times lower resolution)
#         print('Reading, clipping & subsampling GeoTIFF file...')   

    regrid_geotiff(in_file=in_file, out_file=out_file, 
           out_bounds=out_bounds,
           out_xres_sec=out_xres_sec,
           out_yres_sec=out_yres_sec,
           RESAMPLE_ALGO='bilinear', REPORT=True)

#   test3()                   
#-------------------------------------------------------------------
class conflict():
    #---------------------------------------------------------------
    def initialize( self, cfg_file=None ):
    
        home_dir   = os.path.expanduser('~') + os.sep
        
        if (cfg_file is None):
            self.in_dir   = home_dir + 'Conflict/Data/GPW-v4/'
            self.out_dir  = home_dir + 'Conflict/Output/' 
            cfg_file      = self.in_dir  + 'conflict.cfg'
            self.out_file = self.out_dir + 'conflicts.rts'
            self.IDs_file = self.out_dir + 'conflict_IDs.rts'
            self.C1_file  = ''
            self.C2_file  = ''
            #---------------------------
            # Was good for pop count U
            #---------------------------
#             self.U_file   = 'Horn_of_Africa_GPW-v4_pop_count_2020_300sec.tif'
#             self.nx       = 360
#             self.ny       = 360
#             self.c_emerge = 0.01  # (must be in (0,1])
#             self.c_spread = 0.1
#             ## self.c_spread = 0.03
#             ## self.c_spread = 0.05
#             ## self.p_geom   = 0.2
#             self.p_geom   = 0.4
            #--------------------------
            # Is good for pop count U
            #--------------------------
#             self.U_file   = 'Horn_of_Africa_GPW-v4_pop_count_2020_450sec.tif'
#             self.nx       = 240
#             self.ny       = 240
#             self.c_emerge = 0.5  # (must be in (0,1])
#             ## self.c_emerge = 0.1  # (must be in (0,1])
#             ## self.c_emerge = 0.01  # (must be in (0,1])
#             self.c_spread = 0.5
#             ## self.c_spread = 0.1
#             self.p_resolve = 0.4
#             self.p_geom    = 0.4   # (not used now)
            #-------------------------
            # Was good for uniform U
            #-------------------------
            self.U_file   = ''     # (To use uniform U)
            self.nx       = 240
            self.ny       = 240
            self.c_emerge = 0.001     ####
            ## self.c_emerge = 0.2
            ## self.c_emerge = 0.001  # (must be in (0,1])
            ## self.c_spread = 0.1  ####
            self.c_spread = 0.4  ####
            ## self.c_spread = 0.03
            ## self.c_spread = 0.05
            ## self.p_geom   = 0.2
            self.p_resolve = 0.4
            self.p_geom    = 0.4
            self.spread_method = 1
            #--------------------------
            self.time_lag = 1    # (not used yet)
            self.n_steps  = 100
            self.REPORT   = True
        else:
            #-----------------------------------
            # Read params from the config file
            #-----------------------------------
            dum = 0
    
        self.cfg_file     = cfg_file
        self.time_index   = 0
        self.n_conflict_cells = 0
        self.grid_shape  = (self.ny, self.nx)
        ## self.start_time  = time.time()
        self.start_ID    = 1
        self.start_index = 0

        #----------------------------
        # Change to input directory
        #----------------------------
        os.chdir( self.in_dir )
     
        #-----------------------------
        # Open output files to write
        #-----------------------------
        self.out_unit = open( self.out_file, 'wb')
        self.IDs_unit = open( self.IDs_file, 'wb')

        #--------------------------------------       
        # Make grids with col and row numbers
        #--------------------------------------
        cols   = np.arange( self.nx )
        rows   = np.arange( self.ny )
        cg, rg = np.meshgrid( cols, rows )
        self.col_grid = cg
        self.row_grid = rg
            
        self.initialize_U()
        self.initialize_C1()
        self.initialize_C2()
        
        #------------------------------------------------------    
        # Initialize to no conflicts
        # S will later contain 1s in grid cells with conflict
        # Initialize durations to zero also.
        # IDs will contain a unique ID for each conflict.
        # Using 'float32' for IDs now for viewing the RTS.
        #------------------------------------------------------
        self.S    = np.zeros( self.grid_shape, dtype='uint8' )
        self.durs = np.zeros( self.grid_shape, dtype='uint32')
        self.IDs  = np.zeros( self.grid_shape, dtype='float32')
        #----------------------------------------------------------
        # Create a set of random integer IDs, without replacement
        # so when we colorize, it will look better.
        #----------------------------------------------------------
        self.ran_IDs = rn.choice( 10000000, 10000000, replace=False)
        # This next method used built-in random & and problems.
        ### self.ran_IDs = rn.sample( range(1000000), 500000)
        self.start_time  = time.time()
        
    #   initialize()
    #---------------------------------------------------------------
    def initialize_U( self ):

        #-----------------------------------    
        # Start with U = a population grid
        #-----------------------------------
        if (self.U_file != ''):
            self.U = read_geotiff(in_file=self.U_file,
                                  REPORT=True)
            # In case of negative nodata value
            np.maximum(self.U, 0.0, self.U)   # (in place)
        else:        
            #---------------------    
            # Use a grid of ones
            #---------------------
            self.U = np.ones( self.grid_shape, dtype='float32' )
 
        #-----------------------------------
        # Disallow conflict on the 4 edges
        #-----------------------------------
        self.U[0,:]           = 0.0
        self.U[self.ny - 1,:] = 0.0
        self.U[:,0]           = 0.0
        self.U[:,self.nx - 1] = 0.0
                      
    #   initialize_U()
    #---------------------------------------------------------------
    def initialize_C1( self ):

        if (self.C1_file != ''):
            self.C1 = read_geotiff(in_file=self.C1_file,
                                   REPORT=True)
        else:        
            #---------------------    
            # Use a grid of ones
            #---------------------
            self.C1 = np.ones( self.grid_shape, dtype='float32' )
        
    #   initialize_C1()
    #---------------------------------------------------------------
    def initialize_C2( self ):
    
        if (self.C2_file != ''):
            self.C2 = read_geotiff(in_file=self.C2_file,
                                   REPORT=True)
        else:        
            #---------------------    
            # Use a grid of ones
            #---------------------
            self.C2 = np.ones( self.grid_shape, dtype='float32' )
        
    #   initialize_C2()
    #---------------------------------------------------------------
    def update( self ):

        self.update_U()
        self.update_C1()
        self.update_C2()
        #-------------------  
        self.update_p()
        self.update_S()      # also updates IDs
        ## self.update_S1()  # same speed as update_S()
        ## self.update_S2()
        self.update_time()

    #   update()
    #---------------------------------------------------------------
    def update_U( self ):
    
        pass

    #   update_U()
    #---------------------------------------------------------------
    def update_C1( self ):
    
        pass

    #   update_C1()
    #---------------------------------------------------------------
    def update_C2( self ):
    
        pass

    #   update_C2()    
    #---------------------------------------------------------------
    def update_p( self ):

        #----------------------------------------------------------
        # Note:  p is the probability that a conflict emerges in
        #        a grid cell, and is a function of the unrest, U.
        #        In order for p to be a probability, in (0,1],
        #        we need 0 < c_emerge <= 1.
        #----------------------------------------------------------
        self.p_emerge = (self.c_emerge / self.U.max()) * self.U
        
    #   update_p()
    #---------------------------------------------------------------
    def update_S( self ):

        #-----------------------------------------------------------
        # Note:  The previous version of this method generated
        #        Geometric random variables to model conflict
        #        durations.  This new version does not track
        #        durations explicitly, but should produce exactly
        #        the same result.  Here, any conflict ends in the
        #        kth time interval with fixed probability, p.
        #        This is modeled with a Bernoulli random variable.
        #-----------------------------------------------------------
        # Note:  A Bernoulli random variable takes the value 1
        #        with probability p and 0 with probability (1-p).
        #        It is a special case of a Binomial r.v. (n=1).
        #        np.random.binomial() allows p to be an array.
        #----------------------------------------------------------- 
                
        #--------------------------------------------------
        # Initiate new conflicts in cells with no conflict
        #-----------------------------------------------------
        # Generate Bernoulli random variables with parameter
        # p_emerge, and initiate conflicts where B1=1. 
        #-----------------------------------------------------  
        # Convert b from dtype='int64' to dtype='uint8' ?
        # This requires 8 times less memory.
        #-----------------------------------------------------
        B1 = np.random.binomial(1, self.p_emerge)
        w2 = np.logical_and( self.S == 0, B1 == 1 )
        n2 = w2.sum()

        #------------------------------------------------
        # Resolve some conflicts in cells with conflict
        #-----------------------------------------------------        
        # Generate Bernoulli random variables with parameter
        # p_resolve and terminate conflicts that get B2=1.
        # Conflict durations will then turn out to be
        # Geometric random variables, same parameter.
        #-----------------------------------------------------
        B2 = np.random.binomial(1, self.p_resolve, size=self.grid_shape)
        w3 = np.logical_and( self.S == 1, B2 == 1 )
        n3 = w3.sum()        
 
        #------------------------------------       
        # Perform the required updates to S
        #------------------------------------
        i = self.start_index
        self.S[ w2 ]   = B1[ w2 ]
        self.IDs[ w2 ] = self.ran_IDs[i:i + n2]
        self.start_index += n2
        self.n_conflict_cells += n2
        #---------------------------------------------
        self.S[ w3 ] = (1 - B2[ w3 ])
        self.n_conflict_cells -= n3        
         # Reset IDs to zero where resolved (i.e. B2 = 1).
        self.IDs[ w3 ] = self.IDs[w3] * (1 - B2[w3])        
        #---------------------------------------------       
        if (self.REPORT):
            print('time_index =', self.time_index)
            print('Number of new conflict cells =', n2)
            print('Number of resolved conflicts =', n3)
            if (self.spread_method == 0):   
                print()
 
        #----------------------------------   
        # Attempt to spread the conflicts
        #------------------------------------------------
        # Set spread_method == 0 to turn off spreading,
        # e.g. to test against theoretical results.
        #------------------------------------------------
        if (self.spread_method == 1):
            self.spread_conflicts1()
#         elif (self.spread_method == 2):
#             self.spread_conflicts2()
#         elif (self.spread_method == 3):
#             self.spread_conflicts3()
#         else:
#             pass
         
        SAVE_S = True
        if (SAVE_S):
            #---------------------------------
            # Write grid as binary to file
            # (could use .astype('float32'))
            #---------------------------------
            S2 = np.float32(self.S)
            S2.tofile( self.out_unit )
   
        SAVE_IDs = True
        if (SAVE_IDs):
            self.IDs.tofile( self.IDs_unit )
   
    #   update_S()
    #---------------------------------------------------------------
    def update_S1( self ):

        #-----------------------------------------------------------
        # Note:  The previous version of this method generated
        #        Geometric random variables to model conflict
        #        durations.  This new version does not track
        #        durations explicitly, but should produce exactly
        #        the same result.  Here, any conflict ends in the
        #        kth time interval with fixed probability, p.
        #        This is modeled with a Bernoulli random variable.
        #-----------------------------------------------------------
        # Note:  A Bernoulli random variable takes the value 1
        #        with probability p and 0 with probability (1-p).
        #        It is a special case of a Binomial r.v. (n=1).
        #        np.random.binomial() allows p to be an array.
        #----------------------------------------------------------- 

        #----------------------------------------------
        # Make a copy of self.S, i.e. S(k) vs. S(k+1)
        #----------------------------------------------
        # Otherwise, conflicts may be resolved in the
        # same time step as they were initiated.
        # Could also apply self.S updates at end.
        #----------------------------------------------
        S = self.S.copy()   # This may be costly
                
        #--------------------------------------------------
        # Initiate new conflicts in cells with no conflict
        #-----------------------------------------------------
        # Generate Bernoulli random variables with parameter
        # p_emerge, and initiate conflicts where B1=1. 
        #-----------------------------------------------------  
        # Convert b from dtype='int64' to dtype='uint8' ?
        # This requires 8 times less memory.
        #-----------------------------------------------------
        B1 = np.random.binomial(1, self.p_emerge)
        w2 = np.logical_and( S == 0, B1 == 1 )
        n2 = w2.sum()
        i = self.start_index
        self.S[ w2 ]   = B1[ w2 ]
        self.IDs[ w2 ] = self.ran_IDs[i:i + n2]
        self.start_index += n2
        self.n_conflict_cells += n2

        #------------------------------------------------
        # Resolve some conflicts in cells with conflict
        #-----------------------------------------------------        
        # Generate Bernoulli random variables with parameter
        # p_resolve and terminate conflicts that get B2=1.
        # Conflict durations will then turn out to be
        # Geometric random variables, same parameter.
        #-----------------------------------------------------
        B2 = np.random.binomial(1, self.p_resolve, size=self.grid_shape)
        w3 = np.logical_and( S == 1, B2 == 1 )
        n3 = w3.sum()
        self.S[ w3 ] = (1 - B2[ w3 ])
        self.n_conflict_cells -= n3        
         # Reset IDs to zero where resolved (i.e. B2 = 1).
        self.IDs[ w3 ] = self.IDs[w3] * (1 - B2[w3])
 
        if (self.REPORT):
            print('time_index =', self.time_index)
            print('Number of new conflicts =', n2)             
            print('Number of resolved conflicts =', n3)   
            print()

        #----------------------------------   
        # Attempt to spread the conflicts
        #------------------------------------------------
        # Set spread_method == 0 to turn off spreading,
        # e.g. to test against theoretical results.
        #------------------------------------------------
        if (self.spread_method == 1):
            self.spread_conflicts1()
#         elif (self.spread_method == 2):
#             self.spread_conflicts2()
#         elif (self.spread_method == 3):
#             self.spread_conflicts3()
#         else:
#             pass
         
        SAVE_S = True
        if (SAVE_S):
            #---------------------------------
            # Write grid as binary to file
            # (could use .astype('float32'))
            #---------------------------------
            S2 = np.float32(self.S)
            S2.tofile( self.out_unit )
   
        SAVE_IDs = True
        if (SAVE_IDs):
            self.IDs.tofile( self.IDs_unit )
   
    #   update_S1()
    #---------------------------------------------------------------
    def update_S2( self ):

        #-----------------------------------------------------------
        # Note:  The previous version of this method generated
        #        Geometric random variables to model conflict
        #        durations.  This new version does not track
        #        durations explicitly, but should produce exactly
        #        the same result.  Here, any conflict ends in the
        #        kth time interval with fixed probability, p.
        #        This is modeled with a Bernoulli random variable.
        #-----------------------------------------------------------
        # Note:  A Bernoulli random variable takes the value 1
        #        with probability p and 0 with probability (1-p).
        #        It is a special case of a Binomial r.v. (n=1).
        #        np.random.binomial() allows p to be an array.
        #----------------------------------------------------------- 

        #---------------------------------------         
        # Find cells with and without conflict
        #---------------------------------------
        w1 = (self.S == 1)
        w0 = np.invert( w1 )
        n1 = w1.sum()
        n0 = w0.sum()
           
        ## S = self.S.copy()  ########
             
        #--------------------------------------------------
        # Initiate new conflicts in cells with no conflict
        #--------------------------------------------------    
        # Convert b from dtype='int64' to dtype='uint8' ?
        # This requires 8 times less memory.
        #--------------------------------------------------
        B1 = np.random.binomial(1, self.p_emerge[w0])
        w2 = (B1 == 1)
        n2 = w2.sum()
        i = self.start_index
        self.S.flat[ w2 ]   = B1[w2]
        self.IDs.flat[ w2 ] = self.ran_IDs[i:i + n2]
        self.start_index += n2
        self.n_conflict_cells += n2
       
        if (self.REPORT):
            print('Number of new conflicts =', n2)

        #------------------------------   
        # Update S with new conflicts
        #------------------------------
#         if (n3 > 0):
#             i = self.start_index
#             self.S[ w3 ]   = 1          # (for method 1)
#             self.IDs[ w3 ] = self.ran_IDs[i:i + n3]
#             #---------------------------------------------
#             #### self.S[ w0 ]   = B1   # (for method 2)
#             #### self.IDs[ w0 ] = self.ran_IDs[i:i + n3]
#             #---------------------------------------------
#             self.start_index += n3
#             self.n_conflict_cells += n3

        #------------------------------------------------
        # Resolve some conflicts in cells with conflict
        #-----------------------------------------------------        
        # Generate Bernoulli random variables, and terminate
        # conflicts that get B2=1.  Conflict durations will
        # then turn out to be Geometric random variables.
        #-----------------------------------------------------
        B2 = np.random.binomial(1, self.p_resolve, size=n1)
        w3 = (B2 == 1)
        n3 = w3.sum()
        self.S.flat[ w3 ] = (1 - B2[ w3 ])
        self.n_conflict_cells -= n3        
         # Reset IDs to zero where resolved (i.e. B2 = 1).
        self.IDs.flat[ w3 ] *= (1 - B2[w3])
 
        #-----------------------------------    
        # Update S with resolved conflicts
        #-----------------------------------                
        if (self.REPORT):
            print('time_index =', self.time_index)
            print('Number of resolved conflicts =', n3)   
            
        #----------------------------------   
        # Attempt to spread the conflicts
        #------------------------------------------------
        # Set spread_method == 0 to turn off spreading,
        # e.g. to test against theoretical results.
        #------------------------------------------------
        if (self.spread_method == 1):
            self.spread_conflicts1()
#         elif (self.spread_method == 2):
#             self.spread_conflicts2()
#         elif (self.spread_method == 3):
#             self.spread_conflicts3()
#         else:
#             pass
         
        SAVE_S = True
        if (SAVE_S):
            #---------------------------------
            # Write grid as binary to file
            # (could use .astype('float32'))
            #---------------------------------
            S2 = np.float32(self.S)
            S2.tofile( self.out_unit )
   
        SAVE_IDs = True
        if (SAVE_IDs):
            self.IDs.tofile( self.IDs_unit )
   
    #   update_S2()
    #---------------------------------------------------------------
    def update_S_old( self ):

        #-----------------------------------------------------------
        # Notes:  A Bernoulli random variable takes the value 1
        #         with probability p and 0 with probability (1-p).
        #         It is a special case of a Binomial r.v. (n=1).
        #         np.random.binomial() allows p to be an array.
        #----------------------------------------------------------- 

        #------------------------------------------         
        # Reduce the existing durations by 1
        # Durations are integers (# of timesteps)
        #------------------------------------------
        w1 = (self.S == 1)
        # n1 = w1.sum()
        self.durs[ w1 ] -= 1

        #---------------------------------------------
        # Have any conflicts reached their duration?
        # If so, set their S value to 0.
        #---------------------------------------------
        # S[w1][w2] works for retrieving values, but
        # it doesn't work for assignments.
        #---------------------------------------------
#         w2 = (self.durs[ w1 ] == 0)  # can be empty
#         self.S[ w1 ][ w2 ] = 0  # conflict is over


        #-----------------------------------------------------
        # METHOD 1:  Works, but seems no faster than METHOD 2
        #-----------------------------------------------------
#         w2 = (self.durs[ w1 ] == 0)  # can be empty
#         self.S[ w1 ] = (1 - w2.astype('uint8'))
        #----------------------------------------------------  
        # METHOD 2
        #----------------------------------------------------        
        w2 = np.logical_and( (self.S == 1),(self.durs == 0) )
        self.S[ w2 ]   = 0
        self.IDs[ w2 ] = 0   # reset the IDs
        n2 = w2.sum()
        self.n_conflict_cells -= n2
                        
        if (self.REPORT):
            print('time_index =', self.time_index)
            print('Number of resolved conflicts =', n2)

        #--------------------------------------------------
        # Initiate new conflicts;  inherit size from p
        #--------------------------------------------------    
        # Convert b from dtype='int64' to dtype='uint8' ?
        # This requires 8 times less memory.
        #--------------------------------------------------
        dS = np.random.binomial(1, self.p_emerge)
        dS = dS.astype('uint8')
        w3 = np.logical_and(self.S == 0, dS == 1)
        #-------------------------------------------- 
        # This would allow existing conflicts to be
        # "reseeded" and causes error in counting.
        #--------------------------------------------    
        ### w3 = (dS == 1)
        n3 = w3.sum()
        if (self.REPORT):
            print('Number of new conflicts =', n3)

        if (n3 > 0):
            #------------------------------   
            # Update S with new conflicts
            #------------------------------
            self.S[ w3 ]   = 1
            ## self.IDs[ w3 ] = np.arange( n3 ) + self.start_ID
            ## self.start_ID += n3
            i = self.start_index
            self.IDs[ w3 ] = self.ran_IDs[i:i + n3]
            self.start_index += n3

            ### np.maximum( self.S, dS, self.S)    # in place
            #------------------------------------------         
            # New durations are Geometric random vars
            #------------------------------------------           
            g  = np.random.geometric( self.p_geom, size=n3 )
            self.durs[ w3 ] = g
            self.n_conflict_cells += n3

        #----------------------------------   
        # Attempt to spread the conflicts
        #------------------------------------------------
        # Set spread_method == 0 to turn off spreading,
        # e.g. to test against theoretical results.
        #------------------------------------------------
        if (self.spread_method == 1):
            self.spread_conflicts1()
#         elif (self.spread_method == 2):
#             self.spread_conflicts2()
#         elif (self.spread_method == 3):
#             self.spread_conflicts3()
#         else:
#             pass
         
        SAVE_S = True
        if (SAVE_S):
            #---------------------------------
            # Write grid as binary to file
            # (could use .astype('float32'))
            #---------------------------------
            S2 = np.float32(self.S)
            S2.tofile( self.out_unit )
   
        SAVE_IDs = True
        if (SAVE_IDs):
            self.IDs.tofile( self.IDs_unit )
   
    #   update_S_old()
    #---------------------------------------------------------------
    def update_time( self ):

        self.time_index += 1
        
    #   update_time()
    #---------------------------------------------------------------    
    def get_neighbor_cols_and_rows( self, w1, n1 ):

        cols = self.col_grid[ w1 ]
        rows = self.row_grid[ w1 ]
        #--------------------------------------------------
        # 1st index is over grid cells that have conflict.
        # 2nd index is over the 8 nearest neighbors.
        #--------------------------------------------------
        cn = np.zeros( (n1, 8), dtype='int32')
        cn[:,0] = cols-1
        cn[:,1] = cols
        cn[:,2] = cols+1
        cn[:,3] = cols-1
        cn[:,4] = cols+1
        cn[:,5] = cols-1
        cn[:,6] = cols
        cn[:,7] = cols+1
        #--------------------------------------- 
        rn = np.zeros( (n1, 8), dtype='int32')
        rn[:,0] = rows-1
        rn[:,1] = rows-1
        rn[:,2] = rows-1
        rn[:,3] = rows
        rn[:,4] = rows
        rn[:,5] = rows+1
        rn[:,6] = rows+1 
        rn[:,7] = rows+1 
        #------------------
        self.cn = cn
        self.rn = rn
                    
    #   get_neighbor_cols_and_rows()
    #---------------------------------------------------------------
    def get_neighbor_values( self, var, n1 ):
    
        #----------------------------------------
        # Get values of 8 nearest neighbors
        # vals[k,:] = neighbor values of cell k
        #----------------------------------------
        cn = self.cn
        rn = self.rn
        
        vals = np.zeros( (n1, 8), dtype='float32')
        vals[:,0] = var[rn[:,0], cn[:,0]]    # (top left)
        vals[:,1] = var[rn[:,1], cn[:,1]]    # (top center)
        vals[:,2] = var[rn[:,2], cn[:,2]]    # (top right)
        vals[:,3] = var[rn[:,3], cn[:,3]]    # (left center)
        vals[:,4] = var[rn[:,4], cn[:,4]]    # (right center)
        vals[:,5] = var[rn[:,5], cn[:,5]]    # (bottom left)
        vals[:,6] = var[rn[:,6], cn[:,6]]    # (bottom center)
        vals[:,7] = var[rn[:,7], cn[:,7]]    # (bottom right)
        # vals[:,8] = var[rn[:,8], cn[:,8]]  # (center)
        
        return vals
        
    #   get_neighbor_values()    
    #---------------------------------------------------------------
    def spread_conflicts1( self, USE_LOOP=False ):
 
        #-------------------------------------------------   
        # Note:  Can only spread to cells that have S=0.
        #-------------------------------------------------
        w1 = (self.S == 1)
        n1 = w1.sum()
        if (n1 == 0):
            print('No conflicts to spread at time:', self.time_index)
            return

        if (USE_LOOP):
            ID_vals = self.IDs[ w1 ]  #(for the for loop version)
        else:
            ID_vals = np.tile( np.array([self.IDs[w1]]).transpose(), (1,8))
            #------------------
            # This also works
            #------------------
    #         ID_vals = np.zeros((n1,8), dtype='int64')
    #         ID_vals[:,0] = self.IDs[w1]
    #         ID_vals[:,1] = self.IDs[w1]
    #         ID_vals[:,2] = self.IDs[w1]
    #         ID_vals[:,3] = self.IDs[w1]
    #         ID_vals[:,4] = self.IDs[w1]
    #         ID_vals[:,5] = self.IDs[w1] 
    #         ID_vals[:,6] = self.IDs[w1]
    #         ID_vals[:,7] = self.IDs[w1]
      
        #---------------------------------------------
        # Get nearest neighbor values for U, S, & C1
        #---------------------------------------------
        self.get_neighbor_cols_and_rows( w1, n1 )
        #---------------------------------------------
        ## Sn  = self.get_neighbor_values( self.S,  n1) 
        Un  = self.get_neighbor_values( self.U,  n1)
        ## C1n = self.get_neighbor_values( self.C1, n1)

        #------------------------------------------------        
        # Compute probability of spreading to neighbors
        #------------------------------------------------
        # The "None trick" shown here allows us to do
        # the following for all k at once:
        # pn[k,:]  = Un[k,:] * (c2 / Un[k,:].max() )
        # Need c2 = c_spread to be in (0,1].
        # np.amax lets us take the max along an axis.
        #------------------------------------------------
        # NOTE: Un and pn have shape = (n1, 8)
        # NOTE: pn is initialized & defaults to 0.
        #------------------------------------------------    
        Un_max = np.amax( Un, axis=1 )  # a 1D array
        wg = (Un_max > 0)
        pn = np.zeros(Un.shape, dtype='float32')
        pn[ wg,: ] = self.c_spread * Un[wg,:] / (Un_max[wg,None])
       
        #--------------------------------------
        # Alternate method that uses U and C1
        #--------------------------------------
        # Rn = Un * C1n 
        # Rn_max = np.amax( Rn, axis=1 )  # a 1D array
        # wg = (Rn_max > 0)
        # pn = np.zeros(Rn.shape, dtype='float32')
        # pn[ wg,: ] = self.c_spread * Rn[wg,:] / (Rn_max[wg,None])

        #---------------------------------------------
        # Use Bernoulli r.v.s to determine spreading
        #---------------------------------------------
        cn = self.cn
        rn = self.rn
        n_start = self.n_conflict_cells
        
        if (USE_LOOP):        
            for k in range(n1):
                B  = np.random.binomial(1, pn[k,:])  # (8 r.v.s)
                ## B  = B.astype( 'uint8' )  ######
                w2 = np.logical_and( (self.S[rn[k,:], cn[k,:]] == 0),
                                     (B == 1) )
                n2 = w2.sum()
                #-----------------------------------------
                # Spread conflict to some neighbor cells
                #-----------------------------------------
                # w2 is boolean array, but this is okay
                #----------------------------------------------------
                # Set duration to be a geometric r.v.
                ## g1 = np.random.geometric( self.p_geom, size=n2 )
                ## self.durs[ rn[k,w2], cn[k,w2] ] = g1
                #----------------------------------------------------
                self.S[    rn[k,w2], cn[k,w2] ] = 1
                self.IDs[ rn[k,w2], cn[k,w2] ]  = ID_vals[k]
                self.n_conflict_cells += n2
        else:
            #-------------------------------------
            # Spread conflict without a for loop
            # Much faster, and results look similar
            # See notes below re: overcounting.
            #-------------------------------------
            B  = np.random.binomial(1, pn)  # (8 r.v.s)
            ## B  = B.astype( 'uint8' )  ######
            w2 = np.logical_and( (self.S[rn, cn] == 0), (B == 1) )
            n2 = w2.sum()
            #-----------------------------------------
            # Spread conflict to some neighbor cells
            #-----------------------------------------
            # w2 is boolean array, but this is okay
            #-----------------------------------------
            self.S[ rn[w2], cn[w2] ]   = 1
            self.IDs[ rn[w2], cn[w2] ] = ID_vals[w2]
            #--------------------------------------------
            # Without the for loop, several cells can 
            # spread conflict to the same cell and this
            # next line results in over-counting:
            #    self.n_conflict_cells += n2
            #--------------------------------------------
            w_new = (self.S == 1)
            n_new = (w_new.sum() - n1)
            self.n_conflict_cells += n_new
                     
        if (self.REPORT):
            n_spread = (self.n_conflict_cells - n_start)
            print('Number of spread conflicts =', n_spread)
            print()
        
    #   spread_conflicts1()
    #---------------------------------------------------------------
    def finalize( self ):

        #-------------------------    
        # Close the output files
        #-------------------------
        self.out_unit.close()
        self.IDs_unit.close()
        
        if (self.REPORT):
            print()
            run_time = (time.time() - self.start_time)
            run_time = (run_time / 60.0)
            print('run_time  =', run_time, ' [minutes]')
            print('n_steps   =', self.n_steps)
            print('c_emerge  =', self.c_emerge, 'in (0,1)')
            print('c_spread  =', self.c_spread, 'in (0,1)')
            ## print('p_geom   =', self.p_geom)
            print('p_resolve =', self.p_resolve, 'in (0,1)')
            print('time_lag  =', self.time_lag)
            #------------------------------------------
            # S has type 'uint8', so sum may not work
            #------------------------------------------
            ### n_conflict_cells = self.S.sum()
            w = (self.S == 1)
            n_conflict_cells = w.sum()
            #------------------------------------------
            print('n_conflict_cells =', n_conflict_cells)
            print('SELF CHECK...')
            print('n_conflict_cells =', self.n_conflict_cells)
            n_cells = (self.nx - 1) * (self.ny - 1)  # (exclude borders)
            f_conflict_cells = (self.n_conflict_cells / n_cells)
            print('fraction_conflict_cells =', f_conflict_cells)
            print('Finished.')

    #   finalize()
    #---------------------------------------------------------------
    def run_model( self, cfg_file=None ):
    
        self.initialize( cfg_file )
        for k in range( self.n_steps ):
            self.update()
        self.finalize()
        
    #   run_model()
#-------------------------------------------------------------------
#-------------------------------------------------------------------
def get_raster_cellsize( gdal_unit ):

    geotransform = gdal_unit.GetGeoTransform()
    # ulx  = geotransform[0]
    xres = geotransform[1]
    # xrtn = geotransform[2]
    #-----------------------
    # uly  = geotransform[3]
    # yrtn = geotransform[4]  # (not yres !!)
    yres = geotransform[5]  # (not yrtn !!)
    
    return (xres, abs(yres))

#   get_raster_cellsize()
#-------------------------------------------------------------------
def get_raster_bounds( ds, VERBOSE=False):

    #-------------------------------------------------------------
    # Note:  The bounds depend on the map projection and are not
    # necessarily a Geographic bounding box of lons and lats.    
    #-------------------------------------------------------------
    # See:
    # https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html
    # and search on "geotransform".  An example of gdal.SetGeoTransform
    # gives: [xmin, pixel_size, 0, ymax, 0, -pixel_size].
    # Also says args are:
    # [ulx, xDist, rtnX, uly, yDist, rtnY]
    # This is consistent with information below.
    #-------------------------------------------------------------    
    # ulx = upper left x  = xmin
    # uly = upper left y  = ymax
    # lrx = lower right x = xmax
    # lry = lower right y = ymin
    #-----------------------------

    #----------------------------------------------------------  
    # Notice the strange order or parameters here is CORRECT.
    # It is not:  ulx, xres, xskew, uly, yres, yskew
    #----------------------------------------------------------
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    lrx = ulx + (ds.RasterXSize * xres)
    lry = uly + (ds.RasterYSize * yres)

    if (VERBOSE):
        print('ulx, uly   =', ulx, uly)
        print('lrx, lry   =', lrx, lry)
        print('xres, yres = ', xres, yres)
        print('xskew, yskew =', xskew, yskew)
        print('----------------------------------')

    #########################################################
    # Bounding box reported by gdal.info does not match
    # what the GES DISC website is saying.  The result is
    # that gdal.Warp gives all nodata values in output. 
    #########################################################
    return [ulx, lry, lrx, uly]  # [xmin, ymin, xmax, ymax]
 
    #########################################################
    # Bounding box reported by gdal.info does not match
    # what the GES DISC website is saying.  Reversing lats
    # and lons like this doesn't fix the problem.
    #########################################################    
    ## return [lry, ulx, uly, lrx]
    
#   get_raster_bounds()
#------------------------------------------------------------------- 
def bounds_disjoint( bounds1, bounds2, VERBOSE=False):
 
    #-----------------------------------------------------------
    # Note.  Assume both bounds are in same spatial reference
    #        system (SRS), e.g. Geographic lons and lats.
    #------------------------------------------------------------------
    # https://gamedev.stackexchange.com/questions/586/
    # what-is-the-fastest-way-to-work-out-2d-bounding-box-intersection
    #------------------------------------------------------------------    
    b1_xmin = bounds1[0]
    b1_xmax = bounds1[2]
    b2_xmin = bounds2[0]
    b2_xmax = bounds2[2]
#     x_overlap1 = (b1_xmin < b2_xmin) and (b2_xmin < b1_xmax)
#     x_overlap2 = (b2_xmin < b1_xmin) and (b1_xmin < b2_xmax)
#     x_overlap  = (x_overlap1 or x_overlap2)
    
    b1_ymin = bounds1[1]
    b1_ymax = bounds1[3]
    b2_ymin = bounds2[1]
    b2_ymax = bounds2[3]
#     y_overlap1 = (b1_ymin < b2_ymin) and (b2_ymin < b1_ymax) 
#     y_overlap2 = (b2_ymin < b1_ymin) and (b1_ymin < b2_ymax)
#     y_overlap  = (y_overlap1 or y_overlap2)
#     return not(x_overlap and y_overlap)

    disjoint = (b2_xmin > b1_xmax) or (b2_xmax < b1_xmin) or \
               (b2_ymax < b1_ymin) or (b2_ymin > b1_ymax)

    return disjoint
    
#   bounds_disjoint()
#-------------------------------------------------------------------
def read_geotiff(in_file=None, REPORT=True,
                 MAKE_RTG=False, rtg_file=None ):
 
    # Bounds = [ minlon, minlat, maxlon, maxlat ]'

    if (in_file is None):
        in_dir  = '/Users/peckhams/Conflict/Data/GPW-v4/'  
        in_file = 'gpw_v4_population_count_rev11_2020_30_sec.tif'
        in_file = in_dir + in_file

    if (rtg_file is None):
        in_dir   = '/Users/peckhams/Conflict/Data/GPW-v4/'
        rtg_file = 'GPW-v4_global_pop_count.rtg'
        rtg_file = in_dir + rtg_file
          
    #-----------------------------------------    
    # Open the input GeoTIFF file & get info
    #-----------------------------------------
    print('Reading grid from GeoTIFF file...')
    in_unit     = gdal.Open( in_file, gdal.GA_ReadOnly )
    (dx, dy)    = get_raster_cellsize( in_unit )
    in_xres_deg = dx
    in_yres_deg = dy
    in_xres_sec = (in_xres_deg * 3600.0)
    in_yres_sec = (in_yres_deg * 3600.0)
    in_ncols    = in_unit.RasterXSize
    in_nrows    = in_unit.RasterYSize
    in_bounds   = get_raster_bounds( in_unit )   ######

    #------------------------------
    # Get min & max of input grid
    #------------------------------
    in_grid = in_unit.ReadAsArray()
    in_gmin = in_grid.min()
    in_gmax = in_grid.max()
     
    #----------------       
    # Close in_file
    #----------------
    in_unit = None   # Close in_file

    #-----------------------------
    # Option to save as RTG file
    #-----------------------------
    if (MAKE_RTG):
        rtg_path = in_dir + rtg_file   #####
        #-------------------------------
        # Option to create an RTI file
        #-------------------------------
        rti = rti_files.make_info(grid_file=rtg_path,
                  ncols=in_ncols, nrows=in_nrows,
                  xres=in_xres_sec, yres=in_yres_sec,
                  #--------------------------------------
                  data_source='SEDAC',
                  data_type='FLOAT',
                  byte_order=rti_files.get_rti_byte_order(),
                  pixel_geom=0,
                  zres=0.001, z_units='unknown',
                  y_south_edge=in_bounds[1],
                  y_north_edge=in_bounds[3],
                  x_west_edge=in_bounds[0],
                  x_east_edge=in_bounds[2],
                  box_units='DEGREES')
        rti_files.write_info( rtg_path, rti )    
        rtg = rtg_files.rtg_file() 
        OK  = rtg.open_new_file( rtg_path, rti )
        if not(OK):
            print('ERROR during open_new_file().')
            return       
        rtg.write_grid( in_grid, VERBOSE=True )
        rtg.close_file()
         
    #------------------
    # Optional report
    #------------------
    if (REPORT):
        print('GeoTIFF grid info:')
        print('   ' + in_file )
        print('   ncols  =', in_ncols )
        print('   nrows  =', in_nrows )
        print('   xres   =', in_xres_sec, ' [arcsecs]' )
        print('   yres   =', in_yres_sec, ' [arcsecs]' )
        print('   bounds =', in_bounds )
        print('   gmin   =', in_gmin )
        print('   gmax   =', in_gmax ) 
        print('Finished.')
        print()

    #------------------------------
    # Return the grid as an array
    #------------------------------
    return in_grid

#   read_geotiff()
#-------------------------------------------------------------------
def regrid_geotiff(in_file=None, out_file=None, 
                   out_bounds=None,
                   out_xres_sec=None, out_yres_sec=None,
                   RESAMPLE_ALGO='bilinear', REPORT=True):
                   ##### in_nodata=None, out_nodata=None):
                   
    #--------------------------------------------------------
    # Notes:  Read grid from GeoTIFF file, optionally clip
    #         and resample it, then save result as GeoTIFF.
    #--------------------------------------------------------

    #-----------------------------------   
    # Specify the resampling algorithm
    #-----------------------------------
    algo_dict = {
    'nearest'     : gdal.GRA_NearestNeighbour,
    'bilinear'    : gdal.GRA_Bilinear,
    'cubic'       : gdal.GRA_Cubic,
    'cubicspline' : gdal.GRA_CubicSpline,
    'lanczos'     : gdal.GRA_Lanczos,
    'average'     : gdal.GRA_Average,
    'min'         : gdal.GRA_Min,
    'max'         : gdal.GRA_Max,
    'mode'        : gdal.GRA_Mode,
    'med'         : gdal.GRA_Med }
    
    resample_algo = algo_dict[ RESAMPLE_ALGO ]
    
    #---------------------------------------------------------------
    # Note:  DEM_bounds = [dem_xmin, dem_ymin, dem_xmax, dem_ymax]
    #        Give xres, yres in *arcseconds* for Geographic.
    #        They will then be converted to decimal degrees.
    #        gdal.Warp() clips to a bounding box, and can also
    #        resample to a different resolution.
    #        gdal.Translate() is faster for simple clipping.
    #---------------------------------------------------------------
    if (in_file == None):
        #--------------------------------
        # Use Horn of Africa as a test.
        #--------------------------------
        test_dir   = '/Users/peckhams/Conflict/Data/GPW-v4/'
        in_file    = test_dir + 'gpw_v4_population_count_rev11_2020_30_sec.tif'
        out_file   = test_dir + 'Horn_of_Africa_pop_count.tif'
        out_xres_sec = None    # will default to in_xres_sec
        out_yres_sec = None    # will default to in_yres_sec
        # Bounds = [ minlon, minlat, maxlon, maxlat ]
        out_bounds = [ 25.0, -5.0, 55.0, 25.0]
  
    #-----------------------------------------    
    # Open the input GeoTIFF file & get info
    #-----------------------------------------
    in_unit     = gdal.Open( in_file, gdal.GA_ReadOnly )
    (dx, dy)    = get_raster_cellsize( in_unit )
    in_xres_deg = dx
    in_yres_deg = dy
    in_xres_sec = (in_xres_deg * 3600.0)
    in_yres_sec = (in_yres_deg * 3600.0)
    in_ncols    = in_unit.RasterXSize
    in_nrows    = in_unit.RasterYSize
    in_bounds   = get_raster_bounds( in_unit )   ######

    #------------------------------
    # Get min & max of input grid
    #------------------------------
    in_grid = in_unit.ReadAsArray()
    in_gmin = in_grid.min()
    in_gmax = in_grid.max()
    
    #-------------------------------------------------------------
    # If a spatial resolution has not been specified for output,
    # then assume it is the same as the resolution of the input.
    #-------------------------------------------------------------
    if (out_xres_sec is not None):
        out_xres_deg = (out_xres_sec / 3600.0)  # [arcsec] -> [degrees]
    else:
        out_xres_deg = None    # (will default to in_xres_deg)
        ## out_xres_sec = in_xres_sec
        ## out_xres_deg = (out_xres_sec / 3600.0)  # [arcsec] -> [degrees]       
    #----------------------------------------------------------------------
    if (out_yres_sec is not None):
        out_yres_deg = (out_yres_sec / 3600.0)  # [arcsec] -> [degrees]
    else:
        out_yres_deg = None  # (will default to in_yres_deg)
        ## out_yres_sec = in_yres_sec
        ## out_yres_deg = (out_yres_sec / 3600.0)  # [arcsec] -> [degrees]     
                
    #------------------------------------------- 
    # Are out_bounds disjoint from in_bounds ?
    #-------------------------------------------
    # Specifying out_bounds is not required.
    #-------------------------------------------
    if (out_bounds is not None):   
        DISJOINT = bounds_disjoint( in_bounds, out_bounds, VERBOSE=False)
        if (DISJOINT):
            # print('ERROR:  Output bounds do not overlap input bounds.')
            print('ERROR:  Input & output bounding boxes are disjoint.')
            print( '       New grid would contain only nodata.')
            print('  in_file  =', in_file )
            print('  out_file =', out_file )
            in_unit  = None   # Close in_file
            return

    #------------------------------------------------  
    # Resample & clip and write new grid to GeoTIFF
    #------------------------------------------------
    out_unit = gdal.Warp( out_file, in_unit,
        format = 'GTiff',  # (output format string)
        outputBounds=out_bounds,
        xRes=out_xres_deg, yRes=out_yres_deg,
        # srcNodata = in_nodata,      ########  FUTURE
        # dstNodata = out_nodata,     ########  FUTURE
        resampleAlg = resample_algo )

    #-----------------------
    # Get info on new grid
    #-----------------------
    out_ncols  = out_unit.RasterXSize
    out_nrows  = out_unit.RasterYSize
    out_bounds = get_raster_bounds( out_unit )   ######

    #-------------------------------
    # Get min & max of output grid
    #-------------------------------
    out_grid = out_unit.ReadAsArray()
    out_gmin = out_grid.min()
    out_gmax = out_grid.max()
        
    #----------------------------------        
    # Close both in_file and out_file
    #----------------------------------
    in_unit  = None   # Close in_file
    out_unit = None   # Close out_file

    #------------------
    # Optional report
    #------------------
    if (REPORT):
        print('Input grid file:')
        print('   ' + in_file )
        print('   ncols  =', in_ncols )
        print('   nrows  =', in_nrows )
        print('   xres   =', in_xres_sec, ' [arcsecs]' )
        print('   yres   =', in_yres_sec, ' [arcsecs]' )
        print('   bounds =', in_bounds )
        print('   gmin   =', in_gmin )
        print('   gmax   =', in_gmax )
        print()
        print('Output grid file:')
        print('   ' + out_file )
        print('   ncols  =', out_ncols )
        print('   nrows  =', out_nrows )
        print('   xres   =', out_xres_sec, ' [arcsecs]' )
        print('   yres   =', out_yres_sec, ' [arcsecs]' ) 
        print('   bounds =', out_bounds )
        print('   gmin   =', out_gmin )
        print('   gmax   =', out_gmax )      
        print('Finished regridding.')
        print()

    #--------------------------------------------------------  
    # This shows some of the other keywords to gdal.Warp.
    #--------------------------------------------------------      
    # WarpOptions(options=[], format=None, outputBounds=None,
    # outputBoundsSRS=None, xRes=None, yRes=None, targetAlignedPixels=False,
    # width=0, height=0, srcSRS=None, dstSRS=None, srcAlpha=False,
    # dstAlpha=False, warpOptions=None, errorThreshold=None,
    # warpMemoryLimit=None, creationOptions=None, outputType=GDT_Unknown,
    # workingType=GDT_Unknown, resampleAlg=None, srcNodata=None,
    # dstNodata=None, multithread=False, tps=False, rpc=False,
    # geoloc=False, polynomialOrder=None, transformerOptions=None,
    # cutlineDSName=None, cutlineLayer=None, cutlineWhere=None,
    # cutlineSQL=None, cutlineBlend=None, cropToCutline=False,
    # copyMetadata=True, metadataConflictValue=None,
    # setColorInterpretation=False, callback=None, callback_data=None)
    # 
    # Create a WarpOptions() object that can be passed to gdal.Warp()
    # Keyword arguments are : options --- can be be an array of strings,
    # a string or let empty and filled from other keywords.

#    regrid_geotiff()
#-------------------------------------------------------------------
def read_acled_data( excel_file='acled.xlsx', data_dir=None,
                     bounds=None, year_range=[1997,2019],
                     dlat=0.5, dlon=0.5, 
                     rts_file='Horn_of_Africa_deaths_by_year.rts',
                     rtg_file='Horn_of_Africa_deaths_in_year1.rtg'):

    #-----------------------------------------------------
    # Note:  To use this, need:
    #        conda install pandas
    #        conda install openpyxl (for excel files)
    #        DON'T: conda install xlrd  (deprecated)
    #-----------------------------------------------------
    # Note:  Add ability to restrict to bounding box.
    #        ACLED data is for all of Africa, from 1997.
    #-----------------------------------------------------
    # Note:  bounds = [ minlon, minlat, maxlon, maxlat ]
    #        For the Horn of Africa:
    #        bounds = [ 25.0, -5.0, 55.0, 25.0]
    #        if (dlat == 1) and (dlon == 1), then
    #            ncols = 30, nrows = 30
    #-----------------------------------------------------
    if (data_dir is None):
        home_dir = os.path.expanduser('~') + os.sep
        data_dir = home_dir + 'Data/ACLED_Data/'
    excel_filepath = data_dir + excel_file

    df = pd.read_excel( excel_filepath )

    # print(df.columns)  # (print the column header labels)
    # print(df.dtypes)   # (print data type of each column)
    
    lats   = df['LATITUDE'].to_numpy()       # (float64)
    lons   = df['LONGITUDE'].to_numpy()      # (float64)
    deaths = df['FATALITIES'].to_numpy()   # (int64)
    dates  = df['EVENT_DATE']
    years  = df['YEAR'].to_numpy()  # (int64)
  
    #-----------------------------------------------------
    # Note:  bounds = [ minlon, minlat, maxlon, maxlat ]
    #-----------------------------------------------------  
    if (bounds is None):
        minlon = np.floor( lons.min() )
        maxlon = np.ceil( lons.max() )
        minlat = np.floor( lats.min() )
        maxlat = np.ceil( lats.max() )
        IN_BOX = np.ones( lats.size, dtype='bool_' )
    else:
        minlon  = bounds[0]
        minlat  = bounds[1]
        maxlon  = bounds[2]
        maxlat  = bounds[3]
        w1 = np.logical_and( lats >= minlat, lats <= maxlat )
        w2 = np.logical_and( lons >= minlon, lons <= maxlon )
        IN_BOX = np.logical_and( w1, w2 )

    nlons = np.int32( (maxlon - minlon) / dlon )
    nlats = np.int32( (maxlat - minlat) / dlat )
    cols  = np.linspace(minlon, maxlon, nlons)
    rows  = np.linspace(minlat, maxlat, nlats)

    #---------------------------------------------------   
    # Create an RTS file, where each grid in the stack
    # has the number of deaths per cell for one year
    #---------------------------------------------------
    rts_unit = open(rts_file, 'wb')
    rtg_unit = open(rtg_file, 'wb')
    minyear = year_range[0]
    maxyear = year_range[1]
    n_years = (maxyear - minyear) + 1
    year_nums = np.arange(n_years) + minyear
    for year in year_nums:
        print('Working on year:', year)
        w3 = np.logical_and( IN_BOX, years==year )
        year_lons   = lons[w3]
        year_lats   = lats[w3]
        year_deaths = deaths[w3] 
        # year_dates = dates[w3]
        #------------------------------------------------- 
        # Get grid cell row,col for a row in spreadsheet
        #-------------------------------------------------        
        cell_cols   = np.searchsorted( cols, year_lons )
        cell_rows   = np.searchsorted( rows, year_lats )
        death_grid  = np.zeros( (nlats, nlons), dtype='float32')
        n_conflicts = year_lons.size
        for k in range(n_conflicts):
            death_grid[ cell_rows[k], cell_cols[k]] += year_deaths[k]
        # write death_grid to RTS file
        death_grid.tofile( rts_unit )

        #----------------------------------        
        # Make an RTG file for first year
        #----------------------------------
        if (year == minyear):
            death_grid.tofile( rtg_unit)
            rtg_unit.close()
            
    #---------------------       
    # Close the RTS file
    #---------------------
    rts_unit.close()
    print('Finished writing RTS file:')
    print(rts_file)
    return
    #----------------------------------------------------
    # Remaining code is for TOTAL deaths in year_range.
    #----------------------------------------------------
   

    #---------------------------------------------------- 
    # Get grid cell row,col for each row in spreadsheet
    #----------------------------------------------------
    # https://stackoverflow.com/questions/51573164/
    # dividing-geographical-region-into-equal-sized-
    # grid-and-retrieving-indexing-posit
    #----------------------------------------------------
    cell_cols = np.searchsorted( cols, lons )
    cell_rows = np.searchsorted( rows, lats )

    #--------------------------------------------------    
    # Create a geospatial grid of all conflict deaths
    #---------------------------------------------------
    # This for loop is slow, but we need to accumulate
    # the death count from different events.
    # Need a numpy trick to avoid the for loop.
    #-------------------------------------------------------
    # NOTE!  Rows 1 to 74 in the ACLED data for different
    #        months in 1999 all have the same number of
    #        deaths, 1369, which is the max over all rows.
    #        Rows 81 to 198 all have exactly 1000 deaths.
    #        Rows 256-258 are identical except for dates
    #        (3 consecutive days) and each has 400 deaths.
    #        The cause of this overcounting is not clear.
    #-------------------------------------------------------   
    n_conflicts = lons.size
    for k in range(n_conflicts):
        death_grid[ cell_rows[k], cell_cols[k]] += deaths[k]
    return death_grid

    #--------------------------------------------------    
    # Create a geospatial grid of all conflict deaths
    # that occurred in a given year
    #--------------------------------------------------
#     death_grid  = np.zeros( (nlats, nlons), dtype='float64')
#     n_conflicts = lons.size
#     for k in range(n_conflicts):
#         death_grid[ cell_rows[k], cell_cols[k]] += deaths[k]
            
#   read_acled_data()
#-------------------------------------------------------------------

