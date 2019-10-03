#
# Copyright (c) 2001-2014, Scott D. Peckham
#
# Sept 2014.  Moved some functions into outlets.py to avoid cyclic
#             dependencies between BMI_base.py and basins.py.
# Nov. 2016.  This file may be obsolete now, replaced by outlets.py.  ########
#
# January, August 2009
# May 2010 (changes to unit_test(), initialize(), etc.)
#
################################################################
#
# NB!  The update_volume_in() method ONLY tracks precip now.
#      "channels_base.py" now has update_volume_out() also.
#
################################################################

#-----------------------------------------------------------------------
#
#  unit_test()
#
#  class basins_component  (inherits from BMI_base.py)
#
#      initialize()
#      update()          # (non-OpenMI arguments)
#      finalize()
#      read_config_file()
#      -----------------------
#      update_volume_in()    # (commented out)
#      update_volume_out()   # (commented out)
#
#-----------------------------------------------------------------------

## import numpy as np  # (no longer needed)

import os
import os.path

from . import BMI_base
from . import outlets
from . import tf_utils

## from topoflow.utils import BMI_base
## from topoflow.utils import outlets
## from topoflow.utils import tf_utils

#-----------------------------------------------------------------------
def unit_test():

    b = basins_component()
    b.CCA   = False
    b.DEBUG = False
    
    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    cfg_directory = tf_utils.TF_Test_Directory()
    ## os.chdir( cfg_directory )
    b.cfg_directory = cfg_directory
    b.site_prefix = 'Treynor'  ###########

    ## cfg_prefix = tf_utils.TF_Test_case_prefix()     
    cfg_file = None  #####
    
    b.initialize( cfg_file=cfg_file, mode='driver' )

    print('outlet_ID    =', b.outlet_ID)
    print('basin_area   =', b.basin_area)
    print('basin_relief =', b.basin_relief)
    print(' ')
    print('n_outlets    =', b.n_outlets)
    print('outlet_cols  =', b.outlet_cols)
    print('outlet_rows  =', b.outlet_rows)
    print('reliefs      =', b.basin_reliefs)
    print('areas        =', b.basin_areas)
    print(' ')
    print('nx           =', b.nx)
    print('ny           =', b.ny)
    print(' ')
    print("get_status()            = ", b.get_status())
    print("is_scalar('n_outlets')  = ", b.is_scalar('n_outlets'))
    print("is_grid('n_outlets')    = ", b.is_grid('n_outlets'))
    # Next one has double size, since its really a tuple.
    print('b.outlet_IDs =', b.outlet_IDs)
    print('b.basin_area =', b.basin_area)
    print('Finished with unit_test().')
    print(' ')
    
#   unit_test()   
#-----------------------------------------------------------------------
class basins_component( BMI_base.BMI_component ):
 
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        if not(SILENT):
            print('Basins component: Initializing...')

        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file
        
        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        ## self.set_constants()
        self.initialize_config_vars() 
        self.read_grid_info()
 
        #---------------------------------------------
        # Read outlet IDs (IDs of monitored pixels)
        # and their attributes like area and relief.
        # Then read IDs of all cells in the first
        # (or main) basin, i.e. above first outlet.
        #---------------------------------------------
        outlets.read_outlet_file( self )   # (uses nx and ny)
        
        # outlets.read_main_basin_IDs( self )

        #-------------------------------------------
        # Prepare to track total water in and out
        # of the main basin (using basin RTM file)
        #-------------------------------------------
        self.get_pvolume = False   ####
        TRACK_VOLUME = (self.get_pvolume and (self.basin_RTM_file != ''))
        self.TRACK_VOLUME = TRACK_VOLUME
        if (TRACK_VOLUME):
            #-------------------------------------------
            # Prepare to track total water in and out
            # of the main basin (using basin RTM file)
            #----------------------------------------------
            # This requires knowing the IDs of the pixels
            # that lie within the basin (basin_IDs).
            #----------------------------------------------
            outlets.read_main_basin_IDs( self )
            self.volume_in  = self.initialize_scalar( 0, dtype='float64')
            self.volume_out = self.initialize_scalar( 0, dtype='float64')

        self.status = 'initialized'
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, Q, time, dt, da, pv):

        self.status = 'updating'  # (OpenMI)
        
        if (self.TRACK_VOLUME):
            self.update_volume_out(Q, dt)
            self.update_volume_in(time, dt, da, pv)

        #------------------------
        # Update internal clock
        #------------------------
        # self.update_time()

        self.status = 'updated'  # (OpenMI)
        
    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalized'  # (OpenMI)

    #   finalize()
    #-------------------------------------------------------------------
    def read_config_file(self):

        #---------------------------------------------------
        # 6/28/10.  Need this, since there is no CFG file.
        #---------------------------------------------------
        pass
    
    #   read_config_file()
    #-------------------------------------------------------------------

    #-------------------------------------------------------------------
    # NOTE:  May be better to move this into "precip_base.py".
    #-------------------------------------------------------------------
##    def update_volume_in(self, time, dt, da, pv):
##       
##        #----------------------------------------------------------
##        # Notes:  This procedure integrates precip. over the main
##        #         basin using the model vs. sampling timestep.
##
##        #         Recall that da is a grid [km^2].
##        #----------------------------------------------------------
##        if (pv.method == 0): return
##        
##        if (pv.method == 1):
##            #------------------------------------------------
##            # In this case, pv.rates and pv.durations are
##            # 1D vectors (vs. scalar or grid), which does
##            # not conform to the general approach now used
##            # throughout TopoFlow and by PRECIP_METHOD 2.
##            #------------------------------------------------
##            wd  = where(time < pv.duration_sums)
##            nwd = size(wd[0])
##            if (nwd != 0):
##                # rate = pv.rates[wd[0]]     ########
##                rate = pv.rates[wd[0][0]]     ######################
##                dvol = dt * rate * pv.basin_area * float64(1000000)
##                self.volume_in += dvol
##        else:    
##            #----------------------------------------------------
##            # If pv.durations is a scalar, then duration_sums
##            # is equal to the same scalar, but this still works
##            # as written (3/20/07)
##            #----------------------------------------------------
##            n_rates = size(pv.rates)
##            if (n_rates == 1):    
##                P_rates = pv.rates
##            else:    
##                P_rates = pv.rates[self.basin_IDs]
##            #-------------------------------------------------------
##            n_durs = size(pv.duration_sums)
##            if (time <= pv.duration_sums[n_durs - 1]):    
##                if (size(da) == 1):    
##                    nb   = size(self.basin_IDs[0]) ### BUG FIX.
##                    dvol = dt * sum(double(P_rates * da * nb))
##                else:    
##                    dvol = dt * sum(double(P_rates * da[self.basin_IDs]))
##                self.volume_in += dvol
##        
##    #   update_volume_in()
##    #-------------------------------------------------------------------
##    def update_volume_out(self, Q, dt):
##        
##        self.volume_out += (Q[self.outlet_ID] * dt)
##          
##    #   update_volume_out()               
    #-------------------------------------------------------------------
        


