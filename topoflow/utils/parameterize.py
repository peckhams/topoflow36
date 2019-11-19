#------------------------------------------------------------------------
#  Copyright (c) 2019, Scott D. Peckham
#
#  Sep 2019.  First version to put in utils folder.
#------------------------------------------------------------------------

#  test1()
#  test2()
#  get_grid_from_TCA()

#------------------------------------------------------------------------

import numpy as np
from . import rtg_files   # (also in utils directory)
from . import rti_files   # (also in utils directory)
# from topoflow.utils import rtg_files
# from topoflow.utils import rti_files

#------------------------------------------------------------------------
def test1():
    get_grid_from_TCA()

#   test1()
#------------------------------------------------------------------------
def test2():

    site_prefix = 'Baro_Gam_1min'
    cfg_dir = '/Users/peckhams/Desktop/TopoFlow_2019/__Regions/'
    cfg_dir += 'Ethiopia/DEMs/MERIT/Baro_Gam_1min/'

    #--------------------------------------------    
    # Compute channel bottom width grid [meters]
    #---------------------------------------------------------------
    # Baro River at Gambela (lon, lat) = (34.5833333, 8.25)
    # Birbir River mouth (lon, lat)    = (34.9609, 8.2411)
    # Gebba River mouth (lon, lat)     = (34.9609, 8.2411; almost)
    #---------------------------------------------------------------
    # Baro River at Gambela (col, row) = (21,75)
    # Birbir River mouth (col, row)    = (44,74)
    # Gebba River mouth (col, row)     = (45,76)
    #----------------------------------------------
    # Widths estimated with Google Maps.
    #----------------------------------------------       
    get_grid_from_TCA(site_prefix=site_prefix, cfg_dir=cfg_dir,
        area_file='Baro_Gam_1min_area.rtg',
        out_file='Baro_Gam_1min_chan-w.rtg',
        p=0.5, g1=140.0 )
        #--------------------------------------
        # This method gives tiny c and p > 5.
        #--------------------------------------        
        # A1=23567.7, g1=140.0,    # Baro River at Gambella
        # A2=16391.1, g2=30.0 )    # mouth of Birbir River.
        ## A2=4692.35, g2=18.0 )   # mouth of Gebba River

    #------------------------------------    
    # Compute Manning's n grid [...]
    #------------------------------------
    get_grid_from_TCA(site_prefix=site_prefix, cfg_dir=cfg_dir,
        area_file='Baro_Gam_1min_area.rtg',
        out_file='Baro_Gam_1min_chan-n.rtg',
        g1=0.03, g2=0.2 )
        # A1=23567.7, g1=0.03,      # Baro River at Gambella
        # A2=16391.1, g2=0.035 )    # mouth of Birbir River.
        ## A2=4692.35, g2=0.04 )   # mouth of Gebba River
 
    #----------------------------------------   
    # Compute absolute sinuosity grid [m/m]
    #-----------------------------------------------------------
    # NOTE:  The minimum possible value is 1.0, by definition.
    #-----------------------------------------------------------
    get_grid_from_TCA(site_prefix=site_prefix, cfg_dir=cfg_dir,
        area_file='Baro_Gam_1min_area.rtg',
        out_file='Baro_Gam_1min_sinu.rtg',
        g1=1.3, g2=1.0 )
        #-----------------------------------------------------
        #A1=23567.7, g1=1.3,      # Baro River at Gambella
        # A2=16391.1, g2=1.2 )    # mouth of Birbir River.
        ## A2=4692.35, g2=1.1 )   # mouth of Gebba River

    #----------------------------------------   
    # Estimate bankfull depth grid [meters]
    #----------------------------------------
    get_grid_from_TCA(site_prefix=site_prefix, cfg_dir=cfg_dir,
        area_file='Baro_Gam_1min_area.rtg',
        out_file='Baro_Gam_1min_d-bank.rtg',
        g1=8.0, p=0.4 )  # Baro River at A_max.
        ## A1=23567.7, g1=8,      # Baro River at Gambella
        ## A2=16391.1, g2= )    # mouth of Birbir River.
                           
#   test2()   
#------------------------------------------------------------------------
def get_grid_from_TCA( site_prefix=None, cfg_dir=None, 
                       area_file=None, out_file=None,
                       A1=None, g1=None,
                       A2=None, g2=None,
                       p=None, REPORT=True):

    if (site_prefix is None):
        site_prefix = 'Treynor'
    if (cfg_dir is None):
        cfg_dir = '/Users/peckhams/Dropbox/TopoFlow_3.6/topoflow/'
        cfg_dir += 'examples/Treynor_Iowa_30m/'
    if (area_file is None):
        area_file = site_prefix + '_area.rtg'
    if (out_file is None):
        out_file = site_prefix + 'chan-w.rtg'

    #-------------------------------------------------------
    # Read header_file info for both area_file and out_file
    #-------------------------------------------------------
    header_file = (cfg_dir + site_prefix + '.rti')
    grid_info = rti_files.read_info( header_file, REPORT=False)
    
    #-------------------------------------------------------
    # Read the area grid (i.e. D8 total contributing area)
    #-------------------------------------------------------
    # Also available as d8.A
    #-------------------------------------------------------------
    area_file2 = cfg_dir + area_file
    A = rtg_files.read_grid( area_file2, grid_info, RTG_type='FLOAT' )

    #----------------------------------- 
    # Convert units from km^2 to m^2 ?
    #-----------------------------------  
    # A = A * 1000000.0  [km^2 -> m^2]  # (not needed)

    #------------------------------
    # Set c and p based on inputs
    #------------------------------
    if (g1 is not None) and (g2 is not None) and \
       (A1 is not None) and (A2 is not None):
        #-------------------------------------------------
        # Option 1: Use (A1,g1) and (A2,g2) to set (c,p)
        #-------------------------------------------------
        p = np.log( g1 / g2 ) / np.log( A1 / A2 )
        c = g1 / (A1**p)
    elif (p is not None) and (g1 is not None):
        #-------------------------------------------
        # Option 2:  p and g1=g(A_max) are given.
        #-------------------------------------------
        A1 = A.max()
        c  = g1 / (A1**p)
    elif (g1 is not None) and (g2 is not None):
        #----------------------------------------------------
        # Option 3:  g1=g(A_max) and g2=g(A_min) are given.
        # But can have A=0 on the edges, so use min pos.
        #----------------------------------------------------
        A1 = A.max()
        w  = (A > 0)  # (array of True or False)
        A2 = A[ w ].min()
        p  = np.log( g1 / g2 ) / np.log( A1 / A2 )
        c  = g1 / (A1**p)    
    else:
        print('Sorry, not enough input parameters to compute c & p.')
        return
     
    if (REPORT):
        print('Power-law parameters are:')
        print('c = ' + str(c) )
        print('p = ' + str(p) )

    #-------------------------------------------- 
    # Find cells that have TCA = 0 (e.g. edges)
    #--------------------------------------------
    w1 = (A <= 0)  # (array of True or False)
            
    #-------------------------------------------
    # Compute the new grid as power law of TCA
    #-------------------------------------------
    # The TCA area grid has zeros on the edges
    #-------------------------------------------
    if (p < 0):
        #--------------------------------------
        # Need to avoid NaNs in computed grid
        #--------------------------------------
        Ac = A.copy()
        Ac[ w1 ] = 1.0    # (temporary, to avoid warnings)
        grid = c * Ac**p
    else:
        grid = c * A**p

    #-------------------------------------------------------------
    # Zero is okay for some variables, but is not okay for width
    # or Manning's n because we'll get a divide-by-zero error.
    #-------------------------------------------------------------
    ## grid[ w1 ] = np.nan    # (TopoFlow doesn't like these.)
    ## grid[ w1 ] = -1.0
    ## grid[ w1 ] = 0.0
    grid[ w1 ] = 1.0
    
    if (REPORT):
        nbad = w1.sum()
        if (nbad > 0):
            print('Values set to 1.0 where A <= 0.')
            print('  This occurred at', nbad, 'grid cells.')
        print('grid min =', grid.min() )
        print('grid_max =', grid.max() )
        
    #----------------------------------
    # Write computed grid to out_file
    #----------------------------------
    out_file2 = (cfg_dir + out_file)
    rtg_files.write_grid( grid, out_file2, grid_info, RTG_type='FLOAT')
    print( 'Finished writing file: ')
    print( out_file2 )
    print()

#   get_grid_from_TCA()
#------------------------------------------------------------------------
