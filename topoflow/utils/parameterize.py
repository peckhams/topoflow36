#------------------------------------------------------------------------
#  Copyright (c) 2019, Scott D. Peckham
#
#  Sep 2019.  First version to put in utils folder.
#------------------------------------------------------------------------

#  test1()
#  get_grid_from_TCA()

#------------------------------------------------------------------------

from . import rtg_files   # (also in utils directory)
from . import rti_files   # (also in utils directory)
# from topoflow.utils import rtg_files
# from topoflow.utils import rti_files

#------------------------------------------------------------------------
def test1():
    get_grid_from_TCA()

#   test1()    
#------------------------------------------------------------------------
def get_grid_from_TCA( site_prefix=None, cfg_dir=None, 
                       area_file=None, out_file=None,
                       max_value=None, p=None, g0=None):

    if (site_prefix is None):
        site_prefix = 'Treynor'
    if (cfg_dir is None):
        cfg_dir = '/Users/peckhams/Dropbox/TopoFlow_3.6/topoflow/'
        cfg_dir += 'examples/Treynor_Iowa_30m/'
    if (area_file is None):
        area_file = 'Treynor_area.rtg'
    if (out_file is None):
        out_file = 'Treynor_d-bank.rtg'
    if (p is None):
        p = 0.4   # empirical hydraulic geometry for d(Q).
    if (max_value is None):
        max_value = 0.3  # [meters]

    #-------------------------------------------------------
    # Read header_file info for both area_file and out_file
    #-------------------------------------------------------
    header_file = (cfg_dir + site_prefix + '.rti')
    grid_info = rti_files.read_info( header_file, REPORT=True)
    
    #-------------------------------------------------------
    # Read the area grid (i.e. D8 total contributing area)
    #-------------------------------------------------------
    # Also available as d8.A
    #-------------------------------------------------------------
    area_file2 = cfg_dir + area_file
    A = rtg_files.read_grid( area_file2, grid_info, RTG_type='FLOAT' )
    #---------------------------------------------- 
#     rtg = rtg_files.rtg_file()
#     OK  = rtg.open_file( file_name=area_file )   
#     if not(OK):
#         return
#     A = rtg.read_grid( RTG_type='FLOAT' )
    #----------------------------------------------
    A_max = A.max()

    #----------------------------------- 
    # Convert units from km^2 to m^2 ?
    #-----------------------------------  
    # A = A * 1000000.0  [km^2 -> m^2]
    
    #------------------------------------
    # Method 1:  Assume that p is given
    #------------------------------------
    c = max_value / (A_max**p)

    #---------------------------------------------------------
    # Method 2:  Assume that another pair of values is given
    #            as (g0, A0), then compute p
    #---------------------------------------------------------
#     A_min = A.min()
#     A0    = 10.0 * A_min  # (then use d_bankfull = 0.05)
#     p = np.log( max_value / g0 ) / np.log( A_max / A0 )
        
    #-------------------------------------------
    # Compute the new grid as power law of TCA
    #-------------------------------------------
    grid = c * A**p

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
