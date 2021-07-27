import numpy as np

#--------------------------------------------------------------------------
# nx = 29
# ny = 44
# R  = (30.0 / (1000.0 * 3600.0)) # [m/s]
# da = 900  # [m2]
# dur = (75 * 60) # [sec]  (75 minutes or 4500 sec)
# A = (nx * ny * da)
# vol_tot = (A * R * dur)
# print( vol_tot )

#--------------------------------------------------------------------------
# dtype = 'float32'
# nx = 29
# ny = 44
# n_values = nx * ny
# grid_file = 'June_20_67_rain_uniform_30.rtg'
# grid_unit = open( grid_file, 'rb' )
# grid = np.fromfile( grid_unit, count=n_values, dtype=dtype )
# grid.byteswap( True )
# grid = grid.reshape( ny, nx )
# print('min(grid) =', grid.min())
# print('max(grid) =', grid.max())
# # Try to read again from grid_unit
# grid = np.fromfile( grid_unit, count=n_values, dtype=dtype )
# print( grid.dtype )
# print( grid.size )
# print( grid )
# # print('min(grid) =', grid.min())
# # print('max(grid) =', grid.max())
# grid_unit.close()
#--------------------------------------------------------------------------
def create_rainfall_grid():

    # Create a grid of uniform rainfall
    nx = 29
    ny = 44
    grid = np.zeros((ny,nx), dtype='float32') + 60.0  # [mmph]
    new_rtg_file = 'June_20_67_rain_uniform_60.rtg'
    rtg_unit = open(new_rtg_file, 'wb')
    grid.byteswap( True )
    grid.tofile( rtg_unit )
    rtg_unit.close()
    
#   create_rainfall_grid()
#--------------------------------------------------------------------------
def create_rainfall_grid_stack():

    # Create a grid stack of uniform rainfall
    nx = 29
    ny = 44
    grid = np.zeros((ny,nx), dtype='float32') + 30.0  # [mmph]
    new_rts_file = 'June_20_67_rain_uniform_30_75min.rts'
    rts_unit = open(new_rts_file, 'wb')
    grid.byteswap( True )
    for k in range(75):   # 75 grids for 75 minutes
        grid.tofile( rts_unit )
    rts_unit.close()
    
#   create_rainfall_grid_stack()
#--------------------------------------------------------------------------
def read_rainfall_grid():

    # Create a grid of uniform rainfall
    dtype = 'float32'
    nx = 29
    ny = 44
    n_values = nx * ny
    grid_file = 'June_20_67_rain_uniform_60.rtg'
    grid_unit = open( grid_file, 'rb' )
    grid = np.fromfile( grid_unit, count=n_values, dtype=dtype )
    grid.byteswap( True )
    grid = grid.reshape( ny, nx )
    grid_unit.close()
    print('min(grid) =', grid.min())
    print('max(grid) =', grid.max())
    return grid
    
#   read_rainfall_grid()
#--------------------------------------------------------------------------



