#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created by S.D. Peckham, August 2020.
A script to create movies and images from TopoFlow output files.
This script uses topoflow36/topoflow/utils/visualize.py.
It only creates a few of the many, many possible visualizations.
Many other color tables and contrast-enhancement stretches are
   supported by changing the default keywords.
This script does not yet check for existence of required files.
Directories are created for the images needed to make each movie,
   but they are not deleted afterwards.
See the Jupyter notebook for TopoFlow visualization at:
https://github.com/peckhams/topoflow36/blob/master/TopoFlow_Visualization_v2.ipynb

Required arguments:
site_prefix:  e.g. "Baro-Gam_60sec"
case_prefix:  e.g. "Test1"
output_dir:   (path to the TopoFlow output files)
"""

#---------------------------
# Import required packages
#---------------------------
from topoflow.utils import visualize as tfvis
import os, os.path, sys

#---------------------------
# Check required arguments
#---------------------------
if (len(sys.argv) < 3):
    print('Sorry, this script requires 3 arguments:')
    print('output_dir, site_prefix, and case_prefix.')
    return
output_dir  = sys.argv[1]   # (name of script is sys.argv[0])
site_prefix = sys.argv[2]   # e.g. 'Baro-Gam_60sec'
case_prefix = sys.argv[3]   # e.g. 'Test1'
#---------------------------------------------
if (output_dir[-1] != os.sep):
    output_dir += os.sep
   
#------------------------------- 
# Get names of all directories
#-------------------------------
basin_dir  = output_dir
topo_dir   = basin_dir + '__topo/'
met_dir    = basin_dir + '__met/'
soil_dir   = basin_dir + '__soil/'
print('Output directory =', output_dir)

#-----------------------------------
# Create movie for River Discharge
#-----------------------------------
Q_MOVIE = True
if (Q_MOVIE):
    os.chdir(output_dir)
    rts_png_dir = output_dir + 'rts_png_files_Q/'
    if not(os.path.exists( rts_png_dir )): os.mkdir( rts_png_dir )
    rts_file  = case_prefix + '_2D-Q.rts'
    mp4_file  = case_prefix + '_2D-Q.mp4'
    long_name = 'River Discharge'
    stretch   = 'power3'
    #--------------------------------------------------
    # If extent=None, it is determined automatically.
    #--------------------------------------------------
    tfvis.save_rts_as_images( rts_file, rts_png_dir, extent=None,
                              long_name=long_name,
                              stretch=stretch, a=1, b=2, p=0.5,
                              BLACK_ZERO=False,
                              cmap='rainbow', REPORT=True,
                              xsize=6, ysize=6, dpi=192)
    os.chdir(output_dir)
    fps = 10  # frames per second
    tfvis.create_movie_from_images( mp4_file, rts_png_dir, fps=fps, REPORT=True)

#-------------------------------------
# Create movie for River Flood Depth
#-------------------------------------
DF_MOVIE = True
if (DF_MOVIE):
    os.chdir(output_dir)
    rts_png_dir = output_dir + 'rts_png_files_d_flood/'
    if not(os.path.exists( rts_png_dir )): os.mkdir( rts_png_dir )
    rts_file  = case_prefix + '_2D-d-flood.rts'
    mp4_file  = case_prefix + '_2D-d-flood.mp4'
    long_name = 'River Flood Depth'
    stretch   = 'power3'
    tfvis.save_rts_as_images( rts_file, rts_png_dir, extent=None,
                              long_name=long_name,
                              stretch=stretch, a=1, b=2, p=0.5,
                              BLACK_ZERO=False,
                              cmap='rainbow', REPORT=True,
                              xsize=6, ysize=6, dpi=192)
    os.chdir(output_dir)
    fps = 10  # frames per second
    tfvis.create_movie_from_images( mp4_file, rts_png_dir, fps=fps, REPORT=True)

#-------------------------------------------
# Create movie for CHIRPS rainfall (input)
#-------------------------------------------
CHIRPS_RAIN_MOVIE = False
if (CHIRPS_RAIN_MOVIE):
    os.chdir(met_dir)
    rts_png_dir = output_dir + 'rts_chirps_png_files/'
    if not(os.path.exists( rts_png_dir )): os.mkdir( rts_png_dir )
    rts_file  = 'CHIRPS_Rain.rts'
    mp4_file  = 'CHIRPS_Rain.mp4'
    mp4_file  = output_dir + mp4_file
    long_name = 'Rainfall Rate [mmph]'
    stretch   = 'linear'
    tfvis.save_rts_as_images( rts_file, rts_png_dir, extent=None,
                              long_name=long_name,
                              stretch=stretch, a=1, b=2, p=0.5,
                              BLACK_ZERO=False,
                              cmap='rainbow', REPORT=True,
                              xsize=6, ysize=6, dpi=192)
    os.chdir(output_dir)
    fps = 10  # frames per second
    tfvis.create_movie_from_images( mp4_file, rts_png_dir, fps=fps, REPORT=True)

#------------------------------------------
# Create movie for GLDAS rainfall (input)
#------------------------------------------
GLDAS_RAIN_MOVIE = False
if (GLDAS_RAIN_MOVIE):
    os.chdir(met_dir)
    rts_png_dir = output_dir + 'rts_gldas_png_files/'
    if not(os.path.exists( rts_png_dir )): os.mkdir( rts_png_dir )
    rts_file  = 'GLDAS_Rain.rts'
    mp4_file  = 'GLDAS_Rain.mp4'
    mp4_file  = output_dir + mp4_file
    long_name = 'Rainfall Rate [mmph]'
    stretch   = 'linear'
    tfvis.save_rts_as_images( rts_file, rts_png_dir, extent=None,
                              long_name=long_name,
                              stretch=stretch, a=1, b=2, p=0.5,
                              BLACK_ZERO=False,
                              cmap='rainbow', REPORT=True,
                              xsize=6, ysize=6, dpi=192)
    os.chdir(output_dir)
    fps = 10  # frames per second
    tfvis.create_movie_from_images( mp4_file, rts_png_dir, fps=fps, REPORT=True)

#----------------------------------------
# Create movie for GPM rainfall (input)
#----------------------------------------
GPM_RAIN_MOVIE = False
if (GPM_RAIN_MOVIE):
    os.chdir(met_dir)
    rts_png_dir = output_dir + 'rts_gpm_png_files/'
    if not(os.path.exists( rts_png_dir )): os.mkdir( rts_png_dir )
    rts_file  = 'GPM_Rain.rts'
    mp4_file  = 'GPM_Rain.mp4'
    mp4_file  = output_dir + mp4_file
    long_name = 'Rainfall Rate [mmph]'
    stretch   = 'linear'
    tfvis.save_rts_as_images( rts_file, rts_png_dir, extent=None,
                              long_name=long_name,
                              stretch=stretch, a=1, b=2, p=0.5,
                              BLACK_ZERO=False,
                              cmap='rainbow', REPORT=True,
                              xsize=6, ysize=6, dpi=192)
    os.chdir(output_dir)
    fps = 10  # frames per second
    tfvis.create_movie_from_images( mp4_file, rts_png_dir, fps=fps, REPORT=True)

#--------------------------------------------------
# Create image for time series of River Discharge
# at one of the monitored grid cells (var_index)
#--------------------------------------------------
# if (im_file == None), the plot is displayed.
# Change im_file extension for other formats.
#--------------------------------------------------
Q_SERIES_IMAGE = True
if (Q_SERIES_IMAGE):
    os.chdir(output_dir)
    var_index = 0  # (1st monitored grid cell)
    nc_file = case_prefix + '_0D-Q.nc'
    im_file = case_prefix + '_0D-Q.png'
    tfvis.plot_time_series(nc_file, output_dir=output_dir, var_index=var_index,
                           marker=',', REPORT=True, xsize=10, ysize=4,
                           im_file=im_file)

#--------------------------------------------------
# Create image for time series of River Velocity
# at one of the monitored grid cells (var_index)
#--------------------------------------------------
# if (im_file == None), the plot is displayed.
# Change im_file extension for other formats.
#--------------------------------------------------
U_SERIES_IMAGE = True
if (U_SERIES_IMAGE):
    os.chdir(output_dir)
    var_index = 0  # (1st monitored grid cell)
    nc_file = case_prefix + '_0D-u.nc'
    im_file = case_prefix + '_0D-u.png'
    tfvis.plot_time_series(nc_file, output_dir=output_dir, var_index=var_index,
                           marker=',', REPORT=True, xsize=10, ysize=4,
                           im_file=im_file)

#--------------------------------------------------
# Create image for time series of River Depth
# at one of the monitored grid cells (var_index)
#--------------------------------------------------
# if (im_file == None), the plot is displayed.
# Change im_file extension for other formats.
#--------------------------------------------------
D_SERIES_IMAGE = True
if (D_SERIES_IMAGE):
    os.chdir(output_dir)
    var_index = 0  # (1st monitored grid cell)
    nc_file = case_prefix + '_0D-d.nc'
    im_file = case_prefix + '_0D-d.png'
    tfvis.plot_time_series(nc_file, output_dir=output_dir, var_index=var_index,
                           marker=',', REPORT=True, xsize=10, ysize=4,
                           im_file=im_file)
                                                                                                     
#----------------------------------------------------
# Create image for time series of River Flood Depth
# at one of the monitored grid cells (var_index)
#----------------------------------------------------
# if (im_file == None), the plot is displayed.
# Change im_file extension for other formats.
#--------------------------------------------------
DF_SERIES_IMAGE = True
if (DF_SERIES_IMAGE):
    os.chdir(output_dir)
    var_index = 0  # (1st monitored grid cell)
    nc_file = case_prefix + '_0D-d-flood.nc'
    im_file = case_prefix + '-0D-d-flood.png'
    tfvis.plot_time_series(nc_file, output_dir=output_dir, var_index=var_index,
                           marker=',', REPORT=True, xsize=10, ysize=4,
                           im_file=im_file)

#-----------------------------------------------------
# Create image for the DEM (Digital Elevation Model)
#-----------------------------------------------------
DEM_IMAGE = True
if (DEM_IMAGE):
    os.chdir(output_dir)
    long_name = 'Land Surface Elevation'
    rtg_file  = topo_dir + site_prefix + '_DEM.rtg'
    im_file   = topo_dir + site_prefix + '_DEM.png'
    stretch   = 'hist_equal'
    tfvis.read_and_show_rtg( rtg_file, long_name, cmap='jet',
                             ### BLACK_ZERO=False,
                             im_file=im_file,
                             stretch=stretch, VERBOSE=True,
                             xsize=6, ysize=6, dpi=None)
                                               
#-------------------------------------
# Create image for the D8 Slope grid
#-------------------------------------
SLOPE_IMAGE = True
if (SLOPE_IMAGE):
    os.chdir(output_dir)
    long_name = 'Land Surface Slope'
    rtg_file  = topo_dir + site_prefix + '_slope.rtg'
    im_file   = topo_dir + site_prefix + '_slope.png'
    stretch   = 'hist_equal'
    tfvis.read_and_show_rtg( rtg_file, long_name, cmap='jet',
                             ### BLACK_ZERO=False,
                             im_file=im_file,
                             stretch=stretch, VERBOSE=True,
                             xsize=6, ysize=6, dpi=None)

#------------------------------------
# Create image for the D8 Area grid
#------------------------------------
AREA_IMAGE = True
if (AREA_IMAGE):
    os.chdir(output_dir)
    long_name = 'Total Contributing Area'
    rtg_file  = topo_dir + site_prefix + '_d8-area.rtg'
    im_file   = topo_dir + site_prefix + '_d8-area.png'
    stretch   = 'power3'
    # stretch = 'log'
    # stretch = 'hist_equal'
    tfvis.read_and_show_rtg( rtg_file, long_name, cmap='jet',
                             ### BLACK_ZERO=False,
                             im_file=im_file,
                             stretch=stretch, VERBOSE=True,
                             xsize=6, ysize=6, dpi=None)

#------------------------------------
# Create image for the soil Ks grid
#------------------------------------
KS_IMAGE = True
if (KS_IMAGE):
    os.chdir(output_dir)
    long_name = 'Soil Hydraulic Conductivity (Layer 1)'
    rtg_file  = soil_dir + site_prefix + '_sl1_2D-Ks.bin'
    im_file   = soil_dir + site_prefix + '_sl1_2D-Ks.png'
    stretch   = 'hist_equal'
    tfvis.read_and_show_rtg( rtg_file, long_name, cmap='jet',
                             ### BLACK_ZERO=False,
                             im_file=im_file,
                             stretch=stretch, VERBOSE=True,
                             xsize=6, ysize=6, dpi=None)
                                           
#------------------------------------
# Create image for the soil qs grid
#------------------------------------
QS_IMAGE = True
if (QS_IMAGE):
    os.chdir(output_dir)
    long_name = 'Soil Saturated Water Content (Layer 1)'
    rtg_file  = soil_dir + site_prefix + '_sl1_2D-qs.bin'
    im_file   = soil_dir + site_prefix + '_sl1_2D-qs.png'
    stretch   = 'hist_equal'
    tfvis.read_and_show_rtg( rtg_file, long_name, cmap='jet',
                             ### BLACK_ZERO=False,
                             im_file=im_file,
                             stretch=stretch, VERBOSE=True,
                             xsize=6, ysize=6, dpi=None)                                                                  

#----------------------------------------------------------------
# Create movie for a soil z-profile time series:  water content
# at one of the monitored grid cells (var_index)
#----------------------------------------------------------------
Q_PROFILE_MOVIE = False
if (Q_PROFILE_MOVIE):
    os.chdir(output_dir)
    nc1D_png_dir = output_dir + 'q1D_png_files/'
    if not(os.path.exists( nc1D_png_dir )): os.mkdir( nc1D_png_dir)

    nc_file = case_prefix + '_1D-q.nc'
    ymin = 0.0
    ymax = 0.5
    # Set ymin and ymax to make them same for all plots.
    tfvis.save_profile_series_as_images(nc_file, png_dir=nc1D_png_dir,
               ymin=ymin, ymax=ymax,    # (same for all plots)
               # ymin=None, ymax=None,  # (auto for each plot)
               marker=',', REPORT=True,
               xsize=8, ysize=5, dpi=192)
    os.chdir(output_dir)                   
    fps = 20  # frames per second
    mp4_file = case_prefix + '_Profile_Movie_q.mp4'
    tfvis.create_movie_from_images( mp4_file, nc1D_png_dir, fps=fps, REPORT=True)

#----------------------------------------------------------------
# Create movie for a soil z-profile time series:  pressure head
# at one of the monitored grid cells (var_index)
#----------------------------------------------------------------
P_PROFILE_MOVIE = False
if (P_PROFILE_MOVIE):
    os.chdir(output_dir)
    nc1D_png_dir = output_dir + 'p1D_png_files/'
    if not(os.path.exists( nc1D_png_dir )): os.mkdir( nc1D_png_dir)

    nc_file = case_prefix + '_1D-p.nc'
    ymin = None
    ymax = None

    # Set ymin and ymax to make them same for all plots.
    tfvis.save_profile_series_as_images(nc_file, png_dir=nc1D_png_dir,
               ymin=ymin, ymax=ymax,    # (same for all plots)
               # ymin=None, ymax=None,  # (auto for each plot)
               marker=',', REPORT=True,
               xsize=8, ysize=5, dpi=192)
    os.chdir(output_dir)                   
    fps = 20  # frames per second
    mp4_file = case_prefix + '_Profile_Movie_p.mp4'
    tfvis.create_movie_from_images( mp4_file, nc1D_png_dir, fps=fps, REPORT=True)
      
#----------------------------------------------------------------
# Create movie for a soil z-profile time series:  hydr. conductivity
# at one of the monitored grid cells (var_index)
#----------------------------------------------------------------
K_PROFILE_MOVIE = False
if (K_PROFILE_MOVIE):
    os.chdir(output_dir)
    nc1D_png_dir = output_dir + 'K1D_png_files/'
    if not(os.path.exists( nc1D_png_dir )): os.mkdir( nc1D_png_dir)

    nc_file = case_prefix + '_1D-K.nc'
    ymin = None
    ymax = None
    # Set ymin and ymax to make them same for all plots.
    tfvis.save_profile_series_as_images(nc_file, png_dir=nc1D_png_dir,
               ymin=ymin, ymax=ymax,    # (same for all plots)
               # ymin=None, ymax=None,  # (auto for each plot)
               marker=',', REPORT=True,
               xsize=8, ysize=5, dpi=192)
    os.chdir(output_dir)                   
    fps = 20  # frames per second
    mp4_file = case_prefix + '_Profile_Movie_K.mp4'
    tfvis.create_movie_from_images( mp4_file, nc1D_png_dir, fps=fps, REPORT=True)
      
#----------------------------------------------------------------
# Create movie for a soil z-profile time series:  Darcy velocity
# at one of the monitored grid cells (var_index)
#----------------------------------------------------------------
V_PROFILE_MOVIE = False
if (V_PROFILE_MOVIE):
    os.chdir(output_dir)
    nc1D_png_dir = output_dir + 'v1D_png_files/'
    if not(os.path.exists( nc1D_png_dir )): os.mkdir( nc1D_png_dir)

    nc_file = case_prefix + '_1D-v.nc'
    ymin = None
    ymax = None
    # Set ymin and ymax to make them same for all plots.
    tfvis.save_profile_series_as_images(nc_file, png_dir=nc1D_png_dir,
               ymin=ymin, ymax=ymax,    # (same for all plots)
               # ymin=None, ymax=None,  # (auto for each plot)
               marker=',', REPORT=True,
               xsize=8, ysize=5, dpi=192)
    os.chdir(output_dir)                   
    fps = 20  # frames per second
    mp4_file = case_prefix + '_Profile_Movie_v.mp4'
    tfvis.create_movie_from_images( mp4_file, nc1D_png_dir, fps=fps, REPORT=True)
