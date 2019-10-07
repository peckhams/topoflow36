
S.D. Peckham
October 6, 2019

====================================================
 How to Prepare Input Files for the TopoFlow model
====================================================
Step 1.
Obtain a DEM for the region to be modeled.
   -- elevation grid                 [site_prefix]_DEM0.rtg

Step 1b.
(This is optional, but is recommended for mature landscapes.)
Use topoflow/components/smooth_DEM.py to create "profile-smoothed" DEM.
   -- elevation grid                 [site_prefix]_DEM.rtg
Alternately, use topoflow/utils/new_slopes.py to create a new slope grid.

Step 2.
Use topoflow/components/d8_global.py to create set of D8 files from DEM.
   -- D8 flow direction grid         [site_prefix]_flow.rtg
   -- channel slope grid             [site_prefix]_slope.rtg
   -- total contributing area grid   [site_prefix]_area.rtg
   -- D8 flow width grid             [site_prefix]_dw.rtg
   -- monitored basin outlets        [case_prefix]_outlets.txt  May need GIS.

Step 3.
- Copy all of the component "configuration files" (CFG files) from the folder
  "component_cfg_files" into the model's working directory.
- Edit the CFG file called: [case_prefix]_path_info.cfg to change site_prefix,
  case_prefix and the path to TopoFlow's input and output directories.
- Edit the CFG files to choose options.
- Each CFG file also has a setting that allows the component to be Disabled.
- When TopoFlow is used with a GUI, the GUI collects the necessary
  information from the user and writes a set of CFG files like this.

Step 4.
Use topoflow/utils/parameterize.py to create:
   -- channel bankfull depth grid    [site_prefix]_d-bank.rtg
   -- channel bottom width grid      [site_prefix]_chan-w.rtg
   -- channel manning n grid         [site_prefix]_chan-n.rtg
   -- channel flare angle grid       [site_prefix]_chan-a.rtg  (or use constant)
   -- channel abs. sinuosity grid    [site_prefix]_sinu.rtg

Step 5.
Use topoflow/utils/init_depth.py to create:
   -- initial channel depth grid     [site_prefix]_d0.rtg

Step 6.
Obtain meteorological data input grids (e.g. from GES DISC website)
   -- rainfall rates                 [case_prefix]_rainrate.rts
   -- relative humidity              [case_prefix]_RH.rts
   -- air temperature                [case_prefix]_T_air.rts
   -- atmos. air pressure            [case_prefix]_p0.rts

Step 7.
Obtain 4 soil property grids (e.g. from ISRIC website) at various depths.
   -- soil_bulk__mass-per-volume_density
   -- soil_clay__volume_fraction  (percent clay)
   -- soil_matter~organic__volume_fraction  (percent organic matter)
   -- soil_silt__volume_fraction  (percent silt)

Step 8.
Use topoflow/utils/pedotransfer.py to create soil hydraulic property grids.
   -- saturated hydraulic conductivity        [site_prefix]_Ks.rtg
   -- saturated soil water volume fraction    [site_prefix]_qs.rtg
   -- Brooks-Corey lambda parameter           [site_prefix]_lam.rtg
   -- Brooks-Corey-Smith c parameter          [site_prefix]_c.rtg
   -- bubbling pressure head                  [site_prefix]_pB.rtg
   -- Green-Ampt capillary length             [site_prefix]_G.rtg
   --------------------------------------------------------------------
   -- initial soil hydraulic conductivity     [site_prefix]_Ki.rtg
   -- initial soil water volume fraction      [site_prefix]_qi.rtg
   -- residual soil water volume fraction     [site_prefix]_qr.rtg

Step 9.
Create a grid of initial elevation of the groundwater table.
This could be estimated by substracting a constant from the DEM.
   -- [case_prefix]_h0-table.rtg

Step 10.
Edit the file "[case_prefix]_providers.txt" to select which components to use for
each hydrologic process.  This file is read by the EMELI framework.  You should
use a different case_prefix whenever you change the input files or parameters,
or when you plan to run TopoFlow with a different set of model components.


