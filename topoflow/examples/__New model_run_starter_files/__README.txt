
S.D. Peckham
October 5, 2019

====================================================
 How to Prepare Input Files for the TopoFlow model
====================================================
Step 1.
Obtain a DEM for the region to be modeled.
   -- elevation grid                 [site_prefix]_DEM0.rtg

Step 1b. (Optional, but recommended)
Use topoflow/components/smooth_DEM.py to create "profile-smoothed" DEM.
   -- elevation grid                 [site_prefix]_DEM.rtg

Step 2.
Use topoflow/components/d8_global.py to create set of D8 files from DEM.
   -- D8 flow direction grid         [site_prefix]_flow.rtg
   -- channel slope grid             [site_prefix]_slope.rtg
   -- total contributing area grid   [site_prefix]_area.rtg
   -- D8 flow width grid             [site_prefix]_dw.rtg
   -- monitored basin outlets        [case_prefix]_outlets.txt  May need GIS.

Step 3.
Use topoflow/utils/parameterize.py to create:
   -- channel bankfull depth grid    [site_prefix]_d-bank.rtg
   -- channel bottom width grid      [site_prefix]_chan-w.rtg
   -- channel manning n grid         [site_prefix]_chan-n.rtg
   -- channel flare angle grid       [site_prefix]_chan-a.rtg
   -- channel abs. sinuosity grid    [site_prefix]_sinu.rtg

Step 4.
Use topoflow/utils/init_depth.py to create:
   -- initial channel depth grid     [site_prefix]_d0.rtg

Step 4.
Obtain meteorological data input grids (e.g. from GES DISC website)
   -- rainfall rates                 [case_prefix]_rainrate.rts
   -- relative humidity              [case_prefix]_RH.rts
   -- air temperature                [case_prefix]_T_air.rts
   -- atmos. air pressure            [case_prefix]_p0.rts

Step 5.
Obtain 4 soil property grids (e.g. from ISRIC website) at various depths.
   -- soil_bulk__mass-per-volume_density
   -- soil_clay__volume_fraction  (percent clay)
   -- soil_matter~organic__volume_fraction  (percent organic matter)
   -- soil_silt__volume_fraction  (percent silt)

Step 6.
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

Step 7.
- Copy all of the component "configuration files" (CFG files) from the folder
  "component_cfg_files" into the model's working directory.
- Edit the CFG files o choose options.
- Each CFG file also has a setting that allows the component to be Disabled.
- The CFG file called: [case_prefix]_path_info.cfg can be edited to change
  the path to TopoFlow's input and output directories.

Step 8.
Edit the file "Test1_providers.txt" to select which components to use for
each hydrologic process.  This file is read by the EMELI framework.


