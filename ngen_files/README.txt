S.D. Peckham
September 13, 2022
July 3, 2024 (updated)

For use in NextGen, the topoflow36 project folder should be copied into
the NextGen project folder, at ngen/extern/topoflow36.  Please see the
text file:  "How_to_Run_TopoFlow_in_NextGen.txt" in the docs folder for
much more information.

This subfolder, "ngen_files", contains files that should be copied into
the ngen repo tree (or project folder), into the indicated subfolders.
This will allow the topoflow36 Python package to be run from within NextGen.

The "ngen_files/data/topoflow" folder contains:
(1) a "realization config" file for NextGen in the "rc_files" folder
(2) catchment and nexus data GeoJSON files in the "spatial" folder,
    including ones for HUCO1, CAMELS-test, etc.
(3) all of the input files required by TopoFlow_3.6, for several catchments
    in the "input_files" folder (DEM, soil data, met data, CFG files, etc.) 

---------------------------------------------------
NOTES REGARDING THE METEOROLOGICAL FORCING FILES:
---------------------------------------------------
As indicated in the file: "How_to_Run_TopoFlow_in_NextGen.txt", you can
begin by using the example files provided in the "ngen_files" folder.
However, to perform runs for other NextGen catchments you will need to
download meteorological forcing files in CSV format and copy them into
the directory:  ngen/data/topoflow/forcing/huc01. 

Meteorological forcing files are retrieved from this "source" directory by the
prepare_inputs.py module when it creates all of the input files for a given
catchment and places them in a directory like:
   ngen/data/topoflow/input_files/cat-67/__met.

Note:  In a previous version, all of the input files required by TopoFlow_3.6
for several catchments were retrieved from:  ngen/extern/topoflow36/data.
However, the folder: ngen/extern/topoflow36 should only contain files that
are downloaded as part of the topoflow36 repo, so data files should now be
retrieved from:  ngen/data/topoflow/input_files.


