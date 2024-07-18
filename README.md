# TopoFlow

**Description**:  
TopoFlow provides a plug-and-play, component-based modeling framework for spatial hydrologic modeling.  Both use the EMELI framework to couple a set of user-selected components into a functioning model.  The EMELI framework is included in the topoflow "framework" folder. [A paper on EMELI](https://github.com/NOAA-OWP/topoflow/blob/master/docs/Peckham_2014_EMELI_FINAL.pdf) is included in the "docs" folder.

TopoFlow was adapted to model combined snow and ice melt contributions to streamflow in glaciated catchments as part of NextGen modeling efforts, but TopoFlow also provides several other modules related to channel routing, potential evapotranspiration and infiltration estimates, and more. 

Model components are in the topoflow "components" folder and each have a Basic Model Interface (BMI).  Details on BMI can be found here:  https://csdms.colorado.edu/wiki/BMI_Description.  They all inherit from BMI_base.py in the "utils" folder.  Each model component has its own configuration file, a text file with extension ".cfg".  Components may get their input data from other components in a component set (passed by reference), or directly from input files.  By editing model CFG files, users specify which output variables are to be saved to files (and at what "sampling timestep) and values of variables are then written to netCDF files (and optionally to some other file formats). 
Detailed information on each TopoFlow model component --- including its variables and the equations used --- can be found on the CSDMS wiki website at:  https://csdms.colorado.edu/wiki/Model:TopoFlow.  

The "utils" folder has a large collection of utilities that are used to create and read input files, write output files and provide other low-level functionality to all components.  It also includes a complete toolkit for computing a D8 flow direction grid and all related variables (like slope and total contributing area) from a DEM.  The "examples" folder contains several data sets for testing, including the Treynor watershed in Iowa with input data for two historic rainfall events both from June 1967.  Included in these examples are folders showing correct output files that can be used for testing.

  - **Technology stack**: TopoFlow that runs in Python 3.7 or later and contains over 98,000 lines of Python code. The model can be run either in standalone or as part of the [NextGen Framework](https://github.com/NOAA-OWP/ngen). 
  - **Status**:  TopoFlow is currently able to run in standalone. TopoFlow has been run in the NextGen Framework, but a standarization process and associated documentation for users is currently in development. 

## Dependencies

TopoFlow requires a Python version of 3.7 or later. The latest version of Python the model has been tested with is 3.9.

## Installation

Detailed instructions on how to install TopoFlow can be found at [INSTALL](INSTALL.md).

## Configuration

If the software is configurable, describe it in detail, either here or in other documentation to which you link.

## Usage

Show users how to use the software.
Be specific.
Use appropriate formatting when showing code snippets.

## How to test the software

If the software includes automated tests, detail how to run those tests.

## Known issues

Document any known significant shortcomings with the software.

## Getting help

Instruct users how to get help with this software; this might include links to an issue tracker, wiki, mailing list, etc.

**Example**

If you have questions, concerns, bug reports, etc, please file an issue in this repository's Issue Tracker.

## Getting involved

This section should detail why people should get involved and describe key areas you are
currently focusing on; e.g., trying to get feedback on features, fixing certain bugs, building
important pieces, etc.

General instructions on _how_ to contribute should be stated with a link to [CONTRIBUTING](CONTRIBUTING.md).


----

## Open source licensing info
1. [TERMS](TERMS.md)
2. [LICENSE](LICENSE)


----

## Credits and references

1. Related projects:

CSDMS (csdms.colorado.edu) uses an older version of TopoFlow and also provides a graphical user interface via its Web Modeling Tool.  See:  https://csdms.colorado.edu/wiki/CSDMS_Web_Modeling_Tool  (and https://csdms.colorado.edu/wiki/HydrologyExecutorBlanca).  Incomplete efforts to create a wizard-style GUI for TopoFlow in Python --- similar to the GUI used in the old IDL version of TopoFlow 1.6b --- can be found in the "gui_old" folder.
2. Books, papers, talks, or other sources that have meaningful impact or influence on this project:

The best way to learn about TopoFlow --- including its history, capabilities, input data preparation and an example application --- is to read the paper [Peckham et al. (2017)](https://github.com/NOAA-OWP/topoflow/blob/master/docs/Peckham_et_al_2017_GPF.pdf) in the "docs" folder, along with its [appendices](https://github.com/NOAA-OWP/topoflow/blob/master/docs/Peckham_et_al_2017_GPF_Appendices.pdf).  One of the appendices explains how to install the Python package and run a test.  It is convenient to first install the Anaconda Python distribution, create a conda environment like "tf36", and then install the netCDF4 package dependency in that environment.