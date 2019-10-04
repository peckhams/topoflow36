# topoflow36
TopoFlow version 3.6 is an update to version 3.5 that runs in Python 3.7.
TopoFlow 3.5 is in the "topoflow" repository and runs in Python 2.7. 
Both are Python packages that consist of over 72,000 lines of code.

The best way to learn about TopoFlow --- including its history, capabilities, input data preparation and an example application --- is to read the paper in the "docs" folder called:  "Peckham_et_al_2017_GPF.pdf", along with its appendices "Peckham_et_al_2017_GPF_Appendices.pdf".

Detailed information on each TopoFlow model component --- including its variables and the equations used --- can be found on the CSDMS wiki website at:  https://csdms.colorado.edu/wiki/Model:TopoFlow.  Each component has its own very detailed HTML Help Page.  Links to these can be found near the bottom of the main CSDMS page for TopoFlow.

TopoFlow versions 3.5 and 3.6 provide a plug-and-play, component-based modeling framework for spatial hydrologic modeling.  Both use the EMELI framework to couple a set of user-selected components into a functioning model.  The EMELI framework is included in the topoflow "framework" folder.  A paper on EMELI is included in the "docs" folder, called:  "Peckham_2014_EMELI_FINAL.pdf".

Model components are in the topoflow "components" folder and each have a Basic Model Interface (BMI).  They all inherit from BMI_base.py in the "utils" folder.  The "utils" folder has a large collection of utilities that are used by all components.
