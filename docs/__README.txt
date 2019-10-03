
June 19, 2013
November 4, 2013
S.D. Peckham

===================================
 BMI Interface and Standard Names
===================================
This version of TopoFlow has fully-compliant BMI interfaces
for all of the components, with the latest CSDMS standard names.
Each component inherits from BMI_base.py, in "utils" folder.

==========================================
 Installing TopoFlow as a Python Package
==========================================
(Developer version (editable) vs. regular version) ######
To install TopoFlow as a Python package:

#########

If you want to run the unit tests that are included in the new TopoFlow
package, you should also create a folder in your home directory called:
"TopoFlow_Tests" with a subfolder called "Test1".

========================================================
 Running TopoFlow in the Light-weight Python Framework
========================================================
TopoFlow can be run outside of the CSDMS modeling framework by using
the light-weight Python framework called framework.py.  This framework
is now included as part of the TopoFlow package, in the "framework"
folder. This requires creating a "provider_file", based on the example
in the framework folder.

===============
 D8 Component
===============
Right now, TopoFlow uses tf_d8_base.py, instead of d8_base.py.
The latter is located in the "utils" folder.  Need to find
the most powerful and well-tested d8_base.py, possibly in the
Erode3 (Global) project, and make sure TopoFlow can use it.

===================
 Remaining Issues
===================
HIS_base.py needs to be updated with BMI and new standard names.



