
S.D. Peckham
November 8, 2013

The Python code in this "gui" directory provides a partial reproductionof the TopoFlow IDL GUI, using wxPython.  All of the wizard panels and
input/output variable panels are implemented but don't do anything yet.
Additional dialogs in the Create menu, etc. are not yet implemented.To launch the GUI, cd to this directory and type:    % python Main_Dialog.py
It must be launched this way in order for the Exit button and File >
Close to work.

Note that the menu bar appears at the top of the Mac screen, as with
other Mac applications.Note that the bitmaps, dialogs and help folders in this directory contain
the required files, and this version looks for them in /Applications/TopoFlow.
However, they should be stored and found in the gui subpackage instead,
similar to how framework.get_package_path() works.
