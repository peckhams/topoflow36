
SILENT = False
if not(SILENT):
    print('Importing TopoFlow 3.6 packages:')
    print('   topoflow.utils')
    print('   topoflow.utils.tests')
    print('   topoflow.components')
    print('   topoflow.components.tests')
    print('   topoflow.framework')
    print('   topoflow.framework.tests')
    # print '   topoflow.gui (unfinished)'
    print(' ')

import topoflow.utils
import topoflow.utils.tests
#-----------------------------------
import topoflow.components
import topoflow.components.tests
#-----------------------------------
import topoflow.framework
import topoflow.framework.tests

#------------------------------------------
# This imports EMELI and prints more info
#------------------------------------------
import topoflow.framework.tests.test_framework

#--------------------------------------
# The Python GUI is not yet finished.
#--------------------------------------
# import topoflow.gui

#----------------------------------------------
# Idea:  Import the topoflow_test() function,
#        so it can be executed with:
#    >>> import topoflow
#    >>> topoflow.topoflow_test().
#----------------------------------------------
# This doesn't work.
# from topoflow.framework.tests.test_framework import topoflow_test

