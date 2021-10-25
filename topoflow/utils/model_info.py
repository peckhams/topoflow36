
#-------------------------------------------------------------------
#  Copyright (c) 2021, Scott D. Peckham
#
#  Oct. 2021.  Added to bypass old tf_utils.py for version, etc.
#
#-------------------------------------------------------------------

#   name()
#   version_number()
#   build_date()
#   version_string()
#
#-------------------------------------------------------------------   
def name():
 
    return 'TopoFlow'

#   name()
#-------------------------------------------------------------------
def version_number():

    return 3.64

#   version_number()
#-------------------------------------------------------------------   
def build_date():

    #---------------------------------------------------------
    # Notes:  Update this whenever a new version is released
    #---------------------------------------------------------
    return '2021-25-13'

#   build_date()
#-------------------------------------------------------------------
def version_string():

    num_str    = str( version_number() )
    date_str   = ' (' + build_date() + ')'

    ver_string = name() + ' Version ' + num_str + date_str

    return ver_string

#   version_string()
#-------------------------------------------------------------------


