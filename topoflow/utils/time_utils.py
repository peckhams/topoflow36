
#------------------------------------------------------------------------
#  Copyright (c) 2020, Scott D. Peckham
#
#  Jun 2020.  Created for TopoFlow calibration notebooks.
#             Started from code in: balto_gui.py and calibrate.py.
#
#------------------------------------------------------------------------

import datetime
import numpy as np
import sys

#------------------------------------------------------------------------
#
#  standardize_datetime_str()
#  get_datetime_obj_from_str()
#  get_datetime_obj_from_one_str()
#  pad_with_zeros()
#  split_datetime_str()
#  split_date_str()
#  split_time_str()
#  get_time_since_from_datetime()
#  get_datetime_from_time_since()
#  get_month_difference()
#
#  convert_times_from_hhmm_to_minutes()
#  convert_times_from_minutes_to_hhmm()
#  convert_times_from_datetime_to_minutes()
#  convert_times_from_minutes_to_datetime()
#
#  COMMENTED OUT FOR NOW.
#  get_actual_time_units()
#  get_time_delta_str()
#  get_start_datetime_obj()
#  get_end_datetime_obj()
#  get_dt_from_datetime_str()
#
#--------------------------------------------------------------------
def standardize_datetime_str( datetime_str ):

	#---------------------------------------------------
	# Note: Handle any delim between date & time, such
	#       as '', ' ' or 'T' (from Ankush).
	#---------------------------------------------------  
    date_str     = datetime_str[0:10]   # (2015-10-01)
    time_str     = datetime_str[-8:]    # (00:00:00)
    datetime_str = date_str + ' ' + time_str
    return datetime_str
    
#   standardize_datetime_str()
#--------------------------------------------------------------------
def get_datetime_obj_from_str( date_str, time_str='00:00:00' ):

    #---------------------------------------------------
    # date_str = 'YYYY-MM-DD', time_str = 'HH:MM:SS'
    #---------------------------------------------------
    ## e.g. d1 = str(self.datetime_end_date.value)
    ## e.g. t1 = self.datetime_end_time.value

    (y, m1, d) = split_date_str(date_str)
    (h, m2, s) = split_time_str(time_str)
    if( y <= 0 ):
        # msg  = 'Year cannot be < 1 in start date.\n'
        # msg += 'Changed year from ' + str(y) + ' to 1.'
        # self.datetime_notes.value = msg
        print('Year cannot be < 1 in start date.')
        print('Changed year from ' + str(y) + ' to 1.')
        print()
        y = 1
    datetime_obj = datetime.datetime(y, m1, d, h, m2, s) 
    return datetime_obj
    
#   get_datetime_obj_from_str()
#--------------------------------------------------------------------                       
def get_datetime_obj_from_one_str( datetime_str ):

    (date, time) = split_datetime_str( datetime_str )
    (y, m1,  d)  = split_date_str( date )
    (h, m2, s)   = split_time_str( time )
    datetime_obj = datetime.datetime(y, m1, d, h, m2, s)
    return datetime_obj

#   get_datetime_obj_from_one_str()
#------------------------------------------------------------------------
def pad_with_zeros(num, target_len):
  
    num_string = str( int(num) )  # int removes decimal part
    n = len( num_string )
    m = (target_len - n)
    num_string = ('0'*m) + num_string
    return num_string

#   pad_with_zeros()
#--------------------------------------------------------------------  
def split_datetime_str(datetime_obj, datetime_sep=' ',
                       ALL=False):

    #-----------------------------------------------     
    # Note: Still works if datetime_obj is string.
    #-----------------------------------------------
    datetime_str = str(datetime_obj)
    parts = datetime_str.split( datetime_sep )
    ## print('## datetime_str =', datetime_str )
    ## print('## parts =', str(parts) )
    
    date_str = parts[0]
    time_str = parts[1]
    if not(ALL):
        return (date_str, time_str)
    else:
        (y,m1,d) = split_date_str( date_str )
        (h,m2,s) = split_time_str( time_str )
        return (y,m1,d,h,m2,s)

#   split_datetime_str()
#--------------------------------------------------------------------  
def split_date_str(date_str, date_sep='-'):

    date_parts = date_str.split( date_sep )
    year  = int(date_parts[0])
    month = int(date_parts[1])   # NOTE:  int('08') = 8
    day   = int(date_parts[2])
    return (year, month, day)
     
#   split_date_str()
#--------------------------------------------------------------------  
def split_time_str(time_str, time_sep=':'):

    time_parts = time_str.split( time_sep )
    hour   = int(time_parts[0])
    minute = int(time_parts[1])
    second = int(time_parts[2])
    return (hour, minute, second)

#   split_time_str()
#--------------------------------------------------------------------                       
def get_time_since_from_datetime(origin_datetime_obj,
                                 datetime_obj, units='days'):

    #-------------------------------------------------
    # Compute time duration between datetime objects
    #-------------------------------------------------
    origin_obj    = origin_datetime_obj
    duration_obj  = (datetime_obj - origin_obj)
    duration_secs = duration_obj.total_seconds()
    #---------------------------------------------------
    # There is not a fixed number of seconds per month
    # Also 52 (weeks/year) * 7 (days/week) = 364.
    #---------------------------------------------------
    secs_per_unit_map = {
    'years':31536000.0, 'weeks':604800.0, 'days':86400.0,
    'hours':3600.0, 'minutes':60.0, 'seconds':1 }          
    secs_per_unit = secs_per_unit_map[ units ]      
    duration = (duration_secs / secs_per_unit )
    time_since = duration  # (in units provided)

    return time_since
        
#   get_time_since_from_datetime()
#--------------------------------------------------------------------                       
def get_datetime_from_time_since(origin_datetime_obj,
                           time_since, units='days'):
                                         
    # For testing
#         print('## type(times_since) =', type(time_since) )
#         print('## time_since =', time_since )
#         print('## int(time_since) =', int(time_since) )
    
    #---------------------------------------------------   
    # Note: datetime.timedelta() can take integer or
    #       float arguments, and the arguments can be
    #       very large numbers.  However, it does not
    #       accept any numpy types, whether float or
    #       int (e.g. np.int16, np.float32).
    #  https://docs.python.org/3/library/datetime.html
    #---------------------------------------------------
    delta = None
    time_since2 = float(time_since)  ## No numpy types
    #------------------------------------------------------    
    if (units == 'days'):
        delta = datetime.timedelta( days=time_since2 )
    if (units == 'hours'):
        delta = datetime.timedelta( hours=time_since2 )
    if (units == 'minutes'):
        delta = datetime.timedelta( minutes=time_since2 )
    if (units == 'seconds'):
        delta = datetime.timedelta( seconds=time_since2 )
    #------------------------------------------------------
    if (delta is None):
        msg = 'ERROR: Units: ' + units + ' not supported.'
        return

    # For testing
    ## print('#### delta =', delta)
    
    #---------------------------------------------        
    # Create new datetime object from time_since
    #---------------------------------------------
    origin_obj = origin_datetime_obj
    new_dt_obj = (origin_obj + delta)
    return new_dt_obj

#   get_datetime_from_time_since()
#--------------------------------------------------------------------  
def get_month_difference(start_datetime_obj, end_datetime_obj ):

    #-------------------------------------------
    # Example 0: 2017-09 to 2017-09
    # months = (2017-2017)*12 = 0
    # months = (months - 9) = (0-9) = -0
    # months = (months + 9) = 0   (as index)
    #-------------------------------------------           
    # Example 1: 2017-09 to 2018-02
    # 9:10, 10:11, 11:12, 12:1, 1:2 = 5 (if same days)
    # months = (2018-2017)*12  = 12
    # months = (months - 9) = 3
    # months = (months + 2) = 3 + 2 = 5
    #-------------------------------------------
    start_year = start_datetime_obj.year
    end_year   = end_datetime_obj.year
    months     = (end_year - start_year) * 12
    #-------------------------------------------
    start_month = start_datetime_obj.month
    end_month   = end_datetime_obj.month
    months = months - start_month
    months = months + end_month
    ## months = months + 1  # (no: get 1 if dates same)
    ## print('month difference =', months)
    return months

#   get_month_difference()
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def convert_times_from_hhmm_to_minutes( times_hhmm ):

    #----------------------------------------------------      
    # Notes:  For simple time strings with only 'hhmm',
    #         return the time in minutes.
    #----------------------------------------------------    
    n_times   = len( times_hhmm )
    times_min = np.zeros( n_times )
    for k in range(n_times):
        hhmm = times_hhmm[k]
        hour = np.int16( hhmm[:2] )
        min  = np.int16( hhmm[2:] )
        times_min[k] = (hour * 60) + min

    return times_min

#   convert_times_from_hhmm_to_minutes()
#---------------------------------------------------------------------
def convert_times_from_minutes_to_hhmm( times_min ):

    n_times   = len( times_min )
    times_hhmm = np.zeros( n_times )
    for k in range(n_times):
        hour = int( times_min[k] / 60 )
        min  = int( times_min[k] % 60 )
        #-----------------------------------        
        hh   = str(hour)
        if (len(hh) == 1):  hh = ('0' + hh)
        #-----------------------------------
        mm   = str(min)
        if (len(mm) == 1):  mm = ('0' + mm)
        #-----------------------------------
        times_hhmm[k] = int(hh + mm)

    return times_hhmm

#   convert_times_from_minutes_to_hhmm()
#---------------------------------------------------------------------
def convert_times_from_datetime_to_minutes( times_datetime,
                  origin_datetime_obj=None ):
    
    if (origin_datetime_obj is None):
        print('ERROR: origin_datetime_obj is required.')
        return 

    n_times   = len( times_datetime )
    times_min = np.zeros( n_times )
    for k in range(n_times):
        datetime_str = times_datetime[k]
        #---------------------------------------------------
        # Note: Handle any delim between date & time, such
        #       as '', ' ' or 'T' (from Ankush).
        #--------------------------------------------------- 
        datetime_str = standardize_datetime_str( datetime_str )        
        datetime_obj = get_datetime_obj_from_one_str( datetime_str )
        times_min[k] = get_time_since_from_datetime(origin_datetime_obj,
                                datetime_obj, units='minutes')

    return times_min

#   convert_times_from_datetime_to_minutes()
#---------------------------------------------------------------------
def convert_times_from_minutes_to_datetime( times_min,
                  origin_datetime_obj=None ):
    
    if (origin_datetime_obj is None):
        print('ERROR: origin_datetime_obj is required.')
        return 

    n_times   = len( times_min )
    times_datetime = np.zeros( n_times, dtype='U30' )   # string array
    for k in range(n_times):
        time_since = times_min[k]
        times_datetime[k] = get_datetime_from_time_since(origin_datetime_obj,
                                         time_since, units='minutes')

    return times_datetime

#   convert_times_from_minutes_to_datetime()
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#     def get_actual_time_units(self):
# 
# #         secs_per_unit_list = [1, 60.0, 3600.0, 86400, 31536000.0, -1]
# #         next_unit_factor = [60.0, 60.0, 24.0, 365.0, -1, -1]
# 
#         units_list = ['second', 'minute', 'hour',
#                        'day', 'year', 'None']   # ascending, skip month
# 
#         for units in units_list:
#              if (self.time_units_str.startswith(units)):
#                  break
#         if (units != None):
#             units += 's'   # (make units plural now; not before)
#         else:
#             print('ERROR: No match found for units.')
#             return
#         self.time_units = units
# 
#     #   get_actual_time_units()
#     #--------------------------------------------------------------------
#     def get_time_delta_str(self):
# 
#         ## print('### self.time_var.size =', self.time_var.size )
#         ## print('###')
#         
#         #-----------------------------------
#         # Check size of the time_var array
#         #-----------------------------------
#         if (self.time_var.size == 1):
#             dt = 0
#             self.time_delta = '0000-00-00 00:00:00'
#             # print('At top of get_time_delta_str():')
#             # print('self.time_var.size =', self.time_var.size )
#             # print('self.time_delta =', self.time_delta )
#             return
#         if (self.time_var.size > 1):  
#             dt  = (self.time_var[1] - self.time_var[0])
#             print('dt1 =', dt)
#         if (self.time_var.size > 3):
#             dt2 = (self.time_var[2] - self.time_var[1])  ###
#             dt3 = (self.time_var[3] - self.time_var[2])  ###
#             print('dt2 =', dt2)  # check if evenly spaced
#             print('dt3 =', dt3)
#                 
#         #---------------------------------------------------        
#         # Note: Actual time units were stripped from units
#         #       string and saved as self.time_units.
#         #       A full units attribute string may be:
#         #        'hour since 0000-00-00 00:00:00'
#         #---------------------------------------------------
#         units_list = ['seconds', 'minutes', 'hours',
#                       'days', 'years', 'None']  # ascending, skip month
#         secs_per_unit_list = [1, 60.0, 3600.0, 86400, 31536000.0, -1]
#         next_unit_factor   = [60.0, 60.0, 24.0, 365.0, -1, -1]
#         units       = self.time_units
#         units_index = units_list.index( units )
#         #----------------------------------------
#         if (units == 'years'):
#             s = self.pad_with_zeros(dt,4)
#         else:
#             if (len(str(dt)) <= 2):
#                 s = self.pad_with_zeros(dt,2)
#             else:
#                 #-------------------------------
#                 # Must convert units to get dt
#                 # down to 1 or 2 digits.
#                 #-------------------------------
#                 old_dt    = dt
#                 old_units = units
#                 k = units_index
#                 n = len( str(int(dt)) )
#                 while (n > 2) and (units != 'None'):
#                     k     = k + 1
#                     dt    = (dt / next_unit_factor[k-1])
#                     units = units_list[k]
#                     n     = len( str(int(dt)) )
#                 if (units == 'None'):
#                     print('#####################################')
#                     print('ERROR in get_time_delta_str():')
#                     print('      dt has too many digits.')
#                     print('#####################################')
#                     return
#                 else:
#                     # Note that any remainder has been dropped.
#                     s = self.pad_with_zeros(dt,2)
#                     print('Old dt and units =', old_dt, old_units)
#                     print('New dt and units =', dt, units)
#                     print('Remainder not retained yet.')
#         #----------------------------------------------
#         if (units == 'years'):
#             td = (s + '-00-00 00:00:00')
# #         if (units == 'months'):
# #             td= ('0000-' + s + '-00 00:00:00')
#         if (units == 'days'):
#             td = ('0000-00-' + s + ' 00:00:00')
#         if (units == 'hours'):
#             td = ('0000-00-00 ' + s + ':00:00')
#         if (units == 'minutes'):
#             td = ('0000-00-00 00:' + s + ':00')
#         if (units == 'seconds'):
#             td = ('0000-00-00 00:00:' + s)
#         #------------------------------------------------
#         self.time_delta = td
#         # print('At bottom of get_time_delta_str():')
#         # print('self.time_delta =', td)
#         # print()
# 
#     #   get_time_delta_str()
#     #--------------------------------------------------------------------
#     def get_start_datetime_obj(self):
# 
#         #---------------------------------------
#         # d1.value is a datetime "date object"
#         # t1.value is a time string: 00:00:00
#         #---------------------------------------
#         d1 = self.datetime_start_date
#         t1 = self.datetime_start_time
#         if (d1.value is None):
#             return None
#         date_str = str(d1.value)
#         time_str = t1.value   # (already string)
#         ## print('In get_start_datetime_obj():')
#         ## print('date_str =', date_str)
#         ## print('time_str =', time_str)
#         
#         datetime_obj = self.get_datetime_obj_from_str(date_str, time_str)
#         return datetime_obj
#     
#     #   get_start_datetime_obj()
#     #--------------------------------------------------------------------
#     def get_end_datetime_obj(self):
# 
#         #---------------------------------------
#         # d1.value is a datetime "date object"
#         # t1.value is a time string: 00:00:00
#         #---------------------------------------
#         d1 = self.datetime_end_date
#         t1 = self.datetime_end_time
#         if (d1.value is None):
#             return None
#         date_str = str(d1.value)
#         time_str = t1.value   # (already string)
#         ## print('In get_end_datetime_obj():')
#         ## print('date_str =', date_str)
#         ## print('time_str =', time_str)
#         
#         datetime_obj = self.get_datetime_obj_from_str(date_str, time_str)
#         return datetime_obj
# 
#     #   get_end_datetime_obj()
#---------------------------------------------------------------------

