
## Copyright (c) 2009, Scott D. Peckham
## January 16, 2009

## from numpy import *
import numpy as np

from .tf_utils import *

from . import idl_func
import sys, os, time, glob, platform

#******************************************************************
#   save_vars.pro

#   Copyright (c) 2007-2008, Scott D. Peckham
#   Created:  (moved from getvars.pro, extended) Mar 2007
#   Modified:  Mar 2008  (Any number of layers for Richards, etc.)

#******************************************************************

#   Save_Var            (7/24/06)
#   Save_Var2           (7/25/06)
#   Save_Var3           (7/25/06)
#   Save_All_TF_Vars    (7/25/06)
#---------------------------------
#   Read_Vars           (7/25/06)    ;(Name similar to Read_Var )
#   Load_Var
#   Type_Code

#******************************************************************
def Save_Var(LUN, name, var, _TYPE=None, units=None):

    n_params = 5 - [_TYPE, units].count(None)
    fF = '(F12.6)'   #(WILL '(F0)' PRINT ALL DIGITS ?? **************)
    fD = '(D18.8)'
    
    if (idl_func.n_elements(_TYPE) == 0):    
        _TYPE = 'STRING'
    
    #--------------------------
    #Convert "var" to a string
    #--------------------------
    I2PY_expr = _TYPE.upper()
    if I2PY_expr == 'BYTE':    
        str2 = TF_String(np.int16(var))
    elif I2PY_expr == 'FLOAT':    
        str2 = TF_String(var, FORMAT=fF)
    elif I2PY_expr == 'DOUBLE':    
        str2 = TF_String(var, FORMAT=fD)
    elif I2PY_expr == 'STRING':    
        str2 = var
    else:    
        str2 = TF_String(var)
    
    
    #------------------------
    #Left justify the phrase
    #------------------------
    n = np.int16(24)                 #(width of name section)
    _len = len(name)
    m = np.maximum((n - _len), 1)
    str1 = name + str(np.repeat(np.np.uint8(32),m))
    
    #-------------------------
    #Add optional unit string
    #-------------------------
    if (idl_func.n_elements(units) != 0):    
        str3 = units
    else:    
        str3 = ''
    
    #-------------------------------
    #Print phrase and value to file
    #-------------------------------
    file_LUN.write((str1 + str2 + str3))
    
#  Save_Var
#******************************************************************
def Save_Var2(LUN, name, _TYPE, code, ptr, file, units, \
              NOT_PTR=None, FACTOR=None):

#------------------------------------------------------------------
#Notes:  LUN   = logical unit number of output file
#        name  = variable name or description as string ('T_air:')
#        type  = data type ('BYTE', 'FLOAT', 'DOUBLE', etc.)
#        code  = input type code (0=Scalar, 1=Time series, etc.)
#        ptr   = pointer to the heap variable (for scalars)
#        file  = name of file that contains variable (for others)
#        units = units of the variable, as a string

#        The NOT_PTR keyword can be set to indicate that the
#        "ptr" argument is a scalar value (e.g. pv.dt).

#        (2/12/07) Added FACTOR keyword for rainrates.
#------------------------------------------------------------------
    n_params = 7
    NOT_PTR = (NOT_PTR not in [0,None])
    if logical_not((FACTOR not in [0,None])):    
        FACTOR = np.float64(1)
    
    fF = '(F12.6)'   #(WILL '(F0)' PRINT ALL DIGITS ?? ************)
    fD = '(D18.8)'
    
    types = array(['Scalar         ', 'Time_series    ', 'Grid           ', 'Grid_Sequence  '])
    
    #-------------------------------
    #Left justify the variable name
    #-------------------------------
    n = np.int16(24)                 #(width of name section)
    _len = len(name)
    m = np.maximum((n - _len), 1)
    str1 = name + str(np.repeat(np.np.uint8(32),m))
    
    #--------------------------
    #Get the input type string
    #--------------------------
    str2 = types[code]
    
    #------------------------------
    #Get scalar string or filename
    #------------------------------
    if (code == np.uint8(0)):    
        if (NOT_PTR):    
            var = ptr
        else:    
            var = ptr
        
        #-----------------------------
        #Option to change units, etc.
        #-----------------------------
        if (FACTOR != np.float64(1)):    
            var = (var * FACTOR)
        
        #------------------------------------------
        #Make sure var is a scalar (for rainrates)
        #------------------------------------------
        #nv = n_elements(var)
        #if (nv gt 1) then begin
        #    var  = var[0]
        #msg  = ['WARNING: ', ' ',$
        #       'The list of rainrates and durations is not saved when ',$
        #       'the "Uniform in space" method for precipitation is used.',$
        #       'Only the first value in the list will be saved. ',$
        #       ' ',$
        #       'If you use the "General data types" method instead, ',$
        #       'then File > Load Input Vars can restore all settings.',$
        #       ' ']
        #result = GUI_Message(msg, /INFO)
        #str2 = 'Time_series    '
        #endif
        
        #---------------------------
        #Convert scalar to a string
        #---------------------------
        I2PY_expr = _TYPE.upper()
        if I2PY_expr == 'BYTE':    
            str3 = TF_String(np.int16(var))
        elif I2PY_expr == 'FLOAT':    
            str3 = TF_String(var, FORMAT=fF)
        elif I2PY_expr == 'DOUBLE':    
            str3 = TF_String(var, FORMAT=fD)
        else:    
            str3 = TF_String(var)
        
    else:    
        str3 = file
    
    #------------------
    #Left justify str3
    #------------------
    n3 = np.int16(24)               #(width of str3 section)
    len3 = len(str3)
    m3 = np.maximum((n3 - len3), 1)  #********
    str3 = str3 + str(np.repeat(np.np.uint8(32),m3))
    
    #-------------------------------
    #Print name, type, etc. to file
    #-------------------------------
    file_LUN.write((str1 + str2 + str3 + units))
    
#  Save_Var2
#******************************************************************
def Save_Var3(LUN, name, flag, file, units):

    #---------------------------------------------------------------
    # Notes: LUN   = logical unit number of output file
    #        name  = var. name or description as string ('T_air:')
    #        flag  = 1 = save this var, 0 = don't save
    #        file  = name of output file
    #        units = units of the variable, as a string
    #---------------------------------------------------------------
    n_params = 5
    if (file == ''):    
        file = 'n/a'
    
    #-----------------------------
    #Left justify the option name
    #-----------------------------
    n = np.int16(24)                 #(width of name section)
    _len = len(name)
    m = (n - _len)
    
    str1 = name + str(np.repeat(np.np.uint8(32),m))
    
    #--------------------------
    #Get the option flag string
    #---------------------------
    str2 = TF_String(np.int16(flag)) + str(np.repeat(np.np.uint8(32),14))
    
    #-------------------------------
    #Left-justified filename string
    #-------------------------------
    n3 = np.int16(24)               #(width of str3 section)
    len3 = len(file)
    m3 = np.maximum((n3 - len3), 1)
    str3 = file + str(np.repeat(np.np.uint8(32),m3))
    
    #--------------------------------------
    #Print name, flag and filename to file
    #--------------------------------------
    file_LUN.write((str1 + str2 + str3 + units))
    
#  Save_Var3
#******************************************************************
##def Save_All_TF_Vars(mstate):
##
##    #-----------------------------
##    #Get name of output text file
##    #-----------------------------
##    n_params = 1
##    run_prefix = mstate.run_vars.run_prefix
##    filename = ('00_' + run_prefix + '_INPUT.txt')
##    out_file = I2PY_filepath = []
##    app = wx.PySimpleApp()
##    I2PY_dialog = wx.FileDialog(parent=None, defaultDir=os.getcwd(), defaultFile=filename, style=wx.SAVE)
##    if (I2PY_dialog.ShowModal() == wx.ID_OK):
##        I2PY_filepath.append(I2PY_dialog.GetPath())
##    I2PY_dialog.Destroy()
##    #***** /OVERWRITE_PROMPT)
##    if (out_file == ''):    
##        return _ret()
##    Check_Overwrite(out_file, OK)
##    if logical_not(OK):    
##        return _ret()
##    
##    #------------------------
##    #Open text file to write
##    #------------------------
##    unit = TF_Get_LUN(out_file)
##    file_unit = open(out_file, 'w')
##    I2PY_SWAP_ENDIAN = False
##    
##    #--------------------
##    #Local abbreviations
##    #--------------------
##    av = mstate.grid_vars
##    rv = mstate.run_vars
##    pv = mstate.precip_vars
##    cv = mstate.channel_vars
##    mv = mstate.met_vars        #(new, 3/13/07)
##    sv = mstate.snow_vars
##    ev = mstate.ET_vars
##    gv = mstate.GW_vars
##    iv = mstate.infil_vars
##    dv = mstate.diversion_vars
##    fv = mstate.stop_vars
##    
##    #---------------------
##    #Prepare for printing
##    #---------------------
##    line = str(np.repeat(np.np.uint8(45),50))
##    #*** shortline = string(replicate(45b, 22))
##    
##    #-------------------
##    #Print out a header
##    #-------------------
##    file_unit.write('TopoFlow 1.5 Input File')
##    file_unit.write(' ')
##    
##    #---------------------------------------
##    #Read run vars from the "run var" panel
##    #----------------------------------------------
##    #Run vars are not saved until the user goes to
##    #the next panel, and they may not have done so
##    #----------------------------------------------
##    widget_control(rv.prefix_ID, get_value=prefix)
##    widget_control(rv.run_prefix_ID, get_value=run_prefix)
##    widget_control(rv.directory_ID, get_value=directory)
##    widget_control(rv.log_file_ID, get_value=log_file)
##    widget_control(rv.comment_file_ID, get_value=comment_file)
##    #------------------------------------
##    #Alternate approach (see note above)
##    #------------------------------------
##    #;prefix       = rv.prefix
##    #;run_prefix   = rv.run_prefix
##    #;directory    = rv.directory
##    #;log_file     = rv.log_file
##    
##    #;comment_file = rv.comment_file
##    
##    #-----------------------------
##    #Print out the model run vars
##    #-----------------------------
##    file_unit.write(line)
##    file_unit.write('  Model Run Variables')
##    file_unit.write(line)
##    Save_Var(unit, 'Prefix:', prefix)
##    Save_Var(unit, 'Run prefix:', run_prefix)
##    Save_Var(unit, 'Directory:', directory)
##    Save_Var(unit, 'Log file:', log_file)
##    Save_Var(unit, 'Comment file:', comment_file)
##    file_unit.write(' ')
##    
##    #------------------------
##    #Print out the grid vars
##    #------------------------
##    file_unit.write(line)
##    file_unit.write('  Grid Variables')
##    file_unit.write(line)
##    Save_Var(unit, 'RTI file:', av.RTI_file)
##    Save_Var(unit, 'Number of columns:', av.ncols, 'LONG')
##    Save_Var(unit, 'Number of rows:', av.nrows, 'LONG')
##    Save_Var(unit, 'Pixel xsize:', av.xres, 'FLOAT')
##    Save_Var(unit, 'Pixel ysize:', av.yres, 'FLOAT')
##    Save_Var(unit, 'DEM data type:', av.data_type)
##    Save_Var(unit, 'Pixel geom code:', av.pixel_geom, 'BYTE')
##    Save_Var(unit, 'Byte order:', av.byte_order)
##    Save_Var(unit, 'Minimum xsize:', av.min_dx, 'FLOAT')
##    Save_Var(unit, 'Basin RTM file:', av.RTM_file)
##    Save_Var(unit, 'Compute precip vol:', av.get_pvolume, 'BYTE')
##    file_unit.write(' ')
##    
##    #------------------------------------
##    #Print out info for monitored pixels
##    #------------------------------------
##    file_unit.write(line)
##    file_unit.write('  Monitored Pixels')
##    file_unit.write(line)
##    w = I2PY_w = where(av.outlet_IDs > np.int32(0))
##    nw = size(I2PY_w[0])
##    Save_Var(unit, 'Number monitored:', nw, 'LONG')
##    #-------------------------------------------------
##    for k in arange(0, nw):
##        kstr = 'Pixel ' + TF_String(k + 1) + ' '
##        wk = w[k]
##        xk = (av.outlet_IDs)[wk] % av.ncols
##        yk = (av.outlet_IDs)[wk] / av.ncols
##        ak = (av.basin_areas)[wk]
##        rk = (av.basin_reliefs)[wk]
##        #-------------------------------------------------
##        Save_Var(unit, (kstr + 'column:'), xk, 'LONG')
##        Save_Var(unit, (kstr + 'row:'), yk, 'LONG')
##        Save_Var(unit, (kstr + 'area:'), ak, 'DOUBLE', '   [km^2]')
##        Save_Var(unit, (kstr + 'relief:'), rk, 'DOUBLE', '  [m]')
##    file_unit.write(' ')
##    
##    #--------------------------
##    #Print out the precip vars
##    #--------------------------
##    file_unit.write(line)
##    file_unit.write('  Precipitation Process Variables')
##    file_unit.write(line)
##    #*** Save_Var, unit, 'Time step [sec]:', pv.dt, 'DOUBLE'
##    if (pv.method == 1):    
##        #---------------------------------------------------------
##        #Method is "Uniform in space", which does not follow the
##        #standard pattern used for all other methods in TopoFlow.
##        #---------------------------------------------------------
##        #Save rates and durations to files, change the type to
##        #"Time series" and the method to "General data types"
##        #--------------------------------------------------------
##        new_rate_type = np.uint8(1)   #(time series)
##        new_dur_type = np.uint8(1)   #(time series)
##        new_rate_file = '00_' + prefix + '_rainrates.txt'
##        new_dur_file = '00_' + prefix + '_durations.txt'
##        #-------------------------------------------------
##        Check_Overwrite(new_rate_file, OK)
##        if logical_not(OK):    
##            msg = array([' ', 'Please use the next dialog to select a filename', 'for saving the time series of RAINRATES.', ' '])
##            result = GUI_Message(msg, INFO=True)
##            new_rate_file = I2PY_filepath = []
##            app = wx.PySimpleApp()
##            I2PY_dialog = wx.FileDialog(parent=None, defaultDir=os.getcwd(), defaultFile=new_rate_file, style=wx.SAVE | wx.OVERWRITE_PROMPT)
##            if (I2PY_dialog.ShowModal() == wx.ID_OK):
##                I2PY_filepath.append(I2PY_dialog.GetPath())
##            I2PY_dialog.Destroy()
##        #----------------------------------------------------------------
##        Check_Overwrite(new_dur_file, OK)
##        if logical_not(OK):    
##            msg = array([' ', 'Please use the next dialog to select a filename', 'for saving the time series of DURATIONS.', ' '])
##            result = GUI_Message(msg, INFO=True)
##            new_dur_file = I2PY_filepath = []
##            app = wx.PySimpleApp()
##            I2PY_dialog = wx.FileDialog(parent=None, defaultDir=os.getcwd(), defaultFile=new_dur_file, style=wx.SAVE | wx.OVERWRITE_PROMPT)
##            if (I2PY_dialog.ShowModal() == wx.ID_OK):
##                I2PY_filepath.append(I2PY_dialog.GetPath())
##            I2PY_dialog.Destroy()
##        #------------------------------------
##        rate_unit = TF_Get_LUN(new_rate_file)
##        file_rate_unit = open(new_rate_file, 'w')
##        I2PY_SWAP_ENDIAN = False
##        #------------------------------------
##        dur_unit = TF_Get_LUN(new_dur_file)
##        file_dur_unit = open(new_dur_file, 'w')
##        I2PY_SWAP_ENDIAN = False
##        #------------------------------------
##        np = idl_func.n_elements(pv.rates)
##        for k in arange(np.int32(0), np):
##            file_rate_unit.write(((pv.rates)[k] * np.float64(3600000)))   #[mm/hr]
##            file_dur_unit.write((pv.durations)[k])
##            #-----------------------------------------------
##            #** printf, rate_unit, 100.0
##            #** printf, rate_unit, (*pv.rates)[k], (*pv.rates)[k] * 3600000d
##            #** printf, rate_unit, (*pv.rates)[k]   ;[mm/hr]
##        file_rate_unit.close()
##        file_dur_unit.close()
##        #--------------------------------------------------------
##        #NB!  Use FACTOR to convert units from [m/s] to [mm/hr]
##        #--------------------------------------------------------
##        Save_Var(unit, 'Method code:', np.uint8(2), 'BYTE')
##        Save_Var2(unit, 'Time step:', 'DOUBLE', np.uint8(0), pv.dt, '', '[sec]', NOT_PTR=True)
##        Save_Var2(unit, 'Rate:', 'DOUBLE', new_rate_type, pv.rates, new_rate_file, '[mm/hr]', FACTOR=np.float64(3600000))
##        Save_Var2(unit, 'Duration:', 'DOUBLE', new_dur_type, pv.durations, new_dur_file, '[min]')
##        #------------------------------------------------------
##        #Save_Var2, unit, 'T_air:', 'DOUBLE', pv.T_air_type, $  ;(now with met vars)
##        #           pv.T_air, pv.T_air_file, '[deg C]'
##    else:    
##        #--------------------------------------------------------
##        #NB!  Use FACTOR to convert units from [m/s] to [mm/hr]
##        #--------------------------------------------------------
##        Save_Var(unit, 'Method code:', pv.method, 'BYTE')
##        Save_Var2(unit, 'Time step:', 'DOUBLE', np.uint8(0), pv.dt, '', '[sec]', NOT_PTR=True)
##        Save_Var2(unit, 'Rate:', 'DOUBLE', pv.rate_type, pv.rates, pv.rate_file, '[mm/hr]', FACTOR=np.float64(3600000))
##        Save_Var2(unit, 'Duration:', 'DOUBLE', pv.duration_type, pv.durations, pv.duration_file, '[min]')
##        #------------------------------------------------------
##        #Save_Var2, unit, 'T_air:', 'DOUBLE', pv.T_air_type, $  ;(now with met vars)
##        #           pv.T_air, pv.T_air_file, '[deg C]'
##    file_unit.write(' ')
##    
##    #---------------------------
##    #Print out the channel vars
##    #---------------------------
##    file_unit.write(line)
##    file_unit.write('  Channel Process Variables')
##    file_unit.write(line)
##    Save_Var(unit, 'Method code:', cv.method, 'BYTE')
##    #-----------------------------------------------------------------
##    Save_Var(unit, 'Kinematic wave flag:', cv.KINEMATIC_WAVE, 'BYTE')
##    Save_Var(unit, 'Diffusive wave flag:', cv.DIFFUSIVE_WAVE, 'BYTE')
##    Save_Var(unit, 'Dynamic wave flag:  ', cv.DYNAMIC_WAVE, 'BYTE')
##    Save_Var(unit, 'Manning flag:       ', cv.MANNING, 'BYTE')
##    Save_Var(unit, 'Law of Wall flag:   ', cv.LAW_OF_WALL, 'BYTE')
##    #-----------------------------------------------------------------
##    Save_Var2(unit, 'Time step:', 'DOUBLE', np.uint8(0), cv.dt, '', '[sec]', NOT_PTR=True)
##    Save_Var2(unit, 'D8 flow code:', 'BYTE', cv.code_type, cv.codes, cv.code_file, '[none]')
##    Save_Var2(unit, 'D8 slope:', 'DOUBLE', cv.slope_type, cv.slopes, cv.slope_file, '[m/m]')
##    #-----------------------------------------------------------------
##    if (cv.MANNING == np.uint8(1)):    
##        Save_Var2(unit, 'Manning N:', 'DOUBLE', cv.nval_type, cv.nvals, cv.nval_file, '[s/m^(1/3)]')
##    #-----------------------------------------------------------------
##    if (cv.LAW_OF_WALL == np.uint8(1)):    
##        Save_Var2(unit, 'Roughness z0:', 'DOUBLE', cv.z0val_type, cv.z0vals, cv.z0val_file, '[m]')
##    #-----------------------------------------------------------------
##    Save_Var2(unit, 'Bed width:', 'DOUBLE', cv.width_type, cv.widths, cv.width_file, '[m]')
##    Save_Var2(unit, 'Bank angle:', 'DOUBLE', cv.angle_type, cv.angles, cv.angle_file, '[deg]')
##    Save_Var2(unit, 'Init. depth:', 'DOUBLE', cv.d0_type, cv.d0, cv.d0_file, '[m]')
##    Save_Var2(unit, 'Sinuosity:', 'DOUBLE', cv.sinu_type, cv.sinu, cv.sinu_file, '[m/m]')
##    #-----------------------------------------------------------------------------
##    Save_Var2(unit, 'Save grid timestep:', 'DOUBLE', np.uint8(0), cv.save_grid_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save Q grids:', cv.SAVE_Q_GRIDS, cv.Q_gs_file, '[m^3/s]')
##    Save_Var3(unit, 'Save u grids:', cv.SAVE_U_GRIDS, cv.u_gs_file, '[m/s]')
##    Save_Var3(unit, 'Save d grids:', cv.SAVE_D_GRIDS, cv.d_gs_file, '[m]')
##    Save_Var3(unit, 'Save f grids:', cv.SAVE_F_GRIDS, cv.f_gs_file, '[none]')
##    #-----------------------------------------------------------------------------
##    Save_Var2(unit, 'Save pixels timestep:', 'DOUBLE', np.uint8(0), cv.save_pixels_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save Q pixels:', cv.SAVE_Q_PIXELS, cv.Q_ts_file, '[m^3/s]')
##    Save_Var3(unit, 'Save u pixels:', cv.SAVE_U_PIXELS, cv.u_ts_file, '[m/s]')
##    Save_Var3(unit, 'Save d pixels:', cv.SAVE_D_PIXELS, cv.d_ts_file, '[m]')
##    Save_Var3(unit, 'Save f pixels:', cv.SAVE_F_PIXELS, cv.f_ts_file, '[none]')
##    #-----------------------------------------------------------------------------
##    file_unit.write(' ')
##    
##    #-----------------------
##    #Print out the met vars
##    #-----------------------
##    file_unit.write(line)
##    file_unit.write('  Meteorological Variables')
##    file_unit.write(line)
##    #----------------------------------------------
##    Save_Var2(unit, 'rho_H2O:', 'DOUBLE', np.uint8(0), mv.rho_H2O, '', '[kg/m^3]')
##    Save_Var2(unit, 'Cp_air:', 'DOUBLE', np.uint8(0), mv.Cp_air, '', '[J/kg/K]')
##    Save_Var2(unit, 'rho_air:', 'DOUBLE', np.uint8(0), mv.rho_air, '', '[kg/m^3]')
##    #--------------------------------------------------------
##    Save_Var2(unit, 'Qn_SW:', 'DOUBLE', mv.Qn_SW_type, mv.Qn_SW, mv.Qn_SW_file, '[W/m^2]')
##    Save_Var2(unit, 'Qn_LW:', 'DOUBLE', mv.Qn_LW_type, mv.Qn_LW, mv.Qn_LW_file, '[W/m^2]')
##    Save_Var2(unit, 'T_air:', 'DOUBLE', mv.T_air_type, mv.T_air, mv.T_air_file, '[deg C]')
##    Save_Var2(unit, 'T_surf:', 'DOUBLE', mv.T_surf_type, mv.T_surf, mv.T_surf_file, '[deg C]')
##    Save_Var2(unit, 'RH:', 'DOUBLE', mv.RH_type, mv.RH, mv.RH_file, '[none]')
##    Save_Var2(unit, 'p0:', 'DOUBLE', mv.p0_type, mv.p0, mv.p0_file, '[mbar]')
##    Save_Var2(unit, 'uz:', 'DOUBLE', mv.uz_type, mv.uz, mv.uz_file, '[m/s]')
##    Save_Var2(unit, 'z:', 'DOUBLE', mv.z_type, mv.z, mv.z_file, '[m]')
##    Save_Var2(unit, 'z0_air:', 'DOUBLE', mv.z0_air_type, mv.z0_air, mv.z0_air_file, '[m]')
##    file_unit.write(' ')
##    
##    #------------------------
##    #Print out the snow vars
##    #------------------------
##    file_unit.write(line)
##    file_unit.write('  Snowmelt Process Variables')
##    file_unit.write(line)
##    Save_Var(unit, 'Method code:', sv.method, 'BYTE')
##    #** Save_Var, unit, 'Time step [sec]:', sv.dt, 'DOUBLE'
##    Save_Var2(unit, 'Time step:', 'DOUBLE', np.uint8(0), sv.dt, '', '[sec]', NOT_PTR=True)
##    #---------------------------------------------------------------
##    Save_Var2(unit, 'Cp_snow:', 'DOUBLE', np.uint8(0), sv.Cp_snow, '', '[J/kg/K]')
##    Save_Var2(unit, 'rho_snow:', 'DOUBLE', sv.rho_snow_type, sv.rho_snow, sv.rho_snow_file, '[kg/m^3]')
##    Save_Var2(unit, 'c0:', 'DOUBLE', sv.c0_type, sv.c0, sv.c0_file, '[mm/day/deg C]')
##    Save_Var2(unit, 'T0:', 'DOUBLE', sv.T0_type, sv.T0, sv.T0_file, '[deg C]')
##    Save_Var2(unit, 'h0_snow:', 'DOUBLE', sv.h0_snow_type, sv.h0_snow, sv.h0_snow_file, '[m]')
##    Save_Var2(unit, 'h0_swe:', 'DOUBLE', sv.h0_swe_type, sv.h0_swe, sv.h0_swe_file, '[m]')
##    #------------------------------------------------------------------------------
##    Save_Var2(unit, 'Save grid timestep:', 'DOUBLE', np.uint8(0), sv.save_grid_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save mr grids:', sv.SAVE_MR_GRIDS, sv.mr_gs_file, '[m/s]')
##    Save_Var3(unit, 'Save hs grids:', sv.SAVE_HS_GRIDS, sv.hs_gs_file, '[m]')
##    Save_Var3(unit, 'Save sw grids:', sv.SAVE_SW_GRIDS, sv.sw_gs_file, '[m]')
##    Save_Var3(unit, 'Save cc grids:', sv.SAVE_CC_GRIDS, sv.cc_gs_file, '[J/m^2]')
##    Save_Var3(unit, 'Save ea grids:', sv.SAVE_EA_GRIDS, sv.ea_gs_file, '[mbar]')
##    Save_Var3(unit, 'Save es grids:', sv.SAVE_ES_GRIDS, sv.es_gs_file, '[mbar]')
##    #------------------------------------------------------------------------------
##    Save_Var2(unit, 'Save pixels timestep:', 'DOUBLE', np.uint8(0), sv.save_pixels_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save mr pixels:', sv.SAVE_MR_PIXELS, sv.mr_ts_file, '[m/s]')
##    Save_Var3(unit, 'Save hs pixels:', sv.SAVE_HS_PIXELS, sv.hs_ts_file, '[m]')
##    Save_Var3(unit, 'Save sw pixels:', sv.SAVE_SW_PIXELS, sv.sw_ts_file, '[m]')
##    Save_Var3(unit, 'Save cc pixels:', sv.SAVE_CC_PIXELS, sv.cc_ts_file, '[J/m^2]')
##    Save_Var3(unit, 'Save ea pixels:', sv.SAVE_EA_PIXELS, sv.ea_ts_file, '[mbar]')
##    Save_Var3(unit, 'Save es pixels:', sv.SAVE_ES_PIXELS, sv.es_ts_file, '[mbar]')
##    file_unit.write(' ')
##    
##    #----------------------
##    #Print out the ET vars
##    #----------------------
##    file_unit.write(line)
##    file_unit.write('  Evapotranspiration Process Variables')
##    file_unit.write(line)
##    Save_Var(unit, 'Method code:', ev.method, 'BYTE')
##    #** Save_Var, unit, 'Time step [sec]:', ev.dt, 'DOUBLE'
##    Save_Var2(unit, 'Time step:', 'DOUBLE', np.uint8(0), ev.dt, '', '[sec]', NOT_PTR=True)
##    #---------------------------------------------------------------
##    Save_Var2(unit, 'alpha:', 'DOUBLE', ev.alpha_type, ev.alpha, ev.alpha_file, '[none]')
##    Save_Var2(unit, 'Ks:', 'DOUBLE', ev.Ks_type, ev.Ks, ev.Ks_file, '[W/m/deg_C]')
##    Save_Var2(unit, 'soil_x:', 'DOUBLE', ev.soil_x_type, ev.soil_x, ev.soil_x_file, '[m]')
##    Save_Var2(unit, 'T_soil_x:', 'DOUBLE', ev.T_soil_x_type, ev.T_soil_x, ev.T_soil_x_file, '[deg C]')
##    #-----------------------------------------------------------------------------
##    Save_Var2(unit, 'Save grid timestep:', 'DOUBLE', np.uint8(0), ev.save_grid_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save er grids:', ev.SAVE_ER_GRIDS, ev.er_gs_file, '[m/s]')
##    #-----------------------------------------------------------------------------
##    Save_Var2(unit, 'Save pixels timestep:', 'DOUBLE', np.uint8(0), ev.save_pixels_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save er pixels:', ev.SAVE_ER_PIXELS, ev.er_ts_file, '[m/s]')
##    file_unit.write(' ')
##    
##    #-------------------------
##    #Print out the infil vars
##    #-------------------------
##    file_unit.write(line)
##    file_unit.write('  Infiltration Process Variables')
##    file_unit.write(line)
##    Save_Var(unit, 'Method code:', iv.method, 'BYTE')
##    Save_Var(unit, 'Number of layers:', iv.n_layers, 'LONG')
##    #** Save_Var, unit, 'Time step [sec]:', iv.dt, 'DOUBLE'
##    #------------------------------------------------------------
##    Save_Var2(unit, 'Time step:', 'DOUBLE', np.uint8(0), iv.dt, '', '[sec]', NOT_PTR=True)
##    #------------------------------------------------------------
##    if (iv.method < 4):    
##        #---------------------------------------
##        #Vars for Green-Ampt and Smith-Parlange
##        #---------------------------------------
##        Save_Var2(unit, 'Ks:', 'DOUBLE', iv.Ks_type[0], iv.Ks, iv.Ks_file[0], '[m/s]')
##        Save_Var2(unit, 'Ki:', 'DOUBLE', iv.Ki_type[0], iv.Ki, iv.Ki_file[0], '[m/s]')
##        Save_Var2(unit, 'qs:', 'DOUBLE', iv.qs_type[0], iv.qs, iv.qs_file[0], '[none]')
##        Save_Var2(unit, 'qi:', 'DOUBLE', iv.qi_type[0], iv.qi, iv.qi_file[0], '[none]')
##        Save_Var2(unit, 'G:', 'DOUBLE', iv.G_type, iv.G, iv.G_file, '[m]')
##        Save_Var2(unit, 'gamma:', 'DOUBLE', iv.gam_type, iv.gam, iv.gam_file, '[none]')
##        Save_Var(unit, 'Closest soil_type:', iv.soil_type[0], 'STRING')
##    else:    
##        #------------------------------------
##        #Vars for Richards' equation method
##        #------------------------------------
##        for k in arange(0, iv.n_layers):
##            str1 = 'Layer ' + TF_String(k + 1) + ' '
##            Save_Var2(unit, str1 + 'Ks:', 'DOUBLE', iv.Ks_type[k], iv.Ks_val[k], iv.Ks_file[k], '[m/s]')
##            Save_Var2(unit, str1 + 'Ki:', 'DOUBLE', iv.Ki_type[k], iv.Ki_val[k], iv.Ki_file[k], '[m/s]')
##            Save_Var2(unit, str1 + 'qs:', 'DOUBLE', iv.qs_type[k], iv.qs_val[k], iv.qs_file[k], '[unitless]')
##            Save_Var2(unit, str1 + 'qi:', 'DOUBLE', iv.qi_type[k], iv.qi_val[k], iv.qi_file[k], '[unitless]')
##            Save_Var2(unit, str1 + 'qr:', 'DOUBLE', iv.qr_type[k], iv.qr_val[k], iv.qr_file[k], '[unitless]')
##            Save_Var2(unit, str1 + 'pB:', 'DOUBLE', iv.pB_type[k], iv.pB_val[k], iv.pB_file[k], '[m]')
##            Save_Var2(unit, str1 + 'pA:', 'DOUBLE', iv.pA_type[k], iv.pA_val[k], iv.pA_file[k], '[m]')
##            Save_Var2(unit, str1 + 'lam:', 'DOUBLE', iv.lam_type[k], iv.lam_val[k], iv.lam_file[k], '[unitless]')
##            Save_Var2(unit, str1 + 'c:', 'DOUBLE', iv.c_type[k], iv.c_val[k], iv.c_file[k], '[unitless]')
##            Save_Var2(unit, str1 + 'dz:', 'DOUBLE', np.uint8(0), iv.dz_val[k], '', '[m]', NOT_PTR=True)
##            Save_Var(unit, str1 + 'nz:', iv.nz_val[k], 'LONG')
##            Save_Var(unit, str1 + 'soil_type:', iv.soil_type[k], 'STRING')
##    #-----------------------------------------------------------------------------
##    Save_Var2(unit, 'Save grid timestep:', 'DOUBLE', np.uint8(0), iv.save_grid_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save v0 grids:', iv.SAVE_V0_GRIDS, iv.v0_gs_file, '[m/s]')
##    Save_Var3(unit, 'Save q0 grids:', iv.SAVE_Q0_GRIDS, iv.q0_gs_file, '[none]')
##    Save_Var3(unit, 'Save I  grids:', iv.SAVE_I_GRIDS, iv.I_gs_file, '[m]')
##    Save_Var3(unit, 'Save Zw grids:', iv.SAVE_ZW_GRIDS, iv.Zw_gs_file, '[m]')
##    #-----------------------------------------------------------------------------
##    Save_Var2(unit, 'Save pixels timestep:', 'DOUBLE', np.uint8(0), iv.save_pixels_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save v0 pixels:', iv.SAVE_V0_PIXELS, iv.v0_ts_file, '[m/s]')
##    Save_Var3(unit, 'Save q0 pixels:', iv.SAVE_Q0_PIXELS, iv.q0_ts_file, '[none]')
##    Save_Var3(unit, 'Save I  pixels:', iv.SAVE_I_PIXELS, iv.I_ts_file, '[m]')
##    Save_Var3(unit, 'Save Zw pixels:', iv.SAVE_ZW_PIXELS, iv.Zw_ts_file, '[m]')
##    #-----------------------------------------------------------------------------
##    Save_Var2(unit, 'Save stack timestep:', 'DOUBLE', np.uint8(0), iv.save_stack_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save q stacks:', iv.SAVE_Q_STACKS, iv.q_stack_file, '[none]')
##    Save_Var3(unit, 'Save p stacks:', iv.SAVE_P_STACKS, iv.p_stack_file, '[m]')
##    Save_Var3(unit, 'Save K stacks:', iv.SAVE_K_STACKS, iv.K_stack_file, '[m/s]')
##    Save_Var3(unit, 'Save v stacks:', iv.SAVE_V_STACKS, iv.v_stack_file, '[m/s]')
##    #-----------------------------------------------------------------------------
##    Save_Var2(unit, 'Save profile timestep:', 'DOUBLE', np.uint8(0), iv.save_profile_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save q profiles:', iv.SAVE_Q_PROFILES, iv.q_profile_file, '[none]')
##    Save_Var3(unit, 'Save p profiles:', iv.SAVE_P_PROFILES, iv.p_profile_file, '[m]')
##    Save_Var3(unit, 'Save K profiles:', iv.SAVE_K_PROFILES, iv.K_profile_file, '[m/s]')
##    Save_Var3(unit, 'Save v profiles:', iv.SAVE_V_PROFILES, iv.v_profile_file, '[m/s]')
##    #-----------------------------------------------------------------------------
##    file_unit.write(' ')
##    
##    
##    #----------------------
##    #Print out the GW vars
##    #----------------------
##    file_unit.write(line)
##    file_unit.write('  Subsurface Flow Process Variables')
##    file_unit.write(line)
##    Save_Var(unit, 'Method code:', gv.method, 'BYTE')
##    Save_Var(unit, 'Number of layers:', gv.nlayers, 'LONG')
##    #** Save_Var, unit, 'Time step [sec]:', gv.dt, 'DOUBLE'
##    #---------------------------------------------------------------
##    Save_Var2(unit, 'Time step:', 'DOUBLE', np.uint8(0), gv.dt, '', '[sec]', NOT_PTR=True)
##    Save_Var2(unit, 'z_surface:', 'DOUBLE', gv.elev_type, gv.elev, gv.elev_file, '[m]')
##    Save_Var2(unit, 'z0_table:', 'DOUBLE', gv.h0_table_type, gv.h0_table, gv.h0_table_file, '[m]')
##    Save_Var2(unit, 'd_freeze:', 'DOUBLE', gv.d_freeze_type, gv.d_freeze, '', '[m]')
##    Save_Var2(unit, 'd_thaw:', 'DOUBLE', gv.d_thaw_type, gv.d_thaw, '', '[m]')
##    #---------------------------------------------------------------
##    for k in arange(0, gv.nlayers):
##        str1 = 'Layer ' + TF_String(k + 1) + ' '
##        Save_Var2(unit, str1 + 'Ks:', 'DOUBLE', gv.Ks_type[k], gv.Ks[k], gv.Ks_file[k], '[m/s]')
##        Save_Var2(unit, str1 + 'qs:', 'DOUBLE', gv.qs_type[k], gv.qs[k], gv.qs_file[k], '[none]')
##        Save_Var2(unit, str1 + 'thickness:', 'DOUBLE', gv.th_type[k], gv.th[k], gv.th_file[k], '[m]')
##    #-----------------------------------------------------------------------------
##    Save_Var2(unit, 'Save grid timestep:', 'DOUBLE', np.uint8(0), gv.save_grid_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save ht grids:', gv.SAVE_HT_GRIDS, gv.ht_gs_file, '[m]')
##    Save_Var3(unit, 'Save df grids:', gv.SAVE_DF_GRIDS, gv.df_gs_file, '[m]')
##    Save_Var3(unit, 'Save dt grids:', gv.SAVE_DT_GRIDS, gv.dt_gs_file, '[m]')
##    #-----------------------------------------------------------------------------
##    Save_Var2(unit, 'Save pixels timestep:', 'DOUBLE', np.uint8(0), gv.save_pixels_dt, '', '[sec]', NOT_PTR=True)
##    Save_Var3(unit, 'Save ht pixels:', gv.SAVE_HT_PIXELS, gv.ht_ts_file, '[m]')
##    Save_Var3(unit, 'Save df pixels:', gv.SAVE_DF_PIXELS, gv.df_ts_file, '[m]')
##    Save_Var3(unit, 'Save dt pixels:', gv.SAVE_DT_PIXELS, gv.dt_ts_file, '[m]')
##    file_unit.write(' ')
##    
##    #-----------------------------
##    #Print out the diversion vars
##    #-----------------------------
##    file_unit.write(line)
##    file_unit.write('  Diversion Process Variables')
##    file_unit.write(line)
##    Save_Var(unit, 'Method code:', dv.method, 'BYTE')
##    #-----------------------------------------------------------
##    Save_Var3(unit, 'Use sources:', dv.use_sources, dv.source_file, '[N/A]')
##    Save_Var3(unit, 'Use sinks:', dv.use_sinks, dv.sink_file, '[N/A]')
##    Save_Var3(unit, 'Use canals:', dv.use_canals, dv.canal_file, '[N/A]')
##    file_unit.write(' ')
##    
##    #------------------------
##    #Print out the stop vars
##    #------------------------
##    file_unit.write(line)
##    file_unit.write('  Stopping Criterion Variables')
##    file_unit.write(line)
##    Save_Var(unit, 'Method code:', fv.method, 'BYTE')
##    #--------------------------------------------------------------------
##    Save_Var(unit, 'Qp_fraction:', fv.Qp_fraction, 'DOUBLE')
##    Save_Var(unit, 'T_stop_model:', fv.T_stop_model, 'DOUBLE', '   [min]')
##    Save_Var(unit, 'T_stop_real:', fv.T_stop_real, 'DOUBLE', '   [min]')
##    Save_Var(unit, 'n_steps:', fv.n_steps, 'LONG')
##    file_unit.write(' ')
##    
##    #----------------------
##    #Close the output file
##    #----------------------
##    file_unit.close()
##    
##    #---------------------------
##    #Display information dialog
##    #---------------------------
##    msg = array(['All input variables saved to the file:', out_file, ' '])
##    result = GUI_Message(msg, INFO=True, TITLE='Save Settings Status')
##    
###   Save_All_TF_Vars
#******************************************************************
def Read_Vars(file_unit, var2=None, var3=None, var4=None, \
              data_type='STRING', name=None):

    #-------------------------------------------------------------
    # Notes: By default, the STRSPLIT routine uses both blank
    #        spaces and tabs as delimiters.  Others can be
    #        specified with the optional "pattern" argument.

    # NB!   If (type eq 'STRING') (the default), then we want
    #        to read everything after the ":", which may contain
    #        space characters in the interior (e.g a directory),
    #        but with leading/trailing spaces removed.
    #-------------------------------------------------------------
    n_params = 5 - [var2, var3, var4].count(None)
    var1 = None
    _opt = (var2, var3, var4)
    def _ret():
        _opt_rv = list(zip(_opt, [var2, var3, var4]))
        _rv = [var1]
        _rv += [_o[1] for _o in _opt_rv if _o[0] is not None]
        if (len(_rv) == 1):
            return _rv[0]
        else:
            return tuple(_rv)
        
    line = ''
    line = idl_func.readf(file_unit, line)
    
    #-----------------------------------
    # Extract the variable name, which
    # may contain blank spaces
    #-----------------------------------
    p = line.find(':')
    if (p == -1):    
        return _ret()
    name = line[0:0+p]
    line = line[p + 1:]
    
    #-------------------------------
    # Extract variables as strings
    #-------------------------------
    words = line.split()
    count = len(words)

##    print 'words =', words
##    print 'count =', count
    
    #################################    
    try:
        value = eval(words[0])
    except:
        value = words[0]

    #--------------
    # For testing
    #--------------
##    print 'VALUE =', value
##    print 'data_type =', data_type
##    print '-------------------------------------------------'
        
    #------------------------
    # Convert var1 string ?
    #------------------------
    if (count >= 1):    
        dtype = data_type.upper()

        if   (dtype == 'BYTE'):
            var1 = np.int16(value)
            ## var1 = np.uint8(np.int16(value))
        elif (dtype  == 'INTEGER'):    
            var1 = np.int16(value)
        elif (dtype  == 'LONG'):    
            var1 = np.int32(value)
        elif (dtype  == 'FLOAT'):    
            var1 = np.float32(value)
        elif (dtype  == 'DOUBLE'):    
            var1 = np.float64(value)
        elif (dtype  == 'STRING'):    
            var1 = value
        elif (dtype  == 'FILE'):    
            var1 = line.strip()     #(may contain blanks)
        else:    
            var1 = value
 
    #-------------------------------------
    # Extract additional vars as strings
    #-------------------------------------
    if (count >= 2):    
        var2 = words[1]
    if (count >= 3):    
        var3 = words[2]
    if (count >= 4):    
        var4 = words[3]
        
    return _ret()
 
#   Read_Vars
#******************************************************************
def Load_Var(value_str, var_type, factor=np.float64(1)):

    var = None  #####
    
    if (var_type.upper() == 'SCALAR'):    
        if (factor != 1):    
            var = (np.float64(value_str) * factor)
        else:    
            var = np.float64(value_str)
        filename = ''
    else:
        var = None  #######
        filename = value_str
    
    return (var, filename)

#   Load_Var
#******************************************************************
def Type_Code(var_type):

    cmap = {'SCALAR':        0, \
            'TIME_SERIES':   1, \
            'TIME SERIES':   1, \
            'GRID'       :   2, \
            'GRID_STACK' :   3, \
            'GRID STACK' :   3, \
            'GRID_SEQUENCE': 3, \
            'GRID SEQUENCE': 3 }

    return cmap[var_type.upper()]
    
#   Type_Code
#******************************************************************
