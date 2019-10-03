#!/usr/bin/env python

#  S.D. Peckham
#  Nov. 8, 2013. Moved into TopoFlow package; experimental)
#  May 1, 2009.
#  April 28-30, 2009.
#  April 21-24, 2009.

#---------------------------------------------------------------------------
# Notes:  The files in this directory provide a partial reproduction
#         of the TopoFlow IDL GUI, using wxPython.  All of the wizard
#         panels and input/output variable panels are implemented but
#         don't do anything yet.  Additional dialogs in the Create
#         menu, etc. are not yet implemented.
#
#         To launch the GUI, type:
#            % python Main_Dialog.py
#         It must be launched this way in order for the Exit button
#         and File > Close to work.
#
#         Note that the bitmaps, dialogs and help folders in this
#         directory contain the required files, and this version
#         looks for them in /Applications/TopoFlow.  However, they
#         should be stored and found in the gui subpackage instead.
#
#---------------------------------------------------------------------------

import wx
import wx.grid
import webbrowser  ##  standard Python module
import sys
import time

from . import GUI_utils
from .Input_Dialog  import *
from .Output_Dialog import *

## import wx.html
## import xml.dom.minidom

#---------------------------------------------------------------------------

## Note:  widget sensitivity can be set with Disable() and Enable().

#--------------------------------------------------------------------------- 
##def launch():   # (11/8/13)
##
##    #---------------------------------------------------------------
##    # With this approach, Exit button and File > Close don't work.
##    #---------------------------------------------------------------
##    app = wx.App()
##    frame  = TF_Main_Dialog()
##    frame.Show()
##    app.MainLoop()
##
#---------------------------------------------------------------------------    
class TF_Main_Dialog(wx.Frame):

    #-------------------------------------------------------
    # Notes:  This class is for creating the main TopoFlow
    #         wizard dialog.
    
    #         Default is wx.ALIGN_LEFT for labels & text
    #         Can add keyboard shortcuts (accelerator keys
    #         with "&File", for example.
    #-----------------------------------------------------
        
    #---------------------------
    #  Create top-level dialog
    #---------------------------
    def __init__(self, parent=None, ID=-1):
   
        #---------------------------
        # Set resource directories
        #---------------------------
        if (sys.platform == 'win32'):
            self.home_dir = "C:\Program Files\TopoFlow\\"
        elif (sys.platform == 'darwin'):
            self.home_dir = '/Applications/TopoFlow/'
        else:
            self.home_dir = '/usr/local/topoflow/'
        #----------------------------------------------------
        self.help_dir       = self.home_dir + 'help/'
        self.bmp_dir        = self.home_dir + 'bitmaps/'
        self.dialog_in_dir  = self.home_dir + 'dialogs/in/'
        self.dialog_out_dir = self.home_dir + 'dialogs/out/'
        
##        self.dialog_in_dir  = self.home_dir + 'xml/in/'
##        self.dialog_out_dir = self.home_dir + 'xml/out/'
        
        #---------------------------------
        # Embed process module objects ??
        #---------------------------------
        self.p_method = 1
        self.s_method = 0
        self.e_method = 0
        self.i_method = 0
        self.g_method = 0
        self.c_method = 1
        self.d_method = 0
        
        #----------------------------------
        # Used for spacing in some panels
        #----------------------------------
        self.vgap = 8  # (vertical)
        self.hgap = 5  # (horizontal)

        #-------------------------------------------
        # Create the "main frame" for the TopoFlow
        # wizard with a hard-wired size, for now
        #-------------------------------------------
        self.title = "TopoFlow 2.0 beta Dialog (4/23/09)"
        wx.Frame.__init__(self, parent, ID, self.title, \
                          pos=(100,100), size=(640,480))

        #-------------------------------
        # Add TopoFlow menu bar at top
        #-------------------------------
        self.Add_Menu_Bar()

        #----------------------------------------------------
        # Create one panel and sizer that will hold the
        # different wizard panels as sizers (for Hide/Show)
        #----------------------------------------------------
        panel = wx.Panel(self)
        panel.SetBackgroundColour('Light Blue')
        self.main_panel = panel
        #-----------------------------------------
        sizer = wx.BoxSizer(wx.VERTICAL)        
        panel.SetSizer(sizer)       ########
        self.main_sizer = sizer
        
        #-------------------------------------
        # Create the wizard panels, but hide
        # all but the navigation panel
        #-------------------------------------
        self.Add_Navigation_Panel()
        self.Add_Run_Info_Panel()  
        self.Add_Methods_Panel()
        self.Add_Basin_Info_Panel()
        self.Add_Output_Log_Panel()
        self.Add_Preprocessing_Panel()
        self.Add_Plotting_Panel()
        #----------------------------------------
        sizer.Show(self.navigation_sizer)  ####
        self.current_panel = 'navigation'
        
        # self.Plotting_Test()  ### (doesn't work here ?)
        
        #-----------------------------------------------
        # Auto-fitting doesn't work because the hidden
        # panels make the y-size too big
        #-----------------------------------------------
        # sizer.Fit(self)
        # self.Show(True)  # (the frame, not a panel)
        
    #   __init__()
    #----------------------------------------------------------------
    def Add_Menu_Bar(self):

        #--------------------------------
        # Create menus for the menu bar
        #--------------------------------        
        File_menu = wx.Menu()
        open_item  = File_menu.Append(-1, "Open Data Set", "help")
        load_item  = File_menu.Append(-1, "Load Input Vars", "help")
        save_item  = File_menu.Append(-1, "Save Input Vars", "help")
        close_item = File_menu.Append(-1, "Close Any Open File", "help")
        File_menu.AppendSeparator()
        exit_item  = File_menu.Append(-1, "Exit TopoFlow", "help")
        self.Bind(wx.EVT_MENU, self.On_File_Open,  open_item)
        self.Bind(wx.EVT_MENU, self.On_File_Load,  load_item)
        self.Bind(wx.EVT_MENU, self.On_File_Save,  save_item)
        self.Bind(wx.EVT_MENU, self.On_File_Close, close_item)
        self.Bind(wx.EVT_MENU, self.On_File_Exit,  exit_item)
        #-----------------------------------------------------------------
        Goto_menu = wx.Menu()
        nav_item  = Goto_menu.Append(-1, "Main Level", "help")
        pre_item  = Goto_menu.Append(-1, "Preprocessing", "help")
        olog_item = Goto_menu.Append(-1, "Output Log", "help")
        plot_item = Goto_menu.Append(-1, "Plot Results", "help")
        self.Bind(wx.EVT_MENU, self.On_Goto_Navigation,    nav_item)
        self.Bind(wx.EVT_MENU, self.On_Goto_Preprocessing, pre_item)
        self.Bind(wx.EVT_MENU, self.On_Goto_Output_Log,    olog_item)
        self.Bind(wx.EVT_MENU, self.On_Goto_Plotting,      plot_item)
        #-----------------------------------------------------------------
        Create_menu = wx.Menu()
        create1 = Create_menu.Append(-1, "Channel Geometry Grids", "help")
        create2 = Create_menu.Append(-1, "Profile-smoothed DEM",   "help")
        create3 = Create_menu.Append(-1, "RTG file for Initial Depth", "help")
        create4 = Create_menu.Append(-1, "RTS file for Station Data",  "help")
        create5 = Create_menu.Append(-1, "RTS file for Qnet Shortwave Flux", "help")
        create6 = Create_menu.Append(-1, "RTS file for Qnet Longwave Flux", "help")
        create7 = Create_menu.Append(-1, "RTS file for Fractal Rain", "help")
        #-----------------------------------------------------------------
        Plot_menu = wx.Menu()
        plot1 = Plot_menu.Append(-1, "Function",   "help")
        plot2 = Plot_menu.Append(-1, "RTS File",   "help")
        plot3 = Plot_menu.Append(-1, "RTS to MPG", "help")        
        #-----------------------------------------------------------------
        Help_menu = wx.Menu()
        about_item = Help_menu.Append(-1, "About TopoFlow",    "help")
        lic_item   = Help_menu.Append(-1, "License Agreement", "help")
        Help_menu.AppendSeparator()
        new_item   = Help_menu.Append(-1, "What's New in 2.0", "help")
        tut_item   = Help_menu.Append(-1, "Short Tutorial", "help")
        self.Bind(wx.EVT_MENU, self.On_Help_About,     about_item)
        self.Bind(wx.EVT_MENU, self.On_Help_License,   lic_item)
        self.Bind(wx.EVT_MENU, self.On_Help_Whats_New, new_item)
        self.Bind(wx.EVT_MENU, self.On_Help_Tutorial,  tut_item)
        
        #------------------------------------
        # Create menu bar and add the menus
        #------------------------------------
        self.CreateStatusBar()
        TF_menu_bar = wx.MenuBar()
        ## TF_menu_bar = wx.MenuBar(wx.MB_DOCKABLE)
        TF_menu_bar.Append(File_menu,   "File")
        TF_menu_bar.Append(Goto_menu,   "Goto")
        TF_menu_bar.Append(Create_menu, "Create")
        TF_menu_bar.Append(Plot_menu,   "Plot")
        TF_menu_bar.Append(Help_menu,   "Help")
        self.SetMenuBar(TF_menu_bar)

    #   Add_Menu_Bar()       
    #----------------------------------------------------------------
    def Add_Navigation_Panel(self):

        #----------------------------------
        # Create a new panel with a sizer
        #----------------------------------
        panel = self.main_panel
        hgap  = 20
        vgap  = 15
        sizer = wx.GridSizer(2, 2, hgap, vgap)
        self.navigation_sizer = sizer
        
        #-----------------------------------
        # Can we find the button bitmaps ?
        #-----------------------------------
        
        #-------------------------
        # Get the button bitmaps
        #-------------------------
        bmp_nx    = 168
        bmp_ny    = 128
        bmp_file1 = self.bmp_dir + "TF_Button1.bmp"
        bmp_file2 = self.bmp_dir + "TF_Button2.bmp"
        bmp_file3 = self.bmp_dir + "TF_Button3.bmp"
        bmp_file4 = self.bmp_dir + "TF_Button4.bmp"
        bmp1 = wx.Image(bmp_file1, wx.BITMAP_TYPE_BMP).ConvertToBitmap()
        bmp2 = wx.Image(bmp_file2, wx.BITMAP_TYPE_BMP).ConvertToBitmap()
        bmp3 = wx.Image(bmp_file3, wx.BITMAP_TYPE_BMP).ConvertToBitmap()
        bmp4 = wx.Image(bmp_file4, wx.BITMAP_TYPE_BMP).ConvertToBitmap()
        #-----------------------------------------------------------------
        new_button  = wx.BitmapButton(panel, -1, bmp1)
        pre_button  = wx.BitmapButton(panel, -1, bmp2)
        plot_button = wx.BitmapButton(panel, -1, bmp3)
        exit_button = wx.BitmapButton(panel, -1, bmp4)
        #-----------------------------------------------------------------
        proportion = 0
        # flag   = wx.ALIGN_CENTER  # (affects resizing by user)
        flag   = wx.ALL
        border = 20
        sizer.Add(new_button,  proportion, flag, border)
        sizer.Add(pre_button,  proportion, flag, border)
        sizer.Add(plot_button, proportion, flag, border)
        sizer.Add(exit_button, proportion, flag, border)
        #------------------------------------------------------------------
        self.Bind(wx.EVT_BUTTON, self.On_Goto_Run_Info,      new_button)
        self.Bind(wx.EVT_BUTTON, self.On_Goto_Preprocessing, pre_button)
        self.Bind(wx.EVT_BUTTON, self.On_Goto_Plotting,      plot_button)
        self.Bind(wx.EVT_BUTTON, self.On_File_Exit,          exit_button)

        #---------------------------------------
        # Add run_info_sizer to the main_sizer
        #---------------------------------------
        self.main_sizer.Add(sizer, 0, wx.ALL, self.vgap)
        self.main_sizer.Hide(sizer)
        
    #   Add_Navigation_Panel()
    #----------------------------------------------------------------
    def Add_Run_Info_Panel(self):

        #-----------------------------------------------------
        def Browse_Run_Directory(event):
            directory = GUI_utils.Get_Directory()
            if (directory is not None):
                self.run_directory = directory
                self.run_directory_box.SetValue(directory)
        #-----------------------------------------------------
        def Apply_Run_Prefix(event):
            run_prefix   = self.run_prefix_box.GetValue()
            log_file     = (run_prefix + '_LOG.txt')
            comment_file = (run_prefix + '_README.txt')
            self.log_file_box.SetValue(log_file)
            self.comment_file_box.SetValue(comment_file)
        #-----------------------------------------------------
        self.run_directory  = self.home_dir
        self.data_prefix    = 'Treynor'
        self.run_prefix     = 'Case1'
        self.log_file       = 'Case1_LOG.txt'
        self.comment_file   = 'Case1_README.txt'
        self.comments       = 'None'
        self.basin_RTM_file = self.data_prefix + '_basin.rtm'
        ########################################################
    
        #---------------------------------------------------
        # Create a "static box" and associated sizer for
        # everything but the bottom buttons.  Static boxes
        # provide a frame and a title for grouping.
        #---------------------------------------------------
        panel = self.main_panel
        box   = wx.StaticBox(panel, -1, "Model run information:")
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        self.run_info_sizer = sizer
        
        #---------------------------------------------
        #  Create another sizer box for rows of info
        #---------------------------------------------
        fg_sizer = wx.FlexGridSizer(cols=3, hgap=self.hgap, vgap=0)

        #------------------------
        # Labels and text boxes
        #------------------------
        self.run_directory_box = \
            GUI_utils.Labeled_Text_Box(self, panel, fg_sizer, \
                                  "Model run directory:", \
                                  self.run_directory,
                                  button="Browse...", \
                                  function=Browse_Run_Directory)
        self.data_prefix_box = \
            GUI_utils.Labeled_Text_Box(self, panel, fg_sizer, \
                                  "Data set prefix:", \
                                  self.data_prefix)
        self.run_prefix_box = \
            GUI_utils.Labeled_Text_Box(self, panel, fg_sizer, \
                                  "Model run prefix:", \
                                  self.run_prefix, \
                                  button="Apply", \
                                  function=Apply_Run_Prefix)
        self.log_file_box = \
            GUI_utils.Labeled_Text_Box(self, panel, fg_sizer, \
                                  "Log file name:", \
                                  self.log_file)
        self.comment_file_box = \
            GUI_utils.Labeled_Text_Box(self, panel, fg_sizer, \
                                  "Comment file name:", \
                                  self.comment_file)
        self.comments_box = \
            GUI_utils.Labeled_Text_Box(self, panel, fg_sizer, \
                                  "Run comments:", \
                                  self.comments, MULTI=True)
        
        #-------------------------------------
        # Add fg_sizer to the run_info_sizer
        #------------------------------------------------
        # Next 2 lines adds padding inside of StaticBox
        #------------------------------------------------
        # pad_row = wx.StaticText(panel, -1, " ")
        # sizer.Add(pad_row,  0, wx.ALL, self.vgap)
        sizer.Add(fg_sizer, 0, wx.ALL, self.vgap)
    
        #----------------------------------------------
        # Add Back and Next buttons to run_info_sizer
        #----------------------------------------------
        back_fcn = self.On_Goto_Navigation
        next_fcn = self.On_Goto_Methods
        GUI_utils.Add_Back_Next_Buttons(self, panel, sizer, \
                                        back_fcn, next_fcn)

        #---------------------------------------
        # Add run_info_sizer to the main_sizer
        # and hide it initially
        #---------------------------------------
        self.main_sizer.Add(sizer, 0, wx.ALL, self.vgap)
        self.main_sizer.Hide(sizer)
        
    #   Add_Run_Info_Panel() 
    #----------------------------------------------------------------
    def Add_Methods_Panel(self):

        #  Could move down and call "On_Precip_Choice", etc.
        ######################################################
        #-----------------------------------------------------------
        def Set_Precip_Method(event):
            self.p_method = self.precip_droplist.GetSelection()
        #-----------------------------------------------------------
        def Set_Snowmelt_Method(event):
            self.s_method = self.snowmelt_droplist.GetSelection()
        #-----------------------------------------------------------
        def Set_ET_Method(event):
            self.e_method = self.ET_droplist.GetSelection()
        #-----------------------------------------------------------
        def Set_Infil_Method(event):
            self.i_method = self.infil_droplist.GetSelection()
        #-----------------------------------------------------------
        def Set_GW_Method(event):
            self.g_method = self.GW_droplist.GetSelection()
        #-----------------------------------------------------------
        def Set_Channel_Method(event):
            self.c_method = self.channel_droplist.GetSelection()
        #-----------------------------------------------------------
        def Set_Diversion_Method(event):
            self.d_method = self.diversion_droplist.GetSelection()
            
        #---------------------------------------------------
        # Create a "static box" and associated sizer for
        # everything but the bottom buttons.  Static boxes
        # provide a frame and a title for grouping.
        #---------------------------------------------------
        panel = self.main_panel
        box   = wx.StaticBox(panel, -1, "Modeling methods:")
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        self.methods_sizer = sizer

        #---------------------------------------------
        #  Create another sizer box for rows of info
        #---------------------------------------------
        fg_sizer = wx.FlexGridSizer(cols=4, hgap=self.hgap, vgap=0)

        #-------------------------
        # Add the column headers
        #-------------------------
        process_header = wx.StaticText(panel, -1, "Physical process:")
        method_header  = wx.StaticText(panel, -1, "Method to model process:")
        null_header1   = wx.StaticText(panel, -1, "")
        null_header2   = wx.StaticText(panel, -1, "")  # (must be separate)
        #---------------------------------------------------------------------
##        font_size   = 16
##        font_family = wx.DEFAULT
##        font_style  = wx.BOLD | wx.ITALIC
##        font_weight = wx.NORMAL
##        font = wx.Font(font_size, font_family, font_style, font_weight)
##        process_header.SetFont(font)
##        method_header.SetFont(font)
        #---------------------------------------------------------------------        
        fg_sizer.Add(process_header, 0, wx.ALL, self.hgap)
        fg_sizer.Add(method_header,  0, wx.ALL, self.hgap)
        fg_sizer.Add(null_header1,   0, wx.ALL, self.hgap)
        fg_sizer.Add(null_header2,   0, wx.ALL, self.hgap)
        #---------------------------------------------------------------------
        nh1   = wx.StaticText(panel, -1, "")
        nh2   = wx.StaticText(panel, -1, "")
        nh3   = wx.StaticText(panel, -1, "")
        nh4   = wx.StaticText(panel, -1, "")
        fg_sizer.Add(nh1, 0, wx.ALL, self.hgap)
        fg_sizer.Add(nh2, 0, wx.ALL, self.hgap)
        fg_sizer.Add(nh3, 0, wx.ALL, self.hgap)
        fg_sizer.Add(nh4, 0, wx.ALL, self.hgap)
        
        #-----------------------
        # Labels and droplists
        #-----------------------
        self.precip_droplist = \
            GUI_utils.Labeled_Droplist(self, panel, fg_sizer, \
                                  "Precipitation:", \
                                  GUI_utils.Precip_Methods(),
                                  Set_Precip_Method, \
                                  button1="In...", \
                                  b1_function=self.On_Precip_In_Button, \
                                  choice_number=1)
        self.snowmelt_droplist = \
            GUI_utils.Labeled_Droplist(self, panel, fg_sizer, \
                                  "Snowmelt:", \
                                  GUI_utils.Snowmelt_Methods(), \
                                  Set_Snowmelt_Method, \
                                  button1="In...", \
                                  b1_function=self.On_Snowmelt_In_Button, \
                                  button2="Out...", \
                                  b2_function=self.On_Snowmelt_Out_Button)
        self.ET_droplist = \
            GUI_utils.Labeled_Droplist(self, panel, fg_sizer, \
                                  "Evapotranspiration:", \
                                  GUI_utils.ET_Methods(),
                                  Set_ET_Method, \
                                  button1="In...", \
                                  b1_function=self.On_ET_In_Button, \
                                  button2="Out...", \
                                  b2_function=self.On_ET_Out_Button)
        self.infil_droplist = \
            GUI_utils.Labeled_Droplist(self, panel, fg_sizer, \
                                  "Infiltration:", \
                                  GUI_utils.Infil_Methods(), \
                                  Set_Infil_Method, \
                                  button1="In...", \
                                  b1_function=self.On_Infil_In_Button, \
                                  button2="Out...", \
                                  b2_function=self.On_Infil_Out_Button)
        self.GW_droplist = \
            GUI_utils.Labeled_Droplist(self, panel, fg_sizer, \
                                  "Subsurface flow:", \
                                  GUI_utils.GW_Methods(),
                                  Set_GW_Method, \
                                  button1="In...", \
                                  b1_function=self.On_GW_In_Button, \
                                  button2="Out...", \
                                  b2_function=self.On_GW_Out_Button)
        self.channel_droplist = \
            GUI_utils.Labeled_Droplist(self, panel, fg_sizer, \
                                  "Channel flow:", \
                                  GUI_utils.Channel_Methods(), \
                                  Set_Channel_Method, \
                                  button1="In...", \
                                  b1_function=self.On_Channel_In_Button, \
                                  button2="Out...", \
                                  b2_function=self.On_Channel_Out_Button, \
                                  choice_number=1)
        self.diversion_droplist = \
            GUI_utils.Labeled_Droplist(self, panel, fg_sizer, \
                                  "Diversions:", \
                                  GUI_utils.Diversion_Methods(), \
                                  Set_Diversion_Method, \
                                  button1="In...", \
                                  b1_function=self.On_Diversion_In_Button)
        
        #------------------------------------
        # Add fg_sizer to the methods_sizer
        #------------------------------------
        sizer.Add(fg_sizer, 0, wx.ALL, self.vgap)

        #---------------------------------------------
        # Add Back and Next buttons to methods_sizer
        #---------------------------------------------
        back_fcn = self.On_Goto_Run_Info
        next_fcn = self.On_Goto_Basin_Info
        GUI_utils.Add_Back_Next_Buttons(self, panel, sizer, \
                                        back_fcn, next_fcn)
        
       
        #---------------------------------------
        # Add run_info_sizer to the main_sizer
        # and then hide it
        #---------------------------------------
        self.main_sizer.Add(sizer, 0, wx.ALL, self.vgap)
        self.main_sizer.Hide(sizer)
        
    #   Add_Methods_Panel()
    #----------------------------------------------------------------
    def Add_Basin_Info_Panel(self):

        # Move these down and rename to On_Read_Basin_Info, etc. ?
        #-----------------------------------------------------------
        def Read_Table(event):
            print('Not ready to read table yet.')
        #-----------------------------------------------------------
        def Save_Table(event):
            print('Not ready to save table yet.')
        #-----------------------------------------------------------
        def Check_Mass_Balance(event):
            self.check_basin_mass_balance = \
                      not(self.check_basin_mass_balance)
            # print 'choice =', self.check_basin_mass_balance
            ## self.basin_info_sizer.Hide(rtm_sizer)
            ## change sensitivity of text box
        #-----------------------------------------------------------
        self.check_basin_mass_balance = False   ####
        
        #---------------------------------------------------
        # Create a "static box" and associated sizer for
        # everything but the bottom buttons.  Static boxes
        # provide a frame and a title for grouping.
        #---------------------------------------------------
        panel = self.main_panel
        box   = wx.StaticBox(panel, -1, "Info for monitored basins:")
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        self.basin_info_sizer = sizer

        #---------------------------------------------------
        # Create grid (table) with info for several basins
        #---------------------------------------------------
        row_label_width = 60
        column_width    = 120
        table_width     = row_label_width + 4*(column_width)
        grid = wx.grid.Grid(panel, size=(table_width, 190))
        grid.CreateGrid(20, 4)
        grid.SetRowLabelSize(row_label_width)
        col_labels = ['Outlet Col', 'Outlet Row', \
                      'Area [km^2]', 'Relief [m]' ]
        for col in range(4):
            grid.SetColLabelValue(col, col_labels[col])
            grid.SetColSize(col, column_width)
            for row in range(20):
                grid.SetCellValue(row, col, "0.0")

        #---------------------------------
        # Create Read/Save table buttons
        #---------------------------------
        read_button = wx.Button(panel, -1, "Read table from file...")
        save_button = wx.Button(panel, -1, "Save table to file...")
        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        btn_sizer.Add(read_button)
        gap_size = (2*self.hgap, 2*self.hgap)
        btn_sizer.Add(gap_size)
        btn_sizer.Add(save_button) 
        self.Bind(wx.EVT_BUTTON, Read_Table, read_button)
        self.Bind(wx.EVT_BUTTON, Save_Table, save_button) 

        #-------------------------------------
        # Create the "RTM mask file" widgets
        #-------------------------------------
        checkbox = wx.CheckBox(panel, -1, \
                               "Check mass balance for basin 1?")
        self.basin_mass_balance_checkbox = checkbox   #########
        self.Bind(wx.EVT_CHECKBOX, Check_Mass_Balance, checkbox)
        #----------------------------------------------------------
        rtm_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.basin_RTM_file_box = \
             GUI_utils.Labeled_Text_Box(self, panel, rtm_sizer, \
                              "RTM mask file for basin 1:", \
                              self.basin_RTM_file, \
                              box_width=180)
        
        #--------------------------------------
        # Add objects to the basin_info_sizer
        #--------------------------------------
        border = self.vgap
        sizer.Add(grid,      0, wx.ALL, border)
        sizer.Add(btn_sizer, 0, wx.ALIGN_CENTER, border)
        pad_row = wx.StaticText(panel, -1, "")    
        sizer.Add(pad_row,   0, wx.ALL, border)
        sizer.Add(checkbox,  0, wx.ALL, border)
        sizer.Add(rtm_sizer, 0, wx.ALL, border)
        
        #------------------------------------------------
        # Add Back and Next buttons to basin_info_sizer
        #------------------------------------------------
        back_fcn = self.On_Goto_Methods
        next_fcn = self.On_Goto_Output_Log
        GUI_utils.Add_Back_Next_Buttons(self, panel, sizer, \
                                        back_fcn, next_fcn, \
                                        PAD_TOP=False)

        #-----------------------------------------
        # Add basin_info_sizer to the main_sizer
        # and then hide it
        #-----------------------------------------
        self.main_sizer.Add(sizer, 0, wx.ALL, self.vgap)
        self.main_sizer.Hide(sizer)
        
    #   Add_Basin_Info_Panel()
    #----------------------------------------------------------------
    def Add_Output_Log_Panel(self):

        #-----------------------------------------------------------
        def Set_Stop_Method(event):
            self.stop_method = self.stop_droplist.GetSelection()
        #-----------------------------------------------------------
        def Back_One(event):
            self.Show_Panel('basin_info')
        #-----------------------------------------------------------
        def Start_Model(event):
            self.Plotting_Test()  # (draw lines to plot window)
        #-----------------------------------------------------------
        def Clear_Window(event):
            self.output_log_box.Clear()
        #-----------------------------------------------------------
        def Get_Outfile_Size(event):
            print('Not ready yet.')
        #-----------------------------------------------------------
        def Stop_Model(event):
            print('Not ready yet.')
            
        #---------------------------------------------------
        # Create a "static box" and associated sizer for
        # everything but the bottom buttons.  Static boxes
        # provide a frame and a title for grouping.
        #---------------------------------------------------
        panel = self.main_panel
        box   = wx.StaticBox(panel, -1, "Output Log Window:")
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        self.output_log_sizer = sizer

        #------------------------------------------
        # Add stopping criterion options to sizer
        #------------------------------------------
        fg_sizer = wx.FlexGridSizer(cols=3, hgap=self.hgap, vgap=0)
        stop_choices = GUI_utils.Stop_Methods()
        self.stop_droplist = \
             GUI_utils.Labeled_Droplist(self, panel, fg_sizer, \
                               " Stopping criterion:", \
                               stop_choices, \
                               Set_Stop_Method, \
                               button1="Options...", \
                               b1_function=None)   ############
        sizer.Add(fg_sizer, 0, wx.ALL, 2*self.hgap)

        #----------------------------------------------------
        # Add some labels and hydrograph window in a column
        # to the left of a big text box for the output log
        #----------------------------------------------------
        border = 5
        mid_sizer = wx.BoxSizer(wx.HORIZONTAL)
        col_sizer = wx.BoxSizer(wx.VERTICAL)
        #----------------------------------------------------
        L1 = wx.StaticText(panel, -1, "Output log window:")
        col_sizer.Add(L1, 0, wx.ALL, border)        
        L2 = wx.StaticText(panel, -1, "")
        col_sizer.Add(L2, 0, wx.ALL, border)
        L3 = wx.StaticText(panel, -1, "")
        col_sizer.Add(L3, 0, wx.ALL, border)
##        L4 = wx.StaticText(panel, -1, "")
##        col_sizer.Add(L4, 0, wx.ALL, border)
        L5 = wx.StaticText(panel, -1, "Hydrograph:")
        col_sizer.Add(L5, 0, wx.ALL, border)

        #---------------------------------------
        # Add a plot window for the hydrograph
        #---------------------------------------
        # Google wx.Window for various options
        #---------------------------------------
        nx = 110
        ny = 110
        # win_style = (wx.SIMPLE_BORDER | wx.WS_EX_BLOCK_EVENTS)
        win_style = (wx.SIMPLE_BORDER)
        window = wx.Window(panel, -1, size=(nx,ny), \
                           style=win_style)
        window.SetBackgroundColour("white")
##        print 'window.IsDoubleBuffered =', window.IsDoubleBuffered
        self.Bind(wx.EVT_PAINT, None)   ############
        self.plot_window = window
        self.window_nx   = nx
        self.window_ny   = ny
        col_sizer.Add(window, 0, wx.ALL, border)

        #-------------------------------------------
        # Add the text box called the "output log"
        #-------------------------------------------        
        log_box = wx.TextCtrl(panel, -1, 'Clear Test', \
                              size=(360, 240))
        mid_sizer.Add(col_sizer, 0, wx.ALL, self.hgap)
        mid_sizer.Add(log_box,   0, wx.ALL, self.hgap)
        sizer.Add(mid_sizer, 0, wx.ALL, self.vgap)
        self.output_log_box = log_box
        
##        self.output_log_box = \
##             GUI_utils.Labeled_Text_Box(self, panel, fg_sizer, \
##                               "Output log window:", \
##                               box_width=360, \
##                               box_height=200, MULTI=True)         

        
        #------------------------
        # Add buttons at bottom
        #------------------------
        back_button  = wx.Button(panel, -1, "< Back")
        start_button = wx.Button(panel, -1, "Start Model Run")
        clear_button = wx.Button(panel, -1, "Clear Window")
        size_button  = wx.Button(panel, -1, "Get Outfile Size")
        stop_button  = wx.Button(panel, -1, "Stop")
        #--------------------------------------------------
        hgap = 5
        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        btn_sizer.Add(back_button)
        btn_sizer.Add((hgap, hgap), 1)
        btn_sizer.Add(start_button)
        btn_sizer.Add((hgap, hgap), 1)
        btn_sizer.Add(clear_button)
        btn_sizer.Add((hgap, hgap), 1)
        btn_sizer.Add(size_button)
        btn_sizer.Add((hgap, hgap), 1)
        btn_sizer.Add(stop_button)
        btn_sizer.Add((hgap, hgap), 1)
        #--------------------------------------------------------
        self.Bind(wx.EVT_BUTTON, Back_One,         back_button)
        self.Bind(wx.EVT_BUTTON, Start_Model,      start_button)
        self.Bind(wx.EVT_BUTTON, Clear_Window,     clear_button)
        self.Bind(wx.EVT_BUTTON, Get_Outfile_Size, size_button)
        self.Bind(wx.EVT_BUTTON, Stop_Model,       stop_button)
        #--------------------------------------------------------
        pad_row = wx.StaticText(panel, -1, "")
        sizer.Add(pad_row,   0, wx.ALL, hgap)
        sizer.Add(btn_sizer, 0, wx.ALL, hgap)
        
        #---------------------------------------
        # Add run_info_sizer to the main_sizer
        # and then hide it
        #---------------------------------------
        self.main_sizer.Add(sizer, 0, wx.ALL, self.vgap)
        self.main_sizer.Hide(sizer)
        
    #   Add_Output_Log_Panel()
    #----------------------------------------------------------------
    def Plotting_Test(self):

        #--------------------------------------------------------
        # Note:  We may need to save the window's bitmap in a
        #        buffer and refresh it on certain events.
        #        As it stands now, moving the cursor to another
        #        window causes window contents to be lost,
        #        even if we use the Freeze() method.
        #--------------------------------------------------------
        window = self.plot_window
        nx  = self.window_nx
        ny  = self.window_ny
        # win_buffer = self.plot_buffer
        win_buffer = wx.EmptyBitmap(nx,ny)
        #---------------------------------------------
        # Create a device context (don't store them)
        #---------------------------------------------
        dc = wx.BufferedDC(wx.ClientDC(window))
        ## dc = wx.BufferedDC(wx.ClientDC(window), win_buffer)
        ## dc = wx.ClientDC(window)  # (also works)
        ## dc = wx.WindowDC(window)  # (also works)
        pen   = wx.Pen("black", 2, wx.SOLID)
        brush = wx.Brush("white", wx.SOLID)  # (for filling in areas)
        dc.SetPen(pen)
        dc.SetBrush(brush)
        dc.SetBackground(brush)
        dc.Clear()
        #------------------------------------------
        dc.DrawRectangle(0,0,nx,ny)
        dc.DrawLine(0,0,nx,ny)
        dc.DrawCircle(nx/2,ny/2, nx/3)
        # print 'dc.GetSize() =', dc.GetSize()

        ## window.Freeze()   # (see also window.Thaw() )
        ## window.Disable()
        
    #   Plotting_Test()       
    #----------------------------------------------------------------
    def Add_Preprocessing_Panel(self):

        #---------------------------------------------------
        # Create a "static box" and associated sizer for
        # everything but the bottom buttons.  Static boxes
        # provide a frame and a title for grouping.
        #---------------------------------------------------
        panel = self.main_panel
        box   = wx.StaticBox(panel, -1, "Preprocessing Tools:")
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        self.preprocessing_sizer = sizer

        #---------------------------------------
        # Add run_info_sizer to the main_sizer
        # and then hide it
        #---------------------------------------
        self.main_sizer.Add(sizer, 0, wx.ALL, self.vgap)
        self.main_sizer.Hide(sizer)
        
    #   Add_Preprocessing_Panel()
    #----------------------------------------------------------------
    def Add_Plotting_Panel(self):

        #---------------------------------------------------
        # Create a "static box" and associated sizer for
        # everything but the bottom buttons.  Static boxes
        # provide a frame and a title for grouping.
        #---------------------------------------------------
        panel = self.main_panel
        box   = wx.StaticBox(panel, -1, "Plotting Options:")
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        self.plotting_sizer = sizer

        #---------------------------------------
        # Add run_info_sizer to the main_sizer
        # and then hide it
        #---------------------------------------
        self.main_sizer.Add(sizer, 0, wx.ALL, self.vgap)
        self.main_sizer.Hide(sizer)
        
    #----------------------------------------------------------------
    def Show_Panel(self, panel_name):

        #--------------------------------------------
        # Hide whichever panel is currently visible
        #--------------------------------------------
        arg1 = 'self.' + self.current_panel + '_sizer'
        exec('self.main_sizer.Hide(' + arg1 + ')')

        #---------------------
        # Show the new panel
        #---------------------
        arg2 = 'self.' + panel_name + '_sizer'
        exec('self.main_sizer.Show(' + arg2 + ')')
        self.current_panel = panel_name   ####
        
        #----------------------------------
        # Readjust the layout (need this)
        #----------------------------------
        self.main_sizer.Layout()

        #-------------------------------------
        # We're using a fixed panel size now
        #-------------------------------------
        ## self.main_sizer.Fit(self)
             
    #   Show_Panel()
    #----------------------------------------------------------------
    #   EVENT HANDLERS BELOW THIS POINT
    #----------------------------------------------------------------    
    def On_File_Open(self, event):

        wx.MessageBox("You selected: File > Open Data Set")

    #   On_File_Open()
    #----------------------------------------------------------------
    def On_File_Load(self, event):

        wx.MessageBox("You selected: File > Load Input Vars")

    #   On_File_Load()
    #----------------------------------------------------------------
    def On_File_Save(self, event):

        wx.MessageBox("You selected: File > Save Input Vars")

    #   On_File_Save()
    #----------------------------------------------------------------
    def On_File_Close(self, event):

        wx.MessageBox("You selected: File > Close Any Open Files")

    #   On_File_Close()
    #----------------------------------------------------------------
    def On_File_Exit(self, event):

        dlg = wx.MessageDialog(self,
                "Do you really want to close this application?",
                "Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
        result = dlg.ShowModal()
        dlg.Destroy()
        if (result == wx.ID_OK):
            self.Destroy()
  
        ## self.Destroy()

    #   On_File_Exit()
    #----------------------------------------------------------------
    def On_Goto_Navigation(self, event):

        self.Show_Panel('navigation')

    #   On_Goto_Navigation()
    #----------------------------------------------------------------
    def On_Goto_Run_Info(self, event):

        self.Show_Panel('run_info')

    #   On_Goto_Run_Info()
    #----------------------------------------------------------------
    def On_Goto_Methods(self, event):

        self.Show_Panel('methods')

    #   On_Goto_Methods()
    #----------------------------------------------------------------
    def On_Goto_Basin_Info(self, event):

        self.Show_Panel('basin_info')

    #   On_Goto_Basin_Info()
    #----------------------------------------------------------------
    def On_Goto_Output_Log(self, event):

        self.Show_Panel('output_log')

    #   On_Goto_Output_Log()
    #----------------------------------------------------------------
    def On_Goto_Preprocessing(self, event):

        self.Show_Panel('preprocessing')
        
        # wx.MessageBox("You selected: Goto > Preprocessing")

    #   On_Goto_Preprocessing()
    #----------------------------------------------------------------
    def On_Goto_Plotting(self, event):

        self.Show_Panel('plotting')
        
        # wx.MessageBox("You selected: Goto > Plot Results")

    #   On_Goto_Plotting()
    #----------------------------------------------------------------
    def On_Help_About(self, event):

        help_file = self.help_dir + 'about_TF.htm'
        result = webbrowser.open('file://' + help_file)
        
    #   On_Help_About()
    #----------------------------------------------------------------
    def On_Help_License(self, event):

        help_file = self.help_dir + 'TF_License_Agreement.htm'
        result = webbrowser.open('file://' + help_file)
        
    #   On_Help_License()
    #----------------------------------------------------------------
    def On_Help_Whats_New(self, event):

        help_file = self.help_dir + 'TF15b_Release_Notes.htm'
        result = webbrowser.open('file://' + help_file)
        
    #   On_Help_Whats_New()
    #----------------------------------------------------------------
    def On_Help_Tutorial(self, event):

        help_file = self.help_dir + 'TF_Tutorial.htm'
        result = webbrowser.open('file://' + help_file)
        
    #   On_Help_Tutorial()
    #----------------------------------------------------------------
    def On_Precip_In_Button(self, event):
        
        if (self.p_method == 0):
            msg = "This method has no input variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.p_method == 1):
            xml_file = "precip_uniform_space.xml"
        elif (self.p_method == 2):
            xml_file = "precip_general.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return
        
        xml_file = (self.dialog_in_dir + xml_file)
        dialog = TF_Input_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_Snowmelt_In_Button(self, event):
        
        if (self.s_method == 0):
            msg = "This method has no input variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.s_method == 1):
            xml_file = "snowmelt_degree_day.xml"
        elif (self.s_method == 2):
            xml_file = "snowmelt_energy_balance.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return
        
        xml_file = (self.dialog_in_dir + xml_file)
        dialog = TF_Input_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_ET_In_Button(self, event):
        
        if (self.e_method == 0):
            msg = "This method has no input variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.e_method == 1):
            xml_file = "ET_priestley_taylor.xml"
        elif (self.e_method == 2):
            xml_file = "ET_energy_balance.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return
        
        xml_file = (self.dialog_in_dir + xml_file)
        dialog = TF_Input_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_Infil_In_Button(self, event):
        
        if (self.i_method == 0):
            msg = "This method has no input variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.i_method == 1):
            msg = "This method has no input variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.i_method == 2):
            xml_file = "infil_green_ampt.xml"
        elif (self.i_method == 3):
            xml_file = "infil_smith_parlange.xml"
        elif (self.i_method == 4):
            xml_file = "infil_richards_1D.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return
        
        xml_file = (self.dialog_in_dir + xml_file)
        dialog = TF_Input_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_GW_In_Button(self, event):
        
        if (self.g_method == 0):
            msg = "This method has no input variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.g_method == 1):
            xml_file = "GW_darcy_layers.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            
        xml_file = (self.dialog_in_dir + xml_file)
        dialog = TF_Input_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_Channel_In_Button(self, event):
        
        if (self.c_method == 0):
            msg = "This method has no input variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.c_method == 1):
            xml_file = "channel_kinematic_manning.xml"
        elif (self.c_method == 2):
            xml_file = "channel_kinematic_law_of_wall.xml"
        elif (self.c_method == 3):
            xml_file = "channel_diffusive_manning.xml"
        elif (self.c_method == 4):
            xml_file = "channel_diffusive_law_of_wall.xml"
        elif (self.c_method == 5):
            xml_file = "channel_dynamic_manning.xml"
        elif (self.c_method == 6):
            xml_file = "channel_dynamic_law_of_wall.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return
        
        xml_file = (self.dialog_in_dir + xml_file)
        dialog = TF_Input_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_Diversion_In_Button(self, event):
        
        if (self.d_method == 0):
            msg = "This method has no input variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.d_method == 1):
            xml_file = "diversions.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return
        
        xml_file = (self.dialog_in_dir + xml_file)
        dialog = TF_Input_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_Snowmelt_Out_Button(self, event):
        
        if (self.s_method == 0):
            msg = "This method has no output variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.s_method == 1):
            xml_file = "snowmelt_degree_day_out.xml"
        elif (self.s_method == 2):
            xml_file = "snowmelt_energy_balance_out.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return

        xml_file = (self.dialog_out_dir + xml_file)        
        dialog = TF_Output_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_ET_Out_Button(self, event):
        
        if (self.e_method == 0):
            msg = "This method has no output variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.e_method == 1):
            xml_file = "ET_priestley_taylor_out.xml"
        elif (self.e_method == 2):
            xml_file = "ET_energy_balance_out.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return

        xml_file = (self.dialog_out_dir + xml_file)          
        dialog = TF_Output_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_Infil_Out_Button(self, event):
        
        if (self.i_method == 0):
            msg = "This method has no output variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.i_method == 1):
            msg = "This method has no output variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.i_method == 2):
            xml_file = "infil_green_ampt_out.xml"
        elif (self.i_method == 3):
            xml_file = "infil_smith_parlange_out.xml"
        elif (self.i_method == 4):
            xml_file = "infil_richards_1D_out.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return

        xml_file = (self.dialog_out_dir + xml_file)  
        dialog = TF_Output_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_GW_Out_Button(self, event):
        
        if (self.g_method == 0):
            msg = "This method has no output variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.g_method == 1):
            xml_file = "GW_darcy_layers_out.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return

        xml_file = (self.dialog_out_dir + xml_file)  
        dialog = TF_Output_Dialog(main_frame=self, xml_file=xml_file)
        
    #-----------------------------------------------------------------
    def On_Channel_Out_Button(self, event):
        
        if (self.c_method == 0):
            msg = "This method has no output variables."
            wx.MessageBox(msg, caption='SORRY,')
            return
        elif (self.c_method == 1):
            xml_file = "channel_out.xml"
        elif (self.c_method == 2):
            xml_file = "channel_out.xml"
        elif (self.c_method == 3):
            xml_file = "channel_out.xml"
        elif (self.c_method == 4):
            xml_file = "channel_out.xml"
        elif (self.c_method == 5):
            xml_file = "channel_out.xml"
        elif (self.c_method == 6):
            xml_file = "channel_out.xml"
        else:
            msg = "Method number out of range."
            wx.MessageBox(msg, caption='SORRY,')
            return

        xml_file = (self.dialog_out_dir + xml_file)  
        dialog = TF_Output_Dialog(main_frame=self, xml_file=xml_file)
        
#------------------------------------------------------------------------  
class TF_GUI(wx.App):

    def OnInit(self):
    ## def __init__(self):
        frame = TF_Main_Dialog()
        self.SetTopWindow(frame)
        frame.Show()
        return True
    
#---------------------------------------
#  Support two different usage options
#---------------------------------------
if (__name__ == '__main__'):
    app = wx.App()
    frame  = TF_Main_Dialog()
    frame.Show()
    app.MainLoop()

    #----------------------------------
    # Alternate approach, not working
    #----------------------------------
    ## app = TF_GUI()
    ## app.MainLoop()
    
    ## app = wx.PySimpleApp()
    ## app = wx.App()
    #--------------------------------------
    # Didn't help with Mac menu bar issue
    #--------------------------------------
##    frame0 = TF_Hidden_Dialog()
##    frame0.Show()
##    frame  = TF_Main_Dialog(parent=frame0
    #----------------------------------------
    ## frame  = TF_Main_Dialog()
    ## frame.Show()
    ## app.MainLoop()

