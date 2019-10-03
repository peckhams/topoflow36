#!/usr/bin/env python

#  April 21-22, 2009
#  S.D. Peckham

import wx
import webbrowser  ##  standard Python module
import sys
import time
    
## import wx.html
## import xml.dom.minidom

#  For wizard-style panels, use:
#     self.panel1.Show(False)
#     self.panel2.Show(True)
#  or perhaps instead use:
#     self.sizer1.Show(False)
#     self.sizer2.Show(True)

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
    def __init__(self, parent=None, id=-1):
   
        #---------------------------
        # Set resource directories
        #---------------------------
        if (sys.platform == 'win32'):
            self.home_dir = "C:\Program Files\TopoFlow\\"
        elif (sys.platform == 'darwin'):
            self.home_dir = '/Applications/TopoFlow/'
        else:
            self.home_dir = '/usr/local/topoflow/'
        self.help_dir = self.home_dir + 'help/'
        self.bmp_dir  = self.home_dir + 'Images/'

        #----------------------------------
        # Used for spacing in some panels
        #----------------------------------
        self.vgap = 10
        self.hgap = 6
        
        self.title = "TopoFlow 2.0 beta Dialog (4/21/09)"
        wx.Frame.__init__(self, parent, id, self.title)

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

        #-------------------------------------
        # Create the wizard panels, but hide
        # all but the navigation panel
        #-------------------------------------
        self.Get_Navigation_Panel()
        self.Get_Run_Info_Panel()
        #------------------------------    
        # self.Get_Methods_Panel()
        # self.Get_Basin_Info_Panel()
        # self.Get_Output_Log_Panel()
        # self.Get_Preprocessing_Panel()
        #---------------------------------------------
        self.run_info_panel.Hide()
##        self.methods_panel.Hide()
##        self.basin_info_panel.Hide()
##        self.output_log_panel.Hide()
##        self.preprocessing_panel.Hide()
        #---------------------------------------------
        # self.Show(True)  # (the frame, not a panel)
        
    #   __init__()
    #----------------------------------------------------------------
    def Get_Navigation_Panel(self):

        #----------------------------------
        # Create a new panel with a sizer
        #----------------------------------
        panel = wx.Panel(self)
        panel.SetBackgroundColour('Light Blue')
        self.navigation_panel = panel
        sizer = wx.GridSizer(2, 2, self.hgap, self.vgap)
        panel.SetSizer(sizer)
        
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
        flag   = wx.ALIGN_CENTER  # (affects resizing by user)
        border = 10
        sizer.Add(new_button,  proportion, flag, border)
        sizer.Add(pre_button,  proportion, flag, border)
        sizer.Add(plot_button, proportion, flag, border)
        sizer.Add(exit_button, proportion, flag, border)
        #------------------------------------------------------------------
        self.Bind(wx.EVT_BUTTON, self.On_Goto_New_Run,       new_button)
        self.Bind(wx.EVT_BUTTON, self.On_Goto_Preprocessing, pre_button)
        self.Bind(wx.EVT_BUTTON, self.On_Goto_Plotting,      plot_button)
        self.Bind(wx.EVT_BUTTON, self.On_File_Exit,          exit_button)

        #------------------------------------
        # Main panel is displayed at start,
        # but other panels are hidden
        #------------------------------------
        sizer.Fit(self)
        panel.Show()
        ## sizer.Show(True)    # (also works)
        
    #   Get_Navigation_Panel()
    #----------------------------------------------------------------
    def Get_Run_Info_Panel(self):

        #####################################
        self.run_directory = self.home_dir
        self.data_prefix   = 'Treynor'
        #####################################

        #---------------------
        # Create a new panel
        #---------------------
        panel = wx.Panel(self)
        panel.SetBackgroundColour('Light Blue')
        self.run_info_panel = panel
        
        #---------------------------------------------------
        # Create a "static box" and associated sizer for
        # everything but the bottom buttons.  Static boxes
        # provide a frame and a title for grouping.
        #---------------------------------------------------
        box   = wx.StaticBox(panel, -1, "Model run information:")
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        panel.SetSizer(sizer)
        
        #---------------------------------------------
        #  Create another sizer box for rows of info
        #---------------------------------------------
        fg_sizer = wx.FlexGridSizer(cols=3, hgap=self.hgap, vgap=0)

        self.run_info_text_boxes = []
        tbox_width = 250
        
        #------------------------
        # Labels and text boxes
        #------------------------
        L1 = wx.StaticText(panel, -1, "Model run directory:")
        T1 = wx.TextCtrl(panel, -1, self.run_directory, \
                         size=(tbox_width,-1))
        B1 = wx.Button(panel, -1, "Browse...")
        # self.Bind(wx.EVT_BUTTON, self.On_Run_Info_Browse, B1)  ####
        fg_sizer.Add(L1, 0, wx.ALL, self.hgap)
        fg_sizer.Add(T1, 0, wx.ALL, self.hgap)
        fg_sizer.Add(B1, 0, wx.ALL, self.hgap)
        #--------------------------------------
        self.run_info_text_boxes.append( T1 )
        #----------------------------------------------------------
        L2 = wx.StaticText(panel, -1, "Data set prefix:")
        T2 = wx.TextCtrl(panel, -1, self.data_prefix, \
                         size=(tbox_width,-1))
        B2 = wx.StaticText(panel, -1, "")
        fg_sizer.Add(L2, 0, wx.ALL, self.hgap)
        fg_sizer.Add(T2, 0, wx.ALL, self.hgap)
        fg_sizer.Add(B2, 0, wx.ALL, self.hgap)
        #--------------------------------------
        self.run_info_text_boxes.append( T2 )

        #---------------------------------
        # Add fg_sizer to the main sizer
        #------------------------------------------------
        # Next 2 lines adds padding inside of StaticBox
        #------------------------------------------------
        # pad_row = wx.StaticText(panel, -1, " ")
        # sizer.Add(pad_row,  0, wx.ALL, self.vgap)
        sizer.Add(fg_sizer, 0, wx.ALL, self.vgap)
        
        #-----------------------------------
        # Nav panel is displayed at start,
        # but all other panels are hidden
        #-----------------------------------
        sizer.Fit(self)
        panel.Show()
        # panel.Hide()
        
    #   Get_Run_Info_Panel()
    #----------------------------------------------------------------
    def Get_Methods_Panel(self):

        panel = wx.Panel(self)
        panel.SetBackgroundColour('Light Blue')
        self.methods_panel = panel

        #-----------------------------------
        # Nav panel is displayed at start,
        # but all other panels are hidden
        #-----------------------------------
        panel.Show(False)
        
    #   Get_Methods_Panel()
    #----------------------------------------------------------------
    def Get_Basin_Info_Panel(self):

        panel = wx.Panel(self)
        panel.SetBackgroundColour('Light Blue')
        self.basin_info_panel = panel

        #-----------------------------------
        # Nav panel is displayed at start,
        # but all other panels are hidden
        #-----------------------------------
        panel.Show(False)
        
    #   Get_Basin_Info_Panel()
    #----------------------------------------------------------------
    def Get_Output_Log_Panel(self):

        panel = wx.Panel(self)
        panel.SetBackgroundColour('Light Blue')
        self.output_log_panel = panel

        #-----------------------------------
        # Nav panel is displayed at start,
        # but all other panels are hidden
        #-----------------------------------
        panel.Show(False)
        
    #   Get_Output_Log_Panel()
    #----------------------------------------------------------------
    def Get_Preprocessing_Panel(self):

        panel = wx.Panel(self)
        panel.SetBackgroundColour('Light Blue')
        self.preprocessing_panel = panel

        #-----------------------------------
        # Nav panel is displayed at start,
        # but all other panels are hidden
        #-----------------------------------
        panel.Show(False)
        
    #   Get_Preprocessing_Panel()
    #----------------------------------------------------------------
    def Get_Plotting_Panel(self):

        panel = wx.Panel(self)
        panel.SetBackgroundColour('Light Blue')
        self.plotting_panel = panel

        #-----------------------------------
        # Nav panel is displayed at start,
        # but all other panels are hidden
        #-----------------------------------
        panel.Show(False)
        
    #   Get_Plotting_Panel()
    #----------------------------------------------------------------
    #   EVENT HANDLERS BELOW THIS POINT
    #----------------------------------------------------------------    
    def On_File_Open(self, event):

        self.navigation_panel.Show(False)  #### for testing ####
        wx.MessageBox("You selected: File > Open Data Set")

    #   On_File_Open()
    #----------------------------------------------------------------
    def On_File_Load(self, event):

        self.navigation_panel.Show(True)   #### for testing ####
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

        self.Destroy()

    #   On_File_Exit()
    #----------------------------------------------------------------
    def On_Goto_Navigation(self, event):

        self.run_info_panel.Show(False)
        self.navigation_panel.Show(True)
        
        # wx.MessageBox("You selected: Goto > New Model Run")

    #   On_Goto_Navigation()
    #----------------------------------------------------------------
    def On_Goto_New_Run(self, event):

        self.navigation_panel.Hide()
        self.Get_Run_Info_Panel()   #######
        self.run_info_panel.Show()
        
        #self.navigation_panel.Hide()
        #self.run_info_panel.Show()
        #self.Layout()
        
        ## self.run_info_panel.Layout()  ##### (didn't help)
        
        # self.Show(False)  # (this hides entire frame)
        
        # wx.MessageBox("You selected: Goto > New Model Run")

    #   On_Goto_New_Run()
    #----------------------------------------------------------------
    def On_Goto_Output_Log(self, event):

        self.navigation_panel.Hide()
        # self.navigation_panel.Show(False)
        self.output_log_panel.Show(True)
        
        # wx.MessageBox("You selected: Goto > New Model Run")

    #   On_Goto_Output_Log()
    #----------------------------------------------------------------
    def On_Goto_Preprocessing(self, event):

        self.navigation_panel.Show(False)
        self.run_info_panel.Show(False)
        self.plotting_panel.Show(False)
        self.preprocessing_panel.Show(True)

        # wx.MessageBox("You selected: Goto > Preprocessing")

    #   On_Goto_Preprocessing()
    #----------------------------------------------------------------
    def On_Goto_Plotting(self, event):

        self.navigation_panel.Show(False)
        self.preprocessing_panel.Show(False)
        self.plotting_panel.Show(True)
        
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

#-------------------------------------------------------------------------    
##class TF_Hidden_Dialog(wx.Frame):
##    def __init__(self, parent=None, id=-1):
##   
##        wx.Frame.__init__(self, parent, id, "")
##        self.panel = wx.Panel(self)
##        # self.Show(True)

#-------------------------------------------------------------------------  
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

