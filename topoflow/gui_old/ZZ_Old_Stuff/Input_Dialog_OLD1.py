#!/usr/bin/env python
#  August 8, 2008
#  S.D. Peckham

import wx
import wx.html

# class TF_Input_Dialog(wx.Frame, title="Model Input Dialog"):
    
class TF_Input_Dialog(wx.Frame):

    #-----------------------------------------------------
    # Notes:  Default is wx.ALIGN_LEFT for labels & text
    #-----------------------------------------------------
    # title = "Model Input Dialog"
    # model_input_types = ["Scalar", "Time series", "Grid", "Grid sequence"]
    # timestep_str   = "1.0"
    # timestep_units = "[seconds / timestep]"
    # timestep_msg   = "Snowmelt process timestep:"
    # timestep = state.timestep
        
    #---------------------------
    #  Create top-level dialog
    #---------------------------
    def __init__(self, parent, id):
        self.title     = "Model Input Dialog"
        self.var_file  = "my_var_file.xml"
        self.help_file = "my_help_file.html"
        self.pad       = 10
        # self.xsize = 300
        # self.ysize = 200
        wx.Frame.__init__(self, parent, id, self.title)
                          # size=(self.xsize, self.ysize))
        self.panel = wx.Panel(self, -1)
        self.panel.SetBackgroundColour('White')
        #------------------------------------
        self.read_var_info()
        blank1 = wx.StaticLine(self.panel)
        hbox = self.header_box()
        base = wx.BoxSizer(wx.VERTICAL)
        base.Add((10,10), 1)
        base.Add(hbox, 0, wx.GROW)
        base.Add((10,10), 1)
        
        for row in range(0,3):
            row_box = self.row_box(row)
            base.Add(row_box)
            base.Add((10,10), 1)
                    
        time_box   = self.timestep_box()
        button_bar = self.button_bar()
        blank2 = wx.StaticLine(self.panel)
        #------------------------------------------       
        # base = wx.BoxSizer(wx.VERTICAL)
        # base.Add(blank1, 0, wx.GROW|wx.TOP|wx.BOTTOM, 5)
        # base.Add((10,10), 1)
        # base.Add(hbox, 0, wx.GROW)
        ## for row in range(0,3):
        ##     base.Add(row_box)
        base.Add((10,10), 1)
        base.Add(time_box)
        base.Add((10,10), 1)
        base.Add(button_bar, 0, wx.GROW)
        base.Add((10,10), 1)
        # base.Add(blank2, 0, wx.GROW|wx.TOP|wx.BOTTOM, 5)
        
##        base.Add(self.create_header_box(panel))
##        for row in range(0,3):
##            base.Add(self.create_row_box(panel, row))
##        base.Add(self.create_timestep_box(panel))
##        base.Add(self.create_button_bar(panel))

        self.panel.SetSizer(base)
        base.Fit(self)
        base.SetSizeHints(self)

        # panel.Layout()

##        self.SetSizer(base)
##        self.Layout()
        
    #--------------------------------------------------
    #  Read descriptions of input variables from file
    #--------------------------------------------------
    def read_var_info(self):
        self.var_names = ["c0:", "T0:", "T_air:"]
        self.var_units = ["[mm/day/deg_C]", "[deg_C]", "[deg_C]"]
        self.var_setting = ["1.0", "20.0", "25.0"]
        self.var_types = [0, 0, 0]
        #self.var_types = [state.snow_vars.c0_type,
        #                  state.snow_vars.T0_type,
        #                  state.met_vars.T_air_type]

    #-------------------------------------------
    #  Create sizer box for the column headers
    #-------------------------------------------
    def header_box(self):
        L1  = wx.StaticText(self.panel, -1, "Variable: ")
        L2  = wx.StaticText(self.panel, -1, "Type: ")
        L3  = wx.StaticText(self.panel, -1, "Scalar or Grid Filename: ")
        L4  = wx.StaticText(self.panel, -1, "Units: ")

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        # sizer.Add((10,10), 1)
        sizer.Add(L1, 0, wx.GROW, border=self.pad)
        sizer.Add(L2, 0, wx.GROW, border=self.pad)
        sizer.Add(L3, 0, wx.GROW, border=self.pad)
        sizer.Add(L4, 0, wx.GROW, border=self.pad)

        # Next line is bad for some reason.
        # self.SetSizer(sizer)
        return sizer

##        box.AddMany([L1,L2,L3,L4]) # , 1, wx.EXPAND)
        
##        box.Add(L1, 1, wx.EXPAND)
##        box.Add(L2, 1, wx.EXPAND)
##        box.Add(L3, 1, wx.EXPAND)
##        box.Add(L4, 1, wx.EXPAND)

        
    #------------------------------------------------
    #  Create sizer box for a single input variable
    #------------------------------------------------
    def row_box(self, row):

        label = wx.TextCtrl(self.panel, -1, value=self.var_names[row])
        model_input_types = ["Scalar", "Time series", "Grid", "Grid sequence"]
        dlist = wx.Choice(self.panel, -1, choices=model_input_types)
        #  Set to current input type here
        text  = wx.TextCtrl(self.panel, -1, self.var_setting[row])
        ustr  = wx.StaticText(self.panel, -1, self.var_units[row])
        # self.Bind(wx.EVT_CHOICE, self.onChoice, self.var_type)

        box = wx.BoxSizer(wx.HORIZONTAL)
        box.Add(label, 0, wx.GROW, border=self.pad)
        box.Add(dlist, 0, wx.GROW, border=self.pad)
        box.Add(text,  0, wx.GROW, border=self.pad)
        box.Add(ustr,  0, wx.GROW, border=self.pad)
        # box.AddMany([label, dlist, text, ustr])  #, 1, wx.EXPAND) 
        # self.SetSizer(box)
        return box
    
    #---------------------------------------------
    #  Create sizer box for the process timestep
    #---------------------------------------------
    def timestep_box(self):
        timestep_str   = "1.0"
        timestep_units = "[seconds / timestep]"
        timestep_msg   = "Snowmelt process timestep:"
        #-----------------------------------------------
        L1    = wx.StaticText(self.panel, -1, timestep_msg)
        text  = wx.TextCtrl(self.panel, -1, timestep_str)
        L2    = wx.StaticText(self.panel, -1, timestep_units)
        box   = wx.BoxSizer(wx.HORIZONTAL)
        box.Add((10,10), 1)
        box.Add(L1)
        box.Add((5,5), 1)
        box.Add(text)
        box.Add((5,5), 1)
        box.Add(L2)
        box.Add((10,10), 1)
        # box.AddMany([L1, text, L2]) #, 1, wx.EXPAND)

        # self.SetSizer(box)
        return box
    
    #----------------------------
    #  Create bottom button bar
    #----------------------------
    def button_bar(self):
        start_btn  = wx.Button(self.panel, -1, "Start")
        help_btn   = wx.Button(self.panel, -1, "Help")
        cancel_btn = wx.Button(self.panel, -1, "Cancel")
        
        box = wx.BoxSizer(wx.HORIZONTAL)
        box.Add((10,10), 1)
        box.Add(start_btn)
        box.Add((10,10), 1)
        box.Add(help_btn)
        box.Add((10,10), 1)
        box.Add(cancel_btn)
        box.Add((10,10), 1)
                
        # box.AddMany([start_btn, help_btn, cancel_btn]) #, 1, wx.EXPAND)
##        box.Add(start_btn,  1, wx.EXPAND)
##        box.Add(help_btn,   1, wx.EXPAND)
##        box.Add(cancel_btn, 1, wx.EXPAND)
        
##        bar = wx.Panel(parent)
##        self.start_btn  = wx.Button(bar, -1, "Start")
##        self.help_btn   = wx.Button(bar, -1, "Help")
##        self.cancel_btn = wx.Button(bar, -1, "Cancel")
        #-------------------------------------------------
        self.Bind(wx.EVT_BUTTON, self.on_Start,  start_btn)
        self.Bind(wx.EVT_BUTTON, self.on_Help,   help_btn)
        self.Bind(wx.EVT_BUTTON, self.on_Cancel, cancel_btn)
        
##        self.Bind(wx.EVT_BUTTON, self.on_Start,  self.start)
##        self.Bind(wx.EVT_BUTTON, self.on_Help,   self.help)
##        self.Bind(wx.EVT_BUTTON, self.on_Cancel, self.cancel)

        # self.SetSizer(box)
        return box
    
    #-------------------------
    #  Event handler methods
    #-------------------------
    def on_Start(self, event):
        pass

    def on_Help(self, event):
        pass

    def on_Cancel(self, event):
        self.Destroy()
    
#---------------------------------------
#  Support two different usage options
#---------------------------------------
if __name__ == '__main__':
    print("Creating dialog...")
    app = wx.PySimpleApp()
    frame = TF_Input_Dialog(parent=None, id=-1)
    frame.Show()
    app.MainLoop()
