#!/usr/bin/env python

#  April 27, 2009
#  S.D. Peckham

import wx
from TF_Input_Var_Box import *

#-----------------------------------------------------------------------
class Layer_Data():
    def __init__(self):
        tlist = ['Scalar', 'Time series', 'Grid', 'Grid sequence']
        #-----------------------------------------------------------
        self.var_names   = ['', '', '', '']
        self.var_labels  = ['K_s', 'K_i', 'theta_s', 'theta_i']
        self.var_types   = ['Scalar', 'Scalar', 'Scalar', 'Scalar']
        self.var_values  = ['7.20e-06', '9.849e-08', '0.485', '0.3758']
        self.var_units   = ['m/s', 'm/s', 'none', 'none']
        self.var_type_choices = [tlist, tlist, tlist, tlist]
        
#-----------------------------------------------------------------------
class TF_Notebook_Test(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, title="TF Notebook Test")
                          ## style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        self.SetBackgroundColour('Light Blue')
        
        #-----------------------------
        # Create a "tabbed" notebook
        #-----------------------------
##        notebook = wx.Notebook(self, style=wx.BK_DEFAULT)
        notebook = wx.Notebook(self, style=wx.NB_TOP)
##        notebook = wx.Notebook(self, style=wx.BK_TOP)
##        notebook = wx.Notebook(self, style=wx.BK_BOTTOM)
##        notebook = wx.Notebook(self, style=wx.BK_LEFT)
##        notebook = wx.Notebook(self, style=wx.BK_RIGHT)
        
        #-----------------------------
        # Create data for each layer
        #-----------------------------
        layer1_data = Layer_Data()
        layer2_data = Layer_Data()
        layer3_data = Layer_Data()
        
        #--------------------------------------
        # Create notebook "pages" as children
        #--------------------------------------
        page1 = TF_Input_Var_Box(parent=notebook, data=layer1_data, \
                                 box_label="Layer 1 variables:")
        page2 = TF_Input_Var_Box(parent=notebook, data=layer2_data, \
                                 box_label="Layer 2 variables:")
        page3 = TF_Input_Var_Box(parent=notebook, data=layer3_data, \
                                 box_label="Layer 3 variables:")
        
        #------------------------------------------
        # Add pages to notebook with "tab labels"
        #------------------------------------------
        notebook.AddPage(page1, "Layer 1")
        notebook.AddPage(page2, "Layer 2")
        notebook.AddPage(page3, "Layer 3")

        #----------------------------------------------
        # Create sizer box for the process timestep
        # that is *not* part of the "tabbed notebook"
        #----------------------------------------------
        L1    = wx.StaticText(self, -1, 'Infiltration time step:')
        text  = wx.TextCtrl(self,   -1, '60.0')
        L2    = wx.StaticText(self, -1, '[ seconds / timestep ]')
        #------------------------------------------------------------
        hgap = 10
        vgap = 6
        ts_box   = wx.BoxSizer(wx.HORIZONTAL)
        ts_box.Add((3*hgap, 3*hgap), 1)
        ts_box.Add(L1)
        ts_box.Add((hgap, hgap), 1)
        ts_box.Add(text)
        ts_box.Add((hgap, hgap), 1)
        ts_box.Add(L2)
        
        #-----------------------------------------------------
        # Put notebook in a sizer for panel to manage layout
        #-----------------------------------------------------
        frame_sizer = wx.BoxSizer(wx.VERTICAL)
        frame_sizer.Add(notebook, 0, wx.EXPAND|wx.ALL, 5)
        frame_sizer.Add(ts_box, 0, wx.EXPAND|wx.ALL, 5)  ##########
        self.SetSizer(frame_sizer)
        self.Fit()    # (need this)
        ### self.Centre()

#-------------------------------------------------------------
if __name__ == "__main__":
    app = wx.App()
    frame = TF_Notebook_Test()
    frame.Show()
    app.MainLoop()
    

