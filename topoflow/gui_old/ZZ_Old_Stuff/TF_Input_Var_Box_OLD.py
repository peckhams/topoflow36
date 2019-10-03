#!/usr/bin/env python

#  April 27, 2009
#  S.D. Peckham

import wx

#-------------------------------------------------------------
class TF_Input_Var_Box(wx.Panel):

    #-----------------------------------------------------
    # Notes:  This class is for creating "input variable
    #         panels" within TopoFlow input dialogs.
    #         Initial settings, labels, etc. are read
    #         from the XML file provided
    #         Default is wx.ALIGN_LEFT for labels & text
    #-----------------------------------------------------
    def __init__(self, parent=None, id=-1, \
                 box_label="Input variables:", data=None):
 
        wx.Panel.__init__(self, parent)

        #------------------------------------------------
        # Saving parent allows collected values to be
        # stored in parent frame before this one closes.
        #------------------------------------------------
        self.parent     = parent
        self.data       = data       ###########
        #--------------------------
        self.vgap       = 10
        self.hgap       = 6
        self.text_width = 160
        self.type_code = {'Scalar':0, 'Time series':1, \
                          'Grid':2, 'Grid sequence':3}
        self.type_name = {0:'Scalar', 1:'Time series', \
                          2:'Grid', 3:'Grid sequence'}
        
        #-------------------------------------------
        #  Create sizer box for all input var info
        #  (This provides a frame and "box label".)
        #-------------------------------------------
        vbox  = wx.StaticBox(self, -1, box_label)
        sizer = wx.StaticBoxSizer(vbox, wx.VERTICAL)
        
        #---------------------------------------------
        #  Create another sizer box for rows of info
        #---------------------------------------------
        #  Use "vgap=0" for most compact dialog
        #---------------------------------------------
        header = ["Variable:", "Type:", "Scalar or Grid Filename:",
                  "Units:"]
        nh = len(header)
        fg_sizer = wx.FlexGridSizer(cols=nh, hgap=self.hgap, vgap=0)
        
        #----------------------------------------------       
        #  Specify which columns can expand on resize
        #----------------------------------------------
        # fg_sizer.AddGrowableCol(1)

        ####################################        
        self.text_boxes = []
        ####################################

        #-------------------------
        # Add the column headers
        #-------------------------
        for row in range(nh):
            L1 = wx.StaticText(self, -1, header[row])
            fg_sizer.Add(L1, 0, wx.ALL, self.hgap)

        if (data == None):
            return  #################################
        
        #-----------------------------------------------
        # Create a row in the dialog for each variable
        # using names, types, values and units.
        #-----------------------------------------------
        for row in range(len(data.names)):
            vstr  = data.names[row] + ": "
            label = wx.StaticText(self, -1, vstr)
            #----------------------------------------------------------
            row_ID = (5000 + row)   #####
            ### row_ID = 'tf_row_' + str(row)  # (ID can't be string.)
            dlist = wx.Choice(self, row_ID, choices=data.typelists[row])
            dlist.Select(self.type_code[data.types[row]])
            self.Bind(wx.EVT_CHOICE, self.On_Type_Choice, dlist)      
            #----------------------------------------------------------            
            text  = wx.TextCtrl(self, -1, data.values[row],
                                size=(self.text_width,-1))
            ####################################
            self.text_boxes.append( text )
            ####################################
            #----------------------------------------------------------
            ustr  = wx.StaticText(self, -1, data.units[row])
            #----------------------------------------------------------
            fg_sizer.Add(label, 1, wx.EXPAND|wx.ALL, self.hgap)
            fg_sizer.Add(dlist, 1, wx.EXPAND|wx.ALL, self.hgap)
            fg_sizer.Add(text,  1, wx.EXPAND|wx.ALL, self.hgap)
            fg_sizer.Add(ustr,  1, wx.EXPAND|wx.ALL, self.hgap)

        #---------------------------------
        # Add fg_sizer to the main sizer
        #---------------------------------
        sizer.Add(fg_sizer, 1, wx.EXPAND|wx.ALL, 5)
        self.SetSizer(sizer)
        
    #   __init__()
    #----------------------------------------------------------------
    def On_Type_Choice(self, event):
        
        #-----------------------------------------
        #  Event handler for the Type droplists.
        #-----------------------------------------
        print('droplist index  =', event.GetSelection())
        print('droplist string =', event.GetString())
        # print 'event.Id      =', event.Id
        # print 'event.GetId() =', event.GetId()

        #---------------------------------------
        # Need to know the row of the droplist
        #---------------------------------------
        row = (event.Id - 5000)
        # Idea:  row_ID = 'tf_row_' + str(row)
        # row = int(event.Id[7:])  # (ID can't be string ?)
        print('row =', row)
        index = event.GetSelection()
        self.data.types[ row ] = self.type_name[ index ]
        print(' ')
        
    #   On_Type_Choice()
    #----------------------------------------------------------------

#-----------------------------------------------------------------------
class Layer_Data():
    def __init__(self):
        typelist = ['Scalar', 'Time series', \
                    'Grid', 'Grid sequence']
        self.typelist = typelist
        #-----------------------------------------------------------
        self.names  = ['K_s', 'K_i', 'theta_s', 'theta_i']
        self.types  = ['Scalar', 'Scalar', 'Scalar', 'Scalar']
        self.values = ['7.20e-06', '9.849e-08', '0.485', '0.3758']
        self.units  = ['m/s', 'm/s', 'none', 'none']
        self.typelists = [typelist, typelist, typelist, typelist]
        
#-----------------------------------------------------------------------
class TF_Notebook_Test(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, title="TF Notebook Test")
                          ## style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        
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
        ts_box.Add((24,24), 1)
        ts_box.Add(L1)
        ts_box.Add((hgap, hgap), 1)
        ts_box.Add(text)
        ts_box.Add((hgap, hgap), 1)
        ts_box.Add(L2)
        ts_box.Add((vgap, vgap), 1)
        
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
    
