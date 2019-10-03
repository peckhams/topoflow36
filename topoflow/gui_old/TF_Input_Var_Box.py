#!/usr/bin/env python

#  April 27, 2009
#  S.D. Peckham

import wx

#-------------------------------------------------------------
#   class TF_Input_Var_Box
#         __init__
#         On_Type_Choice()

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
                 main_frame=None, \
                 box_label="Input variables:", data=None):
 
        wx.Panel.__init__(self, parent)
        self.SetBackgroundColour('Light Blue')
        
        #------------------------------------------------
        # Saving parent allows collected values to be
        # stored in parent frame before this one closes.
        #------------------------------------------------
        self.main_frame = main_frame
        self.parent     = parent
        self.data       = data       ###########
        #-----------------------------
        self.vgap       = 10
        self.hgap       = 6
        self.text_width = 240
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

        ####################################        
        parent.text_boxes = []
        ####################################

        #-------------------------
        # Add the column headers
        #-------------------------
        for row in range(nh):
            L1 = wx.StaticText(self, -1, header[row])
            fg_sizer.Add(L1, 0, wx.ALL, self.hgap)

        if (data == None):
            msg = 'No data passed to TF_Input_Var_Box!'
            wx.MesageBox(msg, caption='SORRY,')
            return  #################################
        
        #-----------------------------------------------
        # Create a row in the dialog for each variable
        # using labels, types, values and units.
        #-----------------------------------------------
        for row in range(len(data.var_names)):
            vstr  = data.var_labels[row] + ": "
            label = wx.StaticText(self, -1, vstr)
            #----------------------------------------------------------
            row_ID = (5000 + row)   #####
            ### row_ID = 'tf_row_' + str(row)  # (ID can't be string.)
            dlist = wx.Choice(self, row_ID, choices=data.var_type_choices[row])
            dlist.Select(self.type_code[data.var_types[row]])
            self.Bind(wx.EVT_CHOICE, self.On_Type_Choice, dlist)      
            #----------------------------------------------------------            
            text  = wx.TextCtrl(self, -1, data.var_values[row],
                                size=(self.text_width,-1))
            ####################################
            parent.text_boxes.append( text )
            ####################################
            #----------------------------------------------------------
            units_str = "[" + data.var_units[row] + "]"
            # units_str = data.var_units[row]
            ustr  = wx.StaticText(self, -1, units_str)
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

        #########################################################
        # We need the state of the "top-most" parent here !!
        #########################################################
        if (self.main_frame == None):
            print('Sorry:  The TopoFlow main dialog is missing,')
            print('so droplist choices cannot be saved.')
            print(' ')
            return
        
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
        print(' ')

        #----------------------------------
        # Change choice in internal state
        #----------------------------------
        index = event.GetSelection()
        self.data.var_types[ row ] = self.type_name[ index ]

        #--------------------------------
        # Save choice in parent's state
        #--------------------------------
        # self.main_frame....
        
    #   On_Type_Choice()
    #----------------------------------------------------------------

        
