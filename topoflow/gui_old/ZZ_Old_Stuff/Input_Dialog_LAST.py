#!/usr/bin/env python

#  August 8 & 11, 2008
#  February 10, 2009
#  April 21, 2009
#  S.D. Peckham

##  This needs to set TF home directory, etc.  #################

import wx
import wx.html
import xml.dom.minidom
import time
import webbrowser  ##  standard Python module

#  Check these out later.
#  import urllib
#  import xmllib (deprecated: use xml.sax instead)
#  import htmllib

#-------------------------------------------------------------
class TF_Input_Dialog(wx.Frame):

    #-----------------------------------------------------
    # Notes:  This class is for creating TopoFlow input
    #         dialogs, which all use a similar template.
    #         Initial settings, labels, etc. are read
    #         from the XML file provided
    #         Default is wx.ALIGN_LEFT for labels & text
    #-----------------------------------------------------
        
    #---------------------------
    #  Create top-level dialog
    #---------------------------
    def __init__(self, parent=None, id=-1, \
                 xml_file="xml/snowmelt_degree_day.xml", \
                 title="Snowmelt: Degree-Day Input Dialog"):

        #-------------------------------------------
        # Initialize a wxPython frame, add a panel
        #-------------------------------------------
        wx.Frame.__init__(self, parent, id, title)
        panel = wx.Panel(self, -1)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetBackgroundColour('Light Blue')
        self.panel = panel
        self.main_sizer = sizer
        # self.panel.SetBackgroundColour('White')
        # self.panel.SetForegroundColour('White')
        
        #------------------------------------------------
        # Saving parent allows collected values to be
        # stored in parent frame before this one closes.
        #------------------------------------------------
        self.parent = parent
        self.title  = title
        self.vgap   = 10
        self.hgap   = 6
        self.type_code = {'Scalar':0, 'Time series':1, \
                          'Grid':2, 'Grid sequence':3}
        self.type_name = {0:'Scalar', 1:'Time series', \
                          2:'Grid', 3:'Grid sequence'}

        #--------------------------------------------
        # Set XML file to read info from, including
        # the name of the HTML help file
        #--------------------------------------------
        self.xml_file = xml_file
        self.read_var_info()
        
        #-----------------------------------------
        # Create objects to appear in the dialog
        #-----------------------------------------
        var_box    = self.variable_box()  # (returns a sizer)
        pad_row1   = wx.StaticLine(panel)
        time_box   = self.timestep_box()
        pad_row2   = wx.StaticLine(panel)
        button_bar = self.button_bar()
        pad_row3   = wx.StaticLine(panel)

        self.var_box = var_box  ##########  EXPERIMENT
        
        #--------------------------------
        # Add objects to the main sizer
        #--------------------------------
        vpad = 0
        sizer.Add(var_box,    0, wx.ALL, self.vgap)
        sizer.Add(pad_row1,   0, wx.ALL, self.vgap)
        sizer.Add(time_box,   0, wx.ALL, vpad)
        sizer.Add(pad_row2,   0, wx.ALL, self.vgap)
        sizer.Add(button_bar, 0, wx.ALL, vpad)
        sizer.Add(pad_row3,   0, wx.ALL, self.vgap)

        panel.SetSizer(sizer)
        sizer.Fit(self)
        ## sizer.SetSizeHints(self)   # (not needed)
        
        #--------------------------------------------
        # Doesn't work for some reason (see p. 328)
        #--------------------------------------------
        # self.SetSizer(sizer)
        # self.Fit()
        
    #   __init__()
    #----------------------------------------------------------------
    def get_xml_tag_data(self, var, tag_name):

        #-------------------------------------------
        #  Get data string for an XML variable tag
        #-------------------------------------------
        vstr = var.getElementsByTagName(tag_name)[0].firstChild.data
        return vstr.strip()

    #   get_xml_tag_data()    
    #----------------------------------------------------------------
    def read_var_info(self):

        #--------------------------------------------------
        #  Read descriptions of input variables from file
        #--------------------------------------------------
        self.var_names    = []
        self.var_units    = []
        self.var_values   = []
        self.var_types    = []
        self.var_typelist = []
        
        #---------------------------------------------
        #  Read variable info from an XML file
        #---------------------------------------------
        file_obj = open(self.xml_file, 'r')
        #  Read entire file into a string
        doc_string = file_obj.read()
        dom = xml.dom.minidom.parseString( doc_string )
        variables = dom.firstChild.getElementsByTagName("variable")
        # print "variables   =", variables
        
        for var in variables:
            symbol   = self.get_xml_tag_data(var, "symbol")
            units    = self.get_xml_tag_data(var, "units")
            value    = self.get_xml_tag_data(var, "value")
            vtype    = self.get_xml_tag_data(var, "type")
            typelist = self.get_xml_tag_data(var, "typelist")
            
            self.var_names.append( symbol )
            self.var_values.append( value )
            self.var_types.append( vtype )
            typelist = typelist.split(",")
            self.var_typelist.append( typelist )
            self.var_units.append( units )

        self.n_vars = len(self.var_names)

##        This doesn't work 
##        symbol = variables.getElementsByTagName("symbol")[0].firstChild.data
##        for symbol in symbols:
##            print "symbol =", symbol
        
        #-----------------------------
        #  Read timestep information
        #-----------------------------
        timesteps = dom.firstChild.getElementsByTagName("timestep")
        if (len(timesteps) > 0):
            var = timesteps[0]
            self.timestep_label = self.get_xml_tag_data(var, "label")
            self.timestep_value = self.get_xml_tag_data(var, "value")
            self.timestep_units = self.get_xml_tag_data(var, "units")
            
        #-------------------------------
        #  Read name of HTML help file
        #-------------------------------
        help_files = dom.firstChild.getElementsByTagName("help_file")
        if (len(help_files) > 0):
            self.help_file = help_files[0].firstChild.data.strip()
            # print "help file =", self.help_file      
        
    #   read_var_info()
    #----------------------------------------------------------------
    def variable_box(self):

        #-------------------------------------------
        #  Create sizer box for all input var info
        #  (This provides a frame and title.)
        #-------------------------------------------
        vbox  = wx.StaticBox(self.panel, -1, "Input variables:")
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
        #fg_sizer = wx.FlexGridSizer(cols=nh, hgap=self.hgap, vgap=self.vgap)
        
        #----------------------------------------------       
        #  Specify which columns can expand on resize
        #----------------------------------------------
        # fg_sizer.AddGrowableCol(1)
        self.text_boxes = []

        #-------------------------
        # Add the column headers
        #-------------------------
        for row in range(nh):
            L1 = wx.StaticText(self.panel, -1, header[row])
            fg_sizer.Add(L1, 0, wx.ALL, self.hgap)

        #------------------------------------------------
        #  Create a row in the dialog for each variable
        #------------------------------------------------
        for row in range(self.n_vars):
            vstr  = self.var_names[row] + ": "
            label = wx.StaticText(self.panel, -1, vstr)
            #----------------------------------------------------------
            row_ID = (5000 + row)   #####
            ### row_ID = 'tf_row_' + str(row)  # (ID can't be string.)
            dlist = wx.Choice(self.panel, row_ID, choices=self.var_typelist[row])
            dlist.Select(self.type_code[self.var_types[row]])
            self.Bind(wx.EVT_CHOICE, self.on_Type_Choice, dlist)      
            #----------------------------------------------------------            
            text  = wx.TextCtrl(self.panel, -1, self.var_values[row],
                                size=(160,-1))
            self.text_boxes.append( text )  #######
            #----------------------------------------------------------
            ustr  = wx.StaticText(self.panel, -1, self.var_units[row])
            #----------------------------------------------------------
            fg_sizer.Add(label, 0, wx.ALL, self.hgap)
            fg_sizer.Add(dlist, 0, wx.ALL, self.hgap)
            fg_sizer.Add(text,  0, wx.ALL, self.hgap)
            fg_sizer.Add(ustr,  0, wx.ALL, self.hgap)

        #---------------------------------
        # Add fg_sizer to the main sizer
        #---------------------------------
        sizer.Add(fg_sizer, 0, wx.ALL, self.vgap)   ######
        return sizer

        #-----------------------------------------
        #  Leave this to "main" to avoid trouble
        #-----------------------------------------
        # self.panel.SetSizer(sizer)
        # sizer.Fit(self)
        # sizer.SetSizeHints(self)
        
    #   variable_box()   
    #----------------------------------------------------------------
    def timestep_box(self):

        #---------------------------------------------
        #  Create sizer box for the process timestep
        #---------------------------------------------
        L1    = wx.StaticText(self.panel, -1, self.timestep_label)
        text  = wx.TextCtrl(self.panel,   -1, self.timestep_value)
        L2    = wx.StaticText(self.panel, -1, self.timestep_units)
        #------------------------------------------------------------        
        box   = wx.BoxSizer(wx.HORIZONTAL)
        box.Add((24,24), 1)
        box.Add(L1)
        box.Add((self.hgap, self.hgap), 1)
        box.Add(text)
        box.Add((self.hgap, self.hgap), 1)
        box.Add(L2)
        box.Add((self.vgap, self.vgap), 1)
        # box.AddMany([L1, text, L2]) #, 1, wx.EXPAND)

        # self.SetSizer(box)  ## (Bad to do this here.)
        return box
    
    #   timestep_box()
    #----------------------------------------------------------------
    def button_bar(self):

        #----------------------------
        #  Create bottom button bar
        #----------------------------
        start_btn  = wx.Button(self.panel, -1, "Start")
        help_btn   = wx.Button(self.panel, -1, "Help")
        cancel_btn = wx.Button(self.panel, -1, "Cancel")
        #-------------------------------------------------        
        box = wx.BoxSizer(wx.HORIZONTAL)
        box.Add((20,20), 1)
        box.Add(start_btn)
        box.Add((10,10), 1)
        box.Add(help_btn)
        box.Add((10,10), 1)
        box.Add(cancel_btn)
        box.Add((10,10), 1)        
        # box.AddMany([start_btn, help_btn, cancel_btn])
        #-----------------------------------------------------
        self.Bind(wx.EVT_BUTTON, self.on_Start,  start_btn)
        self.Bind(wx.EVT_BUTTON, self.on_Help,   help_btn)
        self.Bind(wx.EVT_BUTTON, self.on_Cancel, cancel_btn)

        # self.SetSizer(box)  ## (Bad to do this here.)
        return box

    #   button_bar()
    #----------------------------------------------------------------
    def on_Type_Choice(self, event):
        
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
        self.var_types[ row ] = self.type_name[ index ]
        print(' ')
        
    #   on_Type_Choice()       
    #----------------------------------------------------------------
    def on_Start(self, event):

        ################################################
        #  Make sure that timestep gets saved also.
        ################################################
        
        #-----------------------------------------
        #  Event handler for the Start button.
        #  Read values from text boxes and save
        #  them somewhere with droplist settings
        #--------------------------------------------
        #  Can use validator with TextCtrl to check
        #  that values entered are valid (p. 282)
        #--------------------------------------------
        k = 0
        for box in self.text_boxes:
            val_string = box.GetValue()
            if (len(val_string) == 0):
                wx.MessageBox("Each text box must contain a number or a filename.", "Sorry,")
                box.SetBackgroundColour("pink")
                time.sleep(0.5)
                box.SetFocus()
                box.Refresh()
                return
            else:
                #-----------------------------------------------
                #  Save values from GUI (numerical or string)
                #  into the state of the dialog's parent frame.
                #  (use var_types, var_units, var_values)
                #-----------------------------------------------
                self.var_values[k] = val_string  # (update self)
                print(' ')
                print('Saving values into parent...')
                value  = val_string
                name   = self.var_names[k]
                tcode  = self.type_code[ self.var_types[k] ]
                field_str = 'self.parent.' + name
                print("value =", value)
                
                exec(field_str + "_type = " + str(tcode))
                if (tcode == 0):
                    exec(field_str + '= ' + value)  # (scalar)
                    exec(field_str + "_file = ''")
                else:
                    exec(field_str + '= None')
                    exec(field_str + '_file = "' + value + '"')
                k += 1

        #-----------------------------------------
        # At this point, values in self are lost
        # but values in parent TLB are retained.
        #-----------------------------------------
        self.Destroy()

    #   on_Start()
    #----------------------------------------------------------------
    def on_Help(self, event):

        #------------------------------------
        # EXPERIMENT (4/22/09)  This works.
        #------------------------------------
##        self.main_sizer.Hide(self.var_box)
##        self.main_sizer.Layout()
##        return
    
        #------------------------------------------
        #  Event handler for the Help button.
        #  Open the default browser (e.g. Safari)
        #------------------------------------------
        result = webbrowser.open('file://' + self.help_file)

        # For testing only
##        print 'For testing:  Some values saved in parent.'
##        print 'self.c0      =', self.c0
##        print 'self.c0_type =', self.c0_type
##        print 'self.c0_file =', self.c0_file
        
    #   on_Help()
    #----------------------------------------------------------------
    def on_Help2(self, event):
        
        #-----------------------------------------------
        #  Event handler for the Help button.
        #  Alternate method that uses HTML_Help_Window
        #  class, defined below.
        #-----------------------------------------------
        app   = wx.PySimpleApp()
        frame = HTML_Help_Window(parent=None,
                                 title="HTML Help System",
                                 html_file=self.help_file)
        frame.Show()
        app.MainLoop()

    #   on_Help2()       
    #----------------------------------------------------------------
    def on_Cancel(self, event):

        #----------------------------------------
        #  Event handler for the Cancel button.
        #----------------------------------------
        self.Destroy()

    #   on_Cancel()
    #----------------------------------------------------------------
        
#-----------------------------------------
#  Class for displaying HTML help
#  (now using webbrowser module instead).
#-----------------------------------------        
##class HTML_Help_Window(wx.Frame):
##    def __init__(self, parent, title, html_file):
##        wx.Frame.__init__(self, parent, -1, title,
##                          size=(700,800), pos=(600,50))
##        html = wx.html.HtmlWindow(self)
##        if "gtk2" in wx.PlatformInfo:
##            html.SetStandardFonts()
##
##        html.LoadPage(html_file)

#-------------------------------------------------------------
def Get_Input_Vars():

    #------------------------------------------------------
    # Note:  This is for testing.  It creates two dialogs
    #        that are identical, on top of each other.
    #        Change values in the top one (the child) and
    #        then click on its Start button.  The child
    #        dialog closes.  Now click on the Help button
    #        of the remaining (parent) to print out some
    #        of the stored values.  Later on, the parent
    #        will be the TopoFlow main wizard.
    #------------------------------------------------------
    # Note:  You need to comment out the last few lines
    #        in this file before using this function.
    #------------------------------------------------------
    
    #----------------------------------
    # Open a TopoFlow input dialog as
    # the "main program window"
    #----------------------------------
    app = wx.PySimpleApp()
    frame1 = TF_Input_Dialog(xml_file="xml/snowmelt_degree_day.xml")
    frame1.Show()

    #--------------------------------------------------
    # Open a 2nd dialog that has frame1 as its parent
    # and which will store its values into its parent
    # before it is destroyed.
    #--------------------------------------------------
    frame2 = TF_Input_Dialog(parent=frame1, \
                             xml_file="xml/snowmelt_degree_day.xml")
    frame2.Show()

    #-----------------------
    # Start the event loop
    #-----------------------
    app.MainLoop()
    
#   Get_Input_Vars()
#-------------------------------------------------------------    
        
#---------------------------------------
#  Support two different usage options
#---------------------------------------
if (__name__ == '__main__'):
    app = wx.PySimpleApp()
    # frame = TF_Input_Dialog(parent=None, id=-1)
    frame = TF_Input_Dialog(parent=None, id=-1, \
                 xml_file="xml/infil_richards_1D.xml", \
                 title="Infiltration: Richards 1D Input Dialog")
    frame.Show()
    app.MainLoop()

    
