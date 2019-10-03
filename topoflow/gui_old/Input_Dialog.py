
# Add ability to pass TF_Input_Dialog a data structure
# instead of an XML filename.  XML file is only for
# initial creation, but for updates or showing user-
# modified settings, need this other option.
#---------------------------------------------------------

#!/usr/bin/env python

#  August 8 & 11, 2008
#  February 10, 2009
#  April 21-28, 2009
#  S.D. Peckham

import wx
import xml.dom.minidom
import time
import webbrowser  ##  standard Python module
## import wx.html
import glob

from .TF_Input_Var_Box import *   ################

#  Check these out later.
#  import urllib
#  import xmllib (deprecated: use xml.sax instead)
#  import htmllib

#-------------------------------------------------------------

#  class TF_In_Variable_Info
#  class TF_In_Variable_List
#  class TF_Process_Info
#  class TF_Input_Dialog
#        In_Variable_Notebook
#        In_Variable_Wizard    #### (not ready yet)
#        Timestep_Box
#        Button_Box
#        On_Type_Choice
#        On_OK
#        On_Help
#        On_Cancel
#        Read_XML_Info

#----------------------------------------------------------------
class TF_In_Variable_Info():
    def __init__(self, name='', label='', vtype='Scalar',
                 value='', units='', type_choices=None):
        
        self.name   = name
        self.label  = label
        self.type   = vtype
        self.value  = value
        self.units  = units

        if (type_choices == None):
            self.type_choices = ['Scalar', 'Time series', \
                                  'Grid', 'Grid sequence']
        
#   __init__()
#----------------------------------------------------------------
class TF_In_Variable_List():
    def __init__(self, n_variables=1):

        #-----------------
        # First approach
        #-----------------
        ## self.variables = []
        ## for each in range(n_variables):
        ##     self.variables.append( TF_Variable_Info() )

        #------------------------------------------
        # Alternate approach with some advantages
        #------------------------------------------
        self.var_names        = []
        self.var_labels       = []
        self.var_types        = []
        self.var_values       = []
        self.var_units        = []
        self.var_type_choices = []
        
#   __init__()
#----------------------------------------------------------------
class TF_Process_Info():
    def __init__(self, name='Infiltration',
                 n_layers=1, n_variables=1):
        #----------------------------------------
        # Each layer has a list of TF variables
        #----------------------------------------
        self.n_layers    = n_layers
        self.n_variables = n_variables  # (per layer)
        self.layers      = []
        for each in range(n_layers): 
            self.layers.append( TF_In_Variable_List(n_variables=n_variables) )
            
        self.timestep  = TF_In_Variable_Info()
        self.help_file = ''
             
#   __init__()
#-------------------------------------------------------------
class TF_Input_Dialog(wx.Frame):

    #-----------------------------------------------------
    # Notes:  This class is for creating TopoFlow input
    #         dialogs, which all use a similar template.
    #         Initial settings, labels, etc. are read
    #         from the XML file provided
    #         Default is wx.ALIGN_LEFT for labels & text
    #-----------------------------------------------------
    def __init__(self, parent=None, id=-1, 
                 main_frame=None,
                 xml_file="infil_richards_1D.xml"):

        #-----------------------
        # Try to find XML file
        #-----------------------
        files = glob.glob(xml_file)
        if (len(files) == 0):
            msg = "Can't find the XML file:\n\n"
            msg += (xml_file + "\n")
            dialog = wx.MessageBox(msg, caption='SORRY,')
            return
        
        #--------------------------------------------
        # Set XML file to read info from, including
        # the name of the HTML help file
        #--------------------------------------------
        self.xml_file = xml_file
        self.Read_XML_Info()

        #----------------------------------
        # Get title for this input dialog
        #----------------------------------
        title = self.proc_info.proc_name + " : " + \
                self.proc_info.method_name + \
                " Input Dialog"
        
        #-------------------------------------------
        # Initialize a wxPython frame, add a panel
        #-------------------------------------------
        wx.Frame.__init__(self, parent, id, title)
        panel = wx.Panel(self, -1)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetBackgroundColour('Light Blue')
        self.panel = panel
        self.main_sizer = sizer
  
        #--------------------------------------------------
        # Saving main_frame allows collected values to be
        # stored in its state before this one closes.
        #--------------------------------------------------
        self.main_frame = main_frame
        self.parent     = parent
        self.title      = title
        self.vgap       = 10
        self.hgap       = 6

        #-------------------------------------------------
        # Later move these into Variable_Box() or remove
        #-------------------------------------------------
        self.type_code = {'Scalar':0, 'Time series':1, \
                          'Grid':2, 'Grid sequence':3}
        self.type_name = {0:'Scalar', 1:'Time series', \
                          2:'Grid', 3:'Grid sequence'}
       
        #-----------------------------------------
        # Create objects to appear in the dialog
        #-----------------------------------------
        if (self.proc_info.n_layers == 1):
            data    = self.proc_info.layers[0]
            var_box = TF_Input_Var_Box(parent=self.panel,
                                       main_frame=self.main_frame,
                                       data=data)
        else:
            var_box = self.In_Variable_Notebook()
            ## var_box = self.In_Variable_Wizard()

        ADD_TIME_BOX = (self.proc_info.timestep.value != '')
        if (ADD_TIME_BOX):
            time_box   = self.Timestep_Box()
        button_bar = self.Button_Box()
        
        #--------------------------------
        # Add objects to the main sizer
        #--------------------------------
        border1 = self.vgap
        border2 = 2 * self.vgap
        sizer.Add(var_box,    0, wx.ALL, border1)
        if (ADD_TIME_BOX):
            sizer.Add(time_box,   0, wx.ALL, border1)
        sizer.Add(button_bar, 0, wx.ALL, border2)

        panel.SetSizer(sizer)
        sizer.Fit(self)
        self.Show()        # (need here if called)
        self.Centre()
        ## sizer.SetSizeHints(self)   # (not needed)
        
        #--------------------------------------------
        # Doesn't work for some reason (see p. 328)
        #--------------------------------------------
        # self.SetSizer(sizer)
        # self.Fit()
        
    #   __init__()
    #----------------------------------------------------------------
    def In_Variable_Notebook(self):
        
        notebook = wx.Notebook(self.panel, style=wx.BK_TOP)

        k = 0
        n_layers = self.proc_info.n_layers

        for k in range(n_layers):
            data  = self.proc_info.layers[k]
            kstr  = str(k+1)
            label = "Layer " + kstr + " variables"
            page  = TF_Input_Var_Box(parent=notebook, \
                                     data=data, \
                                     box_label=label)
            notebook.AddPage(page, "Layer " + kstr)

        return notebook
    
    #   In_Variable_Notebook()
    #----------------------------------------------------------------
    def In_Variable_Wizard(self):
        pass
    
    #   In_Variable_Wizard()  
    #----------------------------------------------------------------
    def Timestep_Box(self):

        #---------------------------------------------
        #  Create sizer box for the process timestep
        #---------------------------------------------
        time_info = self.proc_info.timestep
        unit_str  = "[" + time_info.units + "]"
        ## unit_str  = "[" + time_info.units + " / timestep]"
        L1    = wx.StaticText(self.panel, -1, time_info.label + ":")
        text  = wx.TextCtrl(self.panel,   -1, time_info.value)
        L2    = wx.StaticText(self.panel, -1, unit_str)
        #-------------------------------------------------------        
        box   = wx.BoxSizer(wx.HORIZONTAL)
        box.Add((3*self.hgap, 3*self.hgap), 1)
        box.Add(L1)
        box.Add((self.hgap, self.hgap), 1)
        box.Add(text)
        box.Add((self.hgap, self.hgap), 1)
        box.Add(L2)

        # self.SetSizer(box)  ## (Bad to do this here.)
        return box
    
    #   Timestep_Box()
    #----------------------------------------------------------------
    def Button_Box(self):

        #----------------------------
        #  Create bottom button bar
        #----------------------------
        okay_btn   = wx.Button(self.panel, -1, "OK")
        help_btn   = wx.Button(self.panel, -1, "Help")
        cancel_btn = wx.Button(self.panel, -1, "Cancel")
        #-------------------------------------------------        
        box = wx.BoxSizer(wx.HORIZONTAL)
        proportion = 0
        border = 5
        box.Add((20,20),    proportion, border)
        box.Add(okay_btn,   proportion, border)
        box.Add((10,10),    proportion, border)
        box.Add(help_btn,   proportion, border)
        box.Add((10,10),    proportion, border)
        box.Add(cancel_btn, proportion, border)
        box.Add((10,10),    proportion, border)        
        # box.AddMany([okay_btn, help_btn, cancel_btn])
        #-----------------------------------------------------
        self.Bind(wx.EVT_BUTTON, self.On_OK,     okay_btn)
        self.Bind(wx.EVT_BUTTON, self.On_Help,   help_btn)
        self.Bind(wx.EVT_BUTTON, self.On_Cancel, cancel_btn)

        return box

    #   Button_Box()
    #----------------------------------------------------------------
    # Note:  On_Type_Choice() is located in TF_Input_Var_Box.py
    #----------------------------------------------------------------
    #----------------------------------------------------------------
    def On_OK(self, event):

        ################################################
        #  Need to now include "layer" as well when
        #  saving values from text boxes to parent.
        #
        #  Make sure that timestep gets saved also.
        ################################################

        if (self.parent == None):
            print('Sorry:  This dialog has no parent dialog,')
            print('so collected values cannot be saved.')
            print(' ')
            return
        
        #---------------------------------------------
        #  Read values from text boxes and save them
        #  in parent's state with droplist settings
        #---------------------------------------------
        #  Can use validator with TextCtrl to check
        #  that values entered are valid (p. 282)
        #---------------------------------------------
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
                label = self.var_labels[k]
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

    #   On_OK()
    #----------------------------------------------------------------
    def On_Help(self, event):
  
        #---------------------------
        # Open the default browser
        #---------------------------
        help_file = self.proc_info.help_file
        result = webbrowser.open('file://' + help_file)
        
    #   On_Help()
    #----------------------------------------------------------------
    def On_Cancel(self, event):

        self.Destroy()

    #   On_Cancel()
    #----------------------------------------------------------------
    def Read_XML_Info(self):

        #---------------------------------------------------------
        # Notes:  It is helpful to have a look at the DOM specs,
        #         which can be found online at:
        #         http://www.w3.org/TR/1998/REC-DOM-Level-1-
        #                19981001/introduction.html
        #                (see the tree diagram)
        #         http://www.w3.org/TR/1998/REC-DOM-Level-1-
        #                19981001/level-one-core.html
        #                (search for "firstChild")
        #---------------------------------------------------------
        
        #--------------------------------------------
        # Read input variable info from an XML file
        #--------------------------------------------
        file_obj = open(self.xml_file, 'r')
        #  Read entire file into a string
        doc_string = file_obj.read()
        dom = xml.dom.minidom.parseString( doc_string )

        #----------------------------------------------
        # Count all tags in XML file of various types
        #----------------------------------------------
        P_elements = dom.firstChild.getElementsByTagName("process") 
        L_elements = dom.firstChild.getElementsByTagName("layer")
        V_elements = dom.firstChild.getElementsByTagName("variable")
        T_elements = dom.firstChild.getElementsByTagName("timestep")
        D_elements = dom.firstChild.getElementsByTagName("droplist")
        H_elements = dom.firstChild.getElementsByTagName("help_file")
        #-------------------------------------------------------------
        nP = len(P_elements)
        nL = len(L_elements)
        nV = len(V_elements)
        nT = len(T_elements)
        nD = len(D_elements)
        nH = len(H_elements)
    ##    print 'n_layers              =', nL
    ##    print 'n_variables per layer =', (nV / nL)
        
        #-------------------
        # Issue warnings ?
        #-------------------
        msg = 'TopoFlow Input Dialog XML files \n'
        if (nP > 1):
            msg += 'may only contain 1 "process" tag. \n'
            wx.MessageBox(msg, caption='ERROR ')
            return
        if (nL == 0):
            msg += 'must contain 1 or more "layer" tags. \n'
            wx.MessageBox(msg, caption='ERROR ')
            return
##        if (nT > 1) or (nT == 0):
##            msg += 'must contain 1 "timestep" tag. \n'
##            wx.MessageBox(msg, caption='ERROR ')
##            return
        
        #-----------------------------------------
        # There should be 1 top-level child node
        # that has the tag name "process"
        #-----------------------------------------
        # process_node = dom.childNodes[0]
        process_node = dom.firstChild

        #-------------------------------------------------
        # Extract process name and store it in proc_info
        #-------------------------------------------------
        
        #---------------------------------------------
        # Prepare to save process info from XML file
        #---------------------------------------------
        proc_info = TF_Process_Info(n_layers=nL, n_variables=nV/nL)

        #---------------------------------
        # Read proc_name and method_name
        #---------------------------------
        PN_nodes = process_node.getElementsByTagName("proc_name")
        if (len(PN_nodes) > 0):
            proc_info.proc_name = PN_nodes[0].firstChild.data.strip()
        else:
            proc_info.proc_name = "None"
        MN_nodes = process_node.getElementsByTagName("method_name")
        if (len(MN_nodes) > 0):
            proc_info.method_name = MN_nodes[0].firstChild.data.strip()
        else:
            proc_info.method_name = "None"

        #-------------------------------------------------------
        # For each variable in each layer, save the associated
        # information.  The process node should have at least
        # one layer tag.
        #------------------------------------------------------
        layer_nodes = process_node.getElementsByTagName("layer")
        k = 0
        for L_node in layer_nodes:
            var_nodes = L_node.getElementsByTagName("variable")
            j = 0
            for V_node in var_nodes:
                N_nodes = V_node.getElementsByTagName("name")
                name    = N_nodes[0].firstChild.data.strip()
                #---------------------------------------------------
                S_nodes = V_node.getElementsByTagName("label")
                label  = S_nodes[0].firstChild.data.strip()
                #---------------------------------------------------
                vtypes = V_node.getElementsByTagName("type")
                vtype  = vtypes[0].firstChild.data.strip()
                #---------------------------------------------------
                values = V_node.getElementsByTagName("value")
                value  = values[0].firstChild.data.strip()
                #---------------------------------------------------
                units  = V_node.getElementsByTagName("units")
                unit   = units[0].firstChild.data.strip()
                #---------------------------------------------------
                tlists = V_node.getElementsByTagName("type_choices")
                tlist  = tlists[0].firstChild.data.strip().split(",")
                #---------------------------------------------------
##                proc_info.layers[k].variables[j].name   = name
##                proc_info.layers[k].variables[j].label = label
##                proc_info.layers[k].variables[j].type   = vtype
##                proc_info.layers[k].variables[j].value  = value
##                proc_info.layers[k].variables[j].units  = unit
##                proc_info.layers[k].variables[j].type_choices = tlist
                #-----------------------------------------------------
                proc_info.layers[k].var_names.append( name )
                proc_info.layers[k].var_labels.append( label )
                proc_info.layers[k].var_types.append( vtype )
                proc_info.layers[k].var_values.append( value )
                proc_info.layers[k].var_units.append( unit )
                proc_info.layers[k].var_type_choices.append( tlist )
                #-----------------------------------------------------
                j += 1
            k += 1
        
        #----------------------------
        # Read timestep information
        #----------------------------
        if (nT > 0):
            T_nodes = process_node.getElementsByTagName("timestep")
            L_nodes = T_nodes[0].getElementsByTagName("label")
            V_nodes = T_nodes[0].getElementsByTagName("value")
            U_nodes = T_nodes[0].getElementsByTagName("units")
            #-------------------------------------------------------
            proc_info.timestep.label = L_nodes[0].firstChild.data
            proc_info.timestep.value = V_nodes[0].firstChild.data
            proc_info.timestep.units = U_nodes[0].firstChild.data

        #------------------------------
        # Read name of HTML help file
        #------------------------------
        H_nodes = process_node.getElementsByTagName("help_file")
        if (len(H_nodes) > 0):
            proc_info.help_file = H_nodes[0].firstChild.data.strip()     

        self.proc_info = proc_info

    #   Read_XML_Info()
    #----------------------------------------------------------------
        
#---------------------------------------
#  Support two different usage options
#---------------------------------------
if (__name__ == '__main__'):
    app = wx.PySimpleApp()
    #-------------------------------------
    # Example with one layer and no tabs
    #-------------------------------------
##    frame = TF_Input_Dialog(parent=None, id=-1, \
##                 xml_file="xml/snowmelt_degree_day.xml")
    #---------------------------------
    # Example with 3 layers and tabs
    #---------------------------------
    frame = TF_Input_Dialog(parent=None, id=-1, \
                 xml_file="xml/infil_richards_1D.xml")
    app.MainLoop()

    
