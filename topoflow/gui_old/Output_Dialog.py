
# Add ability to pass TF_Output_Dialog a data structure
# instead of an XML filename.  XML file is only for
# initial creation, but for updates or showing user-
# modified settings, need this other option.

# Read_XML_File() should return True or False depending
# on whether XML file is found and formatted correctly.

#---------------------------------------------------------

#!/usr/bin/env python

#  April 30, 2009
#  S.D. Peckham

import wx
import xml.dom.minidom
import time
import webbrowser  ##  standard Python module
## import wx.html
import glob

from .TF_Output_Var_Box import *   ################

#  Check these out later.
#  import urllib
#  import xmllib (deprecated: use xml.sax instead)
#  import htmllib

#-------------------------------------------------------------

#  class TF_Out_Variable_Info
#  class TF_Out_Variable_List
#  class TF_out_info
#  class TF_Output_Dialog
#        XML_File_Found
#        Out_Variable_Notebook
#        Out_Variable_Wizard    #### (not ready yet)
#        Button_Box
#        On_Type_Choice
#        On_OK
#        On_Help
#        On_Cancel
#        Read_XML_Info

#----------------------------------------------------------------
class TF_Out_Variable_Info():
    def __init__(self, name='', label='', value='', units=''):
        
        self.name   = name
        self.label  = label
        self.value  = value
        self.units  = units
        
#   __init__()
#----------------------------------------------------------------
class TF_Out_Variable_List():
    def __init__(self, n_variables=1):

        self.var_names  = []
        self.var_labels = []
        self.var_values = []
        self.var_units  = []
        self.timestep   = TF_Out_Variable_Info()
        
#   __init__()
#----------------------------------------------------------------
class TF_Out_Info():
    def __init__(self, n_boxes=2, n_variables=1):
        
        #----------------------------------------
        # Each box has a list of TF variables
        #----------------------------------------
        self.n_boxes     = n_boxes
        self.n_variables = n_variables  # (per box)
        self.boxes       = []
        for each in range(n_boxes):
            self.boxes.append( TF_Out_Variable_List() )
            ## self.boxes.append( TF_Out_Variable_List(n_variables=n_variables) )
            
        self.help_file = ''
             
#   __init__()
#-------------------------------------------------------------
class TF_Output_Dialog(wx.Frame):

    #-----------------------------------------------------
    # Notes:  This class is for creating TopoFlow output
    #         dialogs, which all use a similar template.
    #         Initial settings, labels, etc. are read
    #         from the XML file provided
    #         Default is wx.ALIGN_LEFT for labels & text
    #-----------------------------------------------------
    def __init__(self, parent=None, id=-1, \
                 main_frame=None, \
                 xml_file="xml/infil_richards_1d_out_options.xml"):

        #-----------------------
        # Try to find XML file
        #-----------------------
        FOUND = self.XML_File_Found(xml_file)
        if not(FOUND):  return
             
        #--------------------------------------------
        # Set XML file to read info from, including
        # the name of the HTML help file
        #--------------------------------------------
        self.xml_file = xml_file
        self.Read_XML_Info()
            
        #-----------------------------------
        # Get title for this output dialog
        #-----------------------------------
        proc_name   = self.out_info.proc_name
        method_name = self.out_info.method_name
        if (proc_name != "None"):
            title = proc_name
        else:
            title = ""
        if (method_name != "None"):
            title += ": " + method_name
        title += " Output Dialog"
        
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
       
        #-----------------------------------------
        # Create objects to appear in the dialog
        #-----------------------------------------
        if (self.out_info.n_boxes == 1):
            data       = self.out_info.boxes[0]
            var_box    = TF_Output_Var_Box(parent=self.panel,
                                           main_frame=self.main_frame,
                                           data=data)
        else:
            var_box    = self.Out_Variable_Notebook()
            ## var_box = self.Out_Variable_Wizard()
        button_bar = self.Button_Box()
        
        #--------------------------------------------
        # Add objects to the main sizer
        # wx.ALL adds border space on all 4 sides.
        #-------------------------------------------
        border = 2 * self.vgap
        sizer.Add(var_box,    0, wx.ALL, 0)
        sizer.Add(button_bar, 0, wx.ALL, border)

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
    def XML_File_Found(self, xml_file):
        
        files = glob.glob(xml_file)
        if (len(files) == 0):
            msg = "Can't find the XML file:\n\n"
            msg += (xml_file + "\n")
            dialog = wx.MessageBox(msg, caption='SORRY,')
            return False
        else:
            return True
        
    #   XML_File_Found()      
    #----------------------------------------------------------------
    def Out_Variable_Notebook(self):
        
        notebook = wx.Notebook(self.panel, style=wx.BK_TOP)

        k = 0
        n_boxes = self.out_info.n_boxes

        labels = ['Grids', 'Values at Pixels', 'Stacks', \
                  'Z-profiles at Pixels']
        for k in range(n_boxes):
            data  = self.out_info.boxes[k]
            page  = TF_Output_Var_Box(parent=notebook, \
                                      data=data)
            notebook.AddPage(page, labels[k])

        return notebook
    
    #   Out_Variable_Notebook()
    #----------------------------------------------------------------
    def Out_Variable_Wizard(self):
        pass
    
    #   Out_Variable_Wizard()  
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
    # Note:  On_Check_Box() is located in TF_Output_Var_Box.py
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
    
        #------------------------
        # Try to find HTML file
        #------------------------
        help_file = self.out_info.help_file
        files = glob.glob(help_file)
        if (len(files) == 0):
            msg = "Can't find the HTML help file:\n\n"
            msg += (help_file + "\n")
            dialog = wx.MessageBox(msg, caption='SORRY,')
            return
        
        #---------------------------
        # Open the default browser
        #---------------------------
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
        B_elements = dom.firstChild.getElementsByTagName("box")
        V_elements = dom.firstChild.getElementsByTagName("variable")
        T_elements = dom.firstChild.getElementsByTagName("timestep")
        H_elements = dom.firstChild.getElementsByTagName("help_file")
        #-------------------------------------------------------------
        nP = len(P_elements)
        nB = len(B_elements)
        nV = len(V_elements)
        nT = len(T_elements)
        nH = len(H_elements)
    ##    print 'n_boxes               =', nB
    ##    print 'n_variables per layer =', (nV / nB)
        
        #-------------------
        # Issue warnings ?
        #-------------------
        msg = 'TopoFlow Output Dialog XML files \n'
        if (nP > 1):
            msg += 'may only contain 1 "process" tag. \n'
            wx.MessageBox(msg, caption='XML FILE ERROR ')
            return
        if (nB == 0):
            msg += 'must contain 1 or more "box" tags. \n'
            wx.MessageBox(msg, caption='XML FILE ERROR ')
            return
        if (nT == 0) or (nT != nB):
            msg += 'must contain 1 "timestep" tag per box. \n'
            wx.MessageBox(msg, caption='XML FILE ERROR ')
            return
        
        #-----------------------------------------
        # There should be 1 top-level child node
        # that has the tag name "process"
        #-----------------------------------------
        # top_node = dom.childNodes[0]
        top_node = dom.firstChild

        #-------------------------------------------------
        # Extract process name and store it in proc_info
        #-------------------------------------------------
        
        #---------------------------------------------
        # Prepare to save process info from XML file
        #---------------------------------------------
        out_info = TF_Out_Info(n_boxes=nB, n_variables=nV/nB)

        #---------------------------------
        # Read proc_name and method_name
        #---------------------------------
        PN_nodes = top_node.getElementsByTagName("proc_name")
        if (len(PN_nodes) > 0):
            out_info.proc_name = PN_nodes[0].firstChild.data.strip()
        else:
            out_info.proc_name = "None"
        MN_nodes = top_node.getElementsByTagName("method_name")
        if (len(MN_nodes) > 0):
            out_info.method_name = MN_nodes[0].firstChild.data.strip()
        else:
            out_info.method_name = "None"

        #-------------------------------------------------------
        # For each variable in each box, save the associated
        # information.  The process node should have at least
        # one "box" tag.
        #------------------------------------------------------
        box_nodes = top_node.getElementsByTagName("box")
        k = 0
        for B_node in box_nodes:
            var_nodes = B_node.getElementsByTagName("variable")
            j = 0
            for V_node in var_nodes:
                N_nodes = V_node.getElementsByTagName("name")
                name    = N_nodes[0].firstChild.data.strip()
                #---------------------------------------------------
                S_nodes = V_node.getElementsByTagName("label")
                label  = S_nodes[0].firstChild.data.strip()
                #---------------------------------------------------
                values = V_node.getElementsByTagName("value")
                value  = values[0].firstChild.data.strip()
                #---------------------------------------------------
                units  = V_node.getElementsByTagName("units")
                unit   = units[0].firstChild.data.strip()
                #---------------------------------------------------
##                out_info.boxes[k].variables[j].name   = name
##                out_info.boxes[k].variables[j].label  = label
##                out_info.boxes[k].variables[j].value  = value
##                out_info.boxes[k].variables[j].units  = unit
                #-----------------------------------------------------
                out_info.boxes[k].var_names.append( name )
                out_info.boxes[k].var_labels.append( label )
                out_info.boxes[k].var_values.append( value )
                out_info.boxes[k].var_units.append( unit )
                #-----------------------------------------------------
                j += 1

            #-----------------------------------------
            # Read timestep information for this box
            #-----------------------------------------
            T_nodes  = B_node.getElementsByTagName("timestep")
            TL_nodes = T_nodes[0].getElementsByTagName("label")
            TV_nodes = T_nodes[0].getElementsByTagName("value")
            TU_nodes = T_nodes[0].getElementsByTagName("units")
            #----------------------------------------------------------------
            out_info.boxes[k].timestep.label = TL_nodes[0].firstChild.data
            out_info.boxes[k].timestep.value = TV_nodes[0].firstChild.data
            out_info.boxes[k].timestep.units = TU_nodes[0].firstChild.data
            #------------------
            k += 1
        
        #------------------------------
        # Read name of HTML help file
        #------------------------------
        H_nodes = top_node.getElementsByTagName("help_file")
        if (len(H_nodes) > 0):
            out_info.help_file = H_nodes[0].firstChild.data.strip()     

        self.out_info = out_info

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
##    frame = TF_Output_Dialog(parent=None, id=-1, \
##                 xml_file="xml/snowmelt_output_choices.xml")
    #---------------------------------
    # Example with 3 layers and tabs
    #---------------------------------
    frame = TF_Output_Dialog(parent=None, id=-1, \
                 xml_file="xml/infil_richards_1D_out_options.xml")
    app.MainLoop()

    
