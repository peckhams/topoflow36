
#!/usr/bin/env python

#  April 28, 2009
#  S.D. Peckham

import wx
import xml.dom.minidom

#   class TF_Variable_Info
#   class TF_Variable_List
#   class TF_Process_Info
#   def   Read_XML_Info

#----------------------------------------------------------------
class TF_Variable_Info():
    def __init__(self, name='', label='', vtype='Scalar',
                 value='', units='', type_choices=None):
        
        self.name   = name
        self.label = label
        self.type   = vtype
        self.value  = value
        self.units  = units

        if (type_choices == None):
            self.type_choices = ['Scalar', 'Time series', \
                                  'Grid', 'Grid sequence']
        
#   __init__()
#----------------------------------------------------------------
class TF_Variable_List():
    def __init__(self, n_variables=1):

        self.variables = []
        for each in range(n_variables):
            self.variables.append( TF_Variable_Info() )

        #------------------------------------------
        # Alternate approach with some advantages
        #------------------------------------------
##        self.var_names         = []
##        self.var_labels       = []
##        self.var_types         = []
##        self.var_values        = []
##        self.var_units         = []
##        self.var_type_choices = []
        
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
            self.layers.append( TF_Variable_List(n_variables=n_variables) )
            
        self.timestep  = TF_Variable_Info()
        self.help_file = ''
             
#   __init__()
#----------------------------------------------------------------
def Read_XML_Info(xml_file="xml/infil_green_ampt.xml"):

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
    file_obj = open(xml_file, 'r')
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
    TF_XML_str = 'ERROR: TopoFlow XML files '
    if (nP > 1):
        print(TF_XML_str + 'may only contain 1 process tag.')
        return
    if (nL == 0):
        print(TF_XML_str + 'must contain 1 or more layer tags.')
        return
    if (nT > 1) or (nT == 0):
        print(TF_XML_str + 'must contain 1 timestep tag.')
        return
    
    #-----------------------------------------
    # There should be 1 top-level child node
    # that has the tag name "process"
    #-----------------------------------------
    # process_node = dom.childNodes[0]
    process_node = dom.firstChild

    #----------------------------------------------
    # Extract process name and store it in p_info
    #----------------------------------------------
    
    #---------------------------------------------
    # Prepare to save process info from XML file
    #---------------------------------------------
    p_info = TF_Process_Info(n_layers=nL, n_variables=nV/nL)
        
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
            #-----------------------------------------------------
            S_nodes = V_node.getElementsByTagName("label")
            label  = S_nodes[0].firstChild.data.strip()
            #-----------------------------------------------------
            vtypes = V_node.getElementsByTagName("type")
            vtype  = vtypes[0].firstChild.data.strip()
            #-----------------------------------------------------
            values = V_node.getElementsByTagName("value")
            value  = values[0].firstChild.data.strip()
            #-----------------------------------------------------
            units  = V_node.getElementsByTagName("units")
            unit   = units[0].firstChild.data.strip()
            #-----------------------------------------------------
            tlists = V_node.getElementsByTagName("typelist")
            tlist  = tlists[0].firstChild.data.strip().split(",")
            #-----------------------------------------------------
            p_info.layers[k].variables[j].name   = name
            p_info.layers[k].variables[j].label = label
            p_info.layers[k].variables[j].type   = vtype
            p_info.layers[k].variables[j].value  = value
            p_info.layers[k].variables[j].units  = unit
            p_info.layers[k].variables[j].type_choices = tlist
            #-----------------------------------------------------
##            p_info.layers[k].var_names.append( name )
##            p_info.layers[k].var_labels.append( label )
##            p_info.layers[k].var_types.append( vtype )
##            p_info.layers[k].var_values.append( value )
##            p_info.layers[k].var_units.append( unit )
##            p_info.layers[k].var_type_choices.append( tlist )
            #-----------------------------------------------------
            j += 1
        k += 1
    
    #----------------------------
    # Read timestep information
    #----------------------------
    T_nodes = process_node.getElementsByTagName("timestep")
    L_nodes = T_nodes[0].getElementsByTagName("label")
    V_nodes = T_nodes[0].getElementsByTagName("value")
    U_nodes = T_nodes[0].getElementsByTagName("units")
    #-------------------------------------------------------
    p_info.timestep.name  = L_nodes[0].firstChild.data
    p_info.timestep.value = V_nodes[0].firstChild.data
    p_info.timestep.units = U_nodes[0].firstChild.data

    #------------------------------
    # Read name of HTML help file
    #------------------------------
    H_nodes = process_node.getElementsByTagName("help_file")
    if (len(H_nodes) > 0):
        p_info.help_file = H_nodes[0].firstChild.data.strip()
        # print "help file =", p_info.help_file      

    return p_info

#   Read_XML_Info()
#----------------------------------------------------------------
        
#---------------------------------------
#  Support two different usage options
#---------------------------------------
##if (__name__ == '__main__'):
