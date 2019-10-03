
#  April 23, 2009
#  S.D. Peckham

import wx

#----------------------------------------------------------------

# Get_Directory
# Labeled_Text_Box
# Labeled_Droplist
# Add_Back_Next_Buttons

# Precip_Methods
# Snowmelt_Methods
# ET_Methods
# Infil_Methods
# GW_Methods
# Channel_Methods
# Diversion_Methods
# Stop_Methods

#----------------------------------------------------------------    
def Get_Directory():

    dialog = wx.DirDialog(None, "Choose a directory:", \
                          style=wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON)
    if (dialog.ShowModal() == wx.ID_OK):
        directory = dialog.GetPath()
    else:
        directory = None
    dialog.Destroy()
    
    return directory
    
#   Get_Directory()
#----------------------------------------------------------------
def Labeled_Text_Box(frame, panel, sizer, label, text='', \
                     MULTI=False, \
                     button=None, function=None, \
                     border=5, box_width=300, box_height=-1):

    L1 = wx.StaticText(panel, -1, label)
    sizer.Add(L1, 0, wx.ALL, border)

    if (MULTI):
        T1 = wx.TextCtrl(panel, -1, text, \
                         size=(box_width, box_height), \
                         style=wx.TE_MULTILINE)
    else:
        T1 = wx.TextCtrl(panel, -1, text, \
                         size=(box_width, box_height))
    sizer.Add(T1, 0, wx.ALL, border)

    #-------------------------
    # Option to add a button
    #-------------------------
    if (button is None):
        null_item = wx.StaticText(panel, -1, "")
        sizer.Add(null_item, 0, wx.ALL, border)
    else:
        B1 = wx.Button(panel, -1, button)
        ### B1.Disable()  ###  (this works)
        sizer.Add(B1, 0, wx.ALL, border)
        frame.Bind(wx.EVT_BUTTON, function, B1)  ####
        
    return T1   # (text box object)

#   Labeled_Text_Box
#----------------------------------------------------------------
def Labeled_Droplist(frame, panel, sizer, \
                     label, choice_list, function, \
                     button1=None, b1_function=None, \
                     button2=None, b2_function=None, \
                     box_width=280, border=5, \
                     choice_number=0):

    L1 = wx.StaticText(panel, -1, label)
    sizer.Add(L1, 0, wx.ALL, border)
    
    D1 = wx.Choice(panel, -1, choices=choice_list, \
                   size=(box_width, -1))
    D1.Select(choice_number)
    sizer.Add(D1, 0, wx.ALL, border)
    frame.Bind(wx.EVT_CHOICE, function, D1)

    #---------------------------
    # Option to add 1st button
    #---------------------------
    if (button1 is None):
        null_item1 = wx.StaticText(panel, -1, "")
        sizer.Add(null_item1, 0, wx.ALL, border)
    else:
        B1 = wx.Button(panel, -1, button1)
        sizer.Add(B1, 0, wx.ALL, border)
        frame.Bind(wx.EVT_BUTTON, b1_function, B1)
        
    #---------------------------
    # Option to add 2nd button
    #---------------------------
    if (button2 is None):
        null_item2 = wx.StaticText(panel, -1, "")
        sizer.Add(null_item2, 0, wx.ALL, border)
    else:
        B2 = wx.Button(panel, -1, button2)
        sizer.Add(B2, 0, wx.ALL, border)
        frame.Bind(wx.EVT_BUTTON, b2_function, B2)
        
    return D1   # (droplist object)

#   Labeled_Droplist
#----------------------------------------------------------------
def Add_Back_Next_Buttons(frame, panel, sizer, \
                          back_fcn=None, next_fcn=None, \
                          hgap=20, vgap=10, PAD_TOP=True):
 
    back_button = wx.Button(panel, -1, "< Back")
    next_button = wx.Button(panel, -1, "Next >")
    btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
    btn_sizer.Add(back_button)
    btn_sizer.Add((hgap, hgap), 1)
    btn_sizer.Add(next_button)
    
    frame.Bind(wx.EVT_BUTTON, back_fcn, back_button)
    frame.Bind(wx.EVT_BUTTON, next_fcn, next_button)

    if (PAD_TOP):
        pad_row = wx.StaticText(panel, -1, "")    
        sizer.Add(pad_row,   0, wx.ALL, vgap)
    sizer.Add(btn_sizer, 0, wx.ALL, vgap)

#   Add_Back_Next_Buttons
#----------------------------------------------------------------
def Precip_Methods():

    mlist = ['None', \
             'Uniform in space, variable durations',\
             'Variable in space and time']
    return mlist

#   Precip_Methods()
#----------------------------------------------------------------
def Snowmelt_Methods():

    mlist = ['None', \
             'Degree-Day',\
             'Energy Balance']
    return mlist

#   Snowmelt_Methods()
#----------------------------------------------------------------
def ET_Methods():

    mlist = ['None', \
             'Priestley-Taylor',\
             'Energy Balance']
    return mlist

#   ET_Methods()
#----------------------------------------------------------------
def Infil_Methods():

    mlist = ['None', \
             '100% infiltration until (h >= z)',\
             'Simple Green-Ampt, one storm', \
             'Smith-Parlange 3-param., one storm', \
             'Richards 1D Equation 3 layers']
             # 'Beven Exponential K, one storm']  # (no GUI yet)
    return mlist

#   Infil_Methods()
#----------------------------------------------------------------
def GW_Methods():

    mlist = ['None', \
             "Darcy's Law, Surface-parallel layers" ]
    return mlist

#   GW_Methods()
#----------------------------------------------------------------
def Channel_Methods():

    mlist = ['None', \
             'Kinematic Wave, Manning friction', \
             'Kinematic Wave, Law of Wall friction', \
             'Diffusive Wave, Manning friction', \
             'Diffusive Wave, Law of Wall friction', \
             'Dynamic Wave, Manning friction', \
             'Dynamic Wave, Law of Wall friction' ]
    return mlist

#   Channel_Methods()
#----------------------------------------------------------------
def Diversion_Methods():

    mlist = ['None', \
             'Sources, Sinks or Canals']
    return mlist

#   Diversion_Methods()
#----------------------------------------------------------------
def Stop_Methods():

    mlist = ["Run until Q drops to P% of Q_peak", \
             "Run for a specified time (model)", \
             "Run for a specified number of steps"]
    return mlist

#   Stop_Methods()
#----------------------------------------------------------------

