import wx

class PanelOne(wx.Panel):
    def __init__(self,parent):
            print """panel One Constructor"""
            wx.Panel.__init__(self, parent, -1)
            self.parent = parent
            self.label = wx.StaticText(self, -1, "panel one")
            self.switch_panel = wx.Button(self, -1, "Switch to panel two",
            (50,50))
            self.Bind(wx.EVT_BUTTON, self.__OnButton, self.switch_panel)
            self.__do_layout()

    def __do_layout(self):
            self.grid_sizer = wx.GridSizer(5, 5, 4, 4)
            self.grid_sizer.Add(self.label, 0, wx.FIXED_MINSIZE, 0)
            self.grid_sizer.Add(self.switch_panel, 0, wx.FIXED_MINSIZE, 0)

            self.SetAutoLayout(True)
            self.SetSizer(self.grid_sizer)

            self.grid_sizer.Fit(self)
            self.grid_sizer.SetSizeHints(self)

    def __OnButton(self,event):
            print "OnButton"
            self.parent.loadNewPanel(self)

class PanelTwo(wx.Panel):
    def __init__(self,parent):
            print """panel Two Constructor"""
            wx.Panel.__init__(self, parent, -1)
            self.parent = parent
            self.label = wx.StaticText(self, -1, "panel two")
            self.switch_panel = wx.Button(self, -1, "Switch to panel one",
            (50,50))
            self.Bind(wx.EVT_BUTTON, self.__OnButton, self.switch_panel)
            self.__do_layout()

    def __do_layout(self):
            self.grid_sizer = wx.GridSizer(5, 5, 4, 4)
            self.grid_sizer.Add(self.label, 0, wx.FIXED_MINSIZE, 0)
            self.grid_sizer.Add(self.switch_panel, 0, wx.FIXED_MINSIZE, 0)

            self.SetAutoLayout(True)
            self.SetSizer(self.grid_sizer)

            self.grid_sizer.Fit(self)
            self.grid_sizer.SetSizeHints(self)

    def __OnButton(self,event):
            print "OnButton"
            self.parent.loadNewPanel(self)
