import wx
from panels import *

active_frame = 1

class MyFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.panel = PanelOne(self)
        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle("frame_1")

    def __do_layout(self):
        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_1.Add(self.panel, 1, wx.EXPAND, 0)
        self.SetAutoLayout(True)
        self.SetSizer(self.sizer_1)
        self.sizer_1.Fit(self)
        self.sizer_1.SetSizeHints(self)
        self.Layout()

    def loadNewPanel(self,invoker):
        if isinstance(invoker,PanelOne):
            print "loading panel two"
            self.panel = PanelTwo(self)
        else:
            print "loading panel one"
            self.panel = PanelOne(self)
            self.sizer_1.Fit(self)
            self.sizer_1.SetSizeHints(self)

class MyApp(wx.App):
    def OnInit(self):
        wx.InitAllImageHandlers()
        frame = MyFrame(None, -1, "This is a wx.Frame", \
                        pos=(0,0), size=(640,480), \
                        style = wx.DEFAULT_FRAME_STYLE)
        self.SetTopWindow(frame)
        frame.Show()
        return 1

if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()
