
import wx

#----------------------------------------------------------------------
class MyTabbedDlg(wx.Dialog):
    def __init__(self, parent):
        title = "Resize the dialog and see how controls adapt!"
        wx.Dialog.__init__(self, parent, -1, title,
                           style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)

        notebook = wx.Notebook(self, -1, size=(450,300))
        panel1 = wx.Panel(notebook)
        panel2 = wx.Panel(notebook)
        notebook.AddPage(panel1, "Panel 1")
        notebook.AddPage(panel2, "Panel 2")

        dialog_sizer = wx.BoxSizer(wx.VERTICAL)
        dialog_sizer.Add(notebook, 1, wx.EXPAND|wx.ALL, 5)

        panel1_sizer = wx.BoxSizer(wx.VERTICAL)
        text = wx.TextCtrl(panel1, -1, "Hi!", size=(400,90), style=wx.TE_MULTILINE)
        button1 = wx.Button(panel1, -1, "I only resize horizontally...")
        panel1_sizer.Add(text, 1, wx.EXPAND|wx.ALL, 10)
        panel1_sizer.Add(button1, 0, wx.EXPAND|wx.ALL, 10)
        panel1.SetSizer(panel1_sizer)

        panel2_sizer = wx.BoxSizer(wx.HORIZONTAL)
        button2 = wx.Button(panel2, -1, "I resize vertically")
        button3 = wx.Button(panel2, -1, "I don't like resizing!")
        panel2_sizer.Add(button2, 0, wx.EXPAND|wx.ALL, 20)
        panel2_sizer.Add(button3, 0, wx.ALL, 100)
        panel2.SetSizer(panel2_sizer)

        if "__WXMAC__" in wx.PlatformInfo:
           self.SetSizer(dialog_sizer)
        else:
           self.SetSizerAndFit(dialog_sizer)
        self.Centre()

        self.Bind(wx.EVT_BUTTON, self.OnButton)


    def OnButton(self, evt):
        self.EndModal(0)


#----------------------------------------------------------------------
class MyApp(wx.App):
    def OnInit(self):
        dlg = MyTabbedDlg(None)
        dlg.ShowModal()
        dlg.Destroy()
        return True

myapp = MyApp(redirect=False)
myapp.MainLoop()
