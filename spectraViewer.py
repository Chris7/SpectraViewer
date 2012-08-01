import wx, MegaGrid, fileIterators, operator, os, re, wx.aui, ConfigParser, random, matplotlib as mpl, figureIons
import numpy as np

#some dangerous globals for us
searchPaths = set([])
fileNames = {}
titleParse = re.compile('(.+?)\.\d+\.\d+\.\d+|\w+')

def loadFolder(path):
    global searchPaths,fileNames
    for root,dir,files in os.walk(path):
        for fileName in files:
            if '.mgf' in fileName and '.mgfi' not in fileName:
                fileNames[fileName[:fileName.find('.mgf')]] = os.path.join(root,fileName)
    searchPaths.add(path)

def loadConfig():
    global searchPaths
    config = ConfigParser.RawConfigParser()
    config.read('spectra.cfg')
    try:
        paths = config.get('File Paths', 'pathList')
    except ConfigParser.NoSectionError:
        return
    searchPaths = set(paths.split(','))
    for i in searchPaths:
        loadFolder(i)
    
def saveConfig():
    global searchPaths
    config = ConfigParser.RawConfigParser()
    config.add_section('File Paths')
    paths = ','.join(searchPaths)
    config.set('File Paths', 'pathList', paths)
    with open('spectra.cfg', 'wb') as configfile:
        config.write(configfile)

# Define File Drop Target class
#from http://wiki.wxpython.org/DragAndDrop
class FileDropTarget(wx.FileDropTarget):
   """ This object implements Drop Target functionality for Files """
   def __init__(self, parent):
      """ Initialize the Drop Target, passing in the Object Reference to
          indicate what should receive the dropped files """
      # Initialize the wxFileDropTarget Object
      wx.FileDropTarget.__init__(self)
      # Store the Object Reference for dropped files
      self.parent = parent

   def OnDropFiles(self, x, y, filenames):
      """ Implement File Drop """
      for file in filenames:
         self.parent.addPage(file)

class ViewerGrid(MegaGrid.MegaGrid):
    def __init__(self, *args, **kwrds):
        MegaGrid.MegaGrid.__init__(self, *args, **kwrds)
        self.gridType = None
        self.dataSet = {}#used for storing groups of protein
        self.groupBy = None
        self.newGroup = False
        self.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.onLabelLeftClick)
        
    def onLabelLeftClick(self, evt):
        row, col = evt.GetRow(), evt.GetCol()
        if self.groupBy:
            if col == -1: self.rowUnGroup(row, evt)
            
    def setRowSize(self,row, size):
        self.SetRowMinimalHeight(row,size)
        self.SetRowSize(row,size)
                
    def rowUnGroup(self, row, event):
        grouper = self.dataSet[self._table.GetValue(row,self.groupBy)]
        pid = grouper[0][0]
        if isinstance(pid,int):
            return
        if pid[:-3] != str(row):
            return
        if pid[-2] == '+':
            txt = pid.replace('+','-')
            self._table.SetRowLabel(int(pid[:-3]),txt)
            rsize = 20
        else:
            txt = pid.replace('-','+')
            self._table.SetRowLabel(int(pid[:-3]),txt)
            rsize = 0
        grouper[0] = (txt,grouper[0][1])
        for i in grouper[1:]:
            self.setRowSize(i[0],rsize)
        self.Reset()
        
    def colPopup(self, col, evt):
        """(col, evt) -> display a popup menu when a column label is
        right clicked"""
        x = self.GetColSize(col)/2
        menu = wx.Menu()
        id1 = wx.NewId()
        sortID = wx.NewId()
        groupID = wx.NewId()

        xo, yo = evt.GetPosition()
        self.SelectCol(col)
        cols = self.GetSelectedCols()
        self.Refresh()
        menu.Append(sortID, "Sort Column (disabled after grouping)")
        menu.Append(groupID, "Group By Column")

        def sort(event, self=self, col=col):
            self.SortColumn(col, False)
            self.Reset()
            
        def group(event, self=self, col=col):
            self.groupBy=col
            self.newGroup = True
            self.regroup()

        self.Bind(wx.EVT_MENU, group, id=groupID)

        if len(cols) == 1:
            self.Bind(wx.EVT_MENU, sort, id=sortID)

        self.PopupMenu(menu)
        menu.Destroy()
        return
    
    def SortColumn(self, col, reindex):
        """
        col -> sort the data based on the column indexed by col
        """
        name = self._table.colIndexes[col]
        _data = []

        for row in self._table.data:
            rowname, entry = row
            _data.append((entry.get(name, None), row))

        _data.sort()
        if self.newGroup:
            self._table.data = []
            self.newGroup = False
            for index,row in enumerate(_data):
                    if reindex: #changed so the row index is changed with sorting
                        self._table.data.append((index,row[1][1]))
                    else:
                        self._table.data.append(row[1])
        elif not self.groupBy:
            self._table.data = []
            for row in _data:
                self._table.data.append(row)
    
    def regroup(self):
        #have we grouped already?
        if self.dataSet:
            for i in self.dataSet:
                if isinstance(self.dataSet[i][0][0],int):
                    continue
                rid = int(self.dataSet[i][0][0][:-3])
                self._table.SetRowLabel(rid,rid)
                for j in self.dataSet[i][1:]:
                    self.setRowSize(j[0],20)
        #first we sort so things get nested properly
        self.SortColumn(self.groupBy, True)
        self.SetRowMinimalAcceptableHeight(0)
        self.dataSet = {}
        gset = self.dataSet
        col = self._table.GetColLabelValue(self.groupBy)
        toDelete = []
        for i in self._table.data:
            groupby = i[1][col]
            try:
                gset[groupby]#we have it already? 
                gset[groupby].append(i)
                toDelete.append(i[0])
            except KeyError:
                #nope, we don't have it
                gset[groupby] = [i]
        for i in gset.keys():
            if len(gset[i])>1:
                row = gset[i][0]
                self._table.SetRowLabel(row[0],'%d(+)'%row[0])
                gset[i][0] = ('%d(+)'%row[0],row[1]) 
        if toDelete:
            for i in toDelete:
                self.setRowSize(i, 0)
        self.Reset()
                
    def appendRow(self, inputData, **kwrds):
        row=self._table.GetNumberRows()
        self.AppendRow(row,inputData)
        
    def setType(self, gridType):
        self.gridType = gridType.lower()
        
    def cellPopup(self, row, col, event):
        """(row, col, evt) -> display a spectra when cell is right clicked """
        if not self.gridType:
            return
        title = self._table.GetValue(row, 0)
        self.parent.parent.loadScan(title)
        return

class MainFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "MGF Viewer", size=(1000,700))
        loadConfig()
        self._mgr = wx.aui.AuiManager()
        self.nb = wx.aui.AuiNotebook(self, size=self.GetSize())
        self._mgr.AddPane(self.nb, wx.aui.AuiPaneInfo().Name("notebook_content").CenterPane().PaneBorder(False))
        self.tabs = {}
        self.mgfFiles = {}
        self.gffFiles = {}
        self.pepFiles = {}
        dt3 = FileDropTarget(self)
       # Link the Drop Target Object to the Text Control
        self.SetDropTarget(dt3)
        #menu bar
        menuBar = wx.MenuBar()
        menu = wx.Menu()
        settings  = wx.MenuItem(menu, -1, "Search Path", "Path to search for mgf files")
        menu.AppendItem(settings)
        menuBar.Append(menu, "Settings")
        self.Bind(wx.EVT_CLOSE, self.onClose)
        self.SetMenuBar(menuBar)
        self.Bind(wx.EVT_MENU, self.onSettings, settings)
        self.Show()
        
    def onClose(self, event):
        saveConfig()
        event.Skip()
        
    def onSettings(self, event):
        dialog = wx.DirDialog(None, "Choose a directory:",style=wx.DD_DEFAULT_STYLE)
        if dialog.ShowModal() == wx.ID_OK:
            loadFolder(dialog.GetPath())
        
    def addPage(self, path):
        self.tabs[path] = ViewerPanel(self,path,size=(1000,700))#SplitterWindow(self.nb, -1, size=self.GetSize())
        tabName = os.path.split(os.path.splitext(path)[0])[1]
        self.nb.AddPage(self.tabs[path], tabName)

class PeptidePanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.parent = parent
        self.peptide = ""
        self.ions = {}
        self.SetMaxSize((10000,100))
        self.SetMinSize((0,100))
        self.Bind(wx.EVT_PAINT, self.onPaint)
        self.Bind(wx.EVT_SIZE, self.onSize)
        
    def onSize(self, event):
        event.Skip()
        self.Refresh()
        
    def plotPeptide(self, pep, ions):
        self.SetSize((self.parent.GetSize()[0],100))
        self.peptide = pep
        self.sLen = len(pep)
        self.revTypes = set([])
        self.ions = ions
        self.Refresh()
        
    def onPaint(self, event):
        dc = wx.PaintDC(self)
        if dc is None:
            dc = wx.ClientDC(self)
        dc.BeginDrawing()
        w = dc.GetSize()[0]
        h = dc.GetSize()[1]
        dc.SetBrush(wx.RED_BRUSH)
        tpos=0
        lx=0
        ly=0
        tw=0
        th=0
        fsize = 1
        dc.SetFont(wx.Font(fsize, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
        try:
            pp = self.peptide[0]
        except IndexError:
            return
        for i in self.peptide:
            pw,ph = dc.GetTextExtent(pp)
            if pw > tw:
                tw = pw
                pp = i
            if ph > th:
                th = pw
                pp = i
        isize=1
        dc.SetFont(wx.Font(isize, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
        tw,th = dc.GetTextExtent("10")
        while th < h*1/10:
            isize+=1
            dc.SetFont(wx.Font(isize, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
            tw,th = dc.GetTextExtent("10")
        ih=th
        while (th*1.25+ih*2 <= h*1/2) and (tw*len(self.peptide) < w*1/2):
            dc.SetFont(wx.Font(fsize, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
            try:
                tw,th = dc.GetTextExtent(pp)
            except IndexError:
                break
            fsize+=1
        totalw = len(self.peptide)*tw
        dc.SetFont(wx.Font(isize, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
        ixw,ixh = dc.GetTextExtent('x')
        iyw,iyh = dc.GetTextExtent('y')
        izw,izh = dc.GetTextExtent('z')
        iaw,iah = dc.GetTextExtent('a')
        ibw,ibh = dc.GetTextExtent('b')
        icw,ich = dc.GetTextExtent('c')
        cy = h/2
        lx = w/2-totalw/2
        toDraw = {}
        sLen = self.sLen
        for entry in self.ions:
            m = entry[4].split('|')
            fType = m[0][0]
            index = int(m[0][1:])
            try:
                toDraw[index].append(fType)
            except KeyError:
                toDraw[index] = [fType]
        for i,v in enumerate(self.peptide):
            dc.SetFont(wx.Font(fsize, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
            tw,th = dc.GetTextExtent(v)
            dc.SetTextForeground((0,0,0))
            dc.SetPen(wx.Pen((0,0,0)))
            dc.DrawText(v, lx-tw/2,cy-th/2)
            if i==sLen-1:
                break
            dc.SetFont(wx.Font(isize, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
            dc.DrawLine(lx+tw/2,cy-th/2-ixh/2-iyh-izh,lx+tw/2,cy+th/2+iah+ibh+ich/2)
            dc.SetTextForeground((255,0,255))
            ind = str(sLen-i-1)
            iw,ih = dc.GetTextExtent(ind)
            dc.DrawText(ind, lx+tw/2,cy-th/2-ixh/2-iyh-izh-ih)
            dc.DrawText(str(i+1), lx+tw/2,cy+th/2+iah+ibh+ich/2)
            try:
                todraw = toDraw[i+1]
                for fragType in todraw:
                    if fragType == 'A':
                        dc.SetPen(wx.Pen((255,0,0)))
                        ay = cy+th/2
                        dc.DrawLine(lx,ay,lx+tw/2,ay)
                        dc.SetTextForeground((255,0,0))
                        dc.DrawText("a", lx-iw,ay-iah/2)
                    elif fragType == 'B':
                        by = cy+th/2+iah
                        dc.SetPen(wx.Pen((255,0,0)))
                        dc.DrawLine(lx,by,lx+tw/2,by)
                        dc.SetTextForeground((255,0,0))
                        dc.DrawText("b", lx-iw,by-ibh/2)
                    elif fragType == 'C':
                        cy2 = cy+th/2+iah+ibh
                        dc.SetPen(wx.Pen((255,0,0)))
                        dc.DrawLine(lx,cy2,lx+tw/2,cy2)
                        dc.SetTextForeground((255,0,0))
                        dc.DrawText("c", lx-iw,cy2-ich/2)
            except KeyError:
                pass
            try:
                todraw = toDraw[int(ind)]
                for fragType in todraw:
                    if fragType == 'X':
                        #we're x/y/z ions, we go in reverse
                        xy = cy-th/2-izh-iyh
                        dc.SetPen(wx.Pen((0,0,255)))
                        dc.DrawLine(lx+tw,xy,lx+tw/2,xy)
                        dc.SetTextForeground((0,0,255))
                        dc.DrawText("x", lx+tw,xy-ixh/2)
                    elif fragType == 'Y':
                        #we're x/y/z ions, we go in reverse
                        yy = cy-th/2-ixh
                        dc.SetPen(wx.Pen((0,0,255)))
                        dc.DrawLine(lx+tw,yy,lx+tw/2,yy)
                        dc.SetTextForeground((0,0,255))
                        dc.DrawText("y", lx+tw,yy-iyh/2)
                    elif fragType == 'Z':
                        #we're x/y/z ions, we go in reverse
                        zy = cy-th/2
                        dc.SetPen(wx.Pen((0,0,255)))
                        dc.DrawLine(lx+tw,zy,lx+tw/2,zy)
                        dc.SetTextForeground((0,0,255))
                        dc.DrawText("z", lx+tw,zy-izh/2)
                    
            except KeyError:
                pass
            lx+=tw+5
        dc.EndDrawing()
        
class ViewerPanel(wx.SplitterWindow):
    def __init__(self, parent, path, **kwrds):
        wx.SplitterWindow.__init__(self, parent, -1,**kwrds)
        self.parent = parent
        self.gridPanel = wx.Panel(self, -1)
        self.gridPanel.parent = self
        self.path = path
        pepSpecType = ('gff', 'xml')
        specType = ('mgf',)
        if [i for i in pepSpecType if i in path]:
            self.fileType = [i for i in pepSpecType if i in path][0]#this changes based on the file input somewhat
            colnames = ["Scan Title", "Peptide", "Modifications", "Charge", "Accession"]
            data = []
            plugins = {}
            self.dataGrid = ViewerGrid(self.gridPanel, data, colnames, plugins)
            self.dataGrid.setType('pepspectra')
            self.dataType = 'pepspectra'
            if self.fileType == 'gff':
                gffIterator = fileIterators.GFFIterator(path, random=['SpectraId1', 'SpectraId2'])
                self.parent.gffFiles[path] = gffIterator
    #            rnum=0
                for i in gffIterator:
                    if not i:
                        continue
    #                rnum+=1
                    sid = i.getAttribute('SpectraId1')
                    if sid:
                        self.dataGrid.appendRow([sid,i.getAttribute('Sequence'),i.getAttribute('Modifications'), i.getAttribute('Charge'), i.getAttribute('Name')])
    #                if rnum>50:
    #                    break
            elif self.fileType == 'xml':
                pepParser = fileIterators.spectraXML(path)
                self.parent.pepFiles[path] = pepParser
                for i in pepParser.getScans():
                    self.dataGrid.appendRow([i.getId(), i.getPeptide(), i.getModifications(), i.getCharge(),i.getAccession()])
        elif [i for i in specType if i in path]:
            self.fileType = 'spectra'#these are all generic more or less, so spectra works
            self.dataType = 'spectra'
            colnames = ["Scan Title", "Charge", "RT", "Precursor Mass"]
            data = []
            plugins = {}
            self.dataGrid = ViewerGrid(self.gridPanel, data, colnames, plugins)
            self.dataGrid.setType('spectra')
            mgf = fileIterators.mgfIterator(path, random=True)
            self.parent.mgfFiles[path] = mgf
            for i in mgf:
                if not i:
                    continue
                self.dataGrid.appendRow([i.getTitle(),i.getCharge(),i.getRT(), i.getPrecursor()])
        else:
            data = []
            plugins = {}
            colnames = ["none"]
            self.fileType = 'none'
            self.dataType = 'none'
            self.dataGrid = ViewerGrid(self.gridPanel, data, colnames, plugins)
            self.dataGrid.setType('none')
        megasizer = wx.BoxSizer(wx.VERTICAL)
        megasizer.Add(self.dataGrid, 1, wx.EXPAND)
        self.gridPanel.SetSizer(megasizer)
        self.draw = DrawFrame(self)
        self.SplitHorizontally(self.draw,self.gridPanel, sashPosition=self.GetSize()[1]*2/3)
        
    def getTolerance(self):
        return float(self.draw.etolerance.GetValue())
        
    def reloadScan(self):
        try:
            self.loadScan(self.title)
        except AttributeError:
            #no scans loaded yet
            return
        
    def loadScan(self, title):
        path = self.path
        self.title=title
        if self.fileType == 'gff':
#        load a scan from a gff3 of mascot/X!Tandem output
            try:
                it = self.parent.gffFiles[path]
            except KeyError:
                loadGFF(path)
                it = self.parent.gffFiles[path]
            gob = it.getGFF('SpectraId1',title)
            mods = gob.getAttribute('Modifications')
            fileName = titleParse.match(title).group(1)
            if not fileNames.has_key(fileName):
                return
            path = fileNames[fileName]
            try:
                scan = self.parent.mgfFiles[path].getScan(title)
            except KeyError:
                mgf = fileIterators.mgfIterator(path)
                self.parent.mgfFiles[path] = mgf
                try:
                    scan = mgf.getScan(title)
                except StopIteration:
                    scan = None
            if not scan:
                title = gob.getAttribute('SpectraId2')
                scan = self.parent.mgfFiles[path].getScan(title)
                if not scan:
                    return
            self.precursor = scan.getPrecursor()
            mz = scan.getMZ()
            x = []
            y = []
            for i in mz:
                x.append(float(i[0]))
                y.append(float(i[1]))
#            y = np.array(y)/np.max(y)*100
            self.pepSequence = gob.getAttribute('Sequence')
            a = figureIons.figureIons(self.pepSequence,gob.getAttribute('Charge'),mods, self.getTolerance())
            self.draw.cleanup()
            self.draw.setTitle(title)
            self.plotIons(x,y,a)
            if self.draw.ionView['all']:
                self.draw.plotXY(x,y)
        elif self.fileType == 'xml':
            scan = self.parent.pepFiles[path].getScan(title)
            if not scan:
                return
            mods = scan.getModifications()
            self.precursor = scan.getPrecursor()
            mz = scan.getMZ()
            x = []
            y = []
            for i in mz:
                x.append(float(i[0]))
                y.append(float(i[1]))
#            y = np.array(y)/np.max(y)*100
            self.pepSequence = scan.getPeptide()
            a = figureIons.figureIons(self.pepSequence,scan.getCharge(),mods, self.getTolerance())
            self.draw.cleanup()
            self.draw.setTitle(title)
            self.plotIons(x,y,a)
            if self.draw.ionView['all']:
                self.draw.plotXY(x,y)
        elif self.fileType == 'spectra':
            try:
                scan = self.parent.mgfFiles[path].getScan(title)
            except:
                mgf = fileIterators.mgfIterator(path)
                self.parent.mgfFiles[path] = mgf
                scan = self.parent.mgfFiles[path].getScan(title)
            if not scan:
                return
            mz = scan.getMZ()
            x = []
            y = []
            for i in mz:
                x.append(float(i[0]))
                y.append(float(i[1]))
            self.draw.cleanup()
            self.draw.plotXY(x,y)
        
    def plotIons(self, x,y,a):
        ionList = a.assignPeaks(x,y)
        self.draw.plotIons(ionList)
        self.draw.peptidePanel.plotPeptide(self.pepSequence,ionList)        
               
class PlotPanel(wx.Panel):
    #modified from this sourcE:
    #http://www.scipy.org/Matplotlib_figure_in_a_wx_panel
    def __init__( self, parent, color=None, dpi=None, **kwargs ):
        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
        from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
        from matplotlib.figure import Figure
        self.parent = parent
        # initialize Panel
        if 'id' not in kwargs.keys():
            kwargs['id'] = wx.ID_ANY
        if 'style' not in kwargs.keys():
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__( self, parent, **kwargs )

        # initialize matplotlib stuff
        self.figure = Figure( None, dpi )
        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )
        self.canvas.mpl_connect('button_press_event', self.onMouseDown)
        self.canvas.mpl_connect('button_release_event', self.onMouseRelease)
        self.canvas.SetMinSize((20,20))
        self.toolbar = NavigationToolbar(self.canvas)
        self.SetColor( color )
        self.vbox = wx.FlexGridSizer(1,1)
        self.vbox.AddGrowableCol(0)
        self.vbox.AddGrowableRow(2)
        self.peptidePanel = PeptidePanel(self)
        self.vbox.Add(self.toolbar,0,wx.EXPAND)
        self.vbox.Add(self.peptidePanel, 0, wx.EXPAND)
        self.vbox.Add(self.canvas,10,wx.EXPAND)
        self.SetSizer(self.vbox)
        self.SetAutoLayout(True)
        self.mouse = 0
        self.vbox.Fit(self)

    def SetColor( self, rgbtuple=None ):
        """Set figure and canvas colours to be the same."""
        if rgbtuple is None:
            rgbtuple = wx.SystemSettings.GetColour( wx.SYS_COLOUR_BTNFACE ).Get()
        clr = [c/255. for c in rgbtuple]
        self.figure.set_facecolor( clr )
        self.figure.set_edgecolor( clr )
        self.canvas.SetBackgroundColour( wx.Colour( *rgbtuple ) )
        
    def onMouseDown(self, event):
        self.mouse = 1
        
    def onMouseRelease(self, event):
        self.mouse = 0

class DrawFrame(PlotPanel):
    def __init__(self, parent, *args,**kwargs):
        self.parent = parent
        PlotPanel.__init__(self, parent, **kwargs)        
        self.cleanup()
        self.addFlag('y')
        self.addFlag('b')
        self.ionView = {'x': False, 'y': True, 'z': False, 'a': False, 'b': True, 'c': False, 'all': False}
        self.bbox = dict(boxstyle="round", fc="0.8")
        self.arrowprops = dict(
            arrowstyle = "->",
            connectionstyle = "angle,angleA=0,angleB=90,rad=10")
        self.annotate = None
        if self.parent.fileType == 'spectra':
            self.peptidePanel.Hide()
        self.canvas.mpl_connect('motion_notify_event', self.onMouseMove)
        #toolbar stuff
        #draw bitmaps for labels
        tsize = self.toolbar.GetSize()[1]
        idmap = {}
        for ion,desc in zip(('x','y','z','a','b','c'),('X ions', 'Y Ions', 'Z Ions', 'A Ions', 'B Ions', 'C Ions')):
            image = wx.EmptyBitmap(30,tsize*3/4)
            lid = wx.NewId()
            idmap[ion] = lid
            dc = wx.MemoryDC()
            dc.SetFont(wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
            dc.SelectObject(image)
            dc.SetBrush(wx.Brush(self.toolbar.BackgroundColour))
            dc.DrawRectangle(-1,-1,35,tsize)
            dc.SetTextForeground((255,0,0))
            dc.SetPen(wx.Pen((255,0,0)))
            dc.DrawText(ion,7,3)
            dc.SelectObject(wx.NullBitmap)
            self.toolbar.AddCheckLabelTool(lid, ion, bitmap=image,shortHelp = 'Show %s'%desc)
            if ion is 'y' or ion is 'b':
                self.toolbar.ToggleTool(lid, True)
        #a little icon to see all spectra
        image = wx.EmptyBitmap(30,self.toolbar.GetSize()[1]*3/4)
        dc = wx.MemoryDC()
        dc.SetFont(wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
        dc.SelectObject(image)
        dc.SetBrush(wx.Brush(self.toolbar.BackgroundColour))
        dc.DrawRectangle(-1,-1,35,tsize)
        dc.SetTextForeground((105,105,105))
        dc.SetPen(wx.Pen((105,105,105)))
        for i in xrange(0,30,3):
            rn = float('%0.3f'%random.uniform(0.35,1.0))
            dc.DrawRectangle(i,tsize*rn,2,tsize-tsize*rn)
        dc.SelectObject(wx.NullBitmap)
        lid = wx.NewId()
        self.toolbar.AddCheckLabelTool(lid, "Spectra", bitmap=image, shortHelp="Plot Entire spectra")
        self.Bind(wx.EVT_TOOL, self.onViewSpectra, id=lid)
        self.Bind(wx.EVT_TOOL, self.onXIons, id=idmap['x'])
        self.Bind(wx.EVT_TOOL, self.onYIons, id=idmap['y'])
        self.Bind(wx.EVT_TOOL, self.onZIons, id=idmap['z'])
        self.Bind(wx.EVT_TOOL, self.onAIons, id=idmap['a'])
        self.Bind(wx.EVT_TOOL, self.onBIons, id=idmap['b'])
        self.Bind(wx.EVT_TOOL, self.onCIons, id=idmap['c'])
        #error tolerance
        self.etolerance = wx.TextCtrl(self.toolbar, -1, "0.01",size=(50,tsize))
        self.etolerance.SetToolTip(wx.ToolTip("Mass Error Tolerance (da)"))
        self.toolbar.AddControl(self.etolerance)
        self.toolbar.Realize()
        
    def onViewSpectra(self, event):
        self.ionView['all'] = event.Checked()
        self.parent.reloadScan()
        
    def onXIons(self, event):
        self.ionView['x'] = event.Checked()
        self.parent.reloadScan()
        
    def onYIons(self, event):
        self.ionView['y'] = event.Checked()
        self.parent.reloadScan()
        
    def onZIons(self, event):
        self.ionView['z'] = event.Checked()
        self.parent.reloadScan()
        
    def onAIons(self, event):
        self.ionView['a'] = event.Checked()
        self.parent.reloadScan()
        
    def onBIons(self, event):
        self.ionView['b'] = event.Checked()
        self.parent.reloadScan()
        
    def onCIons(self, event):
        self.ionView['c'] = event.Checked()
        self.parent.reloadScan()
        
    def addFlag(self, flag):
        if not self.testFlag(flag):
            self.flag+=self.flags[flag]
        
    def removeFlag(self, flag):
        if self.testFlag(flag):
            self.flag-=self.flags[flag]
        
    def toggleFlag(self, flag):
        if self.testFlag(flag):
            self.removeFlag(flag)
        else:
            self.addFlag(flag)
        
    def testFlag(self, flag):
        return self.flag & self.flags[flag]
    
    def plotIons(self, peaks):
        itypes = {'x':[],
                  'y':[],
                  'z':[],
                  'a':[],
                  'b':[],
                  'c':[]
                  }
        for i in peaks:
            if not i:
                continue
            mz,inten,aa,itype,desc = i
            if self.ionView[itype[0].lower()]:
                self.points.append((mz,inten,aa,itype,desc))
        self.draw()
            
    
    def _plotIons(self, ions, ionType):
        x = []
        y = []
        aas = []
        for i in ions:
            if i == False:
                x.append(0)
                y.append(0)
                aas.append(False)
            else:
                x.append(i[0])
                y.append(i[1])
                aas.append(i[2])
        col = [0.0,0.0,0.0]
        if ionType.lower() in set(['x','y','z']):
            nums = [i for i in reversed(xrange(len(aas)))]
            col = [0.0,0.0,1.0]
        else:
            nums = [i for i in xrange(len(aas))]
            col = [1.0,0.0,0.0]
        self.points.append(([x,y],ionType))
        self.colors.append([1.0,0.0,0.0])
        txt = []
        for x,y,aa,ind in zip(x,y,aas,nums):
            if not aa:
                continue
            self.hitMapXF[x] = (aa,y)
            try:
                self.hitMapX[int(x)].append((aa,y))
            except:
                self.hitMapX[int(x)] = [(aa,y)]
            txt.append((x,y+2,ionType+str(ind+1),col))
        self.text.append(txt)
        self.draw()
   
    def plotXY(self, xco, yco):
#        self.Canvas.ClearAll()
        self.x = xco
        self.y = yco
        for x,y in zip(xco,yco):
            self.points.append((x,y,'','spectra',''))
        self.draw()

    def cleanup(self):
        self.figure.clear()
        self.subplot = self.figure.add_subplot( 111 )
        self.ions = {}
        self.points = []
        self.colors = []
        self.text = []
        self.flags = {'y':1,
                      'b':2,
                      'a':4,
                      'c':8,
                      'x':16,
                      'z':32,
                      'all':64}
        self.flag = 0
        self.hitMapX = {}
        self.hitMapXF = {}
        
    def setTitle(self, title):
        self.figure.suptitle(title)
        
    def draw(self):
        xmax = 10
        ymax = 10
        self.hitMapX = {}
        blues = set(('x','y','z'))
        for i,pt_list in enumerate(self.points):
            x,y,aa,fragtype,desc = pt_list
            fragtype = fragtype[0].lower()
            if fragtype in blues:
                col=[0.0,0.0,1.0]
                tcolor='b'
            elif fragtype == 'spectra':
                t='k'
                col=[0.0,0.0,0.0] 
            else:
                col=[1.0,0.0,0.0]
                tcolor='r'
            if (x > xmax):
                xmax = x
            if (y > ymax):
                ymax = y
            try:
                self.hitMapX[x].add((y,desc.replace('|',' ')))
            except KeyError:
                self.hitMapX[x] = set([(y,desc.replace('|', ' '))])
            bar = self.subplot.bar(x,y,color=col,align='center')
            if fragtype == 'spectra':
                bar[0].set_zorder(0)
            else:
                if not desc:
                    continue
                desc = desc.split('|')
                hlen=len(self.hitMapX[x])*5
                if len(desc) == 3:
                    txt = '$\mathtt{%s^{%s^{%s+}}}$'%(desc[0].lower(),desc[1],desc[2][:-1])
                else:
                    txt = '$\mathtt{%s^{%s+}}$'%(desc[0].lower(),desc[1][:-1])
                self.subplot.text(x,y+hlen,txt,color=tcolor)
        self.subplot.axes.set_xbound(lower=0, upper=xmax+10)
        self.subplot.axes.set_ybound(lower=0, upper=ymax+20)
        self.canvas.draw()
        
    def onMouseMove(self,event):
        if self.mouse:
            return
        inv = self.subplot.axes.transData.inverted()
        x,y = inv.transform((event.x,event.y))
        txt = ""
        hList = []
        for i in self.hitMapX:
            if x-3 < i < x+3:
                yco = self.hitMapX[i]
                for j in yco:
                    if y<=j[0]:
                        hList.append((abs(i-x),i,j))
        if self.annotate:
            try:
                self.subplot.axes.texts.remove(self.annotate)
            except ValueError:
                pass
            self.annotate = None
        if hList:
            sorted(hList,key=operator.itemgetter(0))
            txt+='%s m/z = %0.3f (Int: %0.3f)\n'%(hList[0][2][1],hList[0][1],hList[0][2][0])
            self.annotate = self.subplot.axes.annotate(txt,
            (x,y), xytext=(-2*20, 5), textcoords='offset points',
            bbox=self.bbox, arrowprops=self.arrowprops)
        self.canvas.draw()


app = wx.App(False)
frame = MainFrame()
#frame.addPage('A1.2012_06_07_12_20_00.t.xml')
#import wx.lib.inspection
#wx.lib.inspection.InspectionTool().Show()
app.MainLoop()