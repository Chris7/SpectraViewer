import wx, MegaGrid, fileIterators, operator, os, re, wx.aui, ConfigParser, random
from Bio.Data import IUPACData

import numpy as np

# Define File Drop Target class
#from http://wiki.wxpython.org/DragAndDrop
class FileDropTarget(wx.FileDropTarget):
   """ This object implements Drop Target functionality for Files """
   def __init__(self, obj, pFrame):
      """ Initialize the Drop Target, passing in the Object Reference to
          indicate what should receive the dropped files """
      # Initialize the wxFileDropTarget Object
      wx.FileDropTarget.__init__(self)
      # Store the Object Reference for dropped files
      self.pFrame = pFrame

   def OnDropFiles(self, x, y, filenames):
      """ Implement File Drop """
      for file in filenames:
#         print file
         self.pFrame.addPage(file)

class ViewerGrid(MegaGrid.MegaGrid):
    def __init__(self, *args, **kwrds):
        MegaGrid.MegaGrid.__init__(self, *args, **kwrds)
        self.gridType = None
        self.dataSet = {}
        
    def AppendRow(self, inputData, **kwrds):
        row=self._table.GetNumberRows()
        if kwrds:
            if kwrds.has_key('unique'):
                tup = tuple([inputData[i] for i in kwrds['unique']])
                if self.dataSet.has_key(tup):
                    return
                self.dataSet[tup] = row
            if kwrds.has_key('group'):
                tup = tuple([inputData[i] for i in kwrds['group']])
                try:
                    index = self.dataSet[tup]
                    for i in xrange(0,len(inputData)):
                        val = self._table.GetValue(index,i)
                        if val and val != inputData[i] and inputData[i] not in val:
                            self._table.SetValue(index,i,val+','+inputData[i])
                    #del self.dataSet[tup]
                    #MegaGrid.MegaGrid.SetRowValue(self,index,val)
                    return
                except KeyError:
                    pass
                self.dataSet[tup] = row
        MegaGrid.MegaGrid.AppendRow(self,row,inputData)
        
    def setType(self, gridType):
        self.gridType = gridType.lower()
        
    def cellPopup(self, row, col, event):
        """(row, col, evt) -> display a popup menu when a cell is right clicked"""
        if not self.gridType:
            return
        viewSequence = wx.NewId()
        x = self.GetRowSize(row)/2
        y = self.GetColSize(col)/2

        #menu = wx.Menu()
        xo, yo = event.GetPosition()
        #menu.Append(viewSequence, "View Sequence")
        tidList = {}
        
        def onViewMultiSequence(event, self=self, row=row):
            self.parent.parent.loadScan(tidList[event.GetId()])
            return
        
        #def onViewSequence(event, self=self, row=row):
        if self.gridType == 'spectra':
            title = self._table.GetValue(row, 0)
            #path = self._table.GetHiddenValue(row, 0)
            self.parent.parent.loadScan(title)
            return
        elif self.gridType == 'pepspectra':
            title = self._table.GetValue(row, 0)
            titles = title.split(',')
#            if len(titles) > 1:
#                print type(self.parent)
            menua = wx.Menu()
            for i in titles:
                tid = wx.NewId()
                menua.Append(tid, i)
                tidList[tid] = i
                self.Bind(wx.EVT_MENU,onViewMultiSequence,id=tid)
#            else:
#                self.parent.parent.loadGFFScan(title)
            self.PopupMenu(menua)
            menua.Destroy()
            return
        else:
            return
            
#        self.Bind(wx.EVT_MENU, onViewSequence, id=viewSequence)
#        self.PopupMenu(menu)
#        menu.Destroy()
#        return

#masses from: http://www.weddslist.com/ms/tables.html#tm4
#protein weights with some mascot custom mods too
protein_weights =  {'G': 57.021464,
                    'A': 71.037114,
                    'S': 87.032028,
                    'P': 97.052764,
                    'V': 99.068414,
                    'T': 101.047678,
                    'C': 103.009184,
                    'I': 113.084064,
                    'L': 113.084064,
                    'N': 114.042927,
                    'D': 115.026943,
                    'Q': 128.058578,
                    'K': 128.094963,
                    'E': 129.042593,
                    'M': 131.040485,
                    'H': 137.058912,
                    'F': 147.068414,
                    'R': 156.101111,
                    'Y': 163.063329,
                    'W': 186.079313
                    }

mod_weights = {'h': 1.007825,
               'h2o': 18.010565,
               'nh3': 17.026549,
               'ch2': 14.015650,
               'methylation': 14.015650,
               'o': 15.994915,
               'oxidation': 15.994915,
               'acetylation': 42.010565,
               'carbamidation': 57.021464,
               'carboxylation': 58.005479,
               'phosphorylation': 79.966330,
               'amidation': 0.984016,
               'formylation': 27.994915,
               'cho': 29.002739665,
               'nh2': 16.01872407,
               'co': 27.99491463,
               'oh': 17.00274
               }
modParse = re.compile('([A-Z])(\d+)\((.+)\)')
class figureIons(object):
    def __init__(self,seq,charge,mods, tolerance):
        self.sequence=seq.upper()
        mass = 0
        self.a=[]
        self.b=[]
        self.c=[]
        self.x=[]
        self.y=[]
        self.z=[]
        modList = {}
        charge = int(charge)
        hcharge = mod_weights['h']*charge
        if mods:
            mods = mods.split('|')
            for mod in mods:
                modification, start, modType = modParse.search(mod).groups()
                modList[int(start)] = modType.lower()
        for i,v in enumerate(self.sequence):
            mass+=protein_weights[v]
            try:
                mass+=mod_weights[modList[i+1]]
#                print v,'b ion modified by',modList[i+1]
            except KeyError:
                pass
            self.a.append((mass-mod_weights['cho']+mod_weights['h']+hcharge)/charge)
            self.b.append((mass+mod_weights['h']-mod_weights['h']+hcharge)/charge)
            self.c.append((mass+mod_weights['nh2']+mod_weights['h']+hcharge)/charge)
        self.c.pop()
        self.c.append(9999999)
#        mass+=18.010565
        sLen = len(self.sequence)
        for i,v in enumerate(self.sequence):
            if i == 0:
                self.x.append(0)
            else:
                self.x.append((mass+mod_weights['co']+mod_weights['oh']-mod_weights['h']+hcharge)/charge)
            self.y.append((mass+mod_weights['h']+mod_weights['oh']+hcharge)/charge)
            self.z.append((mass-mod_weights['nh2']+mod_weights['oh']+hcharge)/charge)
            mass-=protein_weights[v]
#            print i,sLen,v
            try:
                mass-=mod_weights[modList[i+1]]
#                print v, 'y ion modified by',modList[i+1]
            except KeyError:
                pass
        self.tolerance = tolerance
        #self.y.reverse()
#        print self.a
#        print self.b
#        print self.c
#        print self.x
#        print self.y
#        print self.z

    def assignABC(self, x, y, atype):
        """
        given a list of masses, assigns ions from a sequence
        """
        start = 0
        self.bSeries = []
        if atype == 'a':
            iseq = self.a
        if atype == 'b':
            iseq = self.b
        if atype == 'c':
            iseq = self.c
        for seqindex,b in enumerate(iseq):
            candidates = []
            seen = False
            for index,m in enumerate(x[start:]):
                if b > m-self.tolerance and b < m+self.tolerance:
                    seen = True
                    candidates.append((x[start+index],y[start+index],self.sequence[seqindex], index))
                elif b > m+self.tolerance and seen:
                    start=index
                    break
            if candidates:
                candidates.sort(key=operator.itemgetter(1))
                self.bSeries.append(candidates[0])
            else:
                self.bSeries.append(False)
        return self.bSeries
        
    def assignXYZ(self, x,y, atype):
        """
        given a list of masses, assigns ions from a sequence
        """
        #a = figureIons('GGGFGGGSSFGGGSGFSGGGFGGGGFGGGR',2)
        start = len(x)
        self.ySeries = []
#        print self.sequence
        seq = list(self.sequence)
        seq.reverse()
        seq = ''.join(seq)
        if atype == 'x':
            iseq = self.x
        if atype == 'y':
            iseq = self.y
        if atype == 'z':
            iseq = self.z
        for seqindex,yi in enumerate(iseq):
#            print yi
            candidates = []
            seen = False
#            print seqindex
            for index,m in enumerate(x[:start]):
                if yi > m-self.tolerance and yi < m+self.tolerance:
                    seen = True
                    candidates.append((x[index],y[index],seq[seqindex], index))
                elif yi > m+self.tolerance and seen:
                    start=index
                    break
            if candidates:
#                print len(candidates)
                candidates.sort(key=operator.itemgetter(1))
                self.ySeries.append(candidates[0])
            else:
                self.ySeries.append(False)
        return self.ySeries

#some dangerous globals for us
searchPaths = set([])
fileNames = {}
titleParse = re.compile('(.+?)\.\d+\.\d+\.\d+|\w+')

def loadFolder(path):
    global searchPaths,fileNames
    for root,dir,files in os.walk(path):
        for fileName in files:
            pos = fileName.find('.mgf')
            pos2 = fileName.find('.mgfi')
            if pos != -1 and pos2 == -1:
                fileNames[fileName[:pos]] = os.path.join(root,fileName)
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

class MainFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "MGF Viewer", size=(1000,700))
        loadConfig()
#        self.tabPanel = wx.Panel(self, -1, size=self.GetSize())
        self._mgr = wx.aui.AuiManager()
        self.nb = wx.aui.AuiNotebook(self, size=self.GetSize())
        self._mgr.AddPane(self.nb, wx.aui.AuiPaneInfo().Name("notebook_content").CenterPane().PaneBorder(False))
        self.tabs = {}
        #plugins = {"Gene Id":MegaGrid.MegaFontRendererFactory("red", "ARIAL", 8),}
        self.mgfFiles = {}
        self.gffFiles = {}
        self.pepFiles = {}
        dt3 = FileDropTarget(self.nb, self)
       # Link the Drop Target Object to the Text Control
        self.nb.SetDropTarget(dt3)
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
        self.ions = ions
        self.Refresh()
        
    def onPaint(self, event):
        dc = wx.PaintDC(self)
        if dc is None:
            dc = wx.ClientDC(self)
        dc.BeginDrawing()
        w = dc.GetSize()[0]
        h = dc.GetSize()[1]
        
#            self.SetSize((w,200))
#            h = 200
        #write the peptide
        dc.SetBrush(wx.RED_BRUSH)
        tpos=0
#        print self.ions
        try:
            x = self.ions['x']
        except KeyError:
            x = []
        try:
            y = self.ions['y']
        except KeyError:
            y = []
        try:
            z = self.ions['z']
        except KeyError:
            z = []
        try:
            a = self.ions['a']
        except KeyError:
            a = []
        try:
            b = self.ions['b']
        except KeyError:
            b = []
        try:
            c = self.ions['c']
        except KeyError:
            c = []
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
#            print th, h, tw, w
            dc.SetFont(wx.Font(fsize, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
            try:
                tw,th = dc.GetTextExtent(pp)
            except IndexError:
                break
            fsize+=1
        totalw = len(self.peptide)*tw
        
        cy = h/2
        lx = w/2-totalw/2
        for i,v in enumerate(self.peptide):
            dc.SetFont(wx.Font(fsize, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
            tw,th = dc.GetTextExtent(v)
            dc.SetTextForeground((0,0,0))
            dc.SetPen(wx.Pen((0,0,0)))
            dc.DrawText(v, lx-tw/2,cy-th/2)
            dc.SetFont(wx.Font(isize, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
#            print i,y[i]
            ixw,ixh = dc.GetTextExtent('x')
            iyw,iyh = dc.GetTextExtent('y')
            izw,izh = dc.GetTextExtent('z')
            iaw,iah = dc.GetTextExtent('a')
            ibw,ibh = dc.GetTextExtent('b')
            icw,ich = dc.GetTextExtent('c')
            dc.DrawLine(lx+tw/2,cy-th/2-ixh/2-iyh-izh,lx+tw/2,cy+th/2+iah+ibh+ich/2)
            dc.SetTextForeground((255,0,255))
            ind = str(len(self.peptide)-i)
            iw,ih = dc.GetTextExtent(ind)
            dc.DrawText(ind, lx+tw/2,cy-th/2-ixh/2-iyh-izh-ih)
            dc.DrawText(str(i+1), lx+tw/2,cy+th/2+iah+ibh+ich/2)
            if x and x[i]:
                #we're x/y/z ions, we go in reverse
                xy = cy-th/2-izh-iyh
                dc.SetPen(wx.Pen((0,0,255)))
                dc.DrawLine(lx+tw,xy,lx+tw/2,xy)
                dc.SetTextForeground((0,0,255))
                dc.DrawText("x", lx+tw,xy-ixh/2)
            if y and y[i]:
#                print y[i]
                #we're x/y/z ions, we go in reverse
                yy = cy-th/2-ixh
                dc.SetPen(wx.Pen((0,0,255)))
                dc.DrawLine(lx+tw,yy,lx+tw/2,yy)
                dc.SetTextForeground((0,0,255))
                dc.DrawText("y", lx+tw,yy-iyh/2)
            if z and z[i]:
#                print y[i]
                #we're x/y/z ions, we go in reverse
                zy = cy-th/2
                dc.SetPen(wx.Pen((0,0,255)))
                dc.DrawLine(lx+tw,zy,lx+tw/2,zy)
                dc.SetTextForeground((0,0,255))
                dc.DrawText("z", lx+tw,zy-izh/2)
            if a and a[i]:
                dc.SetPen(wx.Pen((255,0,0)))
                ay = cy+th/2
                dc.DrawLine(lx,ay,lx+tw/2,ay)
                dc.SetTextForeground((255,0,0))
                dc.DrawText("a", lx-iw,ay-iah/2)
            if b and b[i]:
                by = cy+th/2+iah
                dc.SetPen(wx.Pen((255,0,0)))
                dc.DrawLine(lx,by,lx+tw/2,by)
                dc.SetTextForeground((255,0,0))
                dc.DrawText("b", lx-iw,by-ibh/2)
            if c and c[i]:
                cy2 = cy+th/2+iah+ibh
                dc.SetPen(wx.Pen((255,0,0)))
                dc.DrawLine(lx,cy2,lx+tw/2,cy2)
                dc.SetTextForeground((255,0,0))
                dc.DrawText("c", lx-iw,cy2-ich/2)
                #dc.DrawLine()
            lx+=tw+5
            #draw any needed ions
            
#        print 'done drawing',self.peptide,w,h
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
                        self.dataGrid.AppendRow([sid,i.getAttribute('Sequence'),i.getAttribute('Modifications'), i.getAttribute('Charge'), i.getAttribute('Name')], group=[4])
    #                if rnum>50:
    #                    break
            elif self.fileType == 'xml':
                pepParser = fileIterators.spectraXML(path)
                self.parent.pepFiles[path] = pepParser
                for i in pepParser.getScans():
                    self.dataGrid.AppendRow([i.getId(), i.getPeptide(), i.getModifications(), i.getCharge(),i.getAccession()], group=[4])
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
                self.dataGrid.AppendRow([i.getTitle(),i.getCharge(),i.getRT(), i.getPrecursor()])
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
        self.loadScan(self.title)
        
    def loadScan(self, title):
#        load a scan from a gff3 of mascot/X!Tandem output
        path = self.path
        self.title=title
#        print 'loading',title
        if self.fileType == 'gff':
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
    #        print title,scan
            if not scan:
                title = gob.getAttribute('SpectraId2')
                scan = self.parent.mgfFiles[path].getScan(title)
    #            print title,scan
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
            a = figureIons(self.pepSequence,gob.getAttribute('Charge'),mods, self.getTolerance())
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
            a = figureIons(self.pepSequence,1,mods, self.getTolerance())#for some reason x!Tandem treats everything as a single charge
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
        ionList = {}
        if self.draw.ionView['x']:
            ionList['x'] = a.assignXYZ(x,y,'x')
            self.draw.plotIons(ionList['x'], 'x')
        if self.draw.ionView['y']:
            ionList['y'] = a.assignXYZ(x,y,'y')
#            print 'y ions',ionList['y']
            self.draw.plotIons(ionList['y'], 'y')
        if self.draw.ionView['z']:
            ionList['z'] = a.assignXYZ(x,y,'z')
            self.draw.plotIons(ionList['z'], 'z')
        if self.draw.ionView['a']:
            ionList['a'] = a.assignABC(x,y,'a')
            self.draw.plotIons(ionList['a'], 'a')
        if self.draw.ionView['b']:
            ionList['b'] = a.assignABC(x,y,'b')
#            print 'b ions', ionList['b']
            self.draw.plotIons(ionList['b'], 'b')
        if self.draw.ionView['c']:
            ionList['c'] = a.assignABC(x,y,'c')
            self.draw.plotIons(ionList['c'], 'c')
        self.draw.peptidePanel.plotPeptide(self.pepSequence,ionList)
        
#class SettingsPanel(wx.Panel):
#    def __init__(self, parent):
#        wx.Panel.__init__(self, parent, -1)
#        self.parent = parent
#        #ion selection:
#        self.ions = wx.CheckListBox(self, -1, choices=['x', 'y', 'z', 'unmatched'])
#        self.ions.SetChecked([1])
#        self.ions.Bind(wx.EVT_CHECKLISTBOX, self.onChoice)
#        self.ions2 = wx.CheckListBox(self, -1, choices=['a', 'b', 'c'])
#        self.ions2.SetChecked([1])
#        self.ionLabel = wx.StaticText(self,-1,'Ions to Show:')
#        self.error = wx.TextCtrl(self, -1, "0.01")
#        self.errorLabel = wx.StaticText(self,-1,'Mass Error Tolerance (da)')
#        sizer = wx.FlexGridSizer()
#        sizer.SetFlexibleDirection(wx.HORIZONTAL)
#        sizer.Add(self.ionLabel)
#        sizer.Add(self.ions)
#        sizer.Add(self.ions2)
#        sizer.Add(self.errorLabel)
#        sizer.Add(self.error) 
#        self.SetSizer(sizer)
#        
#    def onChoice(self, event):
#        self.parent.reloadScan()
        
               
class PlotPanel(wx.Panel):
    #modified from this sourcE:
    #http://www.scipy.org/Matplotlib_figure_in_a_wx_panel
    """The PlotPanel has a Figure and a Canvas. OnSize events simply set a 
    flag, and the actual resizing of the figure is triggered by an Idle event."""
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
#        self.vbox.AddGrowableRow(1)
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
        for i in xrange(0,30,random.randint(2,7)):
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
        self.etolerance = wx.TextCtrl(self.toolbar, -1, "0.01",size=(30,tsize))
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
    
    def plotIons(self, ions, ionType):
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
#        print 'adding text'
        for x,y,aa,ind in zip(x,y,aas,nums):
            if not aa:
                continue
            self.hitMapXF[x] = (aa,y)
            try:
                self.hitMapX[int(x)].append((aa,y))
            except:
                self.hitMapX[int(x)] = [(aa,y)]
#            print ionType,x,y+5,ionType+str(ind+1)
            txt.append((x,y+2,ionType+str(ind+1),col))
#            txt.append((x,y+5,aa,col))
        self.text.append(txt)
        self.draw()
   
    def plotXY(self, xco, yco):
#        self.Canvas.ClearAll()
        self.x = xco
        self.y = yco
        self.points.append(([xco,yco],'spectra'))
        self.colors.append([0.0,0.0,0.0])
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
        for i,pt_list in enumerate(self.points):
            plot_pts = np.array( pt_list[0] )
            pType = pt_list[1]
#            print plot_pts[0,:]
            xm = np.max(plot_pts[0,:])
            if (xm > xmax):
                xmax = xm
            ym = np.max(plot_pts[1,:])
            if (ym > ymax):
                ymax = ym
            for xco,yco in zip(plot_pts[0,:],plot_pts[1,:]):
                try:
                    self.hitMapX[xco].add(yco)
                except KeyError:
                    self.hitMapX[xco] = set([yco])
            for barEntry in self.subplot.bar( plot_pts[0,:], plot_pts[1,:], color=self.colors[i],align='center'):
#                barEntry.set_linewidth(0)
                barEntry.set_edgecolor(self.colors[i])
                if pType == 'spectra':
                    barEntry.set_zorder(0)
                else:
                    barEntry.set_zorder(1)
        for i in self.text:
            for x,y,text,c in i:
                self.subplot.text(x,y,text,color=c)
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
            if x-1 < i < x+1:
                yco = self.hitMapX[i]
                for j in yco:
                    if y<=j:
                        hList.append((i-x,i,j))
        if self.annotate:
            self.subplot.axes.texts.remove(self.annotate)
            self.annotate = None
        if hList:
            sorted(hList,key=operator.itemgetter(0))
            txt+='m/z = %0.3f (Int: %0.3f)\n'%(hList[0][1],hList[0][2])
            self.annotate = self.subplot.axes.annotate(txt,
            (x,y), xytext=(-2*20, 5), textcoords='offset points',
            bbox=self.bbox, arrowprops=self.arrowprops)
        self.canvas.draw()


app = wx.App(False)
frame = MainFrame()
#frame.addPage('C:\Users\Chris\Desktop\A1.2012_06_07_12_20_00.t.xml')
#import wx.lib.inspection
#wx.lib.inspection.InspectionTool().Show()
app.MainLoop()