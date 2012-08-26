import fileIterators, operator, os, re, ConfigParser, random, matplotlib as mpl, figureIons, collections, sys
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import numpy as np, time
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import SpecView

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


class PeptidePanel(QWidget):
    def __init__(self, parent):
        QWidget.__init__(self)
        self.parent = parent
        self.peptide = "ABCDEF"
        self.sLen = len(self.peptide)
        self.ions = {}
        
    def plotPeptide(self, pep, ions):
        self.peptide = pep
        self.sLen = len(pep)
        self.revTypes = set([])
        self.ions = ions
        self.Update()
        
    def paintEvent(self, event):
        if not self.peptide:
            return
        qp = QPainter()
        w,h = self.width(), self.height()
        qp.begin(self)
        qp.setPen(QColor(255,0,0))
        font = QFont()
        fm = QFontMetrics(font)
        tpos=0
        lx=0
        ly=0
        tw=0
        th=0
        fsize = 1
        for i in self.peptide:
            tsize = fm.size(Qt.TextSingleLine,i)
            pw = tsize.width()
            ph = tsize.height()
            if pw > tw:
                tw = pw
                pp = i
            if ph > th:
                th = pw
                pp = i
        isize=1
        font.setPointSize(isize)
        fm = QFontMetrics(font)
        tsize = fm.size(Qt.TextSingleLine,"10")
        tw, th = tsize.width(),tsize.height()
        while th < h*1/10:
            isize+=1
            font.setPointSize(isize)
            fm = QFontMetrics(font)
            tsize = fm.size(Qt.TextSingleLine,"10")
            tw,th = tsize.width(),tsize.height()
        ih=th
        fsize=1
        while (th*1.25+ih*2 <= h*1/2) and (tw*len(self.peptide) < w*1/2):
            font.setPointSize(fsize)
            fm = QFontMetrics(font)
            try:
                tsize = fm.size(Qt.TextSingleLine,pp)
                tw,th = tsize.width(),tsize.height()
            except IndexError:
                break
            fsize+=1
        totalw = len(self.peptide)*tw
        font.setPointSize(isize)
        fm = QFontMetrics(font)
        tsize = fm.size(Qt.TextSingleLine,"x")
        ixw,ixh = tsize.width(),tsize.height()
        tsize = fm.size(Qt.TextSingleLine,"y")
        iyw,iyh = tsize.width(),tsize.height()
        tsize = fm.size(Qt.TextSingleLine,"z")
        izw,izh = tsize.width(),tsize.height()
        tsize = fm.size(Qt.TextSingleLine,"a")
        iaw,iah = tsize.width(),tsize.height()
        tsize = fm.size(Qt.TextSingleLine,"b")
        ibw,ibh = tsize.width(),tsize.height()
        tsize = fm.size(Qt.TextSingleLine,"c")
        icw,ich = tsize.width(),tsize.height()
        cy = h/2
        lx = w/2-totalw/2
        toDraw = {}
        sLen = self.sLen
        for entry in self.ions:
            fType = entry[2]
            index = entry[3]
            try:
                toDraw[index].append(fType)
            except KeyError:
                toDraw[index] = [fType]
        for i,v in enumerate(self.peptide):
            font.setPointSize(fsize)
            qp.setFont(font)
            fm = QFontMetrics(font)
            tsize = fm.size(Qt.TextSingleLine,v)
            tw,th = tsize.width(),tsize.height()
            qp.setPen(QColor((0,0,0)))
            qp.drawText(QPoint(lx-tw/2,cy), v)
            if i==sLen-1:
                break
            font.setPointSize(isize)
            qp.setFont(font)
            qp.drawLine(QLineF(lx+tw/2,cy-th/2-ixh/2-iyh-izh,lx+tw/2,cy+th/2+iah+ibh+ich/2))
            qp.setPen(QColor((255,0,255)))
            ind = str(sLen-i-1)
            fm = QFontMetrics(font)
            tsize = fm.size(Qt.TextSingleLine,ind)
            iw,ih = tsize.width(),tsize.height()
            qp.drawText(QPoint(lx+tw/2,cy-ixh/2-iyh-izh-ih),ind)
            qp.drawText(QPoint(lx+tw/2,cy+th/2+iah+ibh+ich/2),str(i+1))
            try:
                todraw = toDraw[i+1]
                for fragType in todraw:
                    if fragType == 'a':
                        qp.setPen(QColor((255,0,0)))
                        ay = cy+th/2
                        qp.drawLine(QLineF(lx,ay,lx+tw/2,ay))
                        qp.setPen(QColor(255,0,0))
                        qp.drawText(QPoint(lx-iw,ay-iah/2),"a")
                    elif fragType == 'b':
                        by = cy+th/2+iah
                        qp.setPen(QColor((255,0,0)))
                        qp.drawLine(QLineF(lx,by,lx+tw/2,by))
                        qp.setPen(QColor(255,0,0))
                        qp.drawText(QPoint(lx-iw,by-ibh/2), "b")
                    elif fragType == 'c':
                        cy2 = cy+th/2+iah+ibh
                        qp.setPen(QColor((255,0,0)))
                        qp.drawLine(QLineF(lx,cy2,lx+tw/2,cy2))
                        qp.setPen(QColor(255,0,0))
                        qp.drawText(QPoint(lx-iw,cy2-ich/2), "c")
            except KeyError:
                pass
            try:
                todraw = toDraw[int(ind)]
                for fragType in todraw:
                    if fragType == 'x':
                        #we're x/y/z ions, we go in reverse
                        xy = cy-th/2-izh-iyh
                        qp.setPen(QColor((0,0,255)))
                        qp.drawLine(QLineF(lx+tw,xy,lx+tw/2,xy))
                        qp.setPen(QColor(0,0,255))
                        qp.drawText(QPoint(lx+tw,xy-ixh/2),"x")
                    elif fragType == 'y':
                        #we're x/y/z ions, we go in reverse
                        yy = cy-th/2-ixh
                        qp.setPen(QColor((0,0,255)))
                        qp.drawLine(QLineF(lx+tw,yy,lx+tw/2,yy))
                        qp.setPen(QColor(0,0,255))
                        qp.drawText(QPoint(lx+tw,yy-iyh/2), "y")
                    elif fragType == 'z':
                        #we're x/y/z ions, we go in reverse
                        zy = cy-th/2
                        qp.setPen(QColor((0,0,255)))
                        qp.drawLine(QLineF(lx+tw,zy,lx+tw/2,zy))
                        qp.setPen(QColor(0,0,255))
                        qp.drawText(QPoint(lx+tw,zy-izh/2), "z")
            except KeyError:
                pass
            lx+=tw+5
        qp.end() 

class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setWindowTitle("Spectrum Viewer")
        self._form = SpecView.Ui_MainWindow()
        self._form.setupUi(self)
        self._form.retranslateUi(self)
        self.fileType = ""
        splitter = QSplitter()
        self._form.verticalLayout_3.addWidget(splitter)
        splitter.setOrientation(Qt.Vertical)
        splitter.setParent(self._form.tab)
        self.peptidePanel = PeptidePanel(self)
        self.draw = DrawFrame(self)
        splitter.addWidget(self.peptidePanel)
        splitter.addWidget(self.draw)
        self.pepFiles = {}
        self.mgfFiles = {}
        self.gffFiles = {}
        #tree view
        self.tree = QTreeWidget()
        splitter.addWidget(self.tree)
        self.show()
        
    def addPage(self, path):
        self.data = collections.OrderedDict()
        self.objMap = {}
        pepSpecType = ('gff', 'xml')
        specType = ('mgf',)
        if [i for i in pepSpecType if i in path]:
            self.fileType = [i for i in pepSpecType if i in path][0]#this changes based on the file input somewhat
            colnames = ["Scan Title", "Peptide", "Modifications", "Charge", "Accession"]
            self.tree.setColumnCount(len(colnames))
            self.tree.setHeaderLabels(QStringList(colnames))
            self.dataType = 'pepspectra'
            if self.fileType == 'gff':
                gffIterator = fileIterators.GFFIterator(path, random=['SpectraId1', 'SpectraId2'])
                self.parent.gffFiles[path] = gffIterator
                for i in gffIterator:
                    if not i:
                        continue
                    sid = i.getAttribute('SpectraId1')
                    if sid:
                        toAdd = [sid,i.getAttribute('Sequence'),i.getAttribute('Modifications'), i.getAttribute('Charge'), i.getAttribute('Name')]
                        nid = toAdd[1]#group on peptide by default
                        node = self.objMap.get(nid)
                        if node is None:
                            newNode = QTreeWidgetItem(QStringList(toAdd))
                            newNode.subnodes = []
                            self.data[newNode] = newNode
                            self.objMap[nid] = newNode
                        else:
                            newNode = QTreeWidgetItem(QStringList(toAdd))
                            self.data[node].subnodes.append(newNode)
            elif self.fileType == 'xml':
                pepParser = fileIterators.spectraXML(path)
                self.pepFiles[path] = pepParser
                for i in pepParser.getScans():
                    toAdd = [i.getId(), i.getPeptide(), i.getModifications(), i.getCharge(),i.getAccession()]
                    nid = toAdd[1]#group on peptide by default
                    node = self.objMap.get(nid)
                    if node is None:
                        newNode = QTreeWidgetItem(QStringList(toAdd))
                        newNode.subnodes = []
                        self.data[newNode] = newNode
                        self.objMap[nid] = newNode
                    else:
                        newNode = QTreeWidgetItem(QStringList(toAdd))
                        self.data[node].subnodes.append(newNode)
        elif [i for i in specType if i in path]:
            self.fileType = 'spectra'#these are all generic more or less, so spectra works
            self.dataType = 'spectra'
            colnames = ["Scan Title", "Charge", "RT", "Precursor Mass"]
            mgf = fileIterators.mgfIterator(path, random=True)
            self.mgfFiles[path] = mgf
            for i in mgf:
                if not i:
                    continue
                toAdd = [i.getTitle(),i.getCharge(),i.getRT(), i.getPrecursor()]
                nid = toAdd[0]#group on title by default (should be unique)
                newNode = QTreeWidgetItem(QStringList(toAdd))
                newNode.subnodes = []
                node = self.objMap.get(nid)
                if node is None:
                    self.data[newNode] = newNode
                    self.objMap[nid] = newNode
                else:
                    self.data[node].subnodes.append(newNode)
        else:
            colnames = ["none"]
            self.fileType = 'none'
            self.dataType = 'none'
        self.tree.setColumnCount(len(colnames))
        self.tree.setHeaderLabels(colnames)
        if self.data:
            self.tree.addTopLevelItems(self.data.keys())
            for i in self.data:
                i.addChildren(self.data[i].subnodes)
                
    def Group(self, col):
        newData = collections.OrderedDict()
        objMap = {}
        for i in self.data:
            i.subnodes.insert(0,i)
            for subnode in i.subnodes:
                toAdd = subnode.data
                nid = toAdd[col]#group on peptide by default
                newNode = QTreeWidgetItem(toAdd)
                newNode.subnodes = []
                node = objMap.get(nid)
                if node is None:
                    newData[newNode] = newNode
                    objMap[nid] = newNode
                else:
                    newData[node].subnodes.append(newNode)
            i.subnodes.pop(0)
        self.data.clear()
        self.objMap.clear()
        self.data = newData
        self.objMap = objMap
        self.tree.addTopLevelItems(self.data.keys())
        for i in self.data:
            i.addChildren(self.data[i])
        
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
                it = self.gffFiles[path]
            except KeyError:
                loadGFF(path)
                it = self.gffFiles[path]
            gob = it.getGFF('SpectraId1',title)
            mods = gob.getAttribute('Modifications')
            fileName = titleParse.match(title).group(1)
            if not fileNames.has_key(fileName):
                return
            path = fileNames[fileName]
            try:
                scan = self.mgfFiles[path].getScan(title)
            except KeyError:
                mgf = fileIterators.mgfIterator(path)
                self.mgfFiles[path] = mgf
                try:
                    scan = mgf.getScan(title)
                except StopIteration:
                    scan = None
            if not scan:
                title = gob.getAttribute('SpectraId2')
                scan = self.mgfFiles[path].getScan(title)
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
            except KeyError:
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


class PlotPanel(QWidget):
    def __init__( self, parent):
        QWidget.__init__(self)
        self.dpi = 100
        self.figure = mpl.figure.Figure((3.0, 2.0), dpi = self.dpi)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(parent)
        self.parent = parent
        # initialize matplotlib stuff
        self.canvas.mpl_connect('button_press_event', self.onMouseDown)
        self.canvas.mpl_connect('button_release_event', self.onMouseRelease)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.vbox = QVBoxLayout(self)
        self.peptidePanel = self.parent.peptidePanel
        self.vbox.addWidget(self.toolbar)
        self.vbox.addWidget(self.canvas)
        self.canvas.setGeometry(QRect(10, 150, 490, 390))
        self.axes = self.figure.add_subplot(111)
        self.mouse = 0

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
        height = self.toolbar.height()
        idmap = {}
        for ion,desc in zip(('x','y','z','a','b','c'),('X ions', 'Y Ions', 'Z Ions', 'A Ions', 'B Ions', 'C Ions')):
#            
#            dc = QPainter()
#            im = QImage(QSize(35,height))
#            dc.begin(im)
#            font = QFont()
#            font.setPointSize(15)
#            dc.setFont(font)
#            dc.setPen(QColor(255,0,0))
#            dc.drawText(QPoint(7,3), ion)
#            dc.end()
            self.toolbar.addAction(ion)
            #self.toolbar.AddAction(ion, bitmap=im,shortHelp = 'Show %s'%desc)
#            if ion is 'y' or ion is 'b':
#                self.toolbar.ToggleTool(lid, True)
        #a little icon to see all spectra
#        im = QImage(QSize(30,height*3/4))
#        dc = QPainter()
#        dc.begin(im)
#        dc.setFont(font)
#        print self.toolbar.BackgroundColour
#        dc.setBrush(self.toolbar.BackgroundColour)
#        dc.drawRectangle(-1,-1,35,tsize)
#        qp.setPen(QColor(105,105,105))
#        for i in xrange(0,30,3):
#            rn = float('%0.3f'%random.uniform(0.35,1.0))
#            dc.drawRectangle(i,tsize*rn,2,tsize-tsize*rn)
#        lid = wx.NewId()
#        dc.end()
#        self.toolbar.AddCheckLabelTool(lid, "Spectra", bitmap=im, shortHelp="Plot Entire spectra")
#        self.Bind(wx.EVT_TOOL, self.onViewSpectra, id=lid)
#        self.Bind(wx.EVT_TOOL, self.onXIons, id=idmap['x'])
#        self.Bind(wx.EVT_TOOL, self.onYIons, id=idmap['y'])
#        self.Bind(wx.EVT_TOOL, self.onZIons, id=idmap['z'])
#        self.Bind(wx.EVT_TOOL, self.onAIons, id=idmap['a'])
#        self.Bind(wx.EVT_TOOL, self.onBIons, id=idmap['b'])
#        self.Bind(wx.EVT_TOOL, self.onCIons, id=idmap['c'])
        #error tolerance
        self.toolbar.addAction("unmatched")
        self.etolerance = QLineEdit("0.01", self.toolbar)
        self.etolerance.setToolTip("Mass Error Tolerance (da)")
        self.toolbar.addWidget(self.etolerance)
        
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
        for i in peaks:
            if not i:
                continue
            mz,inten,fragType, fragNum,charge,loss,aa = i
            if self.ionView[fragType]:
                self.points.append((mz,inten,fragType, fragNum,charge,loss,aa))
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
            self.points.append((x,y,'spectra', None,None,None,None))
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
            x,y,fragType, fragNum,charge,loss,aa = pt_list
            if fragType in blues:
                col=[0.0,0.0,1.0]
                tcolor='b'
            elif fragType == 'spectra':
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
                self.hitMapX[x].add((y,'m/z: %d Int: %d'%(x,y)))
            except KeyError:
                self.hitMapX[x] = set([(y,'m/z: %d Int: %d'%(x,y))])
            bar = self.subplot.bar(x,y,color=col,align='center')
            if fragType == 'spectra':
                bar[0].set_zorder(0)
            else:
                if not fragType:
                    continue
                hlen=len(self.hitMapX[x])*5
                if loss:
                    txt = '$\mathtt{%s%d^{%s^{%s}}}$'%(fragType,fragNum,loss,''.join(['+' for i in xrange(charge)]))
                else:
                    txt = '$\mathtt{%s%d^{%s}}$'%(fragType,fragNum,''.join(['+' for i in xrange(charge)]))
                self.subplot.text(x,y+hlen,txt,color=tcolor)
        self.subplot.axes.set_xbound(lower=0, upper=xmax+10)
        self.subplot.axes.set_ybound(lower=0, upper=ymax+20)
        self.canvas.draw()
        
    def onMouseMove(self,event):
        if self.mouse:
            return
        inv = self.subplot.axes.transData.inverted()
        x,y = inv.transform((event.x,event.y))
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
        txt = ""
        if hList:
            sorted(hList,key=operator.itemgetter(0))
            txt+='m/z = %0.3f (Int: %0.3f)\n'%(hList[0][1],hList[0][2][0])
            self.annotate = self.subplot.axes.annotate(txt,
            (x,y), xytext=(-2*20, 5), textcoords='offset points',
            bbox=self.bbox, arrowprops=self.arrowprops)
        self.canvas.draw()

app = QApplication(sys.argv)
w = MainWindow()
w.addPage('A1.2012_06_07_12_20_00.t.xml')
app.exec_()

#some dangerous globals for us
searchPaths = set([])
fileNames = {}
titleParse = re.compile('(.+?)\.\d+\.\d+\.\d+|\w+')



