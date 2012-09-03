"""
Author: Chris Mitchell (chris.mit7@gmail.com)
Copyright (C) 2012 Chris Mitchell

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import fileIterators, operator, os, re, ConfigParser, random, matplotlib as mpl, figureIons, collections, sys
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import numpy as np, time, thread
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import SpecView

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
        
class FileObject(object):
    def __init__(self, path):
        self.path = path

class LoaderThread(QThread):
    """
    thread to load files
    """
    def __init__(self, parent, files):
        QThread.__init__(self, parent)
        self.parent = parent
        self.files = files
        self.data = {}
        self.objMap = {}
        self.connect(self, SIGNAL('finished()'), self.parent.threadReturn)
        
    def run(self):
        self.its = []
        for path in self.files:
            iterObj = FileObject(path)
            self.its.append(iterObj)
            pepSpecType = ('gff', 'xml', 'msf')
            specType = ('mgf',)
            if [i for i in pepSpecType if i in path]:
                iterObj.fileType = [i for i in pepSpecType if i in path][0]#this changes based on the file input somewhat
                self.colnames = ["Scan Title", "Peptide", "Modifications", "Charge", "Accession"]
                self.groupBy = 1
                iterObj.dataType = 'pepspectra'
                if iterObj.fileType == 'gff':
                    iterObj.iObj = fileIterators.GFFIterator(path, random=['SpectraId1', 'SpectraId2'])
                    for i in iterObj.iObj:
                        if not i:
                            continue
                        self.emit(SIGNAL('updateProgress'),iterObj.iObj.getProgress())
                        sid = i.getAttribute('SpectraId1')
                        if sid:
                            toAdd = [sid,i.getAttribute('Sequence'),i.getAttribute('Modifications'), i.getAttribute('Charge'), i.getAttribute('Name')]
                            nid = toAdd[self.groupBy]#group on peptide by default
                            node = self.objMap.get(nid)
                            if node is None:
                                newNode = QTreeWidgetItem(QStringList(toAdd))
                                newNode.fileName = path
                                newNode.subnodes = []
                                newNode.data = toAdd
                                self.data[newNode] = newNode
                                self.objMap[nid] = newNode
                            else:
                                newNode = QTreeWidgetItem(QStringList(toAdd))
                                newNode.fileName = path
                                newNode.data = toAdd
                                self.data[node].subnodes.append(newNode)
                elif iterObj.fileType == 'xml' or iterObj.fileType == 'msf':
                    if iterObj.fileType == 'xml':
                        iterObj.iObj = fileIterators.XTandemXML(path)
                    elif iterObj.fileType == 'msf':
                        iterObj.iObj = fileIterators.ThermoMSFIterator(path)
                        self.colnames.append('Spectrum ID')
                        self.colnames.append('Confidence Level')
                        self.colnames.append('Search Engine Rank')
                    for i in iterObj.iObj:
                        self.emit(SIGNAL('updateProgress'),iterObj.iObj.getProgress())
                        if iterObj.fileType == 'msf':
                            toAdd = [i.getId(), i.getPeptide(), i.getModifications(), str(i.getCharge()),i.getAccession(),str(i.spectrumId), str(i.confidence), str(i.rank)]
                        else:
                            toAdd = [i.getId(), i.getPeptide(), i.getModifications(), str(i.getCharge()),i.getAccession()]
                        nid = toAdd[self.groupBy]#group on peptide by default
                        node = self.objMap.get(nid)
                        if node is None:
                            newNode = QTreeWidgetItem(QStringList(toAdd))
                            newNode.fileName = path
                            newNode.data = toAdd
                            newNode.subnodes = []
                            self.data[newNode] = newNode
                            self.objMap[nid] = newNode
                        else:
                            newNode = QTreeWidgetItem(QStringList(toAdd))
                            newNode.fileName = path
                            newNode.data = toAdd
                            self.data[node].subnodes.append(newNode)
            elif [i for i in specType if i in path]:
                iterObj.fileType = 'spectra'#these are all generic more or less, so spectra works
                iterObj.dataType = 'spectra'
                self.colnames = ["Scan Title", "Charge", "RT", "Precursor Mass"]
                iterObj.iObj = fileIterators.mgfIterator(path, random=True)
                self.groupBy = 0
                for i in iterObj.iObj:
                    if not i:
                        continue
                    self.emit(SIGNAL('updateProgress'),iterObj.iObj.getProgress())
                    toAdd = [i.getTitle(),i.getCharge(),i.getRT(), i.getPrecursor()]
                    nid = toAdd[self.groupBy]#group on title by default (should be unique)
                    newNode = QTreeWidgetItem(QStringList(toAdd))
                    newNode.fileName = path
                    newNode.data = toAdd
                    newNode.subnodes = []
                    node = self.objMap.get(nid)
                    if node is None:
                        self.data[newNode] = newNode
                        self.objMap[nid] = newNode
                    else:
                        self.data[node].subnodes.append(newNode)
            else:
                self.colnames = ["none"]
                iterObj.fileType = 'none'
                iterObj.dataType = 'none'

class PeptidePanel(QWidget):
    def __init__(self, parent):
        QWidget.__init__(self)
        self.parent = parent
        self.peptide = ""
        self.sLen = len(self.peptide)
        self.ions = {}
        
    def sizeHint(self):
        return QSize(self.width(), 100)
        
    def plotPeptide(self, pep, ions):
        self.peptide = pep
        self.sLen = len(pep)
        self.revTypes = set([])
        self.ions = ions
        self.update()
        
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
                        qp.drawText(QPoint(lx-iw/2,ay),"a")
                    elif fragType == 'b':
                        by = cy+th/2+iah
                        qp.setPen(QColor((255,0,0)))
                        qp.drawLine(QLineF(lx,by,lx+tw/2,by))
                        qp.setPen(QColor(255,0,0))
                        qp.drawText(QPoint(lx-iw/2,by), "b")
                    elif fragType == 'c':
                        cy2 = cy+th/2+iah+ibh
                        qp.setPen(QColor((255,0,0)))
                        qp.drawLine(QLineF(lx,cy2,lx+tw/2,cy2))
                        qp.setPen(QColor(255,0,0))
                        qp.drawText(QPoint(lx-iw/2,cy2), "c")
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
                        qp.drawText(QPoint(lx+tw,xy),"x")
                    elif fragType == 'y':
                        #we're x/y/z ions, we go in reverse
                        yy = cy-th/2-ixh
                        qp.setPen(QColor((0,0,255)))
                        qp.drawLine(QLineF(lx+tw,yy,lx+tw/2,yy))
                        qp.setPen(QColor(0,0,255))
                        qp.drawText(QPoint(lx+tw,yy), "y")
                    elif fragType == 'z':
                        #we're x/y/z ions, we go in reverse
                        zy = cy-th/2
                        qp.setPen(QColor((0,0,255)))
                        qp.drawLine(QLineF(lx+tw,zy,lx+tw/2,zy))
                        qp.setPen(QColor(0,0,255))
                        qp.drawText(QPoint(lx+tw,zy), "z")
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
        self.pepFiles = {}
        self.mgfFiles = {}
        self.gffFiles = {}
        self.tabs = []
        self.resize(1000,600)
        self.setAcceptDrops(True)
        #our menu
        menu = self.menuBar()
        self.filemenu = menu.addMenu('&File')
        exitAction = QAction('Exit', self)
        exitAction.triggered.connect(qApp.quit)
        openAction = QAction('&Open', self)
        openAction.triggered.connect(self.onOpen)
        self.filemenu.addAction(openAction)
        self.filemenu.addAction(exitAction)
        self.settings = menu.addMenu('&Settings')
        self.validExtensions = set(['.xml', '.msf', '.mgf'])
        self.show()
        
    def isValidFile(self,name):
        if os.path.splitext(name)[1].lower() in self.validExtensions:
            return True
        return False
        
    def onOpen(self):
        fileNames = list(QFileDialog.getOpenFileNames(self, 'Select File(s) To Open'))
        files = [str(fileName) for fileName in fileNames if self.isValidFile(fileName)]
        
    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()
        
    def dropEvent(self, event):
        if event.mimeData().hasUrls:
            event.setDropAction(Qt.CopyAction)
            event.accept()
            files = [str(fileName.toLocalFile()) for fileName in event.mimeData().urls() if self.isValidFile(str(fileName.toLocalFile()))]
            self.addPage(files)
        else:
            event.ignore()
    
    def addPage(self, files):
        if not self._form.tabWidget.count():
            self._form.instructions.deleteLater()
        vt = ViewerTab(self, files)
        self.tabs.append(vt)
        if len(files)>1:
            self._form.tabWidget.addTab(vt, '%d Merged Files'%len(files))
        else:
            self._form.tabWidget.addTab(vt, os.path.split(files[0])[1])
        self._form.tabWidget.setCurrentWidget(vt)
        vt.threadPage()
        
class ViewerTab(QWidget):
    def __init__(self, parent, files):
        QWidget.__init__(self)
        self.parent = parent
        self.data = collections.OrderedDict()
        self.objMap = {}
        self.fileType = ""
        self.files = files
        
    def onClick(self, item, col):
        self.item = item
        scanTitle = item.text(0)
        if self.getFileType(item.fileName) == 'msf':
            self.loadScan(item.fileName, str(scanTitle),specId=item.text(5),peptide=item.text(1))
        else:
            self.loadScan(item.fileName, str(scanTitle))
        
    def onHeaderRC(self, pos):
        gpos = self.mapToGlobal(pos)
        menu = QMenu()
        menu.addAction("Group By Column")
        ty = self.tree.pos().y()
        gpos.setY(gpos.y()+ty)
        selected = menu.exec_(self.mapToParent(gpos))
#        connect(selected, SIGNAL(triggered()), this, SLOT(Group()))
        if selected:
            self.Group(self.tree.header().logicalIndexAt(pos))
            
    def threadPage(self):
        self.progress = QProgressBar()
        self.vl = QVBoxLayout(self.parent._form.tabWidget.currentWidget())
        self.vl.setObjectName('vl')
        self.progressText = QLabel()
        self.progressText.setText("Loading Spectra File...MSF Files may take longer to start showing load status")
        self.vl.addWidget(self.progressText, 0, Qt.AlignTop)
        self.vl.addWidget(self.progress, 15)
        self.LoadThread = LoaderThread(self, self.files)
        self.connect(self.LoadThread, SIGNAL('updateProgress'), self.updateProgress)
        self.LoadThread.start()
        
    def updateProgress(self, perc):
        self.progress.setValue(perc)
    
    def threadReturn(self):
        #set up our page layout
        splitter = QSplitter()
        splitter.setChildrenCollapsible(True)
        self.vl.removeWidget(self.progress)
        self.vl.removeWidget(self.progressText)
        self.progress.deleteLater()
        self.progressText.deleteLater()
        #self.vl = QVBoxLayout(self.parent._form.tabWidget.currentWidget())
        #self.vl.setObjectName('vl')
        self.vl.addWidget(splitter)
        splitter.setOrientation(Qt.Vertical)
        splitter.setParent(self)
        self.peptidePanel = PeptidePanel(self)
        self.draw = DrawFrame(self)
        splitter.addWidget(self.peptidePanel)
        splitter.addWidget(self.draw)
        self.tree = QTreeWidget()
        splitter.addWidget(self.tree)
        self.tree.header().setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.header().customContextMenuRequested.connect(self.onHeaderRC)
        self.tree.itemDoubleClicked.connect(self.onClick)
#        self.dataType = self.LoadThread.dataType
#        self.fileType = self.LoadThread.fileType
#        self.groupBy = self.LoadThread.groupBy
        iterObjects = self.LoadThread.its
        self.colnames = self.LoadThread.colnames
        self.objMap = self.LoadThread.objMap
        self.data = self.LoadThread.data
        for it in iterObjects:
            if it.fileType == 'gff':
                self.parent.gffFiles[it.path] = it.iObj
            elif it.fileType == 'xml' or it.fileType == 'msf':
                self.parent.pepFiles[it.path] = it.iObj
            elif it.fileType == 'spectra':
                self.parent.mgfFiles[it.path] = it.iObj
        self.tree.setColumnCount(len(self.colnames))
        self.tree.setHeaderLabels(self.colnames)
        if self.data:
            self.tree.addTopLevelItems(self.data.keys())
            for i in self.data:
                i.addChildren(self.data[i].subnodes)
        self.tree.setSortingEnabled(True)
        
    def Group(self, col):
        self.tree.setSortingEnabled(False)
        newData = collections.OrderedDict()
        objMap = {}
        self.groupBy = col
        for i in self.data:
            i.subnodes.insert(0,i)
            for subnode in i.subnodes:
                toAdd = subnode.data
                nid = toAdd[col]#group on peptide by default
                newNode = QTreeWidgetItem(toAdd)
                newNode.subnodes = []
                node = objMap.get(nid)
                newNode.data = toAdd
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
        self.tree.clear()
        self.tree.addTopLevelItems(self.data.keys())
        for i in self.data:
            sn = self.data[i].subnodes
            if sn:
                i.addChildren(sn)
        self.tree.setSortingEnabled(True)
        
    def getTolerance(self):
        return float(self.draw.etolerance.text())
        
    def reloadScan(self):
        try:
            self.onClick(self.item, 0)
        except AttributeError:
            return
        
    def getFileType(self, path):
        if os.path.splitext(path)[1][1:] == 'mgf':
            return 'spectra'
        else:
            return os.path.splitext(path)[1][1:]
        
    def loadScan(self, fileO, title, **kwrds):
#        print fileO,title,kwrds
#        print self.parent.pepFiles, self.getFileType(fileO)
        path = fileO
        self.title=title
        self.kwrds = kwrds
        if self.getFileType(path) == 'gff':
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
        elif self.getFileType(path) == 'xml' or self.getFileType(path) == 'msf':
            if self.getFileType(path) == 'xml':
                scan = self.parent.pepFiles[path].getScan(title)
            elif self.getFileType(path) == 'msf':
                scan = self.parent.pepFiles[path].getScan(title,kwrds['specId'],kwrds['peptide'])
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
        elif self.getFileType(path) == 'spectra':
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
        self.ionView = {'x': False, 'y': True, 'z': False, 'a': False, 'b': True, 'c': False, 'all': False, '++': False, '>2': False}
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
        for ion,desc,index in zip(('x','y','z','a','b','c', '++', '>2'),('X ions', 'Y Ions', 'Z Ions', 'A Ions', 'B Ions', 'C Ions', 'Doubly Charged Ions', 'All Possible Fragments'),xrange(8)):
            qpixmap = QPixmap(35, height)
            qpixmap.fill(QColor(255,255,255))
            qp = QPainter()
            qp.begin(qpixmap)
            if index < 3:
                qp.setPen(QColor(0,0,255))
            else:
                qp.setPen(QColor(255,0,0))
            font = QFont()
            font.setPointSize(20)
            qp.setFont(font)
            qp.drawText(QPointF(0,height*3/4), ion)
            qp.end()
            icon = QIcon(qpixmap)
            a = self.toolbar.addAction(icon, ion)
            a.setToolTip('Show %s'%desc)
            a.setCheckable(True)
            if ion is 'y' or ion is 'b':
                a.setChecked(True)
        qpixmap = QPixmap(35, height)
        qpixmap.fill(QColor(255,255,255))
        qp = QPainter()
        qp.begin(qpixmap)
        for i,j in zip([1,5,7,9,10,17,20,27,35],[25,5,30,16,7,10,4,35,28]):
            qp.drawLine(i,height,i,height-j)
        qp.end()
        icon = QIcon(qpixmap)
        a = self.toolbar.addAction(icon, 'all')
        a.setCheckable(True)
        a.setToolTip('Show All Spectra')
        self.toolbar.actionTriggered.connect(self.onAction)
        #error tolerance
        self.etolerance = QLineEdit("0.01", self.toolbar)
        self.etolerance.setToolTip("Mass Error Tolerance (da)")
        self.etolerance.editingFinished.connect(self.onToleranceEdit)
        self.toolbar.addWidget(self.etolerance)
        
    def onToleranceEdit(self):
        self.parent.reloadScan()
        
    def onAction(self, action):
        txt = str(action.iconText())
        if txt in set(['x','y','z','a','b','c', 'all', '++', '>2']):
            self.ionView[txt] = action.isChecked()
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
        self.figure.axes[0].set_xlabel('M/Z')
        self.figure.axes[0].set_ylabel('Intensity')
        
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
#w.addPage(['A1.2012_06_07_12_20_00.t.xml'])
#w.addPage('sampleMgf.mgf')
app.exec_()