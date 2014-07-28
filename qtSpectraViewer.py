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

import operator, os, re, xml.etree.cElementTree as etree, random, figureIons, sys, masses, itertools, time, math
from collections import OrderedDict
import pythomics.proteomics.parsers as fileIterators
from pythomics.proteomics.structures import ScanObject, PeptideObject
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import pyqtgraph as pg
import SpecView, SettingsPanel
import numpy as np

#some globals for us
fileNames = {}
titleParse = re.compile('(.+?)\.\d+\.\d+\.\d+|\w+')
SCAN_MAP = {}
FILE_MAP = {}

def loadConfig():
    try:
        config = etree.ElementTree(file='spectraConf.xml')
    except IOError:#doesn't exist
        return
    root = config.getroot()
    losses = root.find('LossMasses')
    for child in losses:
        loss = []
        for mass, frags, display in zip(child.findall('Mass'),child.findall('Fragments'), child.findall('Display')):
            if mass.text and frags.text and display.text:
                loss.append((float(mass.text),tuple(frags.text.split(',')),display.text))
        loss = tuple(loss)
        masses.lossMasses[child.text] = loss


def saveConfig():
    root = etree.Element('SpectraViewerConfig')
    losses = etree.SubElement(root,'LossMasses')
    for aa in masses.lossMasses:
        info = masses.lossMasses[aa]
        nl = etree.SubElement(losses,'AminoAcid')
        nl.text = aa
        for mass,fragments,display in info:
            nm = etree.SubElement(nl,'Mass')
            nm.text = str(mass)
            nf = etree.SubElement(nl,'Fragments')
            nf.text = ','.join(fragments)
            nd = etree.SubElement(nl,'Display')
            nd.text = display
    elTree = etree.ElementTree(root)
    elTree.write('spectraConf.xml', xml_declaration=True)

# class SortableTreeWidgetItem(QTreeWidgetItem):
#     def __lt__(self, other):
#         return self < other

class FileObject(object):
    def __init__(self, path):
        self.path = path
        self.scanMap = {}
        self.titleRegex = None
        pepSpecType = ('pep.xml', 'xml', 'msf', 'dat')
        specType = ('mgf', 'mzml')
        self.groupBy = 1
        self.idColumn = 0
        if [i for i in pepSpecType if i in path.lower()]:
            fileType = [i for i in pepSpecType if i in path][0]#this changes based on the file input somewhat
            self.colnames = OrderedDict([("Scan Title",(str,'id',1)), ("Peptide",(str,'peptide',1)), ("Modifications",(str,'getModifications',0)), ("Charge",(int,'charge',1)), ("Accession",(str, 'acc',1))])
            if fileType == 'pep.xml':
                self.iObj = fileIterators.PepXMLIterator(path)
                self.colnames['Expect'] = (float,'expect',1)
            elif fileType == 'xml':
                self.iObj = fileIterators.XTandemXMLIterator(path)
                self.colnames['Expect'] = (float,'expect',1)
            elif fileType == 'dat':
                self.iObj = fileIterators.MascotDATIterator(path)
                self.colnames['Hit Id'] = (int,'hit',1)
                self.colnames['Rank'] = (int,'rank',1)
            elif fileType == 'msf':
                self.iObj = fileIterators.ThermoMSFIterator(path)
                self.colnames['Spectrum ID'] = (int, 'spectrumId',1)
                self.idColumn = 5
                self.colnames['Confidence Level'] = (int, 'confidence',1)
                self.colnames['Search Engine Rank'] = (int, 'rank',1)
                self.colnames["RT"] = (float,'rt',1)
        elif [i for i in specType if i in path.lower()]:
            fileType = [i for i in specType if i in path.lower()][0]
            if fileType == 'mzml':
                self.colnames = OrderedDict([("Scan Title",(str,'title',1)), ("MS Level",(int,'ms_level',1)), ("Charge",(int,'charge',1)), ("RT",(float,'rt',1)), ("Precursor Mass",(str,'mass',1))])
                self.iObj = fileIterators.MZMLIterator(path)
            else:
                self.colnames = OrderedDict([("Scan Title",(str,'title',1)), ("Charge",(int,'charge',1)), ("RT",(float,'rt',1)), ("Precursor Mass",(str,'mass',1))])
                self.iObj = fileIterators.MGFIterator(path, random=True)
            self.groupBy = 0
        else:
            self.colnames = ["none"]

    def getChromatogram(self):
        return self.iObj.getChromatogram()

    def getBaseTrace(self):
        return self.iObj.getChromatogram()

    def getScan(self, treeItem):
        try:
            scan = self.scanMap[treeItem]
        except:
            scanTitle = treeItem.data[self.idColumn]
            if self.titleRegex:
                title = self.titleRegex.search(scanTitle)
                if title:
                    scanTitle = title.group(1)
            peptide = treeItem.data[1]
            # print scanTitle, peptide, type(self.iObj)
            scan = self.iObj.getScan(scanTitle, peptide=peptide)
            # convert the scans to numpy arrays
            scan.scans = np.array(scan.scans)
            # try:
            #     scan.ms1_scan.scans = np.array(scan.ms1_scan.scans)
            # except:
            #     pass
            self.scanMap[treeItem] = scan
            # print scan
        return scan

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
        for fileindex,path in enumerate(self.files):
            iterObj = FileObject(path)
            if path not in FILE_MAP:
                FILE_MAP[path] = iterObj
            self.its.append(iterObj)
            self.colnames = iterObj.colnames
            if not self.colnames:
                self.emit(SIGNAL('fileDone'),fileindex,len(self.files))
                return
            self.groupBy = iterObj.groupBy
            added = set([])
            for i in iterObj.iObj:
                if not i or i in added:
                    continue
                added.add(i)
                try:
                    self.emit(SIGNAL('updateProgress'), iterObj.iObj.getProgress())
                except AttributeError:
                    pass
                toAdd = tuple(getattr(i,j[1]) if j[2] else getattr(i,j[1])() for j in self.colnames.values())
                if toAdd in added:
                    continue
                added.add(toAdd)
                nid = toAdd[self.groupBy]#group on peptide by default
                node = self.objMap.get(nid)
                if node is None:
                    newNode = QTreeWidgetItem()
                    [newNode.setText(k,str(v)) for k,v in enumerate(toAdd)]
                    newNode.fileName = path
                    newNode.data = toAdd
                    newNode.subnodes = []
                    self.data[newNode] = newNode
                    self.objMap[nid] = newNode
                else:
                    newNode = QTreeWidgetItem(node)
                    [newNode.setText(k,str(v)) for k,v in enumerate(toAdd)]
                    newNode.fileName = path
                    newNode.data = toAdd
                    self.data[node].subnodes.append(newNode)
            self.emit(SIGNAL('fileDone'),fileindex,len(self.files))

class PeptidePanel(QDockWidget):
    def __init__(self, parent):
        super(PeptidePanel, self).__init__()
        self.setWindowTitle("Peptide")
        self.parent = parent
        self.peptide = ""
        self.sLen = len(self.peptide)
        self.ions = []
        self.errors = []
        self.errorboxes = []
        self.setMouseTracking(True)

    def mouseMoveEvent(self, event):
        x,y = int(event.x()),int(event.y())
        for i in self.errorboxes:
            if i.contains(x,y):
                QToolTip.showText(self.mapToGlobal(i.center().toPoint()), i.ion)

    def sizeHint(self):
        return QSize(self.width(), 100)

    def plotPeptide(self, pep, ions):
        self.peptide = pep
        self.sLen = len(pep)
        self.revTypes = set([])
        self.ions = [i[0] for i in ions]
        self.errors = [i[1] for i in ions]
        self.update()

    def paintEvent(self, event):
        """
        3/4 of left side is devoted to peptide, 1/4 of right is error measurement
        """
        if not self.peptide:
            return
        qp = QPainter()
        w,h = self.width()*9/10, self.height()
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
        for entry,error in zip(self.ions,self.errors):
            fType = entry[2]
            index = entry[3]
            idesc = '%s%d%s'%(fType,index,entry[5]) if entry[5] else '%s%d'%(fType,index)
            try:
                toDraw[index].append((fType,error,idesc))
            except KeyError:
                toDraw[index] = [(fType,error,idesc)]

        #draw our error line on the right (a black line up and down)
        esize = self.width()/20
        #"normalize" our error to the width so a delta of 0.0001 actually has pixels
        escale = esize/self.parent.getTolerance()
        elx = self.width()-2*esize
        ely = (self.height()*9/10)/len(self.peptide)
        qp.drawLine(QLineF(elx,0,elx,h))
        #draw some axes for the top
        qp.drawLine(QLineF(elx-esize,0,elx+esize,0))
        font.setPointSize(isize)
        qp.setFont(font)
        tsize = fm.size(Qt.TextSingleLine,str(self.parent.getTolerance()))
        tw,th = tsize.width(),tsize.height()
        qp.drawText(QPointF(elx-esize-tw,th),str(self.parent.getTolerance()))
        qp.drawText(QPointF(elx+esize,th),str(self.parent.getTolerance()))
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
                for fragType,error,errortype in todraw:
                    esize = error*escale
                    if fragType == 'a':
                        qp.setPen(QColor((255,0,0)))
                        ay = cy+th/2
                        qp.drawLine(QLineF(lx,ay,lx+tw/2,ay))
                        qp.setPen(QColor(255,0,0))
                        errorBox = QRectF(elx-esize,ely*(i+1),3,3)
                        qp.drawRect(errorBox)
                        errorBox.ion = errortype
                        self.errorboxes.append(errorBox)
                        qp.drawText(QPoint(lx-iw/2,ay),"a")
                    elif fragType == 'b':
                        by = cy+th/2+iah
                        qp.setPen(QColor((255,0,0)))
                        qp.drawLine(QLineF(lx,by,lx+tw/2,by))
                        qp.setPen(QColor(255,0,0))
                        qp.drawText(QPoint(lx-iw/2,by), "b")
                        errorBox = QRectF(elx-esize,ely*(i+1),3,3)
                        qp.drawRect(errorBox)
                        errorBox.ion = errortype
                        self.errorboxes.append(errorBox)
                    elif fragType == 'c':
                        cy2 = cy+th/2+iah+ibh
                        qp.setPen(QColor((255,0,0)))
                        qp.drawLine(QLineF(lx,cy2,lx+tw/2,cy2))
                        qp.setPen(QColor(255,0,0))
                        qp.drawText(QPoint(lx-iw/2,cy2), "c")
                        errorBox = QRectF(elx-esize,ely*(i+1),3,3)
                        qp.drawRect(errorBox)
                        errorBox.ion = errortype
                        self.errorboxes.append(errorBox)
            except KeyError:
                pass
            try:
                todraw = toDraw[int(ind)]
                for fragType,error,errortype in todraw:
                    esize = error*escale
                    if fragType == 'x':
                        #we're x/y/z ions, we go in reverse
                        xy = cy-th/2-izh-iyh
                        qp.setPen(QColor((0,0,255)))
                        qp.drawLine(QLineF(lx+tw,xy,lx+tw/2,xy))
                        qp.setPen(QColor(0,0,255))
                        qp.drawText(QPoint(lx+tw,xy),"x")
                        errorBox = QRectF(elx-esize,ely*(i+1),3,3)
                        qp.drawRect(errorBox)
                        errorBox.ion = errortype
                        self.errorboxes.append(errorBox)
                    elif fragType == 'y':
                        #we're x/y/z ions, we go in reverse
                        yy = cy-th/2-ixh
                        qp.setPen(QColor((0,0,255)))
                        qp.drawLine(QLineF(lx+tw,yy,lx+tw/2,yy))
                        qp.setPen(QColor(0,0,255))
                        qp.drawText(QPoint(lx+tw,yy), "y")
                        errorBox = QRectF(elx-esize,ely*(i+1),3,3)
                        qp.drawRect(errorBox)
                        errorBox.ion = errortype
                        self.errorboxes.append(errorBox)
                    elif fragType == 'z':
                        #we're x/y/z ions, we go in reverse
                        zy = cy-th/2
                        qp.setPen(QColor((0,0,255)))
                        qp.drawLine(QLineF(lx+tw,zy,lx+tw/2,zy))
                        qp.setPen(QColor(0,0,255))
                        qp.drawText(QPoint(lx+tw,zy), "z")
                        errorBox = QRectF(elx-esize,ely*(i+1),3,3)
                        qp.drawRect(errorBox)
                        errorBox.ion = errortype
                        self.errorboxes.append(errorBox)
            except KeyError:
                pass
            lx+=tw+5
        qp.end()

class CustomClick(object):
    def mouseClickEvent(self, ev):
        pass


class AUCText(pg.TextItem, CustomClick):
    def __init__(self, *args, **kwargs):
        super(AUCText, self).__init__(*args, **kwargs)
        self.setAcceptHoverEvents(True)
        self.text = kwargs.get('text', '')
        self.color = kwargs.get('color', [1,1,1])
        self.contextMenu = [QAction("Copy Value", self)]
        self.contextMenu[0].triggered.connect(self.copyText)

    def copyText(self):
        print 'do it'

    def hoverEnterEvent(self, ev):
        self.setText(text=self.text, color='r')
        ev.ignore()

    def hoverLeaveEvent(self, ev):
        self.setText(text=self.text, color=self.color)
        ev.ignore()

    def mouseClickEvent(self, ev):
        super(AUCText, self).mouseClickEvent(ev)
        if ev.double():
            self.line.hide()
            self.hide()



class ChromatogramPanel(pg.PlotWidget):
    def __init__(self, *args, **kwargs):
        super(ChromatogramPanel, self).__init__(*args, **kwargs)
        self.setBackground([255, 255, 255, 255])
        self.chromatogram = None
        self.controlDown = False
        self.startPos = False
        self.appended = False
        self.plot = pg.PlotDataItem()
        self.vb = self.getViewBox()
        self.addItem(self.plot)
        self.sceneObj.sigMouseMoved.connect(self.onMouseMove)
        self.sceneObj.sigMouseClicked.connect(self.onMouseClick)

    def keyPressEvent(self, event):
        mods = event.modifiers()
        if mods == Qt.ControlModifier:
            self.controlDown = True
        else:
            self.controlDown = False
        event.ignore()

    def plotChromatogram(self, **kwargs):
        x = np.array(kwargs.get('x', []))
        y = np.array(kwargs.pop('y', []))
        append = kwargs.pop('append', False)
        print 'append status',append
        if x.any() and y.any():
            # self.clear()
            # FIXME: make a failsafe against different x values
            if append:
                if self.appended:
                    y = self.y+y
                else:
                    y = y
                    self.appended = True
            self.plot.setData(x, y, pen=[0,0,0])
            self.x, self.y = x,y
            return
        self.appended = False
        if self.chromatogram is None:
            return
        # self.clear()
        self.x = self.chromatogram.times
        self.y = self.chromatogram.intensities
        self.plot.setData(self.x, self.y, pen=[0,0,0])
        # print self.chromatogram.times

    def integrate(self, y_vals, h):
        """
        Integration using simpson's rule
        http://stackoverflow.com/questions/13320262/calculating-the-area-under-a-curve-given-a-set-of-coordinates-without-knowing-t
        """
        i=1
        total=y_vals[0]+y_vals[-1]
        for y in y_vals[1:-1]:
            if i%2 == 0:
                total+=2*y
            else:
                total+=4*y
            i+=1
        return total*(h/3.0)

    def onMouseClick(self, click):
        if click.double():
            self.plotChromatogram()
        view_click = click.pos()
        vcx, vcy = view_click.x(), view_click.y()
        clicks = []
        for i in self.items():
            if not issubclass(type(i), CustomClick):
                continue
            pos = i.pos()
            pos = self.vb.mapFromView(pos)
            br = i.boundingRect()
            br_w, br_h = br.width(), br.height()
            x = pos.x()
            if x < vcx < x+br_w:
                y = pos.y()
                if y < vcy < y+br_h:
                    clicks.append(i)
        for i in clicks:
            i.mouseClickEvent(click)
        if self.controlDown:
            if self.startPos:
                #we have a start position and clicked again, so we're ending our line
                self.controlDown = False
                start_xy = self.mapToView(self.startPos.scenePos())
                end_xy = self.mapToView(click.scenePos())
                #we don't care about y
                x1,y1 = start_xy.x(), start_xy.y()
                x2, y2 = end_xy.x(), end_xy.y()
                if x1 > x2:
                    x1, x2, y1, y2 = x2, x1, y2, y1
                self.startPos = False
                # get the area
                x = []
                y = []
                overlay_y = []
                slope = (y2-y1)/(x2-x1)
                b = y1-slope*x1
                for i,j in zip(self.x, self.y):
                    if x1 <= i <= x2:
                        x.append(i)
                        y_pos = slope*i+b
                        overlay_y.append(j)
                        if y_pos < 0:
                            y_pos = 0
                        if j >= y_pos:
                            y.append(j-y_pos)
                        else:
                            y.append(0)
                # average delta
                deltas = [a - x[i - 1] for i, a in enumerate(x)][1:]
                gap = sum(deltas)/len(deltas)
                if not gap:
                    area = sum(y)
                else:
                    area = self.integrate(y, gap)
                txt='AUC:%s'%area
                self.annotate = AUCText(text = txt, color=[1,1,1])
                self.annotate.line = self.chromLine
                self.addItem(self.annotate)
                self.annotate.setZValue(100)
                self.annotate.setPos(end_xy.x(), end_xy.y())
                # fill the curve
                # self.overlay = pg.PlotDataItem()
                # self.overlay.setData(x, overlay_y, pen=[1,1,1])
                # self.addItem(self.overlay)
                # self.fill = pg.FillBetweenItem(self.overlay, self.chromLine, brush='k')
                # self.fill.setZValue(100)
                # self.addItem(self.fill)
            else:
                self.chromLine = pg.PlotDataItem()
                self.startPos = click
                self.addItem(self.chromLine)

    def onMouseMove(self,event):
        if self.controlDown and self.startPos:
            # draw a line from start to here
            plot_xy = self.mapToView(self.startPos.scenePos())
            x = [plot_xy.x()]
            y = [plot_xy.y()]
            plot_xy = self.mapToView(event)
            x.append(plot_xy.x())
            y.append(plot_xy.y())
            self.chromLine.setData(x, y, pen=[1,1,1])

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
        settingsAction = QAction('&Settings', self)
        settingsAction.triggered.connect(self.onSettings)
        self.filemenu.addAction(openAction)
        self.filemenu.addAction(settingsAction)
        self.filemenu.addAction(exitAction)
        self.validExtensions = set(['.xml', '.msf', '.mgf', '.dat', '.mzml', '.gz'])
        self._form.tabWidget.setTabsClosable(True)
        self._form.tabWidget.tabCloseRequested.connect(self.onTabClose)
        loadConfig()
        self.connect(self,SIGNAL('triggered()'),self.closeEvent)
        self.show()


    def closeEvent(self, event):
        saveConfig()

    def settingsSave(self):
        global dError
        txt = self.dlg.defaultError.text()
        if txt:
            try:
                dError = int(txt)
            except ValueError:
                pass

    def addLoss(self):
        aa = str(self.dlg.aaBox.currentText())
        ions = []
        if self.dlg.checkBoxa.checkState() == Qt.Checked: ions.append('a')
        if self.dlg.checkBoxb.checkState() == Qt.Checked: ions.append('b')
        if self.dlg.checkBoxc.checkState() == Qt.Checked: ions.append('c')
        if self.dlg.checkBoxx.checkState() == Qt.Checked: ions.append('x')
        if self.dlg.checkBoxy.checkState() == Qt.Checked: ions.append('y')
        if self.dlg.checkBoxz.checkState() == Qt.Checked: ions.append('z')
        loss = float(self.dlg.massLoss.text())
        display = str(self.dlg.displayName.text())
        if masses.lossMasses.has_key(aa):
            closses = list(masses.lossMasses[aa])
        else:
            closses = set([])
        if not closses.count((loss,tuple(ions),display)):
            closses.append((loss,tuple(ions),display))
            masses.lossMasses[aa] = tuple(closses)
        # print masses.lossMasses[aa]


    def onSettings(self, args):
        self.dtmp = QDialog()
        self.dlg = SettingsPanel.Ui_Dialog()
        self.dlg.setupUi(self.dtmp)
        self.connect(self.dlg.buttonBox, SIGNAL('accepted()'),self.settingsSave)
        self.dlg.addLoss.pressed.connect(self.addLoss)
        self.dtmp.show()

    def isValidFile(self,name):
        if os.path.splitext(str(name))[1].lower() in self.validExtensions:
            return True
        return False

    def onOpen(self):
        fileNames = list(QFileDialog.getOpenFileNames(self, 'Select File(s) To Open'))
        files = [str(fileName) for fileName in fileNames if self.isValidFile(fileName)]
        self.addPage(files)

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
            if (event.keyboardModifiers()&Qt.ControlModifier):
                self.appendPage(files)
            else:
                self.addPage(files)
        else:
            event.ignore()

    def appendPage(self, files):
        #get current tab
        self._form.tabWidget.currentWidget().appendFiles(files)

    def addPage(self, files):
        if not self._form.instructions.isHidden():
            self._form.instructions.setHidden(True)
        vt = ViewerTab(self, files)
        self.tabs.append(vt)
        if len(files)>1:
            self._form.tabWidget.addTab(vt, '%d Merged Files'%len(files))
        elif files:
            self._form.tabWidget.addTab(vt, os.path.split(files[0])[1])
        self._form.tabWidget.setCurrentWidget(vt)
        vt.threadPage()

    def onTabClose(self, tab):
        self.tabs[tab].prepareDestruct()
        self._form.tabWidget.removeTab(tab)
        self.tabs.pop(tab)
        if not self._form.tabWidget.count():
            self._form.instructions.setHidden(False)

class LossPanel(QWidget):
    """
    a panel to add neutral losses
    """
    def __init__(self, parent):
        pass

class ViewerTab(QMainWindow):
    def __init__(self, parent, files):
        super(ViewerTab, self).__init__()
        self.parent = parent
        self.data = {}
        self.objMap = {}
        self.skipLosses = {}
        self.fileType = ""
        self.files = files

    def prepareDestruct(self):
        del self.data
        del self.objMap
#        self.deleteLater()

    def onClick(self, item):
        self.item = item
        self.loadScan(FILE_MAP[item.fileName].getScan(item))

    def onHeaderRC(self, pos):
        gpos = self.mapToGlobal(pos)
        menu = QMenu()
        menu.addAction("Group By Column")
        ty = self.tree.pos().y()
        gpos.setY(gpos.y()+ty+self.tree.parent().pos().y())
        selected = menu.exec_(gpos)#self.mapToParent(gpos))
        if selected:
            self.Group(self.tree.header().logicalIndexAt(pos))

    def appendFiles(self, files):
        self.LoadThread = LoaderThread(self, files)
        self.files+=files
        self.LoadThread.data = self.data
        self.LoadThread.objMap = self.objMap
        self.progress.show()
        self.progressText.show()
        self.onFileText.show()
        self.connect(self.LoadThread, SIGNAL('updateProgress'), self.updateProgress)
        self.connect(self.LoadThread, SIGNAL('fileDone'), self.fileDone)
        self.LoadThread.start()

    def threadPage(self):
        self.progress = QProgressBar()
        self.vl = QVBoxLayout(self.parent._form.tabWidget.currentWidget())
        self.vl.setObjectName('vl')
        self.progressText = QLabel()
        self.progressText.setText("Loading Spectra File...MSF & DAT Files may take longer to start showing load status")
        self.onFileText = QLabel()
        self.onFileText.setText("Loading file 1 of %d"%len(self.files))
        self.vl.addWidget(self.progressText, 0, Qt.AlignTop)
        self.vl.addWidget(self.onFileText, 0, Qt.AlignTop)
        self.vl.addWidget(self.progress, 15)
        self.LoadThread = LoaderThread(self, self.files)
        self.connect(self.LoadThread, SIGNAL('updateProgress'), self.updateProgress)
        self.connect(self.LoadThread, SIGNAL('fileDone'), self.fileDone)
        self.LoadThread.start()

    def updateProgress(self, perc):
        self.progress.setValue(perc)

    def fileDone(self, i,j):
        self.onFileText.setText('Loading file %d of %d'%(i+1,j))

    def threadReturn(self):
        #set up our page layout
        # if hasattr(self, 'splitter'):
        #     iterObjects = self.LoadThread.its
        #     self.data = self.LoadThread.data
        #     self.objMap = self.LoadThread.objMap
        #     for it in iterObjects:
        #         if it.fileType == 'gff':
        #             self.parent.gffFiles[it.path] = it.iObj
        #         elif it.dataType == 'pepspectra' or it.fileType == 'xml' or it.fileType == 'msf' or it.fileType == 'dat':
        #             self.parent.pepFiles[it.path] = it.iObj
        #         elif it.fileType == 'spectra':
        #             self.parent.mgfFiles[it.path] = it.iObj
        #     if self.data:
        #         for i in self.data.keys():
        #             if not i.treeWidget():
        #                 self.tree.addTopLevelItem(i)
        #                 i.addChildren(self.data[i].subnodes)
        #     self.tree.setSortingEnabled(True)
        #     self.onFileText.hide()
        #     self.progress.hide()
        #     self.progressText.hide()
        #     return
        # self.splitter = QSplitter()
        # self.splitter.setChildrenCollapsible(True)
        self.vl.removeWidget(self.progress)
        self.vl.removeWidget(self.progressText)
        # self.vl.addWidget(self.splitter)
        self.progress.hide()
        self.progressText.hide()
        self.onFileText.hide()
        self.vl.addWidget(self.onFileText)
        self.vl.addWidget(self.progressText)
        self.vl.addWidget(self.progress)
        # self.splitter.setOrientation(Qt.Vertical)
        # self.splitter.setParent(self)
        # self.peptidePanelDock = QDockWidget()
        self.peptidePanel = PeptidePanel(self)
        # self.peptidePanelDock.setWidget(self.peptidePanel)
        self.chromaPanel = ChromatogramPanel(self)
        self.baseTracePanel = ChromatogramPanel(self)
        self.draw = DrawFrame(self)
        self.ms1_draw = DrawFrame(self)
        self.searchBox = QLineEdit()
        self.searchBox.editingFinished.connect(self.onSearch)
        #search button
        self.sbutton = QToolButton(self.searchBox)
        self.sicon = QPixmap('icons/search.ico')
        self.sbutton.setIcon(QIcon(self.sicon))
        self.sbutton.setIconSize(self.sicon.size())
        self.sbutton.setCursor(Qt.ArrowCursor)
        self.sbutton.setStyleSheet("border: none; padding: 0px; ")
        self.searchBox.setStyleSheet("padding-left: %dpx;"%(int(self.sicon.width())))
        self.searchCols = QComboBox()
        self.searchGroup = QWidget()
        #filter box
        self.filterBox = QLineEdit()
        self.filterBox.editingFinished.connect(self.onFilter)
        self.filterList = OrderedDict()
        self.filterBoxes = []
        #filter button
        self.fbutton = QToolButton(self.filterBox)
        self.ficon = QPixmap('icons/filter.ico').scaled(self.sicon.width(), self.sicon.height())
        self.fbutton.setIcon(QIcon(self.ficon))
        self.fbutton.setIconSize(self.ficon.size())
        self.fbutton.setCursor(Qt.ArrowCursor)
        self.fbutton.setStyleSheet("border: none; padding: 0px;")
        self.filterBox.setStyleSheet("padding-left: %dpx;"%(int(self.ficon.width())))
        #our close icon for filters
        self.cicon = QPixmap('icons/close.ico')
        self.searchGroupLayout = QGridLayout()
        self.searchGroupLayout.addWidget(self.searchBox,0,1)
        self.searchGroupLayout.addWidget(self.filterBox,0,2)
        self.searchGroupLayout.addWidget(self.searchCols,0,3)
        self.filterBoxLayout = QHBoxLayout()
        self.searchGroupLayout.addLayout(self.filterBoxLayout,1,1,1,4,Qt.AlignLeft)
        self.searchGroup.setLayout(self.searchGroupLayout)

        dock = QDockWidget("Chromatography")
        dock.setWidget(self.chromaPanel)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.peptidePanel)
        self.tabifyDockWidget(self.peptidePanel, dock)
        base_dock = QDockWidget("Base Trace")
        base_dock.setWidget(self.baseTracePanel)
        self.addDockWidget(Qt.RightDockWidgetArea, base_dock)
        self.tabifyDockWidget(dock, base_dock)
        dock = QDockWidget("MS2")
        dock.setWidget(self.draw)
        self.draw.dock = dock
        ms1dock = QDockWidget("MS1")
        ms1dock.setWidget(self.ms1_draw)
        self.ms1_draw.dock = ms1dock
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        self.addDockWidget(Qt.RightDockWidgetArea, ms1dock)
        self.tabifyDockWidget(dock, ms1dock)
        dock = QDockWidget("Searching and Filtering")
        dock.setWidget(self.searchGroup)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        # self.splitter.addWidget(self.chromaPanel)
        # self.splitter.addWidget(self.peptidePanelDock)
        # self.splitter.addWidget(self.draw)
        # self.splitter.addWidget(self.searchGroup)#search box

        self.tree = QTreeWidget()
        dock = QDockWidget("Mass Spectra Table")
        dock.setWidget(self.tree)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        # self.splitter.addWidget(self.tree)
        self.tree.header().setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.header().customContextMenuRequested.connect(self.onHeaderRC)
        self.tree.itemDoubleClicked.connect(self.onClick)
        self.groupBy = self.LoadThread.groupBy
        self.colnames = self.LoadThread.colnames
        self.searchCols.addItems(self.colnames.keys())
        if self.chromaPanel.chromatogram is None:
            for iterObj in self.LoadThread.its:
                chroma = iterObj.getChromatogram()
                if chroma:
                    self.chromaPanel.chromatogram = chroma
                    self.chromaPanel.plotChromatogram()
                    break
                bt = iterObj.getBaseTrace()
                if bt:
                    self.baseTracePanel.chromatogram = bt
                    self.baseTracePanel.plotChromatogram()
                    break
        self.searchCols.setCurrentIndex(self.groupBy)
        self.objMap = self.LoadThread.objMap
        self.data = self.LoadThread.data
        self.tree.setColumnCount(len(self.colnames))
        self.tree.setHeaderLabels(self.colnames.keys())
        if self.data:
            self.tree.addTopLevelItems(self.data.keys())
            for i in self.data:
                i.addChildren(self.data[i].subnodes)
        self.tree.setSortingEnabled(True)

    def recurseTree(self, itemList, item):
        itemList.add(item)
        for i in xrange(item.childCount()):
            self.recurseTree(itemList, item.child(i))

    def onFilter(self):
        """
        the stored filter list contains items that were filtered out
        """
        filterTerm = self.filterBox.text()
        if not filterTerm:
            return
        col = self.searchCols.currentIndex()
        items = set(self.tree.findItems(filterTerm,Qt.MatchRecursive|Qt.MatchContains,column=col))
        allItems = set([])
        for i in xrange(self.tree.topLevelItemCount()):
            self.recurseTree(allItems, self.tree.topLevelItem(i))
        #we have all items, and we want to remove those matching our filter -- the resulting list is who to hide
        allItems = allItems-items
        if not allItems:
            return
        if self.filterList.has_key((filterTerm,col)):
            return
        self.filterList[(filterTerm,col)] = allItems
        for i in self.data:
            i.subnodes.insert(0,i)
            for subnode in i.subnodes:
                if not subnode.isHidden() and subnode in allItems:
                    subnode.setHidden(True)
            i.subnodes.pop(0)
        newFilter = QPushButton()
        newFilter.setText('%s on %s'%(filterTerm, self.searchCols.itemText(col)))
        newFilter.setIcon(QIcon(self.cicon))
        newFilter.setStyleSheet("image-position: right")
        newFilter.setLayoutDirection(Qt.RightToLeft)
        newFilter.connect(newFilter, SIGNAL("clicked()"), lambda fbutton=newFilter,fterm=filterTerm,fcol=col: self.onFilterClick(fbutton,fterm,fcol))
        spol = QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
        newFilter.setSizePolicy(spol)
        self.filterBoxes.append(newFilter)
        self.filterBoxLayout.addWidget(newFilter)

    def onFilterClick(self, fbut,fterm,fcol):
        fbut.deleteLater()
        index = self.filterList.keys().index((fterm,fcol))
        forder = []
        for i,v in enumerate(self.filterList.values()):
            if i!=index:
                forder.append(v)
        #unlist forder
        forder = set([i for sublist in forder for i in sublist])
        #forder is now a list of items who are still hidden
        for i in self.data:
            i.subnodes.insert(0,i)
            for subnode in i.subnodes:
                if subnode.isHidden():#it's an item that has been filtered out from a prior step
                    if forder:
                        if subnode not in forder:#forder are the items we are not keeping hidden
                            subnode.setHidden(False)
                    else:
                        subnode.setHidden(False)
            i.subnodes.pop(0)
        del self.filterList[(fterm,fcol)]

    def onSearch(self):
        #what column do we search, if none -- do group by column
        searchTerm = self.searchBox.text()
        if not searchTerm:
            return
        #current index we're at
        col = self.searchCols.currentIndex()
        cindex = self.tree.currentItem()
        #first we find our current index in the widgets and go from there
        items = self.tree.findItems(searchTerm, Qt.MatchRecursive|Qt.MatchContains, column=col)
        if not items:
            return
        try:
            nindex = items.index(cindex)
            if nindex<len(items)-1:
                self.tree.setCurrentItem(items[nindex+1])
            else:
                self.tree.setCurrentItem(items[0])
        except ValueError:
            self.tree.setCurrentItem(items[0])
        except IndexError:
            #doesn't exist
            self.tree.setCurrentItem(items[0])

    def Group(self, col):
        self.tree.setSortingEnabled(False)
        newData = {}
        objMap = {}
        self.groupBy = col
        for i in self.data:
            i.subnodes.insert(0,i)
            for subnode in i.subnodes:
                toAdd = subnode.data
                nid = toAdd[col]#group on peptide by default
                newNode = QTreeWidgetItem()
                [newNode.setText(ind,str(val)) for ind,val in enumerate(toAdd)]
                newNode.fileName = subnode.fileName
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
            self.onClick(self.item)
        except AttributeError:
            return

    def getFileType(self, path):
        if os.path.splitext(path)[1][1:].lower() in {'mgf', 'mzml'}:
            return 'spectra'
        else:
            return os.path.splitext(path)[1][1:]

    def loadScan(self, scan):
        if isinstance(scan, PeptideObject):
            a = figureIons.figureIons(scan, self.getTolerance())
            if scan.ms_level != 1:
                canvas = self.draw
            else:
                canvas = self.ms1_draw
            canvas.cleanup()
            canvas.setTitle(scan.title)
            canvas.dock.raise_()
            self.plotIons(a)
            if canvas.ionView['all']:
                mz = scan.scans
                x,y = zip(*[(float(i), float(j)) for i,j in mz if int(j)])
                canvas.plotXY(x,y)
            if hasattr(scan, 'ms1_scan'):
                canvas = self.ms1_draw
                canvas.cleanup()
                mz = scan.ms1_scan.scans
                x,y = zip(*[(float(i), float(j)) for i,j in mz if int(j)])
                canvas.plotXY(x,y, xRange=(min(x), max(x)))
        elif isinstance(scan, ScanObject):
            self.draw.cleanup()
            mz = scan.scans
            x,y = zip(*[(float(i), float(j)) for i,j in mz if int(j)])
            if scan.ms_level != 1:
                canvas = self.draw
            else:
                canvas = self.ms1_draw
            canvas.plotXY(x,y)
            canvas.dock.raise_()

    def plotIons(self, a, canvas=None):
        if canvas is None:
            canvas = self.draw
        ionList = a.assignPeaks()
        self.pepSequence = a.scan.peptide
        canvas.plotIons(ionList)
        canvas.peptidePanel.plotPeptide(self.pepSequence, ionList)


class PlotPanel(QWidget):
    def __init__( self, parent):
        super(PlotPanel, self).__init__()
        self.dpi = 100
        self.parent = parent
        self.toolbar = QToolBar()
        self.vbox = QVBoxLayout(self)
        self.peptidePanel = self.parent.peptidePanel
        self.chromaPanel = self.parent.chromaPanel
        self.baseTracePanel = self.parent.baseTracePanel
        self.pw = pg.PlotWidget(background=[255,255,255,255])
        self.canvas = self.pw.getViewBox()
        # self.canvas.enableAutoRange()
        # self.canvas.mouseEnabled = True
        self.vbox.addWidget(self.toolbar)
        self.vbox.addWidget(self.pw)
        self.controlDown = False
        self.startPos = False
        self.pw.sceneObj.sigMouseMoved.connect(self.onMouseMove)
        self.pw.sceneObj.sigMouseClicked.connect(self.onMouseClick)

class CustomLossWidget(QWidget):
    def __init__(self, parent):
        QWidget.__init__(self)
        self.parent = parent#parent is reference to viewertab
        self.layout = QGridLayout()
        self.setLayout(self.layout)
        self.setWindowFlags(Qt.Popup)
        mpos = QCursor().pos()
        self.move(mpos.x(),mpos.y())
        index=0
        self.enabled = []
        def makeCallback(ind,ami,l):
            return lambda: self.onToggle(ind,ami,l)
        for aa in masses.lossMasses:
            for loss in masses.lossMasses[aa]:
                self.enabled.append(QCheckBox())
                self.enabled[index].stateChanged.connect(makeCallback(index,aa,loss))
                if not self.parent.skipLosses.has_key(aa) or not loss in self.parent.skipLosses[aa]:
                    self.enabled[index].setChecked(True)
                label = QLabel('%s:%s'%(aa,loss[0]))
                display = QLabel('%s'%loss[2])
                self.layout.addWidget(self.enabled[index],index,0)
                self.layout.addWidget(label,index,1)
                self.layout.addWidget(display,index,2)
                index+=1
        self.show()

    def onToggle(self, index,aa,loss):
        obj = self.enabled[index]
        if obj.isChecked():
            if self.parent.skipLosses.has_key(aa):
                self.parent.skipLosses[aa].discard(loss)
                self.parent.reloadScan()
        else:
            try:
                self.parent.skipLosses[aa].add(loss)
            except KeyError:
                self.parent.skipLosses[aa] = set([loss])
            self.parent.reloadScan()
        # print self.parent.skipLosses

    def leaveEvent(self, event):
        self.deleteLater()
        self.destroy()

class SpectraPlot(pg.PlotDataItem):
    def __init__(self, *args, **kwargs):
        super(SpectraPlot, self).__init__(*args, **kwargs)
        self.setAcceptHoverEvents(True)
        pg.GraphicsScene.registerObject(self)

    def hoverEnterEvent(self, ev):
        print 'enter'
    def hoverLeaveEvent(self, ev):
        print 'leave'


class DrawFrame(PlotPanel):
    def __init__(self, parent, *args,**kwargs):
        self.parent = parent
        super(DrawFrame, self).__init__(parent, *args, **kwargs)
        self.cleanup()
        self.ionView = {'x': False, 'y': True, 'z': False, 'a': False, 'b': True, 'c': False, 'all': False, '++': False, '>2': False, 'pairs': True, 'signature': True}

        #for mouse drags
        self.startPos = False
        self.endPos = False
        self.controlDown = False
        self.lastClick = time.time()

        self.annotate = None
        if self.parent.fileType == 'spectra':
            self.peptidePanel.Hide()
        #toolbar stuff
        #draw bitmaps for labels
        height = self.toolbar.height()
        self.lossButton = QPushButton('Losses')
        self.lossButton.pressed.connect(self.customLosses)
        self.toolbar.addWidget(self.lossButton)
        idmap = {}
        icon_height = 35
        for ion,desc,index in zip(('x','y','z','a','b','c', '++', '>2', 'sig'),('X ions', 'Y Ions', 'Z Ions', 'A Ions', 'B Ions', 'C Ions', 'Doubly Charged Ions', 'All Charges', 'Signature ions'),xrange(9)):
            qpixmap = QPixmap(35, icon_height)
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
            qp.drawText(QPointF(0,icon_height*3/4), ion)
            qp.end()
            icon = QIcon(qpixmap)
            a = self.toolbar.addAction(icon, ion)
            a.setToolTip('Show %s'%desc)
            a.setCheckable(True)
            if ion is 'y' or ion is 'b':
                a.setChecked(True)
        qpixmap = QPixmap(35, icon_height)
        qpixmap.fill(QColor(255,255,255))
        qp = QPainter()
        qp.begin(qpixmap)
        for i,j in zip([1,5,7,9,10,17,20,27,35],[25,5,30,16,7,10,4,35,28]):
            qp.drawLine(i,icon_height,i,icon_height-j)
        qp.end()
        icon = QIcon(qpixmap)
        a = self.toolbar.addAction(icon, 'all')
        a.setCheckable(True)
        a.setChecked(True)
        a.setToolTip('Show All Spectra')

        #whether to connect or not connect dots
        qpixmap = QPixmap(35, icon_height)
        qpixmap.fill(QColor(255,255,255))
        qp = QPainter()
        qp.begin(qpixmap)
        for i,j in zip([1,5,7,9,10,17,20,27,35],[1,5,7,9,10,17,20,27,35]):
            qp.drawLine(i,icon_height,i,icon_height-j)
        qp.end()
        icon = QIcon(qpixmap)
        a = self.toolbar.addAction(icon, 'pairs')
        a.setCheckable(True)
        a.setChecked(True)
        a.setToolTip('Plot bars instead of a line plot')

        #whether to connect or not connect dots
        self.clear_quant_button = QPushButton('Clear Quant')
        self.clear_quant_button.setToolTip('This clears the chromaogram quantification selections')
        self.clear_quant_button.connect(self.clear_quant_button, SIGNAL("clicked()"),
                                        self.slot_clear_quant)
        self.toolbar.addWidget(self.clear_quant_button)



        self.toolbar.actionTriggered.connect(self.onAction)
        #error tolerance
        self.etolerance = QLineEdit("40", self.toolbar)
        self.etolerance.setToolTip("Mass Error Tolerance (ppm)")
        self.etolerance.editingFinished.connect(self.onToleranceEdit)
        self.toolbar.addWidget(self.etolerance)

    def customLosses(self):
        CustomLossWidget(self.parent)

    def onToleranceEdit(self):
        self.parent.reloadScan()

    def onAction(self, action):
        txt = str(action.iconText())
        if txt in set(['x','y','z','a','b','c', 'all', '++', '>2', 'pairs']):
            self.ionView[txt] = action.isChecked()
            self.parent.reloadScan()

    def slot_clear_quant(self):
        self.chromaPanel.plotChromatogram()
        self.baseTracePanel.plotChromatogram()

    def plotIons(self, peaks):
        for i in peaks:
            if not i:
                continue
            mz,inten,fragType, fragNum,charge,loss,aa = i[0]
            if self.ionView[fragType]:
                if not self.ionView['>2'] and charge > 2:
                    continue
                if not self.ionView['++'] and charge == 2:
                    continue
                self.points.append((mz,inten,fragType, fragNum,charge,loss,aa))
        self.draw()

    def plotXY(self, x, y, xRange=None):
        self.x = x
        self.y = y
        for i,j in zip(self.x,self.y):
            self.points.append((i,j,'spectra', None,None,None,None))
        self.draw()
        if xRange is not None:
            self.canvas.setRange(xRange=xRange)

    def cleanup(self):
        self.canvas.clear()
        self.ions = {}
        self.points = []
        self.colors = []
        self.text = []
        self.hitMapX = {}
        self.hitMapXF = {}

    def setTitle(self, title):
        self.pw.setTitle(title)
        self.pw.setLabel('left', '<font color="black">M/Z</black>')
        self.pw.setLabel('bottom', '<font color="black">Intensity</black>')

    def draw(self):
        xmax = 10
        ymax = 10
        self.hitMapX = {}
        blues = set(('x','y','z'))
        bars = []
        z_change = []
        self.points.sort(key=operator.itemgetter(2))
        for group_id,points in itertools.groupby(self.points, key=lambda x: x[2]):
            _x = []
            _y = []
            for pt_list in points:
                x,y,fragType, fragNum,charge,loss,aa = pt_list
                if fragType in blues:
                    col=[0,0,255]
                    tcolor='blue'
                elif fragType == 'spectra':
                    t='black'
                    col=[0,0,0]
                else:
                    col=[255,0,0]
                    tcolor='red'
                if (x > xmax):
                    xmax = x
                if (y > ymax):
                    ymax = y
                try:
                    self.hitMapX[x].add((y,'m/z: %d Int: %d'%(x,y)))
                except KeyError:
                    self.hitMapX[x] = set([(y,'m/z: %d Int: %d'%(x,y))])
                if self.ionView['pairs']:
                    _x += [x,x] #were connecting pairs so we do x,x
                    _y += [0, y] # and 0,y to make a vertical line
                else:
                    _x.append(x)
                    _y.append(y)
                if fragType == 'spectra':
                    z_change.append(0)
                else:
                    if not fragType:
                        continue
                    hlen=len(self.hitMapX[x])*5
                    if loss:
                        txt = '<font color="{0}">{1}{2}<sup>{3}<sup>{4}</sup></sup></font>'.format(tcolor, fragType,fragNum,loss,''.join(['+' for i in xrange(charge)]))
                    else:
                        txt = '<font color="{0}">{1}{2}<sup>{3}</sup></font>'.format(tcolor, fragType,fragNum,''.join(['+' for i in xrange(charge)]))
                    ti = pg.TextItem(html=txt)
                    # self.canvas.connect
                    self.canvas.addItem(ti)
                    plot_xy = self.pw.mapToScene(QPoint(x,y))
                    ti.setPos(plot_xy.x(), plot_xy.y())
            plot = pg.PlotDataItem()
            if group_id == 'spectra':
                _z = 0
            else:
                _z = 1
            plot.setData(_x,_y, pen=col, connect='pairs' if self.ionView['pairs'] else 'all')
            plot.setZValue(_z)
            self.canvas.addItem(plot)
        # x, y, colors = zip(*bars)
        self.canvas.getViewBox().setRange(xRange=(0, xmax), yRange=(0, ymax))
        self.pw.addLine(y=0, pen=[0,0,0])

    def mousePressEvent(self, event):
        mods = event.modifiers()
        if mods == Qt.ControlModifier:
            self.controlDown = True
        else:
            self.controlDown = False
        event.ignore()

    def onMouseClick(self, click):
        if self.controlDown:
            if self.startPos:
                #we have a start position and clicked again, so we're ending our line
                self.controlDown = False
                start_xy = self.pw.mapToView(self.startPos.scenePos())
                end_xy = self.pw.mapToView(click.scenePos())
                #we don't care about y
                x1 = start_xy.x()
                x2 = end_xy.x()
                if x1 > x2:
                    x1, x2 = x2, x1
                self.selectedRange = [x1, x2]
                self.startPos = False
                # now we update our chromatagram over the range selected
                # select just ms1 scans
                # FIXME: Make this more generic
                items = set(self.parent.tree.findItems('1',Qt.MatchRecursive|Qt.MatchContains,column=1))
                if not items:
                    items = set(self.parent.tree.findItems('*',Qt.MatchRecursive|Qt.MatchWildcard,column=1))
                new_plot = []
                for row in items:
                    scan = FILE_MAP[row.fileName].getScan(row)
                    try:
                        scans = np.array(scan.ms1_scan.scans)
                    except AttributeError:
                        scans = scan.scans
                    # get the relevant m/z
                    # c = np.array([np.array([(i,i+5) for i in range(20)]) for j in range(50)])
                    # c[np.where((c[:,:,0]>=5)&(c[:,:,0]<=15))].sum()
                    intensity = scans[np.where((scans[:,0]>x1)&(scans[:,0]<=x2)),1].sum()
                    new_plot.append((float(row.data[3]),intensity))
                # update our chromatogram
                new_plot.sort(key=operator.itemgetter(0))
                if new_plot:
                    rx,ry = zip(*new_plot)
                    self.chromaPanel.plotChromatogram(x=rx, y=ry,append=True)
            else:
                try:
                    self.canvas.removeItem(self.chromLine)
                except AttributeError:
                    pass
                self.chromLine = pg.PlotDataItem()
                self.startPos = click
                self.canvas.addItem(self.chromLine)

    def onMouseMove(self,event):
        if self.controlDown and self.startPos:
            # draw a line from start to here
            plot_xy = self.pw.mapToView(self.startPos.scenePos())
            x = [plot_xy.x()]
            y = [plot_xy.y()]
            plot_xy = self.pw.mapToView(event)
            x.append(plot_xy.x())
            y.append(plot_xy.y())
            self.chromLine.setData(x, y, pen=[1,1,1])
            return
        x = event.x()
        y = event.y()
        hList = []
        plot_xy = self.pw.mapToView(event)
        x = plot_xy.x()
        y = plot_xy.y()
        for i in self.hitMapX:
            if x-3 < i < x+3:
                yco = self.hitMapX[i]
                for j in yco:
                    if y<=j[0]:
                        hList.append((abs(i-x),i,j))
        if self.annotate:
            try:
                self.canvas.removeItem(self.annotate)
            except ValueError:
                pass
            self.annotate = None
        txt = ""
        if hList:
            hList.sort(key=operator.itemgetter(0))
            txt+='m/z = %0.3f (Int: %0.3f)\n'%(hList[0][1],hList[0][2][0])
            self.annotate = pg.TextItem(text = txt, color=[1,1,1])
            self.canvas.addItem(self.annotate)
            self.annotate.setPos(plot_xy.x(), plot_xy.y())

app = QApplication(sys.argv)
w = MainWindow()
if len(sys.argv) > 1 and sys.argv[1] == 'test':
    import os
    sample_dir = 'samples'
    files = [os.path.join(sample_dir, filename) for filename in os.listdir(sample_dir)]
    w.addPage(files)
    for i in files:
        w.addPage([i])
    import time
    while len(FILE_MAP) != len(files):
        # wait until they all load
        time.sleep(1)
    app.exit()
else:
    app.exec_()