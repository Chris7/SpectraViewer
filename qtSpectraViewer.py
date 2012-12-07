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

import fileIterators, operator, os, re, xml.etree.cElementTree as etree, random, matplotlib as mpl, figureIons, sys, masses
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import SpecView, SettingsPanel

try:
    from collections import OrderedDict
except ImportError:
    #compatibility
    ## {{{ http://code.activestate.com/recipes/576693/ (r9)
    # Backport of OrderedDict() class that runs on Python 2.4, 2.5, 2.6, 2.7 and pypy.
    # Passes Python2.7's test suite and incorporates all the latest updates.
    
    try:
        from thread import get_ident as _get_ident
    except ImportError:
        from dummy_thread import get_ident as _get_ident
    
    try:
        from _abcoll import KeysView, ValuesView, ItemsView
    except ImportError:
        pass
    
    
    class OrderedDict(dict):
        'Dictionary that remembers insertion order'
        # An inherited dict maps keys to values.
        # The inherited dict provides __getitem__, __len__, __contains__, and get.
        # The remaining methods are order-aware.
        # Big-O running times for all methods are the same as for regular dictionaries.
    
        # The internal self.__map dictionary maps keys to links in a doubly linked list.
        # The circular doubly linked list starts and ends with a sentinel element.
        # The sentinel element never gets deleted (this simplifies the algorithm).
        # Each link is stored as a list of length three:  [PREV, NEXT, KEY].
    
        def __init__(self, *args, **kwds):
            '''Initialize an ordered dictionary.  Signature is the same as for
            regular dictionaries, but keyword arguments are not recommended
            because their insertion order is arbitrary.
    
            '''
            if len(args) > 1:
                raise TypeError('expected at most 1 arguments, got %d' % len(args))
            try:
                self.__root
            except AttributeError:
                self.__root = root = []                     # sentinel node
                root[:] = [root, root, None]
                self.__map = {}
            self.__update(*args, **kwds)
    
        def __setitem__(self, key, value, dict_setitem=dict.__setitem__):
            'od.__setitem__(i, y) <==> od[i]=y'
            # Setting a new item creates a new link which goes at the end of the linked
            # list, and the inherited dictionary is updated with the new key/value pair.
            if key not in self:
                root = self.__root
                last = root[0]
                last[1] = root[0] = self.__map[key] = [last, root, key]
            dict_setitem(self, key, value)
    
        def __delitem__(self, key, dict_delitem=dict.__delitem__):
            'od.__delitem__(y) <==> del od[y]'
            # Deleting an existing item uses self.__map to find the link which is
            # then removed by updating the links in the predecessor and successor nodes.
            dict_delitem(self, key)
            link_prev, link_next, key = self.__map.pop(key)
            link_prev[1] = link_next
            link_next[0] = link_prev
    
        def __iter__(self):
            'od.__iter__() <==> iter(od)'
            root = self.__root
            curr = root[1]
            while curr is not root:
                yield curr[2]
                curr = curr[1]
    
        def __reversed__(self):
            'od.__reversed__() <==> reversed(od)'
            root = self.__root
            curr = root[0]
            while curr is not root:
                yield curr[2]
                curr = curr[0]
    
        def clear(self):
            'od.clear() -> None.  Remove all items from od.'
            try:
                for node in self.__map.itervalues():
                    del node[:]
                root = self.__root
                root[:] = [root, root, None]
                self.__map.clear()
            except AttributeError:
                pass
            dict.clear(self)
    
        def popitem(self, last=True):
            '''od.popitem() -> (k, v), return and remove a (key, value) pair.
            Pairs are returned in LIFO order if last is true or FIFO order if false.
    
            '''
            if not self:
                raise KeyError('dictionary is empty')
            root = self.__root
            if last:
                link = root[0]
                link_prev = link[0]
                link_prev[1] = root
                root[0] = link_prev
            else:
                link = root[1]
                link_next = link[1]
                root[1] = link_next
                link_next[0] = root
            key = link[2]
            del self.__map[key]
            value = dict.pop(self, key)
            return key, value
    
        # -- the following methods do not depend on the internal structure --
    
        def keys(self):
            'od.keys() -> list of keys in od'
            return list(self)
    
        def values(self):
            'od.values() -> list of values in od'
            return [self[key] for key in self]
    
        def items(self):
            'od.items() -> list of (key, value) pairs in od'
            return [(key, self[key]) for key in self]
    
        def iterkeys(self):
            'od.iterkeys() -> an iterator over the keys in od'
            return iter(self)
    
        def itervalues(self):
            'od.itervalues -> an iterator over the values in od'
            for k in self:
                yield self[k]
    
        def iteritems(self):
            'od.iteritems -> an iterator over the (key, value) items in od'
            for k in self:
                yield (k, self[k])
    
        def update(*args, **kwds):
            '''od.update(E, **F) -> None.  Update od from dict/iterable E and F.
    
            If E is a dict instance, does:           for k in E: od[k] = E[k]
            If E has a .keys() method, does:         for k in E.keys(): od[k] = E[k]
            Or if E is an iterable of items, does:   for k, v in E: od[k] = v
            In either case, this is followed by:     for k, v in F.items(): od[k] = v
    
            '''
            if len(args) > 2:
                raise TypeError('update() takes at most 2 positional '
                                'arguments (%d given)' % (len(args),))
            elif not args:
                raise TypeError('update() takes at least 1 argument (0 given)')
            self = args[0]
            # Make progressively weaker assumptions about "other"
            other = ()
            if len(args) == 2:
                other = args[1]
            if isinstance(other, dict):
                for key in other:
                    self[key] = other[key]
            elif hasattr(other, 'keys'):
                for key in other.keys():
                    self[key] = other[key]
            else:
                for key, value in other:
                    self[key] = value
            for key, value in kwds.items():
                self[key] = value
    
        __update = update  # let subclasses override update without breaking __init__
    
        __marker = object()
    
        def pop(self, key, default=__marker):
            '''od.pop(k[,d]) -> v, remove specified key and return the corresponding value.
            If key is not found, d is returned if given, otherwise KeyError is raised.
    
            '''
            if key in self:
                result = self[key]
                del self[key]
                return result
            if default is self.__marker:
                raise KeyError(key)
            return default
    
        def setdefault(self, key, default=None):
            'od.setdefault(k[,d]) -> od.get(k,d), also set od[k]=d if k not in od'
            if key in self:
                return self[key]
            self[key] = default
            return default
    
        def __repr__(self, _repr_running={}):
            'od.__repr__() <==> repr(od)'
            call_key = id(self), _get_ident()
            if call_key in _repr_running:
                return '...'
            _repr_running[call_key] = 1
            try:
                if not self:
                    return '%s()' % (self.__class__.__name__,)
                return '%s(%r)' % (self.__class__.__name__, self.items())
            finally:
                del _repr_running[call_key]
    
        def __reduce__(self):
            'Return state information for pickling'
            items = [[k, self[k]] for k in self]
            inst_dict = vars(self).copy()
            for k in vars(OrderedDict()):
                inst_dict.pop(k, None)
            if inst_dict:
                return (self.__class__, (items,), inst_dict)
            return self.__class__, (items,)
    
        def copy(self):
            'od.copy() -> a shallow copy of od'
            return self.__class__(self)
    
        @classmethod
        def fromkeys(cls, iterable, value=None):
            '''OD.fromkeys(S[, v]) -> New ordered dictionary with keys from S
            and values equal to v (which defaults to None).
    
            '''
            d = cls()
            for key in iterable:
                d[key] = value
            return d
    
        def __eq__(self, other):
            '''od.__eq__(y) <==> od==y.  Comparison to another OD is order-sensitive
            while comparison to a regular mapping is order-insensitive.
    
            '''
            if isinstance(other, OrderedDict):
                return len(self)==len(other) and self.items() == other.items()
            return dict.__eq__(self, other)
    
        def __ne__(self, other):
            return not self == other
    
        # -- the following methods are only used in Python 2.7 --
    
        def viewkeys(self):
            "od.viewkeys() -> a set-like object providing a view on od's keys"
            return KeysView(self)
    
        def viewvalues(self):
            "od.viewvalues() -> an object providing a view on od's values"
            return ValuesView(self)
    
        def viewitems(self):
            "od.viewitems() -> a set-like object providing a view on od's items"
            return ItemsView(self)
    ## end of http://code.activestate.com/recipes/576693/ }}}


#some dangerous globals for us
searchPaths = set([])
fileNames = {}
titleParse = re.compile('(.+?)\.\d+\.\d+\.\d+|\w+')

#filename mapping -- we keep it global to save memory
fileMapping = {}

def loadFolder(path):
    global searchPaths,fileNames
    for root,dir,files in os.walk(path):
        for fileName in files:
            if '.mgf' in fileName and '.mgfi' not in fileName:
                fileNames[fileName[:fileName.find('.mgf')]] = os.path.join(root,fileName)
    searchPaths.add(path)

def loadConfig():
    global searchPaths
    try:
        config = etree.ElementTree(file='spectraConf.xml')
    except IOError:#doesn't exist
        return
    root = config.getroot()
    searchPaths = root.find('GFFFilePaths').text
    if searchPaths:
        searchPaths = searchPaths.split(',')
    losses = root.find('LossMasses')
    for child in losses:
        loss = []
        for mass, frags, display in zip(child.findall('Mass'),child.findall('Fragments'), child.findall('Display')):
            if mass.text and frags.text and display.text:
                loss.append((float(mass.text),tuple(frags.text.split(',')),display.text))
        loss = tuple(loss)
        masses.lossMasses[child.text] = loss
    
    
def saveConfig():
    global searchPaths
    root = etree.Element('SpectraViewerConfig')
    sub = etree.SubElement(root,'GFFFilePaths')
    if searchPaths:
        sub.text = ','.join(searchPaths)
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
        for fileindex,path in enumerate(self.files):
            if fileMapping.values().count(path):
                for i in fileMapping:
                    if fileMapping[i] == path:
                        fileIndex = i
                        break
            else:
                fileIndex = len(fileMapping)
                fileMapping[fileIndex] = path
            iterObj = FileObject(path)
            self.its.append(iterObj)
            pepSpecType = ('gff', 'xml', 'msf', 'dat')
            specType = ('mgf',)
            if [i for i in pepSpecType if i in path]:
                iterObj.fileType = [i for i in pepSpecType if i in path][0]#this changes based on the file input somewhat
                self.colnames = OrderedDict([("Scan Title",(str,'id',1)), ("Peptide",(str,'peptide',1)), ("Modifications",(str,'getModifications',0)), ("Charge",(int,'charge',1)), ("Accession",(str, 'acc',1))])
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
                            toAdd = [getattr(i,j[1]) if j[2] else getattr(i,j[1])() for j in self.colnames.values()]
                            nid = toAdd[self.groupBy]#group on peptide by default
                            node = self.objMap.get(nid)
                            if node is None:
                                newNode = QTreeWidgetItem()
                                [newNode.setText(i,str(v)) for i,v in enumerate(toAdd)]
                                newNode.fileName = fileIndex
                                newNode.subnodes = []
                                newNode.data = toAdd
                                self.data[newNode] = newNode
                                self.objMap[nid] = newNode
                            else:
                                newNode = QTreeWidgetItem(node)
                                [newNode.setText(i,str(v)) for i,v in enumerate(toAdd)]
                                newNode.fileName = fileIndex
                                newNode.data = toAdd
                                self.data[node].subnodes.append(newNode)
                elif iterObj.fileType == 'xml' or iterObj.fileType == 'msf' or iterObj.fileType == 'dat':
                    if iterObj.fileType == 'xml':
                        iterObj.iObj = fileIterators.XTandemXML(path)
                        self.colnames['Expect'] = (float,'expect',1)
                    elif iterObj.fileType == 'dat':
                        iterObj.iObj = fileIterators.MascotDATIterator(path)
                        self.colnames['Hit Id'] = (int,'hit',1)
                        self.colnames['Rank'] = (int,'rank',1)
                    elif iterObj.fileType == 'msf':
                        iterObj.iObj = fileIterators.ThermoMSFIterator(path)
                        self.colnames['Spectrum ID'] = (int, 'hit',1)
                        self.colnames['Confidence Level'] = (int, 'confidence',1)
                        self.colnames['Search Engine Rank'] = (int, 'rank',1)
                    for i in iterObj.iObj:
                        self.emit(SIGNAL('updateProgress'),iterObj.iObj.getProgress())
                        toAdd = [getattr(i,j[1]) if j[2] else getattr(i,j[1])() for j in self.colnames.values()]
                        nid = toAdd[self.groupBy]#group on peptide by default
                        node = self.objMap.get(nid)
                        if node is None:
                            newNode = QTreeWidgetItem()
                            [newNode.setText(i,str(v)) for i,v in enumerate(toAdd)]
                            newNode.fileName = fileIndex
                            newNode.data = toAdd
                            newNode.subnodes = []
                            self.data[newNode] = newNode
                            self.objMap[nid] = newNode
                        else:
                            newNode = QTreeWidgetItem(node)
                            [newNode.setText(i,str(v)) for i,v in enumerate(toAdd)]
                            newNode.fileName = fileIndex
                            newNode.data = toAdd
                            self.data[node].subnodes.append(newNode)
            elif [i for i in specType if i in path]:
                iterObj.fileType = 'spectra'#these are all generic more or less, so spectra works
                iterObj.dataType = 'spectra'
                self.colnames = OrderedDict([("Scan Title",(str,'title',1)), ("Charge",(int,'charge',1)), ("RT",(str,'rt',1)), ("Precursor Mass",(str,'mass',1))])
                iterObj.iObj = fileIterators.mgfIterator(path, random=True)
                self.groupBy = 0
                for i in iterObj.iObj:
                    if not i:
                        continue
                    self.emit(SIGNAL('updateProgress'),iterObj.iObj.getProgress())
                    toAdd = [getattr(i,j[1]) if j[2] else getattr(i,j[1])() for j in self.colnames.values()]
                    nid = toAdd[self.groupBy]#group on title by default (should be unique)
                    newNode = QTreeWidgetItem()
                    [newNode.setText(i,str(v)) for i,v in enumerate(toAdd)]
                    newNode.fileName = fileIndex
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
            self.emit(SIGNAL('fileDone'),fileindex,len(self.files))

class PeptidePanel(QWidget):
    def __init__(self, parent):
        QWidget.__init__(self)
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
        self.validExtensions = set(['.xml', '.msf', '.mgf', '.dat'])
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
        print masses.lossMasses[aa]  
        
        
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
        else:
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
        
class ViewerTab(QWidget):
    def __init__(self, parent, files):
        QWidget.__init__(self)
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
        
    def onClick(self, item, col):
        self.item = item
        scanTitle = item.text(0)
        if self.getFileType(fileMapping[item.fileName]) == 'msf':
            self.loadScan(fileMapping[item.fileName], str(scanTitle),specId=item.text(5),peptide=item.text(1))
        elif self.getFileType(fileMapping[item.fileName]) == 'dat':
            self.loadScan(fileMapping[item.fileName], str(scanTitle),hitId=item.text(5),rank=item.text(6))
        else:
            self.loadScan(fileMapping[item.fileName], str(scanTitle))
        
    def onHeaderRC(self, pos):
        gpos = self.mapToGlobal(pos)
        menu = QMenu()
        menu.addAction("Group By Column")
        ty = self.tree.pos().y()
        gpos.setY(gpos.y()+ty)
        selected = menu.exec_(self.mapToParent(gpos))
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
        if hasattr(self, 'splitter'):
            iterObjects = self.LoadThread.its
            self.data = self.LoadThread.data
            self.objMap = self.LoadThread.objMap
            for it in iterObjects:
                if it.fileType == 'gff':
                    self.parent.gffFiles[it.path] = it.iObj
                elif it.fileType == 'xml' or it.fileType == 'msf' or it.fileType == 'dat':
                    self.parent.pepFiles[it.path] = it.iObj
                elif it.fileType == 'spectra':
                    self.parent.mgfFiles[it.path] = it.iObj
            if self.data:
                for i in self.data.keys():
                    if not i.treeWidget():
                        self.tree.addTopLevelItem(i)
                        i.addChildren(self.data[i].subnodes)
            self.tree.setSortingEnabled(True)
            self.onFileText.hide()
            self.progress.hide()
            self.progressText.hide()
            return
        self.splitter = QSplitter()
        self.splitter.setChildrenCollapsible(True)
        self.vl.removeWidget(self.progress)
        self.vl.removeWidget(self.progressText)
        self.vl.addWidget(self.splitter)
        self.progress.hide()
        self.progressText.hide()
        self.onFileText.hide()
        self.vl.addWidget(self.onFileText)
        self.vl.addWidget(self.progressText)
        self.vl.addWidget(self.progress)
        self.splitter.setOrientation(Qt.Vertical)
        self.splitter.setParent(self)
        self.peptidePanel = PeptidePanel(self)
        self.draw = DrawFrame(self)
        self.searchBox = QLineEdit()
        self.searchBox.editingFinished.connect(self.onSearch)
        #search button
        self.sbutton = QToolButton(self.searchBox)
        self.sicon = QPixmap('icons/search.ico')
        self.sbutton.setIcon(QIcon(self.sicon))
        self.sbutton.setIconSize(self.sicon.size())
        self.sbutton.setCursor(Qt.ArrowCursor)
        self.sbutton.setStyleSheet("border: none; padding: 0px; ");
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
        self.fbutton.setStyleSheet("border: none; padding: 0px; ");
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
        
        self.splitter.addWidget(self.peptidePanel)
        self.splitter.addWidget(self.draw)
        self.splitter.addWidget(self.searchGroup)#search box
        
        self.tree = QTreeWidget()
        self.splitter.addWidget(self.tree)
        self.tree.header().setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.header().customContextMenuRequested.connect(self.onHeaderRC)
        self.tree.itemDoubleClicked.connect(self.onClick)
        self.groupBy = self.LoadThread.groupBy
        iterObjects = self.LoadThread.its
        self.colnames = self.LoadThread.colnames
        self.searchCols.addItems(self.colnames.keys())
        self.searchCols.setCurrentIndex(self.groupBy)
        self.objMap = self.LoadThread.objMap
        self.data = self.LoadThread.data
        for it in iterObjects:
            if it.fileType == 'gff':
                self.parent.gffFiles[it.path] = it.iObj
            elif it.fileType == 'xml' or it.fileType == 'msf' or it.fileType == 'dat':
                self.parent.pepFiles[it.path] = it.iObj
            elif it.fileType == 'spectra':
                self.parent.mgfFiles[it.path] = it.iObj
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
                [newNode.setText(i,str(v)) for i,v in enumerate(toAdd)]
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
            self.onClick(self.item, 0)
        except AttributeError:
            return
        
    def getFileType(self, path):
        if os.path.splitext(path)[1][1:] == 'mgf':
            return 'spectra'
        else:
            return os.path.splitext(path)[1][1:]
        
    def loadScan(self, fileO, title, **kwrds):
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
            self.precursor = scan.mass
            mz = scan.scans
            x = []
            y = []
            for i in mz:
                x.append(float(i[0]))
                y.append(float(i[1]))
            self.pepSequence = gob.getAttribute('Sequence')
            a = figureIons.figureIons(self.pepSequence,gob.getAttribute('Charge'),mods, self.getTolerance(),skipLoss=self.skipLosses)
            self.draw.cleanup()
            self.draw.setTitle(title)
            self.plotIons(x,y,a)
            if self.draw.ionView['all']:
                self.draw.plotXY(x,y)
        elif self.getFileType(path) == 'xml' or self.getFileType(path) == 'msf' or self.getFileType(path) == 'dat':
            if self.getFileType(path) == 'xml':
                scan = self.parent.pepFiles[path].getScan(title)
            elif self.getFileType(path) == 'msf':
                scan = self.parent.pepFiles[path].getScan(title,kwrds['specId'],kwrds['peptide'])
            elif self.getFileType(path) == 'dat':
                scan = self.parent.pepFiles[path].getScan(title,kwrds['hitId'], kwrds['rank'])
            if not scan:
                return
            mods = scan.getModifications()
            self.precursor = scan.mass
            mz = scan.scans
            x = []
            y = []
            for i in mz:
                x.append(float(i[0]))
                y.append(float(i[1]))
            self.pepSequence = scan.peptide
            a = figureIons.figureIons(self.pepSequence,scan.charge,mods, self.getTolerance())
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
            mz = scan.scans
            x = []
            y = []
            for i in mz:
                x.append(float(i[0]))
                y.append(float(i[1]))
            self.draw.cleanup()
            self.draw.plotXY(x,y)
        
    def plotIons(self, x,y,a):
        ionList = a.assignPeaks(x,y)
        self.draw.plotIons([i[0] for i in ionList])
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
        self.toolbar.coordinates = None
        self.vbox = QVBoxLayout(self)
        self.peptidePanel = self.parent.peptidePanel
        self.vbox.addWidget(self.toolbar)
        self.vbox.addWidget(self.canvas)
        self.canvas.setGeometry(QRect(10, 150, 490, 390))
        self.axes = self.figure.add_subplot(111)
        self.mouse = 0

    def onMouseDown(self, event):
        self.mouse = 1
        
    def onMouseRelease(self, event):
        self.mouse = 0

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
        print self.parent.skipLosses
            
    def leaveEvent(self, event):
        self.deleteLater()
        self.destroy()

class DrawFrame(PlotPanel):
    def __init__(self, parent, *args,**kwargs):
        self.parent = parent
        PlotPanel.__init__(self, parent, **kwargs)        
        self.cleanup()
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
        self.lossButton = QPushButton('Losses')
        self.lossButton.pressed.connect(self.customLosses)
        self.toolbar.addWidget(self.lossButton)
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
        
    def customLosses(self):
        CustomLossWidget(self.parent)
        
    def onToleranceEdit(self):
        self.parent.reloadScan()
        
    def onAction(self, action):
        txt = str(action.iconText())
        if txt in set(['x','y','z','a','b','c', 'all', '++', '>2']):
            self.ionView[txt] = action.isChecked()
            self.parent.reloadScan()
        
    def plotIons(self, peaks):
        for i in peaks:
            if not i:
                continue
            mz,inten,fragType, fragNum,charge,loss,aa = i
            if self.ionView[fragType]:
                if not self.ionView['>2'] and charge > 2:
                    continue
                if not self.ionView['++'] and charge == 2:
                    continue
                self.points.append((mz,inten,fragType, fragNum,charge,loss,aa))
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
#w.addPage(['sampleMgf.mgf'])
app.exec_()
