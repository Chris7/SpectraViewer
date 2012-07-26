import re, os

class GFFObject(object):
    def __init__(self, infoList, filters, filterOnly,keydelim,exclude):
        seqid, source, gtype, start, end, score, strand, phase, info = infoList
        self.keydelim = keydelim
        self.children = []
        found=False
        if len(filters):
            for i in filters:
                if info.find(i) != -1:
                    found=True
        else:
            found=True
        for i in exclude:
            if info.find(i) != -1:
                found=False
        if not found:
            return None
        if type(score) == float:
            score = float(score)
        elif type(score) == int:
            score = int(score)
        if type(phase) == int:
            phase = int(phase)
        self.GFFInfo = {'seqid': seqid, 'source': source, 'gtype': gtype, 'start': int(start), 'end': int(end),
                        'score': score, 'strand': strand, 'phase': phase, 'attributes': {}}
#        print info
        for entry in info.split(';'):
            if not entry:
                continue
            i,v = entry.split(keydelim)
#            print filterOnly
            if not filterOnly or (filterOnly and i.strip() in filters):
#                print i,v
                self.GFFInfo['attributes'][i.strip()] = v.strip().replace('"','')
        self.id = self.getAttribute('ID')
        self.parent = self.getAttribute('Parent')
            
    def getType(self):
        return self.GFFInfo['gtype']
            
    def getId(self):
        return self.id
                
    def addChild(self, child):
        cid = child.getId()
        if cid:
            self.children[cid] = child
        
    def getParent(self):
        return self.parent
    
    def getChildren(self):
        return self.children
    
    def addParent(self, parent):
        self.parent = parent
        
    def getStart(self):
        return self.GFFInfo['start']
    
    def getEnd(self):
        return self.GFFInfo['end'] 
            
    def getSeqId(self):
        return self.GFFInfo['seqid']
    
    def getStrand(self):
        return self.GFFInfo['strand']
        
    def getAttribute(self, attribute):
        try:
            return self.GFFInfo['attributes'][attribute]
        except KeyError:
            return None
        except AttributeError:
            return None
    
    def listAttributes(self):
        return self.GFFInfo['attributes'].keys()
    
    def addAttribute(self, key, value):
        self.GFFInfo['attributes'][key]=value
    
    def writeRecord(self, handle):
        t = self.GFFInfo
        aInfo = []
        for i in t['attributes']:
            aInfo.append(str(i)+self.keydelim+str(t['attributes'][i]))
        handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(t['seqid'],t['source'],t['gtype'], str(t['start']),
                                                             str(t['end']),str(t['score']),t['strand'],str(t['phase']),
                                                             ';'.join(aInfo)))
        
    def writeFasta(self, handle, header, fasta):
        handle.write('>%s\n%s\n'%(header, fasta))
    
class GFFIterator():
    def __init__(self, filename, **kwrds):
        """
        Keywords:
        filter=[] Key Entry must be in this list to be parsed
        fitlerOnly=Boolean, True to only keep key=value pairs in filter
        keydelim= how keys are separated, defautl is =
        fast=[] List of key entries we want fast access to
        exclude=[] List of key entries that if found, we skip the record
        cols=(0,8) columns of the gff (for cases with side by side gff records made by bedTools 
        """
        self.filters=[]
        self.filterOnly=False
        self.keydelim='='
        self.fastAcc = []#keeps some keywords at the topmost level
        self.exclude = []
        self.fast = {}
        self.cols=False
        self.rand = False
        if kwrds:
            if kwrds.has_key('filter'):
                self.filters=kwrds['filter']
            if kwrds.has_key('filterOnly'):
                self.filterOnly=kwrds['filterOnly']
            if kwrds.has_key('keydelim'):
                self.keydelim=kwrds['keydelim']
            if kwrds.has_key('exclude'):
                self.exclude=kwrds['exclude']
            if kwrds.has_key('cols'):
                self.cols=kwrds['cols']
            if kwrds.has_key('random'):
                self.randAtt=kwrds['random']
                self.rand = True 
                self.ra = {}
                for i in self.randAtt:
                    self.ra[i] = {}
        if type(filename) is str or type(filename) is unicode:
            self.f = open(filename)
        else:
            self.f = filename
            
    def loadAll(self):
        try:
            while self.next():
                continue
        except StopIteration:
            return
        
    def readGFF(self, row):
        entry = row.strip().split('\t')
        if self.cols:
            entry = entry[self.cols[0]:self.cols[1]]
        ob = GFFObject(entry, self.filters,self.filterOnly,self.keydelim,self.exclude)
        return ob
        
    def getGFF(self, att, value):
        if not self.rand:
            self.rand = True
            self.ra[att] = {}
        try:
            pos = self.ra[att][value]
            self.f.seek(pos[0])
            info = self.f.read(pos[1]-pos[0])
            return self.readGFF(info)
        except KeyError:
            while not self.ra[att].has_key(value):
                ob = self.next()
            return ob
        
    def __iter__(self):
        return self
    
    def next(self):
        pos = self.f.tell()
        row = self.f.readline()
        if row == '':
            raise StopIteration
        ob = self.readGFF(row)
        if self.rand:
            for i in self.randAtt:
                self.ra[i][ob.getAttribute(i)] = (pos,self.f.tell())
        return ob

            
    def getAttribute(self, attribute, *args):
        """
        Optional keyword to search for GFF entries with a given value
        """
        value = False
        if args:
            value = args[0]
        out = {}
        if self.fast.has_key(attribute):
            if value:
                if self.fast[attribute].has_key(value):
                    for i in self.fast[attribute][value]:
                        out[i] = value
                else:
                    return None
            else:
                for i in self.fast[attribute]:
                    for ob in self.fast[attribute][i]:
                        out[ob] = i 
        else:
            for i in self.gffObjects:
                res = i.getAttribute(attribute)
                if res:
                    if not value or (value and res == value):
                        out[i] = res     
        return out
    
class GFFParser():
    def __init__(self, filename, **kwrds):
        """
        Keywords:
        filter=[] Key Entry must be in this list to be parsed
        fitlerOnly=Boolean, True to only keep key=value pairs in filter
        keydelim= how keys are separated, defautl is =
        fast=[] List of key entries we want fast access to
        exclude=[] List of key entries that if found, we skip the record
        cols=(0,8) columns of the gff (for cases with side by side gff records made by bedTools 
        """
        filters=[]
        filterOnly=False
        keydelim='='
        fastAcc = []#keeps some keywords at the topmost level
        exclude = []
        self.fast = {}
        cols=False
        if kwrds:
            if kwrds.has_key('filter'):
                filters=kwrds['filter']
            if kwrds.has_key('filterOnly'):
                filterOnly=kwrds['filterOnly']
            if kwrds.has_key('keydelim'):
                keydelim=kwrds['keydelim']
            if kwrds.has_key('fast'):
                fastAcc=kwrds['fast']
            if kwrds.has_key('exclude'):
                exclude=kwrds['exclude']
            if kwrds.has_key('cols'):
                cols=kwrds['cols']
        if type(filename) == type('string'):
            f = open(filename)
        else:
            f = filename
        self.gffObjects = set([])
        for row in f:
            entry = row.strip().split('\t')
            if cols:
                entry = entry[cols[0]:cols[1]]
#            print row
#            print entry
            ob = GFFObject(entry, filters,filterOnly,keydelim,exclude)
#            print ob.listAttributes()
            if ob:
                if fastAcc:
                    for i in fastAcc:
                        attr = ob.getAttribute(i)
                        if attr:
                            try:
                                self.fast[i][attr].append(ob)
                            except KeyError:
                                try:
                                    self.fast[i][attr] = [ob]
                                except KeyError:
                                    self.fast[i] = {attr: [ob]}
                self.gffObjects.add(ob)
            else:
                print row
            
    def getAttribute(self, attribute, *args):
        """
        Optional keyword to search for GFF entries with a given value
        """
        value = False
        if args:
            value = args[0]
        out = {}
        if self.fast.has_key(attribute):
            if value:
                if self.fast[attribute].has_key(value):
                    for i in self.fast[attribute][value]:
                        out[i] = value
                else:
                    return None
            else:
                for i in self.fast[attribute]:
                    for ob in self.fast[attribute][i]:
                        out[ob] = i 
        else:
            for i in self.gffObjects:
                res = i.getAttribute(attribute)
                if res:
                    if not value or (value and res == value):
                        out[i] = res     
        return out

class scanObject(object):
    """
    A scan object to store peaklist information in
    """
    def __init__(self):
        self.scans = []
        pass
    
    def addTitle(self, title):
        self.title = title
        
    def addCharge(self, charge):
        self.charge = charge
        
    def addMass(self, mass):
        self.mass = mass
        
    def addScan(self, scan):
        s = scan.split(' ')
        if len(s) > 1 and float(s[1]) != 0.0:
            self.scans.append(scan)
        
    def getInfo(self):
        return self.scans
    
    def getCharge(self):
        return self.charge
    
    def addRT(self, rt):
        self.rt = rt
    
    def getTitle(self):
        return self.title
    
    def getRT(self):
        try:
            return self.rt
        except:
            return None
        
    def getMZ(self):
        out = []
        for i in self.scans:
            out.append(i.split(' '))
        return out
    
    def getPrecursor(self):
        return self.mass
    
    def writeScan(self, o):
        o.write('BEGIN IONS\n')
        o.write('TITLE=%s\n'%self.title)
        try:
            o.write('RTINSECONDS=%s\n'%self.rt)
        except:
            pass
        o.write('PEPMASS=%s\n'%self.mass)
        o.write('CHARGE=%s\n'%self.charge)
        for i in self.scans:
            o.write('%s\n'%i)
        o.write('END IONS\n\n')

class mgfIterator(object):
    def __init__(self, filename, **kwrds):
        #load our index      
        tFile = list(filename)
        tFile.reverse()
        tFile = tFile[tFile.index('.')+1:]
        tFile.reverse()
        indexFile=''.join(tFile)+'.mgfi'
        self.rand = True
        self.ra = {}
        if isinstance(filename,(str,unicode)):
            self.f = open(filename, 'rb')
        elif isinstance(filename,file):
            self.f = filename
        else:
            raise Exception(TypeError,"Unknown Type of filename -- must be a file handle or a file path")
        self.tparse = re.compile(r'TITLE=(\d+),(\d+): Scan (\d+) \(rt=(.+)\)')
        self.openMGFIndex(indexFile)
        
    def openMGFIndex(self, path):
        try:
            f = open(path, 'rb')
            for row in f:
                entry = row.strip().split('\t')
                self.ra[entry[0]] = (int(entry[1]),int(entry[2]))
        except IOError:
            print 'building index for:',path
            if os.path.exists(path):
                raise Exception("Index file path: %s found, but appears incomplete"%path)
            f = open(path, 'wb')
            while True:
                try:
                    self.next()
                except StopIteration:
                    break
            for i in self.ra:
                f.write('%s\t%d\t%d\n'%(i,self.ra[i][0],self.ra[i][1]))
            
    def getScan(self, title):
        """
        allows random lookup
        """
        if self.ra.has_key(title):
            self.f.seek(self.ra[title][0],0)
            toRead = self.ra[title][1]-self.ra[title][0]
            info = self.f.read(toRead)
            scan = self.parseScan(info)
        else:
            return None
        return scan
        
    def parseScan(self, scan):
        """
        All input follows the BEGIN IONS row and ends before END IONS
        """
        setupScan = True
        foundCharge = False
        foundMass = False
        foundTitle = False
        scanObj = scanObject()
#        print 'stuff in scan',scan
        for row in scan.split('\n'):
            if not row:
                continue
            entry = row.strip().split('=')
            if len(entry) >= 2:
                if entry[0] == 'PEPMASS':
                    scanObj.addMass(entry[1])
                    foundMass = True
                elif entry[0] == 'CHARGE':
                    scanObj.addCharge(entry[1])
                    foundCharge = True
                elif entry[0] == 'TITLE':
                    title = entry[1]
                    foundTitle = True
                    scanObj.addTitle(entry[1])
                elif entry[0] == 'RTINSECONDS':
                    scanObj.addRT(entry[1])
            else:
                scanObj.addScan(row.strip())
        if foundCharge and foundMass and foundTitle:
            return scanObj
        return None
            
    def __iter__(self):
        return self
    
    def next(self):
        row = self.f.readline()
        if row == '':
            raise StopIteration
        setupScan = False
        scanInfo = ""
        while row:
            if row.find('_DISTILLER') != -1:
                raise StopIteration
            elif row.find('BEGIN IONS') != -1:
                if self.rand:
                    pStart=self.f.tell()
                setupScan=True
#                newScan=True
            elif row.find('END IONS') != -1:
                #scanObj.writeScan(open('/home/chris/test.mgf', 'wb'))
                scan = self.parseScan(scanInfo)
                if scan:
                    if self.rand:
                        self.ra[scan.getTitle()] = (pStart,pos)
                    return scan
                return None
                #newScan = False
            elif setupScan:
                scanInfo+=row
            pos = self.f.tell()
            row = self.f.readline()

class indexFolder():
    def __init__(self, folder):
        for root,dirs,files in os.walk(folder):
            for fileName in files:
                if fileName.find('.mgf') != -1 and fileName.find('.mgfi') == -1:
                    mgfIterator(os.path.join(root,fileName))

class mgfParser(object):
    def __init__(self, filename):
        if type(filename) == type('string'):
            self.f = open(filename)
        else:
            self.f = filename
        self.scans = {}
        newScan=False
        distiller=False
        setupScan = False
        dmap = {}
        dparse = re.compile(r'_DISTILLER_RAWFILE\[(\d+)\]=.+\\(.+)')
        tparse = re.compile(r'TITLE=(\d+),(\d+): Scan (\d+) \(rt=(.+)\)')
        for row in self.f:
            if row.find('_DISTILLER') != -1:
                return None
                distiller=True
                m = dparse.match(row)
                if m:
                    fname = m.group(2)
                    pos = fname.lower().find('.raw')
                    if pos != -1:
                        fname = fname[:pos]
                    dmap[m.group(1)] = fname
                else:
                    continue
            elif row.find('BEGIN IONS') != -1:
                scanObj = scanObject()
                setupScan=True
                newScan=True
            elif row.find('END IONS') != -1:
                self.scans[title] = scanObj
#                scanObj.writeScan(open('/home/chris/test.mgf', 'wb'))
                newScan = False
            elif setupScan:
                entry = row.strip().split('=')
                if len(entry) >= 2:
                    if entry[0] == 'PEPMASS':
                        scanObj.addMass(entry[1])
                    elif entry[0] == 'CHARGE':
                        scanObj.addCharge(entry[1])
                    elif entry[0] == 'TITLE':
                        if distiller:
                            m = tparse.match(row)
                            print '%s.%s.%s'%(dmap[int(m.group(1))-1],m.group(3),m.group(3))
                            print m, m.groups() 
                        else:
                            title = entry[1]
                            scanObj.addTitle(entry[1])
                else:
                    row.strip()
                    scanObj.addScan(row.strip())
                    setupScan=False
            elif newScan and not setupScan:
                scanObj.addScan(row.strip())
        print 'parsed',len(self.scans)
        
    def hasScan(self, scan):
        return self.scans.has_key(scan)
    
    def getScans(self):
        return self.scans
    
    def getScan(self, scan):
        try:
            return self.scans[scan]
        except KeyError:
            return None
        
#f = mgfIterator('/home/chris/cluster/common/src/parallel_tandem_10-12-01-1/bin/NaiveCD4/11BR4_GELFREE-PAGE_AP/Naive_CD4_GelfrEE_B41.mgf')
##    print scan.getTitle()
#f.getScan('Naive_CD4_GelfrEE_B41.548.548.2')
##    print i.getTitle()