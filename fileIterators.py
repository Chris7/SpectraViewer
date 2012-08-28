import re, os, masses, sys
try:
    from lxml import etree
except ImportError:
    print 'lxml is required to parse X!tandem xml files due to the namespaces employed'

#regex for common use
scanSplitter = re.compile(r'[\t\s]')
distillerParse = re.compile(r'_DISTILLER_RAWFILE\[(\d+)\]=\(1\)(.+)')

class scanObject(object):
    """
    A scan object to store peaklist information in.  There are some getter/setter functions to keep variables similiar across file types
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
        
    def addScan(self, mz, intensity):
        self.scans.append((float(mz),float(intensity)))
        
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
            return ''
        
    def getMZ(self):
        return self.scans
    
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
        
class peptideObject(scanObject):
    """
    An enhanced scan object that can store peptide information as well
    """
    def __init__(self):
        scanObject.__init__(self)
        self.mods = set([])
        self.peptide = ""
        
    def setPeptide(self, peptide):
        self.peptide = peptide
        
    def addModification(self, aa,position, modMass):
        """
        Modifications are represented in the mascot manner:
        AA#Type
        such as: M35(o) for an oxidized methionine at residue 35
        """
        #clean up xtandem
        try:
            float(modMass)
            if modMass[0] == '-':
                modMass=modMass[1:]
            modType = ""
            for i in masses.mod_weights:
                if modMass in str(masses.mod_weights[i]):
                    #found it
                    modType = i
            if not modType:
                print 'mod not found',modMass
        except ValueError:
            #string passed, nice
            modType = modMass
        self.mods.add('%s%s(%s)'%(aa,position,modType))
        
    def setExpect(self, expect):
        self.expect = expect
    
    def setId(self, id):
        self.id = id
        
    def setAccession(self, acc):
        self.acc = acc
        
    def getId(self):
        return self.id
    
    def getAccession(self):
        return self.acc
    
    def getModifications(self):
        return '|'.join(self.mods)
    
    def getPeptide(self):
        return self.peptide
    

class XTandemXML(object):
    """
    Parser for X!Tandem XML Files.
    """
    def __init__(self, filename, **kwrds):
        #parse in our X!Tandem xml file
        if kwrds and kwrds['exclude']:
            exclude = set(kwrds['exclude'])
        else:
            exclude = set()
        try:
            dom1 = etree.parse(filename)
            self.lxml = True
        except NameError:
            self.lxml = False 
            print 'XTandem parsing unavailable: lxml is required to parse X!tandem xml files due to the namespaces employed'
            return 
        if self.lxml:
            self.group = dom1.findall("group")
            self.groupMap = {}
            self.index = 0
        else: 
            self.nest = 0
        self.db = None
        
    def __iter__(self):
        return self
    
    def parselxml(self, group):
        subnote = list(group.iter("note"))
        for i in subnote:
            if (i.attrib["label"] == "Description"):
                experiment = i.text.strip()
        try:
            expect = group.attrib["expect"]
        except KeyError:
            self.index+=1
            self.next()
        charge = group.attrib["z"]
        premass = group.attrib["mh"]
        rt = group.attrib["rt"]
        proteins = list(group.iter("protein"))
        fullProtein = ""
        for protein in proteins:
            scanObj = peptideObject()
            scanObj.addCharge(charge)
            scanObj.addMass(premass*int(charge))
            scanObj.addRT(rt)
            sgroup = group.iter("group")
            for i in sgroup:
                #This is horridly inefficient...
                if "fragment ion mass spectrum" in i.attrib["label"]:
                    ab = i.iter('{http://www.bioml.com/gaml/}Xdata')
                    for j in ab:
                        mzIter = j.iter('{http://www.bioml.com/gaml/}values')
                        for k in mzIter:
#                                print k.text
                            mz = [mval for mval in k.text.strip().replace('\n',' ').split(' ')]
                    ab = i.iter('{http://www.bioml.com/gaml/}Ydata')
                    for j in ab:
                        mzIter = j.iter('{http://www.bioml.com/gaml/}values')
                        for k in mzIter:
#                                print k.text
                            inten = [mval for mval in k.text.strip().replace('\n',' ').split(' ')]
                    for j,k in zip(mz,inten):
                        scanObj.addScan(j,k)
            domain = list(protein.iter("domain"))[0]#we only have one domain per protein instance
            note = list(protein.iter("note"))[0]#same here
            mods = list(protein.iter("aa"))#we can have multiple modifications
            if not self.db:
                files = list(protein.iter("file"))
                self.db = files[0].attrib["URL"]
            id = domain.attrib["id"]
            start = domain.attrib["start"]
            end = domain.attrib["end"]
            peptide = domain.attrib["seq"]
            pExpect = domain.attrib["expect"]
            for mod in mods:
                scanObj.addModification(mod.attrib["type"],mod.attrib["at"],mod.attrib["modified"])
            scanObj.setPeptide(peptide)
            scanObj.setExpect(pExpect)
            scanObj.setId(id)
            self.groupMap[id] = self.index
            scanObj.addTitle(id)
            scanObj.setAccession(note.text)
        self.index+=1
        return scanObj
    
    def next(self):
        try:
            group = self.group[self.index]
        except IndexError:
            raise StopIteration
        if self.lxml:
            return self.parselxml(group)
        else:
            return self.parsetext()
                
    def getScan(self, id):
        index = self.groupMap[id]
        return self.parselxml(self.group[index])
        
class GFFObject(object):
    def __init__(self, infoList, filters, filterOnly,keydelim,exclude):
        seqid, source, gtype, start, end, score, strand, phase, info = infoList
        self.keydelim = keydelim
        self.children = []
        found=False
        if len(filters):
            for i in filters:
                if i in info:
                    found=True
        else:
            found=True
        for i in exclude:
            if i in info:
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
        for entry in info.split(';'):
            if not entry:
                continue
            i,v = entry.split(keydelim)
            if not filterOnly or (filterOnly and i.strip() in filters):
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
            ob = GFFObject(entry, filters,filterOnly,keydelim,exclude)
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


class mgfIterator(object):
    def __init__(self, filename, **kwrds):
        #load our index      
        tFile = list(filename)
        tFile.reverse()
        tFile = tFile[tFile.index('.')+1:]
        tFile.reverse()
        indexFile=''.join(tFile)+'.mgfi'
        self.scanSplit = re.compile(r'[\s\t]')
        self.rand = True
        self.titleMap = {}
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
            self.f.seek(0)
            
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
#                    if self.titleMap:
#                        pos = entry[1].find(',')
#                        title = self.titleMap[int(entry[1][:entry[1].find(',')])]
#                    else:
                    title = '='.join(entry[1:])
                    foundTitle = True
                    scanObj.addTitle(title)
                elif entry[0] == 'RTINSECONDS':
                    scanObj.addRT(entry[1])
            else:
                mz,intensity = self.scanSplit.split(row.strip())
                scanObj.addScan(mz,intensity)
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
            if '_DISTILLER' in row:
                if row:
                    m = distillerParse.match(row)
                    if m:
                        self.titleMap[int(m.group(1))+1] = m.group(2)
                pass
            elif 'BEGIN IONS' in row:
                if self.rand:
                    pStart=self.f.tell()
                setupScan=True
#                newScan=True
            elif 'END IONS' in row:
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
                if '.mgf' in fileName and '.mgfi' not in fileName:
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
            if '_DISTILLER' in row:
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
            elif 'BEGIN IONS' in row:
                scanObj = scanObject()
                setupScan=True
                newScan=True
            elif 'END IONS' in row:
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