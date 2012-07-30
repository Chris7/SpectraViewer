import re, masses, operator

modParse = re.compile('([A-Z])(\d+)\((.+)\)')
class figureIons(object):
    def __init__(self,seq,ioncharge,mods, tolerance):
        self.sequence=seq.upper()
        self.ioncharge = ioncharge
        self.modList = {}
        if mods:
            mods = mods.split('|')
            for mod in mods:
                modification, start, modType = modParse.search(mod).groups()
                self.modList[int(start)] = modType.lower()
        self.tolerance = tolerance
        
    def predictPeaks(self):
        mass = 0
        self.peakTable = []
        modList = self.modList
#        ioncharge = self.ioncharge
        hw = masses.mod_weights['h'][0]
        losses = [(0,('A','B','C','X','Y','Z'),'')]
        lossMasses = masses.lossMasses
        charge=1#for starting nh3
        sLen=len(self.sequence)
        for i,v in enumerate(self.sequence):
            ionIndex = i+1
            try:
                if self.sequence[ionIndex] in lossMasses:
                    for k in lossMasses[self.sequence[ionIndex]]:
                        losses.append(k)
            except IndexError:
                pass
            mass+=masses.protein_weights[v][0]
            charge+=masses.protein_weights[v][1]
            try:
                mass+=masses.mod_weights[modList[i+1]][0]
                charge+=masses.mod_weights[modList[i+1]][1]
            except KeyError:
                pass
            for icharge in xrange(1,charge+1):
                for lossType in losses:
                    lossMass = lossType[0]
#                    print 'abc loss type',v,lossType,lossMass,i
                    tcharge = icharge
                    if 0<i<sLen-1 and not lossType[2]:
                        a = (mass-masses.mod_weights['cho'][0]+hw+hw*tcharge)/tcharge
                        self.peakTable.append((a,'A%d'%tcharge,'A%d+%d'%(ionIndex,tcharge),v))
                    elif 0<i<sLen-1 and 'A' in lossType[1]:
                        a = (mass-masses.mod_weights['cho'][0]+hw+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((a,'A%d'%tcharge,'A%d%s+%d'%(ionIndex,lossType[2],tcharge),v))
                    if 0<i<sLen-1 and not lossType[2]:
                        b = (mass+masses.mod_weights['h'][0]-hw+hw*tcharge)/tcharge
                        self.peakTable.append((b,'B%d'%tcharge,'B%d+%d'%(ionIndex,tcharge),v))
                    elif 0<i<sLen-1 and 'B' in lossType[1]:
                        b = (mass+masses.mod_weights['h'][0]-hw+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((b,'B%d'%tcharge,'B%d%s+%d'%(ionIndex,lossType[2],tcharge),v))
                    if i<sLen-1 and not lossType[2]:
                        c = (mass+masses.mod_weights['nh2'][0]+hw+hw*tcharge)/tcharge
                        self.peakTable.append((c,'C%d'%tcharge,'C%d+%d'%(ionIndex,tcharge),v))
                    elif i<sLen-1 and 'C' in lossType[1]:
                        c = (mass+masses.mod_weights['nh2'][0]+hw+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((c,'C%d'%tcharge,'C%d%s+%d'%(ionIndex,lossType[2],tcharge),v))
        charge=int(self.ioncharge)-1#for rightmost charge -- we remove the NH3 charge
        losses = [(0,('A','B','C','X','Y','Z'),'')]
        mass = 0
        for i,v in enumerate(reversed(self.sequence)):
            mass+=masses.protein_weights[v][0]
            charge+=masses.protein_weights[v][1]
            try:
                mass+=masses.mod_weights[modList[i+1]][0]
                charge+=masses.mod_weights[modList[i+1]][1]
            except KeyError:
                pass
            ionIndex = i+1
            #we need to look one ahead
            try:
                if v in lossMasses:
                    for k in lossMasses[v]:
                        losses.append(k)
            except IndexError:
                pass
            for icharge in xrange(1,charge+1):
                for lossType in losses:
                    lossMass = lossType[0]
#                    print 'xyz loss type',v,lossType,lossMass,lossType[1]
#                    tcharge = charge+icharge
                    tcharge=icharge
                    #tcharge = icharge+1
                    if not lossType[2]:
                        x = (mass+masses.mod_weights['co'][0]+masses.mod_weights['oh'][0]-hw+hw*tcharge)/tcharge
                        self.peakTable.append((x,'X%d'%tcharge,'X%d+%d'%(ionIndex,tcharge),v))
                    elif 'X' in lossType[1]:
                        x = (mass+masses.mod_weights['co'][0]+masses.mod_weights['oh'][0]-hw+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((x,'X%d'%tcharge,'X%d%s+%d'%(ionIndex,lossType[2],tcharge),v))
                    if not lossType[2]:
                        y = (mass+masses.mod_weights['h'][0]+masses.mod_weights['oh'][0]+hw*tcharge)/tcharge
                        self.peakTable.append((y,'Y%d'%tcharge,'Y%d+%d'%(ionIndex,tcharge),v))
                    elif 'Y' in lossType[1]:
                        y = (mass+masses.mod_weights['h'][0]+masses.mod_weights['oh'][0]+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((y,'Y%d'%tcharge,'Y%d%s+%d'%(ionIndex,lossType[2],tcharge),v))
                    if not lossType[2]:
                        z = (mass-masses.mod_weights['nh2'][0]+masses.mod_weights['oh'][0]+hw+hw*tcharge)/tcharge
                        self.peakTable.append((z,'Z%d'%tcharge,'Z%d+%d'%(ionIndex,tcharge),v))
                    elif 'Z' in lossType[1]:
                        z = (mass-masses.mod_weights['nh2'][0]+masses.mod_weights['oh'][0]+hw+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((z,'Z%d'%tcharge,'Z%d%s+%d'%(ionIndex,lossType[2],tcharge),v))
                    
                    
#                    self.peakTable.append((y,'Y%d'%icharge,'Y%d%s+%d'%((sLen-i),lossType[2],icharge+1)))
#                    self.peakTable.append((z,'Z%d'%icharge,'Z%d%s+%d'%((sLen-i),lossType[2],icharge+1)))
#            mass-=masses.protein_weights[v][0]
#            charge-=masses.protein_weights[v][1]
#            try:
#                mass-=masses.mod_weights[modList[i+1]][0]
#                charge-=masses.mod_weights[modList[i+1]][1]
#            except KeyError:
#                pass
        self.peakTable.sort(key=operator.itemgetter(0))
        return self.peakTable
    
    def assignPeaks(self, x, y):
        """
        given a list of masses, returns what predicted peaks match up
        return format: [(m/z, intensity, aa, type, type description)]
        """
        self.predictPeaks()
        out = []
        start=0
        for seqindex,peaks in enumerate(self.peakTable):
            mz,pType, pDesc,aa = peaks
            candidates = []
            seen = False
            for index,m in enumerate(x[start:]):
                if m+self.tolerance > mz > m-self.tolerance:
                    seen = True
                    candidates.append((x[start+index],y[start+index],aa, pType, pDesc))
                elif mz > m+self.tolerance and seen:
                    start=index
                    break
            if candidates:
                candidates.sort(key=operator.itemgetter(1))
                out.append(candidates[0])
#        print out
        return out

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
        start = len(x)
        self.ySeries = []
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
            candidates = []
            seen = False
            for index,m in enumerate(x[:start]):
                if yi > m-self.tolerance and yi < m+self.tolerance:
                    seen = True
                    candidates.append((x[index],y[index],seq[seqindex], index))
                elif yi > m+self.tolerance and seen:
                    start=index
                    break
            if candidates:
                candidates.sort(key=operator.itemgetter(1))
                self.ySeries.append(candidates[0])
            else:
                self.ySeries.append(False)
        return self.ySeries
#for i in figureIons('AAAAK',2, '',10).predictPeaks():
#    print i