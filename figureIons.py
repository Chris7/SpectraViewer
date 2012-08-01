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
        """
        returns a list of peaks predicted from a peptide:
        format sorted by low->high mz:
        list of tuples of: [(m/z, fragment type(A,B,C,X,Y,Z), fragment #, charge state,neutral loss, peptide),(...)]
        """
        mass = 0
        self.peakTable = []
        modList = self.modList
        maxcharge = int(self.ioncharge)-1
        hw = masses.mod_weights['h'][0]
        losses = [(0,('a','b','c','x','y','z'),'')]
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
            if charge>maxcharge:
                charge=maxcharge
            for icharge in xrange(1,charge+1):
                for lossType in losses:
                    lossMass = lossType[0]
#                    print 'abc loss type',v,lossType,lossMass,i
                    tcharge = icharge
                    if 0<i<sLen-1 and not lossType[2]:
                        a = (mass-masses.mod_weights['cho'][0]+hw+hw*tcharge)/tcharge
                        self.peakTable.append((a,'a',ionIndex,tcharge,None,v))
                    elif 0<i<sLen-1 and 'a' in lossType[1]:
                        a = (mass-masses.mod_weights['cho'][0]+hw+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((a,'a',ionIndex,tcharge,lossType[2],v))
                    if 0<i<sLen-1 and not lossType[2]:
                        b = (mass+masses.mod_weights['h'][0]-hw+hw*tcharge)/tcharge
                        self.peakTable.append((b,'b',ionIndex,tcharge,None,v))
                    elif 0<i<sLen-1 and 'b' in lossType[1]:
                        b = (mass+masses.mod_weights['h'][0]-hw+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((b,'b',ionIndex,tcharge,lossType[2],v))
                    if i<sLen-1 and not lossType[2]:
                        c = (mass+masses.mod_weights['nh2'][0]+hw+hw*tcharge)/tcharge
                        self.peakTable.append((c,'c',ionIndex,tcharge,None,v))
                    elif i<sLen-1 and 'c' in lossType[1]:
                        c = (mass+masses.mod_weights['nh2'][0]+hw+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((c,'c',ionIndex,tcharge,lossType[2],v))
        charge=int(self.ioncharge)-1#for rightmost charge -- we remove the NH3 charge
        losses = [(0,('a','b','c','x','y','z'),'')]
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
            if charge>maxcharge:
                charge=maxcharge
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
                    tcharge=icharge
                    if not lossType[2]:
                        x = (mass+masses.mod_weights['co'][0]+masses.mod_weights['oh'][0]-hw+hw*tcharge)/tcharge
                        self.peakTable.append((x,'x',ionIndex,tcharge,None,v))
                    elif 'x' in lossType[1]:
                        x = (mass+masses.mod_weights['co'][0]+masses.mod_weights['oh'][0]-hw+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((x,'x',ionIndex,tcharge,lossType[2],v))
                    if not lossType[2]:
                        y = (mass+masses.mod_weights['h'][0]+masses.mod_weights['oh'][0]+hw*tcharge)/tcharge
                        self.peakTable.append((y,'y',ionIndex,tcharge,None,v))
                    elif 'y' in lossType[1]:
                        y = (mass+masses.mod_weights['h'][0]+masses.mod_weights['oh'][0]+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((y,'y',ionIndex,tcharge,lossType[2],v))
                    if not lossType[2]:
                        z = (mass-masses.mod_weights['nh2'][0]+masses.mod_weights['oh'][0]+hw+hw*tcharge)/tcharge
                        self.peakTable.append((z,'z',ionIndex,tcharge,None,v))
                    elif 'z' in lossType[1]:
                        z = (mass-masses.mod_weights['nh2'][0]+masses.mod_weights['oh'][0]+hw+hw*tcharge+lossMass)/tcharge
                        self.peakTable.append((z,'z',ionIndex,tcharge,lossType[2],v))
        self.peakTable.sort(key=operator.itemgetter(0))
        return self.peakTable
    
    def assignPeaks(self, x, y):
        """
        given a list of masses, returns what predicted peaks match up
        return format: [(m/z, intensity, amino acid, fragmentation type, fragment #, charge, neutral loss)]
        """
        self.predictPeaks()
        out = []
        start=0
        primaries = set([])
        for seqindex,peaks in enumerate(self.peakTable):
            mz,fragType, fragNum,charge,loss,aa = peaks
            candidates = []
            seen = False
            for index,m in enumerate(x[start:]):
                if m+self.tolerance > mz > m-self.tolerance:
                    seen = True
                    candidates.append((x[start+index],y[start+index],fragType, fragNum,charge,loss,aa))
                elif mz > m+self.tolerance and seen:
                    start=index
                    break
            if candidates:
                candidates.sort(key=operator.itemgetter(1))
                toadd = candidates[0]
                if not toadd[5]:
                    primaries.add((toadd[2:5]))
                out.append(candidates[0])
        #now we filter our secondary ions if the primary ion is missing
        out[:] = [x for x in out if not x[5] or (x[5] and x[2:5] in primaries)]
        return out