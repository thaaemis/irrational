# Writing for Anaconda (python 2.7) for hdf5? or for JSON? Which is better?
# -- This program should take a string of elements as arguments
#    and output the resulting quantities:
#    -- Value
#    -- Tree of convergents
#    -- Closest convergent; last convergent
#    -- Distance to closest convergent
#    -- parent convergents

import sys, itertools, datetime, json, time
import numpy as np
from pylab import *
sys.path.insert(0, '/u/bkraus/Python/')
import colormaps as cmaps

def extendCF(CF):
    output = [x for x in CF]
    for i in range(30):
        output.append(1)
    return output

def getValue(CF):
    val = 0
    for i in reversed(range(len(CF))):
        val = val + CF[i]
        val = val**-1
    return val

class convergent:
    
    def __init__(self, p, q, alpha):
        self.p = int(p)
        self.q = int(q)
        try:
            self.value = float(p)/float(q)
            self.dCrit = abs(alpha * float(self.q)**2 - float(self.p) * float(self.q)) # for k = 2!
        except ZeroDivisionError:
            self.value = None
            self.dCrit = 1.

    def __lt__(self, other):
        if self.dCrit < other.dCrit:
            return True
        else:
            return False
            
    def __repr__(self):
        return str(self.p)+'/'+str(self.q)

class continuedFraction:
    def __init__(self, alpha, CF, finalConvergent, closestConvergent, dCrit):
        self.alpha = float(alpha)
        self.CF = [int(x) for x in CF]
        self.finalConvergent = finalConvergent
        self.closestConvergent = closestConvergent
        self.dCrit = dCrit

def makeCF(CFinput):

    CForig = [x for x in CFinput]

    if CForig[-1] != 1:
        CForig[-1] -= 1
        CForig.append(1)
        
    CF = [int(x) for x in CForig]
    value = getValue(extendCF(CF))
    
    # Chain of convergents
    convergents = [convergent(0,1,value)]
    if value > 0.5:
        convergents = [convergent(0,1,value), convergent(1,1,value)]
        CF.pop(0)
    else:
        convergents = [convergent(1,0,value), convergent(0,1,value)]
    for i in range(len(CF)):
        pNew = CF[i] * convergents[i+1].p + convergents[i].p
        qNew = CF[i] * convergents[i+1].q + convergents[i].q
        convergents.append(convergent(pNew, qNew, value))
    print([x for x in convergents])
    closestConvergent = min(convergents)
    return continuedFraction(value, CForig, 
           convergents[-1], closestConvergent, closestConvergent.dCrit)

def outputCSV(CFs):
    # antiquated: built for MySQL
    def printItem(item):
        string = '"' + str(item) + '",'
        return string
        
    output = ''
    for x in CFs:
        output += printItem(x.alpha)
        for i in range(30):
            try:
                output += printItem(x.CF[i])
            except IndexError:
                output += printItem(1)
        output += printItem(x.finalConvergent.p)
        output += printItem(x.finalConvergent.q)        
        output += printItem(x.closestConvergent.p)            
        output += printItem(x.closestConvergent.q)
        output += printItem(x.dCrit) + '\n'
    with open('CFs.csv','w') as f:
        f.write(output)

def outputJSON(CFs,length):
    
    JSONs = {} 
    for x in CFs:
        lastInd = 0 # last element that isn't 1
        for i, element in enumerate(x.CF):
            if element != 1:
                lastInd = i
        try:
            CFsave = x.CF[:lastInd+1]
        except IndexError:
            CFsave = []
        CFsave = ''.join([str(int(b)) for b in CFsave])
        for i in range(len(CFsave),length):
            CFsave += '1'
        value = x.alpha
        final = [x.finalConvergent.p, x.finalConvergent.q]
        closest = [x.closestConvergent.p, x.closestConvergent.q]
        dCrit = x.dCrit
        JSON = {"value":value, "finalConvergent":final,
                "closestConvergent":closest, "dCrit":dCrit}
        JSONs[CFsave]=JSON
    with open("CFs.json","w") as f:
        json.dump(JSONs,f,indent=1)

def makeOutput(nmax, length):

    iteration = itertools.product(range(1,nmax+1), repeat = length)
    results = []
    maxElement = 0
    for x in iteration:
        if max([i for i in x]) > maxElement:
            maxElement = max([i for i in x])
            print(datetime.datetime.now().isoformat(), maxElement)
        results.append(makeCF(x))
    outputJSON(results,length)
    
def output():
    nmax, length = [int(x) for x in sys.argv[1:]]
    makeOutput(nmax, length)

def retrieveInfo():

    def sortdCrit(CFs):
        CFlist, dCritList, valList, convergents = [], [], [], []
        for CF,obj in CFs.items():
           CFlist.append(str(CF))
           dCritList.append(obj['dCrit'])
           value = obj['value']
           p, q = obj['finalConvergent']
           valList.append(value)
           convergents.append(convergent(p, q, value))
        sortInd = [x[1] for x in sorted((e,i) for i,e in enumerate(dCritList))]
        sortInd.reverse()
        CFlist = np.asarray(CFlist)[sortInd]
        dCritList = np.asarray(dCritList)[sortInd]
        valList = np.asarray(valList)[sortInd]
        convergents = np.asarray(convergents)[sortInd]
        return CFlist, dCritList, valList, convergents

    with open('CFs66.json','r') as f:
        CFs6 = json.load(f)
    
    filtCF = {k:v for k,v in CFs6.iteritems() if (k[:3] == '564') and \
              (k[-1] == '3')}
    
    # lastly, sort by dCrit and print
    CFlist, dCritList, valList, convergents = sortdCrit(filtCF)
    for i,num in enumerate(CFlist):
        print(num, dCritList[i]) 
    dCritMax, dCritMin = dCritList[0],dCritList[-1]
    cmap = cmaps.plasma
    for i, value in enumerate(valList):
        scale = (dCritList[i]-dCritMin)/(dCritMax-dCritMin)
        plot(value, convergents[i].q, 'o',color=cmap(scale))
    yscale('log')
    show()
retrieveInfo()
