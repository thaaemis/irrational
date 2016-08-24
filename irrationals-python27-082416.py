# Writing for Anaconda (python 2.7) for hdf5? or for JSON? Which is better?
# -- This program should take a string of elements as arguments
#    and output the resulting quantities:
#    -- Value
#    -- Tree of convergents
#    -- Closest convergent; last convergent
#    -- Distance to closest convergent
#    -- parent convergents

import sys, itertools, datetime, json

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

def outputJSON(CFs):
    
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
        value = x.alpha
        final = [x.finalConvergent.p, x.finalConvergent.q]
        closest = [x.closestConvergent.p, x.closestConvergent.q]
        dCrit = x.dCrit
        JSON = {"value":value, "finalConvergent":final,
                "closestConvergent":closest, "dCrit":dCrit}
        JSONs[CFsave]=JSON
    with open("CFs.json","w") as f:
        json.dump(JSONs,f,indent=1)

def main(nmax, length):

    iteration = itertools.product(range(1,nmax+1), repeat = length)
    results = []
    maxElement = 0
    for x in iteration:
        if max([i for i in x]) > maxElement:
            maxElement = max([i for i in x])
            print(datetime.datetime.now().isoformat(), maxElement)
        results.append(makeCF(x))
    outputJSON(results)
    
nmax, length = [int(x) for x in sys.argv[1:]]
main(nmax, length)
