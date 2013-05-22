from HvEDAlib import *
import pickle

import Gnuplot



Xvar = []
test = HvEDA(100,FonctionObject01)


filehandler = open('filename', 'r')   # load Copula from file
test.C  = pickle.load(filehandler) 
filehandler.close()

print "final"
for i in range(m):
    XX, YY = test.C[i].generate_xy(50)
    X1 = XX.tolist()
    Y1 = YY.tolist()
    Xvar.append(X1+Y1)
Zvar = zip(*Xvar)
Indexs = []
for index in range(len(Zvar)):
    if ( any(map(lambda x:True if (x < 0 or x >1) else False ,Zvar[index])) == True ):
        Indexs.append(Zvar[index])
for index in range(len(Indexs)):
    Zvar.remove(Indexs[index])
test.ListNonDominated = deepcopy(Zvar)
print "Filter"
test.filtrate()
test.Tofile("best_final",test.ListNonDominated,2)

