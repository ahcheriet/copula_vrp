
from HvEDA_m import *





test = HvEDA(100,FonctionObject01)
test.Evaluate()
test.Tofile("best",test.ListNonDominated,2)
print "Starting"
gplot = Gnuplot.Gnuplot()
for i in range(500):
    if (i%1 ==0):
        print "Iteration ",i ," ",len(test.ListNonDominated)
        d=Gnuplot.Data(test.ApplyFunctions(test.ListNonDominated))
        gplot.plot(d)
    test.Evaluate_After()
   # test.filtrate()
    
print "Iteration ",i ," ",len(test.ListNonDominated)
test.Tofile("best2",test.ListNonDominated,2)
Xvar = []
print "finale"
for i in range(m):
    XX, YY = test.C[i].generate_xy(500)
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
print "Filtre"
test.filtrate()
test.Tofile("best_final",test.ListNonDominated,2)

##-----------------------------------------------------------------------------
