#
#    Utilisation d Fast Sorting et permutation de X , Y
#
from math import *
from copy import *
from threading import *
import random as rd
from decimal import *
from scipy.misc.common import *
import os
import numpy as np
import Gnuplot

from copulalib.copulalib import Copula

from hv import *



Contrainte01 = []
FonctionObject01 = []
Distance = lambda a,b:sqrt(sum( (a - b)**2 for a, b in zip(a, b))) #distance between a and b
fct1 = lambda x:len([ a for a in x if a<>[]]) # number of routes
fct2 = lambda a:sum([sum([Distance(coordinate[a[j][i]],coordinate[a[j][i+1]]) for i in range(len(a[j])-1)]) for j in range(len(a)) ])
FonctionObject01.append('lambda x:fct1(x)')
FonctionObject01.append('lambda x:fct2(x)')
#Contrainte01.append('lambda x: True')


#Construction of VRP Potential Solutions

m = 30 # nombre des sous-variables
Num_Indv = 10 # number of Individuals
Demand = range(m)
Capacity = 100
Result = [] # Delete this because its the result of lecture
coordinate = [[3,2],[2,3],[4,2],[2,6],[1,3],[3,2],[2,3],[4,2],[2,6],[1,3],[3,2],[2,3],[4,2],[2,6],[1,3],[3,2],[2,3],[4,2],[2,6],[1,3],[3,2],[2,3],[4,2],[2,6],[1,3],[3,2],[2,3],[4,2],[2,6],[1,3]]

class HvEDA:
    def __init__(self,nbits,fonction):
        " testing HV With EDA "
        self.ind = []
        self.C = range(m)
        self.ListDominated = []
        self.fonctionObject = fonction
        self.X = []
        self.Tbit = nbits
    def IsDominate( self,Lfonction, x, y , MinOrMax):
        fonction = eval(Lfonction)
        if MinOrMax == 1:
            if min(fonction(x),fonction(y)) == fonction(x):
                return True
            else:
                return False
        else:
            if max(fonction(x),fonction(y)) == fonction(x):
                return True
            else:
                return False           
                
    def Fast_Sorting(self,X):
        tmpX = deepcopy(X)
        n = [-1]* len(tmpX)
        S = [0]* len(tmpX)
        rank = [-1]* len(tmpX)
        F = [[]]
        for x_i in range(len(tmpX)):
            S[x_i] = []
            n[x_i] = 0
            for x_j in  range(len(tmpX)):
                if (x_j != x_i):
                    k = 0
                    CountI = 0
                    CountJ = 0        
                    while k < len ( self.fonctionObject ):                    
                        if  self.IsDominate( self.fonctionObject[k], tmpX[x_i], tmpX[x_j],1): # minimisation problem
                            CountI = CountI + 1
                        if  self.IsDominate( self.fonctionObject[k], tmpX[x_j], tmpX[x_i],1):
                            CountJ = CountJ + 1
                        k = k +1
                    if CountI == len( self.fonctionObject ):
                        S[x_i].append(x_j)
                    if CountJ == len( self.fonctionObject ):
                        n[x_i] = n[x_i] + 1
            if n[x_i] == 0:
                rank[x_i] = 0
                F[0]=F[0]+[x_i]
        i = 0
        ri = 1
        while F[i] != []:
            Q = []
            for p in F[i]:
                for q in S[p]:
                    n[q] = n[q]-1
                    if n[q] == 0:
                        rank[q] = ri + 1
                        Q.append(q)
            i = i + 1
            ri = ri + 1
            F.append(deepcopy(Q))
        return F
    
    def Load_ListND(self,X,n):
        tmpND = []
        f = self.Fast_Sorting(X)
        for i in range(len(f)):
            for j in range(len(f[i])):
                tmpND.append(X[f[i][j]])
                if len(tmpND) == n:
                    return tmpND ,tmpND[:len(f[0])]
        return tmpND , tmpND[:len(f[0])]

    def Make_First_Pop(self):
        ligne = []
        i = 0
        while i < m*Num_Indv:  # all individuals
            ligne.append((np.random.uniform(0,1,m)).tolist())
            i = i + 1
        self.P = deepcopy(ligne) 
        return ligne
    
    def Lecture(self,X):
        Result = []
        for i in range(len(X)/m):
            tmpList=[[] for k in range(m)]
            Indv = zip(*X[(i*m):(m*(i+1))])  # Extracting of Individuals
            for j in range(m):
                tmpList[np.argmax(Indv[j])].append(j) # insert Routes ?!
            Result.append(deepcopy(tmpList))
        self.before_contrainte = deepcopy(Result)
        return Result
    
    def isContrainte(self,Result):
        Accepted = []
        for i in range(len(Result)):
            boolean = True
            PAccepted = Result[i]
            for j in range(m):
                som = 0
                for x in PAccepted[j]:
                    som = Demand[x]+som
                if som > Capacity:
                    boolean = False
            if boolean :
                Accepted.append(i)
        self.P = self.Load_PModel(Accepted)
        self.Constrained = [ Result[i] for i in Accepted]
        return self.Constrained
 
    def Load_PModel(self,indexes):
        res = []
        for i in indexes:
            res = res + self.P[(i*m):((i+1)*m)]
        return res

    def Make_X_UsingCopula2(self):
        indexes = [self.Constrained.index(b) for b in self.ListNonDominated ]        
        self.ListNonDominatedP = self.Load_PModel(indexes)
        Xvar = []
        xx = rd.sample( self.ListNonDominatedP , len(self.ListNonDominatedP)/2) # demi de meilleur sol and random
        yy = [item for item in self.ListNonDominatedP if item not in xx]
        if len(xx)>len(yy):
            xx.pop()            # test egality
        if len(yy)>len(xx):
            yy.pop()        
        for j in range(m):
            xm = []
            ym = []
            for i in range(len(yy)):   
                xm.append(deepcopy(xx[i][j]))
                ym.append(deepcopy(yy[i][j]))
            x = np.array(xm)
            y = np.array(ym)
            foo = Copula(x, y, family='frank') # gumbel, clayton, frank
            self.C[j] = deepcopy(foo)            
            XX, YY = foo.generate_xy(150)
            X1 = XX.tolist()
            Y1 = YY.tolist()
            Xvar.append(X1+Y1)
        Zvar = zip(*Xvar)
        Zvar = Zvar + self.ListNonDominatedP
        self.P = Zvar
        new = self.isContrainte(self.Lecture(Zvar))
  #      self.P = self.P + self.ListNonDominatedP
        return new


    def Evaluate_After(self):  
        " evaluter avec Fast Sorting et Copula "
        self.ListNonDominated = self.Make_X_UsingCopula2() + self.ListNonDominated
        self.ListNonDominated , self.front_0 = self.Load_ListND(self.ListNonDominated, 10)
        
    def ApplyFunctions(self,ListX):
        tmpResult = []
        for i in range(len( ListX )):
            fonction1 = eval(self.fonctionObject[0])
            fonction2 = eval(self.fonctionObject[1])
            tmpResult.append([fonction1(ListX[i]),fonction2(ListX[i])])
        return tmpResult
    
    def Tofile( self , Str_File , X  , AorW):
        if AorW == 1:
            fileData = open(Str_File+".txt","a")
        else:
            fileData = open(Str_File+".txt","w")            
        i = 0
        while i < len( X ):
            fonction1 = eval(self.fonctionObject[0])
            fonction2 = eval(self.fonctionObject[1])
            fileData.write(str(fonction1(X[i]))+'\t'+str(fonction2(X[i]))+'\n')
            i = i +1
        fileData.close()            


referencePoint = [60.0,60.0]
hv = HyperVolume(referencePoint)      
volume = []
test = HvEDA(100,FonctionObject01)
test.X = test.isContrainte(test.Lecture(test.Make_First_Pop()))
test.ListNonDominated , test.front_0 = test.Load_ListND(test.X, 30)


#########################
# main program running  #
#########################
print "Starting"
gplot = Gnuplot.Gnuplot()
for i in range(30):
    if (i%1 ==0):
        print "Iteration ",i ," ",len(test.ListNonDominated)," ",len(test.front_0)
   #     d=Gnuplot.Data(test.ApplyFunctions(test.front_0))
   #     gplot.plot(d)
    test.Evaluate_After()
    front = test.ApplyFunctions(test.front_0)
    volume.append (hv.compute(front))
   # test.filtrate()

fileData = open("NSGA_eda.txt","w")            
i = 0
while i < len( volume ):
    fileData.write(str(i)+'\t'+str(volume[i])+'\n')
    i = i +1
fileData.close()            




print "Iteration ",i ," ",len(test.ListNonDominated)
test.Tofile("best2",test.front_0,2)
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
#test.filtrate()
test.Tofile("best_final",test.ListNonDominated,2)

##-----------------------------------------------------------------------------
