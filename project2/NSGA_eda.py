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
#from pylab import *

from hv import *


##foo_clayton = Copula(x, y, family='clayton')
##foo_frank = Copula(x, y, family='frank')
##foo_gumbel = Copula(x, y, family='gumbel')


Contrainte01 = []
FonctionObject01 = []
Contrainte01.append('lambda x: True')
FonctionObject01.append('lambda x:x[0]')
f = lambda x:1+9*sum(map(lambda a:a/(m-1),x[1:]))
g = lambda x:(1-sqrt(x[0]/f(x) ))
FonctionObject01.append('lambda x:f(x)*g(x)')



Contrainte02 = []
FonctionObject02 = []
Contrainte02.append('lambda x: True')
FonctionObject02.append('lambda x:x[0]')
g2 = lambda x:(1-pow((x[0]/(f(x))),2))
FonctionObject02.append('lambda x:f(x)*g2(x)')

Contrainte03 = []
FonctionObject03 = []
Contrainte03.append('lambda x: True')
FonctionObject03.append('lambda x:x[0]')
g3 = lambda x:1-sqrt(x[0]/(f(x)))-(x[0]/(f(x)))*sin(10*pi*x[0])
FonctionObject03.append('lambda x:f(x)*g3(x)')



Contrainte04 = []
FonctionObject04 = []
Contrainte04.append('lambda x: True')
FonctionObject04.append('lambda x:x[0]')
FonctionObject04.append('lambda x:(1+10*(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))*(1-sqrt(x[0]/(1+10(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))))')

Contrainte05 = []
FonctionObject05 = []
Contrainte05.append('lambda x: True')
FonctionObject05.append('lambda x:x[0]')
FonctionObject05.append('lambda x:(1+10*(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))*(1-sqrt(x[0]/(1+10(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))))')

Contrainte06 = []
FonctionObject06 = []
Contrainte06.append('lambda x: True')
FonctionObject06.append('lambda x:1-exp(-4*x[0])*pow(sin(6*pi*x[0]),6)')
FonctionObject06.append('lambda x:(1+9*pow(sum(x[1:])/(m-1),0.25))*(1-pow( ( 1-exp(-4*x[0])*pow(sin(6*pi*x[0]),6) )/( 1+9*pow(sum(x[1:])/(m-1),0.25 )) ,2))')



m = 30 # nombre des sous-variables
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
                
    def isContrainte(self,X):
        Zvar = deepcopy(X)
        Indexs = []
        for index in range(len(Zvar)): # test du contraintes
#            if ( any(map(lambda x:True if (x < -5 or x >5) else False ,Zvar[index][1:])) == True ) or (Zvar[index][0]<0 or Zvar[index][0] > 1):     # pour f4
            if ( any(map(lambda x:True if (x < 0 or x >1) else False ,Zvar[index])) == True ):
                Indexs.append(Zvar[index])
        for index in range(len(Indexs)):
            Zvar.remove(Indexs[index])
        return Zvar
    
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
        for i in range(m):
#            ligne.append((np.random.uniform(-5,5,100)).tolist()) # 100 individu pour fonction4
            ligne.append((np.random.uniform(0,1,100)).tolist()) # 100 individu
        return zip(*ligne)
    
    def Make_X_UsingCopula2(self):
        Xvar = range(m)        # pour que la liste puisse contenir 30 elements
#        Varrtmp = self.Choice()
        varr = list(zip(*self.ListNonDominated))
        indexs = rd.sample(range(m),m)
        xx = rd.sample( self.ListNonDominated , len(self.ListNonDominated)/2) # random demi de meilleur sol
        yy = [item for item in self.ListNonDominated if item not in xx]
        if len(xx)>len(yy):
            xx.pop()            # test egality
        if len(yy)>len(xx):
            yy.pop()
        k = 0        
        for j in range(m/2):
            xm = deepcopy(varr[k])
            ym = deepcopy(varr[k+1])
            x = np.array(xm)
            y = np.array(ym)
            foo = Copula(x, y, family='frank') # gumbel, clayton, frank
            self.C[j] = deepcopy(foo)            
            XX, YY = foo.generate_xy(100)
            X1 = XX.tolist()
            Y1 = YY.tolist()
            Xvar[k] = deepcopy(X1)
            Xvar[k+1] = deepcopy(Y1)
            k = k + 2
        Zvar = zip(*Xvar)
        Indexs = []
        for index in range(len(Zvar)): # test du contraintes
      #      if ( any(map(lambda x:True if (x < -5 or x >5) else False ,Zvar[index][1:])) == True ) or (Zvar[index][0]<0 or Zvar[index][0] > 1):     # pour f4
            if ( any(map(lambda x:True if (x < 0 or x >1) else False ,Zvar[index])) == True ):
                Indexs.append(Zvar[index])
        for index in range(len(Indexs)):
            Zvar.remove(Indexs[index])
        return Zvar

    def Evaluate_After(self):  
        " evaluter avec Fast Sorting et Copula "
        self.ListNonDominated = self.Make_X_UsingCopula2() + self.ListNonDominated
        self.ListNonDominated , self.front_0 = self.Load_ListND(self.ListNonDominated, 50)
        
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


referencePoint = [6.0,6.0]
hv = HyperVolume(referencePoint)      
volume = []
test = HvEDA(100,FonctionObject02)
test.X = test.isContrainte(test.Make_First_Pop())
test.ListNonDominated , test.front_0 = test.Load_ListND(test.X, 30)
#test.Evaluate()
print "Starting"
gplot = Gnuplot.Gnuplot()
for i in range(300):
    if (i%1 ==0):
        print "Iteration ",i ," ",len(test.ListNonDominated)," ",len(test.front_0)
        d=Gnuplot.Data(test.ApplyFunctions(test.front_0))
        gplot.plot(d)
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
