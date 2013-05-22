from math import *
from copy import *
from threading import *
import random as rd
from decimal import *
from scipy.misc.common import *
import os
import numpy as np
#import Gnuplot
from openturns import *
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
    def Dominate( self ,i): # les solutions dominer par i
        j = 0
        while j < len(self.X):
            if ( i != j ):
                k = 0
                CountI = 0
                CountJ = 0        
                while k < len ( self.fonctionObject ):                    
                      if  self.IsDominate( self.fonctionObject[k], self.X[i], self.X[j],1): # minimisation problem
                        CountI = CountI + 1
                      if  self.IsDominate( self.fonctionObject[k], self.X[j], self.X[i],1):
                        CountJ = CountJ + 1                               
                      k = k +1
                if CountI == len( self.fonctionObject ):
                    self.ListDominated.append(deepcopy(self.X[j]))
                if CountJ == len( self.fonctionObject ):
                    self.ListDominated.append(deepcopy(self.X[i]))
            j = j + 1
            
    def Choice(self):
        List_Sorted = []
        list_one = map(eval(self.fonctionObject[0]),self.ListNonDominated)#fonction1)
        list_two = sorted(list_one)
        for i in range(len(self.ListNonDominated)):
            index_one = list_one.index(list_two[i])
            element = deepcopy(self.ListNonDominated[index_one])
            List_Sorted.append(element)
        return List_Sorted
    def isContrainte(self,X):
        Zvar = deepcopy(X)
        Indexs = []
        for index in range(len(Zvar)): # test du contraintes
            if ( any(map(lambda x:True if (x < 0 or x >1) else False ,Zvar[index])) == True ):
                Indexs.append(Zvar[index])
        for index in range(len(Indexs)):
            Zvar.remove(Indexs[index])
        return Zvar

    def Evaluate(self):        
        i = 0
        self.X = []
        self.ListNonDominated = []
        self.X = self.isContrainte(self.Make_First_Pop())
        j = len(self.X)-1
        self.Dominate(j)
        i = 0
        while i < len(self.X):
            if self.X[i] not in self.ListDominated:
                self.ListNonDominated.append(deepcopy(self.X[i]))
            i = i +1
    def filtrate(self):
        i = 0
        while i < len(self.ListNonDominated)-1:
                j = i+1
                while j < len(self.ListNonDominated):
                        CountI = 0
                        CountJ = 0
                        k = 0
                        while k < len ( self.fonctionObject ):                    
                            if  self.IsDominate( self.fonctionObject[k], self.ListNonDominated[i], self.ListNonDominated[j],1):
                                CountI = CountI + 1
                            if  self.IsDominate( self.fonctionObject[k], self.ListNonDominated[j], self.ListNonDominated[i],1):
                                CountJ = CountJ + 1                               
                            k = k +1
                        if CountI == len( self.fonctionObject ):
                                x = self.ListNonDominated[j]
                                j = j - 1
                                self.ListNonDominated.remove(x)
                        if CountJ == len( self.fonctionObject ):
                                x = self.ListNonDominated[i]
                                j = j - 1
                                i = i - 1
                                self.ListNonDominated.remove(x)
                                break
                        j = j + 1
                i = i +1
    def Make_First_Pop(self):
        ligne = []
        for i in range(m):
            ligne.append((np.random.normal(0.5,0.2,100)).tolist()) # 100 individu
        return zip(*ligne)
    
    def Make_X_UsingCopula2(self):  # look at the filter
        Xvar = range(30)
        snd = NumericalSample(0,30)  # create a ensemble to enter nondimanated set0
        for i in range(len(self.ListNonDominated)):
            snd.add(deepcopy(self.ListNonDominated[i]))
        factory = NormalCopulaFactory()
        estimatedDistribution = factory.build(snd) #create a normal copula
        list_new = estimatedDistribution.getSample(100) # sampling
        # pour que la liste puisse contenir 30 elements
#        Varrtmp = self.Choice()
        Zvar = []
        for i in range(100):
            Zvar.append(list(list_new[i]))   
        print Zvar 
        Indexs = []
        for index in range(len(Zvar)): # test du contraintes
            if ( any(map(lambda x:True if (x < 0 or x >1) else False ,Zvar[index])) == True ):
                Indexs.append(Zvar[index])
        for index in range(len(Indexs)):
            Zvar.remove(Indexs[index])
        return Zvar

    def Evaluate_After(self):  
        " evaluter avec hypervolume et Copula "
        referencePoint = [6.0,6.0]
        ListNonDominatedtmp = []
        hv = HyperVolume(referencePoint)      
        front = self.ApplyFunctions(self.ListNonDominated)
        volume1 = hv.compute(front)
        self.indexbest = 0
        self.ListNonDominated = self.Make_X_UsingCopula2() + self.ListNonDominated
        self.filtrate()
        if ( len(self.ListNonDominated) > 100 ):
            self.ListNonDominated = rd.sample( self.ListNonDominated , 80)
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
test = HvEDA(100,FonctionObject01)
test.Evaluate()
test.Tofile("best",test.ListNonDominated,2)
print "Starting"
#gplot = Gnuplot.Gnuplot()
for i in range(20):
    if (i%1 ==0):
        print "Iteration ",i ," ",len(test.ListNonDominated)
#        d=Gnuplot.Data(test.ApplyFunctions(test.ListNonDominated))
#        gplot.plot(d)
    test.Evaluate_After()
    front = test.ApplyFunctions(test.ListNonDominated)
    volume.append (hv.compute(front))
   # test.filtrate()

fileData = open("hypervolume1.txt","w")            
i = 0
while i < len( volume ):
    fileData.write(str(i)+'\t'+str(volume[i])+'\n')
    i = i +1
fileData.close()            






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
