from HvEDAlib import *
import pickle

from multiprocessing import Process,Pool


import Gnuplot



Xvar = []
test = HvEDA(100,FonctionObject01)
test.Evaluate()


def Generation(gen):
    print "Starting"
    for i in range(gen):
        print "Iteration ",i ," ",len(test.ListNonDominated)
        test.Evaluate_After()
    
print "Starting"
gplot = Gnuplot.Gnuplot()
for i in range(10):
    if (i%5 ==0):
        print "Iteration ",i ," ",len(test.ListNonDominated)
        d=Gnuplot.Data(test.ApplyFunctions(test.ListNonDominated))
        gplot.plot(d)
    test.Evaluate_After()


#job = []
#job.append(Process(target=Generation, args=(30,)))
#job.append(Process(target=Generation, args=(30,)))
#
#for j in job:
#    j.start()
#for j in job:
#    j.join()


    

test.Tofile("best",test.ListNonDominated,2)
print "Terminated"
filehandler = open('filename', 'w')  # Save Copula to 
pickle.dump(test.C, filehandler) 
filehandler.close()

