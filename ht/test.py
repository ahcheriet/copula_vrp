import numpy as np
from copy import deepcopy
from math import  sqrt
Num_Indv = 15
m = 5
Demand = [2,1,2,1,1]
Capacity = 4
coordinate =  [[3,2],[2,3],[4,2],[2,6],[1,3],[3,2],[2,3],[4,2],[2,6],[1,3],[3,2],[2,3],[4,2],[2,6],[1,3],[3,2],[2,3],[4,2],[2,6],[1,3],[3,2],[2,3],[4,2],[2,6],[1,3],[3,2],[2,3],[4,2],[2,6],[1,3]]

X = []
Result = []
fct1 = lambda x:len([ a for a in x if a<>[]]) # number of routes
Distance = lambda a,b:sqrt(sum( (a - b)**2 for a, b in zip(a, b))) #distance between a and b
fct2 = lambda a:sum([sum([Distance(coordinate[a[j][i]],coordinate[a[j][i+1]]) for i in range(len(a[j])-1)]) for j in range(len(a)) ])

def Make_First_Pop():
    ligne = []
    i = 0
    while i < m*Num_Indv:  # all individuals
        ligne.append((np.random.uniform(0,1,m)).tolist())
        i = i + 1 
    return ligne
def Lecture():
    Result = []
    for i in range(Num_Indv):
        tmpList=[[] for k in range(m)]
        Indv = zip(*X[(i*m):(m*(i+1))])  # Extracting of Individuals
        for j in range(m):
            tmpList[np.argmax(Indv[j])].append(j) # insert Routes ?!
        Result.append(deepcopy(tmpList))
    return Result
def InConstraine():
    Accepted = []
    for i in range(Num_Indv):
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
    return Accepted


X = Make_First_Pop()
Result = Lecture()
a = [[0, 3, 5], [], [13], [], [], [], [], [11], [1, 23, 26], [28], [18], [], [9, 19], [27], [8], [], [], [12], [14], [4, 16, 24], [2, 
 6], [7], [], [15], [21], [], [17, 20, 25], [10], [29], [22]]
print fct2(a)
print InConstraine()
