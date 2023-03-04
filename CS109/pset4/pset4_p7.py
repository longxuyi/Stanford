import sys
import math
import numpy as np
from scipy import stats

fa = "personKeyTimingA.txt"
fb = "personKeyTimingB.txt"
fe = "email.txt"

def mean_var_calculator(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close

    time = []
    time.append(float(lines[0].split(",")[0]))

    for i in range(1, len(lines[1:])):
        time.append(float(lines[i].split(",")[0]) - float(lines[i-1].split(",")[0]))

    mean = np.mean(time) 
    var = np.var(time)
    
    return time, mean, var

timeA, meanA, varA = mean_var_calculator(fa)
timeB, meanB, varB = mean_var_calculator(fb)
timeE, meanE, varE = mean_var_calculator(fe)
print ("User A mean:{}, var: {}".format(round(meanA,3), round(varA,3)))
print ("User B mean:{}, var: {}".format(round(meanB,3), round(varB,3)))

odds =1
for i in range(len(timeE)):
    temp = stats.norm.pdf(timeE[i], meanA, math.sqrt(varA))/stats.norm.pdf(timeE[i], meanB, math.sqrt(varB))
    odds *= temp

print("odds: ", odds)