import matplotlib.pyplot as plt
import numpy as np
import random

filename1 = "peerGrades.csv"

#read learningOutcomes.csv
f = open(filename1,'r')
grades = f.readlines()
f.close()

grades = [int(x.strip()) for x in grades]

def resample(data,n):
    return np.random.choice(data, n, replace = True)

n=5
mean = []
median =[]

for i in range (10000):
    mean.append(np.mean(random.sample(grades,n)))
    median.append(np.median(random.sample(grades,n)))
print("variance of mean: {}".format(np.var(mean)))  
print("variance of median: {}".format(np.var(median)))   

#plot histogram
n_col =2
fig, axes= plt.subplots(n_col)

#print(axes[2])


axes[0].hist(mean, bins = int(max(mean) - min(mean)) )
axes[0].set_title("mean distribution")
axes[1].hist(median, bins = int(max(median) - min(median)) )
axes[1].set_title("median distribution")

plt.show()
