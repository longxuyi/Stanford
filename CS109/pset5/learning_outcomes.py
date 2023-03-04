import math
import numpy as np
import random
import copy

filename1 = "learningOutcomes.csv"
filename2 = "background.csv"

#read learningOutcomes.csv
f = open(filename1,'r')
outcomes = f.readlines()
f.close()

act1, act2 = [],[]
for outcome in outcomes:
    if outcome.split(",")[1] == "activity1":
        act1.append(int(outcome.split(",")[2]))
    if outcome.split(",")[1] == "activity2":
        act2.append(int(outcome.split(",")[2]))


#read background csv
f = open(filename2,'r')
background = f.readlines()
f.close()

less_1, avg_1, more_1 = [],[],[]
less_2, avg_2, more_2 = [],[],[]

count = 0
for line in background:
    if line.split(",")[1].strip() == "less":
        if outcomes[count].split(",")[1] == "activity1":
            less_1.append(int(outcomes[count].split(",")[2]))
        else:
            less_2.append(int(outcomes[count].split(",")[2]))
    if line.split(",")[1].strip() == "average":
        if outcomes[count].split(",")[1] == "activity1":
            avg_1.append(int(outcomes[count].split(",")[2]))
        else:
            avg_2.append(int(outcomes[count].split(",")[2]))
    if line.split(",")[1].strip() == "more":
        if outcomes[count].split(",")[1] == "activity1":
            more_1.append(int(outcomes[count].split(",")[2]))
        else:
            more_2.append(int(outcomes[count].split(",")[2]))
    count +=1

#define function to sample data and calculate p value
def resample(data,n):
    return np.random.choice(data, n, replace = True)

def p_value(activity1, activity2, n):

    sm1 = round(sum(activity1)/len(activity1),3)
    #sv1 = np.var(activity1, ddof=1)
    sm2 = round(sum(activity2)/len(activity2),3)
    #sv2 = np.var(activity2, ddof=1)

    mean1 =[]
    mean2 =[]

    for i in range(10000):
        subsample1 = random.sample(activity1, n)
        mean1.append(np.mean(subsample1))
        subsample2 = random.sample(activity2, n)
        mean2.append(np.mean(subsample2))

    #make universal sample
    all_samples = copy.deepcopy(activity1)
    all_samples.extend(activity2)

    observed_diff = abs(sm1 -sm2)

    #implement bootstraping and calculate p value
    count =0
    iter =10000
    for i in range(iter):
        resample1 = resample(all_samples, len(activity1))
        resample2 = resample(all_samples, len(activity2))
        diff = abs(np.mean(resample1) -np.mean(resample2))

        if diff >= observed_diff: count +=1
    return round(observed_diff,3),count/iter

diff_less, p_less = p_value(less_1,less_2,100)
print("Less experience:")
print("mean difference is {}, p value is {}".format(diff_less, p_less))

diff_avg, p_avg = p_value(avg_1,avg_2,100)
print("Average experience:")
print("mean difference is {}, p value is {}".format(diff_avg, p_avg))

diff_more, p_more = p_value(more_1,more_2,100)
print("More experience:")
print("mean difference is {}, p value is {}".format(diff_more, p_more))
