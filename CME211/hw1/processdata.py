import random
import sys
import time

#read reference file and create a reference string
f = open("reference.txt",'r')
ref = f.readline()
f.close()

#read read file and create a list of reads
f = open("read.txt",'r')
read = f.readlines()
f.close()

#Start the timer
start = time.time()


#First scan list of locations of alignment in reference
location_first = []
for i in range(len(read)):
    reads_loc =ref.find(read[i].strip('\n'), 0, int(len(ref)))
    location_first.append(reads_loc)
    

#Second scan list of locations of alignment in reference

#create second list with same length as first list
location_second = [" "]*len(location_first)

for i in range(len(read)):
    if location_first[i] != -1:
        reads_loc =ref.find(read[i].strip('\n'), location_first[i]+1, int(len(ref)))
        location_second[i] = reads_loc

for i in range(len(read)):
    if location_second[i] == -1:
        location_second[i] = " "



#print alignment positions for each read
# for i in range(len(read)):
#     print("{} {} {}". format(read[i].strip('\n'),location_first[i], location_second[i]))

#calculate percentage of reads with none, one, and mulitple alignments 

#Initialize lists of reads that align once, twice, and zero times
read_one = 0
read_two = 0
read_zero = 0


for i in range(len(read)):
    if location_first[i] != -1:
        read_one += 1
        #read_two += 1
    elif location_first[i] == -1:
        read_zero += 1

for i in range(len(read)):
    if location_second[i] != " ":
         read_two += 1
         read_one -= 1
    if location_second[i] == " ":
         read_one += 0

# print(read_two)
# print(read_one)
# print(read_zero)



#Calculate the eplapsed time
time = time.time() - start


print ("Reference Length: {}".format(len(ref)))
print ("Number Reads: {}".format(len(read)))
print ("Read Length: {}".format(len(read[0])))
print ("Align 0: {}".format(read_zero/len(read)))
print ("Align 1: {}".format(read_one/len(read)))
print ("Align 2: {}".format(read_two/len(read)))
print("elapsed time: {}".format(time))


# Create lists of alignment positions for each read

AlignList =[" "]*len(location_first)
for i in range(len(read)):
     AlignList[i] = "{} {} {}". format(read[i].strip('\n'),location_first[i], location_second[i])
     print(AlignList[i])

#generate a align file
with open('alignfile.txt', 'w') as f:
     for item in AlignList:
         f.write ('{}\n'.format(item))