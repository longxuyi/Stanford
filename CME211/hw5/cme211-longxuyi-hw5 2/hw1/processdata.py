import random
import sys
import time

# print usage message if missing arguments
if __name__ == "__main__":
    if len(sys.argv) <=3:
        print("Usage:")
        print(" $ python3 processdata.py <ref_file> <reads_file> <align_file> ")
        sys.exit(0)

#read reference file and create a reference string
f = open(sys.argv[1],'r')
ref = f.readline()
f.close()

#read read file and create a list of reads
f = open(sys.argv[2],'r')
read = f.readlines()
f.close()

#Start the timer
start = time.time()

#First scan list of locations of alignment in reference
location_first = []
for i in range(len(read)):
    reads_loc =ref.find(read[i].strip('\n'), 0, int(len(ref)))
    location_first.append(reads_loc)
    

#Initialize the second list with same length as the first one
location_second = [" "]*len(location_first)

#Second scan list of locations of alignment in reference
for i in range(len(read)):
    if location_first[i] != -1:
        reads_loc =ref.find(read[i].strip('\n'), location_first[i]+1, int(len(ref)))
        location_second[i] = reads_loc

for i in range(len(read)):
    if location_second[i] == -1:
        location_second[i] = " "



#Initialize lists of reads that align once, twice, and zero times
read_one = 0
read_two = 0
read_zero = 0

# Count number of reads that align zero, one, and two times.
for i in range(len(read)):
    if location_first[i] != -1:
        read_one += 1
    elif location_first[i] == -1:
        read_zero += 1

for i in range(len(read)):
    if location_second[i] != " ":
         read_two += 1
         read_one -= 1
    if location_second[i] == " ":
         read_one += 0


#Calculate the eplapsed time
time = time.time() - start

#output on terminal
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
     #print(AlignList[i])

#generate an align file
with open(sys.argv[3], 'w') as f:
     for item in AlignList:
         f.write ('{}\n'.format(item))
        