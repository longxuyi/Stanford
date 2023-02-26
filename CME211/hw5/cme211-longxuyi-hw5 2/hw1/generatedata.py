import random
import sys

# print usage message if missing arguments
if __name__ == "__main__":
    if len(sys.argv) <= 5:
        print("Usage:")
        print(" $ python3 generatedata.py <ref_length> <nreads> <read_len> <ref_file> <reads_file>")
        sys.exit(0)


# Ask for reference length, number of reads, and read length
ref_length = int(sys.argv[1])
nreads = int(sys.argv[2])
read_len = int(sys.argv[3])

# create a string for reference
ref = ''

# reference length for randomly generated section
non_repeat_reflength =int(ref_length*7.5/10)

# reference length for repeated section
repeat_reflength = int(ref_length*2.5/10)

# randomly generated section of the reference
for i in range (0, non_repeat_reflength):
    ref =ref +str(random.randint(0,3))


# concatenate repeated section to randomly generated section
ref  = ref +ref[-repeat_reflength:]

#Convert reference from numbers to letters
ref = ref.replace('0','A').replace('1','C').replace('2','G').replace('3','T')

#Initialize lists of reads that align once, twice, and zero times
read_one = []
read_two = []
read_zero =[]


# generate lists of reads that align once, twice, and zero times
for i in range (nreads):
    randn = random.random()

     #generate a list of eads align once   
    if randn >0.15 and randn <=0.9:
        rand_start = random.randint(0,int(ref_length/2))
        read =ref[rand_start:rand_start+read_len-1]
        read_one.append(read)
    
    #generate a list of reads with two align
    elif  randn <= 0.1:
        rand_start = random.randint(int(ref_length*0.75), ref_length-read_len)
        read =ref[rand_start:rand_start+read_len-1]
        read_two.append(read) 


    #generate a list of reads align twzero times
    else:
        #initialize a value for noAlign
        noAlign = ref[0:read_len]
        while ref.find(noAlign) != -1:
            noAlign = ''

            for i in range (read_len):
                noAlign = noAlign + str(random.randint(0,3))
            noAlign = noAlign.replace('0','A').replace('1','C').replace('2','G').replace('3','T')

        read_zero.append(noAlign)
            

#print outputs in terminal
print ("Reference Length: {}".format(ref_length))
print ("Number Reads: {}".format(nreads))
print ("Read Length: {}".format(read_len))
print ("Align 0: {}".format(len(read_zero)/nreads))
print ("Align 1: {}".format(len(read_one)/nreads))
print ("Align 2: {}".format(len(read_two)/nreads))

#create a reference.txt file with a line of reference
with open(sys.argv[4], 'w') as f:
    f.write (ref)


#create a r.txt file with a slit of all reads
read = read_one + read_two +read_zero
with open(sys.argv[5], 'w') as f:
    for item in read:
        f.write ('{}\n'.format(item))