ƒimport random
import sys

# print usage message if missing arguments
if len(sys.argv) <=1:
   print("Usage:")
   print(" $ python3 generatedata.py <ref_length> <nreads> <read_len> <ref_file> <reads_file>")
   sys.exit(0)



#---------This section of codes generate reference based on specfied length

# create a string for reference
ref = ''

# Ask for reference length
# input("Enter reference length:")
ref_length = int(sys.argv[1])

# reference length for randomly generated section
non_repeat_reflength =int(ref_length*7.5/10)

# reference length for repeated section
repeat_reflength = int(ref_length*2.5/10)

# randomly generated section of the reference
for i in range (0, non_repeat_reflength):
    ref =ref +str(random.randint(0,3))

# print out references in numbers
# print (ref)

# concatenate repeated section to randomly generated section
ref  = ref +ref[-repeat_reflength:]

#Convert reference from numbers to letters
ref = ref.replace('0','A').replace('1','C').replace('2','G').replace('3','T')

#print out references in letters and length
#print (ref)
#print(len(ref))


#-----This section of codes generate reads  based on specfied read length and number of reads

#read number of reads and read length from arguments
nreads = int(sys.argv[2])
read_len =  int(sys.argv[3])



