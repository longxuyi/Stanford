CME211 hw1 

Part 2 Generate Data

Command line log

$ python3 generatedata.py 1000 600 50 "reference_1.txt" "reads_1.txt"
Reference Length: 1000
Number Reads: 600
Read Length: 50
Align 0: 0.135
Align 1: 0.7716666666666666
Align 2: 0.09333333333333334
$ python3 generatedata.py 10000 6000  50 "reference_2.txt" "reads_2.txt"
Reference Length: 10000
Number Reads: 6000
Read Length: 50
Align 0: 0.14816666666666667
Align 1: 0.7501666666666666
Align 2: 0.10166666666666667

$ python3 generatedata.py 100000  60000  50 "reference_3.txt" "reads_3.txt"
Reference Length: 100000
Number Reads: 60000
Read Length: 50
Align 0: 0.15026666666666666
Align 1: 0.7495333333333334
Align 2: 0.1002


//Describe considerations you used in designing your handwritten data.
 For manually created data, the length of the reference is 10 and there are 5 reads of length 3. To thoroughly check the code, we manually create 3 reads that align once, 1 that aligns twice, and 1 that has zero alignment. The program should work correctly with larger datasets as well.

//Should you expect an exact 15\% 75\% 10\% distribution for the reads that align 0,1,2 times? What other factors will affect the exact distribution of the reads

No. The distribution is never exact, because the distribution is generated using random.randint() function. The Size of the dataset is a major factor that affects the exact disribution. As the size of the dataset increases, reads distribution get closer to the exact 15/75/10 distribution.  

//How much time did you spend writing the code for this part?

The time spent on actual coding for part 2 is about 3 - 4 hours. But it took me additional 4-5 hours in the beginning to properly set up a development environemnt and editor for coding due to the unfamiliarity with python and terminal. Once I had vs_code set up on my laptop, I was able to write code in a much more efficient way.  




Part 3 Process Data

$ python3 processdata.py "reference_1.txt" "reads_1.txt" "alignfile_1.txt"
Reference Length: 1000
Number Reads: 600
Read Length: 50
Align 0: 0.135
Align 1: 0.7716666666666666
Align 2: 0.09333333333333334
elapsed time: 0.006119966506958008

$ python3 processdata.py "reference_2.txt" "reads_2.txt" "alignfile_2.txt"
Reference Length: 10000
Number Reads: 6000
Read Length: 50
Align 0: 0.14816666666666667
Align 1: 0.7498333333333334
Align 2: 0.102
elapsed time: 0.26036620140075684

$ python3 processdata.py "reference_3.txt" "reads_3.txt" "alignfile_3.txt"
Reference Length: 100000
Number Reads: 60000
Read Length: 50
Align 0: 0.15026666666666666
Align 1: 0.7495333333333334
Align 2: 0.1002
elapsed time: 23.887802362442017

// Does the reads distribution exactly match the distributions you computed in part 2?

Yes. Because the processdata program takes exact reads and reference from files generated from generatedata program. 
There is no random function involved in match readings with reference. 

//Discuss scalability of your implementation. Estimate the time to align 
the data for a human at 30x coverage and a read length of 50. Is it feasible to analyze all data
for a human using your program?

We can use the provided data and measured elapsed time to predict how long it takes to process alignments for human gene.By creating a trendline between reference lengths and process times on googlesheet. We can find the estimated relationship is:

Elapsed time ~= 0.114*e^(5.11E-5*ref_length)

When ref_length is 3,000,000, which is just 10% of human reference length, the estimated time to align is 4.3e^65s or 1.37e^58 years. From the data extrapolation, we can conclude that this program is not scalable for analzying human data.



//How much did you spend writing this part?
2-3 hours. Once the part 2 is done, part 3 is fairly easy.