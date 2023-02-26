import sys
import os
import numpy as np
from enum import Enum

if len(sys.argv) < 3:
    print('Usage:')
    print('  python3 checksoln.py <maze_file> <solution_file>')
    sys.exit(0)

maze_file = sys.argv[1]
solution_file = sys.argv[2]

mf = open(maze_file, "r")
first_line = mf.readline()
#read first line of the maze file to get row and col values
r, c = int(first_line.split()[0]),int(first_line.split()[1])

#create a static array
maze = np.zeros([201,201])
for i in range(201):
    for j in range(201):
        if i >= r or j>= c:
             maze[i][j]==1
  
lines = mf.readlines()
for line in lines:
    i,j = int(line.split()[0]),int(line.split()[1])
    maze[i][j] = 1
mf.close()

#read solution file to check if solution is valid
sf = open(solution_file, "r")

#check if the maze enters from the first row
first_line = sf.readline()
entry_row, entry_col = int(first_line.split()[0]),int(first_line.split()[1])
if entry_row !=0:
    raise RuntimeError("Invalid Entry. Must enter from the first row")

lines = sf.readlines()
for line in lines:
    i,j = int(line.split()[0]),int(line.split()[1])
    #raise error when the solution path go through the wall or out of bounds 
    if maze[i][j] != 0 or i >= r or j >= c:
        raise RuntimeError("Invalid solution")

print("valid solution!")
sf.close()
