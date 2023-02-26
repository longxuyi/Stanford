In this assignment, a C++ program is created to compute the solution path for a given maze using the right hand wall following algorithm, and store the solution in a file. A python program is also created to test the solution.

For the C++ code, we first define a static array based on the largest maze provided. Then read the maze file to create walls.
Then we use a switch function to explore the maze and keep track of the direction as we go. Each step is directly written into the solution file until we reach the exit at the last row of the maze.

For python program, it reads the maze file and solution file to check if the solution is valid. For error handling, runtime error will be raised in not entering from the first row, hitting the wall, or moving out of bounds.