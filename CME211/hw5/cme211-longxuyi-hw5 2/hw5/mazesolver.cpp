#include <iostream>
#include <fstream>
#include <string>

#define ni 201
#define nj 201

int main(int argc, char *argv[]) {

    //create a empty 2d static array
    int arr[ni][nj] = {{0}};
    std::ifstream mazefile(argv[1]);

    int nRow, nCol;
    if (mazefile.is_open()) {
        
        //read row and col of the maze
        mazefile >> nRow >> nCol;

        //read maze file to add walls
        int nif, njf;
        while (mazefile >> nif >> njf){
            arr[nif][njf] = 1;
        }
    }
    mazefile.close();
    
    int x = 0;
    int y = 0;
    // find the entry position
    std::ofstream f(argv[2]);
    for(int j = 0;  j < nCol; j++) {
        if (arr[0][j] == 0) {
            if (f.is_open()) {
                x =0; y =j;
                f << x << " " << j << std::endl; 
            }
        }
    }

    enum direction {N,S,E,W};
    direction d = S;

    //continue to explore until the solution path reach the final row
    while (x < nRow - 1) {
        switch (d){
            case S:
                //if (left is open) {y = y-1, direction is now west;}
                if (arr[x][y-1]==0) {y = y-1; d = W;}
                // if down is open
                else if (arr[x+1][y]==0) {x= x+1; d = S;}
                // if right is open
                else if (arr[x][y+1]==0) {y = y+1; d = E;}
                // else turn around
                else {x =x-1; d = N;}
                break;
            case N:
                //if (right is open) {y = y+1, direction is now west;}
                if (arr[x][y+1]==0) {y = y+1; d = E;}
                // if up is open
                else if (arr[x-1][y]==0) {x= x-1; d = N;}
                // if left is open
                else if (arr[x][y-1]==0) {y = y-1; d = W;}
                // else turn around
                else {x = x+1; d = S;}
                break;
            case E:
                //if (down is open) {x = x+1, direction is now west;}
                if (arr[x+1][y]==0) {x = x+1; d = S;}
                // if right is open
                else if (arr[x][y+1]==0) {y= y+1; d = E;}
                // if up is open
                else if (arr[x-1][y]==0) {x = x-1; d = N;}
                // else turn around
                else {y = y-1; d = W;}
                break;
            case W:
                //if (up is open) {x = x+1, direction is now west;}
                if (arr[x-1][y]==0) {x = x-1; d = N;}
                // if left is open
                else if (arr[x][y-1]==0) {y= y-1; d = W;}
                // if down is open
                else if (arr[x+1][y]==0) {x = x+1; d = S;}
                // else turn around
                else {y = y+1; d = E;}
                break;
        }
        //write move coordinates into the solution file
        f << x << " " << y << std::endl; 

     }
    f.close();
    return 0;
}