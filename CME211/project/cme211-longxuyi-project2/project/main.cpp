#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "COO2CSR.hpp"
#include "matvecops.hpp"
#include "CGSolver.hpp"

/* main function: read a matrix file in COO format and write a solution file */
int main(int argc, char *argv[]) {

if (argc < 2) {
    std::cout << "./main <input matrix file name> <output solution file name>\n" << std::endl;
    return 0;
}

//read the matrix file in COO format
std::ifstream COOfile(argv[1]);
std::vector<int> row_idx;
std::vector<int> col_idx;
std::vector<double> a;

if (COOfile.is_open()) {
    int nRow, nCol;
    COOfile >> nRow >> nCol;
    int nif, njf;
    double val;
    while (COOfile >> nif >> njf >> val){
        row_idx.push_back(nif);
        col_idx.push_back(njf);
        a.push_back(val); 
    }
 }
    COOfile.close();

//Convert the COO to CSR using provided COOSCSR() function
COO2CSR(a, row_idx, col_idx);

//initialize a solution vector u0, right hand side b, and tolerance
std::vector<double> x(row_idx.size()-1, 1.0);
std::vector<double> b(row_idx.size()-1, 0.0);
double tol = 1.e-5;

//solve the system using CG solver
int niter = CGSolver(a, row_idx, col_idx, b, x, tol); 

//write the solution file
std::ofstream soln_file(argv[2]);
for(int i = 0;  i < x.size(); i++) {
    if (soln_file.is_open()) {soln_file << x[i] << std::endl; }   
}
soln_file.close();

return 0;
}