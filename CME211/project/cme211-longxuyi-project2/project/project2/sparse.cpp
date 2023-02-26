#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "sparse.hpp"
#include "CGSolver.hpp"
#include "COO2CSR.hpp"
#include "matvecops.hpp"

    /* Method to modify sparse matrix dimensions */
    void SparseMatrix::Resize(int nrows, int ncols) {
        this->nrows = nrows;
        this->ncols = ncols;

    }  


    /* Method to add entry to matrix in COO format */
    void SparseMatrix::AddEntry(int i, int j, double val) {
        this->i_idx.push_back(i);
        this->j_idx.push_back(j);
        this->a.push_back(val);
    }


    /* Method to convert COO matrix to CSR format using provided function */
    void SparseMatrix::ConvertToCSR() {
        COO2CSR(a,i_idx,j_idx);
    }


    /* Method to perform sparse matrix vector multiplication using CSR formatted matrix */
    std::vector<double> SparseMatrix::MulVec(std::vector<double> &vec) {
        std::vector<double> result = mat_vec_mult(a, i_idx, j_idx, vec);

        return result;
    }


    /* Solve matrix equation with provided matrix A */
    int SparseMatrix::SolveCG(std::vector<double> &b,
             std::vector<double> &x,
             double              tol,
             std::string soln_prefix) {

        ConvertToCSR();

        int niter = CGSolver(this->a,this->i_idx,this->j_idx,b,x,tol,soln_prefix);

        return niter;      
    }