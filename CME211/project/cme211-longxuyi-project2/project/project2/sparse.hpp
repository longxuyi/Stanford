#ifndef SPARSE_HPP
#define SPARSE_HPP


#include <string>
#include <vector>


class SparseMatrix
{
  private:
    std::vector<int> i_idx;
    std::vector<int> j_idx;
    std::vector<double> a;
    int ncols;
    int nrows;

    /* TODO: Add any additional private data attributes and/or methods  you need */


  public:
    /* Method to modify sparse matrix dimensions */
    void Resize(int nrows, int ncols);

    /* Method to add entry to matrix in COO format */
    void AddEntry(int i, int j, double val);

    /* Method to convert COO matrix to CSR format using provided function */
    void ConvertToCSR();

    /* Method to perform sparse matrix vector multiplication using CSR formatted matrix */
    std::vector<double> MulVec(std::vector<double> &vec);

    /* TODO: Add any additional public methods you need */

    //add a SolveCG method to solve linear system using matrix A in CSR format

    int SolveCG(std::vector<double> &b,
             std::vector<double> &x,
             double              tol,
             std::string soln_prefix);

};

#endif /* SPARSE_HPP */
