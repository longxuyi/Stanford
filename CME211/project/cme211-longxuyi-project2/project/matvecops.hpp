#ifndef MATVECOPS_HPP
#define MATVECOPS_HPP

#include <vector>

//calculate norm of a vector
double L2norm(std::vector<double> r0);

//vector multiplication
double vect_mult(std::vector<double> v1, std::vector<double> v2);

//vector addition
std::vector<double> vect_add(std::vector<double> v1, 
                          std::vector<double> v2);

//vector substraction
std::vector<double> vect_sub(std::vector<double> v1, 
                          std::vector<double> v2);

//scalar vector multiplication
std::vector<double> scal_vec_mult(std::vector<double> vec, double num);

//vector multiplication
double vect_mult(std::vector<double> v1, std::vector<double> v2);                        

// matrix and vector multiplication
std::vector<double> mat_vec_mult(std::vector<double> &val,
                              std::vector<int>    &row_ptr,
                              std::vector<int>    &col_idx,
                              std::vector<double>    &u);


#endif /* MATVECOPS_HPP */