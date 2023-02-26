#include "matvecops.hpp"
#include <vector>
#include <iostream>
#include <cmath>

//calculate norm of a vector
double L2norm(std::vector<double>  r0) {
    double norm_square = 0;
    double norm = 0;
    for (int i =0; i < r0.size(); i++) {
        norm_square = norm_square + pow(r0[i], 2);
    }
    norm = sqrt(norm_square);
    return norm;
}

//vector addition
std::vector<double> vect_add(std::vector<double> v1, std::vector<double> v2) {
  std::vector<double> result;
  for (int i =0; i < v1.size(); i++) {
    double new_val = v1[i] + v2[i];
    result.push_back(new_val);
  }
  return result;
}

//vector substraction
std::vector<double> vect_sub(std::vector<double> v1, std::vector<double> v2) {
  std::vector<double> result;
  for (int i =0; i < v1.size(); i++) {
    double new_val = v1[i] - v2[i];
    result.push_back(new_val);
  }
  return result;
}

//scalar vector multiplication
std::vector<double> scal_vec_mult(std::vector<double> vec, double num) {
  std::vector<double> result;
  for (int i =0; i < vec.size(); i++) {
  double new_val = vec[i]*num;
  result.push_back(new_val);
  }
  return result; 
}

//vector multiplication
double vect_mult(std::vector<double> v1, std::vector<double> v2) {
  double result =0;
  for (int i =0; i < v1.size(); i++) {
  double new_val = v1[i]*v2[i];
  result +=new_val;
  }
  return result; 
}

// CSR matrix and vector multiplication
std::vector<double> mat_vec_mult(std::vector<double> &val, std::vector<int> &row_ptr,
                              std::vector<int> &col_idx, std::vector<double> &u) {
  
std::vector<double> result(row_ptr.size()-1,0);
  for (int i = 0; i < row_ptr.size()-1; i++) {
    for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) {
      result[i] += val[j]* u[col_idx[j]];
    }
  }
  return result;
}


