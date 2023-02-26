#ifndef HEAT_HPP
#define HEAT_HPP

#include <string>
#include <vector>

#include "sparse.hpp"

class HeatEquation2D
{
  private:
    SparseMatrix A;
    std::vector<double> b, x;

    /* TODO: Add any additional private data attributes and/or methods you need */

  public:
    /* Method to setup Ax=b system */
    int Setup(std::string inputfile);

    /* Method to solve system using CGsolver */
    int Solve(std::string soln_prefix);

    /* TODO: Add any additional public methods you need */

};

#endif /* HEAT_HPP */
