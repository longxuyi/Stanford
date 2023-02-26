#include "CGSolver.hpp"
#include "matvecops.hpp"
#include <vector> 
#include <iostream>
#include <cmath>

int CGSolver(std::vector<double> &val,
             std::vector<int>    &row_ptr,
             std::vector<int>    &col_idx,
             std::vector<double> &b,
             std::vector<double> &x,
             double              tol) {


//solve initial residual vector r0 and p0
std::vector<double> Au0 = mat_vec_mult(val, row_ptr, col_idx, x);
std::vector<double> r0 = vect_sub(b, Au0);
std::vector<double> p0 = r0;

double L2normr0 = L2norm(r0);
double L2normr = L2norm(r0);
double normsq_cur = 0;
double normsq_new = 0;

int niter =0;
int nitermax = 50;
double alpha =0;
double beta = 0;

while (niter < nitermax) {
    niter ++;
    //solve alpha 
    std::vector<double> Ap = mat_vec_mult(val, row_ptr, col_idx, p0);
    double pAp = vect_mult(p0, Ap);
    normsq_cur = pow(L2norm(r0),2);
    alpha  = normsq_cur/pAp;

    //update solution vector x
    std::vector<double> alpha_p;
    alpha_p = scal_vec_mult(p0, alpha);
    x = vect_add(x, alpha_p);

    //update reisdual vector r0
    r0 = vect_sub(r0, scal_vec_mult(Ap, alpha));

    //update norm of r0
    L2normr = L2norm(r0);

    normsq_new = pow(L2norm(r0),2);

    if ((L2normr/L2normr0) < tol ) {
        std::cout<< "SUCCESS: CG Solver converged in " << niter << " iterations." << std::endl;
        break;
    }
    //update beta 
    normsq_new = pow(L2norm(r0),2);

    beta = normsq_new/normsq_cur;
    //update p0
    p0 = vect_add(r0,scal_vec_mult(p0, beta));
}
return niter;
}
