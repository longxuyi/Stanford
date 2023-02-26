#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "CGSolver.hpp"
#include "COO2CSR.hpp"
#include "sparse.hpp"
#include "heat.hpp"


/* Method to setup Ax=b system */
int HeatEquation2D::Setup(std::string inputfile) {
    float l, w, h;
    int Tc, Th;
    //read input file
    std::ifstream input(inputfile);
    if (input.is_open()) {
        input >> l >> w >> h;
        input >> Tc >> Th;
    }
    input.close();

    //determine number of unknown nodes along x and y direction
    int n_x =(int) (l/h);
    int n_y =(int) (w/h -1);

    int num_unknowns = n_x * n_y;
    
    //define the size of matrix A 
    this->A.Resize(num_unknowns, num_unknowns);

    //initialize vector x
    this->x.resize(num_unknowns, 1.0);
    
    //iterate through each row of stencil
    for (int j = 0; j < n_y; j++) {
        //iterate through the column of stencil
        for (int i = 0; i < n_x; i++) {

            int center = i + j*n_x ;
            int up = (i + n_x) + j*n_x;
            int down = (i - n_x) + j*n_x;

            // special case for periodic boundary: center i = 0; left = n_x - 1;
            int left = (i + n_x - 1)%n_x + j*n_x; 

            // special case for periodic boundary: center i = n_x - 1; right = 0;
            int right = (i + n_x + 1)%n_x + j*n_x;

            this->A.AddEntry(center, center, 4.0);
            this->A.AddEntry(center, left, -1.0);
            this->A.AddEntry(center, right, -1.0);

            //if there is no Tc and Th points, set right hand side b element as 0
            if (j != 0 && j != n_y-1) this->b.push_back(0.0);

            //when down node is on Tc boundary, solve for temperature at lower boundary.
            double T_cold = -Tc*(exp(-10*(i*h - l/2)*(i*h - l/2)) - 2); 

            if (j == 0) this->b.push_back(T_cold);
            else this->A.AddEntry(center, down, -1.0);

            //when up node is on Th boundary
            if (j == n_y -1) this->b.push_back(Th);
            else this->A.AddEntry(center, up, -1.0);
        }
    }

    return 0;
}


/* Method to solve system using CGsolver */
int HeatEquation2D::Solve(std::string soln_prefix) {

    //Solve linear system
    double tol = 1.e-5;
    this->A.SolveCG(this->b,this->x,tol, soln_prefix);

    return 0;
}