#pragma once
#include "DGspecific.hpp"
#include "RHS.hpp"
#include<Eigen/Dense>

void calcBoundary(Container& cache, Discretization DG, double t);
double timestep(Discretization DG, ParameterProvider para, double CFL, double tend);
Eigen::VectorXd solve(Container& cache, Discretization& DG, ParameterProvider& para, double tstart, double tend, double CFL, VectorXd start);