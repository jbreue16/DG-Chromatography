#pragma once
#include "DGspecific.hpp"
#include "RHS.hpp"
#include<Eigen/Dense>

double timestep(Discretization DG, ParameterProvider para, double CFL, double tend);
Eigen::VectorXd solveRK3(Container& cache, Discretization& DG, ParameterProvider& para, double tstart, double tend, double CFL, VectorXd start);
Eigen::VectorXd solveRK4(Container& cache, Discretization& DG, ParameterProvider& para, double tstart, double tend, double CFL, VectorXd start);
Eigen::VectorXd solveEuler(Container& cache, Discretization& DG, ParameterProvider& para, double tstart, double tend, double CFL, VectorXd start);
