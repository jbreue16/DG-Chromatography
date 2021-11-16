#pragma once
#include "DGspecific.hpp"
#include "RHS.hpp"
#include<Eigen/Dense>

double timestep(Discretization DG, ParameterProvider para, double CFL);
Eigen::MatrixXd solveRK3(Container& cache, Discretization& DG, ParameterProvider& para, double tStart, double tEnd, double CFL, VectorXd start,
						 int t_save = 1);
Eigen::MatrixXd solveRK4(Container& cache, Discretization& DG, ParameterProvider& para, double tStart, double tEnd, double CFL, VectorXd start,
						 int t_save = 1);
Eigen::VectorXd solveEuler(Container& cache, Discretization& DG, ParameterProvider& para, double tStart, double tEnd, double CFL, VectorXd start,
						   int t_save = 1);
