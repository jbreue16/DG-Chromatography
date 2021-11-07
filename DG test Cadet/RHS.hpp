#pragma once
#include "DGspecific.hpp"



void volumeIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer);
void calcSurfaceFlux(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, bool aux);
void surfaceIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer, bool aux);
void ConvDisp(Container& cache, Discretization& DG, ParameterProvider& para);
void applyJacobian(Discretization DG, VectorXd& state);

Eigen::MatrixXd calcDQDC(Container& cache, ParameterProvider& para, int cell, int node);
void calcC(Container& cache, ParameterProvider& para);