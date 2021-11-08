#pragma once
#include "DGspecific.hpp"



void volumeIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer);
void calcSurfaceFlux(Container& cache, Discretization& DG, ParameterProvider& para);
void calcSurfaceFluxAuxiliary(Container& cache, Discretization& DG, ParameterProvider& para);
void surfaceIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer, bool aux);
void ConvDisp(Container& cache, Discretization& DG, ParameterProvider& para, double t);
void applyJacobian(Discretization DG, VectorXd& state);

Eigen::MatrixXd calcDQDC(Container& cache, ParameterProvider& para, int cell, int node);
void calcDC(Container& cache, ParameterProvider& para);