#pragma once
#include "DGspecific.hpp"



void volumeIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer);
void calcSurfaceFlux(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, bool aux);
void surfaceIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer, bool aux);
void ConvDisp(Container& cache, Discretization& DG, ParameterProvider& para);
void applyJacobian(Discretization DG, VectorXd& state);