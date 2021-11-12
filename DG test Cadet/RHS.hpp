#pragma once
#include "DGspecific.hpp"



void volumeIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer);
void InterfaceFlux(Container& cache, Discretization& DG, ParameterProvider& para);
void InterfaceFluxAuxiliary(Container& cache, Discretization& DG, ParameterProvider& para);
void surfaceIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer, bool aux);
void ConvDisp(Container& cache, Discretization& DG, ParameterProvider& para, double t);
void applyMapping(Discretization DG, VectorXd& state);
void applyMapping_Aux(Discretization DG, ParameterProvider para, VectorXd& state);
void calcH(Container& cache, ParameterProvider& para);

Eigen::MatrixXd calcDQDC(Container& cache, ParameterProvider& para, int cell, int node);
void calcDC(Container& cache, ParameterProvider& para);