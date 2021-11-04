
#define _USE_MATH_DEFINES

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>

#include "DGspecific.hpp"
#include<Eigen/Dense>
using namespace Eigen;


/**
* @brief calculates the volume Integral of the auxiliary equation
* @param [in] current state vector
* @param [in] stateDer vector to be changed
* @param [in] aux true if auxiliary, else main equation
*/
void volumeIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer) {
	// loop over all Cells, Nodes, Components
	for (int Cell = 0; Cell < para.nCells; Cell++) {
		for (int Node = 0; Node < DG.nNodes; Node++) {
			for (int Comp = 0; Comp < para.nComp; Comp++) {
				// strong volume integral -> D c
				for (int ii = 0; ii < DG.nNodes; ii++) {
					stateDer[Cell * para.strideCell() + Node * para.strideNode() + Comp * para.strideComp()]
					+= DG.polyDerM(Node, ii) * state(Cell * para.strideCell() + ii * para.strideNode() + Comp * para.strideComp());
				}
			}
		}
	}
}
/**
* @brief calculates and fills the surface flux values
* @param [in] aux true if auxiliary, else main equation
*/
void calcSurfaceFlux(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, bool aux) {
	// switch to choose between auxiliary equation and main equation
	Flux flux = (aux == 1) ? auxiliaryFlux : advectionDispersionFlux;
	riemannSolver numFlux = (aux == 1) ? centralFlux : DG.numFlux;

	// calculate inner interface fluxes
	for (int Cell = 1; Cell < para.nCells; Cell++) {
		for (int Comp = 0; Comp < para.nComp; Comp++) {
			cache.surfaceFlux[Cell * para.strideNode() + Comp * para.strideComp()] // inner interfaces
				= DG.numFlux(state[Cell * para.strideCell() + Comp * para.strideComp()], // left cells
					state[Cell * para.strideCell() + Comp * para.strideComp() - para.strideNode()], // right cells
					flux, para);
		}
	}
	// calculate boundary interface fluxes
	for (int Comp = 0; Comp < para.nComp; Comp++) {
		cache.surfaceFlux[Comp * para.strideComp()] // left boundary interface
			= DG.numFlux(cache.boundary[Comp * para.strideComp()], // boundary value
				state[Comp * para.strideComp()], // first cell first node
				flux, para);
		cache.surfaceFlux[(para.nCells) * para.strideNode() + Comp * para.strideComp()] // right boundary interface
			= DG.numFlux(state[(para.nCells - 1) * para.strideCell() + (DG.nNodes - 1) * para.strideNode() + Comp * para.strideComp()], // last cell last node
				cache.boundary[para.strideNode() + Comp * para.strideComp()], // boundary value
				flux, para);
	}
}
/**
* @brief calculates the surface Integral
* @param [in] State vector to be changed
* @param [in] aux true if auxiliary, false for main equation
*/
void surfaceIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer, bool aux) {
	// calc numerical flux values c*
	calcSurfaceFlux(cache, DG, para, state, aux);
	// calc surface integral
	for (int Cell = 0; Cell < para.nCells; Cell++) {
		for (int Comp = 0; Comp < para.nComp; Comp++) {
			// strong surface integral -> M^-1 B (c*-c)
			stateDer[Cell * para.strideCell() + Comp * para.strideComp()] // left interfaces
				+= - DG.invWeights[0] * (cache.surfaceFlux(Cell * para.nComp + Comp) 
				- state[Cell * para.strideCell() + Comp * para.strideComp()] );
			stateDer[Cell * para.strideCell() + Comp * para.strideComp() + para.polyDeg * para.strideNode()] // right interfaces
				+= DG.invWeights[para.polyDeg] * (cache.surfaceFlux((Cell + 1) * para.nComp + Comp)
				- state[Cell * para.strideCell() + Comp * para.strideComp() + para.polyDeg * para.strideNode()] );
		}
	}
}
/**
* @brief calculates the substitute h = D_ax S(c) - vc
*/
void calcH(Container& cache, ParameterProvider& para) {
	cache.h = para.dispersion * cache.S - para.velocity * cache.u;
}

/**
* @brief applies the inverse Jacobian of the mapping
*/
void applyJacobian(Discretization DG, VectorXd& state) {
	state = state * (2 / DG.deltaX);
}
/**
* @brief calculates the convection dispersion part of right hand side
*
*/
void ConvDisp(Container& cache, Discretization& DG, ParameterProvider& para) {
	// reset du, auxiliary q, and surface flux
	cache.du.setZero();
	cache.S.setZero();
	cache.surfaceFlux.setZero();
	// TODO: actualize boundary flux values !

	// first solve the auxiliary system for dispersion operator
	volumeIntegral(cache, DG, para, cache.u, cache.S);
	surfaceIntegral(cache, DG, para, cache.u, cache.S, 1);
	applyJacobian(DG, cache.S);
	cache.surfaceFlux.setZero(); // surface flux storage is used twice

	// solve dispersion convection part of main equation
	calcH(cache, para); // calculate the substitute h(S(c), c) = D_ax S(c) - v c
	volumeIntegral(cache, DG, para, cache.u, cache.du);
	surfaceIntegral(cache, DG, para, cache.u, cache.du, 0);
	applyJacobian(DG, cache.du);
}
