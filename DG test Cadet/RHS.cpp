
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
				// strong volume integral -> D h
				for (int ii = 0; ii < DG.nNodes; ii++) {
					stateDer[Cell * para.strideCell() + Node * para.strideNode() + Comp * para.strideComp()]
					+= DG.polyDerM(Node, ii) * state(Cell * para.strideCell() + ii * para.strideNode() + Comp * para.strideComp());
				}
			}
		}
	}
}
/**
* @brief calculates and fills the surface flux values for auxiliary equation
* @param [in] aux true if auxiliary, else main equation
*/
void calcSurfaceFluxAuxiliary(Container& cache, Discretization& DG, ParameterProvider& para) {

	// calculate inner interface fluxes
	for (int Cell = 1; Cell < para.nCells; Cell++) {
		for (int Comp = 0; Comp < para.nComp; Comp++) {
			cache.surfaceFlux[Cell * para.strideNode() + Comp * para.strideComp()] // inner interfaces
				= centralFlux(cache.c[Cell * para.strideCell() + Comp * para.strideComp()], // left cells
					cache.c[Cell * para.strideCell() + Comp * para.strideComp() - para.strideNode()], // right cells
					auxiliaryFlux, para);
		}
	}
	 // calculate boundary interface fluxes
	for (int Comp = 0; Comp < para.nComp; Comp++) {
		cache.surfaceFlux[Comp * para.strideComp()] // left boundary interface
			= centralFlux(cache.boundary[Comp * para.strideComp()], // boundary value
				cache.c[Comp * para.strideComp()], // first cell first node
				auxiliaryFlux, para);
		cache.surfaceFlux[(para.nCells) * para.strideNode() + Comp * para.strideComp()] // right boundary interface
			= centralFlux(cache.c[(para.nCells - 1) * para.strideCell() + (DG.nNodes - 1) * para.strideNode() + Comp * para.strideComp()], // last cell last node
				cache.boundary[para.strideNode() + Comp * para.strideComp()], // boundary value
				auxiliaryFlux, para);
	}
}
/**
* @brief calculates and fills the surface flux values for main equation
* @param [in] aux true if auxiliary, else main equation
*/
void calcSurfaceFlux(Container& cache, Discretization& DG, ParameterProvider& para) {

	// calculate inner interface fluxes
	for (int Cell = 1; Cell < para.nCells; Cell++) {
		for (int Comp = 0; Comp < para.nComp; Comp++) {
			// h* = h*_conv + h*_disp
			cache.surfaceFlux[Cell * para.strideNode() + Comp * para.strideComp()] // inner interfaces
				= - DG.numFlux(cache.c[Cell * para.strideCell() + Comp * para.strideComp()], // left cells
					cache.c[Cell * para.strideCell() + Comp * para.strideComp() - para.strideNode()], // right cells
					convectionFlux, para) // convection part
				+ centralFlux(cache.S[Cell * para.strideCell() + Comp * para.strideComp()], // left cells
					cache.S[Cell * para.strideCell() + Comp * para.strideComp() - para.strideNode()], // right cells
					dispersionFlux, para); // dispersion part
		}
	}
	//calculate boundary interface fluxes
	for (int Comp = 0; Comp < para.nComp; Comp++) {
		// h* = h*_conv + h*_disp
		// left boundary interface
		cache.surfaceFlux[Comp * para.strideComp()]
			= - DG.numFlux(cache.boundary[Comp * para.strideComp()], // boundary value c
						 cache.c[Comp * para.strideComp()], // first cell first node
						 convectionFlux, para) // convection part
			+ centralFlux(cache.boundary[2 * para.nComp + Comp * para.strideComp()], // boundary value S
						  cache.S[Comp * para.strideComp()], // first cell first node
						  dispersionFlux, para); // dispersion part
		 // right boundary interface
		cache.surfaceFlux[(para.nCells) * para.strideNode() + Comp * para.strideComp()]
			= - DG.numFlux(cache.c[(para.nCells - 1) * para.strideCell() + (DG.nNodes - 1) * para.strideNode() + Comp * para.strideComp()], // last cell last node
						  cache.boundary[para.strideNode() + Comp * para.strideComp()], // boundary value c
						  convectionFlux, para) // convection part
			+ centralFlux(cache.S[(para.nCells - 1) * para.strideCell() + (DG.nNodes - 1) * para.strideNode() + Comp * para.strideComp()], // last cell last node
						  cache.boundary[(2* para.nComp-1) + para.strideNode() + Comp * para.strideComp()], // boundary value S
						  dispersionFlux, para); // dispersion part
	}
}
/**
* @brief calculates the surface Integral
* @param [in] State vector to be changed
* @param [in] aux true if auxiliary, false for main equation
*/
void surfaceIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer, bool aux) {
	// calc numerical flux values c* or h* depending on equation switch aux
	(aux == 1) ? calcSurfaceFluxAuxiliary(cache, DG, para) : calcSurfaceFlux(cache, DG, para);
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
	cache.h = sqrt(para.dispersion) * cache.S - para.velocity * cache.c;
}

/**
* @brief returns (F* dq/dc + I)^-1
* @param [in] cell index of the cell
*/
MatrixXd calcDQDC(Container& cache, ParameterProvider& para, int cell, int node) {
	//TODO calculation depending on model isotherm
	MatrixXd Jac = MatrixXd::Zero(para.nComp, para.nComp);
	if (para.isotherm == "Langmuir") {
		// first calc (1 + sum(b_i * c_i))^-1
		double factor = 1.0;
		for (int comp = 0; comp < para.nComp; comp++) {
			factor += para.ADratio[comp] * cache.c[cell * para.strideCell() + node * para.strideNode() + comp];
		}
		factor = 1.0 / factor;
		// estimate Jacobian dq/dc + I
		for (int i = 0; i < para.nComp; i++) {
			for (int j = 0; j < para.nComp; j++) {
				Jac(i, j) = para.porosity * (- (pow(factor, 2)) * para.adsorption[i] * para.ADratio[j] * cache.c[cell * para.strideCell() + node * para.strideNode() + i * para.strideComp()]);
				if (i == j) {
					Jac(i, j) += para.porosity * (para.adsorption[i] * factor) + 1; // +1 for F*Jac + I
				}
			}
		}
	}
	else {
		throw std::invalid_argument("spelling error or this isotherm is not implemented yet");
	}
	return Jac.inverse();
}

/**
* @brief calculates c from w
*/
void calcDC(Container& cache, ParameterProvider& para) {
	if (para.isotherm == "Linear") { // diagonal Jacobian dq/dc
		for (int comp = 0; comp < para.nComp; comp++) {
			double factor = 1.0 / (1.0 + para.porosity * para.adsorption[comp]);
			for (int nodes = 0; nodes < para.nCells * (para.polyDeg + 1); nodes++) {
				cache.dc[nodes * para.strideNode() + comp] = factor * cache.w[nodes * para.strideNode() + comp];
			}
		}
	}
	else { // Dense Jacobian dq/dc
		for (int cell = 0; cell < para.nCells; cell++) {
			for (int node = 0; node < para.polyDeg + 1; node++) {
				// get (I + F* J)^-1 for each node, J=dq/dc being the isotherm Jacobian
				MatrixXd Jac = calcDQDC(cache, para, cell, node);
				// c_components =  (I + F* J)^-1 * w_components
				for (int comp = 0; comp < para.nComp; comp++) {
					for (int comp2 = 0; comp2 < para.nComp; comp2++){
					cache.dc[cell * para.strideCell() + node * para.strideNode() + comp * para.strideComp()]
						+= Jac(comp, comp2) * cache.w[cell * para.strideCell() + node * para.strideNode() + comp2 * para.strideComp()];
					}
				}
			}
		}
	}
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
void ConvDisp(Container& cache, Discretization& DG, ParameterProvider& para, double t) {
	// reset w, h, auxiliary S, dc (rhs), surface flux
	cache.w.setZero();
	cache.h.setZero();
	cache.S.setZero();
	cache.dc.setZero();
	cache.surfaceFlux.setZero();
	// boundary values are calculated in time Integration

	// first solve the auxiliary system for dispersion operator
	DG.BoundCond(t, cache, DG.BoundFunc, para); // update boundary values for c 
	volumeIntegral(cache, DG, para, cache.c, cache.S);
	surfaceIntegral(cache, DG, para, cache.c, cache.S, 1);
	applyJacobian(DG, cache.S);
	//std::cout << "auxiliary S " << cache.S << "auxiliary S " << std::endl;
	cache.surfaceFlux.setZero(); // surface flux storage is used twice

	// solve dispersion convection part of main equation
	DG.BoundCond(t, cache, DG.BoundFunc, para); // update boundary values for S
	calcH(cache, para); // calculate the substitute h(S(c), c) = sqrt(D_ax) S(c) - v c
	//std::cout << "h(S,c)= " << cache.h << "h(S,c) " << std::endl;//bei FSP t=0 korrekt
	volumeIntegral(cache, DG, para, cache.h, cache.w);
	surfaceIntegral(cache, DG, para, cache.h, cache.w, 0);
	//std::cout << "w(c)= " << cache.w << "w(c) " << std::endl;
	//std::cout << "surfFlux= " << cache.surfaceFlux << "surfFlux" << std::endl;
	applyJacobian(DG, cache.w);

	// estimate state vector derivative dc of mobile phase from substitues w, h
	// TODO check calcDC !
	calcDC(cache, para);
	std::cout << "dc= " << cache.dc << "dc " << std::endl;
}
