
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
using namespace std;

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
				// strong volume integral -> [- D* state]
				for (int ii = 0; ii < DG.nNodes; ii++) {
					stateDer[Cell * para.strideCell() + Node * para.strideNode() + Comp * para.strideComp()]
					-= DG.polyDerM(Node, ii) * state(Cell * para.strideCell() + ii * para.strideNode() + Comp * para.strideComp());
				}
			}
		}
	}
}
/**
* @brief calculates and fills the surface flux values for auxiliary equation
* @param [in] aux true if auxiliary, else main equation
*/
void InterfaceFluxAuxiliary(Container& cache, Discretization& DG, ParameterProvider& para) {
	// [c* = 0.5 (c_l + c_r)]
	// calculate inner interface fluxes
	for (int Cell = 1; Cell < para.nCells; Cell++) {
		for (int Comp = 0; Comp < para.nComp; Comp++) {
			cache.surfaceFlux[Cell * para.strideNode() + Comp * para.strideComp()] // inner interfaces
				= centralFlux(cache.c[Cell * para.strideCell() + Comp * para.strideComp() - para.strideNode()], // left node
							  cache.c[Cell * para.strideCell() + Comp * para.strideComp()], // right node
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
* @detail calculates [h* = h*_conv + h*_disp = DG.riemannSolver(v c_l, v c_r) + 0.5 sqrt(D_ax) (S_l + S_r)] and stores it into DG.numFlux
* @param [in] aux true if auxiliary, else main equation
*/
void InterfaceFlux(Container& cache, Discretization& DG, ParameterProvider& para) {
	// [h* = h*_conv + h*_disp = numFlux(v c_l, v c_r) + 0.5 sqrt(D_ax) (S_l + S_r)]
	// calculate inner interface fluxes
	for (int Cell = 1; Cell < para.nCells; Cell++) {
		for (int Comp = 0; Comp < para.nComp; Comp++) {
			// h* = h*_conv + h*_disp
			cache.surfaceFlux[Cell * para.strideNode() + Comp * para.strideComp()] // inner interfaces
				= DG.numFlux(cache.c[Cell * para.strideCell() + Comp * para.strideComp() - para.strideNode()], // left cell
							 cache.c[Cell * para.strideCell() + Comp * para.strideComp()], // right cell 
							 convectionFlux, para) // convection part
				- centralFlux(cache.S[Cell * para.strideCell() + Comp * para.strideComp() - para.strideNode()], // left cell
							  cache.S[Cell * para.strideCell() + Comp * para.strideComp()], // right cell
							  dispersionFlux, para); // dispersion part
		}
	}
	//calculate boundary interface fluxes
	for (int Comp = 0; Comp < para.nComp; Comp++) {
		// h* = h*_conv + h*_disp
		// left boundary interface
		cache.surfaceFlux[Comp * para.strideComp()]
			= DG.numFlux(cache.boundary[Comp * para.strideComp()], // left boundary value c
						 cache.c[Comp * para.strideComp()], // first cell first node
						 convectionFlux, para) // convection part
			- centralFlux(cache.boundary[2 * para.nComp + Comp * para.strideComp()], // left boundary value S
						  cache.S[Comp * para.strideComp()], // first cell first node
						  dispersionFlux, para); // dispersion part
		 // right boundary interface
		cache.surfaceFlux[(para.nCells) * para.strideNode() + Comp * para.strideComp()]
			= DG.numFlux(cache.c[(para.nCells - 1) * para.strideCell() + (DG.nNodes - 1) * para.strideNode() + Comp * para.strideComp()], // last cell last node
						 cache.boundary[para.nComp + Comp * para.strideComp()], // right boundary value c
						 convectionFlux, para) // convection part
			- centralFlux(cache.S[(para.nCells - 1) * para.strideCell() + (DG.nNodes - 1) * para.strideNode() + Comp * para.strideComp()], // last cell last node
						  cache.boundary[3 * para.nComp + Comp * para.strideComp()], // right boundary value S
						  dispersionFlux, para); // dispersion part
	}
}
/**
* @brief calculates the surface Integral
* @param [in] state relevant state vector
* @param [in] stateDer state derivative vector the solution is added to
* @param [in] aux true for auxiliary, false for main equation
*/
void surfaceIntegral(Container& cache, Discretization& DG, ParameterProvider& para, VectorXd& state, VectorXd& stateDer, bool aux) {
	// calc numerical flux values c* or h* depending on equation switch aux
	(aux == 1) ? InterfaceFluxAuxiliary(cache, DG, para) : InterfaceFlux(cache, DG, para);
	// calc surface integral
	for (int Cell = 0; Cell < para.nCells; Cell++) {
		for (int Comp = 0; Comp < para.nComp; Comp++) {
			// strong surface integral -> M^-1 B [state - state*]
			stateDer[Cell * para.strideCell() + Comp * para.strideComp()] // first node
				-= DG.invWeights[0] * (state[Cell * para.strideCell() + Comp * para.strideComp()] // first node
					- cache.surfaceFlux(Cell * para.nComp + Comp));
			stateDer[Cell * para.strideCell() + Comp * para.strideComp() + para.polyDeg * para.strideNode()] // last node
				+= DG.invWeights[para.polyDeg] * (state[Cell * para.strideCell() + Comp * para.strideComp() + para.polyDeg * para.strideNode()]
					- cache.surfaceFlux((Cell + 1) * para.nComp + Comp));
		}
	}
}

/**
* @brief calculates the substitute h = vc - sqrt(D_ax) S(c)
*/
void calcH(Container& cache, ParameterProvider& para) {
	cache.h = para.velocity * cache.c - sqrt(para.dispersion) * cache.S;
}

/**
* @brief returns (F* dq/dc + I)^-1 for a specified node
* @param [in] cell index
* @param [in] node index
* @param [out] inverse Jacobian isotherm
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
* @brief calculates dc from w
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
void applyMapping(Discretization DG, VectorXd& state) {
	state = state * (2 / DG.deltaX);
}
/**
* @brief applies the inverse Jacobian of the mapping and Auxiliary factor
*/
void applyMapping_Aux(Discretization DG, ParameterProvider para, VectorXd& state) {
	state = state * (-2 / DG.deltaX) * ((para.dispersion==0) ? 1.0:sqrt(para.dispersion));
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
	//std::cout << "boundaryS= " << cache.boundary << "boundaryS" << std::endl;
	volumeIntegral(cache, DG, para, cache.c, cache.S);
	surfaceIntegral(cache, DG, para, cache.c, cache.S, 1);
	applyMapping_Aux(DG, para, cache.S);
	//std::cout << "auxiliary S " << cache.S << "auxiliary S " << std::endl;
	cache.surfaceFlux.setZero(); // surface flux storage is used twice

	// solve dispersion convection part of main equation
	DG.BoundCond(t, cache, DG.BoundFunc, para); // update boundary values for S
	//std::cout << "boundary= " << cache.boundary << "boundary" << std::endl;
	calcH(cache, para); // calculate the substitute h(S(c), c) = sqrt(D_ax) S(c) - v c
	//std::cout << "h(S,c)= " << cache.h << "h(S,c) " << std::endl;//bei FSP t=0 korrekt
	volumeIntegral(cache, DG, para, cache.h, cache.w);
	//std::cout << "(c)= " << cache.c << "(c) " << std::endl;
	surfaceIntegral(cache, DG, para, cache.h, cache.w, 0);
	//std::cout << "surfFlux= " << cache.surfaceFlux << "surfFlux" << std::endl;
	//std::cout << "w(c)= " << cache.w << "w(c) " << std::endl;
	applyMapping(DG, cache.w);

	// estimate state vector derivative dc of mobile phase from substitues w, h
	calcDC(cache, para);
	//std::cout << "dc= " << cache.dc << "dc " << std::endl;
}
