#pragma once

#include <iostream>

#include "DGspecific.hpp"
#include "RHS.hpp"
#include<Eigen/Dense>

/** 
*@file implements the third order explicit RK scheme as used by Kristian
*/


/**
*@brief returns timestep for linear advection
*/
double timestep(Discretization DG, ParameterProvider para, double CFL) {
	// Delta X_eff = (DG.deltaX/DG.nNodes)
	return CFL* (DG.deltaX/DG.nNodes) / abs(para.velocity);
}
/**
*@brief 3rd order expl. RK
* @param [in] t_save how many timesteps should be saved at outlet
* @param [in] outlet Matrix to store the outlet solution
*/
Eigen::MatrixXd solveRK3(Container& cache, Discretization& DG, ParameterProvider& para, double tStart, double tEnd, double CFL, VectorXd start,
	int t_save = 1) {
	double t = tStart;
	double Dt = 0.0;
	cache.c = start;
	int tCount = 0;

	VectorXd v1 = VectorXd::Zero(cache.c.size());
	VectorXd v2 = VectorXd::Zero(cache.c.size());

	while(t < tEnd) {
		// determine current timestep
		Dt = (timestep(DG, para, CFL) + t > tEnd) ? tEnd - t : timestep(DG, para, CFL);
		//std::cout << "Timestep: " << Dt << std::endl;
		// calc convection dispersion rhs
		ConvDisp(cache, DG, para, t); // stores rhs in cache.dc

		// expl. 3-stage RK
		v1 = cache.c + Dt * cache.dc;
		v2 = cache.c + Dt * (0.25 * cache.dc + 0.25 * v1);
		cache.c = cache.c + Dt * ((1.0 / 6.0) * cache.dc + (2.0 / 3.0) * v2 + (1.0 / 6.0) * v1);

		t += Dt;
		tCount++;
		std::cout << "t=" << t << std::endl;
	}
	std::cout << "timestep:" << tCount << std::endl;
	return Eigen::MatrixXd::Zero(1,1);
}

/**
*@brief Carpenter Kennedy 4th order low storage expl. RK
* @detail is freestream preserving (dispersion = 0.0)
* @param [in] t_save how many timesteps should be saved at outlet
* @param [in] outlet Matrix to store the outlet solution
*/
Eigen::MatrixXd solveRK4(Container& cache, Discretization& DG, ParameterProvider& para, double tStart, double tEnd, double CFL, VectorXd start,
	int t_save = 1) {
	double t = tStart;
	double Dt = 0.0;
	int tCount = 0;
	int tSaveCount = 0;
	VectorXd outlet = VectorXd::Zero(para.strideNode() * t_save);
	VectorXd tSave = VectorXd::LinSpaced(t_save, tStart, tEnd);

	// RK coefficients A, B, C for each stage
	VectorXd A(5);
	A << 0.0, -567301805773.0 / 1357537059087.0, -2404267990393.0 / 2016746695238.0,
		-3550918686646.0 / 2091501179385.0, -1275806237668.0 / 842570457699.0;
	VectorXd B(5);
	B << 0.0, 1432997174477.0 / 9575080441755.0, 2526269341429.0 / 6820363962896.0,
		2006345519317.0 / 3224310063776.0, 2802321613138.0 / 2924317926251.0;
	VectorXd C(5);
	C << 1432997174477.0 / 9575080441755.0, 5161836677717.0 / 13612068292357.0,
		1720146321549.0 / 2090206949498.0, 3134564353537.0 / 4481467310338.0, 2277821191437.0 / 14882151754819.0;
	VectorXd g = VectorXd::Zero(cache.c.size());

	cache.c = start;
	while (t < tEnd) {
		// determine current timestep
		Dt = (timestep(DG, para, CFL) + t > tEnd) ? tEnd - t : timestep(DG, para, CFL);
		//std::cout << "Timestep: " << Dt << std::endl;
		// calc convection dispersion rhs
		ConvDisp(cache, DG, para, t); // stores rhs in cache.dc
		g = cache.dc;
		// expl. 5-step RK
		for (int step = 0; step < 5; step++) {
			ConvDisp(cache, DG, para, t + B[step] * Dt);
			g = A[step] * g + cache.dc;
			cache.c += C[step] * Dt * g;
		}
		//save outlet solution
		if (t == tSave[tSaveCount] && t_save > 1) {
			for (int comp = 0; comp < para.nComp; comp++) {
				outlet[(t_save - tSaveCount - 1) * para.strideNode() + comp]
					= cache.c[para.nCells * para.strideCell() - para.nComp + comp];
			}
			tSaveCount++;
		}
		t += Dt;
		tCount++;
		std::cout << "t=" << t << std::endl;
	}
	std::cout << "timesteps: " << tCount << std::endl;
	return outlet;
}
/**
*@brief expl. Euler
* @param [in] t_save how many timesteps should be saved at outlet
* @param [in] outlet Matrix to store the outlet solution
*/
Eigen::VectorXd solveEuler(Container& cache, Discretization& DG, ParameterProvider& para, double tStart, double tEnd, double CFL, VectorXd start,
							int t_save = 1) {
	double t = tStart;
	double Dt = 0.0;
	int tCount = 0;
	int tSaveCount = 0;
	VectorXd outlet = VectorXd::Zero(para.strideNode()*t_save);
	VectorXd tSave = VectorXd::LinSpaced(t_save, tStart, tEnd);

	cache.c = start;
	while (t < tEnd) {
		// determine current timestep
		Dt = ((timestep(DG, para, CFL) + t > tSave[tSaveCount]) ? tSave[tSaveCount] - t : timestep(DG, para, CFL));
		//std::cout << "Timestep: " << Dt << std::endl;
		// calc convection dispersion rhs
		ConvDisp(cache, DG, para, t); // stores rhs in cache.dc
		// expl. Euler
		cache.c += Dt * cache.dc;
		//save outlet solution
		if (t == tSave[tSaveCount] && t_save > 1) {
			for (int comp = 0; comp < para.nComp; comp++) {
				outlet[(t_save - tSaveCount - 1) * para.strideNode() + comp]
					= cache.c[para.nCells * para.strideCell() - para.nComp + comp];
			}
		tSaveCount++;
		}
		t += Dt;
		tCount++;
		std::cout << "t=" << t << std::endl;
	}
	std::cout << "timesteps: " << tCount << std::endl;
	std::cout << "outlet saves: " << tSaveCount + 1 << std::endl;
	return outlet;
}
