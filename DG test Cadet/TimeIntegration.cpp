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
double timestep(Discretization DG, ParameterProvider para, double CFL, double tend) {
	// Delta X_eff = (DG.deltaX/DG.nNodes)
	return CFL* (DG.deltaX/DG.nNodes) / abs(para.velocity);
}
/**
*@brief 3rd order expl. RK
*/
Eigen::VectorXd solveRK3(Container& cache, Discretization& DG, ParameterProvider& para, double tstart, double tend, double CFL, VectorXd start) {
	double t = tstart;
	double Dt = 0.0;
	cache.c = start;

	VectorXd v1 = VectorXd::Zero(cache.c.size());
	VectorXd v2 = VectorXd::Zero(cache.c.size());

	while(t < tend) {
		// determine current timestep
		Dt = (timestep(DG, para, CFL, tend) + t > tend) ? tend - t : timestep(DG, para, CFL, tend);
		//std::cout << "Timestep: " << Dt << std::endl;
		// calc convection dispersion rhs
		ConvDisp(cache, DG, para, t); // stores rhs in cache.dc

		// expl. 3-stage RK
		v1 = cache.c + Dt * cache.dc;
		v2 = cache.c + Dt * (0.25 * cache.dc + 0.25 * v1);
		cache.c = cache.c + Dt * ((1.0 / 6.0) * cache.dc + (2.0 / 3.0) * v2 + (1.0 / 6.0) * v1);

		t += Dt;
	}
	return cache.c;
}

/**
*@brief Carpenter Kennedy 4th order low storage expl. RK
*/
Eigen::VectorXd solveRK4(Container& cache, Discretization& DG, ParameterProvider& para, double tstart, double tend, double CFL, VectorXd start) {
	double t = tstart;
	double Dt = 0.0;
	cache.c = start;

	// RK coefficients A, B, C for each stage
	VectorXd A(5);
	A << 0.0, -567301805773 / 1357537059087, -2404267990393 / 2016746695238,
		-3550918686646 / 2091501179385, -1275806237668 / 842570457699;
	VectorXd B(5);
	B << 0.0, 1432997174477 / 9575080441755, 2526269341429 / 6820363962896,
		2006345519317 / 3224310063776, 2802321613138 / 2924317926251;
	VectorXd C(5);
	C << 1432997174477 / 9575080441755, 5161836677717 / 13612068292357,
		1720146321549 / 2090206949498, 3134564353537 / 4481467310338, 2277821191437 / 14882151754819;
	VectorXd g = VectorXd::Zero(cache.c.size());
	while (t < tend) {
		// determine current timestep
		Dt = (timestep(DG, para, CFL, tend) + t > tend) ? tend - t : timestep(DG, para, CFL, tend);
		//std::cout << "Timestep: " << Dt << std::endl;
		// calc convection dispersion rhs
		ConvDisp(cache, DG, para, t); // stores rhs in cache.dc
		g = cache.dc;
		// expl. 3-stage RK
		for (int stage = 0; stage < 5; stage++) {
			ConvDisp(cache, DG, para, t + B[stage] * Dt);
			g = A[stage] * g + cache.dc;
			cache.c += C[stage] * Dt * g;
		}

		t += Dt;
	}
	return cache.c;
}

