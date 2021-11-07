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
	return CFL * DG.deltaX / abs(para.velocity);
}

/**
*@brief calculates boundary values for current timestep
*/
void calcBoundary(Container& cache, Discretization DG, double t) {
	cache.boundary = DG.BoundFunc(t);
}

Eigen::VectorXd solve(Container& cache, Discretization& DG, ParameterProvider& para, double tstart, double tend, double CFL, VectorXd start) {
	double t = tstart;
	double Dt = 0.0;
	cache.c = start;


	VectorXd v1 = VectorXd::Zero(cache.c.size());
	VectorXd v2 = VectorXd::Zero(cache.c.size());

	while(t < tend) {
		Dt = (timestep(DG, para, CFL, tend) + t > tend) ? t - tend : timestep(DG, para, CFL, tend);
		calcBoundary(cache, DG, t);

		ConvDisp(cache, DG, para); // stores rhs in cache.dc

		// expl. 3-stage RK
		v1 = cache.c + Dt * cache.dc;

		v2 = cache.c + Dt * (0.25 * cache.dc + 0.25 * v1);
		cache.c = cache.c + Dt * ((1.0 / 6.0) * cache.dc + (2.0 / 3.0) * v2 + (1.0 / 6.0) * v1);

		t += Dt;
	}
	return cache.c;
}


