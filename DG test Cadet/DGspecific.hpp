#pragma once
#include<Eigen/Dense>
using namespace Eigen;


/**
* @brief vgl cadet parameterprovider
*/
class ParameterProvider {
public:
	unsigned int nComp;
	unsigned int nCells;
	unsigned int nNodes;
	unsigned int polyDeg;
	double velocity;
	double dispersion;
	double porosity;
	VectorXd adsorption;
	VectorXd ADratio;
	std::string isotherm;
	// strides to switch to next entry in state vector
	inline int strideCell() { return (polyDeg + 1) * nComp; };
	inline int strideComp() { return 1; };
	inline int strideNode() { return nComp; };
	ParameterProvider(int nComp, int nCells, int polyDeg, double velocity, double disp, double porosity = 0.0, std::string isotherm = "Langmuir");
};

class Container {
public:
	VectorXd c; //!< state vector of mobile phase
	VectorXd dc; //!< state derivatove vector of mobile phase
	VectorXd w; //!< mobile phase + solidphase rhs
	VectorXd S; //!< auxiliary variable du/dx
	VectorXd h; //!< substitute h = D_ax S - v u
	VectorXd q; //!< stationary phase q
	VectorXd surfaceFlux; //!< stores the surface flux values
	VectorXd boundary; //!< stores current boundary values for c and S -> [c_1,...c_n, S_1,...,S_n]
	Container(int nCells, int nNodes, int nComp);
};

// physical and numerical fluxes
typedef double (*Flux)(double point, ParameterProvider para);
double auxiliaryFlux(double point, ParameterProvider para);
double advectionDispersionFlux(double point, ParameterProvider para);
double dispersionFlux(double point, ParameterProvider para);
double convectionFlux(double point, ParameterProvider para);
typedef double (*riemannSolver)(double left, double right, Flux flux, ParameterProvider para);
double centralFlux(double left, double right, Flux flux, ParameterProvider para);
/**
* @brief Lax-Friedrichs numerical flux
*/
double laxFriedrichsFlux(double left, double right, Flux flux, ParameterProvider para);
/**
* @brief upwind numerical flux
*/
double upwindFlux(double left, double right, Flux flux, ParameterProvider para);
// basis functions
void lglNodesWeights(int polyDeg, VectorXd& nodes, VectorXd& weights);
void polynomialDerivativeMatrix(int polyDeg, VectorXd& nodes, MatrixXd& D, MatrixXd& DT);
void qAndL(int polyDeg, double x, double* L, double* q, double* qder);
//boundary treatment
typedef double(*boundaryFunction)(double t, int component);
double freestream01(double t, int component);
double pulse1Comp(double t, int component);
double pulse2Comp(double t, int component);
/**
*@brief Boundary condition
*/
typedef void(*boundaryCondition)(double t, Container& cache, boundaryFunction boundFunc, ParameterProvider para);
/**
*@brief Danckwert boundary condition
*/
void Danckwert(double t, Container& cache, boundaryFunction boundFunc, ParameterProvider para);
/**
*@brief Freeflow boundary condition
*/
void Freeflow(double t, Container& cache, boundaryFunction boundFunc, ParameterProvider para);
/**
* @brief implements Periodic boundary conditions
*/
void Periodic(double t, Container& cache, boundaryFunction boundFunc, ParameterProvider para);

class Discretization {  // ~ vgl Discretization in cadet
public:

	unsigned int polyDeg; //!< polynomial degree
	unsigned int nNodes; //!< Number of nodes in each cell = polynomial degree - 1
	Eigen::VectorXd nodes;		//!< Array with positions of nodes in reference element
	Eigen::VectorXd weights; //!< Array with weights for numerical quadrature of size nNodes
	Eigen::VectorXd invWeights; //!< Array with inverse weights for numerical quadrature of size nNodes
	Eigen::MatrixXd polyDerM; //!< Array with polynomial derivative Matrix of size nNodes^2
	Eigen::MatrixXd polyDerMtranspose; //!< Array with D^T
	riemannSolver numFlux; //!< numerical flux to serve as Riemann solver
	boundaryFunction BoundFunc; //!< boundary function
	boundaryCondition BoundCond; //!< boundary condition
	double deltaX; //<! cell spacing

	Discretization(int polyDeg, double dX, riemannSolver numFlux = laxFriedrichsFlux, boundaryCondition BoundCond = Danckwert, boundaryFunction BoundFunc = pulse1Comp);
};

// computation of physical nodes
VectorXd physNodes(double x_start, double x_end, ParameterProvider para, Discretization DG);