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
		unsigned int polyDeg;
		double velocity;
		double dispersion;
		VectorXd porosity;
		VectorXd adsorption;
		VectorXd ADratio;
		std::string isotherm;
		// strides to switch to next entry in state vector
		inline int strideCell() { return (polyDeg + 1) * nComp; };
		inline int strideComp() { return 1; };
		inline int strideNode() { return nComp; };
		ParameterProvider(int nComp, int nCells, int polyDeg, double velocity, double disp, std::string isotherm = "Linear");
	};

class Container {
public:
	VectorXd c; //!< state vector of mobile phase
	VectorXd w; //!< mobile phase + solidphase rhs
	VectorXd S; //!< auxiliary variable du/dx
	VectorXd h; //!< substitute h = D_ax S - v u
	VectorXd surfaceFlux; //!< stores the surface flux values
	VectorXd boundary; //!< stores current boundary values
	Container(int nCells, int nNodes, int nComp);
};

// physical and numerical fluxes
typedef double (*Flux)(double point, ParameterProvider para);
double auxiliaryFlux(double point, ParameterProvider para);
double advectionDispersionFlux(double point, ParameterProvider para);
double dispersionFlux(double point, ParameterProvider para);
double advectionFlux(double point, ParameterProvider para);
typedef double (*riemannSolver)(double left, double right, Flux flux, ParameterProvider para);
double centralFlux(double left, double right, Flux flux, ParameterProvider para);
double laxFriedrichsFlux(double left, double right, Flux flux, ParameterProvider para);
// basis functions
void lglNodesWeights(int polyDeg, VectorXd& nodes, VectorXd& weights);
void polynomialDerivativeMatrix(int polyDeg, VectorXd& nodes, MatrixXd& D, MatrixXd& DT);
void qAndL(int polyDeg, double x, double* L, double* q, double* qder);

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

	double deltaX;
	Discretization(int polyDeg, double dX, riemannSolver numFlux = laxFriedrichsFlux);
};

// computation of physical nodes
VectorXd physNodes(double x_start, double x_end, ParameterProvider para, Discretization DG);