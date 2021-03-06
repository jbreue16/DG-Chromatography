
#define _USE_MATH_DEFINES

#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>

#include<Eigen/Dense>
using namespace Eigen;


/**
 * @brief computes the Legendre polynomial L_N and q = L_N+1 - L_N-2 and q' at point x
 * @param [in] polyDeg polynomial degree of spatial Discretization
 * @param [in] x evaluation point
 * @param [in] L pointer for allocated memory of L(x)
 * @param [in] q pointer for allocated memory of q(x)
 * @param [in] qder pointer for allocated memory of q'(x)
 */
void qAndL(int polyDeg, double x, double* L, double* q, double* qder) {
    // auxiliary variables (Legendre polynomials)
    double L_2 = 1.0;
    double L_1 = x;
    double Lder_2 = 0.0;
    double Lder_1 = 1.0;
    double Lder = 0.0;
    for (double k = 2; k <= polyDeg; k++) {
        *L = ((2 * k - 1) * x * L_1 - (k - 1) * L_2) / k;
        Lder = Lder_2 + (2 * k - 1) * L_1;
        L_2 = L_1;
        L_1 = *L;
        Lder_2 = Lder_1;
        Lder_1 = Lder;
    }
    *q = ((2.0 * polyDeg + 1) * x * *L - polyDeg * L_2) / (polyDeg + 1.0) - L_2;
    *qder = Lder_1 + (2.0 * polyDeg + 1) * L_1 - Lder_2;
}

/**
 * @brief computes and assigns the Legendre-Gauss nodes and weights
 * @param [in] polyDeg polynomial degree of spatial Discretization
 * @param [in] pre-initialized nodes array
 * @param [in] pre-initialized weights array
 */
void lglNodesWeights(int polyDeg, VectorXd& nodes, VectorXd& weights) {
    // tolerance and max #iterations for Newton iteration
    int nIterations = 10;
    double tolerance = 1e-15;
    // Legendre polynomial and derivative
    double L = 0;
    double q = 0;
    double qder = 0;
    switch (polyDeg) {
    case 0:
        throw std::invalid_argument("Polynomial degree must be at least 1 !");
        break;
    case 1:
        nodes[0] = -1;
        weights[0] = 1;
        nodes[1] = 1;
        weights[1] = 1;
        break;
    default:
        nodes[0] = -1;
        nodes[polyDeg] = 1;
        weights[0] = 2.0 / (polyDeg * (polyDeg + 1.0));
        weights[polyDeg] = weights[0];
        // use symmetrie, only compute half of points and weights
        for (int j = 1; j <= floor((polyDeg + 1) / 2) - 1; j++) {
            //  first guess for Newton iteration
            nodes[j] = -cos(M_PI * (j + 0.25) / polyDeg - 3 / (8.0 * polyDeg * M_PI * (j + 0.25)));
            // Newton iteration to find zero points of Legendre Polynomial
            for (int k = 0; k <= nIterations; k++) {
                qAndL(polyDeg, nodes[j], &L, &q, &qder);
                nodes[j] = nodes[j] - q / qder;
                if (abs(q / qder) <= tolerance * abs(nodes[j])) {
                    break;
                }
            }
            // calculate weights
            qAndL(polyDeg, nodes[j], &L, &q, &qder);
            weights[j] = 2 / (polyDeg * (polyDeg + 1.0) * pow(L, 2));
            nodes[polyDeg - j] = -nodes[j]; // copy to second half of points and weights
            weights[polyDeg - j] = weights[j];
        }
    }
    if (polyDeg % 2 == 0) { // for even polyDeg we have an odd number of points which include 0.0
        qAndL(polyDeg, 0.0, &L, &q, &qder);
        nodes[polyDeg / 2] = 0;
        weights[polyDeg / 2] = 2 / (polyDeg * (polyDeg + 1.0) * pow(L, 2));
    }
}

/**
 * @brief computation of barycentric weights for fast polynomial evaluation
 * @param [in] polyDeg polynomial degree of spatial Discretization
 * @param [in] nodes node vector
 * @param [in] weights pre-allocated vector for barycentric weights. Must be set to ones!
 */
void barycentricWeights(int polyDeg, VectorXd& nodes, VectorXd& baryWeights) {
    for (int j = 1; j <= polyDeg; j++) {
        for (int k = 0; k <= j - 1; k++) {
            baryWeights[k] = baryWeights[k] * (nodes[k] - nodes[j]) * 1.0;
            baryWeights[j] = baryWeights[j] * (nodes[j] - nodes[k]) * 1.0;
        }
    }
    for (int j = 0; j <= polyDeg; j++) {
        baryWeights[j] = 1 / baryWeights[j];
    }
}

/**
 * @brief initialization of polynomial derivative matrix D and D^T
 * @param [in] polyDeg polynomial degree of spatial Discretization
 * @param [in] nodes LGL node vector
 * @param [in] D pointer to pre-allocated derivative matrix
 * @param [in] DT pointer to pre-allocated transposed derivative matrix
 */
void polynomialDerivativeMatrix(int polyDeg, VectorXd& nodes, MatrixXd& D, MatrixXd& DT) {
    VectorXd baryWeights = VectorXd::Ones(polyDeg + 1);
    barycentricWeights(polyDeg, nodes, baryWeights);
    for (int i = 0; i <= polyDeg; i++) {
        for (int j = 0; j <= polyDeg; j++) {
            if (i != j) {
                D(i, j) = baryWeights[j] / (baryWeights[i] * (nodes[i] - nodes[j]));
                D(i, i) += - D(i, j);
            }
        }
    }
    DT = D.transpose();
}

/**
* TODO? polynomial evaluation
*/
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
    double porosity; //<! global porosity coefficient
    VectorXd adsorption; //<! adsoprtion for each component
    VectorXd ADratio; //<! adsoprtion/desorption ratio for each component
    std::string isotherm;
    // strides to switch to next entry in state vector
    inline int strideCell() { return (polyDeg + 1) * nComp; };
    inline int strideComp() { return 1; };
    inline int strideNode() { return nComp; };
    ParameterProvider(int nComp, int nCells, int polyDeg, double velocity, double disp,
                      double porosity = 0.0, std::string isotherm = "Langmuir");
};

ParameterProvider::ParameterProvider(int nComp, int nCells, int polyDeg, double v, double disp,
                                     double porosity, std::string isotherm)
    : nComp(nComp),
    nCells(nCells),
    nNodes(polyDeg + 1),
    polyDeg(polyDeg),
    velocity(v),
    dispersion(disp),
    porosity(porosity),
    adsorption(VectorXd::Zero(nComp)),
    ADratio(VectorXd::Zero(nComp)),
    isotherm(isotherm)
{
}
/**
*@brief container for state vector, auxiliary variable and rhs
*/
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
Container::Container(int nCells, int nNodes, int nComp)
    : c(VectorXd::Zero(nCells* nNodes* nComp)),
    dc(VectorXd::Zero(nCells* nNodes* nComp)),
    w(VectorXd::Zero(nCells* nNodes* nComp)),
    S(VectorXd::Zero(nCells* nNodes* nComp)),
    h(VectorXd::Zero(nCells* nNodes* nComp)),
    q(VectorXd::Zero(nCells* nNodes* nComp)),
    surfaceFlux(VectorXd::Zero(nComp* (nCells + 1))),
    boundary(VectorXd::Zero(4 * nComp))
{
}
/**
*@brief physical flux function
*/
typedef double (*Flux)(double point, ParameterProvider para);
double auxiliaryFlux(double point, ParameterProvider para) {
    return point;
}
double convectionFlux(double point, ParameterProvider para) {
    return para.velocity * point;
}
double dispersionFlux(double point, ParameterProvider para) {
    return sqrt(para.dispersion) * point;
}
double advectionDispersionFlux(double point, ParameterProvider para) {
    return dispersionFlux(point, para) - convectionFlux(point, para);
}

/**
* @brief numerical flux function
*/
typedef double (*riemannSolver)(double left, double right, Flux flux, ParameterProvider para);
/**
* @brief central numerical flux
*/
double centralFlux(double left, double right, Flux flux, ParameterProvider para) {
    return 0.5 * (flux(left, para) + flux(right, para));
}
/**
* @brief Lax-Friedrichs numerical flux
*/
double laxFriedrichsFlux(double left, double right, Flux flux, ParameterProvider para) {
    //double lambda = abs(para.velocity); // (local) dissipation parameter
    // or choose lambda = Delta x/ Delta t for global LF flux
    return 0.5 * (flux(left, para) + flux(right, para)) - (abs(para.velocity) /2) * (left - right);
}
/**
* @brief upwind numerical flux
*/
double upwindFlux(double left, double right, Flux flux, ParameterProvider para) {
    return flux(left, para);
}

/**
* @brief boundary functions 1D, returning c at inlet, outlet
*/
typedef double(*boundaryFunction)(double t, int component);
/**
* @brief Freestream boundary with all Component const 0.1
*/
double freestream01(double t, int component) {
    return 0.1;
}
/**
* @brief discontinous input pulse for one Component
*/
double pulse1Comp(double t, int component) {
    double bound;
    (t <= 1.0) ? bound = 1.0 : bound = 0.0;
    return bound;
}
/**
* @brief discontinous input pulse for two Components
*/
double pulse2Comp(double t, int component) {
    double bound = 0.0;
    if (component == 0) {
        if(t <= 12.0) bound = 10.00;
    }
    else {
        if (t <= 12.0) { bound = 10.00; }
        //if (1.0 <= t < 2.0) { bound = 3.0; }
    }
    return bound;
}
/**
* @brief boundary conditions returning the ghost nodes
* @param [in] boundary funtion giving c at inlet, outlet
* @param [in] t time
*/
typedef void(*boundaryCondition)(double t, Container& cache, boundaryFunction boundFunc, ParameterProvider para);
/**
* @brief implements the Danckwert boundary conditions
*/
void Danckwert(double t, Container& cache, boundaryFunction boundFunc, ParameterProvider para) {
    for (int comp = 0; comp < para.nComp; comp++) {
        cache.boundary[comp]                  = cache.c[comp]; // c_l inlet
        cache.boundary[para.nComp + comp]     = cache.c[para.nCells * para.strideCell() 
                                                        - para.strideNode() + comp]; // c_r outlet
        // store S_inlet
        //cache.boundary[2 * para.nComp + comp] = (-sqrt(para.dispersion) * cache.S[comp] - 2.0 * para.velocity *
        //    (cache.c[comp] + boundFunc(t, comp))) / ((para.dispersion == 0) ? 1.0 : sqrt(para.dispersion)); // S_l inlet
        // store h_inlet
        cache.boundary[2 * para.nComp + comp] = - cache.h[comp] + 2* para.velocity* boundFunc(t, comp); // h_l inlet = -h_r + 2vc_in
        cache.boundary[3 * para.nComp + comp] = cache.S[para.nCells * para.strideCell()
                                                        - para.strideNode() + comp]; // S_r outlet
    }
}
/**
* @brief implements Freeflow boundary conditions
*/
void Freeflow(double t, Container& cache, boundaryFunction boundFunc, ParameterProvider para) {
    for (int comp = 0; comp < para.nComp; comp++) {
        cache.boundary[comp] = cache.c[comp]; // c_l inlet
        cache.boundary[para.nComp + comp] = cache.c[para.nCells * para.strideCell()
                                                    - para.strideNode() + comp]; // c_r outlet
        cache.boundary[2 * para.nComp + comp] = cache.S[comp]; // S_l inlet
        cache.boundary[3 * para.nComp + comp] = cache.S[para.nCells * para.strideCell()
                                                        - para.strideNode() + comp]; // S_r outlet
    }
    // std::cout << "boundary values " << std::endl << cache.boundary << std::endl;
}
/**
* @brief implements Periodic boundary conditions
*/
void Periodic(double t, Container& cache, boundaryFunction boundFunc, ParameterProvider para) {
    for (int comp = 0; comp < para.nComp; comp++) {
        cache.boundary[comp] = cache.c[para.nCells * para.strideCell()
            - para.strideNode() + comp]; // c_l inlet
        cache.boundary[para.nComp + comp] = cache.c[comp]; // c_r outlet
        cache.boundary[2 * para.nComp + comp] = cache.S[para.nCells * para.strideCell()
            - para.strideNode() + comp]; // S_l inlet
        cache.boundary[3 * para.nComp + comp] = cache.S[comp]; // S_r outlet
    }
    // std::cout << "boundary values " << std::endl << cache.boundary << std::endl;
}
/**
 * @brief class to create DG objects containing all DG specifics
 * @detail is actually the Discretization _disc in Cadet
 */
class Discretization { // ~ _disc in CADET code
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

    Discretization(int polyDeg, double dX, riemannSolver numFlux = laxFriedrichsFlux,
                   boundaryCondition BoundCond = Danckwert, boundaryFunction BoundFunc = pulse1Comp);
};

void initializeDG(Discretization& dgsem) {
    lglNodesWeights(dgsem.polyDeg, dgsem.nodes, dgsem.weights);
    polynomialDerivativeMatrix(dgsem.polyDeg, dgsem.nodes, dgsem.polyDerM, dgsem.polyDerMtranspose);
    for (int n = 0; n <= dgsem.polyDeg; n++) {
        dgsem.invWeights[n] = 1.0 / dgsem.weights[n];
    }
}

Discretization::Discretization(int degree, double dX, riemannSolver numFlux,
                               boundaryCondition BoundCond, boundaryFunction BoundFunc)
    : polyDeg(degree),
    nNodes(degree + 1),
    nodes(VectorXd::Zero(nNodes)),
    weights(VectorXd::Ones(nNodes)),
    invWeights(VectorXd::Ones(nNodes)),
    polyDerM(MatrixXd::Zero(nNodes, nNodes)),
    polyDerMtranspose(MatrixXd::Zero(nNodes, nNodes)),
    deltaX(dX),
    numFlux(numFlux),
    BoundCond(BoundCond),
    BoundFunc(BoundFunc)
{
    // compute/initialize DG members 
    initializeDG(*this);
}

/**
*@brief returns the physical node coordinates
* @param [in] x_start start of interval
* @param [in] x_end end of interval
* @param [in] para Parameterprovider
* @param [in] DG DG Discretization
*/
VectorXd physNodes(double x_start, double x_end, ParameterProvider para, Discretization DG) {
    VectorXd x_l = VectorXd::LinSpaced(para.nCells + 1, x_start, x_end);
    double DeltaX = (x_l[1] - x_l[0]);
    VectorXd physNodes = VectorXd::Zero(para.nCells * DG.nNodes);
    for (int i = 0; i < para.nCells; i++) {
        for (int j = 0; j <= para.polyDeg; j++) {
            physNodes[i * DG.nNodes + j] = x_l[i] + 0.5 * DeltaX * (1 + DG.nodes[j]);// mapping 
        }
    }
    return physNodes;
}