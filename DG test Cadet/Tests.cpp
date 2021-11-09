#define CATCH_CONFIG_MAIN

#include "C:\Users\jmbr\source\repos\DG test Cadet\DG test Cadet\Catch\catch.hpp"
#include "C:\Users\jmbr\source\repos\DG test Cadet\DG test Cadet\DGspecific.hpp"
#include "RHS.hpp"
#include "TimeIntegration.hpp"
//#define WITHOUT_NUMPY // somehow cant find numpy/arrayobject.h file despite it exists in the fucking folder im including
//#include<matplotlibcpp.h>

using namespace std;

TEST_CASE("Test DG basic functions") {

    for (int POLYDEG = 1; POLYDEG <= 10; POLYDEG++) {
        bool test = true;
        int NCELLS = 1;
        Discretization dgsem(POLYDEG, 1 / NCELLS);

        Eigen::MatrixXd polynomials = Eigen::MatrixXd::Zero(POLYDEG + 1, POLYDEG + 2);
        Eigen::MatrixXd derivative = Eigen::MatrixXd::Zero(POLYDEG + 1, POLYDEG + 1);
        Eigen::MatrixXd num_derivative = Eigen::MatrixXd::Zero(POLYDEG + 1, POLYDEG + 1);

        for (int j = 0; j <= POLYDEG; j++) {
            for (int i = 0; i <= POLYDEG; i++) {
                polynomials(i, j + 1) = pow(dgsem.nodes[i], j);
            }
            num_derivative.col(j) = dgsem.polyDerM * polynomials.col(j + 1);
            derivative.col(j) = (j)*polynomials.col(j);
        }
        if (!num_derivative.isApprox(derivative, 1e-15)) {
            test = false;
        }
        REQUIRE(test);
    }
}

TEST_CASE("Test auxiliary equation solver") {

    for (int POLYDEG = 1; POLYDEG <= 10; POLYDEG++) {
        int NCELLS = 3;
        int NCOMP = 1;
        int NNODES = POLYDEG + 1;
        double deltaX = 1.0 / NCELLS;
        double velocity = 1.0;
        double dispersion = 0.01;

        Discretization DG = Discretization(POLYDEG, deltaX);
        ParameterProvider para = ParameterProvider(NCOMP, NCELLS, POLYDEG, velocity, dispersion);
        Container cache = Container(NCELLS, NNODES, NCOMP);

        VectorXd phNodes = physNodes(0.0, 1.0, para, DG);
        VectorXd derivative = VectorXd::Zero(phNodes.size());

        for (int j = 1; j <= POLYDEG; j++) {
            // fill the state vector u and boundary values in Container
            for (int i = 0; i < cache.c.size(); i++) {
                cache.c[i] = pow(phNodes[i], j);
                derivative[i] = j * pow(phNodes[i], j - 1);
            }
            cache.boundary[0] = pow(phNodes[0], j);
            cache.boundary[1] = pow(phNodes[phNodes.size() - 1], j);
            ConvDisp(cache, DG, para, 0.0);

            //cout << " u: " << endl;
            //for (int i = 0; i < cache.S.size(); i++) {
            //    std::cout << cache.u[i] << std::endl;
            //}
            //std::cout << " surface flux: " << std::endl;
            //for (int i = 0; i < cache.surfaceFlux.size(); i++) {
            //    std::cout << cache.surfaceFlux[i] << std::endl;
            //}
            //cout << " q: " << endl;
            //for (int i = 0; i < cache.S.size(); i++) {
            //    std::cout << cache.S[i] << std::endl;
            //}
            //cout << " derivative: " << endl;
            //for (int i = 0; i < cache.S.size(); i++) {
            //    std::cout << derivative[i] << std::endl;
            //}

            // TODO: somehow tolerance e-14 neccessery despite calculated and displayed values are exact
            REQUIRE(cache.S.isApprox(derivative, 1e-14)); 
        }
    }
}

TEST_CASE("Test isotherm Jacobian") {
    int POLYDEG = 3;
    int NCELLS = 3;
    int NCOMP = 3;
    int NNODES = POLYDEG + 1;
    double deltaX = 1.0 / NCELLS;
    double velocity = 0.1;
    double dispersion = 0.01;
    double porosity = 0.35;

    Discretization DG = Discretization(POLYDEG, deltaX);
    Container cache = Container(NCELLS, NNODES, NCOMP);
    try {
        ParameterProvider para = ParameterProvider(NCOMP, NCELLS, POLYDEG, velocity, dispersion, porosity, "Lineal");
        ConvDisp(cache, DG, para, 0.0);
    }
    catch (std::invalid_argument const& err) {
        REQUIRE(err.what() == std::string("spelling error or this isotherm is not implemented yet"));
    }
    ParameterProvider para = ParameterProvider(NCOMP, NCELLS, POLYDEG, velocity, dispersion, porosity, "Langmuir");
    para.adsorption[0] = 0.2; // "wenig" adsorption
    para.adsorption[1] = 0.5;
    para.adsorption[2] = 0.8; // "viel" adsorption
    para.ADratio[0] = 0.1; // starke desorption zu adsorption
    para.ADratio[1] = 0.5;
    para.ADratio[2] = 0.9; // schwache desorption zu adsorption
    VectorXd phNodes = physNodes(0.0, 1.0, para, DG);
    for (int i = 0; i < DG.nNodes* para.nCells; i++) {
        cache.c[3 * i] = 0.5; // Comp 1
        cache.c[3 * i + 1] = 1; // Comp 2
        cache.c[3 * i + 2] = 1.5; // Comp 3
    }
    MatrixXd Jac = calcDQDC(cache, para, 0, 0);
    MatrixXd testJac = MatrixXd::Zero(para.nComp, para.nComp);
    // testJac = F*(dq/dc) + I)
    testJac << para.porosity*0.067776456 + 1.0, para.porosity * -0.0059453032, para.porosity * -0.010701545,
              para.porosity * -0.0059453032, para.porosity * 0.142687277 + 1.0, para.porosity * -0.053507728,
              para.porosity * -0.014268727, para.porosity * -0.071343638, para.porosity * 0.147443519 + 1.0;
    REQUIRE(Jac.isApprox(testJac.inverse(), 1e-8)); // tolerance somehow digits-1 exact.. 
}

TEST_CASE("Free-Stream Preservation") {
    int POLYDEG = 3;
    int NCELLS = 10;
    int NCOMP = 1;
    int NNODES = POLYDEG + 1;
    double deltaX = 1.0 / NCELLS;
    double velocity = 0.8;
    double dispersion = 0.0; // Freestream only without dispersion !

    Discretization DG = Discretization(POLYDEG, deltaX, laxFriedrichsFlux, Freeflow, freestream01);
    ParameterProvider para = ParameterProvider(NCOMP, NCELLS, POLYDEG, velocity, dispersion);
    Container cache = Container(NCELLS, NNODES, NCOMP);

    VectorXd phNodes = physNodes(0.0, 1.0, para, DG);

    // constant conditions
    VectorXd start = 0.1 * VectorXd::Ones(cache.c.size());

    double tStart = 0.0;
    double tEnd = 10;
    double CFL = 0.95;

    VectorXd solution = solveRK4(cache, DG, para, tStart, tEnd, CFL, start);

    // test freestream preservation
    REQUIRE(start == solution);
}

//TEST_CASE("Spielwiese") {
//        ParameterProvider para = ParameterProvider(1, 2, 3, 1.0, 0.01);
//        Container cache = Container(2, 4, 1);
//        std::cout << " size surface flux: " << cache.surfaceFlux.size() << std::endl;
//}

