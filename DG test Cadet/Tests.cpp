#define CATCH_CONFIG_MAIN
#define _USE_MATH_DEFINES

#include "C:\Users\jmbr\source\repos\DG test Cadet\DG test Cadet\Catch\catch.hpp"
#include "C:\Users\jmbr\source\repos\DG test Cadet\DG test Cadet\DGspecific.hpp"
#include "RHS.hpp"
#include "TimeIntegration.hpp"

#include "gnuplot-iostream.h"
#include "AnalysisTools.hpp"
//#define WITHOUT_NUMPY // somehow cant find numpy/arrayobject.h but it exists
//#include<python.h>
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
            derivative = derivative * (sqrt(para.dispersion)); // derivative independet factor that is applied to account for dispersion flux
            cache.boundary[0] = pow(phNodes[0], j);
            cache.boundary[1] = pow(phNodes[phNodes.size() - 1], j);
            ConvDisp(cache, DG, para, 0.0);

            //cout << " u: " << endl;
            //for (int i = 0; i < cache.S.size(); i++) {
            //    std::cout << cache.c[i] << std::endl;
            //}
            //std::cout << " surface flux: " << std::endl;
            //for (int i = 0; i < cache.surfaceFlux.size(); i++) {
            //    std::cout << cache.surfaceFlux[i] << std::endl;
            //}
            //cout << " S: " << endl;
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

    Discretization DG = Discretization(POLYDEG, deltaX, upwindFlux, Periodic, freestream01);
    ParameterProvider para = ParameterProvider(NCOMP, NCELLS, POLYDEG, velocity, dispersion);
    Container cache = Container(NCELLS, NNODES, NCOMP);

    VectorXd phNodes = physNodes(0.0, 1.0, para, DG);

    // constant conditions
    VectorXd start = 0.1 * VectorXd::Ones(cache.c.size());

    double tStart = 0.0;
    double tEnd = 10;
    double CFL = 0.1;

    VectorXd solution = solveEuler(cache, DG, para, tStart, tEnd, CFL, start);

    // test freestream preservation
    REQUIRE(start==solution);
}

TEST_CASE("RHS functions") {
    int POLYDEG = 4;
    int NCELLS = 20;// dont change!
    int NCOMP = 1;
    int NNODES = POLYDEG + 1;
    double deltaX = 2.0 * M_PI / NCELLS;
    double velocity = M_PI;
    double dispersion = 0.0;

    Discretization DG = Discretization(POLYDEG, deltaX, upwindFlux, Periodic, NULL);// Boundary function is irrelevant for periodic BC
    ParameterProvider para = ParameterProvider(NCOMP, NCELLS, POLYDEG, velocity, dispersion);
    Container cache = Container(NCELLS, NNODES, NCOMP);

    VectorXd phNodes = physNodes(0.0, 2.0 * M_PI, para, DG);
    // initial conditions and exact rhs solution
    VectorXd sinus = VectorXd::Zero(cache.c.size());
    VectorXd cosinus = VectorXd::Zero(cache.c.size());
    for (int i = 0; i < phNodes.size(); i++) {
        sinus[i] = sin(phNodes[i]);
        cosinus[i] = cos(phNodes[i]);
    }

    // TEST RHS ConvDisp function calls one by one 
    // 1) test boundary value estimation (periodic boundary)
    cache.c[0] = 12.0; cache.c[cache.c.size() - 1] = 13.0;
    cache.S[0] = 22.0; cache.S[cache.c.size() - 1] = 23.0;
    DG.BoundCond(0.0, cache, DG.BoundFunc, para);
    VectorXd boundTest(4);
    boundTest << 13.0, 12.0, 23, 22;
    REQUIRE(cache.boundary == boundTest);

    // 2) test calculation of h = vc - sqrt(D_ax) S
    cache.c.Random(NCELLS * NNODES * NCOMP); cache.S.Random(NCELLS * NNODES * NCOMP);
    VectorXd calcHTest = para.velocity * cache.c - sqrt(para.dispersion) * cache.S;
    calcH(cache, para);
    REQUIRE(cache.h == calcHTest);
    // 3) test volume integral [- D * state] and applyMapping [* 2/Dx]
    cache.c = sinus;
    cache.dc.Zero(NCELLS * NNODES * NCOMP);
    volumeIntegral(cache, DG, para, cache.c, cache.dc);
    applyMapping(DG, cache.dc);
    REQUIRE(cache.dc.isApprox(-cosinus, 1e-5)); // tolerance in relation to polydeg

    // 4) test Interface Flux
    cache.c = VectorXd::Random(size(cache.c)); cache.surfaceFlux.Zero(NCELLS + 1);
        // periodic boundary; could also use DG.BoundCond(0.0, cache, DG.BoundFunc, para);
    cache.boundary[0] = cache.c[NCELLS * NNODES * NCOMP - 1]; cache.boundary[1] = cache.c[0];
    InterfaceFlux(cache, DG, para);
    VectorXd upwind(NCELLS + 1);
    for (int cell = 1; cell < NCELLS; cell++) {
        upwind[cell] = para.velocity * cache.c[cell * para.strideCell() - 1];// left interface nodes for upwind scheme
    }
    upwind[0] = para.velocity * cache.boundary[0];
    upwind[NCELLS] = para.velocity * cache.c[NCELLS * NNODES * NCOMP - 1];
    REQUIRE(cache.surfaceFlux == upwind);

    // 5) test surface integral M^-1 B (h* - h)
    cache.c = sinus; cache.surfaceFlux.Zero(NCELLS + 1);
    VectorXd surfIntContribution = VectorXd::Zero(cache.c.size());
    VectorXd test = VectorXd::Zero(cache.c.size());
    // add discontinuities at cell interfaces by shifting nodes
    for (int cell = 0; cell < floor(para.nCells / 4); cell++) {
        cache.c(cell * para.strideCell() + para.polyDeg * para.strideNode()) += 0.2;
    }
    for (int cell = floor(para.nCells / 4); cell < floor(para.nCells / 2); cell++) {
        cache.c(cell * para.strideCell()) += 0.1;
    }
    for (int cell = floor(para.nCells / 2); cell < floor(3*para.nCells / 4); cell++) {
        cache.c(cell * para.strideCell()) += 0.3;
        cache.c(cell * para.strideCell() + para.polyDeg * para.strideNode()) += 0.2;
    }
    for (int cell = floor(3*para.nCells / 4); cell < para.nCells; cell++) {
        cache.c(cell * para.strideCell()) += 0.1;
        cache.c(cell * para.strideCell() + para.polyDeg * para.strideNode()) += 0.1;
    }
    // calculated by hand
    VectorXd eingabe = VectorXd::Zero(20);
    eingabe << 1, 2, 2, 2, 2, 1, -1, -1, -1, -1, -3, -1, -1, -1, -1, 1, 0.0, 0.0, 0.0, 0.0;
    for (int cell = 0; cell < NCELLS; cell++) {
        test[cell * para.strideCell()] = eingabe[cell]*velocity;
    }
    DG.BoundCond(0.0, cache, DG.BoundFunc, para);
    calcH(cache, para);
    surfaceIntegral(cache, DG, para, cache.h, surfIntContribution, 0);
    REQUIRE(test.isApprox(surfIntContribution, 1e-14)); // machine accuracy
}

TEST_CASE("Advection") {
        int POLYDEG = 4;
        int NCELLS = 20;
        int NCOMP = 1;
        int NNODES = POLYDEG + 1;
        double deltaX = 2.0 * M_PI / NCELLS;
        double velocity = M_PI;
        double dispersion = 0.0;
    
        Discretization DG = Discretization(POLYDEG, deltaX, upwindFlux, Periodic, pulse1Comp);// Boundary function is irrelevant for periodic BC
        ParameterProvider para = ParameterProvider(NCOMP, NCELLS, POLYDEG, velocity, dispersion);
        Container cache = Container(NCELLS, NNODES, NCOMP);
    
        VectorXd phNodes = physNodes(0.0, 2.0 * M_PI, para, DG);

        // initial conditions
        VectorXd start = VectorXd::Zero(cache.c.size());
        VectorXd analyticalSol = VectorXd::Zero(cache.c.size());
        for (int i = 0; i < phNodes.size(); i++) {
            start[i] = sin(phNodes[i]);
            analyticalSol[i] = -velocity* cos(phNodes[i]);
        }
        cache.c = start;

    double tStart = 0.0;
    double CFL = 0.1;
    double tEnd = 2.0;

    solveEuler(cache, DG, para, tStart, tEnd, CFL, start);
    //Plotter plot;
    //plot.line_plot(cache.c, phNodes, "bitte", "x", "C");

    REQUIRE(cache.c.isApprox(start, 1e-1));
}

//TEST_CASE("Spielwiese") {
//    int POLYDEG = 4;
//    int NCELLS = 20;
//    int NCOMP = 1;
//    int NNODES = POLYDEG + 1;
//    double deltaX = 2.0 * M_PI / NCELLS;
//    double velocity = M_PI;
//    double dispersion = 0.1;
//
//    Discretization DG = Discretization(POLYDEG, deltaX, upwindFlux, Periodic, pulse1Comp);// Boundary function is irrelevant for periodic BC
//    ParameterProvider para = ParameterProvider(NCOMP, NCELLS, POLYDEG, velocity, dispersion);
//    Container cache = Container(NCELLS, NNODES, NCOMP);
//
//    VectorXd phNodes = physNodes(0.0, 2.0 * M_PI, para, DG);
//
//    // initial conditions
//    VectorXd start = VectorXd::Zero(cache.c.size());
//    VectorXd analyticalSol = VectorXd::Zero(cache.c.size());
//    for (int i = 0; i < phNodes.size(); i++) {
//        start[i] = sin(phNodes[i]);
//        analyticalSol[i] = -velocity * cos(phNodes[i]);
//    }
//    cache.c = start;
//
//    double tStart = 0.0;
//    double CFL = 0.1;
//    double tEnd = 2.0;
//
//    solveEuler(cache, DG, para, tStart, tEnd, CFL, start);
//    Plotter plot;
//    plot.line_plot(cache.c, phNodes, "bitte", "x", "C");
//
//}