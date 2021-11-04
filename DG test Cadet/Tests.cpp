#define CATCH_CONFIG_MAIN
#include "C:\Users\jmbr\source\repos\DG test Cadet\DG test Cadet\Catch\catch.hpp"
#include "C:\Users\jmbr\source\repos\DG test Cadet\DG test Cadet\DGspecific.hpp"
#include "RHS.hpp"

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
            for (int i = 0; i < cache.u.size(); i++) {
                cache.u[i] = pow(phNodes[i], j);
                derivative[i] = j * pow(phNodes[i], j - 1);
            }
            cache.boundary[0] = pow(phNodes[0], j);
            cache.boundary[1] = pow(phNodes[phNodes.size() - 1], j);
            ConvDisp(cache, DG, para);

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

//TEST_CASE("Spielwiese") {
//        ParameterProvider para = ParameterProvider(1, 2, 3, 1.0, 0.01);
//        Container cache = Container(2, 4, 1);
//        std::cout << " size surface flux: " << cache.surfaceFlux.size() << std::endl;
//}

