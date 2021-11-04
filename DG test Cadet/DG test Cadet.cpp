#define _USE_MATH_DEFINES

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>

#include<Eigen/Dense>

using namespace Eigen;

#include "DGspecific.hpp"
// #include "mkl.h" unnoetig -> use eigen

#define min(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;



/*int main() {

    // Test derivative Matrix
    int POLYDEG = 5;
    int nCells = 10;
    Discretization dgsem(POLYDEG, nCells);

        // initialize test functions
    Eigen::VectorXd constant = Eigen::VectorXd::Zero(POLYDEG + 1);
    Eigen::VectorXd linear = Eigen::VectorXd::Zero(POLYDEG + 1);
    Eigen::VectorXd square = Eigen::VectorXd::Zero(POLYDEG + 1);
    Eigen::VectorXd cubic = Eigen::VectorXd::Zero(POLYDEG + 1);
    for (int i = 0; i <= POLYDEG; i++) {
        constant[i] = 1;
        linear[i] = dgsem.nodes[i];
        square[i] = pow(dgsem.nodes[i], 2);
        cubic[i] = pow(dgsem.nodes[i], 3);
    }

    Eigen::VectorXd solConst = dgsem.polyDerM * constant;
    Eigen::VectorXd solLin = dgsem.polyDerM * linear;
    Eigen::VectorXd solSq = dgsem.polyDerM * square;
    Eigen::VectorXd solCubic = dgsem.polyDerM * cubic;
    for(int i = 0; i <= POLYDEG; i++) {
        cout << solConst[i] << "    "<< solLin[i] << "    " << solSq[i] << "    " << solCubic[i] << endl;
    }
    

    /* Test DG initialization
    int POLYDEG = 4;
    Discretization dgsem(POLYDEG, 10);
    for (int i = 0; i <= dgsem.polyDeg; i++) {
        cout << "node " << i+1 << ": " << dgsem.nodes[i] << endl;
    }
    for (int i = 0; i <= dgsem.polyDeg; i++) {
        cout << "weight " << i + 1 << ": " << dgsem.weights[i] << endl;
    }
    for (int i = 0; i <= dgsem.polyDeg; i++) {
        for (int j = 0; j <= dgsem.polyDeg; j++) {
            cout << dgsem.polyDerM(i, j) << " ";
        }
        cout << endl;
    }*/
    
    /* Test lglNodesWeights()->verified !
    const int polyDeg = 4;
    double nodes[polyDeg + 1] = { 0 };
    double weights[polyDeg + 1] = { 0 };
    lglNodesWeights(polyDeg, nodes, weights);
    for (int i = 0; i <= polyDeg; i++) {
        cout << "Knoten " << i << ": " << nodes[i] << " ";
    }
    cout << endl;
    for (int i = 0; i <= polyDeg; i++) {
        cout << "Gewicht " << i << ": " << weights[i] << " ";
    }*/


    /* Test qAndL()->verified !
    int polyDeg = 2; // must be at least 2 ! (for polyDeg = 1 this function is never called)
    int x = 0.0;
    double L = 0;
    double q = 0;
    double qder = 0;

    qAndL(polyDeg, x, &L, &q, &qder);
    cout << "polyDeg = " << polyDeg << " x= " << x << endl;
    cout << "L = " << L << endl;
    cout << "q = " << q << " (L_N+1 - L_N-2)" <<endl;
    cout << "q' = " << qder << endl;
    */

/*}
*/

// OLD CODE DEPOSIT //

/**
 * @brief computes the Legendre polynomial and derivative at point x
 * @param [in] polyDeg polynomial degree
 * @param [in] x evaluation point
 * @param [in] L pointer for allocated memory of L(x)
 * @param [in] Lder pointer for allocated memory of L'(x)
 */
 /*void legendrePolynomialAndDerivative(int polyDeg, double x, double* L, double* Lder) {
     switch (polyDeg) {
     case 0:
         *L = 1;
         *Lder = 0;
         break;
     case 1:
         *L = x;
         *Lder = 1;
         break;
     default:
         double L2 = 1;
         double L1 = x;
         double Lder2 = 0;
         double Lder1 = 1;
         for (double k = 2; k <= polyDeg; k++) {
             *L = ((2 * k - 1) * x * L1 - (k - 1) * L2) / k;
             *Lder = Lder2 + (2 * k - 1) * L1;
             L2 = L1;
             L1 = *L;
             Lder2 = Lder1;
             Lder1 = *Lder;
         }
     }
 }
 */
