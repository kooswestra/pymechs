/* 
 * File:   main.cpp
 * Author: westra
 *
 * Created on November 28, 2017, 5:45 PM
 */

#include <iostream>
#include <eigen3/Eigen/Dense>
#include "Evolution/dna.h"
#include "Evolution/mechanism.h"
#include "Evolution/objective.h"
#include "Evolution/mutator.h"
#include "Evolution/pickandplace.h"

using namespace Eigen;

int main() {

    /** A simple test program for the c++ library while more rigorous testing is still absent */

    Matrix<int,3,3> incidence_matrix;
    incidence_matrix << 0, 1, 1,
                        1, 1, 1,
                        1, 0, 0;


    VectorXi edge_labels(3);
    edge_labels << 1, 0, 2;

    RowVectorXd mass(3);
    mass << -0.5, -0.5, 5;

    std::vector<RowVectorXd> masses = {mass};

    RowVectorXd H(2);
    H << 0, 0;

    RowVectorXd E(3);
    E << 1, -1, 0.2;

    RowVectorXd K(6);
    K << 0, 0, -0.5, -1, 1, 2;

    std::vector<RowVectorXd> parameters = {E, H, K};

    // Test DNA structure
    Mech::DNA dna(incidence_matrix, edge_labels, masses, parameters);

    VectorXi labels(2);
    labels << 5, 3;

    Mech::PickAndPlace testt({1, 2}, {-2, 0});

    Mech::Mechanism mech(dna);
    mech.simulate(10, 200);
    testt.evaluate(&mech);

}

