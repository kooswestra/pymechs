#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "controllerdna.h"

// Controller generates control outputs from the state of the mechanism and a reference trajectory
class Controller
{
public:
    Controller(ControllerDNA dna);
    double controlOutput(double q, double qd);

    // Controller dimensionality
    int dim() const;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    int dim_;
    double control_output_;

    ControllerDNA dna_;
    // Possible reference trajectory, either time based or just a fixed reference
    //MatrixXd reference_trajectory_;

    void initializeStructure();

};

#endif // CONTROLLER_H
