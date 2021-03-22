#include "controller.h"

Controller::Controller(ControllerDNA dna)
    :   dna_(dna)
{
    // Construct the controller structure from the controller dna
    initializeStructure();
}

double Controller::controlOutput(double q, double qd)
{
    control_output_ = dna_.K()(0)*q + dna_.K()(1)*qd;

    return control_output_;
}

int Controller::dim() const
{
    return dim_;
}

void Controller::initializeStructure() 
{

}
