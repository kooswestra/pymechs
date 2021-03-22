#ifndef PID_H
#define PID_H

#include "controller.h"

namespace Mech {

// Controller generates control outputs from the state of the mechanism and a reference (trajectory)
class Pid : public Controller
{
public:
    Pid(double time, Eigen::VectorXd parameters);
    double controlOutput(double t, ConstVectorRef& q, ConstVectorRef& qd) override;

    // Controller dimensionality
    int dim() const;
    Eigen::VectorXd parameters() const;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    int dim_ = 1;

    double time_;
    Eigen::VectorXd parameters_;

};

}

#endif // PID_H
