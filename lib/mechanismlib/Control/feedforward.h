#ifndef FEEDFORWARD_H
#define FEEDFORWARD_H

#include "../common.h"
#include "controller.h"

namespace Mech {

class FeedForward : public Controller
{
public:
    FeedForward(double w0, Eigen::VectorXd amplitudes, Eigen::VectorXd phase);
    double controlOutput(double t, ConstVectorRef& q, ConstVectorRef& qd) override;

    int order() const;
    double controlEffort() const;

private:
    Eigen::VectorXd amplitudes_;
    Eigen::VectorXd phase_;
    double w0_;
    double control_effort_ = 0;
};

}

#endif // FEEDFORWARD_H
