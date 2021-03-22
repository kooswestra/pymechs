#include "feedforward.h"

using namespace Eigen;
using namespace Mech;

FeedForward::FeedForward(double w0, VectorXd amplitudes, VectorXd phase)
    :   amplitudes_(amplitudes),
        phase_(phase),
        w0_(w0)
{
    if(amplitudes.size() != phase.size())
        throw std::invalid_argument("The feed-forward controller is poorly defined, check the order");

}

double FeedForward::controlOutput(double t, ConstVectorRef& /*q*/, ConstVectorRef& /*qd*/)
{
    double output = 0;

    for(int i = 0; i < amplitudes_.size(); i++) {
        output += amplitudes_(i)*cos((i + 1)*w0_*t + phase_(i));
    }

    control_effort_ += abs(output);

    return output;
}

int FeedForward::order() const
{
    return amplitudes_.size();
}

double FeedForward::controlEffort() const
{
    return control_effort_;
}
