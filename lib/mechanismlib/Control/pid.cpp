#include "pid.h"

using namespace Eigen;
using namespace Mech;

Pid::Pid(double time, VectorXd parameters)
    :  time_(time),
       parameters_(parameters)
{

}

double Pid::controlOutput(double t, ConstVectorRef& q, ConstVectorRef& qd)
{
    double Kp = parameters_(0);
    double Kd = parameters_(1);
    double ref1 = parameters_(2);

    // Ref 2 should always be the starting point, i.e. theta = 0
    if (t < time_ / 2) {
        return Kp * (ref1 - q(q.size() - 1)) - Kd * qd(qd.size() - 1);
    }
    else {
        return Kp * (0 - q(q.size() - 1)) - Kd * qd(qd.size() - 1);
    }
}

int Pid::dim() const
{
    return 1;
}

VectorXd Pid::parameters() const
{
    return parameters_;
}