#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "../common.h"
#include <Eigen/Dense>

/** @brief Controller is the abstract interface for controller implementations
 *
 * Each controller has to implement a controlOutput function which provides a torque as a
 * function of the state of the mechanism. Two implementations are provided, a python bridge
 * to define the controlOutput in python, and a simple PID controller.
 */

namespace Mech {

class Controller
{
public:
    Controller();
    virtual ~Controller();
    virtual double controlOutput(double t, ConstVectorRef& q, ConstVectorRef& qd) = 0;
};

}

#endif // CONTROLLER_H
