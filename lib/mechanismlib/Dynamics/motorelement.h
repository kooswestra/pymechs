#ifndef MOTORELEMENT_H
#define MOTORELEMENT_H

#include "../common.h"
#include <memory>
#include <Eigen/Dense>
#include "connection.h"
#include "../Control/controller.h"

/** @brief MotorElement is an interface for active force inducing elements
 *
 * A motor element applies a force that is passed to it from an external source.
 * Note that this actively powers the mechanism and as such adds or removes energy. It still
 * requires the state as input in order to determine moment arms as the
 * force does not have to be applied at the c.o.m. of the connected elements
 */

namespace Mech {

class MotorElement : public Connection
{
public:
    using Connection::Connection;
    virtual void applyForce(VectorRef force_vector, double t, ConstVectorRef& q, ConstVectorRef& qd) = 0;

    void assignController(std::shared_ptr<Controller> controller);
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
    std::shared_ptr<Controller> controller_;
};

}

#endif // MOTORELEMENT_H
