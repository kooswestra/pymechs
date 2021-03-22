#ifndef DYNAMICELEMENT_H
#define DYNAMICELEMENT_H

#include "../common.h"
#include <Eigen/Dense>
#include "connection.h"

/** @brief DynamicElement is an interface for passive force inducing elements such as springs
 *
 * A dynamic element has no active torque that it can apply and as such has a force output
 * that depends on the time and mechanism state only. Note that a dynamic element does
 * NOT actively add energy to the system. It can only store, release and dampen
 */
namespace Mech {

class DynamicElement : public Connection
{
public:
    using Connection::Connection;

    /** evaluateForce() calculates the forces applied by a dynamic element to the elements
    * that this dynamic element connects. It takes as input the time and the state(s) of the
    * connect(ed) link(s), this can be either a 3d vector (ground connection) or a 6d vector
    * (link-link connection). And correspondingly outputs either a 3d or 6d force vector.
    */
    virtual void applyForce(VectorRef force_vector, ConstVectorRef& state, ConstMatrixRef& sin_cos_list) = 0;
private:
};

}

#endif // DYNAMICELEMENT_H
