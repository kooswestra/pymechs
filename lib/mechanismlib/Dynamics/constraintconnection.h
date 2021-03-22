#ifndef DYNAMICCONNECTION_H
#define DYNAMICCONNECTION_H

#include "../common.h"
#include <Eigen/Dense>
#include "connection.h"

/** @brief ConstraintConnection is an interface type for constraint inducing connection elements.
 *
 * A constraint inducing element has to have the constraint equations, the jacobian of the
 * constraint equation w.r.t. the state of the elements that are connected, and the convective
 * terms defined explicitly in it's functions.
 */

namespace Mech {

// ConstraintConnection is an interface type for constraint inducing connection elements
class ConstraintConnection : public Connection
{
public:
    using Connection::Connection;

    /** constraint() implements the calculation of the constraints, this is done to
     * determine numerical errors during the simulation which are consequently corrected.
     * */
    virtual void constraint(VectorRef D, ConstVectorRef& q, ConstMatrixRef& sin_cos_list) = 0;

    /** constraintD1() implements the calculation of the inter-jacobian
     *  \f$J_i\left(\mathbf{D_i}\right)=\frac{\partial\mathbf{D_i}}{\partial\mathbf{x_i}}\f$
     *  defined by a ConstraintConnection such as a hinge joint
     *
     * */
    virtual void constraintD1(MatrixRef Dk, ConstVectorRef& q, ConstMatrixRef& sin_cos_list) = 0;

    /** constraintD2() implements the calculation of the convective accelerations
     *  determined by the constraints and the velocities of the elements.
     * */
    virtual void constraintD2(VectorRef Dkk, ConstVectorRef& q, ConstVectorRef& qd, ConstMatrixRef& sin_cos_list) = 0;

    /** initializeJacobian() implements the initialization of the static elements of the inter-jacobian
     *  defined by this element, i.e. those independent of the state (usually 0 and 1 values)
     * */
    virtual void initializeJacobian(MatrixRef Dk) = 0;
private:
};

}

#endif // DYNAMICCONNECTION_H
