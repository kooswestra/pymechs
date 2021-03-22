#include "motor.h"
#include "pid.h"

using namespace Eigen;
using namespace Mech;

template<>
Motor<Grounded>::Motor(Vector2i element_connections)
    :   MotorElement(element_connections, Grounded)
{
    /** @todo Build a more efficient fix */
    // assign an empty controller on construction to prevent memory issues if the user forgot to define a controller
    controller_.reset(new Pid(10, Vector4d::Zero()));
    is_ground_connection_ = true;
}

template<>
void Motor<Grounded>::applyForce(VectorRef force_vector, double t, ConstVectorRef& q, ConstVectorRef& qd) {

    force_vector(indx2_ + 2) += controller_->controlOutput(t, q, qd);
}

template<>
Motor<Floating>::Motor(Vector2i element_connections)
    :   MotorElement(element_connections, Floating)
{
    /** @todo Build a more efficient fix */
    // assign an empty controller on construction to prevent memory issues if the user forgot to define a controller
    controller_.reset(new Pid(10, Vector4d::Zero()));
    is_ground_connection_ = false;
}

template<>
void Motor<Floating>::applyForce(VectorRef force_vector, double t, ConstVectorRef& q, ConstVectorRef& qd) {

    double torque = controller_->controlOutput(t, q, qd);
    force_vector(indx1_ + 2) -= torque;
    force_vector(indx2_ + 2) += torque;
}
