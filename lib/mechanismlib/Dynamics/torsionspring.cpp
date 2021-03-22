#include "torsionspring.h"

using namespace Eigen;
using namespace Mech;

template<>
TorsionSpring<Grounded>::TorsionSpring(RowVector4d T, Vector2i element_connections)
    :   DynamicElement(element_connections, Grounded),
        t0_(T(2)),
        k_(abs(T(3)))
{
    is_ground_connection_ = true;
}

template<>
void TorsionSpring<Grounded>::applyForce(VectorRef force_vector, ConstVectorRef& state, ConstMatrixRef& /*sin_cos_list*/) {

    double torque = (t0_ - state(indx2_ + 2))*k_;
    force_vector(indx2_ + 2) += torque;
}

template<>
TorsionSpring<Floating>::TorsionSpring(RowVector4d T, Vector2i element_connections)
    :   DynamicElement(element_connections, Floating),
        t0_(T(2)),
        k_(abs(T(3)))
{
    is_ground_connection_ = false;
}

template<>
void TorsionSpring<Floating>::applyForce(VectorRef force_vector, ConstVectorRef& state, ConstMatrixRef& /*sin_cos_list*/) {

    double torque = (t0_ + state(indx1_ + 2) - state(indx2_ + 2))*k_;

    force_vector(indx1_ + 2) -= torque;
    force_vector(indx2_ + 2) += torque;
}
