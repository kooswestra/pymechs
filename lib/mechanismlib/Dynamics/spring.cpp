#include "spring.h"
#include <iostream>

using namespace Eigen;
using namespace Mech;

template<>
Spring<Grounded>::Spring(Matrix<double,1,6> K, Matrix2d moment_arms, Vector2i element_connections)
    :   DynamicElement(element_connections, Grounded),
        k_(abs(K(5))),
        L0_(abs(K(4))),
        K_(K)
{
    moment_arms_ = moment_arms;
    is_ground_connection_ = true;
}

template<>
void Spring<Grounded>::applyForce(VectorRef force_vector, ConstVectorRef& state, ConstMatrixRef& sin_cos_list)
{
    double st = sin_cos_list(indx2_/3, 0);
    double ct = sin_cos_list(indx2_/3, 1);

    double rx = moment_arms_(0,0)*ct - moment_arms_(1,0)*st;
    double ry = moment_arms_(0,0)*st + moment_arms_(1,0)*ct;

    difference_vector_(0) = K_(0) - (state(indx2_ + 0) - rx);
    difference_vector_(1) = K_(1) - (state(indx2_ + 1) - ry);

    // Normalize the difference vector to obtain direction of force
    double magnitude = (difference_vector_.norm() - L0_)*k_;

    force_vector.segment<2>(indx2_ + 0) += magnitude*difference_vector_;
    force_vector(indx2_ + 2) += magnitude*difference_vector_(0)*ry - magnitude*difference_vector_(1)*rx;
}

template<>
Spring<Floating>::Spring(Matrix<double,1,6> K, Matrix2d moment_arms, Vector2i element_connections)
    :   DynamicElement(element_connections, Floating),
        k_(K(5)),
        L0_(K(4)),
        K_(K)
{
    moment_arms_ = moment_arms;
    is_ground_connection_ = false;
}

template<>
void Spring<Floating>::applyForce(VectorRef force_vector, ConstVectorRef& state, ConstMatrixRef& sin_cos_list)
{
    double st1 = sin_cos_list(indx1_/3, 0);
    double st2 = sin_cos_list(indx2_/3, 0);
    double ct1 = sin_cos_list(indx1_/3, 1);
    double ct2 = sin_cos_list(indx2_/3, 1);

    double rx1 = moment_arms_(0,0)*ct1 - moment_arms_(1,0)*st1;
    double ry1 = moment_arms_(0,0)*st1 + moment_arms_(1,0)*ct1;
    double rx2 = moment_arms_(0,1)*ct2 - moment_arms_(1,1)*st2;
    double ry2 = moment_arms_(0,1)*st2 + moment_arms_(1,1)*ct2;

    // (X1 - R(th1)*Moment arm) - (X2 - R(th2)*Moment arm) = difference vector
    difference_vector_(0) = (state(indx1_ + 0) - rx1) - (state(indx2_ + 0) - rx2);
    difference_vector_(1) = (state(indx1_ + 1) - ry1) - (state(indx2_ + 1) - ry2);

    // Get the magnitude of the force exerted by the spring
    double magnitude = (difference_vector_.norm() - L0_)*k_;
    // Normalize the difference vector to obtain direction of force
    difference_vector_.normalize();

    force_vector.segment<2>(indx1_ + 0) -= magnitude*difference_vector_;
    force_vector(indx1_ + 2) -= magnitude*difference_vector_(0)*ry1 - magnitude*difference_vector_(1)*rx1;
    force_vector.segment<2>(indx2_ + 0) += magnitude*difference_vector_;
    force_vector(indx2_ + 2) += magnitude*difference_vector_(0)*ry2 - magnitude*difference_vector_(1)*rx2;
}
