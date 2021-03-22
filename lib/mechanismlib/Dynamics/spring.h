#ifndef SPRING_H
#define SPRING_H

#include "../common.h"
#include <Eigen/Dense>
#include "dynamicelement.h"

namespace Mech {

template <bool type>
class Spring : public DynamicElement
{
public:
    Spring(Eigen::Matrix<double,1,6> K, Eigen::Matrix2d moment_arms, Eigen::Vector2i elementConnections);

    void applyForce(VectorRef force_vector, ConstVectorRef& state, ConstMatrixRef& sin_cos_list) override;

    int edgeLabel() const override {return 2;}

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    double k_;
    double L0_;
    Eigen::Matrix<double,1,6> K_;
    Eigen::Vector2d difference_vector_;

    // These are initialized on construction
    Eigen::Matrix2d moment_arms_;
};

}

#endif // SPRING_H
