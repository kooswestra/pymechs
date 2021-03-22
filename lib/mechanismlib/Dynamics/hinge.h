#ifndef MASSHINGE_H
#define MASSHINGE_H

#include "../common.h"
#include <Eigen/Dense>
#include "constraintconnection.h"

namespace Mech {

// Connection of a mass with a hinge to a ground or other mass
template <bool type>
class Hinge : public ConstraintConnection
{
public:
    Hinge(Eigen::RowVector2d H, Eigen::Matrix2d com_list, Eigen::Vector2i element_Connections);

    // Constraint equations
    void constraint(VectorRef D, ConstVectorRef& q, ConstMatrixRef& sin_cos_list) override;
    void constraintD1(MatrixRef Dk, ConstVectorRef& q, ConstMatrixRef& sin_cos_list) override;                  
    void constraintD2(VectorRef Dkk, ConstVectorRef& q, ConstVectorRef& qd, ConstMatrixRef& sin_cos_list) override;

    int edgeLabel() const override {return 0;}
    void initializeJacobian(MatrixRef Dk) override;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Eigen::RowVector2d H_;
    Eigen::Matrix2d com_list_;
    Eigen::Matrix2d l_;
};

}

#endif // MASSHINGE_H
