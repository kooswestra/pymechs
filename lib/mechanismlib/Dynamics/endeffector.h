#ifndef ENDEFFECTOR_H
#define ENDEFFECTOR_H

#include "../common.h"
#include <Eigen/Dense>
#include "constraintconnection.h"

namespace Mech {

class EndEffector : public ConstraintConnection
{
public:
    EndEffector(Eigen::RowVector2d H, Eigen::RowVector2d com_list, Eigen::Vector2i element_Connections);

    // Constraint equations
    void constraint(VectorRef D, ConstVectorRef& q, ConstMatrixRef& sin_cos_list) override;
    void constraintD1(MatrixRef Dk, ConstVectorRef& q, ConstMatrixRef& sin_cos_list) override;                  
    void constraintD2(VectorRef Dkk, ConstVectorRef& q, ConstVectorRef& qd, ConstMatrixRef& sin_cos_list) override;

    int edgeLabel() const override {return 1;}
    void initializeJacobian(MatrixRef Dk) override;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Eigen::RowVector2d l_;
};

}


#endif // ENDEFFECTOR_H
