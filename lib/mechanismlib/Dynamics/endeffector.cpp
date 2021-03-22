#include "endeffector.h"

using namespace Eigen;
using namespace Mech;

EndEffector::EndEffector(RowVector2d H, RowVector2d com_list, Vector2i element_connections)
    :   ConstraintConnection(element_connections, Floating)
{
    l_ = com_list - H;
}

void EndEffector::constraint(VectorRef D, ConstVectorRef& q, ConstMatrixRef& sin_cos_list)
{
    double st1 = sin_cos_list(indx1_/3, 0);
    double ct1 = sin_cos_list(indx1_/3, 1);

    D(0) = q(indx2_ + 0) - q(indx1_ + 0) + l_(0)*ct1 - l_(1)*st1;
    D(1) = q(indx2_ + 1) - q(indx1_ + 1) + l_(0)*st1 + l_(1)*ct1;
}

void EndEffector::constraintD1(MatrixRef Dk, ConstVectorRef& q, ConstMatrixRef& sin_cos_list)
{
    double st1 = sin_cos_list(indx1_/3, 0);
    double ct1 = sin_cos_list(indx1_/3, 1);

    Dk(0,indx1_ + 2) = -l_(0)*st1 - l_(1)*ct1;
    Dk(1,indx1_ + 2) = l_(0)*ct1 - l_(1)*st1;
}

void EndEffector::constraintD2(VectorRef Dkk, ConstVectorRef& q, ConstVectorRef& qd, ConstMatrixRef& sin_cos_list)
{
    double st1 = sin_cos_list(indx1_/3, 0);
    double ct1 = sin_cos_list(indx1_/3, 1);

    Dkk(0) = qd(indx1_ + 2)*qd(indx1_ + 2)*(-l_(0)*ct1 + l_(1)*st1);
    Dkk(1) = qd(indx1_ + 2)*qd(indx1_ + 2)*(-l_(0)*st1 - l_(1)*ct1);
}

void EndEffector::initializeJacobian(MatrixRef Dk)
{
    Dk(0,indx1_ + 0) = -1;
    Dk(1,indx1_ + 1) = -1;
    Dk(0,indx2_ + 0) = 1;
    Dk(1,indx2_ + 1) = 1;
}
