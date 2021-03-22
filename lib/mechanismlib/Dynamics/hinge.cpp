#include "hinge.h"

using namespace Eigen;
using namespace Mech;

// Mass - ground connections ------------------------------------------------------------

template<>
Hinge<Grounded>::Hinge(RowVector2d H, Matrix2d com_list, Vector2i element_connections)
    :   ConstraintConnection(element_connections, Grounded),
        H_(H)
{
    com_list_ = com_list;
    l_.row(0) = com_list.row(0) - H;

    is_ground_connection_ = true;
}

template<>
void Hinge<Grounded>::constraint(VectorRef D, ConstVectorRef& q, ConstMatrixRef& sin_cos_list)
{
    double st = sin_cos_list(indx2_/3, 0);
    double ct = sin_cos_list(indx2_/3, 1);

    D(0) = q(indx2_ + 0) - H_(0) - l_(0,0)*ct + l_(0,1)*st;
    D(1) = q(indx2_ + 1) - H_(1) - l_(0,0)*st - l_(0,1)*ct;
}

template<>
void Hinge<Grounded>::constraintD1(MatrixRef Dk, ConstVectorRef& q, ConstMatrixRef& sin_cos_list)
{
    double st = sin_cos_list(indx2_/3, 0);
    double ct = sin_cos_list(indx2_/3, 1);

    Dk(0, indx2_ + 2) = l_(0,0)*st + l_(0,1)*ct;
    Dk(1, indx2_ + 2) = -l_(0,0)*ct + l_(0,1)*st;
}

template<>
void Hinge<Grounded>::constraintD2(VectorRef Dkk, ConstVectorRef& q, ConstVectorRef& qd, ConstMatrixRef& sin_cos_list)
{
    double st = sin_cos_list(indx2_/3, 0);
    double ct = sin_cos_list(indx2_/3, 1);

    Dkk(0) = qd(indx2_ + 2)*qd(indx2_ + 2)*(l_(0,0)*ct - l_(0,1)*st);
    Dkk(1) = qd(indx2_ + 2)*qd(indx2_ + 2)*(l_(0,0)*st + l_(0,1)*ct);
}

template<>
void Hinge<Grounded>::initializeJacobian(MatrixRef Dk)
{
        Dk(0,indx2_ + 0) = 1;
        Dk(1,indx2_ + 1) = 1;
}


// Mass - Mass connections -----------------------------------------------------------------

template<>
Hinge<Floating>::Hinge(RowVector2d H, Matrix2d com_list, Vector2i element_connections)
    :   ConstraintConnection(element_connections, Floating),
        H_(H)
{
    com_list_ = com_list;
    l_.row(0) = com_list.row(0) - H;
    l_.row(1) = com_list.row(1) - H;

    is_ground_connection_ = false;
}

template<>
void Hinge<Floating>::constraint(VectorRef D, ConstVectorRef& q, ConstMatrixRef& sin_cos_list)
{
    double st1 = sin_cos_list(indx1_/3, 0);
    double st2 = sin_cos_list(indx2_/3, 0);
    double ct1 = sin_cos_list(indx1_/3, 1);
    double ct2 = sin_cos_list(indx2_/3, 1);

    D(0) = q(indx2_ + 0) - q(indx1_ + 0) + l_(0,0)*ct1 - l_(0,1)*st1 -
            l_(1,0)*ct2 + l_(1,1)*st2;

    D(1) = q(indx2_ + 1) - q(indx1_ + 1) + l_(0,0)*st1 + l_(0,1)*ct1 -
            l_(1,0)*st2 - l_(1,1)*ct2;
}

template<>
void Hinge<Floating>::constraintD1(MatrixRef Dk, ConstVectorRef& q, ConstMatrixRef& sin_cos_list)
{
    double st1 = sin_cos_list(indx1_/3, 0);
    double st2 = sin_cos_list(indx2_/3, 0);
    double ct1 = sin_cos_list(indx1_/3, 1);
    double ct2 = sin_cos_list(indx2_/3, 1);

    Dk(0,indx1_ + 2) = -l_(0,0)*st1 - l_(0,1)*ct1;
    Dk(1,indx1_ + 2) = l_(0,0)*ct1 - l_(0,1)*st1;
    Dk(0,indx2_ + 2) = l_(1,0)*st2 + l_(1,1)*ct2;
    Dk(1,indx2_ + 2) = -l_(1,0)*ct2 + l_(1,1)*st2;
}

template<>
void Hinge<Floating>::constraintD2(VectorRef Dkk, ConstVectorRef& q, ConstVectorRef& qd, ConstMatrixRef& sin_cos_list)
{
    double st1 = sin_cos_list(indx1_/3, 0);
    double st2 = sin_cos_list(indx2_/3, 0);
    double ct1 = sin_cos_list(indx1_/3, 1);
    double ct2 = sin_cos_list(indx2_/3, 1);

    Dkk(0) = qd(indx1_ + 2)*qd(indx1_ + 2)*(-l_(0,0)*ct1 + l_(0,1)*st1) +
            qd(indx2_ + 2)*qd(indx2_ + 2)*(l_(1,0)*ct2 - l_(1,1)*st2);

    Dkk(1) = qd(indx1_ + 2)*qd(indx1_ + 2)*(-l_(0,0)*st1 - l_(0,1)*ct1) +
            qd(indx2_ + 2)*qd(indx2_ + 2)*(l_(1,0)*st2 + l_(1,1)*ct2);
}

template<>
void Hinge<Floating>::initializeJacobian(MatrixRef Dk)
{
    Dk(0,indx1_ + 0) = -1;
    Dk(1,indx1_ + 1) = -1;
    Dk(0,indx2_ + 0) = 1;
    Dk(1,indx2_ + 1) = 1;
}
