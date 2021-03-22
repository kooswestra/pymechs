#ifndef COMMON_H
#define COMMON_H

#include <Eigen/Core>

namespace Mech {

/** Templated specializations are used for floating and grounded connections*/
const bool Grounded = true;
const bool Floating = false;

/** Convenience typedefs for some often used Eigen types */
typedef const Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> ConstRowMatrixRef;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> RowMatrixXd;

typedef const Eigen::Ref<const Eigen::MatrixXd> ConstMatrixRef;
typedef Eigen::Ref<Eigen::MatrixXd> MatrixRef;

typedef const Eigen::Ref<const Eigen::RowVectorXd> ConstRowVectorRef;
typedef Eigen::Ref<Eigen::RowVectorXd> RowVectorRef;

typedef const Eigen::Ref<const Eigen::VectorXd> ConstVectorRef;
typedef Eigen::Ref<Eigen::VectorXd> VectorRef;

}

#endif // COMMON_H