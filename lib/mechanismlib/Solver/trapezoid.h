#ifndef TRAPEZOID_H
#define TRAPEZOID_H

#include <Eigen/Dense>
#include "solver.h"

namespace Mech {

class Trapezoid  : public Solver
{
public:
    using Solver::Solver;
    const Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> solve(StateEquation &state_equation,
                                                                                               Eigen::RowVectorXd x0,
                                                                                               double t_end,
                                                                                               int steps) override;

    Eigen::RowVectorXd k1_, k2_;
};

}

#endif // TRAPEZOID_H
