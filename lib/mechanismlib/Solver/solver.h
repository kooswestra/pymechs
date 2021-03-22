#ifndef SOLVER_H
#define SOLVER_H

#include "../common.h"
#include <Eigen/Dense>

namespace Mech {

class StateEquation;

/** @brief Solver is the abstract class interface for integration methods
  */
class Solver
{
public:
    Solver(double precision, int iter_limit);
    virtual ConstRowMatrixRef solve(StateEquation &state_equation, Eigen::RowVectorXd x0, double t_end, int steps) = 0;
    bool integration_success() const;

protected:
    double precision_;
    int iter_limit_;
    bool integration_success_ = false;

    RowMatrixXd r_;
    Eigen::RowVectorXd state_update_;
    Eigen::RowVectorXd velocity_update_;

    // Constraint solvers
    void solveConstraints(StateEquation& state_equation, RowVectorRef r_n);
    void solveVelocities(StateEquation& state_equation, RowVectorRef r_n);
};

}

#endif // SOLVER_H
