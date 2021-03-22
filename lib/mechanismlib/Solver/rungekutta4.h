#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H

#include <Eigen/Dense>
#include "../common.h"
#include "solver.h"

namespace Mech {

class RungeKutta4 : public Solver
{
public:
    using Solver::Solver;
    ConstRowMatrixRef solve(StateEquation &state_equation,
                        Eigen::RowVectorXd x0,
                        double t_end,
                        int steps) override;
};

}

#endif // RUNGEKUTTA4_H
