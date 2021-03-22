#include "solver.h"
#include <../Dynamics/stateequation.h>

using namespace Eigen;
using namespace Mech;

Solver::Solver(double precision, int iter_limit)
    :   precision_(precision),
        iter_limit_(iter_limit)
{

}

bool Solver::integration_success() const
{
    return integration_success_;
}

// Newton's algorithm to make sure constraints are satisfied within precision at every step
void Solver::solveConstraints(StateEquation &state_equation, Ref<RowVectorXd> r_n) {

    // precision is the l2 norm of the maximum allowed error
    double error = state_equation.constraints().norm();

    // If the error is too big
    if(error > precision_) {
        // r[n+1] = r[n] - pinv(J)*f(r[n])
        state_equation.updateSinCosList(r_n);
        state_equation.updateJacobian(r_n);

        state_update_ = state_equation.jacobianTranspose()*(state_equation.jacobian()*state_equation.jacobianTranspose()).selfadjointView<Eigen::Lower>().llt().solve(state_equation.constraints());
        r_n.segment(0, state_equation.dim()/2) -= state_update_;
    }
}

void Solver::solveVelocities(StateEquation &state_equation, Ref<RowVectorXd> r_n) {

    state_equation.updateSinCosList(r_n);
    state_equation.updateJacobian(r_n);

    VectorXd error = state_equation.jacobian()*(r_n.segment(state_equation.dim()/2, state_equation.dim()/2)).transpose();

    r_n.segment(state_equation.dim()/2, state_equation.dim()/2) -= state_equation.jacobian().fullPivLu().solve(error);
}
