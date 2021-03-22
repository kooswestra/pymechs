#include "trapezoid.h"
#include "../Dynamics/stateequation.h"

using namespace Eigen;
using namespace Mech;

ConstRowMatrixRef Trapezoid::solve(StateEquation &state_equation, RowVectorXd x0, double t_end, int steps) {
    double t, dt;
    int dim = state_equation.dim();
    k1_ = RowVectorXd::Zero(dim);
    k2_ = RowVectorXd::Zero(dim);

    // Initialize matrix and set initial conditions
    dt = t_end/steps;
    t = 0;
    r_ = RowMatrixXd(steps,dim);
    r_.row(0) = x0;
    RowVectorXd xn_1(dim);

    // Running the loop to obtain an integrated result
    for(int i = 0; i < steps-1; i++) {

        /* Use class method of StateEquation to evaluate, StateEquation
         * is a function in the form dx/dt = f(x,t)
         *
         * Updated version that uses the Eigen library, less headache more speed */

        // k1 = f(t,y_t)
        state_equation.evaluate(t, r_.row(i), k1_);
        xn_1 = r_.row(i) + dt*k1_;

        // Solve the constraints at xn+1, using the euler step as first guess
        solveConstraints(state_equation, xn_1);
        // Also correct the velocity, using the corrected state value
        solveVelocities(state_equation, xn_1);

        // Evaluate xdd at the estimated and corrected xn+1
        state_equation.evaluate(t, xn_1, k2_);

        // Take the average of forward and backward integration to get vn+1 more exactly
        r_.row(i+1) = (r_.row(i) + dt*k2_ + xn_1)/2;

        // Use the average velocity get a better estimate of the new state

        // Constraint check again
        solveConstraints(state_equation, r_.row(i+1));

        if(integration_success_ == false) {break;}
        // Update timestep
        t = t + dt;
    }

    // Return the timeseries
    return r_;
}
