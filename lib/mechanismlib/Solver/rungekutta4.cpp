#include "rungekutta4.h"
#include "stateequation.h"

using namespace Eigen;
using namespace Mech;

ConstRowMatrixRef RungeKutta4::solve(StateEquation &state_equation,
                                 RowVectorXd x0,
                                 double t_end,
                                 int steps)
{
    double t, dt;
    int dim = state_equation.dim();

    RowVectorXd k1(dim);
    RowVectorXd k2(dim);
    RowVectorXd k3(dim);
    RowVectorXd k4(dim);

    // Initialize matrix and set initial conditions
    dt = t_end/steps;
    t = 0;
    r_ = RowMatrixXd(steps,dim);
    r_.row(0) = x0;

    // Running the loop to obtain an integrated result
    for(int i = 0; i < steps-1; i++) {

        /* Use class method of StateEquation to evaluate, StateEquation
         * is a function in the form dx/dt = f(x,t)
         *
         * Updated version that uses the Eigen library, less headache more speed */

        // k1 = f(t,y_t)
        state_equation.evaluate(t, r_.row(i), k1);

        // k2 = f(t + dt/2, y_t + dt*k1/2)
        state_equation.evaluate(t + dt/2, r_.row(i) + dt*k1/2, k2);

        // k3 = f(t + dt/2, y_t + dt*k2/2)
        state_equation.evaluate(t + dt/2, r_.row(i) + dt*k2/2, k3);

        // k4 = f(t+dt, y_t + dt*k3)
        state_equation.evaluate(t + dt, r_.row(i) + dt*k3, k4);

        /* Store result in r
         * r[i+1] = r[i] + dt*(k1 + 2*k2 + 2*k3 + k4)/6 */
        r_.row(i+1) = r_.row(i) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;

        // Solve D(r[i+1]) <= precision
        solveConstraints(state_equation, r_.row(i+1));

        // Check for NaN -> this indicates integration has become unstable
        if(r_.row(i+1).hasNaN()) {
            integration_success_ = false;
            return r_;
        }
        // Update timestep for next loop iteration
        t = t + dt;
    }

    integration_success_ = true;
    return r_;
}
