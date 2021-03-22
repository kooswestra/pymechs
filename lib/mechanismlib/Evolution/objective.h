#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <Eigen/Dense>

/** @brief Objective is the interface for fitness evaluation of the mechanism
 *
 * It defines the objective function and provides functions to evaluate the mechanisms for their
 * performance on that objective function. This includes complexity penalties on the DNA.
 */

namespace Mech {
// Forward declare mechanism class
class Mechanism;

class Objective
{
public:
    // Need to think about the structure here, reference trajectory defines points in order that have to be passed
    Objective();
    virtual ~Objective();
    virtual Eigen::VectorXd evaluate(Mechanism* mechanism)  = 0;

};
}

#endif // OBJECTIVE_H
