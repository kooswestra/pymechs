#include "pickandplace.h"
#include "mechanism.h"

using namespace Eigen;
using namespace Mech;

PickAndPlace::PickAndPlace(RowVector2d pick_location, RowVector2d place_location)
    :   pick_location_(pick_location),
        place_location_(place_location)
{
    l0_ = 2*(pick_location - place_location).norm();
}

VectorXd PickAndPlace::evaluate(Mechanism* mechanism)
{
        VectorXd fitness = VectorXd::Zero(4);
        double distance = 0;

        if(!mechanism->isSimulated()) {
            throw std::invalid_argument("Mechanism has to be successfully simulated to be evaluated");
        }

        if(!mechanism->dna().withEndEffector()) {
            throw std::invalid_argument("A mechanism without end-effector cannot be evaluated by this objective");
        }

        // End effector is always the last two states
        int nr_of_states = mechanism->nrOfStates();

        MatrixXd q = mechanism->statesTime().middleCols<2>(nr_of_states - 2);
        MatrixXd qd = mechanism->statesTime().rightCols(2);

        // Time dependency makes life easier, be at x at time y
        int step = q.rows()/2;

        fitness(1) += (q.row(0) - pick_location_).squaredNorm() + 0.2*qd.row(0).squaredNorm();
        fitness(1) += (q.row(step - 1) - place_location_).squaredNorm() + 0.2*qd.row(step - 1).squaredNorm();
        fitness(1) += (q.row(2*step - 1) - pick_location_).squaredNorm() + 0.2*qd.row(2*step - 1).squaredNorm();

        // Total traversed path minimization
        for(int i = 1; i < q.rows(); i++) {
            distance += (q.row(i) - q.row(i-1)).norm();
        }

        fitness(2) += abs(distance - l0_);

        // Complexity score
        fitness(3) += mechanism->dna().complexity();

        // Weighted total score (was 0.03)
        fitness(0) = - pow(fitness(1),1) - pow(0.12 * fitness(2), 2) - pow(0.03 * fitness(3), 2);

        return fitness;

}
