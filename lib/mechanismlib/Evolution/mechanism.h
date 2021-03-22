#ifndef MECHANISM_H
#define MECHANISM_H

#include "../common.h"
#include <Eigen/Dense>
#include "dna.h"
#include "objective.h"
#include "../Dynamics/connection.h"
#include "../Solver/solver.h"
#include "../Dynamics/constraintconnection.h"
#include "../Dynamics/dynamicelement.h"
#include "../Dynamics/motorelement.h"
#include "../Dynamics/link.h"
#include "../Control/controller.h"

namespace Mech {
// Forward declare state equation
class StateEquation;

class Mechanism
{
public:
    Mechanism(const DNA dna);
    // Direct, raw constructor that auto-initializes dna for convenience
    Mechanism(Eigen::MatrixXi incidenceMatrix,
              Eigen::VectorXi edge_labels,
              std::vector<Eigen::RowVectorXd> masses,
              std::vector<Eigen::RowVectorXd> parameters);
    ~Mechanism();

    void simulate(double t, int steps);
    void evaluate(Objective& objective, double t, int steps);
    void assignControllers(std::vector<std::shared_ptr<Controller>> controller_list);

    ConstRowMatrixRef statesTime() const;

    // Delete copy operators, mechanisms should only ever be explicitly reinitialized
    Mechanism & operator=(const Mechanism&) = delete;
    Mechanism(const Mechanism&) = delete;

    DNA dna() const;
    // The element lengths are calculated dynamically, make it accessible for animation and plotting purposes
    std::vector<Eigen::MatrixXd> momentarmList() const;

    Eigen::RowVectorXd initialState() const;

    // StateEquation is a friend class that handles the mathematics of calculating the dynamics
    friend class StateEquation;

    bool isSimulated() const;
    bool isEvaluated() const;
    double fitness() const;
    double time() const;
    int nrOfStates() const;
    int nrOfMasses() const;
    int nrOfConnections() const;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    // The DNA of the mechanism, which fully determines the dynamical structure. Does not change during a mechanisms lifetime
    const DNA dna_;

    // The result of the simulation
    RowMatrixXd states_time_;

    // Constructor functions for the mechanism structure
    void initializeLinks(const DNA& dna);
    void initializeMomentarmList(const DNA& dna);
    void initializeConnections(const DNA& dna);

    // The dynamic model of the mechanism that is created for simulation on creation of a mechanism object
    std::unique_ptr<StateEquation> state_equation_;
    std::unique_ptr<Solver> solver_;

    // The elements of the mechanism
    std::vector<std::unique_ptr<ConstraintConnection>> constraint_connections_;
    std::vector<std::unique_ptr<DynamicElement>> dynamic_elements_;
    std::vector<std::unique_ptr<MotorElement>> motor_elements_;
    std::vector<std::unique_ptr<Link>> links_;

    // The controllers
    std::vector<std::shared_ptr<Controller>> controllers_;

    //Structure parameters of the mechanism
    std::vector<Eigen::MatrixXd> moment_arm_list_;
    Eigen::RowVectorXd initial_state_;

    bool valid_mechanism_ = 0;
    bool simulated_ = 0;
    bool evaluated_ = 0;

    // The fitness value of the mechanism in the evolutionary race
    double fitness_ = 0;
    double time_ = 0;

    int nr_of_connections_;
    int nr_of_masses_;
    int nr_of_states_;

    int nr_of_constraints_ = 0;
    int nr_of_constraint_connections_ = 0;
    int nr_of_dynamic_elements_ = 0;
    int nr_of_motor_elements_ = 0;
};
}

#endif // MECHANISM_H
