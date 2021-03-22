#include "mechanism.h"
#include "../Dynamics/hinge.h"
#include "../Dynamics/endeffector.h"
#include "../Dynamics/motor.h"
#include "../Dynamics/spring.h"
#include "../Dynamics/torsionspring.h"
#include "../Dynamics/stateequation.h"
#include "../Control/pid.h"
#include "../Solver/rungekutta4.h"

using namespace Eigen;
using namespace Mech;

Mechanism::Mechanism(const DNA dna)
    :   dna_(dna)
{
    if(!dna.isConnected())
            throw std::invalid_argument("Trying to construct a mech with disconnected elements");
    if(!dna.isGrounded())
            throw std::invalid_argument("Trying to construct a mech that isn't grounded");
    if(!dna.isDynamic())
            throw std::invalid_argument("Trying to construct a mech that isn't dynamic");

    // The number of connections is equal to the number of columns in the incidence matrix
    nr_of_connections_ = dna.incidenceMatrix().cols();

    // Calculate number of states, links have three (x,y,θ) while point-masses have two (x,y)
    if(dna.withEndEffector()) {
        nr_of_masses_ = dna.incidenceMatrix().rows() - 2;
        nr_of_states_ = nr_of_masses_*3 + 2;
    }
    else {
        nr_of_masses_ = dna.incidenceMatrix().rows() - 1;
        nr_of_states_ = nr_of_masses_*3;
    }

    // Build the list of mass elements
    initializeLinks(dna);

    initial_state_ = RowVectorXd::Zero(nr_of_states_);
    // The initial state vector is the initial centre of mass of the links
    for(int i = 0; i < links_.size(); i++) {
        initial_state_.segment<2>(3*i) = links_[i]->centreOfMass();
    }

    // Build the moment arm list using the centre of mass list
    initializeMomentarmList(dna);
    // And then the connection structure of the mechanism
    initializeConnections(dna);

    // If all the tests passed, build the state equation and attach a solver
    state_equation_.reset(new StateEquation(*this));
    solver_.reset(new RungeKutta4(1e-12, 5));
}

Mechanism::Mechanism(MatrixXi incidence_matrix,
                     VectorXi edge_labels,
                     std::vector<RowVectorXd> masses,
                     std::vector<RowVectorXd> parameters)
    :   Mechanism(DNA(incidence_matrix, edge_labels, masses, parameters))
{
}

Mechanism::~Mechanism() 
{
}

void Mechanism::simulate(double t, int steps)
{
    // The initial conditions is the initial state, with velocity of zero
    RowVectorXd initial_conditions = RowVectorXd::Zero(nr_of_states_*2);
    initial_conditions.head(nr_of_states_) = initialState();

    states_time_ = solver_->solve(*state_equation_, initial_conditions, t, steps);
    simulated_ = solver_->integration_success();

    time_ = t;
}

void Mechanism::evaluate(Objective &objective, double t, int steps)
{
    if(evaluated_)
        return;

    time_ = t;
    simulate(t, steps);

    if(simulated_)
        fitness_ = objective.evaluate(this)(0);
    else
        fitness_ = std::numeric_limits<double>::lowest();

    evaluated_ = true;
}

void Mechanism::assignControllers(std::vector<std::shared_ptr<Controller>> controller_list)
{
   if(controller_list.size() != nr_of_motor_elements_)  {
       throw std::invalid_argument("The new controller list must match the amount of motors in the mechanism");
   }

   for(size_t i = 0; i < controller_list.size(); i++)
   {
       motor_elements_[i]->assignController(controller_list[i]);
   }
}

const Ref<const Matrix<double,Dynamic,Dynamic,RowMajor>> Mechanism::statesTime() const
{
    return states_time_;
}

DNA Mechanism::dna() const
{
    return dna_;
}

void Mechanism::initializeLinks(const DNA& dna) {

    // For each mass
    for(int i = 0; i < nr_of_masses_; i++) {

        // What connections are associated with the mass?
        VectorXi element_connections = dna.linkConnections(i+1);
        MatrixXd connection_locations = Matrix<double, Dynamic, 2>::Zero(element_connections.size(), 2);

        // For each connection of the mass
        for(int j = 0; j < element_connections.size(); j++) {
            // Determine connection positions, these are (almost) always the first 2 parameters, except for linear springs
            // These are a special case, as the connection point to element 1 is not equal to that of 2, there is no constraint
            if(dna.edgeLabels(element_connections[j]) == DNA::SPRING) {
                Vector2i massnrs = dna.edgeConnections(element_connections[j]);
                // Check whether this mass is the first or second connected element by the spring
                if(i == massnrs[0]-1) {
                    connection_locations.row(j) = dna.parameters(element_connections[j]).segment<2>(0);
                }
                else {
                    connection_locations.row(j) = dna.parameters(element_connections[j]).segment<2>(2);
                }
            }
            else
            {
                connection_locations.row(j) = dna.parameters(element_connections[j]).segment<2>(0);
            }
        }

        links_.emplace_back(new Link(dna.masses(i), connection_locations, element_connections, i));
    }
}

void Mechanism::initializeMomentarmList(const DNA& dna) {

    for(int i = 0; i < nr_of_connections_; i++) {

        VectorXi index = dna.edgeConnections(i);
        index[0] -= 1;
        index[1] -= 1;

        Matrix2d moment_arms;
        // Check if it is a ground connection, which only has one nonzero moment arm. Also check if the
        // connection is a spring, as a spring defines no constraint and thus has two different connection points
        if(index[0] == -1) {
            if(dna.edgeLabels(i) == DNA::SPRING) {
                moment_arms.col(0) = links_[index[1]]->centreOfMass() - dna.parameters(i).segment<2>(2);
            }
            else {
                moment_arms.col(0) = links_[index[1]]->centreOfMass() - dna.parameters(i).segment<2>(0);
            }
            moment_arm_list_.push_back(moment_arms);
        }
        else {
            if(dna.edgeLabels(i) == DNA::SPRING) {
                moment_arms.col(0) = links_[index[0]]->centreOfMass() - dna.parameters(i).segment<2>(0);
                moment_arms.col(1) = links_[index[1]]->centreOfMass() - dna.parameters(i).segment<2>(2);
            }
            else if(dna.edgeLabels(i) == DNA::END_EFFECTOR) {
                moment_arms.col(0) = links_[index[0]]->centreOfMass() - dna.parameters(i).segment<2>(0);
            }
            else {
                moment_arms.col(0) = links_[index[0]]->centreOfMass() - dna.parameters(i).segment<2>(0);
                moment_arms.col(1) = links_[index[1]]->centreOfMass() - dna.parameters(i).segment<2>(0);
            }
            moment_arm_list_.push_back(moment_arms);
        }
    }

}

void Mechanism::initializeConnections(const DNA& dna)
{
    // Build the connections and elements present in the DNA
    for(int i = 0; i < nr_of_connections_; i++) {

        // Read the incidence matrix column related to this connection &
        // Interpret the incidence matrix, determine connections
        Vector2i element_connections = dna.edgeConnections(i);

        // Check if the element also defines a constraint, and which one. Then add it to the constraint connection list
        if(dna.edgeLabels(i) == DNA::HINGE || dna.edgeLabels(i) == DNA::TORSION_SPRING || dna.edgeLabels(i) == DNA::HINGE_MOTOR )
        {
            nr_of_constraint_connections_++;
            // Hinges add two constraints
            nr_of_constraints_ = nr_of_constraints_+2;
            // Grab the parameters related to the connection
            RowVector2d H = dna.parameters(i).segment<2>(0);
            if((element_connections.array() == 0).any())
            {
                Matrix2d com;
                com.row(0) = links_[element_connections[1]-1]->centreOfMass();
                constraint_connections_.emplace_back(new Hinge<Grounded>(H, com, element_connections));
            }
            else
            {
                Matrix2d com;
                com << links_[element_connections[0]-1]->centreOfMass(), links_[element_connections[1]-1]->centreOfMass();
                constraint_connections_.emplace_back(new Hinge<Floating>(H, com, element_connections));
            }
        }

        if(dna.edgeLabels(i) == DNA::END_EFFECTOR)
        {
            nr_of_constraint_connections_++;
            nr_of_constraints_ = nr_of_constraints_+2;
            RowVector2d H = dna.parameters(i).segment<2>(0);
            constraint_connections_.emplace_back(new EndEffector(H, links_[element_connections[0]-1]->centreOfMass(), element_connections));
            initial_state_.tail(2) = H;
        }

        // Update the dynamic and motor element lists
        switch(dna.edgeLabels(i)) {
        case DNA::HINGE :
        {
            break;
        }
        case DNA::END_EFFECTOR :
        {
            break;
        }
        case DNA::SPRING :
        {
            nr_of_dynamic_elements_++;
            if((element_connections.array() == 0).any()) {
                dynamic_elements_.emplace_back(new Spring<Grounded>(dna.parameters(i), moment_arm_list_[i], element_connections));
            }
            else {
                dynamic_elements_.emplace_back(new Spring<Floating>(dna.parameters(i), moment_arm_list_[i], element_connections));
            }
            break;
        }
        case DNA::TORSION_SPRING :
        {
            nr_of_dynamic_elements_++;
            if((element_connections.array() == 0).any()) {
                dynamic_elements_.emplace_back(new TorsionSpring<Grounded>(dna.parameters(i), element_connections));
            }
            else {
                dynamic_elements_.emplace_back(new TorsionSpring<Floating>(dna.parameters(i), element_connections));
            }
            break;
        }
        case DNA::HINGE_MOTOR :
        {
            nr_of_motor_elements_++;
            if((element_connections.array() == 0).any()) {
                motor_elements_.emplace_back(new Motor<Grounded>(element_connections));
            }
            else {
                motor_elements_.emplace_back(new Motor<Floating>(element_connections));
            }
            break;
        }
        default: {
            throw std::invalid_argument("Invalid or not yet implemented edge label provided to mechanism constructor");
        }
        }
    }
}

int Mechanism::nrOfConnections() const
{
    return nr_of_connections_;
}

int Mechanism::nrOfMasses() const
{
    return nr_of_masses_;
}

int Mechanism::nrOfStates() const
{
    return nr_of_states_;
}

RowVectorXd Mechanism::initialState() const
{
    return initial_state_;
}

std::vector<MatrixXd> Mechanism::momentarmList() const
{
    return moment_arm_list_;
}

bool Mechanism::isSimulated() const
{
    return simulated_;
}

bool Mechanism::isEvaluated() const
{
    return evaluated_;
}

double Mechanism::fitness() const
{
    return fitness_;
}

double Mechanism::time() const
{
    return time_;
}
