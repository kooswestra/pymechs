#include "stateequation.h"
#include "constraintconnection.h"
#include "dynamicelement.h"
#include "motorelement.h"

using namespace Eigen;
using namespace Mech;

StateEquation::StateEquation(Mechanism &mechanism)
    :   mch_(mechanism)
{
    // Initialize state and initial conditions vectors in the right (dynamic) size
    equation_vector_ = VectorXd::Zero(mch_.nr_of_states_);
    gravity_vector_ = VectorXd::Zero(mch_.nr_of_states_);
    force_vector_ = VectorXd::Zero(mch_.nr_of_states_);
    sin_cos_list_ = MatrixXd::Zero(mch_.nr_of_masses_, 2);

    // Initialize mass matrix, connections and centre of mass data from dna
    initializeMassMatrix();
    initializeGravity();
    initializeJacobian();

    // Initialize constraint matrices and force vector in the right size
    D_ = VectorXd::Zero(mch_.nr_of_constraints_);
    Dkk_ = VectorXd::Zero(mch_.nr_of_constraints_);
}

StateEquation::~StateEquation() 
{

}

int StateEquation::dim() const
{
    return mch_.nr_of_states_ * 2;
}

void StateEquation::evaluate(double t, ConstRowVectorRef &state, RowVectorRef stated)
{

    // First  precompute required sines and cosines and store in the list [[sin(th1), cos(th1)], [sin(th2), cos(th2)], ... ]
    updateSinCosList(state);

    // Update constraint equations
    updateConstraints(state);

    // Update the jacobian of the constraints
    updateJacobian(state);

    // Update the convective acceleration terms
    updateConvectiveAccelerations(state);

    // Update the force vector
    updateForces(t, state);

    // Simply copy the velocity
    stated.segment(0, mch_.nr_of_states_) = state.segment(mch_.nr_of_states_,  mch_.nr_of_states_);

    // Solve for the accelerations and constraints using the square, symmetric and invertible properties of the matrix
    Dk_T_ = Dk_.transpose().eval();
    equation_vector_ = Dk_*M_.asDiagonal()*force_vector_ + Dkk_;

    // Solve the matrix inversion to obtain the reaction forces
    (Dk_*M_.asDiagonal()*Dk_T_).selfadjointView<Eigen::Lower>().llt().solveInPlace(equation_vector_);

    // Solve the rest of the equation with the trivial known M^-1 to obtain the accelerations
    stated.segment(mch_.nr_of_states_, mch_.nr_of_states_) = M_.asDiagonal()*(force_vector_ - Dk_T_*equation_vector_);

}

void StateEquation::updateSinCosList(ConstRowVectorRef& state)
{
    for(int i = 0; i < mch_.nr_of_masses_; i++)
    {
        sin_cos_list_(i,0) = sin(state(3*i + 2));
        sin_cos_list_(i,1) = cos(state(3*i + 2));
    }
}

void StateEquation::updateJacobian(ConstRowVectorRef& state)
{

    // Iterate over the constraints as defined by the constraint inducing elements, each of which fully defines two rows of the jacobian matrix
    // -> Updates the jacobian of the complete system in place
    for(int i = 0; i < mch_.nr_of_constraint_connections_; i++)
    {
        // J(D) = f(q)
        mch_.constraint_connections_[i]->constraintD1(Dk_.block(i*2,0,2,Dk_.cols()),
                                                      state.segment(0,mch_.nr_of_states_),
                                                      sin_cos_list_);
    }
}

ConstMatrixRef StateEquation::jacobian() const {
    return Dk_;
}

ConstMatrixRef StateEquation::jacobianTranspose() const {
    return Dk_T_;
}


void StateEquation::updateConvectiveAccelerations(ConstRowVectorRef& state)
{
    // Iterate over the convective accelerations as defined by the constraint inducing elements
    // -> Updates the convective accelerations of the complete system in place
    for(int i = 0; i < mch_.nr_of_constraint_connections_; i++)
    {
        // Dkk = f(q,qd)
        mch_.constraint_connections_[i]->constraintD2(Dkk_.segment<2>(i*2),
                                                      state.segment(0,mch_.nr_of_states_),
                                                      state.segment(mch_.nr_of_states_, mch_.nr_of_states_),
                                                      sin_cos_list_);
    }
}

ConstVectorRef StateEquation::convectiveAccelerations() const
{
    return Dkk_;
}


void StateEquation::updateConstraints(ConstRowVectorRef& state)
{
    // Evaluate the constraint matrices for this specific state
    for(int i = 0; i < mch_.nr_of_constraint_connections_; i++)
    {
        mch_.constraint_connections_[i]->constraint(D_.segment<2>(i*2),
                                                    state.segment(0,mch_.nr_of_states_),
                                                    sin_cos_list_);
    }
}

ConstVectorRef StateEquation::constraints() const
{
    return D_;
}


void StateEquation::updateForces(double t, ConstRowVectorRef& state)
{

    // First apply gravity to the system
    force_vector_ = gravity_vector_;
    //force_vector_.setZero();

    // Compute forces induced by conservative elements (spring-like)
    for(int i = 0; i < mch_.nr_of_dynamic_elements_; i++)
    {
        mch_.dynamic_elements_[i]->applyForce(force_vector_,
                                              state.segment(0, mch_.nr_of_states_),
                                              sin_cos_list_);
    }

    // Compute forces induced by dissipative and active elements (dampers, motors etc.)
    for(int i = 0; i < mch_.nr_of_motor_elements_; i++)
    {
        mch_.motor_elements_[i]->applyForce(force_vector_,
                                            t,
                                            state.segment(0, mch_.nr_of_states_),
                                            state.segment(mch_.nr_of_states_, mch_.nr_of_states_));
    }
}

void StateEquation::initializeMassMatrix()
{
    // Initialize the mass matrix and mass list
    M_ = VectorXd::Zero(mch_.nr_of_states_);

    for(int i = 0; i < mch_.nr_of_masses_; i++)
    {
        // x mass
        M_(i*3) = 1./mch_.links_[i]->mass();
        // y mass
        M_(i*3+1) = 1./mch_.links_[i]->mass();
        // inertia
        M_(i*3+2) = 1./mch_.links_[i]->inertia();
    }

    if(mch_.dna_.withEndEffector()) {
        M_(3*mch_.nr_of_masses_) = 1./mch_.dna_.parameters(0)(2);
        M_(3*mch_.nr_of_masses_ + 1) = 1./mch_.dna_.parameters(0)(2);
    }
}

void StateEquation::initializeGravity()
{
    for(int i = 0; i < mch_.nr_of_masses_; i++)
    {
        gravity_vector_(3*i+1) = -mch_.links_[i]->mass();
    }

    if(mch_.dna_.withEndEffector()) {
        gravity_vector_(3*mch_.nr_of_masses_ + 1) = -mch_.dna_.parameters(0)(2);
    }
}

void StateEquation::initializeJacobian()
{
    Dk_ = MatrixXd::Zero(mch_.nr_of_constraints_, mch_.nr_of_states_);
    Dk_T_ = MatrixXd::Zero(mch_.nr_of_states_, mch_.nr_of_constraints_);

    for(int i = 0; i < mch_.nr_of_constraint_connections_; i++)
    {
        mch_.constraint_connections_[i]->initializeJacobian(Dk_.block(i*2,0,2,Dk_.cols()));
    }

}
