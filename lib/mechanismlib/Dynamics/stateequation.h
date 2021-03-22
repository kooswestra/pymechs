#ifndef STATEEQUATION_H
#define STATEEQUATION_H

#include "../common.h"
#include <Eigen/Dense>
#include "../Evolution/dna.h"
#include "../Evolution/mechanism.h"

/** @brief StateEquation calculates the equations of motion for a corresponding mechanism structure
 *
 * It implements the actual calculation of the equations of motion by iterating over the elements
 * present in the mechanism structure and updating the corresponding elements in the state equation
 * given by:
 * \f[
        \begin{bmatrix}\mathbf{M} & J\left(\mathbf{D}\right)^{T}\\
        J\left(\mathbf{D}\right) & \mathbf{0}
        \end{bmatrix}\begin{bmatrix}\ddot{\mathbf{x}}\\
        \boldsymbol{\lambda}
        \end{bmatrix}=\begin{bmatrix}\sum\mathbf{F}\\
        \frac{\partial\left(J\left(\boldsymbol{D}\right)\dot{\mathbf{x}}\right)}{\partial\mathbf{x}}\dot{\mathbf{x}}
        \end{bmatrix}
  \f]
 * where \f$J\left(D\right)=\frac{\partial\mathbf{D}}{\partial\mathbf{x}}\f$ is the Jacobian
 * of the constraints w.r.t. the mechanism state. It does this by calculating the inter-jacobian
 * corresponding with every ConstraintConnection and then updating the main jacobian matrix with these sub blocks
 * at the positions corresponding with the masses this ConstraintConnection connects. It also calculates
 * the convective accelerations corresponding with this constraint. Additionally the forces
 * implemented by DynamicElement and MotorElement elements are calculated and the force vector updated.
 */

namespace Mech {

class StateEquation
{
public:
    StateEquation(Mechanism& mechanism);
    ~StateEquation();
    int dim() const;

    // Delete copy operators
    StateEquation & operator=(const StateEquation&) = delete;
    StateEquation(const StateEquation&) = delete;

    void evaluate(double t, ConstRowVectorRef& state, RowVectorRef stated);

    void updateSinCosList(ConstRowVectorRef& state);

    void updateJacobian(ConstRowVectorRef& state);
    ConstMatrixRef jacobian() const;
    ConstMatrixRef jacobianTranspose() const;

    void updateConvectiveAccelerations(ConstRowVectorRef& state);
    ConstVectorRef convectiveAccelerations() const;

    void updateConstraints(ConstRowVectorRef& state);
    ConstVectorRef constraints() const;

    void updateForces(double t, ConstRowVectorRef& state);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:

    Mechanism& mch_;

    void initializeMassMatrix();
    void initializeGravity();
    void initializeJacobian();

    // A function that determines element connections from the incidence matrix
    Eigen::Vector2i incMatrixToConnections(Eigen::VectorXi incidence_vector);
    int dim_;

    Eigen::VectorXd M_;
    Eigen::VectorXd D_;
    Eigen::Matrix<double, Eigen::Dynamic, 2> sin_cos_list_;
    Eigen::MatrixXd Dk_;
    Eigen::MatrixXd Dk_T_;
    Eigen::VectorXd Dkk_;
    Eigen::VectorXd force_vector_;

    // State, initial condition, force and gravity vectors
    Eigen::VectorXd equation_vector_;
    Eigen::VectorXd gravity_vector_;

};

}

#endif // STATEEQUATION_H
