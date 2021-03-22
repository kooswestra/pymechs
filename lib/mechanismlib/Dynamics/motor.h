#ifndef MOTOR_H
#define MOTOR_H

#include <Eigen/Dense>
#include "motorelement.h"

namespace Mech {

template <bool type>
class Motor : public MotorElement
{
public:
    // Attach either a torque for sinusoidal application (testing purposes) or a controller?
    Motor(Eigen::Vector2i element_connections);

    void applyForce(VectorRef force_vector, double t, ConstVectorRef& q, ConstVectorRef& qd) override;
    int edgeLabel() const override {return 4;}

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:

};

}

#endif // MOTOR_H
