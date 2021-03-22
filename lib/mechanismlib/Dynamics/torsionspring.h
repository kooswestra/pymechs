#ifndef TORSIONSPRING_H
#define TORSIONSPRING_H

#include "../common.h"
#include <Eigen/Dense>
#include "dynamicelement.h"

namespace Mech {

template <bool type>
class TorsionSpring : public DynamicElement
{
public:
    TorsionSpring(Eigen::Matrix<double,1,4> T, Eigen::Vector2i element_connections);

    void applyForce(VectorRef force_vector, ConstVectorRef& state, ConstMatrixRef& sin_cos_list) override;
                    
    int edgeLabel() const override {return 3;}

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    double t0_;
    double k_;
};

}

#endif // TORSIONSPRING_H
