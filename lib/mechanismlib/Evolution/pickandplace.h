#ifndef PICKANDPLACE_H
#define PICKANDPLACE_H

#include <Eigen/Dense>
#include "objective.h"

namespace Mech {

class PickAndPlace : public Objective
{
public:
    PickAndPlace(Eigen::RowVector2d pick_location, Eigen::RowVector2d place_location);
    Eigen::VectorXd evaluate(Mechanism* mechanism) override;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Eigen::RowVector2d pick_location_;
    Eigen::RowVector2d place_location_;
    double l0_;
};

}


#endif // PICKANDPLACE_H
