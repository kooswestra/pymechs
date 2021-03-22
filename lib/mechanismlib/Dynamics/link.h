#ifndef LINK_H
#define LINK_H

#include <Eigen/Dense>

namespace Mech {

class Link
{
public:
    Link(Eigen::RowVector3d mass, Eigen::Matrix<double, Eigen::Dynamic, 2> connection_locations,
            Eigen::VectorXi element_connections, int link_number);
            
    Eigen::RowVector2d centreOfMass() const;
    Eigen::VectorXi elementConnections() const;

    double inertia() const;
    double mass() const;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    int link_number_;
    double mass_;
    double inertia_;

    Eigen::RowVector2d centre_of_mass_;
    Eigen::VectorXi element_connections_;
    Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> polygon_;

    void buildPolygon(const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 2>> connection_locations);
    void calculatePolygonInertia();
};

}

#endif // LINK_H
