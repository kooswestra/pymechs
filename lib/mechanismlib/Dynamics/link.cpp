#include "link.h"
#include <vector>

using namespace Eigen;
using namespace Mech;

Link::Link(RowVector3d mass, Matrix<double, Dynamic, 2> connection_locations, VectorXi element_connections, int link_number)
    :   link_number_(link_number),
        mass_(abs(mass(2))),
        element_connections_(element_connections)
{
    centre_of_mass_ = mass.segment<2>(0);

    switch(element_connections.size()) {
    case 1 :
        inertia_ = mass_*(connection_locations.row(0) - centre_of_mass_).squaredNorm() / 3;
        break;
    case 2 :
        centre_of_mass_ = (connection_locations.row(1) + connection_locations.row(0)) / 2;
        inertia_ = mass_*(connection_locations.row(0) - centre_of_mass_).squaredNorm() / 3;
        break;
    default :
        buildPolygon(connection_locations);
        calculatePolygonInertia();
        break;
    }
}

void Link::buildPolygon(const Ref<const Matrix<double, Dynamic, 2>> unsorted_polygon)
{
    std::vector<double> angles;
    std::vector<int> indx;
    RowVector2d mean = RowVector2d::Zero();
    polygon_ = Matrix<double, Dynamic, 2>::Zero(unsorted_polygon.rows(), 2);

    // First determine mean
    if(element_connections_.size() > 1) {
        for(int i = 0; i < unsorted_polygon.rows(); i++) {
            mean = mean + unsorted_polygon.row(i)/unsorted_polygon.rows();
         }
    }

    // Then determine  all the incident angles
    for(int i = 0; i < polygon_.rows(); i++) {
        angles.push_back(atan2(unsorted_polygon(i,1) - mean(1), unsorted_polygon(i,0) - mean(0)));
        indx.push_back(i);
    }

    // Finally sort the polygon by incident angle from the mean
    std::sort(indx.begin(), indx.end(), [&](size_t i, size_t j){return angles[i] < angles[j];});

    for(int i = 0; i < polygon_.rows(); i++) {
        polygon_.row(i) = unsorted_polygon.row(indx[i]);
    }

}

void Link::calculatePolygonInertia()
{
    inertia_ = 0;
    double area = 0;

    // First determine the area of the polygon, as the density is given by m/A
    for(int i = 0; i < polygon_.rows(); i++) {
        // j = (i+1) with wraparound such that j = 0 for i = polygon_.rows()
        int j = (i+1) % polygon_.rows();

        area += polygon_(i,0)*polygon_(j,1) - polygon_(i,1)*polygon_(j,0);
    }

    area = area / 2;

    centre_of_mass_ = RowVector2d::Zero();

    // Determine center of mass, which is the centroid of the polygon
    for(int i = 0; i < polygon_.rows(); i++) {
        // j = (i+1) with wraparound such that j = 0 when i = polygon_.rows()
        int j = (i+1) % polygon_.rows();

        centre_of_mass_(0) += (polygon_(j,0) + polygon_(i,0))*(polygon_(i,0)*polygon_(j,1) - polygon_(i,1)*polygon_(j,0));
        centre_of_mass_(1) += (polygon_(j,1) + polygon_(i,1))*(polygon_(i,0)*polygon_(j,1) - polygon_(i,1)*polygon_(j,0));
    }

    centre_of_mass_ = centre_of_mass_/(6*area);

    for(int i = 0; i < polygon_.rows(); i++) {
        polygon_.row(i) = polygon_.row(i) - centre_of_mass_;
    }


    // Polygonal moment of inertia formula, for derivation see thesis appendix
    for(int i = 0; i < polygon_.rows(); i++) {
        // j = (i+1) with wraparound such that j = 0 when i = polygon_.rows()
        int j = (i+1) % polygon_.rows();

        inertia_ += (polygon_(j,1) - polygon_(i,1))*(polygon_(j,0) + polygon_(i,0))*(pow(polygon_(j,0),2) + pow(polygon_(i,0),2))
                    - (polygon_(j,0) - polygon_(i,0))*(polygon_(j,1) + polygon_(i,1))*(pow(polygon_(j,1),2) + pow(polygon_(i,1),2));
    }

    inertia_ = mass_*inertia_ / (12*area);

}

RowVector2d Link::centreOfMass() const
{
    return centre_of_mass_;
}

VectorXi Link::elementConnections() const
{
    return element_connections_;
}

double Link::inertia() const
{
    return inertia_;
}

double Link::mass() const
{
    return mass_;
}
