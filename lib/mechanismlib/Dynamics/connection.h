#ifndef CONNECTION_H
#define CONNECTION_H

#include <Eigen/Dense>

/** @brief Connection is the abstract class that implements connections in the mechanism
 *
 * Each connection has the information of which links it connects, as well as a corresponding
 * edge label that defines which type of connection it is. All connection elements inherit
 * from this class
  */

namespace Mech {

class Connection
{
public:
    Connection(Eigen::Vector2i element_connections, bool is_ground_connection);

    /** Each connection contains a vector of the element numbers it connects which can be read*/
    Eigen::Vector2i elementConnections() const;
    /** Each connection knows if it is connected to the ground */
    bool isGroundConnection() const;
    /** Each connection has a corresponding edge label that can be read*/
    virtual int edgeLabel() const = 0;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
    bool is_ground_connection_;
    Eigen::Vector2i element_connections_;

    int indx1_ = 0;
    int indx2_ = 0;
};

}

#endif // CONNECTION_H
