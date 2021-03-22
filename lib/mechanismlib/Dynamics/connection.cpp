#include "connection.h"

using namespace Eigen;
using namespace Mech;

Connection::Connection(Vector2i element_connections, bool is_ground_connection)
    :   is_ground_connection_(is_ground_connection),
        element_connections_(element_connections)
{
    // Extract the state indices of the connected elements
    indx1_ = (element_connections[0]-1)*3;
    indx2_ = (element_connections[1]-1)*3;
}

Vector2i Connection::elementConnections() const
{
    return element_connections_;
}

bool Connection::isGroundConnection() const
{
    return is_ground_connection_;
}
