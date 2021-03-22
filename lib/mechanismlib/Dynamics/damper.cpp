#include "damper.h"

using namespace Eigen;
using namespace Mech;

template<>
Damper<Grounded>::Damper(Vector2i element_connections)
	: DynamicElement(element_connections, Grounded)
{

}

template<>
void Damper<Grounded>::applyForce(Ref<VectorXd> force_vector, const Ref<const VectorXd>& q, const Ref<const VectorXd>& qd)
{
	double c = 2;
	force_vector(indx1_ + 2) = -c * qd(indx1_ + 2);
}


template<>
Damper<Floating>::Damper(Vector2i element_connections)
	: DynamicElement(element_connections, Floating)
{

}

template<>
void Damper<Floating>::applyForce(Ref<VectorXd> force_vector, const Ref<const VectorXd>& q, const Ref<const VectorXd>& qd)
{
	double c = 2;
}