#ifndef DAMPER_H
#define DAMPER_H
#include "../common.h"
#include <Eigen/Dense>
#include "dynamicelement.h"

namespace Mech {

template <bool type>
class Damper : public DynamicElement
{
public:
	Damper(Eigen::Vector2i element_connections);
	void applyForce(VectorRef force_vector, ConstVectorRef& state, ConstMatrixRef& sin_cos_list);

private:
	int edgeLabel() const override { return 2; }

};

}

#endif //DAMPER_H