#include "motorelement.h"

using namespace Eigen;
using namespace Mech;

void MotorElement::assignController(std::shared_ptr<Controller> controller)
{
    controller_ = controller;
}
