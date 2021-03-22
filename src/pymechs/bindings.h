#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "Evolution/dna.h"
#include "Evolution/mechanism.h"
#include "Evolution/objective.h"
#include "Evolution/mutator.h"
#include "Evolution/pickandplace.h"
#include "Control/pid.h"
#include "Control/feedforward.h"
#include "threader.h"

namespace py = pybind11;

// In order to inherit from the controller, and define the control laws from python this trampoline class has
// to be constructed. It handles the virtual function overloading

namespace Mech {

class PyController : public Controller
{
public:
    using Controller::Controller;

    // Trampoline method for the virtual Controller class
    double controlOutput(double t, const Eigen::Ref<const Eigen::VectorXd>& q, const Eigen::Ref<const Eigen::VectorXd>& qd) override {
        PYBIND11_OVERLOAD_PURE(
                    double,
                    Controller,
                    controlOutput,
                    t,
                    q,
                    qd);
    }
};

class PyObjective : public Objective
{
public:
    using Objective::Objective;

    // Trampoline method for the virtual Objective class
    Eigen::VectorXd evaluate(Mechanism* mechanism) override {
        PYBIND11_OVERLOAD_PURE(
                    Eigen::VectorXd,
                    Objective,
                    evaluate,
                    mechanism);
    }
};

}
