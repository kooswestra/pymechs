#include "bindings.h"

using namespace Eigen;
using namespace Mech;

// Define the bindings for the mechanism, controller, dna and objective classes
PYBIND11_MODULE(cmechlib, m) {
    m.doc() = "A library for 2D mechanism simulation using a graph based DNA Representation";

    // Bindings for the DNA class
    py::class_<DNA> dna(m, "_DNA");
    // Parametrized constructor
    dna.def(py::init<MatrixXi,
        VectorXi,
        std::vector<RowVectorXd>,
        std::vector<RowVectorXd>>())

        // Copy constructor, required to build pydna from a cdna reference that belongs to a mechanism
        .def(py::init<const DNA&>())

        // Accessible properties of the DNA
        .def_property_readonly("incidenceMatrix", &DNA::incidenceMatrix)
        .def("edgeLabels", py::overload_cast<>(&DNA::edgeLabels, py::const_))
        .def("edgeLabels", py::overload_cast<int>(&DNA::edgeLabels, py::const_))
        .def("parameters", py::overload_cast<>(&DNA::parameters, py::const_))
        .def("parameters", py::overload_cast<int>(&DNA::parameters, py::const_))
        .def("masses", py::overload_cast<>(&DNA::masses, py::const_))
        .def("masses", py::overload_cast<int>(&DNA::masses, py::const_))
        .def_static("incMatrixToConnections", &DNA::incMatrixToConnections)
        .def("linkConnections", &DNA::linkConnections)
        .def("edgeConnections", &DNA::edgeConnections)
        .def("isGrounded", &DNA::isGrounded)
        .def("isDynamic", &DNA::isDynamic)
        .def("isConnected", &DNA::isConnected);

    // Bindings for the mechanism class
    py::class_<Mechanism> mechanism(m,"_Mechanism");
    mechanism.def(py::init<const DNA&>())
            .def(py::init<MatrixXi,
                 VectorXi,
                 std::vector<RowVectorXd>,
                 std::vector<RowVectorXd>>())
            .def("simulate", &Mechanism::simulate)
            .def("assignControllers", &Mechanism::assignControllers)
            .def("evaluate", py::overload_cast<Objective&, double, int>(&Mechanism::evaluate))
            .def("statesTime", &Mechanism::statesTime, py::return_value_policy::reference_internal)
            .def_property_readonly("dna", &Mechanism::dna, py::return_value_policy::reference_internal)
            .def_property_readonly("momentarmList", &Mechanism::momentarmList)
            .def_property_readonly("isSimulated", &Mechanism::isSimulated)
            .def_property_readonly("time", &Mechanism::time)
            .def_property_readonly("fitness", &Mechanism::fitness)
            .def_property_readonly("initialState", &Mechanism::initialState);

    // Bindings for the mutator class
    py::class_<Mutator> mutator(m,"_Mutator");
    mutator.def(py::init<VectorXi, std::vector<RowVectorXd>>())
            .def("eval", &Mutator::eval)
            .def("reduce", &Mutator::reduce)
            .def(py::self * py::self)
            .def_property_readonly("structure", &Mutator::structure)
            .def_property_readonly("parameters", &Mutator::parameters);

    // Controller bindings
    py::class_<Controller, std::shared_ptr<Controller>, PyController  /* the trampoline class*/> controller(m, "_Controller");
    controller.def(py::init<>());

    // Pid controller example implementation in C++
    py::class_<Pid, std::shared_ptr<Pid>, Controller> pid(m, "_Pid");
    pid.def(py::init<double, VectorXd>());

    py::class_<FeedForward, std::shared_ptr<FeedForward>, Controller> feedforward(m, "_FeedForward");
    feedforward.def(py::init<double, VectorXd, VectorXd>())
        .def_property_readonly("controlEffort", &FeedForward::controlEffort);

    // Bindings for the virtual objective class
    py::class_<Objective, std::shared_ptr<Objective>, PyObjective> objective(m, "_Objective");
    objective.def(py::init<>());

    // Pick and place task objective implementation in c++
    py::class_<PickAndPlace, std::shared_ptr<PickAndPlace>, Objective> pickandplace(m, "_PickAndPlace");
    pickandplace.def(py::init<RowVector2d, RowVector2d>())
        .def("evaluate", &PickAndPlace::evaluate);

    // Bindings for the threader
    py::class_<Threader> threader(m,"Threader");
    threader.def(py::init<int>())
            .def("simulate", &Threader::simulate)
            .def("evaluate", &Threader::evaluate);

    // Exports the enums
    py::enum_<DNA::Labels>(dna,"Labels")
            .value("HINGE", DNA::Labels::HINGE)
            .value("END_EFFECTOR", DNA::Labels::END_EFFECTOR)
            .value("SPRING", DNA::Labels::SPRING)
            .value("TORSION_SPRING", DNA::Labels::TORSION_SPRING)
            .value("HINGE_MOTOR", DNA::Labels::HINGE_MOTOR)
            .value("PRISMATIC_MOTOR", DNA::Labels::PRISMATIC_MOTOR);

    py::enum_<Mutator::Types>(mutator,"Types")
            .value("PARAMETER", Mutator::Types::PARAMETER)
            .value("RELABEL", Mutator::Types::RELABEL)
            .value("TRANSFORM", Mutator::Types::TRANSFORM)
            .value("MOVE_EFFECTOR", Mutator::Types::MOVE_EFFECTOR)
            .value("ADD_EDGE", Mutator::Types::ADD_EDGE)
            .value("REMOVE_EDGE", Mutator::Types::REMOVE_EDGE)
            .value("ADD_VERTEX", Mutator::Types::ADD_VERTEX)
            .value("REMOVE_VERTEX", Mutator::Types::REMOVE_VERTEX);
}
