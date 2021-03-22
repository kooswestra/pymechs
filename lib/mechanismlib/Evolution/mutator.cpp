#include "mutator.h"

using namespace Eigen;
using namespace Mech;

Mutator::Mutator(VectorXi structure, std::vector<RowVectorXd> parameters)
    :   structure_(structure),
        parameters_(parameters)
{
    // This does not actually activate the operator yet, this is to prevent computation
    // until the full operator is determined from the combination of multiple operators

    if(structure.size() != parameters.size())
        throw std::invalid_argument("Invalid mutator defined, each label must have corresponding parameters");
}

// Overloaded method, combines operators into a new operator
/** @todo change this at some point to reduce copy and allocation overhead for large chains of operators */
Mutator Mutator::operator *(const Mutator &other) const
{
    VectorXi new_structure(structure_.size() + other.structure().size());
    new_structure << other.structure(), structure_;

    std::vector<RowVectorXd> new_parameters(other.parameters());

    for(int i = 0; i < parameters_.size(); i++) {
        new_parameters.push_back(parameters_[i]);
    }

    return Mutator(new_structure, new_parameters);
}

void Mutator::activate(DNA const &dna)
{
    int vertex_adds = 0;
    int vertex_del = 0;
    int edge_adds = 0;
    int edge_del = 0;
    int end_effector_movs = 0;

    // Interpret the complete representation and form the operators, do some basic sanity checks
    for(int i = 0; i < structure_.size(); i++)
    {
        switch(structure_[i]) {
        case ADD_VERTEX:
            vertex_additions_.push_back(dna.incidenceMatrix().rows() + vertex_adds - vertex_del);
            vertex_adds++;
            break;
        case REMOVE_VERTEX:
            if(parameters_[i](0) > dna.incidenceMatrix().rows() + vertex_adds - vertex_del - 1) {
                throw std::invalid_argument("Operator tried to remove a link that doesn't exist");
            }
            if(parameters_[i](0) == 0) {
                throw std::invalid_argument("Operator tried to remove the ground");
            }
            vertex_removals_.push_back(parameters_[i](0));
            vertex_del++;
            break;
        case ADD_EDGE:
            edge_additions_.push_back(dna.incidenceMatrix().cols() + edge_adds - edge_del);
            if(parameters_[i](0) == parameters_[i](1)) {
                throw std::invalid_argument("Cannot connect an element to itself");
            }
            edge_adds++;
            break;
        case REMOVE_EDGE:
            if(parameters_[i](0) > dna.incidenceMatrix().cols() + edge_adds - edge_del - 1) {
                throw std::invalid_argument("Operator tried to remove a connection that doesn't exist");
            }
            edge_removals_.push_back(parameters_[i](0));
            edge_del++;
            break;
        case RELABEL:
            if(parameters_[i](0) > dna.incidenceMatrix().cols() + edge_adds - edge_del) {
                throw std::invalid_argument("Operator tried to modify the label of a connection that doesn't exist");
            }
            label_mods_.push_back(parameters_[i](0));
            break;
        case PARAMETER:
            parameter_mods_.push_back(parameters_[i](0));
            break;
        case MOVE_EFFECTOR:
            end_effector_movs++;
            break;
        case TRANSFORM:
            break;
        default:
            throw std::invalid_argument("invalid operation specified, can't activate operator");
            break;
        }
    }

    // Check if only parameter mutations are defined, if this is the case, it is a reduced complexity parametric mutator
    if(vertex_adds + vertex_del + edge_adds + edge_del + end_effector_movs == 0) {
        for(int i = 0; i < dna.parameters().size(); i++) {
            temp_parameters_.push_back(dna.parameters(i));
        }
        for (int i = 0; i < dna.masses().size(); i++) {
            temp_masses_.push_back(dna.masses(i));
        }
        temp_edge_labels_ = dna.edgeLabels();
        parameter_only_ = true;
        active_ = true;
        return;
    }

    // Possibly apply logical reductions if vertices or edges are both added and removed
    if((vertex_adds > 0 && vertex_del > 0) || (edge_adds > 0 && edge_del > 0)) {
        reduce();
    }

    active_ = true;
}

void Mutator::order()
{

}

void Mutator::reduce()
{
    /** @todo Check if a link or connection is added and then removed and vice versa, and delete these duplicate operations */
}

void Mutator::initializeRepresentation(const DNA &dna)
{
    // Initialize (possibly resized) temporary elements for the new dna structure
    int cols = dna.incidenceMatrix().cols() + edge_additions_.size() - edge_removals_.size();
    int rows = dna.incidenceMatrix().rows() + vertex_additions_.size() - vertex_removals_.size();

    temp_inc_matrix_ = MatrixXi::Zero(rows, cols);
    temp_edge_labels_ = VectorXi(cols);

    // If the incidence matrix downsizes, check which elements are removed
    if(vertex_additions_.size() < vertex_removals_.size() || edge_additions_.size() < edge_removals_.size()) {
        // Calculate preserved column and preserved row numbers
        std::vector<int> preserved_rows;
        std::vector<int> preserved_columns;

        // Inner iterator
        int k = 0;
        // Ground is always preserved
        preserved_rows.push_back(0);
        for(int i = 1; i < dna.incidenceMatrix().rows() - (dna.withEndEffector() ? 1 : 0); i++) {
            // If the row is not in the vertex removal list, store it and copy its parameters
            if(std::find(vertex_removals_.begin(), vertex_removals_.end(), i) == vertex_removals_.end()) {
                preserved_rows.push_back(i);
                temp_masses_.push_back(dna.masses(i-1));
                k++;
            }
        }
        // The end-effector is also always preserved, if it exists
        if(dna.withEndEffector())
            preserved_rows.push_back(dna.incidenceMatrix().rows() - 1);

        // Inner iterator
        int l = 0;
        for(int i = 0; i < dna.incidenceMatrix().cols(); i++) {
            // If the column is not in the edge removal list, store it and copy it's edge label
            if(std::find(edge_removals_.begin(), edge_removals_.end(), i) == edge_removals_.end()) {
                preserved_columns.push_back(i);
                temp_edge_labels_(l) = dna.edgeLabels(i);

                // Copy the parameter belonging to the edge label as well
                temp_parameters_.push_back(dna.parameters(i));
                l++;
            }
        }

        // Copy over all the relevant elements of the incidence matrix
        for(int i = 0; i < temp_inc_matrix_.rows(); i++) {
            for(int j = 0; j < temp_inc_matrix_.cols(); j++) {
                temp_inc_matrix_(i,j) = dna.incidenceMatrix()(preserved_rows[i], preserved_columns[j]);
            }
        }
    }

    else {
        temp_inc_matrix_.block(0 , 0, dna.incidenceMatrix().rows(), dna.incidenceMatrix().cols()) = dna.incidenceMatrix();
        temp_edge_labels_.segment(0, dna.edgeLabels().size()) = dna.edgeLabels();
        for (int i = 0; i < dna.parameters().size(); i++) {
            temp_parameters_.push_back(dna.parameters(i));
        }
        for (int i = 0; i < dna.masses().size(); i++) {
            temp_masses_.push_back(dna.masses(i));
        }
        // Swap the end-effector row to the bottom
        temp_inc_matrix_.bottomRows(1).swap(temp_inc_matrix_.row(dna.incidenceMatrix().rows()-1));
    }
}

// Evaluate the operator on the DNA of a mechanism to generate a new dna
DNA Mutator::eval(const DNA &dna)
{
    // Activate the operator if it is dormant
    if(!active_) {
        activate(dna);
    }

    // Make an internal (properly resized) copy of the dna to work with if it changes size
    if(!parameter_only_)
        initializeRepresentation(dna);

    int vertex_adds = dna.incidenceMatrix().rows();
    int edge_adds = dna.incidenceMatrix().cols();

    for(int i = 0; i < structure_.size(); i++)
    {
        switch(structure_[i]) {
        case ADD_VERTEX:
            temp_masses_.push_back(parameters_[i]);
            vertex_adds++;
            break;
        case ADD_EDGE:
            // Define the new connection in the incidence matrix
            temp_inc_matrix_(parameters_[i](0), edge_adds) = 1;

            // Special case link number -1 refers to the last added element
            if(parameters_[i](1) < 0)
                temp_inc_matrix_(vertex_adds - 2, edge_adds) = 1;
            else
                temp_inc_matrix_(parameters_[i](1), edge_adds) = 1;

            // Set the edge label for this new edge
            temp_edge_labels_(edge_adds) = parameters_[i](2);
            // Add the parameter
            temp_parameters_.push_back(parameters_[i].segment(3, parameters_[i].size()-3));
            edge_adds++;
            break;
        case RELABEL:
            // Redefine the label at the right spot
            temp_edge_labels_(parameters_[i](0)) = parameters_[i](1);
            temp_parameters_[parameters_[i](0)] = parameters_[i].segment(2, parameters_[i].size()-2);
            break;
        case PARAMETER:
            // Modify the parameters of the edge or mass
            if(parameters_[i](0) < temp_masses_.size() - 1)
                temp_masses_[parameters_[i](0)] += parameters_[i].segment(1, parameters_[i].size() - 1);
            else
                temp_parameters_[parameters_[i](0) - temp_masses_.size()] += parameters_[i].segment(1, parameters_[i].size() - 1);
            break;
        case MOVE_EFFECTOR:
            temp_inc_matrix_.col(0).setZero();
            temp_inc_matrix_(parameters_[i](0), 0) = 1;
            temp_inc_matrix_(temp_inc_matrix_.rows()-1, 0) = 1;
            break;
        case TRANSFORM:
            // Rotate, scale and offset the position parameters.
            Matrix2d R;
            R << cos(parameters_[i](0)), -sin(parameters_[i](0)), sin(parameters_[i](0)), cos(parameters_[i](0));

            for(int j = 0; j < temp_parameters_.size(); j++) {
                temp_parameters_[j].segment<2>(0) = temp_parameters_[j].segment<2>(0)*R.transpose()*parameters_[i](1) + parameters_[i].segment<2>(2);
                if(temp_edge_labels_[i] == DNA::SPRING)
                    temp_parameters_[j].segment<2>(2) = temp_parameters_[j].segment<2>(2)*R.transpose()*parameters_[i](1) + parameters_[i].segment<2>(2);
            }

            for(int j = 0; j < temp_masses_.size(); j++) {
                temp_masses_[j].segment<2>(0) = temp_masses_[j].segment<2>(0)*parameters_[i](1)*R.transpose() + parameters_[i].segment<2>(2);
            }
            break;
        }
    }

    // Parametric only results in a reduced operator form, which is more efficient
    if(parameter_only_) {
        DNA modified_dna(dna.incidenceMatrix(), temp_edge_labels_, temp_masses_, temp_parameters_);
        return modified_dna;
    }

    // -> Preferably move rather than copy the data
    DNA modified_dna(temp_inc_matrix_, temp_edge_labels_, temp_masses_, temp_parameters_);
    return modified_dna;
}

VectorXi Mutator::structure() const
{
    return structure_;
}

std::vector<RowVectorXd> Mutator::parameters() const
{
    return parameters_;
}
