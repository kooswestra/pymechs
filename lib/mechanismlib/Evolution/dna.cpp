#include "dna.h"
#include <list>

using namespace Eigen;
using namespace Mech;

// Raw constructor of the dna object, need to clean up parameter handling
DNA::DNA(MatrixXi incidence_matrix,
         VectorXi edge_labels,
         std::vector<RowVectorXd> masses,
         std::vector<RowVectorXd> parameters)
    : incidence_matrix_(incidence_matrix),
      edge_labels_(edge_labels)
{
    // Initialize the parameters
    parameter_vector_ = parameters;
    mass_vector_ = masses;
    complexity_ = 0;

    // Do some sanity checking on the construction of the DNA as this function can be directly interfaced from python
    // and as such has no checks on if the right incidence matrix with matching edge labels is passed

    if(edge_labels.size() != incidence_matrix.cols())
        throw std::invalid_argument("Every edge must have a corresponding label");

    // Check if the parameters match the labels and are sane

    for(int i = 0; i < edge_labels.size(); i++) {
        switch (edge_labels(i)) {
        case HINGE :
            if(parameters[i].size() != 2)
                throw std::invalid_argument("Supplied parameters do not match HINGE edge label");
            complexity_ += 2;
            break;
        case END_EFFECTOR :
            if(parameters[i].size() != 3)
                throw std::invalid_argument("Supplied parameters do not match END_EFFECTOR edge label");
            if(i != 0)
                throw std::invalid_argument("End-effector connection must always be the first edge");
            if(parameters[i](2) < 0)
                throw std::invalid_argument("End-effector mass cannot be negative");
            with_end_effector_ = true;
            break;
        case SPRING :
            if(parameters[i].size() != 6)
                throw std::invalid_argument("Supplied parameters do not match SPRING edge label");
            complexity_ += 6;
            break;
        case TORSION_SPRING :
            if(parameters[i].size() != 4)
                throw std::invalid_argument("Supplied parameters do not match TORSION_SPRING edge label");
            complexity_ += 4;
            break;
        case HINGE_MOTOR :
            if(parameters[i].size() != 3)
                throw std::invalid_argument("Supplied parameters do not match HINGE_MOTOR edge label");
            complexity_ += 3;
            break;
        default :
            throw std::invalid_argument("An invalid or not yet implemented edge label has been passed: " + std::to_string(edge_labels(i)));
        }
    }

    for (int i = 0; i < incidence_matrix.rows() - 1 - with_end_effector_; i++) {
        if (masses[i].size() != 3)
            throw std::invalid_argument("Supplied parameters do not match MASS node label");
        complexity_ += 5;
    }
}

MatrixXi DNA::incidenceMatrix() const
{
    return incidence_matrix_;
}

int DNA::edgeLabels(int index) const
{
    return edge_labels_(index);
}

std::vector<RowVectorXd> DNA::parameters() const
{
    return parameter_vector_;
}

RowVectorXd DNA::parameters(int index) const
{
    return parameter_vector_[index];
}

VectorXi DNA::edgeLabels() const
{
    return edge_labels_;
}

std::vector<RowVectorXd> DNA::masses() const
{
    return mass_vector_;
}

RowVectorXd DNA::masses(int index) const
{
    return mass_vector_[index];
}

bool DNA::withEndEffector() const
{
    return with_end_effector_;
}

int DNA::complexity() const
{
    return complexity_;
}

VectorXi DNA::incMatrixToConnections(VectorXi incidence_vector)
{
    // Calculates element numbers that are connected
    VectorXi element_connections(incidence_vector.sum());
    int j = 0;

    for(int i = 0; i < incidence_vector.size(); i++)
    {
        if(incidence_vector[i] == 1)
        {
            element_connections[j] = i;
            j++;
        }
    }

    return element_connections;
}

Vector2i DNA::edgeConnections(int index) const
{
    return incMatrixToConnections(incidence_matrix_.col(index));
}

VectorXi DNA::linkConnections(int index) const
{
    return incMatrixToConnections(incidence_matrix_.row(index));
}

bool DNA::isGrounded() const
{
    // Check if the graph with spring connections deleted represents a fully connected graph
    // contrary to isConnected() this DOES include the ground.

    // The connectivity matrix is the incidence matrix with all spring connections deleted
    // Determine which columns contain springs, delete those
    std::vector<int> indices;
    for(int i = 0; i < incidence_matrix_.cols(); i++) {
        if(edge_labels_(i) != DNA::SPRING) {
            indices.push_back(i);
        }
    }

    // Build the connectivity matrix using these indices
    MatrixXi connectivity_matrix = MatrixXi::Zero(incidence_matrix_.rows(), indices.size());
    for(int i = 0; i < indices.size(); i++) {
        connectivity_matrix.col(i) = incidence_matrix_.col(indices[i]);
    }

    return graphTraversable(connectivity_matrix);
}

bool DNA::isDynamic() const
{
    int nr_of_states = mass_vector_.size()*3;
    int nr_of_constraints = 0;

    // Check if the incidence matrix describes a mechanism that has a degree of freedom > 0 by counting constraints
    // This does not check explicitly for singularities

    for(int i = 0; i < edge_labels_.size(); i++) {
        switch(edge_labels_(i)) {
            case DNA::HINGE : {
                nr_of_constraints += 2;
                break;
            }
            case DNA::HINGE_MOTOR : {
                nr_of_constraints += 2;
                break;
            }
            case DNA::END_EFFECTOR : {
                nr_of_constraints += 0;
                break;
            }
            case DNA::PRISMATIC_MOTOR : {
                nr_of_constraints += 2;
                break;
            }
            case DNA::TORSION_SPRING : {
                nr_of_constraints += 2;
                break;
            }
        }
    }
    return(nr_of_states > nr_of_constraints);
}

bool DNA::isConnected() const
{
    // Check if the incidence matrix describes a fully connected mechanism,
    // i.e. there exists a path from every link to every other link (graph traversal)

    // Initialize traversal list as false
    // Then start exploring, the ground does NOT count. As it is not technically part of the mechanism structure
    // The functions isGrounded() checks for ground connectivity specifically. Start from the first element
    MatrixXi connectivity_matrix = incidence_matrix_.block(1,0,incidence_matrix_.rows()-1, incidence_matrix_.cols());

    return graphTraversable(connectivity_matrix);
}

bool DNA::graphTraversable(const Ref<const MatrixXi> connectivity_matrix) const
{
    Array<bool,Dynamic,1> traversed = Array<bool,Dynamic,1>::Zero(connectivity_matrix.rows());

    std::list<int> q;
    q.push_back(0);
    traversed[0] = true;

    while(!q.empty()) {

        int v = q.front();
        // Check connected edges
        VectorXi edges = incMatrixToConnections(connectivity_matrix.row(v));
        q.pop_front();

        // Visit the nodes connected by those edges, but only the ones that haven't been visited yet
        for(int i = 0; i < edges.size(); i++) {
            VectorXi nodes = incMatrixToConnections(connectivity_matrix.col(edges(i)));
            for(int j = 0; j < nodes.size(); j++) {
                if(!traversed[nodes(j)]) {
                    traversed[nodes(j)] = true;
                    q.push_back(nodes(j));
                }
            }
        }
    }


   // If any elements are not traversed, return false. Otherwise the graph is connected
    return traversed.all();
}
