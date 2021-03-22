#ifndef MUTATOR_H
#define MUTATOR_H

#include <Eigen/Dense>
#include "dna.h"
#include "mechanism.h"

/** @brief Mutator is the operator class, it represents the variational operators that can be applied to mechanisms
 *
 *  It provides the ability to modify the mechanism structure in preset ways
 *  @todo Needs to be refactored to handle operations more cleanly and add checks, currently quite error prone if
 *  a user tries to apply an operator that results in a bad mechanism
 */

namespace Mech {

class Mutator
{
public:
    Mutator(Eigen::VectorXi structure, std::vector<Eigen::RowVectorXd> parameters);

    // Operator overload, operators can be chained to create a new, larger operator
    Mutator operator *(const Mutator &other) const;
    DNA eval(const DNA &dna);

    enum Types {PARAMETER = 0 ,RELABEL = 1, TRANSFORM = 2, ADD_EDGE = 3, REMOVE_EDGE = 4, ADD_VERTEX = 5, REMOVE_VERTEX = 6, MOVE_EFFECTOR = 7};

    Eigen::VectorXi structure() const;
    std::vector<Eigen::RowVectorXd> parameters() const;

    void reduce();
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    /**Representation of the operator is in this subsection, might move to separate class */
    Eigen::VectorXi structure_;
    std::vector<Eigen::RowVectorXd> parameters_;

    // Keep track of the modifications internally -> part of the representation
    std::vector<int> vertex_additions_;
    std::vector<int> edge_additions_;

    std::vector<int> vertex_removals_;
    std::vector<int> edge_removals_;

    std::vector<int> parameter_mods_;
    std::vector<int> label_mods_;

    // Representation reduction and ordering
    void order();

    /**-----------------------------------------------------------------------------------*/

    bool active_ = false;
    bool parameter_only_ = false;

    // Operator internal DNA structure temporaries
    Eigen::MatrixXi temp_inc_matrix_;
    Eigen::VectorXi temp_edge_labels_;
    std::vector<Eigen::RowVectorXd> temp_parameters_;
    std::vector<Eigen::RowVectorXd> temp_masses_;

    void activate(const DNA &dna);
    void initializeRepresentation(const DNA &dna);
};
}

#endif // MUTATOR_H
