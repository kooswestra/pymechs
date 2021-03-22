#ifndef DNA_H
#define DNA_H
#include <Eigen/Dense>
#include <vector>

/** @brief DNA is a container that contains the full represention of the mechanism, with some added helper functions
 *
 * It holds the graph structure of the mechanism, the edge labels, as well as the parameters corresponding with
 * the Connection elements the edge labels defined in the structure. Additionally, it has checks applied
 * on the graph structure to give certain guarantees on the behaviour of the mechanism the DNA represents. This
 * class is further extended by python
 */

namespace Mech {
class DNA {
public:
  // Direct, raw constructor
  DNA(Eigen::MatrixXi incidenceMatrix,
      Eigen::VectorXi edge_labels,
      std::vector<Eigen::RowVectorXd> masses,
      std::vector<Eigen::RowVectorXd> parameters);

  Eigen::MatrixXi incidenceMatrix() const;
  // Functions for edge labels, all or indexed
  Eigen::VectorXi edgeLabels() const;
  int edgeLabels(int index) const;

  // Function for parameters, all or indexed
  std::vector<Eigen::RowVectorXd> parameters() const;
  Eigen::RowVectorXd parameters(int index) const;

  std::vector<Eigen::RowVectorXd> masses() const;
  Eigen::RowVectorXd masses(int index) const;

  // Edge label enum for more readable code
  enum Labels { HINGE = 0, END_EFFECTOR = 1, SPRING = 2, TORSION_SPRING = 3, HINGE_MOTOR = 4, PRISMATIC_MOTOR = 5 };

  // Calculating connections from the incidence matrix is a static helper function
  static Eigen::VectorXi incMatrixToConnections(Eigen::VectorXi incidence_vector);
  Eigen::Vector2i edgeConnections(int index) const;
  Eigen::VectorXi linkConnections(int index) const;

  // Functions that check if this dna describes a valid mechanism structure
  /** A check if every element of the mechanism is connected to the ground using breadth first graph traversal,
   * specifically ignoring spring connections as these can cause "dangling"*/
  bool isGrounded() const;
  /** A check if the mechanism has a degree of freedom >0, checking if it can actually move */
  bool isDynamic() const;
  /** A check if the mechanism is fully connected using breadth first graph traversal, ignoring ground connections. */
  bool isConnected() const;

  bool withEndEffector() const;

  int complexity() const;

private:
  Eigen::MatrixXi incidence_matrix_;
  Eigen::VectorXi edge_labels_;
  bool with_end_effector_ = false;

  // Stores the parameters
  std::vector<Eigen::RowVectorXd> parameter_vector_;
  std::vector<Eigen::RowVectorXd> mass_vector_;
  int complexity_;

  // Helper function which implements graph traversal for a given incidence matrix
  bool graphTraversable(const Eigen::Ref<const Eigen::MatrixXi> connectivity_matrix) const;
};
}

#endif // DNA_H
