# Import the mechanism library
import pymechlib as mch
# Import the labels for clarity
L = mch.DNA.Labels

# Define a weird springy mech
# Incidence Matrix
inc_matrix = [[0, 1, 0],
              [0, 1, 1],
              [1, 0, 1],
              [1, 0, 0]]

# Labels, two hinges two springs
edge_labels = [L.END_EFFECTOR, L.TORSION_SPRING, L.TORSION_SPRING]

# Define the link parameters
masses = [[0, 0, 1],
          [0, 1, 1]]

# Define hinge location
t1 = [0, 0, 0, 4]
t2 = [0.2, 1, 0, 1]
# Define end-effector location and mass
e1 = [1, 1.4, 1]

# Parameters is the ordered list of the parameters belonging to the edge labels
parameters = [e1, t1, t2]

# Build a dna object from the representation
dna = mch.DNA(inc_matrix, edge_labels, masses, parameters)
# Build a mechanism object using this dna object
mech = mch.Mechanism(dna)

# Simulate calls the simulation function of the C library
mech.simulate(20, 2000)
# Animate the mechanism
mech.animate()
