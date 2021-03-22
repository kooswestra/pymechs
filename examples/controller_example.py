# Mechanisms
import pymechlib as mch

# Edge labels
E = mch.DNA.Labels

# This defines a weird mech with two motors
inc_matrix = [[1, 0, 1],
              [1, 1, 1],
              [0, 1, 0]]

edge_labels = [E.HINGE, E.HINGE_MOTOR, E.SPRING]

masses = [[0.5, 0, 5],
          [1.5, 0, 2]]

h1 = [0, 0]
m1 = [1, 0, 0]
s1 = [-1, 1, 0, 0.7, 1.2, 40]

# Parameters is a list of vectors of the parameters corresponding to the labels
parameters = [h1, m1, s1]

# Build a dna object from the structure, checks parameter validity
dna = mch.DNA(inc_matrix, edge_labels, masses, parameters)
# Build a mechanism using this dna
mech = mch.Mechanism(dna)


# An example PD controller model defined using the library, make a class that inherits mch.Controller
class PD(mch.Controller):

    def __init__(self, ref1, ref2, kd, kp):

        # Make sure to explicitly initialize the base
        mch.Controller.__init__(self)
        # Define any additional parameters
        self.Kd = kd
        self.Kp = kp
        self.ref1 = ref1
        self.ref2 = ref2

    # define the controlOutput(self, t, q, qd) function to define the controller behavior
    def controlOutput(self, t, q, qd):

        if t <= 10:
            return self.Kd * (self.ref1 - q[5]) - self.Kp * qd[5]
        else:
            return self.Kd * (self.ref2 - q[5]) - self.Kp * qd[5]


# Initialize an instance with some example parameters
ctrl = PD(-3.1416 / 2, 3.1416 / 2, 3, 3)
# Assign it to the motor on the mechanism
mech.assignControllers([ctrl])

# Simulate and animate
mech.simulate(20, 2000)
mech.animate()