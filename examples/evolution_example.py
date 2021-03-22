import numpy as np
import pymechlib as mch
from matplotlib import pyplot as plt
from pymechlib.evolution import operators as op
import pymechlib.evolution as evo

#mch.setSeed(80085)
# Alias the mechanism edge labels
L = mch.DNA.Labels

# Defines a simple pendulum mechanism as seed mechanism
inc_matrix = [[0, 1],
              [1, 1],
              [1, 0]]
masses = [[-0.25, -0.25, 5]]
edge_labels = [L.END_EFFECTOR, L.HINGE]
parameters = [[1, 1, 0.1], [0, 0]]

seed_mech = mch.Mechanism(inc_matrix, edge_labels, masses, parameters)

operator_set = [#op.SimpleParameter,
                op.MutateParameter,
                op.Transform,
                op.ReLabel,
                op.AddLink,
                op.AddConnection,
                op.RemoveLink,
                op.RemoveConnection,
                op.MoveEndEffector]

objective = evo.PickAndPlace([1, 2], [-1, 0])

#Tree evolution is the core of the library
tree, fitness = evo.treeEvolution(objective, seed_mech, operator_set, leaves=32, population_size=64,
                                    generations=16, epochs=50, initial_spread=5, simulation_steps=200, simulation_time=10)

# Plots the fitness over epochs
np.save("fitness.npy", fitness)
plt.plot(fitness, '.-')
plt.legend(["Maximum fitness", "Average of all nodes", "Average of active nodes"])
plt.xlabel("Epochs")
plt.ylabel("Fitness")
plt.ylim((-3, 0))
plt.grid('on')
plt.show()

# Plots the mechanism space graph
tree.plot_span(40, "span.pdf")
tree.plot_core(40, "min.pdf")

fitness_list = []
mechanism_list = []
controller_dict = {}

# Visualize some of the mechanisms, using a higher simulation framerate for display
for i in range(10):
    tree.structures[tree.ranking[i]].best_member.simulate(10, 600)
    tree.structures[tree.ranking[i]].best_member.animate(objective=objective, trace=True)

# Store the mechanisms and evolution data
for node in tree.structures:
    mechanism_list.append(node.best_member.dna)
    fitness_list.append(node.best_fitness)

mch.storeDna(mechanism_list, "examples.json")
np.save("examples.npy", fitness_list)
tree.save("tree.gt")