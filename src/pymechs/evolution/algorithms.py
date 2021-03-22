# This module provides helper functions for the evolutionary algorithm
from graph_tool import GraphView
import numpy as np
from typing import List, Type
from os import cpu_count
import time

from pymechlib import DNA, Mechanism, Mutator, Objective, complexity, isomorphism, Controller, rnd
from pymechlib.cmechlib import Threader
from graph_tool.all import Graph, graph_draw, radial_tree_layout, all_shortest_paths
from matplotlib.pyplot import cm

# The threader is shared and global to the evolution module
from pymechlib.control import Pid

print("Initialized pymechlib threader with " + str(cpu_count()*2) + " threads")
threader = Threader(cpu_count()*2)

# Prints iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 30, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)

    if iteration == total:
        print("\n")

# Picks and initializes a random operator from the base set
def randomMutator(dna: DNA, operator_set: List[Type[Mutator]], operator_probabilities=None) -> Type[Mutator]:
    # Draw a random operator to be applied from the list, only if > 1 options to reduce random number generation cost
    if len(operator_set) > 1:
        operator_type = rnd.choice(operator_set)
    else:
        operator_type = operator_set[0]

    # Randomly initialize this operator, using the (generally dna dependent) generation function defined by it
    try:
        operator = operator_type.factory(dna)
    except AttributeError as err:
        err.message = err.message + "An operator can only be randomly generated if it provides a factory function"
        raise
    if not operator:
        return randomMutator(dna, operator_set, operator_probabilities)

    return operator


def sortOperators(operator_set: List[Type[Mutator]]):

    operators = {
        "constructive": [],
        "destructive": [],
        "parametric": []
    }

    # Sort the operator set by type
    for i in range(len(operator_set)):
        if operator_set[i].isStructural:
            if operator_set[i].isConstructive:
                # It is a construction operator
                operators['constructive'].append(operator_set[i])
            else:
                # It is a destruction operator
                operators['destructive'].append(operator_set[i])
        else:
            # It is a parametric operator
            operators['parametric'].append(operator_set[i])
    
    return operators


# The fitness tree is the graph of the mechanism search space
# - empty nodes still have a fitness value
# - population size per node, overall population
# - the evolutionary population 'lives' in this graph

# TODO More deeply tie the fitness tree to the graph-tools graph structure
class FitnessNode:

    def __init__(self, root: Mechanism, population_size: int):

        self.best_member = root.clone()
        self.best_fitness = -1e20

        self.fitness = np.zeros(population_size)
        self.population = [root.clone()]

        for i in range(population_size - 1):
            self.population.append(root.clone())

    def evaluate(self, objective, simulation_time, simulation_steps):
        threader.evaluate(self.population, objective, simulation_time, simulation_steps)
        # Update the fitness list and best member
        for n, mech in enumerate(self.population):
            self.fitness[n] = mech.fitness

        ranking = np.argsort(-self.fitness)
        self.best_member = self.population[ranking[0]].clone()
        self.best_fitness = self.fitness[ranking[0]]

    def localEvolution(self, objective, operators, generations, simulation_time, simulation_steps):

        for j in range(generations):
            # Evaluate the population of the node on the objective (multithreaded)
            self.evaluate(objective, simulation_time, simulation_steps)

            for k in range(1, len(self.population)):
                mutator = randomMutator(self.best_member.dna, operators)
                self.population[k] = Mechanism(mutator.eval(self.best_member.dna))

            # Always keep the best member in the population
            self.population[0] = self.best_member.clone()

    def setExtinct(self):
        self.population.clear()

    def setActive(self, root, population_size: int):
        # Wake the node up again and assign a population to it
        if self.isEmpty():
            self.__init__(self.best_member.clone(), population_size)
        # Insert the crossover member into the population
        self.population[1] = root.clone()

    def isEmpty(self) -> bool:
        return not self.population


class FitnessTree(Graph):

    def __init__(self, root, population_size):

        Graph.__init__(self, directed=True)
        self.structures = self.new_vp("object")
        self.rank_position = self.new_vp("int")
        self.complexity = self.new_vp("int")
        self.color = self.new_vp("vector<double>")
        self.operators = self.new_ep("object")
        self.labels = self.new_ep("string")
        self.edge_width = self.new_ep("double")

        self.core_vertex = self.new_vp("bool")
        self.core_edge = self.new_ep("bool")

        # Define the root vertex
        self.add_vertex()
        self.structures[0] = FitnessNode(root, population_size)
        self.complexity[0] = complexity(root.dna)
        self.best_member = root
        self.ranking = []
        self.core_only = False

    # Wrapper for add_vertex() that checks for isomorphisms and double edges
    def addNode(self, mutator: Type[Mutator], source_node, population_size: int):
        # Generate the new mechanism with the mutator, build a new node representing the new structure
        new_node_dna = mutator.eval(self.structures[source_node].best_member.dna)
        isomorphic = False
        # Check if the new node is isomorphic to an existing node, if the case, the new node is the old node.
        for node in self.vertices():
            if isomorphism(DNA(new_node_dna), DNA(self.structures[node].best_member.dna)):
                self.structures[node].setActive(Mechanism(new_node_dna), population_size)
                isomorphic = True
                break

        if not isomorphic:
            new_node = FitnessNode(Mechanism(new_node_dna), population_size)
            # Add a vertex for the new label
            node = self.add_vertex()
            self.structures[node] = new_node
            self.complexity[node] = complexity(new_node_dna)
            # Store the operator leading to this node as the edge label
            e = self.add_edge(source_node, node)
            self.operators[e] = mutator
            self.labels[e] = mutator.__class__.__name__
            self.edge_width[e] = 0.5

            return node

        # If it's isomorphic it defines just another edge, check for duplicates
        else:
            existing_edge = self.edge(source_node, node)
            if existing_edge:
                if self.labels[existing_edge] != str(mutator.__class__.__name__):
                    e = self.add_edge(source_node, node)
                    self.operators[e] = mutator
                    self.labels[e] = mutator.__class__.__name__
                    self.edge_width[e] = 0.5
            else:
                e = self.add_edge(source_node, node)
                self.operators[e] = mutator
                self.labels[e] = mutator.__class__.__name__
                self.edge_width[e] = 0.5

            return node

    def plot(self, output=None):

        # Map the ranking to colors for easier interpretation
        color = cm.winter(np.linspace(0, 1, self.num_vertices()))

        for n, v in enumerate(self.ranking):
            try:
                self.rank_position[v] = n + 1
                self.color[v] = color[-n-1]
            except OverflowError:
                self.color[v] = [1, 0.2, 0.2, 1]

        self.color[0] = [1, 0.2, 0.2, 1]

        if self.core_only is True:
            self.set_edge_filter(self.core_edge)
            self.set_vertex_filter(self.core_vertex)
            self.core_only = False

        # Plot the graph
        graph_draw(self,
                   pos=radial_tree_layout(self, self.vertex(0)),
                   vertex_color=[0, 0, 0, 1],
                   vertex_fill_color=self.color,
                   vertex_pen_width=1,
                   vertex_text=self.rank_position,
                   edge_pen_width=self.edge_width,
                   #edge_text=self.labels,
                   edge_text_out_color=[1, 0, 0, 1],
                   output=output)

    def plot_span(self, nr_top_nodes, output=None):

        for i in range(nr_top_nodes):
            for elist in all_shortest_paths(self, self.vertex(0), self.vertex(self.ranking[i]), edges=True):
                for e in elist:
                    self.edge_width[e] = 2

        self.plot(output)

    def plot_core(self, nr_top_nodes, output=None):

        for i in range(nr_top_nodes):
            for path in all_shortest_paths(self, self.vertex(0), self.vertex(self.ranking[i])):
                for v in path:
                    self.core_vertex[v] = True
            for path in all_shortest_paths(self, self.vertex(0), self.vertex(self.ranking[i]), edges=True):
                for e in path:
                    self.core_edge[e] = True

        # @TODO this global is garbage, needs a proper fix
        self.core_only = True
        self.plot(output)

        

def treeEvolution(objective: Objective, seed_mechanism: Mechanism, operator_set: List[Type[Mutator]], leaves: int,
                  generations: int, population_size: int, epochs: int, simulation_time=10, simulation_steps=200,
                  initial_spread=10, seed_controller=None) -> FitnessTree:
    # TODO Refactor
    # Local controller optimization on top of parameters
    # Global structure

    # Convert the operator set to sorted dict
    operators = sortOperators(operator_set)

    ##########################################

    # Keep track of time to measure performance
    tic = time.perf_counter()

    # Initialize the fitness tree with the seed mechanism as root
    tree = FitnessTree(seed_mechanism, population_size)
    overall_best_fitness = -9
    overall_simulations = 0
    fitness_over_time = []

    print("Growing initial tree...", end=" ")

    # Generate the initial unique structures
    for i in range(leaves):
        # Generate a series of constructive structural mutators from the set, size depends on the desired initial spread
        # This series of operations defines a new mechanism, according to the operator-mechanism equivalency
        # Add a node for this mechanism structure, taking into account and mapping "skipped" nodes as well

        # first modification connects to root
        mutator = randomMutator(seed_mechanism.dna, operators['constructive'])
        new_node = tree.addNode(mutator, tree.vertex(0), population_size)

        # Then a chain of operations mapping the path taken to the final node from the root
        for j in range(initial_spread - 2):
            mutator = randomMutator(tree.structures[new_node].best_member.dna, operators['constructive'])
            new_node = tree.addNode(mutator, tree.vertex(new_node), population_size)

        mutator = randomMutator(tree.structures[new_node].best_member.dna, operators['constructive'])
        tree.addNode(mutator, tree.vertex(new_node), population_size)

    print("generated " + str(tree.num_vertices()) + " initial nodes \n")

    # Start the core algorithm loop
    for i in range(epochs):
        active_fitness_list = []
        fitness_list = [node.best_fitness for node in tree.structures]
        for vertex in tree.vertices():
            node = tree.structures[vertex]
            # Print a progress bar for realtime info updates in the terminal
            printProgressBar(int(vertex), tree.num_vertices() - 1, 
                             prefix="Epoch " + str(i + 1) + "/" + str(epochs),
                             suffix="|| Max fitness " + "{:.5f}".format(overall_best_fitness))

            if not node.isEmpty():
                node.localEvolution(objective, operators['parametric'], generations, simulation_time, simulation_steps)
                fitness_list[int(vertex)] = node.best_fitness
                # Track statistics
                active_fitness_list.append(node.best_fitness)
                overall_simulations += population_size*generations
                if node.best_fitness > overall_best_fitness:
                    overall_best_fitness = node.best_fitness

        # What N nodes performed best? -> order by fitness
        tree.ranking = np.argsort(-np.array(fitness_list))
        tree.best_member = tree.vertex(tree.ranking[0])
        fitness_over_time.append([overall_best_fitness, np.average(fitness_list), np.average(active_fitness_list)])
        #print("Best node is " + str(tree.ranking[0]) + "/" + str(tree.num_vertices()) + "\n")

        # TODO Clean up this hacky part but it works -> Only remove nodes if more active nodes than specified leaves
        if i == epochs - 1:
            break

        # Every epoch structural mutations occur generating additional nodes, poorly performing nodes are made extinct
        k = 0
        while k < tree.num_vertices() - leaves:
            while tree.structures[tree.ranking[-k - 1]].isEmpty() and k < tree.num_vertices() - leaves:
                k += 1
            tree.structures[tree.ranking[-k - 1]].setExtinct()

        for j in range(leaves):
            mutator = randomMutator(tree.structures[tree.ranking[j]].best_member.dna,
                                    operators['constructive'] + operators['destructive'])

            tree.addNode(mutator, tree.vertex(tree.ranking[j]), population_size)
    

    toc = time.perf_counter()
    
    # Print statistics
    print("\nSimulated " + str(overall_simulations) + " mechanisms")
    print("Processed " + str(tree.num_vertices(ignore_filter=True)) + " unique nodes")
    print("Processing took " + "{:.0f}".format(toc - tic) + " seconds")

    return tree, fitness_over_time


def islandEvolution(objective: Objective, seed_mechanism: Mechanism, operator_set: List[Type[Mutator]], islands: int,
                    generations: int, population_size: int, epochs: int, simulation_time=10, simulation_steps=200,
                    initial_spread=10) -> List[FitnessNode]:

    # Keep track of time to measure performance
    tic = time.perf_counter()

    overall_simulations = 0
    fitness_over_time = []
    overall_best_fitness = -9
    node_list = []

    # Generate the initial islands
    for i in range(islands):
        new_dna = seed_mechanism.dna
        for j in range(rnd.choice(initial_spread)):
            new_dna = randomMutator(new_dna, operator_set).eval(new_dna)

        node_list.append(FitnessNode(Mechanism(new_dna), population_size))

    # Start the core algorithm loop
    for i in range(epochs):
        fitness_list = np.empty(len(node_list))
        fitness_list.fill(-1e20)
        for n, node in enumerate(node_list):
            node.localEvolution(objective, operator_set, generations, simulation_time, simulation_steps)
            fitness_list[n] = node.best_fitness
            overall_simulations += population_size*generations
            if node.best_fitness > overall_best_fitness:
                    overall_best_fitness = node.best_fitness

            # Print a progress bar for realtime info updates in the terminal
            printProgressBar(n, islands - 1, 
                            prefix="Epoch " + str(i + 1) + "/" + str(epochs),
                            suffix="|| Max fitness " + "{:.5f}".format(overall_best_fitness))

        fitness_over_time.append(overall_best_fitness)

        # Handle islands crossover
        if i == epochs - 1:
            break

        # What nodes are considered neighbours?
        for n, node in enumerate(node_list):
            # Copy the best member of neighbouring nodes into the node population
            node.population[1] = node_list[n-1].best_member.clone()

    toc = time.perf_counter()

    print("\nSimulated " + str(overall_simulations) + " mechanisms")
    print("Processing took " + "{:.0f}".format(toc - tic) + " seconds")

    return node_list, fitness_list

def localEvolution(objective: Objective, seed_mechanism: Mechanism, operator_set: List[Type[Mutator]], 
                   generations: int, population_size: int, simulation_time=10, simulation_steps=200,
                   initial_spread=10):

    # Keep track of time to measure performance
    tic = time.perf_counter()

    overall_simulations = 0
    best_member = seed_mechanism
    fitness_over_time = []
    overall_best_fitness = -9
    population = []

    # Generate the initial population
    for i in range(population_size):
        new_dna = seed_mechanism.dna
        for j in range(rnd.choice(initial_spread)):
            new_dna = randomMutator(new_dna, operator_set).eval(new_dna)

        population.append(Mechanism(new_dna))

    # Start the core algorithm loop
    fitness = np.empty(population_size)
    fitness.fill(-1e20)
    for i in range(generations):
        # Evaluate the population on the objective (multithreaded)
        threader.evaluate(population, objective, simulation_time, simulation_steps)
        # Update the fitness list and best member
        for n, mech in enumerate(population):
            fitness[n] = mech.fitness

        ranking = np.argsort(-fitness)
        best_member = population[ranking[0]].clone()
        best_fitness = fitness[ranking[0]]

        overall_simulations += population_size
        if best_fitness > overall_best_fitness:
                    overall_best_fitness = best_fitness

        # Print a progress bar for realtime info updates in the terminal
        printProgressBar(i, generations - 1, 
                        prefix="Generation " + str(i + 1) + "/" + str(generations),
                        suffix="|| Max fitness " + "{:.5f}".format(overall_best_fitness))

        fitness_over_time.append(overall_best_fitness)

        population.clear()

        for k in range(0, population_size):
            mutator = randomMutator(best_member.dna, operator_set)
            population.append(Mechanism(mutator.eval(best_member.dna)))

        population[0] = best_member

    toc = time.perf_counter()

    print("\nSimulated " + str(overall_simulations) + " mechanisms")
    print("Processing took " + "{:.0f}".format(toc - tic) + " seconds")

    return best_member, fitness_over_time
