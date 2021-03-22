from .mechanism import Mechanism
from .dna import DNA
from .mutator import Mutator
from .objective import Objective
from .controller import Controller
from graph_tool.topology import isomorphism as _isomorphism
from numpy.random import RandomState as _RandomState

import json as _json

# The random number generator is global to the library
rnd = _RandomState()


# A function to set the global seed used by the algorithms, to guarantee reproducibility
def setSeed(seed):
    rnd.seed(seed)
    return


def loadDna(filename: str) -> DNA:
    if not filename.endswith('.json'):
        filename = filename + '.json'

    with open(filename, "r") as read_file:
        dna_dict_list = _json.load(read_file)

    dna_list = []

    for dna_dict in dna_dict_list:
        dna_list.append(DNA(dna_dict['incidence_matrix'],
                            dna_dict['edge_labels'],
                            dna_dict['masses'],
                            dna_dict['parameters']))

    if len(dna_list) == 1:
        return dna_list[0]
    else:
        return dna_list


def storeDna(dna: DNA, filename: str):
    if not filename.endswith('.json'):
        filename = filename + '.json'

    try:
        len(dna)
    except TypeError:
        dna = [dna]
    dna_dict_list = []

    for i, d in enumerate(dna):
        # Convert parameters from arrays to lists
        parameters = d.parameters()
        for j in range(len(parameters)):
            parameters[j] = parameters[j].tolist()

        masses = d.masses()
        for j in range(len(masses)):
            masses[j] = masses[j].tolist()

        dna_dict = {"number": i,
                    "incidence_matrix": d.incidenceMatrix.tolist(),
                    "edge_labels": d.edgeLabels().tolist(),
                    "masses": masses,
                    "parameters": parameters}

        dna_dict_list.append(dna_dict)

    with open(filename, "w") as write_file:
        _json.dump(dna_dict_list, write_file, separators=(',', ':'))

    return


# Complexity measure evaluation
def complexity(dna: DNA) -> int:
    rating = 0
    nr_of_labels = 4

    # Allowed edge labels -> possible variations count for the complexity
    for i in range(1, len(dna.edgeLabels())):
        # A label has the complexity rating Nlabels * Nparameters
        rating += nr_of_labels * len(dna.parameters(i))

    # End-effector and ground have complexity zero, guaranteed elements
    for i in range(len(dna.masses())):
        rating += 3

    return rating


# TODO Implementation of the mechanism distance measure
def distance(dna_i: DNA, dna_j: DNA, operator_set) -> float:

    return 0


# checks graph structure isomorphism of the mechanisms
def isomorphism(dna_i: DNA, dna_j: DNA) -> bool:
    # Can only possibly be isomorphic with identical complexity scores -> saves time
    if complexity(dna_i) == complexity(dna_j):
        # If identical complexity, check graph isomorphism
        return _isomorphism(dna_i.toGraph(), dna_j.toGraph())
    else:
        return False