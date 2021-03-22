from pymechlib import Mutator, DNA, rnd
from typing import List

T = Mutator.Types
L = DNA.Labels


# Operators are defined as specific cases of the underlying mutator class


# Always inherit from Mutator, which handles the underlying computations
class AddLink(Mutator):
    isStructural = True
    isConstructive = True

    # Define it's initialization function, as a function of the mutable parameters
    def __init__(self, com, mass, source_element, connection_position):

        # Build the operator using the mutator representation defined in the thesis
        struct = [T.ADD_VERTEX, T.ADD_EDGE]
        pars = [com + [mass], [source_element, -1, L.HINGE] + connection_position]

        # Initialize the underlying mutator c++ class with these parameters
        Mutator.__init__(self, struct, pars)

    # In order to use it with the evolutionary algorithm toolbox it also needs to define a (random) generation function
    @staticmethod
    def factory(dna: DNA):
        all_vertices = range(1, len(dna.masses()) + 1)
        # Randomly choose the element the new link is attached to -> NOT the ground
        source_element = rnd.choice(all_vertices)
        # Determine the connection position
        connection_position = generatePosition()
        # Determine mass and com of the new link
        mass = 1
        com = generatePosition()

        return AddLink(com, mass, source_element, connection_position)

# Rules: Don't allow >1 connection between elements.
#        Don't allow elements to be connected to eachother and both to the ground
#        Filter these out in the generation function


class AddConnection(Mutator):
    isStructural = True
    isConstructive = True

    def __init__(self, vertices, edge_label, parameters):

        struct = [T.ADD_EDGE]
        pars = [vertices + [edge_label] + parameters]
        Mutator.__init__(self, struct, pars)

    @staticmethod
    def factory(dna: DNA):
        all_vertices = range(1, len(dna.masses()) + 1)

        vertex1 = rnd.choice(all_vertices)
        vertex1 = rnd.choice(all_vertices)
        connections1 = dna.linkConnections(vertex1)
        # Determine which vertices are not yet connected to this one
        not_connected = [x for x in all_vertices if x not in connections1 + [vertex1]]

        i = 0

        # Don't overconnect the mechanism -> at least one free state
        if 3*len(dna.masses()) < dna.incidenceMatrix.shape[1]*2 + 3:
            vertex2 = 0
        else:
            if len(not_connected) > 1:
                shared_connection = True
                while shared_connection:
                    # Choose the second vertex from the list of not yet connected vertices
                    vertex2 = rnd.choice(not_connected)
                    connections2 = dna.linkConnections(vertex2)
                    # Check for shared connections, if so restart -> filter out springs somehow
                    shared_connection = any(x in connections1 for x in connections2)
                    i = i + 1
                    if i > 5:
                        vertex2 = 0
                        break
            else:
                try:
                    vertex2 = not_connected[0]
                except IndexError:
                    vertex2 = 0

        # If it is a ground connection it is always a spring connection
        if vertex2 == 0:
            edge_label = L.SPRING
        else:
            edge_label = L.HINGE

        # Randomly generate corresponding parameters
        parameters = randomParameter(edge_label)

        # Safety check, if somehow an element will be connected to itself, halt the generation
        if vertex1 == vertex2:
            return False

        return AddConnection([vertex1, vertex2], edge_label, parameters)


class RemoveLink(Mutator):
    isStructural = True
    isConstructive = False

    def __init__(self, dna: DNA, link_number: int):

        struct = [T.REMOVE_VERTEX]
        pars = [[link_number]]

        connections = dna.linkConnections(link_number)
        for connection in connections:
            struct.append(T.REMOVE_EDGE)
            pars.append([connection])

        Mutator.__init__(self, struct, pars)

    @staticmethod
    def factory(dna: DNA):

        removable_vertex_list = []

        # Only remove links that are free
        for i in range(1, len(dna.masses()) + 1):
            connections = dna.linkConnections(i)
            if len(connections) == 1:
                removable_vertex_list.append(i)

        if not removable_vertex_list:
            return False
        else:
            vertex = rnd.choice(removable_vertex_list)
            return RemoveLink(dna, vertex)


class RemoveConnection(Mutator):
    isStructural = True
    isConstructive = False

    def __init__(self, edgenr: int):

        struct = [T.REMOVE_EDGE]
        pars = [[edgenr]]
        Mutator.__init__(self, struct, pars)

    @staticmethod
    def factory(dna: DNA):

        removable_edge_list = []

        # Only remove connections that can be severed, i.e. excess ground connections
        for i in range(2, len(dna.parameters())):
            connections = dna.edgeConnections(i)
            if connections[0] == 0:
                removable_edge_list.append(i)

        if not removable_edge_list:
            return False
        else:
            edge = rnd.choice(removable_edge_list)
            return RemoveConnection(edge)


class MutateParameter(Mutator):
    isStructural = False
    def __init__(self, element, new_parameters):

        struct = [T.PARAMETER]
        pars = [[element] + new_parameters]
        Mutator.__init__(self, struct, pars)

    @staticmethod
    def factory(dna: DNA):
        element = rnd.choice(len(dna.parameters())) + len(dna.masses())
        if element < len(dna.masses()):
            new_parameters = randomParameter(-1)
        else:
            edge = element - len(dna.masses())
            new_parameters = randomParameter(dna.edgeLabels(edge))

            # Can't modify the original ground connection location
            if edge == 1:
                new_parameters[0:2] = dna.parameters(edge)[0:2].copy()

        # There is a recursive 50% chance a mutation becomes more complex
        if rnd.choice([True, False]):
            return MutateParameter.factory(dna)*MutateParameter(element, new_parameters)

        return MutateParameter(element, new_parameters)


class SimpleParameter(Mutator):
    isStructural = False
    def __init__(self, element, new_parameters):

        struct = [T.PARAMETER]
        pars = [[element] + new_parameters]
        Mutator.__init__(self, struct, pars)

    @staticmethod
    def factory(dna: DNA):
        element = rnd.choice(len(dna.parameters())) + len(dna.masses())
        if element < len(dna.masses()):
            new_parameters = [0]*3
            new_parameters[rnd.choice(3)] = rnd.uniform(-0.2, 0.2)
        else:
            edge = element - len(dna.masses())
            new_parameters = [0]*len(dna.parameters(edge))

            if edge == 0:
                new_parameters[rnd.choice(len(dna.parameters(edge)) - 1)] = rnd.uniform(-0.2, 0.2)
            else:
                new_parameters[rnd.choice(len(dna.parameters(edge)))] = rnd.uniform(-0.2, 0.2)

            # Can't modify the original ground connection location
            if edge == 1:
                new_parameters[0:2] = dna.parameters(edge)[0:2].copy()

        return SimpleParameter(element, new_parameters)


class ReLabel(Mutator):
    isStructural = True
    isConstructive = True
    def __init__(self, edge, new_label, new_parameters):

        struct = [T.RELABEL]
        pars = [[edge, new_label] + new_parameters]
        Mutator.__init__(self, struct, pars)

    @staticmethod
    def factory(dna: DNA):

        try:
            edge = rnd.choice(len(dna.edgeLabels()) - 1) + 1
        except ValueError:
            return False
        drawn_label = dna.edgeLabels(edge)

        if drawn_label == L.SPRING:
            return False

        new_label = rnd.choice([x for x in [L.HINGE, L.TORSION_SPRING] if x not in [drawn_label]])
        new_parameters = randomParameter(new_label)

        # The position of the element remains the same -> needs explict copy
        new_parameters[0:2] = dna.parameters(edge)[0:2].copy()

        return ReLabel(edge, new_label, new_parameters)


# Transform
class Transform(Mutator):
    isStructural = False
    def __init__(self, rotation: float, scale: float, offset: List[float]):

        struct = [T.TRANSFORM]
        pars = [[rotation, scale] + offset]
        Mutator.__init__(self, struct, pars)

    @staticmethod
    def factory(dna: DNA):

        rotation = rnd.normal(0, 1)
        # Rescale from 50% to 150%
        scale = 0.5 + 1*rnd.random()
        offset = [0, 0]

        return Transform(rotation, scale, offset)


# Operator that connects the end-effector to a different link
class MoveEndEffector(Mutator):
    isStructural = True
    isConstructive = False
    def __init__(self, link_number: int):

        # Randomly select
        struct = [T.MOVE_EFFECTOR]
        pars = [[link_number]]
        Mutator.__init__(self, struct, pars)

    @staticmethod
    def factory(dna: DNA):

        connections = dna.edgeConnections(0)
        try:
            link_number = rnd.choice([x for x in range(1, len(dna.masses()) + 1) if x not in [connections[0]]])
        except ValueError:
            return False

        return MoveEndEffector(link_number)


def generatePosition():

    # What distribution to use?
    return [rnd.uniform(-0.25, 0.25), rnd.uniform(-0.25, 0.25)]


def randomParameter(label: L):

    if label == L.HINGE:
        return generatePosition()
    if label == L.SPRING:
        return generatePosition() + generatePosition() + [rnd.uniform(-0.25, 0.25), rnd.uniform(-1, 1)]
    if label == L.TORSION_SPRING:
        return generatePosition() + [rnd.uniform(-0.25, 0.25), rnd.uniform(-1.25, 1.25)]
    if label == L.HINGE_MOTOR:
        return generatePosition() + [0]
    if label == L.END_EFFECTOR:
        return generatePosition() + [0]
    if label == -1:
        return generatePosition() + [rnd.uniform(-0.25, 0.25)]

