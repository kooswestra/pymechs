from pymechlib.cmechlib import _DNA
from graph_tool import Graph

class DNA(_DNA):

    def toGraph(self) -> Graph:
        graph = Graph(directed=False)
        inc_matrix = self.incidenceMatrix

        graph.ep.labels = graph.new_edge_property("int")
        graph.add_vertex(inc_matrix.shape[0])

        for n, label in enumerate(self.edgeLabels()):
            connections = self.edgeConnections(n)
            e = graph.add_edge(connections[0], connections[1])
            graph.ep.labels[e] = label

        return graph