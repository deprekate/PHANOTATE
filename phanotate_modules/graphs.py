#!/usr/bin/python

import random
from .edges import Edge


class Graph(dict):
    """The class defining a graph.
    
    Nodes can be numbers, strings, or any hashable objects.
    We would like to compare nodes.
    
    An exemplary graph structure:
    {"A": {"B": Edge("A", "B", 1), "C": Edge("A", "C", 2)}, 
    "B": {"C": Edge("B", "C", 3), "D": Edge("B", "D", 4)}, 
    "C": {"D": Edge("C", "D", 5)}, 
    "D": {"C": Edge("D", "C", 6)}, 
    "E": {"C": Edge("E", "C", 7)}, 
    "F": {}}
    """

    def __init__(self, n=0, directed=False):
        """Load up a Graph instance.
        
        Parameters
        ----------
        n : int (positive; not used, for compatibility only)
        directed : bool, optional (default=False)
        """
        self.n = n
        self.directed = directed

    def is_directed(self):
        """Test if the graph is directed."""
        return self.directed

    def v(self):
        """Return the number of nodes (the graph order)."""
        return len(self)

    def e(self):
        """Return the number of edges in O(V) time."""
        edges = sum(len(self[node]) for node in self)
        return (edges if self.is_directed() else edges / 2)

    def add_node(self, node):
        """Add a node to the graph."""
        if node not in self:
            self[node] = dict()

    def has_node(self, node):
        """Test if a node exists."""
        return node in self
    
    def del_node(self, node):
        """Remove a node from the graph (with edges)."""
        # The dictionary changes size during iteration.
        for edge in list(self.iterinedges(node)):
            self.del_edge(edge)
        if self.is_directed():
            for edge in list(self.iteroutedges(node)):
                self.del_edge(edge)
        del self[node]

    def add_edge(self, edge):
        """Add an edge to the graph (missing nodes are created)."""
        if edge.source == edge.target:
            raise ValueError("loops are forbidden:"+str(edge.source)+" "+str(edge.target))
        self.add_node(edge.source)
        self.add_node(edge.target)
        if edge.target not in self[edge.source]:
            self[edge.source][edge.target] = edge
        else:
            raise ValueError("parallel edges are forbidden")
        if not self.is_directed():
            if edge.source not in self[edge.target]:
                self[edge.target][edge.source] = ~edge
            else:
                raise ValueError("parallel edges are forbidden")

    def del_edge(self, edge):
        """Remove an edge from the graph."""
        del self[edge.source][edge.target]
        if not self.is_directed():
            del self[edge.target][edge.source]

    def has_edge(self, edge):
        """Test if an edge exists (the weight is not checked)."""
        return edge.source in self and edge.target in self[edge.source]

    def weight(self, edge):
        """Return the edge weight or zero."""
        if edge.source in self and edge.target in self[edge.source]:
            return self[edge.source][edge.target].weight
        else:
            return 0

    def iternodes(self):
        """Generate the nodes from the graph on demand."""
        return self.keys()

    def iteradjacent(self, source):
        """Generate the adjacent nodes from the graph on demand."""
        return self[source].keys()

    def iteroutedges(self, source):
        """Generate the outedges from the graph on demand."""
        for target in self[source]:
            yield self[source][target]

    def iterinedges(self, source):
        """Generate the inedges from the graph on demand."""
        if self.is_directed():   # O(V) time
            for target in self.iternodes():
                if source in self[target]:
                    yield self[target][source]
        else:
            for target in self[source]:
                yield self[target][source]

    def iteredges(self):
        """Generate the edges from the graph on demand."""
        for source in self.iternodes():
            for target in self[source]:
                if self.is_directed() or source < target:
                    yield self[source][target]

    def show(self):
        """The graph presentation."""
        for source in self.iternodes():
            print(source, ":",)
            for edge in self.iteroutedges(source):
                if edge.weight == 1:
                    print(edge.target,)
                else:
                    print("%s(%s)" % (edge.target, edge.weight),)
            print("\n")

    def copy(self):
        """Return the graph copy."""
        new_graph = Graph(n=self.n, directed=self.directed)
        for node in self.iternodes():
            new_graph[node] = dict(self[node])
        return new_graph

    def transpose(self):
        """Return the transpose of the graph."""
        new_graph = Graph(n=self.n, directed=self.directed)
        for node in self.iternodes():
            new_graph.add_node(node)
        for edge in self.iteredges():
            new_graph.add_edge(~edge)
        return new_graph

    def degree(self, source):
        """Return the degree of the node in the undirected graph."""
        if self.is_directed():
            raise ValueError("the graph is directed")
        return len(self[source])

    def outdegree(self, source):
        """Return the outdegree of the node."""
        return len(self[source])

    def indegree(self, source):
        """Return the indegree of the node."""
        if self.is_directed():   # O(V) time
            counter = 0
            for target in self.iternodes():
                if source in self[target]:
                    counter = counter + 1
            return counter
        else:                   # O(1) time
            return len(self[source])

    def __eq__(self, other):
        """Test if the graphs are equal."""
        if self.is_directed() is not other.is_directed():
            #print "directed and undirected graphs"
            return False
        if self.v() != other.v():
            #print "|V1| != |V2|"
            return False
        for node in self.iternodes():   # O(V) time
            if not other.has_node(node):
                #print "V1 != V2"
                return False
        if self.e() != other.e():   # inefficient, O(E) time
            #print "|E1| != |E2|"
            return False
        for edge in self.iteredges():   # O(E) time
            if not other.has_edge(edge):
                #print "E1 != E2"
                return False
            if edge.weight != other.weight(edge):
                return False
        return True

    def __ne__(self, other):
        """Test if the graphs are not equal."""
        return not self == other

    def add_graph(self, other):
        """Add a graph to this graph (the current graph is modified)."""
        if self.is_directed() is not other.is_directed():
            raise ValueError("directed vs undirected")
        for node in other.iternodes():
            self.add_node(node)
        for edge in other.iteredges():
            self.add_edge(edge)

    def save(self, file_name, name="Graph"):
        """Export the graph to the adjacency list format with comments."""
        afile = open(file_name, "w")
        afile.write("# NAME=%s\n" % name)
        afile.write("# DIRECTED=%s\n" % self.directed)
        afile.write("# V=%s\n" % self.v())
        afile.write("# E=%s\n" % self.e())
        for edge in self.iteredges():
            if edge.weight == 1:
                afile.write("%s %s\n" % (edge.source, edge.target))
            else:
                afile.write("%s %s %s\n" % (edge.source, edge.target, edge.weight))
        afile.close()

    @classmethod
    def load(cls, file_name):
        """Import the graph from the adjacency list format with comments."""
        afile = open(file_name, "r")
        n = 1
        is_directed = False
        for line in afile:
            if line[0] == "#":
                if "# NAME=" in line:
                    name = line[7:-1]
                elif line == "# DIRECTED=False\n":
                    is_directed = False
                elif line == "# DIRECTED=True\n":
                    is_directed = True
                elif "# V=" in line:
                    n = int(line[4:-1])
                else:   # ignore other
                    graph = cls(n, is_directed)
            else:
                #alist = [int(x) for x in line.split()]
                #alist = [eval(x) for x in line.split()]
                alist = line.split()
                if len(alist) == 3:
                    alist[-1] = eval(alist[-1])
                graph.add_edge(Edge(*alist))
        afile.close()
        return graph

    def save_lgl(self, file_name="graph.lgl"):
        """Export the graph to the adjacency list format (LGL)."""
        if self.is_directed():
            raise ValueError("the graph is directed")
        afile = open(file_name, "w")
        for edge in self.iteredges():
            if edge.weight == 1:
                afile.write("%s %s\n" % (edge.source, edge.target))
            else:
                afile.write("%s %s %s\n" % (edge.source, edge.target, edge.weight))
        afile.close()

    def save_ncol(self, file_name="graph.ncol"):
        """Export the graph to the labelled edge list format (NCOL)."""
        if self.is_directed():
            raise ValueError("the graph is directed")
        afile = open(file_name, "w")
        for node in self.iternodes():
            afile.write("# %s\n" % str(node))
            for edge in self.iteroutedges(node):
                if edge.source < edge.target and edge.weight == 1:
                    afile.write("%s\n" % str(edge.target))
                elif edge.source < edge.target:
                    afile.write("%s %s\n" % (edge.target, edge.weight))
        afile.close()

# EOF
