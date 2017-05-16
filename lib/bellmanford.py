#!/usr/bin/python

from decimal import Decimal

class BellmanFord:
    """The Bellman-Ford algorithm for the shortest path problem.
    
    Attributes
    ----------
    graph : input directed weighted graph
    parent : dict with nodes (shortest path tree)
    distance : dict with nodes (distances to source node)
    source : node
    """

    def __init__(self, graph):
        """The algorithm initialization.
        
        Parameters
        ----------
        graph : directed weighted graph
        """
        if not graph.is_directed():
            raise ValueError("the graph is not directed")
        self.graph = graph
        # Shortest path tree as a dictionary.
        self.parent = dict(((node, None) for node in self.graph.iternodes()))
        self.distance = dict(((node, Decimal('Infinity')) for node in self.graph.iternodes()))
        self.source = None

    def run(self, source):
        """Finding shortest paths from the source.
        
        Parameters
        ----------
        source : node
        """
        self.source = source
        self.distance[source] = Decimal('0')
        for step in xrange(self.graph.v()-1):   # |V|-1 times
	    flag = 1
            for edge in self.graph.iteredges():   # O(E) time
		if self.distance[edge.target] > self.distance[edge.source] + edge.weight:
		    self.distance[edge.target] = self.distance[edge.source] + edge.weight
                    self.parent[edge.target] = edge.source
                    flag = 0
            if(flag):
		break
        # Check for negative cycles.
        for edge in self.graph.iteredges():   # O(E) time
            if self.distance[edge.target] > self.distance[edge.source] + edge.weight:
		print "target", edge.target
		print "source", edge.source
		print "weight", edge.weight
		print "parent", self.parent[edge.source]
		print "dist_source", self.distance[edge.source]
		print "dist_target", self.distance[edge.target]
                raise ValueError("negative cycle detected")

    def path(self, target):
        """Construct a path from source to target."""
        if self.source == target:
            return [self.source]
        elif self.parent[target] is None:
            raise ValueError("no path to target")
        else:
            return self.path(self.parent[target]) + [target]

# EOF
