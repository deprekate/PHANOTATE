#!/usr/bin/python

class Edge:
    def __init__(self, source, target, weight):
        """Load up a directed edge instance.
        
        Parameters
        ----------
        source : starting node
        target : ending node
        weight : number
        """
        self.source = source
        self.target = target
        self.weight = weight

    def __str__(self):
        """Compute the string representation of the edge."""
        return "%s\t%s\t%s" % (
                repr(self.source),
                repr(self.target),
				str(self.weight*1000)
                )

    def __repr__(self):
        """Compute the string representation of the edge."""
        return "%s(%s, %s, %s)" % (
                self.__class__.__name__,
                repr(self.source),
                repr(self.target),
                repr(self.weight))

    def __cmp__(self, other):
        """Comparing of edges (the weight first)."""
        # Check weights.
        if self.weight > other.weight:
            return 1
        if self.weight < other.weight:
            return -1
        # Check the first node.
        if self.source > other.source:
            return 1
        if self.source < other.source:
            return -1
        # Check the second node.
        if self.target > other.target:
            return 1
        if self.target < other.target:
            return -1
        return 0

    def __hash__(self):
        """Hashable edges."""
        return hash(repr(self))

    def __invert__(self):
        """Return the edge with the opposite direction."""
        return self.__class__(self.target, self.source, self.weight)

    inverted = __invert__


class UndirectedEdge(Edge):
    """The class defining an undirected edge."""

    def __init__(self, source, target, weight=1):
        """Load up an edge instance."""
        if source > target:
            self.source = target
            self.target = source
        else:
            self.source = source
            self.target = target
        self.weight = weight

    def __invert__(self):
        """The edge direction is not defined."""
        return self


# EOF
