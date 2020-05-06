
class Node:
	def __init__(self, gene, type, frame, position):
		self.gene = gene
		self.type = type
		self.frame = frame
		self.position = position
	def __hash__(self):
		return hash(repr(self))
	#def __str__(self):
	#	return self.gene + "-" + self.type + "-" + str(self.position)
	def __eq__(self, other):
		return hash(self) == hash(other)
	def __repr__(self):
		"""Compute the string representation of the edge."""
		return "%s(%s,%s,%s,%s)" % (
			self.__class__.__name__,
			repr(self.gene),
			repr(self.type),
			repr(self.frame),
			repr(self.position))
