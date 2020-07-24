#! /usr/bin/env python3

from phanotate_modules.nodes import Node
from phanotate_modules import functions
import phanotate_connect as pc

class Orf:
	def __init__(self, left, right, direction, pstop):
		self.left = left
		self.right = right
		self.direction = direction
		self.pstop = pstop


f = open("orfs", "r")

# write edges to the graph
for line in f:
	orf = eval(line)
	ret = pc.add_edge(orf.left, orf.right, orf.direction, orf.pstop)

# find the best path from a source node to a target node
for edge in pc.get_connections():
	print(edge)
