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


f = open("orfs.txt", "r")

# write edges to the graph
for line in f:
	#orf = eval(line)
	left, right, frame, pstop = line.rstrip().split("\t")
	ret = pc.add_edge(left, right, int(frame), float(pstop))

# find the best path from a source node to a target node
for edge in pc.get_connections():
	print(edge)
