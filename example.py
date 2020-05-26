#! /usr/bin/env python3

from phanotate_modules.nodes import Node
from phanotate_modules import functions
import phanotate_connect as pc

class Orf:
    def __init__(self, left, right, pstop):
        self.left = left
        self.right = right
        self.pstop = pstop


f = open("orfs.txt", "r")

# write edges to the graph
for line in f:
	orf = eval(line)
	ret = pc.add_edge(orf.left, orf.right, orf.pstop)

# find the best path from a source node to a target node
for edge in pc.get_connected():
	if edge[0] < edge[1]:
		print(edge[0], edge[1], functions.score_gap(
