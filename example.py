#! /usr/bin/env python3

from phanotate_modules.nodes import Node
import phanotate_connect as pc

class Orf:
    def __init__(self, left, right, weight):
        self.left = left
        self.right = right
        self.weight = weight


f = open("orfs.txt", "r")

# write edges to the graph
for line in f:
	orf = eval(line)
	ret = pc.add_edge(orf.left, orf.right, orf.weight)

exit()
# find the best path from a source node to a target node
for edge in pc.get_connected():
	print(edge)
