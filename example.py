#! /usr/bin/env python3

from modules.nodes import Node
import phanotate_connect as pc

class Orf:
    def __init__(self, start, stop, frame, *other ):
        self.start = start
        self.stop = stop
        self.frame = frame


f = open("orfs.txt", "r")

# write edges to the graph
for line in f:
	orf = eval(line)
	if orf.frame > 0:
		ret = pc.add_edge(orf.start, orf.stop)
	else:
		pass #print(orf.stop, orf.start)

# find the best path from a source node to a target node
for edge in pc.get_connected():
	print(edge)
