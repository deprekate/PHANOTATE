#! /usr/bin/env python3

import phanotate_connect as pc


f = open("orfs.txt", "r")

# write edges to the graph
for line in f:
	left, right, frame, pstop = line.rstrip().split("\t")
	ret = pc.add_edge(left, right, int(frame), float(pstop))

# find the best path from a source node to a target node
for edge in pc.get_connections():
	print(edge)
