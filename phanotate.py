#!/usr/bin/env python
import os
import sys
import getopt

#from subprocess import Popen, PIPE, STDOUT
import fastpath as fp

from phanotate_modules import file_handling
from phanotate_modules import functions
from phanotate_modules.nodes import Node


#--------------------------------------------------------------------------------------------------#
#                               ARGUMENTS                                                          #
#--------------------------------------------------------------------------------------------------#

args = file_handling.get_args()

#--------------------------------------------------------------------------------------------------#
#                               FILE INPUT                                                         #
#--------------------------------------------------------------------------------------------------#

my_contigs = file_handling.read_fasta(args.infile);
if not my_contigs:
	sys.stdout.write("Error: no sequences found in infile\n")
	sys.exit()

#--------------------------------------------------------------------------------------------------#
#                               MAIN ROUTINE                                                       #
#--------------------------------------------------------------------------------------------------#
for id, seq in my_contigs.items():

	#-------------------------------Find the ORFs----------------------------------------------#
	my_orfs = functions.get_orfs(seq)



	#-------------------------------Create the Graph-------------------------------------------#
	my_graph = functions.get_graph(my_orfs)



	#-------------------------------Run Bellman-Ford-------------------------------------------#
	source = "Node('source','source',0,0)"
	target = "Node('target','target',0," + str(len(seq)+1) + ")"
	# Write edges to the fastpath program, and multiply the weight to not lose decimal places
	for e in my_graph.iteredges():
		ret = fp.add_edge(str(e))

	shortest_path = fp.get_path(source=source, target=target)

	if args.dump:
		[sys.stdout.write(repr(e.source) + "\t" + repr(e.target) + "\t" + str(e.weight*100000) + "\n") for e in my_graph.iteredges()]
		sys.exit()

	
	#-------------------------------Write Output ----------------------------------------------#
	file_handling.write_output(id, args, shortest_path, my_graph, my_orfs)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#

