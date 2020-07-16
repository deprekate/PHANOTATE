#!/usr/bin/env python
import os
import sys
import getopt

#from subprocess import Popen, PIPE, STDOUT
import fastpathz as fz

from phanotate_modules import file_handling
from phanotate_modules import functions
from phanotate_modules.nodes import Node
from phanotate_modules.orfs import Orfs


#--------------------------------------------------------------------------------------------------#
#                               ARGUMENTS                                                          #
#--------------------------------------------------------------------------------------------------#

args = file_handling.get_args()

#--------------------------------------------------------------------------------------------------#
#                               FILE INPUT                                                         #
#--------------------------------------------------------------------------------------------------#

base_trans = str.maketrans('SBVDEFHIJKLMNOPQRUWXYZ','GGGAAAAAAAAAAAAAAAAAAA')
my_contigs = file_handling.read_fasta(args.infile, base_trans)
if not my_contigs:
	sys.stdout.write("Error: no sequences found in infile\n")
	sys.exit()

#--------------------------------------------------------------------------------------------------#
#                               MAIN ROUTINE                                                       #
#--------------------------------------------------------------------------------------------------#
for id, seq in my_contigs.items():

	contig_orfs = Orfs(**vars(args))

	#-------------------------------Find the ORFs----------------------------------------------#

	contig_orfs.parse_contig(seq)
	
	for orf in contig_orfs.iter_orfs():
		start = orf.start
		stop = orf.stop
		if(orf.frame > 0):
			print(start, stop)
			gcfp = contig_orfs.gc_frame_plot
			for base in range(start, stop, 3):
				print(gcfp.min_frame_at(base), gcfp.max_frame_at(base))
			exit()

	exit()

	contig_orfs.score_rbs_sites()

	contig_orfs.calculate_weights(my_orfs)


	#-------------------------------Create the Graph-------------------------------------------#
	my_graph = functions.get_graph(my_orfs)



	#-------------------------------Run Bellman-Ford-------------------------------------------#
	source = "Node('source','source',0,0)"
	target = "Node('target','target',0," + str(len(seq)+1) + ")"
	# Write edges to the fastpath program, and multiply the weight to not lose decimal places
	fz.empty_graph()
	for e in my_graph.iteredges():
		if args.dump: print(e)
		ret = fz.add_edge(str(e))

	if args.dump: sys.exit()

	shortest_path = fz.get_path(source=source, target=target)


	
	#-------------------------------Write Output ----------------------------------------------#
	file_handling.write_output(id, args, shortest_path, my_graph, my_orfs)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#

