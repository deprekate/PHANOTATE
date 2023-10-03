#!/usr/bin/env python3
import os
import sys
import getopt

#from subprocess import Popen, PIPE, STDOUT
sys.path.pop(0)
from phanotate_modules.file import File
import fastpathz as fz

from phanotate_modules import file_handling
from phanotate_modules import functions
from phanotate_modules.nodes import Node
from phanotate_modules.edges import Edge
from phanotate_modules.file_handling import pairwise

sign = lambda x: (1, -1)[x<0]


#--------------------------------------------------------------------------------------------------#
#                               ARGUMENTS                                                          #
#--------------------------------------------------------------------------------------------------#

args = file_handling.get_args(File)
if args.format == 'fasta':
	args.format = 'fna'
#--------------------------------------------------------------------------------------------------#
#                               FILE INPUT                                                         #
#--------------------------------------------------------------------------------------------------#

base_trans = str.maketrans('sbvdefhijklmnopqruwxyz','gggaaaaaaaaaaaaaaaaaaa')
genbank = File(args.infile);
if not genbank.seq():
	sys.stdout.write("Error: no sequences found in infile\n")
	sys.exit()

#--------------------------------------------------------------------------------------------------#
#                               MAIN ROUTINE                                                       #
#--------------------------------------------------------------------------------------------------#
for locus in genbank:
	#-------------------------------Find the ORFs----------------------------------------------#
	locus.start_codons = args.start_codons
	locus.stop_codons = args.stop_codons
	locus.min_orf_len = args.min_orf_len
	orfs = functions.get_orfs(locus)


	#-------------------------------Create the Graph-------------------------------------------#
	graph = functions.get_graph(orfs)


	#-------------------------------Run Bellman-Ford-------------------------------------------#
	source = "Node('source','source',0,0)"
	target = "Node('target','target',0," + str(locus.length()+1) + ")"
	# Write edges to the fastpath program, and multiply the weight to not lose decimal places
	fz.empty_graph()
	for e in graph.iteredges():
		if args.dump: print(e)
		ret = fz.add_edge(str(e))

	if args.dump: sys.exit()

	shortest_path = fz.get_path(source=source, target=target)


	
	#-------------------------------Write Output ----------------------------------------------#
	#file_handling.write_output(locus.name(), args, shortest_path, graph, orfs)
	shortest_path = shortest_path[1:]
	for source, target in pairwise(shortest_path):
		left,right = eval(source) , eval(target)
		weight = graph.weight(Edge(left,right,0))
		pairs = [[left.position , right.position]]
		feature = locus.add_feature(left.gene, sign(left.frame), pairs, {'note':['score:%E' % weight]})
		feature.weight = '%E' % weight
	locus.write(args)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#

