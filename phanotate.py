#!/usr/bin/env python
import os
import sys
import getopt

#from subprocess import Popen, PIPE, STDOUT
import fastpathz as fz
import phanotate_connect as pc

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
contigs = file_handling.read_fasta(args.infile, base_trans)
if not contigs:
	sys.stdout.write("Error: no sequences found in infile\n")
	sys.exit()

#--------------------------------------------------------------------------------------------------#
#                               MAIN ROUTINE                                                       #
#--------------------------------------------------------------------------------------------------#
for id, dna in contigs.items():
	contig_orfs = Orfs(**vars(args))
	#-------------------------------Find the ORFs----------------------------------------------#

	contig_orfs.parse_contig(id, dna)

	contig_orfs.score()

	'''
	for orfs in contig_orfs.iter_out():
		for orf in orfs:
			print(orf.left(), orf.right(), orf.stop, orf.rbs, orf.score_rbs(), orf.pstop, orf.start_codon(), orf.weight, sep='\t')
	exit()
	'''
	#-------------------------------Create the Graph-------------------------------------------#
#	my_graph = functions.get_graph(contig_orfs)
	
	fz.empty_graph()
	for orf in contig_orfs.iter_orfs():
		# write edges to the graph
		ret = fz.add_edge( orf.as_edge() )
		# write edges to pconnect to get interconnections
		ret = pc.add_edge(orf.left_node(), orf.right_node(), orf.frame, orf.pnots)

	for edge in pc.get_connections():
		ret = fz.add_edge( "\t".join(map(str,edge)) )
	
	for orf in contig_orfs.iter_orfs():
		if(orf.left() <= 2000):
			edge = "\t".join(map(str, ['source', orf.left_node(), 1/(contig_orfs.pnots**orf.left())] ))
			ret = fz.add_edge( edge )
		if(contig_orfs.contig_length() - orf.right() <= 2000):
			l = contig_orfs.contig_length() - orf.right()
			edge = "\t".join(map(str, [orf.right_node(), 'target', 1/(contig_orfs.pnots**l)] ))
			ret = fz.add_edge( edge )


	#-------------------------------Run Bellman-Ford-------------------------------------------#

	shortest_path = fz.get_path(source='source', target='target')[1:-1]
	
	#-------------------------------Write Output ----------------------------------------------#
	file_handling.write_output(contig_orfs, shortest_path)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#

