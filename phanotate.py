#!/usr/bin/env python
import os
import sys
import getopt

#from subprocess import Popen, PIPE, STDOUT
__requires__= 'fastpath==1.4'
import pkg_resources
pkg_resources.require("fastpath==1.4")
import fastpathz as fz
import phanotate_connect as pc

#from phanotate_modules import functions
#from phanotate_modules.nodes import Node
from phanotate_modules import file_handling
from phanotate_modules.features import Features


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
	contig_features = Features(**vars(args))
	#-------------------------------Find the ORFs----------------------------------------------#

	contig_features.parse_contig(id, dna)

	contig_features.score_orfs()

	# find other features
	from phanotate_modules.trnas import tRNAs
	for trna in tRNAs(contig_features.dna):
		contig_features.add_feature( trna )

	'''
	print(contig_features.pnots)
	for orfs in contig_features.iter_out():
		for orf in orfs:
			print(orf.left(), orf.right(), orf.stop, orf.rbs, orf.score_rbs(), orf.pstop, orf.start_codon(), orf.weight, sep='\t')

	scores = list()
	for orf in contig_features.iter_orfs():
			scores.append(orf.score_rbs())

	from statistics import mean
	m = mean(scores)

	counts = {'ATG':0,'GTG':0,'TTG':0}
	for orf in contig_features.iter_orfs():
		if orf.score_rbs() > m:
			counts[orf.start_codon()] += 1

	gc = (contig_features.dna.count('G') + contig_features.dna.count('C')) / contig_features.contig_length()
	print(id, gc, counts['ATG'], counts['GTG'], counts['TTG'], sep='\t')
	exit()
	'''
	#-------------------------------Create the Graph-------------------------------------------#

	fz.empty_graph()

	# write source and target edges to graph
	fz.add_edge( contig_features.source_edge() )
	fz.add_edge( contig_features.target_edge() )

	# write edges to the graph
	for orf in contig_features:
		ret = fz.add_edge( orf.as_scaled_edge() )


	# write edges to pconnect to get interconnections
	for edge in fz.get_edges():
		ret = pc.add_edge( edge )

	# find regions with no features that would break the path
	for edge in pc.get_gaps():
		ret = fz.add_edge( edge )
		ret = pc.add_edge( edge )

	# write connections to the graph
	for edge in pc.get_connections(pnots=contig_features.pnots):
		ret = fz.add_edge( edge )

	if args.dump:
		for edge in fz.get_edges():
			print(edge)
		sys.exit()
	#-------------------------------Run Bellman-Ford-------------------------------------------#

	shortest_path = fz.get_path(source= contig_features.source_node(), target= contig_features.target_node())[1:-1]

	#-------------------------------Write Output ----------------------------------------------#
	file_handling.write_output(contig_features, shortest_path)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#

