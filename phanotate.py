#!/usr/bin/env python3
import os
import sys
import getopt

sys.path.pop(0)

import fastpathz as fz
import phanotate.connect as pc

from phanotate import file_handling
from phanotate.features import Features

from phanotate.trnas import tRNAs

#--------------------------------------------------------------------------------------------------#
#                               ARGUMENTS                                                          #
#--------------------------------------------------------------------------------------------------#

args = file_handling.get_args()

#--------------------------------------------------------------------------------------------------#
#                               FILE INPUT                                                         #
#--------------------------------------------------------------------------------------------------#

import time
start = time.time()
base_trans = str.maketrans('SBVDEFHIJKLMNOPQRUWXYZ','GGGAAAAAAAAAAAAAAAAAAA')
contigs = file_handling.read_fasta(args.infile, base_trans)
if not contigs:
	sys.stdout.write("Error: no sequences found in infile\n")
	sys.exit()

print("Input", time.time() - start)

#--------------------------------------------------------------------------------------------------#
#                               MAIN ROUTINE                                                       #
#--------------------------------------------------------------------------------------------------#
for id, dna in contigs.items():
	start = time.time()
	contig_features = Features(**vars(args))
	#-------------------------------Find the ORFs----------------------------------------------#

	contig_features.parse_contig(id, dna)
	print("Input", time.time() - start)

	start = time.time()
	contig_features.score_orfs()
	print("Score", time.time() - start)


	# find other features
	for trna in tRNAs(contig_features.dna):
		contig_features.add_feature( trna )

	#-------------------------------Create the Graph-------------------------------------------#

	fz.empty_graph()
	fz.scaling = 3
	# write source and target edges to graph
	fz.add_edge( contig_features.source_edge() )
	fz.add_edge( contig_features.target_edge() )

	# write edges to the graph
	for orf in contig_features:
		ret = fz.add_edge( orf.as_edge() )

	start = time.time()
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
	print("Connect", time.time() - start)

	if args.dump:
		for edge in fz.get_edges():
			print(edge)
		sys.exit()
	elif args.orfs:
		for orf in contig_features.iter_orfs():
			print(orf.left(), orf.right(), orf.stop, orf.length(), orf.rbs, orf.rbs_score, orf.pstop, orf.start_codon(), orf.weight, orf.good, sep='\t')
		sys.exit()
	#-------------------------------Run Bellman-Ford-------------------------------------------#

	start = time.time()
	shortest_path = fz.get_path(source= contig_features.source_node(), target= contig_features.target_node())[1:-1]
	print("Bellman", time.time() - start)

	#-------------------------------Write Output ----------------------------------------------#
	start = time.time()
	file_handling.write_output(contig_features, shortest_path)
	print("Output", time.time() - start)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#

