#!/usr/bin/env python
import os
import sys
import getopt
import pkg_resources

#from subprocess import Popen, PIPE, STDOUT
pkg_resources.require("fastpath>=1.4")
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


	#-------------------------------Create the Graph-------------------------------------------#

	fz.empty_graph()

	# write source and target edges to graph
	fz.add_edge( contig_orfs.source_edge() )
	fz.add_edge( contig_orfs.target_edge() )

	# write edges to the graph
	for orf in contig_orfs.iter_orfs():
		ret = fz.add_edge( orf.as_scaled_edge() )


	
	# write edges to pconnect to get interconnections
	for edge in fz.get_edges():
		ret = pc.add_edge( edge )

	# find regions with no features that would break the path
	for edge in pc.get_gaps():
		ret = pc.add_edge( edge )
	exit()

	# write edges to the graph
	for edge in pc.get_connections(pnots=contig_orfs.pnots):
		ret = fz.add_edge( edge )

	if args.dump:
		for edge in fz.get_edges():
			print(edge)
		sys.exit()
	#-------------------------------Run Bellman-Ford-------------------------------------------#

	shortest_path = fz.get_path(source= contig_orfs.source_node(), target= contig_orfs.target_node())[1:-1]
	
	#-------------------------------Write Output ----------------------------------------------#
	file_handling.write_output(contig_orfs, shortest_path)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#

