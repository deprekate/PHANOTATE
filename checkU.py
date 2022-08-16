#!/usr/bin/env python3
import os
import sys
import getopt

sys.path.pop(0)

from genbank.locus import Locus
import fastpathz as fz
import phanotate.connect as pc
sys.path.pop(0)

from phanotate import file_handling
from phanotate.features import Features

from phanotate.trnas import tRNAs


def pairwise(iterable):
	a = iter(iterable)
	return zip(a, a)
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

baseline = None
for notstop,stops in zip([None, 'TGA','TAG','TAA'],['TAA,TAG,TGA', 'TAA,TAG', 'TAA,TGA', 'TAG,TGA']):
	args.stop_codons = stops
	for id, dna in contigs.items():
		locus = Locus(id, dna)
		contig_features = Features(**vars(args))
		#-------------------------------Find the ORFs----------------------------------------------#
	
		contig_features.parse_contig(id, dna)
	
		contig_features.score_orfs()
	
		# find other features
		for trna in tRNAs(contig_features.dna):
			contig_features.add_feature( trna )
	
		#-------------------------------Create the Graph-------------------------------------------#

		pc.empty()

		fz.empty_graph()
		fz.scaling = 3
		# write source and target edges to graph
		fz.add_edge( contig_features.source_edge() )
		fz.add_edge( contig_features.target_edge() )
	
		# write edges to the graph
		for orf in contig_features:
			ret = fz.add_edge( orf.as_edge() )
	
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

		#-------------------------------Run Bellman-Ford-------------------------------------------#

		shortest_path = fz.get_path(source= contig_features.source_node(), target= contig_features.target_node())[1:-1]
		for left, right in pairwise(shortest_path):
			feature = contig_features.get_feature(left, right)
			if feature:
				pairs = tuple(map(str,[feature.left(), feature.right()]))
				locus.add_feature('CDS', feature.strand(), (pairs,) ) 

		if not baseline:
			baseline = locus.gene_coverage()
		elif locus.gene_coverage() / baseline > 1.15:
			print('Stop codon readthrough suspected for', notstop, locus.gene_coverage(), 'vs', baseline)


#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#

