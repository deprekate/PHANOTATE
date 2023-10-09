#!/usr/bin/env python3

import os
import io
import sys
from itertools import zip_longest, chain, tee, islice, tee

from genbank.locus import Locus
from phanotate_modules.feature import Feature

def previous_and_next(some_iterable):
	prevs, items, nexts = tee(some_iterable, 3)
	prevs = chain([None], prevs)
	nexts = chain(islice(nexts, 1, None), [None])
	return zip(prevs, items, nexts)

def pairwise(iterable):
	"s -> (s0,s1), (s1,s2), (s2, s3), ..."
	a, b = tee(iterable)
	next(b, None)
	return zip(a, b)

class Locus(Locus, feature=Feature):
	def init(self, args):
		#self.start_codons = ['atg','gtg','ttg']
		#self.stop_codons  = ['taa','tga','tag']
		pass

	def add_feature(self, key, strand, pairs, tags=dict()):
		pairs[-1][-1] += 2
		pairs = [map(str,pair) for pair in pairs]
		feature = super(Locus,self).add_feature(key, strand, pairs, tags)
		#if(strand > 0 and not feature.has_start()) or (strand < 0 and not feature.has_stop()):
		#	feature.set_left('<1')
		#if(strand > 0 and not feature.has_stop()) or (strand < 0 and not feature.has_start()):
		#	feature.set_right('>' + str(self.length()))	
		return feature

	def tabular(self, outfile=sys.stdout):
		outfile.write("#id:\t" + self.name() + "\n")
		outfile.write("#START\tSTOP\tFRAME\tCONTIG\tSCORE\n")
		for feature in self.features(include=['CDS']):
			pairs = feature.pairs
			left,right = pairs[0][0] , pairs[-1][-1]
			if feature.strand < 0:
				right,left = left,right
			outfile.write(left)
			outfile.write('\t')
			outfile.write(right)
			outfile.write('\t')
			outfile.write(chr(44 - feature.strand))
			outfile.write('\t')
			outfile.write(self.name())
			outfile.write('\t')
			outfile.write(feature.weight)				
			outfile.write('\n')
