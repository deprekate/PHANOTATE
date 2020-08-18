import sys
import itertools
from math import log10
from decimal import Decimal

#from .kmeans import KMeans
from .gc_frame_plot import GCFramePlot
from .orf import CDS
from .functions import *
from score_rbs import ScoreXlationInit


class Features(list):
	"""The class holding the orfs"""
	def __init__(self, n=0, **kwargs):
		self.__dict__.update(kwargs)

		self.n = n
		self.dna = None
		self.pstop = None
		self.gc_frame_plot = None
		self.cds = dict()
		self.feature_at = dict()

	def score_orfs(self):
		pos_max = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]
		pos_min = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]

		self.pos_min = pos_min
		self.pos_max = pos_max
		#self.classify_orfs()
		for orfs in self.iter_orfs('in'):
			for orf in orfs:
				if(orf.start_codon() == 'ATG'):
				#if(orf.good):
					#n = int(orf.length()/10)
					#for base in range(start+n, stop-36, 3):
					for min_frame,max_frame in zip(orf.min_frames[10:-10], orf.max_frames[10:-10]):
							pos_min[min_frame] += 1
							pos_max[max_frame] += 1
					break
		# normalize to one
		y = max(pos_min)
		pos_min[:] = [x / y for x in pos_min]	
		y = max(pos_max)
		pos_max[:] = [x / y for x in pos_max]	
		
		for orf in self.iter_orfs():
			orf.weight = 1
			for min_frame,max_frame in zip(orf.min_frames, orf.max_frames):
				orf.weight = orf.weight * ((orf.pnots**pos_min[min_frame])**pos_max[max_frame])
			orf.weight = 1 / orf.weight
			if(orf.start_codon() in self.start_weight):
				orf.weight = orf.weight * Decimal(self.start_weight[orf.start_codon()])
			orf.weight = orf.weight * self.score_rbs(orf.rbs)
			orf.weight = -orf.weight



	def score_rbs(self, seq):
		return 1 + Decimal(str(self.rbs_scorer.score_init_rbs(seq, 20)[0]))

	def seq(self, a, b):
			return self.dna[ a-1 : b ]

	def add_orf(self, start, stop, frame, seq, rbs):
		""" Adds an orf to the factory"""
		if len(seq) < self.min_orf_len: return

		o = CDS(start, stop, frame, self)

		orfs = self.cds
		if stop not in orfs:
			orfs[stop] = dict()
			orfs[stop][start] = o
		elif start not in orfs[stop]:
			orfs[stop][start] = o
		else:
			raise ValueError("orf already defined")
		self.add_feature(o)

	def add_feature(self, feature):
		self.feature_at[ feature.as_scaled_edge()[:2] ] = feature
		self.append(feature)


	def iter_features(self, type_ = None):
		for feature in self:
			if feature.type != type_:
				continue
			yield feature

	def iter_orfs(self, kind=None):
		orfs = self.cds
		if not kind:
			for stop in orfs:
				for start in orfs[stop]:
					yield orfs[stop][start]
		elif kind == 'in':
			for stop in orfs.keys():
				keylist = list(orfs[stop].keys())
				if(orfs[stop][keylist[0]].frame > 0):
					keylist.sort()
				else:
					keylist.sort(reverse=True)
				yield (orfs[stop][start] for start in keylist)
		elif kind == 'out':
			for stop in orfs.keys():
				keylist = list(orfs[stop].keys())
				if(orfs[stop][keylist[0]].frame > 0):
					keylist.sort(reverse=True)
				else:
					keylist.sort()
				yield (orfs[stop][start] for start in keylist)

	def get_orf(self, start, stop):
		orfs = self.cds
		if stop in orfs:
			if start in orfs[stop]:
				return orfs[stop][start]
			else:
				raise ValueError("orf with start codon not found")
		else:
			raise ValueError(" orf with stop codon not found")

	def get_feature(self, left, right):
		return self.feature_at[ (left, right) ]

	def contig_length(self):
		return len(self.dna)
	
	def end(self, frame):
		return self.contig_length() - ((self.contig_length() - (frame-1))%3)
	
	def classify_orfs(self):
		import numpy as np
		from sklearn.preprocessing import StandardScaler
		from sklearn.cluster import KMeans
		from .kmeans import KMeans as KM
		X = []
		for orf in self.iter_orfs():
			point = []
			for aa in list('ARNDCEQGHILKMFPSTWYV'):
				point.append(orf.amino_acid_frequency(aa))
			X.append(point)
			
		kmeans = KM(n_clusters=3).fit(X)
		val, idx = min((val, idx) for (idx, val) in enumerate(kmeans.withinss_))

		for i, orf in enumerate(self.iter_orfs()):
			if kmeans.labels_[i] == idx:
				orf.good = 1
			else:
				orf.good = 0

	def parse_contig(self, id, dna):
		self.id = id
		self.dna = dna
		self.pstop = pstop(self.dna + rev_comp(self.dna))
		self.pnots = 1 - self.pstop

		self.gc_frame_plot = GCFramePlot(dna)
		self.rbs_scorer = ScoreXlationInit()
	

		#background_rbs = [1.0] * 28
		#training_rbs = [1.0] * 28
		# find nucleotide frequency
	
		#y = sum(background_rbs)
		#background_rbs[:] = [x/y for x in background_rbs]
	
		# The dicts that will hold the start and stop codons
		stops = {1:0, 2:0, 3:0, -1:1, -2:2, -3:3}
		starts = {1:[], 2:[], 3:[], -1:[], -2:[], -3:[]}
		if dna[0:3] not in self.start_codons:
			starts[1].append(1)
		if dna[1:4] not in self.start_codons:
			starts[2].append(2)
		if dna[2:5] not in self.start_codons:
			starts[3].append(3)
	
		# Reset iterator and find all the open reading frames
		states = itertools.cycle([1, 2, 3])
		for i in range(1, (len(dna)-1)):
			codon = dna[i-1:i+2]
			frame = next(states)
			if codon in self.start_codons:
				starts[frame].append(i)
			elif rev_comp(codon) in self.start_codons:
				starts[-frame].append(i+2)
			elif codon in self.stop_codons:
				stop = i+2
				for start in reversed(starts[frame]):
					seq = dna[start-1:stop]
					rbs = dna[start-21:start]
					self.add_orf(start, stop-2, frame, seq, rbs)
		
				starts[frame] = []
				stops[frame] = stop
			elif rev_comp(codon) in self.stop_codons:
				stop = stops[-frame]
				for start in starts[-frame]:
					seq = rev_comp(dna[max(0,stop-1):start])
					rbs = rev_comp(dna[start:start+21])
					self.add_orf(start-2, stop, -frame, seq, rbs)
		
				starts[-frame] = []
				stops[-frame] = i
		# Add in any fragment ORFs at the end of the genome
		for frame in [1, 2, 3]:
			for start in reversed(starts[frame]):
				stop = self.end(frame)
				seq = dna[max(0,start-1):stop]
				rbs = dna[start-21:start]
				self.add_orf(start, stop-2, frame, seq, rbs)
			start = self.end(frame)
			if rev_comp(dna[start-3:start]) not in self.start_codons:
				starts[-frame].append(self.end(frame))	
			for start in starts[-frame]:
				stop = stops[-frame]
				seq = rev_comp(dna[max(0,stop-1):start])
				rbs = rev_comp(dna[start:start+21])
				self.add_orf(start-2, stop, -frame, seq, rbs)

	def source_node(self):
		return '0_source'

	def source_edge(self):
		return (self.source_node(), self.source_node(), '1')

	def target_node(self):
		return '%s_target' % str(self.contig_length())

	def target_edge(self):
		return (self.target_node(), self.target_node(), '1')
	
	def first_start(self, stop):
		if stop in self:
			list = sorted(list(self[stop].keys()))
			if(list[0] < stop):
				return list[0]
			else:
				return list[-1]


