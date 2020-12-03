import sys
from math import log
from decimal import Decimal

#from .kmeans import KMeans
from .functions import *


class CDS:
	def __init__(self, start, stop, frame, parent): #, frame, seq=None, rbs=None, start_codons=None, stop_codons=None):
		self.type = 'CDS'
		self.start = start
		self.stop = stop
		self.frame = frame
		self.weight = 1
		self.parent = parent
		self.dna = self.dna()
		#self.amino_acids = self.amino_acids()
		self.rbs = self.rbs()
		self.pstop = pstop(self.dna)
		self.pnots = 1 - self.pstop
		self.min_frames, self.max_frames = self.gc_frame_plot()
		#self.weight = self.weight()
		self.good = False

	def as_scaled_edge(self):
		if self.frame > 0:
			return ("%s_start" % self.left(), "%s_stop" % (self.right()-2), "%s" % (self.weight*1000))
		else:
			return ("%s_stop" % self.left(), "%s_start" % (self.right()-2), "%s" % (self.weight*1000))

	def gc_frame_plot(self):
		gcfp = self.parent.gc_frame_plot
		rev = lambda x : 0 if not x else 3 - x + 1
		min_frames = []
		max_frames = []
		if(self.frame > 0):
			for base in range(self.start, self.stop+1, 3):
				min_frames.append(gcfp.min_frame_at(base))
				max_frames.append(gcfp.max_frame_at(base))
		else:
			for base in range(self.start, self.stop-1, -3):
				min_frames.append(rev(gcfp.min_frame_at(base)))
				max_frames.append(rev(gcfp.max_frame_at(base)))
		return min_frames, max_frames	

	'''
	def score(self):
		weight = 1
		pos_min = self.parent.pos_min
		pos_max = self.parent.pos_max
		for min_frame,max_frame in zip(self.min_frames, self.max_frames):
			weight = weight * ((self.pnots**pos_min[min_frame])**pos_max[max_frame])

		weight = 1 / weight

		if(self.start_codon() in self.parent.start_weight):
			weight = weight * Decimal(self.start_weight[self.start_codon()])

		weight = weight * self.parent.score_rbs(orf.rbs)

		return -weight
	'''	
	'''
	def score_rbs(self):
		return 1 + self.parent.score_rbs(self.rbs)
	'''

	def dna(self):
		if self.frame > 0:
			return          self.parent.seq(self.left(), self.right())
		else:
			return rev_comp(self.parent.seq(self.left(), self.right()))

	def length(self):
		return len(self.dna)

	def rbs(self):
		if self.frame > 0:
			return          self.parent.seq(self.left()-19, self.left()-1)
		else:
			return rev_comp(self.parent.seq(self.right()+1, self.right()+19))

	def direction(self):
		if self.frame > 0:
			return '+'
		else:
			return '-'

	def begin(self):
		if self.frame > 0:
			c = '' if self.start_codon() in self.parent.start_codons else '<'
			return c + str(self.start)
		else:
			c = '' if self.start_codon() in self.parent.start_codons else '>'
			return c + str(self.start+2)

	def end(self):
		if self.frame > 0:
			c = '' if self.stop_codon() in self.parent.stop_codons else '>'
			return c + str(self.stop+2)
		else:
			c = '' if self.stop_codon() in self.parent.stop_codons else '<'
			return c + str(self.stop)

	def left(self):
		if self.frame > 0:
			return self.start
		else:
			return self.stop

	def right(self):
		if self.frame > 0:
			return self.stop + 2
		else:
			return self.start + 2

	def left_node(self):
		return str(self.left()) + 'o'

	def right_node(self):
		return str(self.right()-2) + 'o'

	def codon_counts(self):
		codons = dict()
		for i in range(0, self.length(), 3):
			codon = self.dna[i:i+3]
			codons[codon] = codons.get(codon, 0) + 1
		return codons
		
	def amino_acids(self):
		"""Calculate the amino acid frequency"""
		aa = []
		for i in range(0, self.length(), 3):
			aa.append(self.parent.translate_codon[self.dna[i:i+3]])
		return "".join(aa)

	def amino_acid_frequencies(self):
		amino_acids = self.amino_acids()
		counts = dict()
		total = 0
		for aa in list('ACDEFGHIKLMNPQRSTVWY#+*'):
			counts[aa] = amino_acids.count(aa)
			total += amino_acids.count(aa)
		for aa in list('ACDEFGHIKLMNPQRSTVWY#+*'):
			counts[aa] = counts[aa] / total
		return counts


	def amino_acid_count(self, aa):
		return self.amino_acids().count(aa) 
			
	def amino_acid_frequency(self, aa):
		return self.amino_acid_count(aa) / len(self.amino_acids)

	def amino_acid_entropies(self):
		"""Calculate entropy"""
		counts = self.amino_acid_frequencies()
		H = 0
		for aa in list('ACDEFGHIKLMNPQRSTVW#+*Y'):
			p = -counts[aa]*log(counts[aa]) if counts[aa] else 0
			counts[aa] = p
			H += p
		for aa in list('ACDEFGHIKLMNPQRSTVWY#+*'):
			counts[aa] /= H
		return counts

	def start_codon(self):
		return self.dna[0:3]

	def stop_codon(self):
		return self.dna[-3:]

	def has_start(self):
		return self.start_codon() in self.parent.start_codons
		
	def has_stop(self):
		return self.stop_codon() in self.parent.stop_codons

	def pp_stop(self):
		""" this one does codon base specific pstop calculation"""
		frequency = [None]*4
		frequency[1] = {'A':0, 'T':0, 'C':0, 'G':0}
		frequency[2] = {'A':0, 'T':0, 'C':0, 'G':0}
		frequency[3] = {'A':0, 'T':0, 'C':0, 'G':0}
		count = [Decimal(0)]*4
		states = itertools.cycle([1, 2, 3])
		frame = states.next()
		for _ in range(0,3-abs(self.stop-self.start)%3):
			frame = states.next()
		for base in self.seq:
			if(base not in ['A', 'C', 'T', 'G']):
				continue
			frequency[frame][base] += 1
			count[frame] +=1
			frame = states.next()
		Ptaa = (frequency[1]['T']/count[1]) * (frequency[2]['A']/count[2]) * (frequency[3]['A']/count[3])
		Ptga = (frequency[1]['T']/count[1]) * (frequency[2]['G']/count[2]) * (frequency[3]['A']/count[3])
		Ptag = (frequency[1]['T']/count[1]) * (frequency[2]['A']/count[2]) * (frequency[3]['G']/count[3])
		return Ptaa+Ptga+Ptag

	def __repr__(self):
		"""Compute the string representation of the orf"""
		return "%s(%s,%s,%s,%s)" % (
			self.__class__.__name__,
			repr(self.start),
			repr(self.stop),
			repr(self.frame),
			repr(self.weight))

	def __eq__(self, other):
		"""Override the default Equals behavior"""
		if isinstance(other, self.__class__):
			return self.__dict__ == other.__dict__
		return NotImplemented

	def __ne__(self, other):
		"""Define a non-equality test"""
		if isinstance(other, self.__class__):
			return not self.__eq__(other)
		return NotImplemented

	def __hash__(self):
		"""Override the default hash behavior (that returns the id or the object)"""
		return hash(tuple(sorted(self.__dict__.items())))

