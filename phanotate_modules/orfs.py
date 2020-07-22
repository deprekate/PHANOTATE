import sys
import itertools
from math import log10
from decimal import Decimal

#from .kmeans import KMeans
from .gc_frame_plot import GCFramePlot
from score_rbs import ScoreXlationInit


def rev_comp(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G',
		    'N':'N',
		    'R':'Y','Y':'R','S':'S','W':'W','K':'M','M':'K',
		    'B':'V','V':'B','D':'H','H':'D'}
	return "".join([seq_dict[base] for base in reversed(seq)])

def pstop(seq):
	length = Decimal(len(seq))
	Pa = seq.count('A') / length
	Pt = seq.count('T') / length
	Pg = seq.count('G') / length
	Pc = seq.count('C') / length
	return (Pt*Pa*Pa + Pt*Pg*Pa + Pt*Pa*Pg)
	

class Orfs(dict):
	"""The class holding the orfs"""
	def __init__(self, n=0, **kwargs):
		self.__dict__.update(kwargs)
		self.n = n
		self.dna = None
		self.pstop = None
		self.min_orf_len = 90
		self.other_end = dict()
		self.start_codons = ['ATG', 'GTG', 'TTG']
		self.stop_codons = ['TAA', 'TGA', 'TAG']
		self.start_weight = {'ATG':Decimal('1.00'), 'CAT':Decimal('1.00'),
							 'GTG':Decimal('0.12'), 'CAC':Decimal('0.12'),
							 'TTG':Decimal('0.05'), 'CAA':Decimal('0.05')}
		self.gc_frame_plot = None
		self.rbs_scorer = ScoreXlationInit()

	def score(self):
		pos_max = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]
		pos_min = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]

		self.classify_orfs()
		for orfs in self.iter_in():
			for orf in orfs:
				#if(orf.start_codon() == 'ATG'):
				if(orf.good):
					n = int(orf.length()/10)
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
		print(pos_min)
		print(pos_max)
		exit()
		s = 1/self.hold
		if(self.start_codon() in self.start_weight):
			s = s * self.start_weight[self.start_codon()]
		s = s * Decimal(str(self.weight_rbs))
		self.weight = -s

	def score_rbs(self, seq):
		return 1 + self.rbs_scorer.score_init_rbs(seq, 20)[0]

	def seq(self, a, b):
			return self.dna[ a-1 : b ]

	def add_orf(self, start, stop, frame, seq, rbs):
		""" Adds an orf to the factory"""
		if len(seq) < self.min_orf_len: return
		
		o = Orf(start, stop, frame, self)

		if stop not in self:
			self[stop] = dict()
			self[stop][start] = o
			self.other_end[stop] = start
			self.other_end[start] = stop
		elif start not in self[stop]:
			self[stop][start] = o
			self.other_end[start] = stop
			if(frame > 0 and start < self.other_end[stop]):
				self.other_end[stop] = start
			elif(frame < 0 and start > self.other_end[stop]):
				self.other_end[stop] = start
		else:
			raise ValueError("orf already defined")

	def iter_orfs(self):
		for stop in self:
			for start in self[stop]:
				yield self[stop][start]
	def iter_in(self):
		for stop in self.keys():
			keylist = list(self[stop].keys())
			if(self[stop][keylist[0]].frame > 0):
				keylist.sort()
			else:
				keylist.sort(reverse=True)
				
			yield (self[stop][start] for start in keylist)
	def iter_out(self):
		for stop in self.keys():
			keylist = list(self[stop].keys())
			if(self[stop][keylist[0]].frame > 0):
				keylist.sort(reverse=True)
			else:
				keylist.sort()
			yield (self[stop][start] for start in keylist)
	def first_start(self, stop):
		if stop in self:
			list = sorted(list(self[stop].keys()))
			if(list[0] < stop):
				return list[0]
			else:
				return list[-1]
	def get_orf(self, start, stop):
		if stop in self:
			if start in self[stop]:
				return self[stop][start]
			else:
				raise ValueError("orf with start codon not found")
		else:
			raise ValueError(" orf with stop codon not found")

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

	def parse_contig(self, dna):
		self.dna = dna
		self.pstop = pstop(self.dna + rev_comp(self.dna))

		self.gc_frame_plot = GCFramePlot(dna)
	

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

	def calculate_weights(self):
		pass
		
	
class Orf:
	def __init__(self, start, stop, frame, parent): #, frame, seq=None, rbs=None, start_codons=None, stop_codons=None):
		self.start = start
		self.stop = stop
		self.frame = frame
		self.parent = parent
		self.dna = self.dna()
		self.amino_acids = self.amino_acids()
		self.rbs = self.rbs()
		self.pstop = pstop(self.dna)
		self.pnots = 1 - self.pstop
		self.min_frames, self.max_frames = self.gc_frame_plot()
		#self.weight = self.weight()

		'''
		self.scores = dict()
		self.aa = dict()
		self.med = dict()
		self.good = 0
		'''
		#self.parse_seq()

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
		
	def weight():
		hold = hold * ((self.pnots**pos_max[ind_max])**pos_min[ind_min])

		s = 1/hold
		if(self.start_codon() in self.start_weight):
			s = s * self.start_weight[self.start_codon()]
		s = s * Decimal(str(self.weight_rbs))
		return -s

	def score_rbs(self):
		return self.parent.score_rbs(self.rbs)

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

	def amino_acids(self):
		#calculate the amino acid frequency
		nucs = ['T', 'C', 'A', 'G']
		codons = [a+b+c for a in nucs for b in nucs for c in nucs]
		amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		codon_table = dict(zip(codons, amino_acids))
		aa = []
		for i in range(0, self.length(), 3):
			aa.append(codon_table[self.dna[i:i+3]])
		return "".join(aa)

	def amino_acid_count(self, aa):
		return self.amino_acids.count(aa) 
			
	def amino_acid_frequency(self, aa):
		return self.amino_acid_count(aa) / (len(self.amino_acids) - self.amino_acids.count('*'))

	def something(self):
		frequencies = dict()
		for aa in self.amino_acids:
			frequencies[aa] = frequencies.get(aa, 0) + 1
		return frequencies
		

	def parse_seq(self):
		amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		self.aa = dict()
		for a in amino_acids:
			self.aa[a] = 0
		for i in range(0, self.length(), 3):
			self.aa[codon_table[self.seq[i:i+3]]] += 1
			self.aa['*'] += 1
		#calculate the MED scores
		H = 0
		for aa in list('ACDEFGHIKLMNPQRSTVWY'):
			p = self.aa[aa]/self.aa['*']
			#p = self.aa[aa]
			if(p):
				H += p * log10(p)
		for aa in list('ACDEFGHIKLMNPQRSTVWY'):
			p = self.aa[aa]/self.aa['*']
			if(p):
				self.med[aa] = (1/H) * p * log10(p)
			else:
				self.med[aa] = 0

	def score(self):
		s = 1/self.hold
		if(self.start_codon() in self.start_weight):
			s = s * self.start_weight[self.start_codon()]
		s = s * Decimal(str(self.weight_rbs))
		self.weight = -s
		
	def start_codon(self):
		return self.dna[0:3]

	def stop_codon(self):
		return self.dna[-3:]

	def has_start(self):
		return self.start_codon() in self.start_codons
		
	def has_stop(self):
		return self.stop_codon() in self.stop_codons

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

	def pstop(self):
		frequency = {'A':0, 'T':0, 'C':0, 'G':0}
		for base in self.dna:
			if(base not in ['A', 'C', 'T', 'G']):
				continue
			frequency[base] += 1
		length = Decimal(len(self.dna))
		Pa = frequency['A']/length
		Pt = frequency['T']/length
		Pg = frequency['G']/length
		Pc = frequency['C']/length
		return (Pt*Pa*Pa + Pt*Pg*Pa + Pt*Pa*Pg)

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

