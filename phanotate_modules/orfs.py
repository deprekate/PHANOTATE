import sys
import itertools
from math import log10
from decimal import Decimal

from .kmeans import KMeans
from .gc_frame_plot import GCFramePlot


def rev_comp(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G',
		    'N':'N',
		    'R':'Y','Y':'R','S':'S','W':'W','K':'M','M':'K',
		    'B':'V','V':'B','D':'H','H':'D'}
	return "".join([seq_dict[base] for base in reversed(seq)])

class Orfs(dict):
	"""The class holding the orfs"""
	def __init__(self, n=0, **kwargs):
		self.__dict__.update(kwargs)
		self.n = n
		self.seq = None
		self.pstop = None
		self.min_orf_len = 90
		self.other_end = dict()
		self.start_codons = ['ATG', 'GTG', 'TTG']
		self.stop_codons = ['TAA', 'TGA', 'TAG']
		self.gcframeplot_consensus_min = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]
		self.gcframeplot_consensus_max = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]

	def add_orf(self, start, stop, frame, seq, rbs):
		if len(seq) < self.min_orf_len: return

		o = Orf(start, stop, frame, seq, rbs, self.start_codons, self.stop_codons)
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
		return len(self.seq)
	
	def end(self, frame):
		return self.contig_length() - ((self.contig_length() - (frame-1))%3)
	
	def classify_orfs(self):
		X = []
		for orf in self.iter_orfs():
			point = []
			for aa in list('ARNDCEQGHILKMFPSTWYV'):
				point.append(orf.med[aa])
			X.append(point)

		kmeans = KMeans(n_clusters=2).fit(X)
		val, idx = min((val, idx) for (idx, val) in enumerate(kmeans.withinss_))

		for i, orf in enumerate(self.iter_orfs()):
			if kmeans.labels_[i] == idx:
				orf.good = 1
		#print kmeans.cluster_centers_

	def parse_contig(self, dna):
		self.seq = dna
	
		# find nucleotide frequency, kmers, and create gc frame plot
		frequency = {'A':Decimal(0), 'T':Decimal(0), 'C':Decimal(0), 'G':Decimal(0)}
		frame_plot = GCframe()
		#background_rbs = [1.0] * 28
		#training_rbs = [1.0] * 28

		for i, base in enumerate(dna):
			# nucleotide frequency
			frequency[base] += 1
			frequency[rev_comp(base)] += 1
			#kmers for rbs
			#background_rbs[score_rbs(dna[i:i+21])] += 1
			#background_rbs[score_rbs(rev_comp(dna[i:i+21]))] += 1
			#gc frame plot
			#frame_plot.add_base(base)
		#gc_pos_freq = frame_plot.get()
		print(dna[1:10])
		if 1: return 
	
		Pa = frequency['A'] / (2*self.contig_length())
		Pt = frequency['T'] / (2*self.contig_length())
		Pg = frequency['G'] / (2*self.contig_length())
		Pc = frequency['C'] / (2*self.contig_length())
		self.pstop = (Pt*Pa*Pa + Pt*Pg*Pa + Pt*Pa*Pg)
	
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
	def __init__(self, start, stop, frame, seq, rbs, start_codons, stop_codons):
		self.start = start
		self.stop = stop
		self.frame = frame
		self.seq = seq
		self.rbs = rbs
		self.rbs_score = None
		self.pstop = self.p_stop()
		self.weight = 1
		self.weight_start = 1
		self.weight_rbs = 1
		self.hold = 1
		self.gcfp_mins = 1
		self.gcfp_maxs = 1
		self.start_codons = start_codons
		self.stop_codons = stop_codons
		self.start_weight = {'ATG':Decimal('1.00'), 'CAT':Decimal('1.00'),
				     'GTG':Decimal('0.12'), 'CAC':Decimal('0.12'),
				     'TTG':Decimal('0.05'), 'CAA':Decimal('0.05')}
		self.aa = dict()
		self.med = dict()
		self.good = 0

		self.parse_seq()

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

	def parse_seq(self):
		#calculate the amino acid frequency
		nucs = ['T', 'C', 'A', 'G']
		codons = [a+b+c for a in nucs for b in nucs for c in nucs]
		amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		codon_table = dict(zip(codons, amino_acids))
		for a in amino_acids:
			self.aa[a] = 0
		for i in range(0, len(self.seq)-5, 3):
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
		return self.seq[0:3]

	def stop_codon(self):
		return self.seq[-3:]

	def has_start(self):
		return self.start_codon() in self.start_codons
		
	def has_stop(self):
		return self.stop_codon() in self.stop_codons

	def pp_stop(self):
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

	def p_stop(self):
		frequency = {'A':0, 'T':0, 'C':0, 'G':0}
		for base in self.seq:
			if(base not in ['A', 'C', 'T', 'G']):
				continue
			frequency[base] += 1
		length = Decimal(len(self.seq))
		Pa = frequency['A']/length
		Pt = frequency['T']/length
		Pg = frequency['G']/length
		Pc = frequency['C']/length
		return (Pt*Pa*Pa + Pt*Pg*Pa + Pt*Pa*Pg)

	def __repr__(self):
		"""Compute the string representation of the orf"""
		return "%s(%s,%s,%s,%s,%s)" % (
			self.__class__.__name__,
			repr(self.start),
			repr(self.stop),
			repr(self.frame),
			repr(self.weight_rbs),
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

