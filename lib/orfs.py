from decimal import Decimal
import sys
import itertools

class Orfs(dict):
	"""The class holding the orfs"""
	def __init__(self, n=0):
		self.n = n
		self.pstop = 0
		self.min_orf_len = 90
		self.contig_length = 0
		self.seq = ''
		self.other_end = dict()
		self.p = dict()

	def add_orf(self, start, stop, length, frame, seq, rbs, rbs_score):
		o = Orf(start, stop, length, frame, seq, rbs, rbs_score)
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
			keylist = self[stop].keys()
			if(self[stop][keylist[0]].frame > 0):
				keylist.sort()
			else:
				keylist.sort(reverse=True)
				
			yield (self[stop][start] for start in keylist)
	def iter_out(self):
		for stop in self.keys():
			keylist = self[stop].keys()
			if(self[stop][keylist[0]].frame > 0):
				keylist.sort(reverse=True)
			else:
				keylist.sort()
			yield (self[stop][start] for start in keylist)
	def first_start(self, stop):
		if stop in self:
			list = sorted(self[stop].keys())
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
			raise ValueError("orf with stop codon not found")

class Orf:
	def __init__(self, start, stop, length, frame, seq, rbs, rbs_score):
		self.start = start
		self.stop = stop
		self.frame = frame
		self.seq = seq
		self.rbs = rbs
		self.rbs_score = rbs_score
		self.length = length
		self.pstop = self.p_stop()
		self.weight = 1
		self.weight_start = 1
		self.weight_rbs = 1
		self.hold = 1
		self.gcfp_mins = 1
		self.gcfp_maxs = 1
		self.start_weight = {'ATG':Decimal('1.00'), 'CAT':Decimal('1.00'),
				     'GTG':Decimal('0.42'), 'CAC':Decimal('0.42'),
				     'TTG':Decimal('0.05'), 'CAA':Decimal('0.05')}
	def score(self):
        	s = 1/self.hold
        	if(self.start_codon() in self.start_weight):
                	s = s * self.start_weight[self.start_codon()]
		s = s * Decimal(str(self.weight_rbs))
        	self.weight = -s
		
	def start_codon(self):
		return self.seq[0:3]
	def end(self):
		if self.frame > 0:
			return self.stop+2
		else:
			return self.stop
	def beg(self):
		if self.frame < 0:
			return self.start+2
		else:
			return self.start

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

