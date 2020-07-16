from __future__ import division
from collections import deque
from decimal import Decimal
import sys
import itertools


def min_idx(a, b, c):
	if a < b and a < c:
		return 1
	if b < a and b < c:
		return 2
	if c < a and c < b:
		return 3
	return 0

def max_idx(a, b, c):
	if a > b and a > c:
		return 1
	if b > a and b > c:
		return 2
	if c > a and c > b:
		return 3
	return 0


class GCframe(list):
	def __init__(self, seq=None, window = 120):
		if seq and len(seq) < window:
			self.window = len(seq) // 3 // 2 * 2
		else:
			self.window = window//3
		self.states = itertools.cycle([1, 2, 3])
		self.bases = [None]*4
		self.bases[1] = deque(['-']*self.window)
		self.bases[2] = deque(['-']*self.window)
		self.bases[3] = deque(['-']*self.window)
		self.frequency = [None]*4
		self.frequency[1] = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
		self.frequency[2] = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
		self.frequency[3] = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
		self.total = [deque([None]),deque([None]),deque([None]),deque([None])]
		self.min_frame = []
		self.max_frame = []
		if seq:
			for base in seq:
				self.add_base(base)
			self._close()

	def add_base(self, base):
		frame = next(self.states)
		self.bases[frame].append(base)
		self.frequency[frame][base] += 1
		item = self.bases[frame].popleft()
		self.frequency[frame][item] -= 1
		self.total[frame].append(self.frequency[frame]['G'] + self.frequency[frame]['C'])

	def min_frame_at(self, i):
		return self.min_frame[i]

	def max_frame_at(self, i):
		return self.max_frame[i]

	def add(self, a, b, c):
			self.append( [a,b,c] )
			self.min_frame.append( min_idx(a,b,c) )
			self.max_frame.append( max_idx(a,b,c) )
		
	def _close(self):
		for _ in range(self.window//2):
			for frame in [1,2,3]:
				self.total[frame].popleft()
				item = self.bases[frame].popleft()
				self.frequency[frame][item] -= 1
				self.total[frame].append(self.frequency[frame]['G'] + self.frequency[frame]['C'])
		# this is to set the arrays for one-based indexing to match fasta/genbank
		self.add( 0, 0, 0 )
		# calculate the index of the min/max gc content frame
		for i in range(len(self.total[3])-1):
			self.add( self.total[1][i], self.total[2][i],   self.total[3][i]   )
			self.add( self.total[2][i], self.total[3][i],   self.total[1][i+1] )
			self.add( self.total[3][i], self.total[1][i+1], self.total[2][i+1] )
		i += 1
		#self.add(     self.total[1][i], self.total[2][i],   self.total[3][i]   )
		if(i < len(self.total[1])-1):
			self.add( self.total[2][i], self.total[3][i],   self.total[1][i+1] )
		if(i < len(self.total[2])-1):
			self.add( self.total[3][i], self.total[1][i+1], self.total[2][i+1] )


