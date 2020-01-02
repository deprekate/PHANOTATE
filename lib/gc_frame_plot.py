from __future__ import division
from collections import deque
from decimal import Decimal
import sys
import itertools

def max_idx(a, b, c):
	if(a > b):
		if(a > c):
			return 1;
		else:
			return 3;
	else:
		if(b > c):
			return 2;
		else:
			return 3;
def min_idx(a, b, c):
	if(a > b):
		if(b > c):
			return 3;
		else:
			return 2;
	else:
		if(a > c):
			return 3;
		else:
			return 1;
class GCframe:
	def __init__(self, window = 120):
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
		self.total = [deque([]),deque([]),deque([]),deque([])]
		self.freq =[]

	def add_base(self, base):
		frame = next(self.states)
		self.bases[frame].append(base)
		self.frequency[frame][base] += 1
		item = self.bases[frame].popleft()
		self.frequency[frame][item] -= 1
		self.total[frame].append(self.frequency[frame]['G'] + self.frequency[frame]['C'])


	def _close(self):
		for _ in range(self.window//2):
			for frame in [1,2,3]:
				self.total[frame].popleft()
				item = self.bases[frame].popleft()
				self.frequency[frame][item] -= 1
				self.total[frame].append(self.frequency[frame]['G'] + self.frequency[frame]['C'])

	def get(self):
		self._close()
		self.freq.append( [20,20,20] )
		for i in range(len(self.total[3])-1):
			self.freq.append( [self.total[1][i],self.total[2][i],self.total[3][i]] )
			self.freq.append( [self.total[2][i],self.total[3][i],self.total[1][i+1]] )
			self.freq.append( [self.total[3][i],self.total[1][i+1],self.total[2][i+1]] )
		i += 1
		self.freq.append( [self.total[1][i],self.total[2][i],self.total[3][i]] )
		if(i < len(self.total[1])-1):
			self.freq.append( [self.total[2][i],self.total[3][i],self.total[1][i+1]] )
		if(i < len(self.total[2])-1):
			self.freq.append( [self.total[3][i],self.total[1][i+1],self.total[2][i+1]] )
		return self.freq

