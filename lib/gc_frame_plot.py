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
def mid_idx(a, b, c):
	if(a > b):
		if(b > c):
			return 2;
		elif(a > c):
			return 3;
		else:
			return 1;
	else:
		if(a > c):
			return 1;
		elif(b > c):
			return 3;
		else:
			return 2;
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
		self.window = window/3
		self.states = itertools.cycle([1, 2, 3])
		self.bases = [None]*4
		#self.bases[0] = deque(['-']*self.window)
		self.bases[1] = deque(['-']*self.window)
		self.bases[2] = deque(['-']*self.window)
		self.bases[3] = deque(['-']*self.window)
		self.frequency = [None]*4
		#self.frequency[0] = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
		self.frequency[1] = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
		self.frequency[2] = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
		self.frequency[3] = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
		self.total = [deque([]),deque([]),deque([]),deque([])]
		self.freq =[]

	def add_base(self, base):
		frame = self.states.next()
		##this is for general gc skew
		#self.bases[0].append(base)
		#self.frequency[0][base] += 1
		#item = self.bases[0].popleft()
		#self.frequency[0][item] -= 1
		#self.total[0].append(self.frequency[frame]['G'] + self.frequency[frame]['C'])
		##these are for gc frame plots
		self.bases[frame].append(base)
		self.frequency[frame][base] += 1
		item = self.bases[frame].popleft()
		self.frequency[frame][item] -= 1
		self.total[frame].append(self.frequency[frame]['G'] + self.frequency[frame]['C'])

	def middle(a, b, c):
		return min(max(a,b),max(b,c),max(a,c))

	def _close(self):
		# get rid of the first half of the window
		for _ in range(self.window/2):
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

