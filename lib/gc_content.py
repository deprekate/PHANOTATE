from __future__ import division
from collections import deque
from decimal import Decimal
import sys
import itertools

class GCcontent:
	def __init__(self, window = 120):
		self.window = window
		self.bases = deque(['-']*(self.window))
		self.frequency = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
		self.total = deque([self.frequency])

	def add_base(self, base):
		#this is for general gc skew
		self.bases.append(base)
		self.frequency[base] += 1
		item = self.bases.popleft()
		self.frequency[item] -= 1
		self.total.append(self.frequency.copy())


	def _close(self):
		for _ in range(self.window//2):
			self.total.popleft()
			item = self.bases.popleft()
			self.frequency[item] -= 1
			self.total.append(self.frequency.copy())

	def get(self):
		self._close()
		return self.total

