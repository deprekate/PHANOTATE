from collections import deque
from decimal import Decimal
class GCskew:
	def __init__(self):
		self.bases = deque([])
		self.frequency = {'A':Decimal(0), 'T':Decimal(0), 'C':Decimal(0), 'G':Decimal(0)}
		#for base in seq:
		#	self.bases.append(base)
		#	self.frequency[base] += 1
	def add_base(self, base):
		self.bases.append(base)
		self.frequency[base] += 1
		if(len(self.bases) > 100):
			item = self.bases.popleft()
		#if item:
			self.frequency[item] -= 1
	def p_stop(self):
		Pa = self.frequency['A']/100
		Pt = self.frequency['T']/100
		Pg = self.frequency['G']/100
		Pc = self.frequency['C']/100
        	return 1-(Pt*Pa*Pa + Pt*Pg*Pa + Pt*Pa*Pg)

