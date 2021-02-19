import sys
import tempfile
from subprocess import Popen, PIPE, STDOUT
from decimal import Decimal

class tRNAs(list):
	def __init__(self, dna):
		#self.__dict__.update(kwargs)
		self.n = 0
		self.dna = dna

		f = tempfile.NamedTemporaryFile(mode='wt')
		f.write(">temp\n")
		f.write(dna)
		f.seek(0)

		try:
			output = Popen(["tRNAscan-SE", "-B", "-q", "--brief", f.name], stdout=PIPE, stdin=PIPE, stderr=PIPE).stdout.read()
		except:
			sys.stderr.write("Warning: tRNAscan not found, proceding without tRNA masking.\n")
			return []

		# Iterate over the trnas
		for line in output.decode().splitlines():
			# Add in trna
			name, num, start, stop, type_, *_ = line.split('\t')
			if(start < stop):
				self.append(tRNA(start, stop, type_))
			else:
				self.append(tRNA(start, stop, type_))



class tRNA: 
	def __init__(self, start, stop, type_):
		self.start  = int(start)
		self.stop    = int(stop)
		self.type   = type_
		self.weight = Decimal(-1)
		self.frame  = 4 if start < stop else -4

	def as_edge(self):
		if self.frame > 0:
			return ("%s_begin" % self.left(), "%s_end" % (self.right()-2), "%s" % self.weight)
		else:
			return ("%s_end" % self.left(), "%s_begin" % (self.right()-2), "%s" % self.weight)

	def left(self):
		if self.frame > 0:
			return self.start
		else:
			return self.stop

	def right(self):
		if self.frame > 0:
			return self.stop
		else:
			return self.start

	def begin(self):
		if self.frame > 0:
			return str(self.start)
		else:
			return str(self.start)

	def end(self):
		if self.frame > 0:
			return str(self.stop)
		else:
			return str(self.stop)

	def __repr__(self):
		"""Compute the string representation of the orf"""
		return "%s(%s,%s,%s,%s)" % (
			self.__class__.__name__,
			repr(self.start),
			repr(self.stop),
			repr(self.frame),
			repr(self.weight))

