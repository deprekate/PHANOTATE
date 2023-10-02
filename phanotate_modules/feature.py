from genbank.feature import Feature
import textwrap


class Feature(Feature):

	def more(self):
		return "mooore"

	def end(self):
		if self.strand > 0:
			return self.right()
		else:
			return self.left()
	
	def start_codon(self):
		return self.seq()[:3]

	def stop_codon(self):
		return self.seq()[-3:]

	def has_start(self):
		return self.start_codon() in self.locus.start_codons
	
	def has_stop(self):
		return self.stop_codon() in self.locus.stop_codons

	def set_left(self, left):
		if left:
			pairs = [list(tup) for tup in self.pairs]
			pairs[0][0] = str(left)
			self.pairs = tuple([tuple(lis) for lis in pairs])

	def set_right(self, right):
		if right:
			pairs = [list(tup) for tup in self.pairs]
			pairs[-1][-1] = str(right)
			self.pairs = tuple([tuple(lis) for lis in pairs])
