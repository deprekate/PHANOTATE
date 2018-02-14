import sys
import os.path
import itertools

from edges import Edge
from nodes import Node
from math import log
import argparse
from argparse import RawTextHelpFormatter
from decimal import Decimal

def pairwise(iterable):
	a = iter(iterable)
	return itertools.izip(a, a)

class Range(object):
	def __init__(self, start, end):
		self.start = start
		self.end = end
	def __repr__(self):
		return '{0}-{1}'.format(self.start, self.end)
	def __eq__(self, other):
		return self.start <= other <= self.end
def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def get_args():
	usage = 'thea.py [-opt1, [-opt2, ...]] infile'
	parser = argparse.ArgumentParser(description='THEA: A phage genome annotator', formatter_class=RawTextHelpFormatter, usage=usage)

	parser.add_argument('infile', type=is_valid_file, help='input file in fasta format')

	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	parser.add_argument('-f', '--outfmt', action="store", default="tabular", dest='outfmt', help='format of the output [tabular]', choices=['tabular','genbank','fasta'])

	args = parser.parse_args()

	return args


def read_fasta(filepath):
	my_contigs = dict()
	name = ''
	seq = ''
	with open(filepath, mode="r") as my_file:
		for line in my_file:
			if(line.startswith(">")):
				my_contigs[name] = seq
				name = line.split()[0]
				seq = ''
			else:
				seq += line.replace("\n", "").upper()
		my_contigs[name] = seq

	if '' in my_contigs: del my_contigs['']
	return my_contigs

def write_output(id, args, my_path, my_graph, G):
	outfmt = args.outfmt
	outfile = args.outfile

	if(outfmt == 'tabular'):
		last_node = eval(my_path[-1])
		outfile.write("#id:\t" + str(id[1:]) + "\n")
		outfile.write("#START\tSTOP\tFRAME\tCONTIG\tSCORE\n")
		cutoff = -1/((1-G.pstop)**30)/3
		#cutoff = -Decimal("1")
		for source, target in pairwise(my_path):
			left = eval(source)
			right = eval(target)
			if(left.gene == 'CDS'):
				weight = my_graph.weight(Edge(left,right,0))
				length = abs(right.position-left.position)/3
				if(weight > cutoff):
					continue
				if(left.type == 'start' and right.type == 'stop'):
					outfile.write(str(left.position) + '\t' + str(right.position+2) + '\t+\t' + id[1:] + '\t' + str(weight) + '\t\n')
				elif(left.type == 'stop' and right.type == 'start'):
					outfile.write(str(right.position+2) + '\t' + str(left.position) + '\t-\t' + id[1:] + '\t' + str(weight) + '\t\n')

	elif(outfmt == 'genbank'):
		last_node = eval(my_path[-1])
		outfile.write('LOCUS       UNKNOWN')
		outfile.write(str(last_node.position-1).rjust(10))
		outfile.write(' bp    DNA             PHG\n')
		outfile.write('DEFINITION  ' + id[1:] + '\n')
		outfile.write('FEATURES             Location/Qualifiers\n')
		outfile.write('     source          1..' + str(last_node.position-1) + '\n')
		for source, target in pairwise(my_path):
			left = eval(source)
			right = eval(target)
			outfile.write('     ' + left.gene.ljust(17))
			if(left.type == 'start' and right.type == 'stop'):
				outfile.write(str(left.position) + '..' + str(right.position+2) + '\n')
			elif(left.type == 'stop' and right.type == 'start'):
				outfile.write('complement(' + str(left.position) + '..' + str(right.position+2) + ')\n')
		outfile.write('//\n')
	elif(outfmt == 'fasta'):
		last_node = eval(my_path[-1])
		cutoff = -1/((1-G.pstop)**30)/3
		for source, target in pairwise(my_path):
			left = eval(source)
			right = eval(target)
			if(left.gene == 'CDS'):
				weight = my_graph.weight(Edge(left,right,0))
				if(weight > cutoff):
					continue
				if(left.frame > 0):
					o = G.get_orf(left.position, right.position)
					if(right.position == last_node.position):
						right.position = left.position+3*int((right.position-left.position)/3)-1
						outfile.write(id + "." + str(right.position) + " START=" + str(o.start) + " SCORE=" + str(weight) + "\n")
					else:
						outfile.write(id + "." + str(o.stop+2) + " START=" + str(o.start) + " SCORE=" + str(weight) + "\n")
				else:
					o = G.get_orf(right.position, left.position)
					if(left.position == 0):
						left.position = ((right.position+2)%3)+1
						outfile.write(id + "." + str(left.position) + " START=" + str(o.start) + " SCORE=" + str(weight) + "\n")
					else:
						outfile.write(id + "." + str(o.stop) + " START=" + str(o.start) + " SCORE=" + str(weight) + "\n")
				outfile.write(o.seq)
				outfile.write("\n")


