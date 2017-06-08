import sys
import os.path
import itertools
import settings
from nodes import Node
from settings import weights
import argparse
from argparse import RawTextHelpFormatter

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
	parser.add_argument('-f', '--outfmt', action="store", default="tabular", dest='outfmt', help='format of the output [tabular]', choices=['tabular','genbank'])
	parser.add_argument('-m', '--min_orf_len', action="store", default=60, dest="min_orf_len", type=int, help='the minimum open reading frame length [60]')

	parser.add_argument('--overlap_penalty', action="store", dest="overlap_penalty", type=float, choices=[Range(0.0,1.0)], help='')
	parser.add_argument('--gap_penalty', action="store", dest="gap_penalty", type=float, choices=[Range(0.0,1.0)], help='')
	parser.add_argument('--switch_penalty', action="store", dest="switch_penalty", type=float, choices=[Range(0.0,1.0)], help='')


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

def write_output(id, args, my_path):
	outfmt = args.outfmt
	outfile = args.outfile

	if(outfmt == 'tabular'):
		outfile.write("#gap:\t" + str(weights['gap']) + "\n")
		outfile.write("#overlap:\t" + str(weights['overlap']) + "\n")
		outfile.write("#switch:\t" + weights['switch'] + "\n")
		outfile.write("#min_orf_length:\t" + weights['min_orf_length'] + "\n")
		for source, target in pairwise(my_path):
			left = eval(source)
			right = eval(target)
			if(left.type == 'start' and right.type == 'stop'):
				outfile.write(str(left.position) + '\t' + str(right.position+2) + '\t+\n')
			elif(left.type == 'stop' and right.type == 'start'):
				outfile.write(str(left.position) + '\t' + str(right.position+2) + '\t-\n')

	elif(outfmt == 'genbank'):
		last_node = eval(my_path[-1])
		outfile.write('LOCUS       UNKNOWN')
		outfile.write(str(last_node.position-1).rjust(10))
		outfile.write(' bp    DNA             PHG\n')
		outfile.write('DEFINITION  ' + id + '\n')
		outfile.write('FEATURES             Location/Qualifiers\n')
		outfile.write('     source          1..' + str(last_node.position-1) + '\n')
		for source, target in pairwise(my_path):
			left = eval(source)
			right = eval(target)
			if(left.position == 0):
				left.position = '<1'
			if(right.position == last_node.position):
				right.position = '>' + str(right.position-1)
			else:
				right.position += 2

			outfile.write('     CDS             ')
			if(left.type == 'start' and right.type == 'stop'):
				outfile.write(str(left.position) + '..' + str(right.position) + '\n')
			elif(left.type == 'stop' and right.type == 'start'):
				outfile.write('complement(' + str(left.position) + '..' + str(right.position) + ')\n')
				#outfile.write('                     /color=100 100 100\n')
		outfile.write('//\n')


