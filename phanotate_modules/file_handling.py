import sys
import os.path
import itertools
import textwrap
import argparse
from argparse import RawTextHelpFormatter
from math import log

from phanotate_modules.edges import Edge
from phanotate_modules.nodes import Node


import pkg_resources

try:
    __version__ = pkg_resources.get_distribution('phanotate').version
except Exception:
    __version__ = 'unknown'



def pairwise(iterable):
	a = iter(iterable)
	return zip(a, a)

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
	usage = 'phanotate.py [-opt1, [-opt2, ...]] infile'
	parser = argparse.ArgumentParser(description='PHANOTATE: A phage genome annotator', formatter_class=RawTextHelpFormatter, usage=usage)

	parser.add_argument('infile', type=is_valid_file, help='input file in fasta format')

	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	parser.add_argument('-f', '--outfmt', action="store", default="tabular", dest='outfmt', help='format of the output [tabular]', choices=['tabular','genbank','fasta'])
	parser.add_argument('-d', '--dump', action="store_true")
	parser.add_argument('-V', '--version', action='version', version=__version__)

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

def write_output(id, args, my_path, my_graph, my_orfs):
	outfmt = args.outfmt
	outfile = args.outfile
	
	try:
		my_path = my_path[1:]
	except:
		sys.stdout.write("Error running fastpathz: " + output + '\n')
	if(not my_path):
		outfile.write("#id:\t" + str(id[1:]) + " NO ORFS FOUND\n")
	elif(outfmt == 'tabular'):
		last_node = eval(my_path[-1])
		outfile.write("#id:\t" + str(id[1:]) + "\n")
		outfile.write("#START\tSTOP\tFRAME\tCONTIG\tSCORE\n")
		for source, target in pairwise(my_path):
			left = eval(source)
			right = eval(target)
			weight = my_graph.weight(Edge(left,right,0))
			if(left.gene == 'tRNA'):
				continue
			if(left.position == 0 and right.position == last_node.position):
				left.position = abs(left.frame)
				right.position = '>' + str(left.position+3*int((right.position-left.position)/3)-1)
				left.position = '<' + str(left.position)
			elif(left.position == 0):
				left.position = '<' + str(((right.position+2)%3)+1)
				right.position += 2
			elif(right.position == last_node.position):
				right.position = '>' + str(left.position+3*int((right.position-left.position)/3)-1)
			else:
				right.position += 2
			if(left.type == 'start' and right.type == 'stop'):
				outfile.write(str(left.position) + '\t' + str(right.position) + '\t+\t' + id[1:] + '\t' + str(weight) + '\t\n')
			elif(left.type == 'stop' and right.type == 'start'):
				outfile.write(str(right.position) + '\t' + str(left.position) + '\t-\t' + id[1:] + '\t' + str(weight) + '\t\n')

	elif(outfmt == 'genbank'):
		last_node = eval(my_path[-1])
		outfile.write('LOCUS       ' + id[1:])
		outfile.write(str(last_node.position-1).rjust(10))
		outfile.write(' bp    DNA             PHG\n')
		outfile.write('DEFINITION  ' + id[1:] + '\n')
		outfile.write('FEATURES             Location/Qualifiers\n')
		outfile.write('     source          1..' + str(last_node.position-1) + '\n')
		for source, target in pairwise(my_path):
			#get the orf
			left = eval(source)
			right = eval(target)
			weight = my_graph.weight(Edge(left,right,0))
			if(left.gene == 'tRNA' or right.gene == 'tRNA'):
				outfile.write('     ' + left.gene.ljust(16))
				if(left.frame > 0):
					outfile.write(str(left.position) + '..' + str(right.position) + '\n')
				else:
					outfile.write('complement(' + str(left.position) + '..' + str(right.position) + ')\n')
				continue
				
			if(left.frame > 0):
				orf = my_orfs.get_orf(left.position, right.position)
			else:
				orf = my_orfs.get_orf(right.position, left.position)
			#properly display the orf
			if(not orf.has_start() and not orf.has_stop()):
				left.position = '<' + str(((right.position+2)%3)+1)
			if(right.position == last_node.position):
				right.position = '>' + str(left.position+3*int((right.position-left.position)/3)-1)
			else:
				right.position += 2
			outfile.write('     ' + left.gene.ljust(16))
			if(left.type == 'start' and right.type == 'stop'):
				outfile.write(str(left.position) + '..' + str(right.position) + '\n')
			elif(left.type == 'stop' and right.type == 'start'):
				outfile.write('complement(' + str(left.position) + '..' + str(right.position) + ')\n')
			outfile.write('                     /note="weight=' + '{:.2E}'.format(weight) + ';"\n')
		outfile.write('ORIGIN')
		i = 0
		dna = textwrap.wrap(my_orfs.seq, 10)
		for block in dna:
			if(i%60 == 0):
				outfile.write('\n')
				outfile.write(str(i+1).rjust(9))
				outfile.write(' ')
				outfile.write(block.lower())
			else:
				outfile.write(' ')
				outfile.write(block.lower())
			i += 10
			
		outfile.write('\n')	
		outfile.write('//')
		outfile.write('\n')
	elif(outfmt == 'fasta'):
		last_node = eval(my_path[-1])
		for source, target in pairwise(my_path):
			left = eval(source)
			right = eval(target)
			if(left.gene == 'tRNA'):
				continue
			if(left.frame > 0):
				orf = my_orfs.get_orf(left.position, right.position)
			else:
				orf = my_orfs.get_orf(right.position, left.position)
			if(left.gene == 'CDS'):
				weight = my_graph.weight(Edge(left,right,0))
				if(left.position == 0 and right.position == last_node.position):
					left.position = abs(left.frame)
					right.position = '>' + str(left.position+3*int((right.position-left.position)/3)-1)
					left.position = '<' + str(left.position)
				elif(left.position == 0):
					left.position = '<' + str(((right.position+2)%3)+1)
					right.position += 2
				elif(right.position == last_node.position):
					right.position = '>' + str(left.position+3*int((right.position-left.position)/3)-1)
				else:
					right.position += 2
				if(left.type == 'start' and right.type == 'stop'):
					#outfile.write(str(left.position) + '\t' + str(right.position) + '\t+\t' + id[1:] + '\t' + str(weight) + '\t\n')
					outfile.write(id + "." + str(right.position) + " [START=" + str(left.position) + "] [SCORE=" + str(weight) + "]\n")
				elif(left.type == 'stop' and right.type == 'start'):
					#outfile.write(str(right.position) + '\t' + str(left.position) + '\t-\t' + id[1:] + '\t' + str(weight) + '\t\n')
					outfile.write(id + "." + str(left.position) + " [START=" + str(right.position) + "] [SCORE=" + str(weight) + "]\n")
				outfile.write(orf.seq)
				outfile.write("\n")


