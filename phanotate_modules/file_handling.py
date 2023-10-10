import sys
import os.path
import itertools
import textwrap
import argparse
from argparse import RawTextHelpFormatter
from math import log
from decimal import Decimal

from phanotate_modules.edges import Edge
from phanotate_modules.nodes import Node
from phanotate_modules.functions import rev_comp


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

def get_args(File):
	usage = 'phanotate.py [-opt1, [-opt2, ...]] infile'
	parser = argparse.ArgumentParser(description='PHANOTATE: A phage genome annotator', formatter_class=RawTextHelpFormatter, usage=usage)

	parser.add_argument('infile', type=is_valid_file, help='input file in fasta format')

	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	#parser.add_argument('--outfmt', action="store", default="tabular", dest='outfmt', help='format of the output [tabular]', choices=['tabular','genbank','fasta'])
	parser.add_argument('-f', '--format', help='Output the features in the specified format [tabular]', type=str, default='tabular', choices=File.formats[:7])
	parser.add_argument('-s', '--start_codons', action="store", default="atg:0.85,gtg:0.10,ttg:0.05", dest='start_codons', help='comma separated list of start codons and frequency [atg:0.85,gtg:0.10,ttg:0.05]')
	parser.add_argument('-e', '--stop_codons', action="store", default="tag,tga,taa", dest='stop_codons', help='comma separated list of stop codons [tag,tga,taa]')
	parser.add_argument('-l', '--minlen', action="store", type=int, default=90, dest='min_orf_len', help='to store a variable')
	parser.add_argument('-d', '--dump', action="store_true")
	parser.add_argument('-V', '--version', action='version', version=__version__)
	args = parser.parse_args()

	start_codons = dict()
	for codon,weight in map(lambda x: tuple(x.split(':')), args.start_codons.split(',')):
		start_codons[codon.lower()] = Decimal(weight)
	start_codons = {k: v / m for m in (max(start_codons.values()),) for k, v in start_codons.items()}
	args.start_codons = start_codons
	stop_codons = []
	for codon in args.stop_codons.split(','):
		stop_codons.append(codon.lower())
	args.stop_codons = stop_codons

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
				seq += line.replace("\n", "").lower()
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
		outfile.write("#id:\t" + str(id) + " NO ORFS FOUND\n")
	elif(outfmt == 'tabular'):
		last_node = eval(my_path[-1])
		outfile.write("#id:\t" + str(id) + "\n")
		outfile.write("#START\tSTOP\tFRAME\tCONTIG\tSCORE\n")
		for source, target in pairwise(my_path):
			left = eval(source)
			right = eval(target)
			weight = my_graph.weight(Edge(left,right,0))
			if(left.gene == 'tRNA'):
				continue
			if(left.frame > 0):
				orf = my_orfs.get_orf(left.position, right.position)
			else:
				orf = my_orfs.get_orf(right.position, left.position)
			if(left.frame > 0 and not orf.has_start()) or (left.frame < 0 and not orf.has_stop()):
				left.position = '<1'
			if(left.frame > 0 and not orf.has_stop()) or (left.frame < 0 and not orf.has_start()):
				right.position = '>' + str(last_node.position-1)
			else:
				right.position += 2
			if(left.type == 'start' and right.type == 'stop'):
				outfile.write(str(left.position) + '\t' + str(right.position) + '\t+\t' + id + '\t' + str(weight) + '\t\n')
			elif(left.type == 'stop' and right.type == 'start'):
				outfile.write(str(right.position) + '\t' + str(left.position) + '\t-\t' + id + '\t' + str(weight) + '\t\n')

	elif(outfmt == 'genbank'):
		last_node = eval(my_path[-1])
		outfile.write('LOCUS       ' + id)
		outfile.write(str(last_node.position-1).rjust(10))
		outfile.write(' bp    DNA             PHG\n')
		outfile.write('DEFINITION  ' + id + '\n')
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
			if(left.frame > 0 and not orf.has_start()) or (left.frame < 0 and not orf.has_stop()):
				left.position = '<1'
			if(left.frame > 0 and not orf.has_stop()) or (left.frame < 0 and not orf.has_start()):
				right.position = '>' + str(last_node.position-1)
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
			if(left.gene != 'CDS'):
				continue
			weight = my_graph.weight(Edge(left,right,0))
			if(left.frame > 0):
				orf = my_orfs.get_orf(left.position, right.position)
			else:
				orf = my_orfs.get_orf(right.position, left.position)
			if(not orf.has_start() and not orf.has_stop()):
				pass
			elif(left.frame > 0 and not orf.has_start()):
				orf.seq = my_orfs.seq[:left.position-1] + orf.seq
				left.position = '<1'
			elif(left.frame < 0 and not orf.has_stop()):
				orf.seq = orf.seq + rev_comp(my_orfs.seq[:left.position-1])
				left.position = '<1'
			elif(left.frame > 0 and not orf.has_stop()):
				orf.seq += my_orfs.seq[right.position+2:]
				right.position = '>' + str(last_node.position-1)
			elif(left.frame < 0 and not orf.has_start()):
				orf.seq = rev_comp(my_orfs.seq[right.position+2:]) + orf.seq
				right.position = '>' + str(last_node.position-1)
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


