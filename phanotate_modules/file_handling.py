import sys
import os.path
import itertools
import textwrap
import argparse
from argparse import RawTextHelpFormatter
from math import log


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

	args = parser.parse_args()

	args.min_orf_len  = 87
	args.start_codons = ['ATG' ,'GTG', 'TTG']
	args.stop_codons  = ['TAG' ,'TGA', 'TAA']
	args.start_weight = {
						 'ATG':1.00, 'CAT':1.00,
						 'GTG':0.12, 'CAC':0.12,
						 'TTG':0.05, 'CAA':0.05
						}

	return args


def read_fasta(filepath, base_trans=str.maketrans('','')):
	contigs_dict = dict()
	name = ''
	seq = ''
	with open(filepath, mode="r") as f:
		for line in f:
			if line.startswith(">"):
				contigs_dict[name] = seq
				name = line[1:].split()[0]
				seq = ''
			else:
				#seq += line.replace("\n", "").upper()
				seq += line[:-1].upper()
		contigs_dict[name] = seq.translate(base_trans)

	if '' in contigs_dict: del contigs_dict['']
	return contigs_dict

def write_output(contig_orfs, shortest_path):
	if   contig_orfs.outfmt == 'tabular':
		write_tabular(contig_orfs, shortest_path)
	elif contig_orfs.outfmt == 'genbank':
		write_genbank(contig_orfs, shortest_path)
	elif contig_orfs.outfmt == 'fasta':
		write_fasta(contig_orfs, shortest_path)


def write_tabular(contig_orfs, shortest_path):
	outfile = contig_orfs.outfile

	outfile.write("#id:\t" + contig_orfs.id + "\n")
	outfile.write("#START\tSTOP\tFRAME\tCONTIG\tSCORE\n")
	for left, right in pairwise(shortest_path):
		feature = contig_orfs.get_feature(left, right)
		if(type(feature).__name__  != 'CDS'):
			continue

		outfile.write('\t'.join(map(str, [feature.begin(), feature.end(), feature.direction(), contig_orfs.id, feature.weight, '\n'] )))

def write_genbank(contig_orfs, shortest_path):
	outfile = contig_orfs.outfile

	outfile.write('LOCUS       ')
	outfile.write(contig_orfs.id)
	outfile.write(str(contig_orfs.contig_length()).rjust(10))
	outfile.write(' bp    DNA             PHG')
	outfile.write('\n')
	outfile.write('DEFINITION  ' + contig_orfs.id + '\n')
	outfile.write('FEATURES             Location/Qualifiers\n')
	outfile.write('     source          1..')
	outfile.write(str(contig_orfs.contig_length()))
	outfile.write('\n')
	for left, right in pairwise(shortest_path):
		feature = contig_orfs.get_feature(left, right)
		outfile.write('     ' + type(feature).__name__.ljust(16))
		if feature.frame > 0:
			outfile.write(feature.begin() + '..' + feature.end() + '\n')
		else:
			outfile.write('complement(' + feature.end() + '..' + feature.begin() + ')\n')
		outfile.write('                     /note="weight=' + '{:.2E}'.format(feature.weight) + ';"\n')
	outfile.write('ORIGIN')
	return
	i = 0
	dna = textwrap.wrap(contig_orfs.dna, 10)
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

def write_fasta(contig_orfs, shortest_path):
	outfile = contig_orfs.outfile
	for left, right in pairwise(shortest_path):
		feature = contig_orfs.get_feature(left, right)

		if(feature.type == 'CDS'):
			#outfile.write('\t'.join(map(str, [feature.begin(), feature.end(), feature.direction(), contig_orfs.id, feature.weight, '\n'] )))
			#outfile.write(str(left.position) + '\t' + str(right.position) + '\t+\t' + id[1:] + '\t' + str(weight) + '\t\n')
			outfile.write(">" + contig_orfs.id + "_" + feature.end().strip("<>"))
			outfile.write(" [START=" + feature.begin().strip("<>") + "]")
			outfile.write(" [STOP=" + feature.end().strip("<>") + "]")
			outfile.write(" [SCORE=" + '{:.2E}'.format(feature.weight) + "]")
			outfile.write("\n")
			outfile.write(feature.dna)
			outfile.write("\n")


