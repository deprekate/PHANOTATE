import os
import sys
import copy
import itertools
import tempfile
from subprocess import Popen, PIPE, STDOUT
from decimal import Decimal
from math import log10

from orfs import Orfs
from nodes import Node
from edges import Edge
from graphs import Graph
from gc_frame_plot import GCframe
from gc_frame_plot import max_idx
from gc_frame_plot import min_idx


def rev_comp(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G',
		    'N':'N',
		    'R':'Y','Y':'R','S':'S','W':'W','K':'M','M':'K',
		    'B':'V','V':'B','D':'H','H':'D'}
	return "".join([seq_dict[base] for base in reversed(seq)])

def score_orf(pstop, startcodon, length):
	score = (1-pstop)**(length/3)
	score = 1/score
	if(startcodon):
		score = score * start_weight[startcodon]
	return -score

def score_overlap(length, direction, pstop):
	o = Decimal(1-pstop)
	s = Decimal('0.05')

	score = Decimal(o)**Decimal(length)
	score = 1/score
	if(direction == 'diff'):
		score = score + (1/s)
	return score

def score_gap(length, direction, pgap):
	g = Decimal(1-pgap)
	s = Decimal('0.05')

	if(length > 300):
			return Decimal(g)**Decimal(100) + length
	score = Decimal(g)**Decimal(length/3)
	score = 1/score
	if(direction == 'diff'):
		score = score + (1/s)
	return score

def push(dict, key, value):
	try:
		dict[key].append(value)
	except:
		dict[key] = [value]
def comp(list1, list2):
	for val in list1:
		if val in list2:
			return True
	return False

def score_rbs(seq):
	score = 0
	kmers = dict()

	for length in range(3,7):
		for p in itertools.product('ACTG', repeat=length):
			kmers["".join(p)] = []
	
	for n in range(0,len(seq)-5):
		for x in range(3, 7):
			push(kmers, seq[n:n+x], len(seq)-n-x)
	push(kmers, seq[-3:len(seq)], 0)
	push(kmers, seq[-4:len(seq)], 0)
	push(kmers, seq[-5:len(seq)], 0)

	#These are motifs/scores taken from Prodigal
	if(comp(kmers['AGGAGG'],range(5,11))):
		score = 27
	elif(comp(kmers['AGGAGG'],range(3,5))):
		score = 26
	elif(comp(kmers['AGGAGG'],range(11,13))):
		score = 25
	elif(comp(kmers['GGAGG'],range(5,11))):
		score = 24
	elif(comp(kmers['GGAGG'],range(3,5))):
		score = 23
	elif(comp(kmers['AGGAG'],range(5,11))):
		score = 22
	elif(comp(kmers['AGGAG'],range(3,5))):
		score = 21
	elif(comp(kmers['AGGAG']+kmers['GGAGG'],range(3,5))):
		score = 20
	elif(comp(kmers['AGAAGG']+kmers['AGTAGG']+kmers['AGCAGG']+kmers['AGGTGG']+kmers['AGGCGG']+kmers['AGGGGG'],range(5,11))):
		score = 19
	elif(comp(kmers['AGAAGG']+kmers['AGTAGG']+kmers['AGCAGG']+kmers['AGGTGG']+kmers['AGGCGG']+kmers['AGGGGG'],range(3,5))):
		score = 18
	elif(comp(kmers['AGAAGG']+kmers['AGTAGG']+kmers['AGCAGG']+kmers['AGGTGG']+kmers['AGGCGG']+kmers['AGGGGG'],range(11,13))):
		score = 17
	elif(comp(kmers['GGAG']+kmers['GAGG'],range(5,11))):
		score = 16
	elif(comp(kmers['AGGA'],range(5,11))):
		score = 15
	elif(comp(kmers['GGTGG']+kmers['GGCGG']+kmers['GGGGG'],range(5,11))):
		score = 14
	elif(comp(kmers['GGA']+kmers['GAG']+kmers['AGG'],range(5,11))):
		score = 13
	elif(comp(kmers['AGGA']+kmers['GGAG']+kmers['GAGG'],range(11,13))):
		score = 12
	elif(comp(kmers['AGGA']+kmers['GGAG']+kmers['GAGG'],range(3,5))):
		score = 11
	elif(comp(kmers['AGGAG']+kmers['GGAGG']+kmers['AGGAGG'],range(13,16))):
		score = 10
	elif(comp(kmers['AGAAG']+kmers['AGTAG']+kmers['AGCAG'],range(5,11))):
		score = 9
	elif(comp(kmers['GGTGG']+kmers['GGCGG']+kmers['GGGGG'],range(3,5))):
		score = 8
	elif(comp(kmers['GGTGG']+kmers['GGCGG']+kmers['GGGGG'],range(11,13))):
		score = 7
	elif(comp(kmers['GGA']+kmers['GAG']+kmers['AGG'],range(11,13))):
		score = 6
	elif(comp(kmers['AGAAG']+kmers['AGTAG']+kmers['AGCAG'],range(3,5))):
		score = 5
	elif(comp(kmers['AGAAG']+kmers['AGTAG']+kmers['AGCAG'],range(12,13))):
		score = 4
	elif(comp(kmers['AGGA']+kmers['GGAG']+kmers['GAGG']+kmers['AGTAGG']+kmers['AGCAGG']+kmers['AGGTGG']+kmers['AGGCGG']+kmers['AGGGGG'],range(13,16))):
		score = 3
	elif(comp(kmers['GGA']+kmers['GAG']+kmers['AGG']+kmers['AGAAG']+kmers['AGTAG']+kmers['AGCAG']+kmers['GGTGG']+kmers['GGCGG']+kmers['GGGGG'],range(13,16))):
		score = 2
	elif(comp(kmers['GGA']+kmers['GAG']+kmers['AGG'],range(3,5))):
		score = 1
	return score

def p_stop(seq):
	frequency = {'A':Decimal(0), 'T':Decimal(0), 'C':Decimal(0), 'G':Decimal(0)}
	for base in seq:
		if(base not in ['A', 'C', 'T', 'G']):
			continue
		frequency[base] += 1
	Pa = frequency['A']/len(seq)
	Pt = frequency['T']/len(seq)
	Pg = frequency['G']/len(seq)
	Pc = frequency['C']/len(seq)
       	return (Pt*Pa*Pa + Pt*Pg*Pa + Pt*Pa*Pg)
def ave(a):
	return Decimal(sum(a)/len(a))

def get_orfs(dna):
	# This is the object that holds all the orfs
	my_orfs = Orfs()
	my_orfs.seq = dna
	my_orfs.contig_length = len(dna)

	# find nucleotide frequency, kmers, and create gc frame plot
	frequency = {'A':Decimal(0), 'T':Decimal(0), 'C':Decimal(0), 'G':Decimal(0)}
	frame_plot = GCframe()
	background_rbs = [1.0] * 28
	training_rbs = [1.0] * 28

	for i, base in enumerate(dna):
		if(base not in ['A', 'C', 'T', 'G']):
			continue
		#nucl frequency
		frequency[base] += 1
		frequency[rev_comp(base)] += 1
		#kmers for rbs
		background_rbs[score_rbs(dna[i:i+21])] += 1
		background_rbs[score_rbs(rev_comp(dna[i:i+21]))] += 1
		#gc frame plot
		frame_plot.add_base(base)


	Pa = frequency['A']/(my_orfs.contig_length*2)
	Pt = frequency['T']/(my_orfs.contig_length*2)
	Pg = frequency['G']/(my_orfs.contig_length*2)
	Pc = frequency['C']/(my_orfs.contig_length*2)
       	my_orfs.pstop = (Pt*Pa*Pa + Pt*Pg*Pa + Pt*Pa*Pg)

	gc_pos_freq = frame_plot.get()

	y = sum(background_rbs)
	background_rbs[:] = [x/y for x in background_rbs]

	#Pstarts = Pa*Pt*Pg + Pg*Pt*Pg + Pc*Pt*Pg
	#start_weight['ATG'] = (Decimal('0.85') * (Pa*Pt*Pg/Pstarts))/(Decimal('0.85') * (Pa*Pt*Pg/Pstarts))
	#start_weight['CAT'] = (Decimal('0.85') * (Pa*Pt*Pg/Pstarts))/(Decimal('0.85') * (Pa*Pt*Pg/Pstarts))
	#start_weight['GTG'] = (Decimal('0.10') * (Pg*Pt*Pg/Pstarts))/(Decimal('0.85') * (Pa*Pt*Pg/Pstarts))
	#start_weight['CAC'] = (Decimal('0.10') * (Pg*Pt*Pg/Pstarts))/(Decimal('0.85') * (Pa*Pt*Pg/Pstarts))
	#start_weight['TTG'] = (Decimal('0.05') * (Pc*Pt*Pg/Pstarts))/(Decimal('0.85') * (Pa*Pt*Pg/Pstarts))
	#start_weight['CAA'] = (Decimal('0.05') * (Pc*Pt*Pg/Pstarts))/(Decimal('0.85') * (Pa*Pt*Pg/Pstarts)b
	
	# The dicts that will hold the start and stop codons
	stops = {1:0, 2:0, 3:0, -1:0, -2:0, -3:0}
	starts = {1:[0], 2:[0], 3:[0], -1:[0], -2:[0], -3:[0]}


	# Reset iterator and find all the open reading frames
	states = itertools.cycle([1, 2, 3])
	for i in range(1, (len(dna)-1)):
		codon = dna[i-1:i+2]
		frame = states.next()

		if codon in ['ATG', 'GTG', 'TTG']:
			starts[frame].append(i)
		elif rev_comp(codon) in ['ATG', 'GTG', 'TTG']:
			starts[-frame].append(i)
		elif codon in ['TAA', 'TGA', 'TAG']:
			for start in reversed(starts[frame]):
				stop = i
				length = stop-start+3
				if(length >= my_orfs.min_orf_len):
					seq = dna[max(0,start-1):stop+2]
					rbs = dna[start-21:start]
					rbs_score = score_rbs(dna[start-21:start])
					my_orfs.add_orf(start, stop, length, frame, seq, rbs, rbs_score)
					training_rbs[rbs_score] += 1
	
			starts[frame] = []
			stops[frame] = i
		elif rev_comp(codon) in ['TAA', 'TGA', 'TAG']:
			for start in starts[-frame]:
				stop = stops[-frame]
				length = start-stop+3
				if(length >= my_orfs.min_orf_len):
					seq = rev_comp(dna[max(0,stop-1):start+2])
					rbs = rev_comp(dna[start+2:start+23])
					rbs_score = score_rbs(rev_comp(dna[start+2:start+23]))
					my_orfs.add_orf(start, stop, length, -frame, seq, rbs, rbs_score)
					training_rbs[rbs_score] += 1
	
			starts[-frame] = []
			stops[-frame] = i

	# Add in any fragment ORFs at the end of the genome
	for frame in [1, 2, 3]:
		for start in reversed(starts[frame]):
			stop = i+3
			length = stop-start
			if(length >= my_orfs.min_orf_len):
				seq = dna[max(0,start-1):stop+2]
				rbs = dna[start-21:start]
				rbs_score = score_rbs(dna[start-21:start])
				my_orfs.add_orf(start, stop, length, frame, seq, rbs, rbs_score)
				training_rbs[rbs_score] += 1

		starts[-frame].append(i+3)	
		for start in starts[-frame]:
			stop = stops[-frame]
			length = start-stop
			if(length >= my_orfs.min_orf_len):
				seq = rev_comp(dna[max(0,stop-1):start+2])
				rbs = rev_comp(dna[start+2:start+23])
				rbs_score = score_rbs(rev_comp(dna[start+2:start+23]))
				my_orfs.add_orf(start, stop, length, -frame, seq, rbs, rbs_score)
				training_rbs[rbs_score] += 1

	#-------------------------------Score ORFs based on RBS motif--------------------------------------#
	y = sum(training_rbs)
	training_rbs[:] = [x/y for x in training_rbs]
	for orf in my_orfs.iter_orfs():
		orf.weight_rbs = log10(10*(training_rbs[orf.rbs_score]/background_rbs[orf.rbs_score]))

	#-------------------------------Score ORFs based on GC frame plot----------------------------------#
	pos_max = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]
	pos_min = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]
	for orfs in my_orfs.iter_in():
		for orf in orfs:
			if(orf.start_codon() == 'ATG'):
				start = orf.start
				stop = orf.stop
				loc_pos_max = [Decimal(0), Decimal(0), Decimal(0), Decimal(0)]
				loc_pos_min = [Decimal(0), Decimal(0), Decimal(0), Decimal(0)]
				if(start < stop and stop-start):
					n = ((stop-start)/8)*3
					if(start == 0):
						start = (stop+2)%3+1
					for base in range(start+3+n, min(stop-33, len(dna)-1), 3):
						pos_max[max_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])] += 1
						pos_min[min_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])] += 1
				elif(stop and stop < start and start-stop):
					n = ((start-stop)/8)*3
					if(start >= len(dna)):
						start = len(dna)-(stop%3)-2
					for base in range(start-n, stop+33, -3):
						pos_max[max_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])] += 1
						pos_min[min_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])] += 1
				break

	# normalize to one
	y = max(pos_max)
	pos_max[:] = [x / y for x in pos_max]	
	y = max(pos_min)
	pos_min[:] = [x / y for x in pos_min]	

	for orf in my_orfs.iter_orfs():
		mins = []
		maxs = []
		start = orf.start
		stop = orf.stop
		if(orf.frame > 0):
			if(start == 0):
				start = orf.frame
			for base in range(start+3+0, min(stop-0,len(dna)-1), 3):
				ind_max = max_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])
				ind_min = min_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])
				orf.hold = orf.hold * ((((1-orf.pstop)**pos_max[ind_max]))**pos_min[ind_min])
				maxs.append(pos_max[ind_max])
				mins.append(pos_min[ind_min])
		else:
			if(start >= len(dna)):
				start = len(dna)-(stop%3)-2
			for base in range(start-0, stop+0, -3):
				ind_max = max_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])
				ind_min = min_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])
				orf.hold = orf.hold * ((((1-orf.pstop)**pos_max[ind_max]))**pos_min[ind_min])
				maxs.append(pos_max[ind_max])
				mins.append(pos_min[ind_min])
		orf.gcfp_min = ave(mins)
		orf.gcfp_max = ave(maxs)

	for orf in my_orfs.iter_orfs():
		orf.score()
		#print orf.start, orf.stop, orf.pstop, 1/orf.hold, "sep", orf.rbs, orf.weight_rbs, orf.gcfp_min, orf.gcfp_max, orf.weight
	#sys.exit()

	return my_orfs



def get_graph(my_orfs):
	G = Graph(directed=True)
	pgap = my_orfs.pstop

	for orf in my_orfs.iter_orfs():
		if(orf.frame > 0):
			source = Node('CDS', 'start', orf.frame, orf.start)
			target = Node('CDS', 'stop', orf.frame, orf.stop)
		else:
			source = Node('CDS', 'stop', orf.frame, orf.stop)
			target = Node('CDS', 'start', orf.frame, orf.start)
		G.add_edge(Edge(source, target, orf.weight))

	#-------------------------------Check for long noncoding regions that would break the path---------#
	bases = [None] * my_orfs.contig_length
	for orfs in my_orfs.iter_in():
		for orf in orfs:
			mi = min(orf.start, orf.stop)
			ma = max(orf.start, orf.stop)
			for n in range(mi, min(ma, my_orfs.contig_length-1)):
				try:
					bases[n] = n
				except:
					print n
			break
	last = 0
	for base in bases:
		if(base):
			if(base-last > 300):
				for right_node in G.iternodes():
					for left_node in G.iternodes():
						l = left_node.position
						r = right_node.position
						if(last+1 >= l > last-300 and base-1 <= r < base+300):
							if(left_node.frame*right_node.frame > 0):
								if(left_node.type == 'stop' and right_node.type =='start' and left_node.frame > 0):
									score = score_gap(r-l-3, 'same', pgap)
									G.add_edge(Edge(left_node, right_node, score))	
								elif(left_node.type == 'start' and right_node.type =='stop' and left_node.frame < 0):
									score = score_gap(r-l-3, 'same', pgap)
									G.add_edge(Edge(left_node, right_node, score ))	
							else:
								if(left_node.type == 'stop' and right_node.type =='stop' and left_node.frame > 0):
									score = score_gap(r-l-3, 'diff', pgap)
									G.add_edge(Edge(left_node, right_node, score ))	
								elif(left_node.type == 'start' and right_node.type =='start' and left_node.frame < 0):
									score = score_gap(r-l-3, 'diff',pgap)
									G.add_edge(Edge(left_node, right_node, score ))	
			last = base
	#-------------------------------Add in tRNA data---------------------------------------------------#
	
	#add_trnas(my_orfs.seq, G, start_to_stop, stop_to_start)
	
	#-------------------------------Connect the open reading frames to each other----------------------#
	for right_node in G.iternodes():
		r = right_node.position
		r_other = my_orfs.other_end[r]
		for left_node in G.iternodes():
			l = left_node.position
			l_other = my_orfs.other_end[l]
			if(0 < r-l < 300):
				if(l in my_orfs and my_orfs.other_end[l] in my_orfs[l]):
					o1 = my_orfs.get_orf(my_orfs.other_end[l], l)
				else:
					o1 = my_orfs.get_orf(l, my_orfs.other_end[l])
				if(r in my_orfs and my_orfs.other_end[r] in my_orfs[r]):
					o2 = my_orfs.get_orf(my_orfs.other_end[r], r)
				else:
					o2 = my_orfs.get_orf(r, my_orfs.other_end[r])
				pstop = ave([o1.pstop, o2.pstop])
			
				# same directions
				if(left_node.frame*right_node.frame > 0):
					if(left_node.type == 'stop' and right_node.type =='start'):
						if(left_node.frame > 0):
							score = score_gap(r-l-3, 'same', pgap)
							G.add_edge(Edge(left_node, right_node, score ))	
						else:
							if(left_node.frame != right_node.frame):
								if(r < l_other and r_other < l):
									score = score_overlap(r-l+3, 'same', pstop)
									G.add_edge(Edge(right_node, left_node, score ))	
					if(left_node.type == 'start' and right_node.type =='stop'):
						if(left_node.frame > 0):
							if(left_node.frame != right_node.frame):
								if(r < l_other and r_other < l):
									score = score_overlap(r-l+3, 'same', pstop)
									G.add_edge(Edge(right_node, left_node, score ))	
						else:
							score = score_gap(r-l-3, 'same', pgap)
							G.add_edge(Edge(left_node, right_node, score ))	
				# different directions
				else:
					if(left_node.type == 'stop' and right_node.type =='stop'):
						if(right_node.frame > 0):
							if(r_other < l and r < l_other):
								score = score_overlap(r-l+3, 'diff', pstop)
								G.add_edge(Edge(right_node, left_node, score ))	
						else:
							score = score_gap(r-l-3, 'diff', pgap)
							G.add_edge(Edge(left_node, right_node, score ))	
					if(left_node.type == 'start' and right_node.type =='start'):
						if(right_node.frame > 0 and r-l > 2):
							score = score_gap(r-l-3, 'diff', pgap)
							G.add_edge(Edge(left_node, right_node, score ))	
						elif(right_node.frame < 0):
							if(r_other < l and r < l_other):
								score = score_overlap(r-l+3, 'diff', pstop)
								G.add_edge(Edge(right_node, left_node, score ))	
	#-------------------------------Connect open reading frames at both ends to a start and stop-------#
	source = Node('source', 'source', 0, 0)
	target = Node('target', 'target', 0, my_orfs.contig_length+1)
	G.add_node(source)
	G.add_node(target)
	for node in G.iternodes():
		if(node.position <= 2000):
			if( (node.type == 'start' and node.frame > 0) or (node.type =='stop' and node.frame < 0) ):
				score = score_gap(node.position, 'same', pgap)
				G.add_edge(Edge(source, node, score))
		if(my_orfs.contig_length - node.position <= 2000):
			if( (node.type == 'start' and node.frame < 0) or (node.type =='stop' and node.frame > 0) ):
				score = score_gap(my_orfs.contig_length-node.position, 'same', pgap)
				G.add_edge(Edge(node, target, score))

	return G
#---------------------------------------END OF LOOP----------------------------------------------------------#

def add_trnas(seq, G, start_to_stop, stop_to_start):
	trna_end_left = {}
	trna_end_right = {}

	f = tempfile.NamedTemporaryFile()
	f.write(">temp\n")
	f.write(seq)
	f.seek(0)

	try:
		output = Popen(["tRNAscan-SE", "-B", "-q", "-b", f.name], stdout=PIPE, stdin=PIPE, stderr=PIPE).stdout.read()
	except:
		sys.stderr.write("Warning: tRNAscan not found, proceding without tRNA masking.\n")
		return []

	# Iterate over the trnas
	for line in output.splitlines():
		# Add in trna
		column = line.split('\t')
		start = int(column[2])
		stop = int(column[3])
		if(start < stop):
			source = Node('tRNA', 'start', 4, start)
			target = Node('tRNA', 'stop', 4, stop)
		else:
			source = Node('tRNA', 'start', -4, start)
			target = Node('tRNA', 'stop', -4, stop)
		G.add_edge(Edge(source, target, -Decimal(100)))
		start_to_stop[start] = stop
		stop_to_start[stop] = start


