import os
import sys
import copy
import itertools
import tempfile
from subprocess import Popen, PIPE, STDOUT
from decimal import Decimal

from .orfs import Orfs
from .nodes import Node
from .edges import Edge
from .graphs import Graph
from .gc_frame_plot import GCframe
from .gc_frame_plot import max_idx
from .gc_frame_plot import min_idx
#from kmeans import kmeans


def rev_comp(seq):
	seq_dict = {'a':'t','t':'a','g':'c','c':'g',
		    'n':'n',
		    'r':'y','y':'r','s':'s','w':'w','k':'m','m':'k',
		    'b':'v','v':'b','d':'h','h':'d'}
	return "".join([seq_dict[base] for base in reversed(seq)])

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

def score_rbs(seq):
	s = seq[::-1]
	score = 0

	if 'ggagga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 27
	elif 'ggagga' in (s[3:9],s[4:10]):
		score = 26
	elif 'ggagga' in (s[11:17],s[12:18]):
		score = 25
	elif 'ggagg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 24
	elif 'ggagg' in (s[3:8],s[4:9]):
		score = 23
	elif 'gagga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 22
	elif 'gagga' in (s[3:8],s[4:9]):
		score = 21
	elif 'gagga' in (s[11:16],s[12:17]) or 'ggagg' in (s[11:16],s[12:17]):
		score = 20
	elif 'ggacga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggatga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggaaga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggcgga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggggga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggtgga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggaaga' in (s[3:9],s[4:10]) or 'ggatga' in (s[3:9],s[4:10]) or 'ggacga' in (s[3:9],s[4:10]):
		score = 18
	elif 'ggtgga' in (s[3:9],s[4:10]) or 'ggggga' in (s[3:9],s[4:10]) or 'ggcgga' in (s[3:9],s[4:10]):
		score = 18
	elif 'ggaaga' in (s[11:17],s[12:18]) or 'ggatga' in (s[11:17],s[12:18]) or 'ggacga' in (s[11:17],s[12:18]):
		score = 17
	elif 'ggtgga' in (s[11:17],s[12:18]) or 'ggggga' in (s[11:17],s[12:18]) or 'ggcgga' in (s[11:17],s[12:18]):
		score = 17
	elif 'ggag' in (s[5:9],s[6:10],s[7:11],s[8:12],s[9:13],s[10:14]):
		score = 16
	elif 'gagg' in (s[5:9],s[6:10],s[7:11],s[8:12],s[9:13],s[10:14]):
		score = 16
	elif 'agga' in (s[5:9],s[6:10],s[7:11],s[8:12],s[9:13],s[10:14]):
		score = 15
	elif 'ggtgg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 14
	elif 'ggggg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 14
	elif 'ggcgg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 14
	elif 'agg' in (s[5:8],s[6:9],s[7:10],s[8:11],s[9:12],s[10:13]):
		score = 13
	elif 'gag' in (s[5:8],s[6:9],s[7:10],s[8:11],s[9:12],s[10:13]):
		score = 13
	elif 'gga' in (s[5:8],s[6:9],s[7:10],s[8:11],s[9:12],s[10:13]):
		score = 13
	elif 'agga' in (s[11:15],s[12:16]) or 'gagg' in (s[11:15],s[12:16]) or 'ggag' in (s[11:15],s[12:16]):
		score = 12
	elif 'agga' in (s[3:7],s[4:8]) or 'gagg' in (s[3:7],s[4:8]) or 'ggag' in (s[3:7],s[4:8]):
		score = 11
	elif 'gagga' in (s[13:18],s[14:19],s[15:20]) or 'ggagg' in (s[13:18],s[14:19],s[15:20]) or 'ggagga' in (s[13:19],s[14:20],s[15:21]):
		score = 10
	elif 'gaaga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 9
	elif 'gatga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 9
	elif 'gacga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 9
	elif 'ggtgg' in (s[3:8],s[4:9]) or 'ggggg' in (s[3:8],s[4:9]) or 'ggcgg' in (s[3:8],s[4:9]):
		score = 8
	elif 'ggtgg' in (s[11:16],s[12:17]) or 'ggggg' in (s[11:16],s[12:17]) or 'ggcgg' in (s[11:16],s[12:17]):
		score = 7
	elif 'agg' in (s[11:14],s[12:15]) or 'gag' in (s[11:14],s[12:15]) or 'gga' in (s[11:14],s[12:15]):
		score = 6
	elif 'gaaga' in (s[3:8],s[4:9]) or 'gatga' in (s[3:8],s[4:9]) or 'gacga' in (s[3:8],s[4:9]):
		score = 5
	elif 'gaaga' in (s[11:16],s[12:17]) or 'gatga' in (s[11:16],s[12:17]) or 'gacga' in (s[11:16],s[12:17]):
		score = 4
	elif 'agga' in (s[13:17],s[14:18],s[15:19]) or 'gagg' in (s[13:17],s[14:18],s[15:19]) or 'ggag' in (s[13:17],s[14:18],s[15:19]):
		score = 3
	elif 'agg' in (s[13:16],s[14:17],s[15:18]) or 'gag' in (s[13:16],s[14:17],s[15:18]) or 'gga' in (s[13:16],s[14:17],s[15:18]):
		score = 2
	elif 'ggaaga' in (s[13:19],s[14:20],s[15:21]) or 'ggatga' in (s[13:19],s[14:20],s[15:21]) or 'ggacga' in (s[13:19],s[14:20],s[15:21]):
		score = 2
	elif 'ggtgg' in (s[13:18],s[14:19],s[15:20]) or 'ggggg' in (s[13:18],s[14:19],s[15:20]) or 'ggcgg' in (s[13:18],s[14:19],s[15:20]):
		score = 2
	elif 'agg' in (s[3:6],s[4:7]) or 'gag' in (s[3:6],s[4:7]) or 'gga' in (s[3:6],s[4:7]):
		score = 1
	return score

def ave(a):
	return Decimal(sum(a)/len(a))

def get_orfs(locus):
	dna = locus.seq().lower()
	#start_codons = ['atg', 'gtg', 'ttg']
	#stop_codons = ['taa', 'tga', 'tag']
	# This is the object that holds all the orfs
	my_orfs = Orfs(locus)
	my_orfs.seq = dna 
	my_orfs.contig_length = len(dna)

	# find nucleotide frequency, kmers, and create gc frame plot
	frequency = {'a':Decimal(0), 't':Decimal(0), 'c':Decimal(0), 'g':Decimal(0)}
	frame_plot = GCframe()
	background_rbs = [1.0] * 28
	training_rbs = [1.0] * 28

	for i, base in enumerate(dna):
		if(base not in ['a', 'c', 't', 'g']):
			if(base in ['s', 'b', 'v']):
				base = 'g'
			else:
				base = 'a'
		#nucl frequency
		frequency[base] += 1
		frequency[rev_comp(base)] += 1
		#kmers for rbs
		background_rbs[score_rbs(dna[i:i+21])] += 1
		background_rbs[score_rbs(rev_comp(dna[i:i+21]))] += 1
		#gc frame plot
		frame_plot.add_base(base)
	gc_pos_freq = frame_plot.get()

	Pa = frequency['a']/(my_orfs.contig_length*2)
	Pt = frequency['t']/(my_orfs.contig_length*2)
	Pg = frequency['g']/(my_orfs.contig_length*2)
	Pc = frequency['c']/(my_orfs.contig_length*2)
	my_orfs.pstop = (Pt*Pa*Pa + Pt*Pg*Pa + Pt*Pa*Pg)

	y = sum(background_rbs)
	background_rbs[:] = [x/y for x in background_rbs]

	# The dicts that will hold the start and stop codons
	stops = {1:0, 2:0, 3:0, -1:1, -2:2, -3:3}
	starts = {1:[], 2:[], 3:[], -1:[], -2:[], -3:[]}
	if dna[0:3] not in locus.start_codons:
		starts[1].append(1)
	if dna[1:4] not in locus.start_codons:
		starts[2].append(2)
	if dna[2:5] not in locus.start_codons:
		starts[3].append(3)

	# Reset iterator and find all the open reading frames
	states = itertools.cycle([1, 2, 3])
	for i in range(1, (len(dna)-1)):
		codon = dna[i-1:i+2]
		frame = next(states)
		if codon in locus.start_codons:
			starts[frame].append(i)
		elif rev_comp(codon) in locus.start_codons:
			starts[-frame].append(i+2)
		elif codon in locus.stop_codons:
			stop = i+2
			for start in reversed(starts[frame]):
				length = stop-start+1
				if(length >= my_orfs.min_orf_len):
					seq = dna[start-1:stop]
					rbs = dna[start-21:start]
					rbs_score = score_rbs(dna[start-21:start])
					my_orfs.add_orf(start, stop-2, length, frame, seq, rbs, rbs_score)
					training_rbs[rbs_score] += 1
	
			starts[frame] = []
			stops[frame] = stop
		elif rev_comp(codon) in locus.stop_codons:
			stop = stops[-frame]
			for start in starts[-frame]:
				length = start-stop+1
				if(length >= my_orfs.min_orf_len):
					seq = rev_comp(dna[max(0,stop-1):start])
					rbs = rev_comp(dna[start:start+21])
					rbs_score = score_rbs(rev_comp(dna[start:start+21]))
					my_orfs.add_orf(start-2, stop, length, -frame, seq, rbs, rbs_score)
					training_rbs[rbs_score] += 1
	
			starts[-frame] = []
			stops[-frame] = i
	# Add in any fragment ORFs at the end of the genome
	for frame in [1, 2, 3]:
		for start in reversed(starts[frame]):
			stop = my_orfs.contig_length-((my_orfs.contig_length-(frame-1))%3)
			length = stop-start+1
			#stop = i+3
			if(length >= my_orfs.min_orf_len):
				seq = dna[max(0,start-1):stop]
				rbs = dna[start-21:start]
				rbs_score = score_rbs(dna[start-21:start])
				my_orfs.add_orf(start, stop-2, length, frame, seq, rbs, rbs_score)
				training_rbs[rbs_score] += 1
		start = my_orfs.contig_length-((my_orfs.contig_length-(frame-1))%3)
		if rev_comp(dna[start-3:start]) not in locus.start_codons:
			starts[-frame].append(my_orfs.contig_length-((my_orfs.contig_length-(frame-1))%3))	
		for start in starts[-frame]:
			stop = stops[-frame]
			length = start-stop+1
			if(length >= my_orfs.min_orf_len):
				seq = rev_comp(dna[max(0,stop-1):start])
				rbs = rev_comp(dna[start:start+21])
				rbs_score = score_rbs(rev_comp(dna[start:start+21]))
				my_orfs.add_orf(start-2, stop, length, -frame, seq, rbs, rbs_score)
				training_rbs[rbs_score] += 1

	#-------------------------------Score ORFs based on RBS motif--------------------------------------#
	y = sum(training_rbs)
	training_rbs[:] = [x/y for x in training_rbs]
	for orf in my_orfs.iter_orfs():
		orf.weight_rbs = training_rbs[orf.rbs_score]/background_rbs[orf.rbs_score]

	#-------------------------------Score ORFs based on GC frame plot----------------------------------#

	x = []
	pos_max = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]
	pos_min = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]
	for orfs in my_orfs.iter_in():
		for orf in orfs:
			if(orf.start_codon() == 'atg'):
				start = orf.start
				stop = orf.stop
				if(start < stop):
					n = int((stop-start)/8)*3
					for base in range(start+n, stop-36, 3):
						pos_max[max_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])] += 1
						pos_min[min_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])] += 1
				elif(stop < start):
					n = int((start-stop)/8)*3
					for base in range(start-n, stop+36, -3):
						pos_max[max_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])] += 1
						pos_min[min_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])] += 1
				break
	# normalize to one
	y = max(pos_max)
	pos_max[:] = [x / y for x in pos_max]	
	y = max(pos_min)
	pos_min[:] = [x / y for x in pos_min]	

	for orf in my_orfs.iter_orfs():
		start = orf.start
		stop = orf.stop
		if(orf.frame > 0):
			for base in range(start, stop, 3):
				ind_max = max_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])
				ind_min = min_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])
				orf.hold = orf.hold * (((1-orf.pstop)**pos_max[ind_max])**pos_min[ind_min])
		else:
			for base in range(start, stop, -3):
				ind_max = max_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])
				ind_min = min_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])
				orf.hold = orf.hold * (((1-orf.pstop)**pos_max[ind_max])**pos_min[ind_min])

	for orf in my_orfs.iter_orfs():
		orf.score()
		#print orf.start, orf.stop, orf.pstop, 1/orf.hold, "sep", orf.rbs, orf.weight_rbs, orf.weight
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
					sys.stderr.write("error in breaking region"+str(n))
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
	
	add_trnas(my_orfs, G)
	
	#-------------------------------Connect the open reading frames to each other----------------------#
	for right_node in G.iternodes():
		r = right_node.position
		if(right_node.gene == 'CDS'):
			r_other = my_orfs.other_end[r]
		else:
			r_other = my_orfs.other_end['t'+str(r)]
		for left_node in G.iternodes():
			l = left_node.position
			if(left_node.gene == 'CDS'):
				l_other = my_orfs.other_end[l]
			else:
				l_other = my_orfs.other_end['t'+str(l)]
			if(0 < r-l < 500):
				if(l in my_orfs and my_orfs.other_end[l] in my_orfs[l]):
					o1 = my_orfs.get_orf(my_orfs.other_end[l], l).pstop
				elif(l in my_orfs):
					o1 = my_orfs.get_orf(l, my_orfs.other_end[l]).pstop
				else:
					o1 = pgap
				if(r in my_orfs and my_orfs.other_end[r] in my_orfs[r]):
					o2 = my_orfs.get_orf(my_orfs.other_end[r], r).pstop
				elif(r in my_orfs):
					o2 = my_orfs.get_orf(r, my_orfs.other_end[r]).pstop
				else:
					o2 = pgap
				pstop = ave([o1, o2])

				#trna
				if(left_node.gene == 'tRNA' or right_node.gene == 'tRNA'):
					if(left_node.frame*right_node.frame > 0 and left_node.type != right_node.type):
						if(left_node.frame > 0 and left_node.type == 'stop') or (left_node.frame < 0 and left_node.type == 'start'):
							if not G.has_edge(Edge(left_node, right_node, 1)):
								score = score_gap(r-l-3, 'same', pgap)
								G.add_edge(Edge(left_node, right_node, score ))	

					elif(left_node.frame*right_node.frame < 0 and left_node.type == right_node.type):
						if(left_node.frame > 0 and left_node.type == 'stop') or (left_node.frame < 0 and left_node.type == 'start'):
							if not G.has_edge(Edge(left_node, right_node, 1)):
								score = score_gap(r-l-3, 'same', pgap)
								G.add_edge(Edge(left_node, right_node, score ))	
				# same directions
				elif(left_node.frame*right_node.frame > 0):
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
							if(r_other+3 < l and r < l_other):
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
							#print(r_other, l, r, l_other)
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

def add_trnas(my_orfs, G):
	f = tempfile.NamedTemporaryFile(mode='wt')
	f.write(">temp\n")
	f.write(my_orfs.seq)
	f.seek(0)

	trnas = list()
	seen = list()
	output1 = output2 = None
	# Iterate over the trnas
	try:
		# no support for tmRNA yet
		#output2 = Popen(["aragorn", "-m", "-t", "-w", f.name], stdout=PIPE, stdin=PIPE, stderr=PIPE).stdout.read()
		output2 = Popen(["aragorn", "-t", "-w", f.name], stdout=PIPE, stdin=PIPE, stderr=PIPE).stdout.read()
		for line in output2.decode().splitlines():
			if not line.startswith('>') and not line.endswith('found'):
				column = line.split()
				pair = eval(column[2].replace('c',''))
				if 'c' not in column[2]:
					trnas.append(pair)
				else:
					trnas.append(pair[::-1])
				seen.extend(range(*pair))
	except:
		pass
	try:
		output1 = Popen(["tRNAscan-SE", "-B", "-q", "--brief", f.name], stdout=PIPE, stdin=PIPE, stderr=PIPE).stdout.read()
		for line in output1.decode().splitlines():
			column = line.split('\t')
			a,b = map(int, column[2:4])
			if a < b and not set(seen) & set(list(range(a,b))):
				trnas.append([a,b])
			elif a > b and not set(seen) & set(list(range(b,a))):
				trnas.append([a,b])
	except:
		pass
	if output1 is None and output2 is None:
		sys.stderr.write("Warning: tRNAscan or Aragorn were not found, proceding without tRNA masking.\n")
		return

	for start,stop in trnas:
		# Add in trna
		if(start < stop):
			source = Node('tRNA', 'start', 4, start)
			target = Node('tRNA', 'stop', 4, stop-2)
			my_orfs.other_end['t'+str(stop-2)] = start
			my_orfs.other_end['t'+str(start)] = stop-2
		else:
			source = Node('tRNA', 'stop', -4, stop)
			target = Node('tRNA', 'start', -4, start-2)
			my_orfs.other_end['t'+str(start-2)] = stop
			my_orfs.other_end['t'+str(stop)] = start-2
		G.add_edge(Edge(source, target, -Decimal(20)))


