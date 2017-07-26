import os
import sys
import itertools
import tempfile
from subprocess import Popen, PIPE, STDOUT
from decimal import Decimal

from settings import weights
from settings import start_weight
from nodes import Node
from edges import Edge
from graphs import Graph
from gc_frame_plot import GCframe
from gc_frame_plot import max_idx
from gc_frame_plot import min_idx
from math import sqrt



def rev_comp(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
	return "".join([seq_dict[base] for base in reversed(seq)])

def score_orf(pstop, startcodon, length, rbs):
	score = (1-pstop)**((length-6)/3)
	score = 1/score
	if(startcodon):
		score = score * start_weight[startcodon]
	if(rbs):
		score = score * score_rbs(rbs)
	return -score

def score_overlap(length, direction, factor=1):
	if(length > 100):
		return False
	o = Decimal(weights['overlap'])
	s = Decimal(weights['switch'])

	score = Decimal(o)**length
	score = 1/score
	#score = score*factor
	if(direction == 'diff'):
		score = abs(score)*factor + (1/s)
	return score

def score_gap(length, direction):
	g = Decimal(weights['gap'])
	s = Decimal(weights['switch'])

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
	score = Decimal(0)
	bases = ['A', 'C', 'T', 'G']
	kmers = {}

	for length in range(3,7):
		for p in itertools.product('ACTG', repeat=length):
			kmers["".join(p)] = []
	
	for n in range(0,len(seq)-5):
		for x in range(3, 7):
			push(kmers, seq[n:n+x], len(seq)-n-x)
	push(kmers, seq[-3:len(seq)], 0)
	push(kmers, seq[-4:len(seq)], 0)
	push(kmers, seq[-5:len(seq)], 0)

	#These are motifs/scores taken directly from Prodigal
	if(comp(kmers['AGGAGG'],range(4,11))):
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
	score = Decimal(score+10)
	return score.log10()

def p_stop(seq):
	frequency = {'A':Decimal(0), 'T':Decimal(0), 'C':Decimal(0), 'G':Decimal(0)}
	for base in seq:
		frequency[base] += 1
		frequency[rev_comp(base)] += 1
	Pa = frequency['A']/(len(seq)*2)
	Pt = frequency['T']/(len(seq)*2)
	Pg = frequency['G']/(len(seq)*2)
	Pc = frequency['C']/(len(seq)*2)
       	return (Pt*Pa*Pa + Pt*Pg*Pa + Pt*Pa*Pg)
def ave(a):
	return Decimal(sum(a)/len(a))

def parse(dna):
	G = Graph(directed=True)

	# Check for ambiguous bases and create gc frame plot
	frequency = {'A':Decimal(0), 'T':Decimal(0), 'C':Decimal(0), 'G':Decimal(0)}
	frame_plot = GCframe()
	for base in dna:
		if(base not in ['A', 'C', 'T', 'G']):
			sys.stderr.write('Error: ambiguous nucleotide base ' + base + ' found\n')
			sys.exit()
		frame_plot.add_base(base)
		frequency[base] += 1
		frequency[rev_comp(base)] += 1
	gc_pos_freq = frame_plot.get()


	Pa = frequency['A']/(len(dna)*2)
	Pt = frequency['T']/(len(dna)*2)
	Pg = frequency['G']/(len(dna)*2)
	Pc = frequency['C']/(len(dna)*2)
       	weights['gap'] = (1-(Pt*Pa*Pa + Pt*Pg*Pa + Pt*Pa*Pg))
	
	# The dicts that will hold the start and stop codons
	stops = {1:0, 2:0, 3:0, -1:0, -2:0, -3:0}
	starts = {1:[0], 2:[0], 3:[0], -1:[0], -2:[0], -3:[0]}

	# This tells us how far back an overlap can go (you cannot go further back than the length of the ORF)
	start_to_stop = {}
	stop_to_start = {}
	# This tells us the probabity of a stop codon for each ORF
	pstops = {}

	#This is just a quick way to hold orf scores
	start_to_score = {}

	#This holds the good ORFs: the longest ones that start with ATG
	training_orfs = {}

	pos_to_codon = {}
	start_to_rbs = {}


	# Reset iterator and find all the open reading frames
	states = itertools.cycle([1, 2, 3])
	for i in range(1, (len(dna)-1)):
		codon = dna[i-1:i+2]
		frame = states.next()

		if codon in ['ATG', 'GTG', 'TTG']:
			starts[frame].append(i)
			pos_to_codon[i] = codon	
		elif rev_comp(codon) in ['ATG', 'GTG', 'TTG']:
			starts[-frame].append(i)
			pos_to_codon[i] = rev_comp(codon)	
		elif codon in ['TAA', 'TGA', 'TAG']:
			for start in reversed(starts[frame]):
				stop = i
				length = stop-start
				if(length >= Decimal(weights['min_orf_length'])):
					nlength = start - starts[frame][0]
					if(start and pos_to_codon[start]=='TTG' and length < 100):
						continue
					pstop = p_stop(dna[max(0,start-1):stop])
					startcodon = dna[start-1:start+2]
					rbs = dna[start-21:start]
					#score
					score = score_orf(pstop, startcodon, length, rbs)
					source = Node('CDS', 'start', frame, start)
					target = Node('CDS', 'stop', frame, stop)
					G.add_edge(Edge(source, target, score))
	
					pstops[start] = pstop
					pstops[stop] = pstop
					start_to_stop[start] = stop
					stop_to_start[stop] = start
					start_to_rbs[start] = score_rbs(rbs)
					if(start and pos_to_codon[start] == 'ATG'):
						training_orfs[stop] = start	
			starts[frame] = []
			stops[frame] = i
		elif rev_comp(codon) in ['TAA', 'TGA', 'TAG']:
			for start in starts[-frame]:
				stop = stops[-frame]
				length = start-stop
				if(length >= Decimal(weights['min_orf_length'])):
					nlength = starts[-frame][-1] - start
					if(pos_to_codon[start]=='TTG' and length < 100):
						continue
					pstop = p_stop(rev_comp(dna[max(0,stop-1):start]))
					startcodon = rev_comp(dna[start-1:start+2])
					rbs = rev_comp(dna[start+2:start+23])
					#score
					score = score_orf(pstop, startcodon, length, rbs)
					source = Node('CDS', 'stop', -frame, stop)
					target = Node('CDS', 'start', -frame, start)
					G.add_edge(Edge(source, target, score))
	
					pstops[start] = pstop
					pstops[stop] = pstop
					start_to_stop[start] = stop
					stop_to_start[stop] = start
					start_to_rbs[start] = score_rbs(rbs)
					if(stop and pos_to_codon[start] == 'ATG'):
						training_orfs[stop] = start	
			starts[-frame] = []
			stops[-frame] = i
	#for e in G.iteredges():
	#	print e
	#sys.exit()

	# Add in any fragment ORFs at the end of the genome
	for frame in [1, 2, 3]:
		for start in reversed(starts[frame]):
			stop = i+3
			length = stop-start
			if(length >= Decimal(weights['min_orf_length'])):
				if(start):
					pstop = p_stop(dna[start-1:])
				else:
					pstop = p_stop(dna[start:])
				startcodon = dna[start-1:start+2]
				nonlength = start-stops[frame]
				rbs = dna[start-21:start]
				#score
				score = score_orf(pstop, startcodon, length, rbs)
				source = Node('CDS', 'start', frame, start)
				target = Node('CDS', 'stop', frame, stop)
				G.add_edge(Edge(source, target, score))

				pstops[start] = pstop
				pstops[stop] = pstop
				start_to_stop[start] = stop
				stop_to_start[stop] = start
				start_to_rbs[start] = score_rbs(rbs)
		starts[-frame].append(i+3)	
		for start in starts[-frame]:
			stop = stops[-frame]
			length = start-stop
			if(length >= Decimal(weights['min_orf_length'])):
				if(stop):
					pstop = p_stop(rev_comp(dna[stop-1:start+2]))
				else:
					pstop = p_stop(rev_comp(dna[stop:start+2]))
				startcodon = rev_comp(dna[start-1:start+2])
				nonlength = 1
				rbs = rev_comp(dna[start+2:start+23])
				#score
				score = score_orf(pstop, startcodon, length, rbs)
				source = Node('CDS', 'stop', -frame, stop)
				target = Node('CDS', 'start', -frame, start)
				G.add_edge(Edge(source, target, score))

				pstops[start] = pstop
				pstops[stop] = pstop
				start_to_stop[start] = stop
				stop_to_start[stop] = start
				start_to_rbs[start] = score_rbs(rbs)

	#-------------------------------Score ORFs based on GC frame plot----------------------------------#
	pos_max = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]
	pos_min = [Decimal(1), Decimal(1), Decimal(1), Decimal(1)]
	for stop in training_orfs:
		start = training_orfs[stop]
		loc_pos_max = [Decimal(0), Decimal(0), Decimal(0), Decimal(0)]
		loc_pos_min = [Decimal(0), Decimal(0), Decimal(0), Decimal(0)]
		if(start < stop and stop-start > 90):
			n = ((stop-start)/8)*3
			if(start == 0):
				start = (stop+2)%3+1
			for base in range(start+3+n, min(stop-30, len(dna)-1), 3):
				pos_max[max_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])] += 1
				pos_min[min_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])] += 1
		elif(stop and stop < start and start-stop > 90):
			n = ((start-stop)/8)*3
			if(start >= len(dna)):
				start = len(dna)-(stop%3)-2
			for base in range(start-n, stop+30, -3):
				pos_max[max_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])] += 1
				pos_min[min_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])] += 1

	# normalize to one
	y = max(pos_max)
	pos_max[:] = [x / y for x in pos_max]	
	y = max(pos_min)
	pos_min[:] = [x / y for x in pos_min]	
	#print pos_min
	#print pos_max
	
	for edge in G.iteredges():
		maxes = []
		mins = []
		if(edge.source.frame > 0):
			start = edge.source.position
			stop = edge.target.position
			s = 1
			if(start == 0):
				start = edge.source.frame
			for base in range(start+3+30, min(stop-30,len(dna)-1), 3):
				ind_max = max_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])
				ind_min = min_idx(gc_pos_freq[base][0],gc_pos_freq[base][1],gc_pos_freq[base][2])
				maxes.append(pos_max[ind_max])
				mins.append(pos_min[ind_min])
				s = s * pos_max[ind_max]
			start = edge.source.position
		else:
			start = edge.target.position
			stop = edge.source.position
			if(start >= len(dna)):
				start = len(dna)-(stop%3)-2
			for base in range(start-30, stop+30, -3):
				ind_max = max_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])
				ind_min = min_idx(gc_pos_freq[base][2],gc_pos_freq[base][1],gc_pos_freq[base][0])
				maxes.append(pos_max[ind_max])
				mins.append(pos_min[ind_min])
			start = edge.target.position
		#print start, stop, ave(mins), ave(maxes)
		edge.weight = -1 * abs(edge.weight)**ave(mins)
		edge.weight = -1 * abs(edge.weight)**ave(maxes)
		start_to_score[start] = edge.weight

	#sys.exit()
	#-------------------------------Check for long noncoding regions that would break the path---------#
	bases = [None] * (len(dna))
	for edge in G.iteredges():
		for n in range(edge.source.position, min(edge.target.position, len(dna)-1)):
			try:
				bases[n] = n
			except:
				print n
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
									score = score_gap(r-l-3, 'same')
									G.add_edge(Edge(left_node, right_node, score))	
								elif(left_node.type == 'start' and right_node.type =='stop' and left_node.frame < 0):
									score = score_gap(r-l-3, 'same')
									G.add_edge(Edge(left_node, right_node, score ))	
							else:
								if(left_node.type == 'stop' and right_node.type =='stop' and left_node.frame > 0):
									score = score_gap(r-l-3, 'diff')
									G.add_edge(Edge(left_node, right_node, score ))	
								elif(left_node.type == 'start' and right_node.type =='start' and left_node.frame < 0):
									score = score_gap(r-l-3, 'diff')
									G.add_edge(Edge(left_node, right_node, score ))	
			last = base

	#-------------------------------Add in tRNA data---------------------------------------------------#
	#add_trnas(G)
	
	#-------------------------------Connect the open reading frames to each other----------------------#
	for right_node in G.iternodes():
		r = right_node.position
		for left_node in G.iternodes():
			l = left_node.position
			if(0 < r-l < 300):
				weights['overlap'] = (1-(pstops[l]+pstops[r])/2)
				# same directions
				if(left_node.frame*right_node.frame > 0):
					if(left_node.type == 'stop' and right_node.type =='start'):
						if(left_node.frame > 0):
							score = score_gap(r-l-3, 'same')
							G.add_edge(Edge(left_node, right_node, score ))	
						else:
							if(left_node.frame != right_node.frame):
								if(r < stop_to_start[l] and start_to_stop[r] < l):
									score = score_overlap(r-l+3, 'same')
									if(score):
										G.add_edge(Edge(right_node, left_node, score ))	
					if(left_node.type == 'start' and right_node.type =='stop'):
						if(left_node.frame > 0):
							if(left_node.frame != right_node.frame):
								if(r < start_to_stop[l] and stop_to_start[r] < l):
									score = score_overlap(r-l+3, 'same')
									if(score):
										G.add_edge(Edge(right_node, left_node, score ))	
						else:
				
							score = score_gap(r-l-3, 'same')
							G.add_edge(Edge(left_node, right_node, score ))	
				# different directions
				else:
					if(left_node.type == 'stop' and right_node.type =='stop'):
						if(right_node.frame > 0):
							if(stop_to_start[r] < l and r < stop_to_start[l]):
								score = score_overlap(r-l+3, 'diff')
								if(score):
									G.add_edge(Edge(right_node, left_node, score ))	
						else:
							score = score_gap(r-l-3, 'diff')
							G.add_edge(Edge(left_node, right_node, score ))	
					if(left_node.type == 'start' and right_node.type =='start'):
						if(right_node.frame > 0 and r-l > 2):
							score = score_gap(r-l-3, 'diff')
							G.add_edge(Edge(left_node, right_node, score ))	
						elif(right_node.frame < 0):
							if(start_to_stop[r] < l and r < start_to_stop[l]):
								factor = ave([start_weight[pos_to_codon[l]]*start_to_rbs[l], start_weight[pos_to_codon[r]]*start_to_rbs[r]])+2
								score = score_overlap(r-l+3, 'diff', factor)
								if(score):
									G.add_edge(Edge(right_node, left_node, score ))	
	#-------------------------------Connect open reading frames at both ends to a start and stop-------#
	source = Node('source', 'source', 0, 0)
	target = Node('target', 'target', 0, len(dna)+1)
	G.add_node(source)
	G.add_node(target)
	for node in G.iternodes():
		if(node.position <= 2000):
			if( (node.type == 'start' and node.frame > 0) or (node.type =='stop' and node.frame < 0) ):
				score = score_gap(node.position, 'same')
				G.add_edge(Edge(source, node, score))
		if(len(dna) - node.position <= 2000):
			if( (node.type == 'start' and node.frame < 0) or (node.type =='stop' and node.frame > 0) ):
				score = score_gap(len(dna)-node.position, 'same')
				G.add_edge(Edge(node, target, score))

	return G
#---------------------------------------END OF LOOP----------------------------------------------------------#

def add_trnas(G):
	trna_end_left = {}
	trna_end_right = {}
	try:
		output = Popen(["../tRNAscan-SE/tRNAscan-SE", "-B", "-q", "-b", infile], stdout=PIPE, stdin=PIPE, stderr=PIPE).stdout.read()
	except:
		sys.stdout.write("Warning: tRNAscan not found, proceding without tRNA masking.\n")
		return
	#Connect the tRNAs to ORFs and the ends
	for line in output.splitlines():
		column = line.split('\t')
		if(int(column[2]) < int(column[3])):
			l = int(column[2])
			r = int(column[3])
			trna_end_left[l] = 'start'
			trna_end_right[r] = 'stop'
		
		else:
			l = int(column[3])
			r = int(column[2])
			trna_end_left[l] = 'stop'
			trna_end_right[r] = 'start'
		G.add_edge(Edge('t'+str(l), 't'+str(r), -Decimal(10**2)))	
		for i in range(max(0, l-400), l):
			codon = codon_at_position[i]
			if not codon: continue
			score = Decimal(g**(l-i))
			score = 1/score
			if(codon.type == 'stop' and codon.frame > 0) or (codon.type == 'start' and codon.frame < 0):
				G.add_edge(Edge(i, 't'+str(l), score))
		for i in range( r, min(len(dna), r+400)):
			codon = codon_at_position[i]
			if not codon: continue
			score = Decimal(g**(i-r))
			score = 1/score
			if(codon.type == 'start' and codon.frame > 0) or (codon.type == 'stop' and codon.frame < 0):
				G.add_edge(Edge( 't'+str(r), i, score))
		if(l < 400):
			G.add_edge(Edge( 0, 't'+l, l))
		if(len(dna)-r < 400):
			G.add_edge(Edge('t'+str(r), len(dna)+1, r))
	#Connect the tRNAs to each other
	for t1 in trna_end_right:
		for t2 in trna_end_left:
			distance = t2-t1
			score = score_gap(distance, 'same')
			score = 1/score
			if( 0 < distance < 400):
				G.add_edge(Edge('t'+str(t2), 't'+str(t1), score))
