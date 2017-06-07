#!/usr/bin/env python
import sys
import getopt

from subprocess import Popen, PIPE, STDOUT

import os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/lib')
import make_graph
import file_handling
from nodes import Node
from settings import weights

# Test if FastPath was compliled
#if(not os.path.dirname(os.path.realpath(__file__))+'/fastpath/fastpathz'):
#		sys.stderr.write("Error: fastpathz binary not found, did you type 'make'?\n")
#		sys.exit()
#--------------------------------------------------------------------------------------------------#
#                               ARGUEMENTS                                                         #
#--------------------------------------------------------------------------------------------------#

args = file_handling.get_args()

if(args.gap_penalty):
	weights['gap'] = str(args.gap_penalty)
if(args.overlap_penalty):
	weights['overlap'] = str(args.overlap_penalty)
if(args.switch_penalty):
	weights['switch'] = str(args.switch_penalty)

#--------------------------------------------------------------------------------------------------#
#                               FILE INPUT                                                         #
#--------------------------------------------------------------------------------------------------#

my_contigs = file_handling.read_fasta(args.infile);

#--------------------------------------------------------------------------------------------------#
#                               MAIN ROUTINE                                                       #
#--------------------------------------------------------------------------------------------------#
for id, seq in my_contigs.items():

	#-------------------------------Create the Graph-------------------------------------------#
	my_graph = make_graph.parse(seq)
	#-------------------------------Run Bellman-Ford-------------------------------------------#
	source = "Node('source','source',0,0)"
	target = "Node('target','target',0," + str(len(seq)+1) + ")"
	try:
		proc = Popen([os.path.dirname(os.path.realpath(__file__))+'/fastpath/fastpathz', '-s', source, '-t', target], stdout=PIPE, stdin=PIPE, stderr=PIPE)
	except:
		sys.stdout.write("Error running fastpathz. Did you run make to compile the binary?\n")
		sys.exit()
	for e in my_graph.iteredges():
		proc.stdin.write(repr(e.source) + "\t" + repr(e.target) + "\t" + str(e.weight*10) + "\n")
		#sys.stdout.write(repr(e.source) + "\t" + repr(e.target) + "\t" + str(e.weight) + "\n")
	output = proc.communicate()[0].rstrip()

	my_path = output.split('\n')
	#for n in my_path:
	#	print n
	#sys.exit()
	#Determine whether the first edge is a gap or an fragment open reading frame
	try:
		node1 = eval(my_path[1])
		node2 = eval(my_path[2])
		if(node1.type == 'stop' and node1.frame < 0):
			if(node1.frame == node2.frame):
				my_path = my_path[1:]
			else:
				my_path = my_path[2:]
		elif(node1.type == 'start' and node1.frame > 0):
			if(node1.frame == node2.frame):
				my_path = my_path[1:]
			else:
				my_path = my_path[2:]
	except:
		sys.stdout.write("Error running fastpathz: " + output + '\n')
 

	#-------------------------------Write Output ----------------------------------------------#
	file_handling.write_output(id, args, my_path)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#






