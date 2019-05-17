#!/usr/bin/env python
import os
import sys
import getopt

from subprocess import Popen, PIPE, STDOUT

sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/lib')
import file_handling
import functions
from nodes import Node

# Test if FastPath was compliled
if(not os.path.dirname(os.path.realpath(__file__))+'/fastpath/fastpathz'):
		sys.stderr.write("Error: fastpathz binary not found, did you type 'make'?\n")
		sys.exit()
#--------------------------------------------------------------------------------------------------#
#                               ARGUMENTS                                                          #
#--------------------------------------------------------------------------------------------------#

args = file_handling.get_args()

#--------------------------------------------------------------------------------------------------#
#                               FILE INPUT                                                         #
#--------------------------------------------------------------------------------------------------#

my_contigs = file_handling.read_fasta(args.infile);
if not my_contigs:
	sys.stdout.write("Error: no sequences found in infile\n")
	sys.exit()

#--------------------------------------------------------------------------------------------------#
#                               MAIN ROUTINE                                                       #
#--------------------------------------------------------------------------------------------------#
for id, seq in my_contigs.items():

	#-------------------------------Create the Graph-------------------------------------------#
	my_orfs = functions.get_orfs(seq)

	my_graph = functions.get_graph(my_orfs)
	#-------------------------------Run Bellman-Ford-------------------------------------------#
	source = "Node('source','source',0,0)"
	target = "Node('target','target',0," + str(len(seq)+1) + ")"
	try:
		proc = Popen([os.path.dirname(os.path.realpath(__file__))+'/fastpath/fastpathz', '-s', source, '-t', target], stdout=PIPE, stdin=PIPE, stderr=PIPE)
	except:
		sys.stdout.write("Error running fastpathz. Did you run make to compile the binary?\n")
		sys.exit()
	# Write edges to the fastpath program, and multiply the weight to not lose decimal places
	if args.dump:
		[sys.stdout.write(repr(e.source) + "\t" + repr(e.target) + "\t" + str(e.weight*100000) + "\n") for e in my_graph.iteredges()]
		 sys.exit()
	for e in my_graph.iteredges():
		my_message = (repr(e.source) + "\t" + repr(e.target) + "\t" + str(e.weight*100000) + "\n")
		proc.stdin.write(my_message.encode('utf-8'))
	raw_output = proc.communicate()[0].rstrip()
	output = raw_output.decode('utf-8')
	
	my_path = output.split('\n')
	
	#-------------------------------Write Output ----------------------------------------------#
	file_handling.write_output(id, args, my_path, my_graph, my_orfs)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#

