import re
import sys
import itertools
from ast import literal_eval as make_tuple

def pairwise(iterable):
    a = iter(iterable)
    return zip(a, a)

def grouper(input_list, n = 2):
    for i in range(len(input_list) - (n - 1)):
        yield input_list[i:i+n]


if len(sys.argv) < 3:
	print("usage: get_distance.py 191,734,720,1900 edges.txt")
	exit()

edges = dict()

with open(sys.argv[2]) as f:
	for line in f:
		line = re.sub('_\w*', '', line)
		tup = make_tuple(line)
		edges[tup[:2]] = tup[2]

distance = 0
path = [n if not i%2 else str(int(n)-2) for i,n in enumerate(sys.argv[1].split(","))]

for left, right in grouper(path):
	
	d = edges[(left, right)]
	print(left, right, d)
	distance += int(d)
print("----------------")
print(distance)

