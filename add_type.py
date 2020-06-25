import fileinput

stop_locations = dict()

for line_stdin in fileinput.input():

	# get genbank annotations
	if line_stdin.startswith("#id:"):
		with open(line_stdin.split()[1] + "/" + line_stdin.split()[1] + ".gb") as f:
			for line in f:
				if line.startswith("     CDS   "):
					stop_locations[ line.replace("(", " ").replace(")", " ").replace("..", " ").split()[2] ] = "coding"

	# continue with rest of PHANOTATE calls and add genbank annotation type
	elif not line_stdin.startswith("#"):
		stop = line_stdin.split()[1]
		print(line_stdin.rstrip(), stop_locations.get(stop, ''), sep='\t')
