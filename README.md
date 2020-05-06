Introduction
------------

PHANOTATE is a tool to annotate phage genomes.  It uses the assumption that non-coding
bases in a phage genome is disadvantageous, and then populates a weighted graph to
find the optimal path through the six frames of the DNA where open reading frames
are beneficial paths, while gaps and overlaps are penalized paths.

To install `PHANOTATE`,
```
pip3 install phanotate
```

or

```sh
 git clone https://github.com/deprekate/PHANOTATE.git
 cd PHANOTATE
 python3 setup.py install
```


PHANOTATE Example
--------------

Run on included sample data:
```sh
./phanotate.py tests/NC_001416.1.fasta 
```
Output is the predicted ORFs, and should look like
```sh
125     187     +
191     736     +
741     2636    +
2633    2839    +
2836    4437    +
4319    5737    +
...
```

