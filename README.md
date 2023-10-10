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
 pip3 install PHANOTATE/.
```


PHANOTATE Example
--------------

Run on included sample data:
```sh
phanotate.py tests/NC_001416.1.fasta 
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

`PHANOTATE` has the ability to output different formats: genbank, gff, gff3, fasta

Output a genbank file that contains the genes and genome:
```sh
$ phanotate.py tests/phiX174.fasta -f genbank | head 
LOCUS       phiX174                 5386 bp 
FEATURES             Location/Qualifiers
     CDS             100..627
                     /note=score:-4.827981E+02
     CDS             687..1622
                     /note=score:-4.857517E+06
     CDS             1686..3227
                     /note=score:-3.785434E+10
     CDS             3224..3484
                     /note=score:-3.779878E+02
```

Output the nucleotide bases of the gene calls in fasta format:
```sh
$ phanotate.py tests/phiX174.fasta -f fna | head -n2
>phiX174_CDS_[100..627] [note=score:-4.827981E+02]
atgtttcagacttttatttctcgccataattcaaactttttttctgataagctggttctcacttctgttactccagcttcttcggcacctgttttacagacacctaaagctacatcgtcaacgttatattttgatagtttgacggttaatgctggtaatggtggttttcttcattgcattcagatggatacatctgtcaacgccgctaatcaggttgtttctgttggtgctgatattgcttttgatgccgaccctaaattttttgcctgtttggttcgctttgagtcttcttcggttccgactaccctcccgactgcctatgatgtttatcctttgaatggtcgccatgatggtggttattataccgtcaaggactgtgtgactattgacgtccttccccgtacgccgggcaataacgtttatgttggtttcatggtttggtctaactttaccgctactaaatgccgcggattggtttcgctgaatcaggttattaaagagattatttgtctccagccacttaagtga
```

Output the amino-acids of the gene calls in fasta format:
```sh
$ phanotate.py tests/phiX174.fasta -f faa | head -n2
>phiX174_CDS_[100..627] [note=score:-4.827981E+02]
MFQTFISRHNSNFFSDKLVLTSVTPASSAPVLQTPKATSSTLYFDSLTVNAGNGGFLHCIQMDTSVNAANQVVSVGADIAFDADPKFFACLVRFESSSVPTTLPTAYDVYPLNGRHDGGYYTVKDCVTIDVLPRTPGNNVYVGFMVWSNFTATKCRGLVSLNQVIKEIICLQPLK*

