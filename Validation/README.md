# Validating Phanotate gene calls

The big question when developing a new ORF caller is simply *"how do you know it is right?"*

We have come up with several ways to validate phanotate, but none of them are perfect!

In an ideal world we would have some proteomics data where we can identify whether the ORFs we have predicted are actually transcribed and translated. There are a few phage proteomics papers (but not many), and in those cases they have already run the *m/z* measurements through software that compares the measurements to ORFs predicted by other software. We were unable to identify a source of raw data - that is data that has not been compared to a protein file. Obviously that is needed to truly validate an ORF caller.

## Comparing protein lengths

We have included a [Jupyter notebook](Phanotate Gene Length Statistics.ipynb) that demonstrates the statistics we performed on the lengths of the genes predicted by phanotate, [Prodigal](https://github.com/hyattpd/Prodigal), [Glimmer](https://ccb.jhu.edu/software/glimmer/), and [GeneMarkS](http://exon.gatech.edu/GeneMark/). This repository includes the raw data for you to repeat the calculations.

## Comparing to the Sequence Read Archive

The [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) contains approximately 10<sup>16</sup> bp of DNA sequences from all domains of life. In [separate work](https://github.com/linsalrob/partie) we have mined this data to identify all the random community metagenomics sequences (i.e. not the 16S microbiome sequences). We hypothesized that phage sequences that truly encode proteins would be found more frequently in this data than random sequences - or than ORFs that do not encode proteins and thus are not under selection!

We used [lastal](http://last.cbrc.jp/) to compare all the protein sequences we predicted, all the protein sequences predicted by other software, and all the ORFs that are not predicted to be proteins against the random metagenomes in the Sequence Read Archive.

We ran that computation on [Amazon Web Services](https://aws.amazon.com) as [described in "RunningSearches"](RunningSearches.md).

This created 94,652 alignment files, one each with an output from lastal, and you can download a [tarball of all the data](https://edwards.sdsu.edu/data/phanotate_lastal_alignments.tgz) (*Note:* this compressed archive file is 12 GB).

We counted the sequence similarities in those files to determine how many predicted proteins were found in the different metagenomes using [count hits per orf](count_hits_per_orf.pl). The output of those [counts are availble](count_types.tsv.gz) as a tab-separated values file. That file has three columns:

1. the source gene caller &mdash; either ANY for one of the standard gene callers ([Prodigal](https://github.com/hyattpd/Prodigal), [Glimmer](https://ccb.jhu.edu/software/glimmer/), and [GeneMarkS](http://exon.gatech.edu/GeneMark/)), NONE for a CDS that was not predicted to be a gene by any of the software, or PHAN for a CDS that was *only* predicted to be a gene by phanotate.
2. The ID of the protein. This is in the format [RefSeq ID].[ORF Number] thus the ID NC\_000871.20269 is from [NC\_000871](https://www.ncbi.nlm.nih.gov/nuccore/9632893) and is ORF number 20,269 that we identified in that genome.
3. The number of times that predicted ORF was seen in the lastal results.

Those counts become the input to our [jupyter notebook](lastal_counts.ipynb) that shows how to calculate the statistics on this data.
