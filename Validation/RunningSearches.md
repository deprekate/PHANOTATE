# Running Searches on Amazon Web Services

We start by creating a list of all the WGS metagenomes in the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) that we have detected in [partie](https://github.com/linsalrob/partie):

After cloning the repository, you can get those IDs by running this command:

```bash
grep WGS partie/SRA_Metagenome_Types.txt | cut -f 1 > wgs.ids
```
We calculate how many AWS instances we are going to run, and create a base AMI image that we can spawn. The AMI image needs the following software installed:

- [lastal](http://last.cbrc.jp/doc/lastal.html) &mdash; needed for the sequence alignment algorithm
- [SRA Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/) &mdash; needed for fastq-dump
- [fastq2fasta](https://raw.githubusercontent.com/linsalrob/EdwardsLab/master/bin/fastq2fasta.cpp) &mdash; needed to convert fastq to fasta. There are other ways to do that, but this is very lightweight standalone code. You will need to compile this with `c++ -o fastq2fasta ./fastq2fasta.cpp` and ensure that it is in your path. 
- [pigz](https://zlib.net/pigz/) &mdash; you can use gzip, zip, or whatever else you prefer, but we like pigz (and not just because of the name!)

We then create this script somewhere on the master image:

```bash
#!/bin/bash
# Run fastqdump and lastal on some data
# and copy the output to an S3 bucket

# num here is how many sequences we want this instance to process.
# It is determined by the number of lines in wgs.ids from partie
# and the number of instances you are running, which is also a factor
# of how impatient you are
NUM=4850
# find our instance ID. When we start in batch and instantiate a few
# of these, this is an automatically assigned integer.
AMI_ID=$(curl http://169.254.169.254/latest/meta-data/ami-launch-index)
AMI_ID=$(($AMI_ID+1))
END=$(($NUM*$AMI_ID))

export IFS=$'\n'
echo "STARTING to get $NUM FROM $END ";

# note that rather than getting from $NUM to $END, we use head and tail, so we end at $END and get the previous $NUM entries
for SRR in $(head -n $END /usr/local/sra/sra_metagenomes.txt | tail -n $NUM);
do
	echo "SRR IS $SRR"
	# download 100,000 sequences, starting at 1,000 for this entry
	# you may need to specify the full path to fastq-dump here
	# See https://edwards.sdsu.edu/research/sra for more information about this command
	fastq-dump --outdir fastq --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip -N 1000 -X 101000 $SRR
	mkdir -p fasta lastal
	# we do this for the left and right reads, although we may only use one of those later...
	for FQ in $(ls fastq/); do 
		# FA is the fasta output file. Convert fastq to fasta
		FA=$(echo $FQ | sed -e 's/fastq/fasta/');
		fastq2fasta fastq/$FQ fasta/$FA;
		
		# LA is the lastal output file. Run the sequence similarity search
		LA=$(echo $FQ | sed -e 's/fastq/lastal/');
		lastal -F15 /usr/local/sra/all_orfs.ldb fasta/$FA -f BlastTab+ -P 6 > lastal/$LA;
		
		# compress the output and copy it to our s3 bucket. My bucket is called lastaloutput
		pigz lastal/$LA;
		aws s3 cp lastal/$LA.gz s3://lastaloutput/lastal/;
	done;
	# clean up these directories so we don't store the data. That's a good way to crash the machine :)
	rm -rf fasta/ fastq/ lastal/ ncbi/;
done
```

Once we have a master image with this bash script, we start 18-20 instances of the image. Each one gets a unique ID that it retrieves using this line: `AMI_ID=$(curl http://169.254.169.254/latest/meta-data/ami-launch-index)`. Rather than make this script start on boot, we prefer to boot the instances and then make a simple file called `ips.txt` with the IP addresses. We then run the bash script as a screen daemon:

```bash
for I in $(cat ips.txt); do echo $I; ssh $I "screen -d -m /usr/local/genome/run_lastalsearch.sh"; done
```

This way, you can always log in to a couple of instances and make sure everything is running is you believe it should.

The computes will run and the compressed results will quickly start appearing in your s3 bucket. Its then up to you to download and process them!




