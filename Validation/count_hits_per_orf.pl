#!/usr/bin/perl
use Data::Dumper;
#use bignum;

#object structure to hold hits to orfs
package Hits;
sub new {
	my $class = shift;
	my $line = shift;
	my @line = split(/\t/, $line);
	my $self = {
		_evalue => $line[10],
		_bitscore  => $line[1],
		_hits       => ($line),
	};
	bless $self, $class;
	return $self;
}
sub addHit {
	my ($self, $hit) = @_;
	push @{$self->{_hits}}, $hit;
}
sub evalue {
	my( $self ) = @_;
	return $self->{_evalue};
}
sub getHits {
	my( $self ) = @_;
	return @{$self->{_hits}};
}

#get a list of all your lastal output files
my $loc = "/home3/redwards/phage/THEA/testing_predictions/protein/lastal_alignments/";
opendir my $dir, $loc or die "Cannot open directory: $!";
my @files = sort readdir $dir;
closedir $dir;
shift @files;
shift @files;

#loop over each lastal file to count the hits to orfs
my %orf_count;
my $last;
my $object = new Hits();
foreach my $file (@files){
	if($file =~ m/gz/){
		open(IN, "gunzip -c $loc/$file |") || die "can’t open pipe to $file";
	}else{
		open(IN, "$loc/$file") || die "can’t open pipe to $file";
		next;
	}

	while(<IN>){
		next if(m/^#/);
		chomp();
		my @line = split(/\t/);

		next if($line[10] > 1e-10);

		if($line[0] ne $last){
			my @hits = $object->getHits();
			my %seen = {};
			foreach my $hit (@hits){
				my @line = split(/\t/, $hit);
				my $from = $line[1];
				$from =~ s/\..*//;
				#here is where you get rid of ones youve seen
				next if $seen{$from};
				#ones you have not seen yet
				$seen{$from} = 1;
				$orf_count{$line[1]} += 1;
			}
			$object = new Hits($_);
		}
		$object->addHit($_);
		$last = $line[0];
	}
	close(IN);
}
foreach my $orf (keys %orf_count){
	print $orf, "\t", $orf_count{$orf}, "\n";
}
