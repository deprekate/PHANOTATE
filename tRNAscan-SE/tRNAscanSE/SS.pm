# tRNAscanSE/SS.pm
# This class contains parameters and functions for secondary structure parsing used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2011  Patricia P. Chan & Todd M. Lowe
# --------------------------------------------------------------
#

package tRNAscanSE::SS;

use strict;
use tRNAscanSE::Utils;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(valid_structure get_acceptor_half);

sub valid_structure {
	
	my ($ss, $canonical_intron_len) = @_;
	my $stem_index = 0;
	
	my %valid = ();
	$valid{tRNA} = 1;
	$valid{acceptor} = 1;
	$valid{darm} = 1;
	$valid{anticodon} = 1;
	$valid{variable} = 1;
	$valid{tstem} = 1;
	
	my $total_mismatches = 0;
	
	my ($r_stems, $r_mismatches) = &get_stems($ss);

	for ($stem_index = 0; $stem_index < scalar(@$r_mismatches); $stem_index++) {
		$total_mismatches += $r_mismatches->[$stem_index];
	}
	if (($total_mismatches > 1) || (scalar(@$r_stems) < 4) || (scalar(@$r_stems) > 5)) {
		$valid{tRNA} = 0;
	}

	if (scalar(@$r_stems) == 5) {
		$valid{acceptor} = 0 if (($r_mismatches->[0] > 1) || (&get_stem_length($r_stems->[0]) != 7));
		$valid{darm} = 0 if (($r_mismatches->[1] > 1) || (&get_stem_length($r_stems->[1]) < 3));
		$valid{anticodon} = 0 if (($r_mismatches->[2] > 1) || (&get_stem_length($r_stems->[2]) != 5));
		$valid{variable} = 0 if (($r_mismatches->[3] > 1) || (&get_stem_length($r_stems->[3]) < 2));
		$valid{tstem} = 0 if (($r_mismatches->[4] > 1) || (&get_stem_length($r_stems->[4]) != 5));
		$valid{tRNA} = $valid{acceptor} && $valid{darm} && $valid{anticodon} && $valid{variable} &&
			$valid{tstem} && (length($ss) - $canonical_intron_len <= 90);
	}
	elsif (scalar(@$r_stems) == 4) {
		$valid{variable} = 0;
		$valid{acceptor} = 0 if (($r_mismatches->[0] > 1) || (&get_stem_length($r_stems->[0]) != 7));
		$valid{darm} = 0 if (($r_mismatches->[1] > 1) || (&get_stem_length($r_stems->[1]) < 3));
		$valid{anticodon} = 0 if (($r_mismatches->[2] > 1) || (&get_stem_length($r_stems->[2]) != 5));
		$valid{tstem} = 0 if (($r_mismatches->[3] > 1) || (&get_stem_length($r_stems->[3]) != 5));
		$valid{tRNA} = $valid{acceptor} && $valid{darm} && $valid{anticodon} && $valid{tstem} && (length($ss) - $canonical_intron_len <= 80);
	}
	elsif (scalar(@$r_stems) == 3) {
		$valid{variable} = 0;
		$valid{acceptor} = 0 if (($r_mismatches->[0] > 1) || (&get_stem_length($r_stems->[0]) != 7));
		if ($r_mismatches->[1] == 0) {
			$valid{darm} = 0 if (($r_mismatches->[1] > 1) || (&get_stem_length($r_stems->[1]) < 3));
			$valid{anticodon} = 0 if (($r_mismatches->[2] > 1) || (&get_stem_length($r_stems->[2]) != 5));
			$valid{tstem} = 0;
		}
		elsif ($r_mismatches->[2] == 0) {
			$valid{darm} = 0;
			$valid{anticodon} = 0 if (($r_mismatches->[1] > 1) || (&get_stem_length($r_stems->[1]) != 5));
			$valid{tstem} = 0 if (($r_mismatches->[2] > 1) || (&get_stem_length($r_stems->[2]) != 5));
		}
		else {
			$valid{acceptor} = 0;
			$valid{darm} = 0;
			$valid{anticodon} = 0;
			$valid{variable} = 0;
			$valid{tstem} = 0;		
		}
	}
	else {
		$valid{acceptor} = 0;
		$valid{darm} = 0;
		$valid{anticodon} = 0;
		$valid{variable} = 0;
		$valid{tstem} = 0;		
	}
	
	return \%valid;
}

sub get_stem_length {
	my ($r_stem) = @_;
	
	return ($r_stem->{end_left} - $r_stem->{start_left} + 1);
}

sub get_stems {
	
	my ($ss) = @_;
	my %pairs = ();
	my @left = ();
	my @right = ();
	my @stems = ();
	my @mismatches = ();
	my $left_index = -1;
	my $right_index = -1;
	
	my $last_right_index = -1;
	my $start_left_index = -1;
	my $end_left_index = -1;
	my $start_right_index = -1;
	my $end_right_index = -1;
	
	for (my $pos = 0; $pos < length($ss); $pos++) {
		if (substr($ss, $pos, 1) eq ">") {
			push(@left, $pos);
			$pairs{$#left} = -1;
		}
		elsif (substr($ss, $pos, 1) eq "<") {
			push(@right, $pos);
			$left_index = scalar(@left) - 1;
			while (($pairs{$left_index} > -1) && ($left_index > -1)) {
				$left_index--;
			}
			if (($left_index > -1) && ($pairs{$left_index} == -1)) {
				$pairs{$left_index} = scalar(@right) - 1;
			}
		}
	}

	foreach $left_index (sort { $a <=> $b } keys %pairs) {
		if ($last_right_index == -1) {
			$start_left_index = $left_index;
			$end_right_index = $pairs{$left_index};						
		}
		elsif ($pairs{$left_index} != ($last_right_index - 1)) {
			$end_left_index = $left_index - 1;
			$start_right_index = $last_right_index;
			push(@stems, {start_left=>$left[$start_left_index], end_left=>$left[$end_left_index],
						  start_right=>$right[$start_right_index], end_right=>$right[$end_right_index]});
			$start_left_index = $left_index;
			$end_right_index = $pairs{$left_index};			
		}
		$last_right_index = $pairs{$left_index};
	}
	if ($last_right_index > -1) {
		$end_left_index = $left_index - 1;
		$start_right_index = $last_right_index;
		push(@stems, {start_left=>$left[$start_left_index], end_left=>$left[$end_left_index],
					  start_right=>$right[$start_right_index], end_right=>$right[$end_right_index]});
	}
		
    # find mismatches in stems
    for (my $ct = 0; $ct < scalar(@stems); $ct++) {
		$mismatches[$ct] = 0;
		$left_index = $stems[$ct]->{end_left};
		$right_index = $stems[$ct]->{start_right};
		while ($left_index >= $stems[$ct]->{start_left} && $right_index <= $stems[$ct]->{end_right}) {
			if (substr($ss, $left_index, 1) eq ".") {
				if (substr($ss, $right_index, 1) eq ".") {
					$mismatches[$ct] += 1;
					$right_index++;
				}
				$left_index--;
			}
			elsif (substr($ss, $right_index, 1) eq ".") {
				$right_index++;
			}
			else {
				$left_index--;
				$right_index++;
			}
		}
	}
	
	return (\@stems, \@mismatches);
}

sub get_acceptor_half {
	my ($seq, $ss, $half) = @_;
	
	if ($half eq "5h") {
		return substr($seq, 0, 7);
	}
	elsif ($half eq "3h") {
		if (substr($ss, length($ss) - 12) eq "<...........") {
			return substr($seq, length($seq) - 11, 7);
		}
		else {
			return substr($seq, length($seq) - 8, 7);
		}
	}
	else {
		return "";
	}
}