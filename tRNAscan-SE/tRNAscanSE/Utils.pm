# tRNAscanSE/Utils.pm
# This module contains utility functions used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2011 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::Utils;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(check_output_file open_for_read open_for_write open_for_append tempname
                 print_filename rev_comp_seq max min seg_overlap error_exit_status trim);

our %comp_map = (                     # Complement map
                'A' => 'T', 'T' => 'A', 'U' => 'A',
                'G' => 'C', 'C' => 'G',
                'Y' => 'R', 'R' => 'Y', 
                'S' => 'W', 'W' => 'S', 
                'M' => 'K', 'K' => 'M', 
                'B' => 'V', 'V' => 'B', 
                'H' => 'D', 'D' => 'H', 
                'N' => 'N', 'X' => 'X',
                '?' => '?');

sub check_output_file {
    my ($fname, $prompt_for_overwrite) = @_;
    my ($ans, $ansline);
   
    if ((-e $fname) && ($prompt_for_overwrite)) {
        print STDERR "\nWARNING: $fname exists already.\n\n",
            " (O)verwrite file, (A)ppend to file, or (Q)uit program? ";
        $ansline = <STDIN>;
        $ans = substr($ansline, 0, 1);
        while ($ans !~ /[AOQaoq]/) {
            print STDERR "\nReply (O)verwrite (A)ppend, or (Q)uit [O/A/Q]: ";
            $ansline = <STDIN>;
            $ans = substr($ansline, 0, 1);
        }
        if (uc($ans) eq 'Q') {
            die "\ntRNAscan-SE aborted.\n\n";
        }
        elsif  (uc($ans) eq 'A') {
            print STDERR "\n Appending to $fname...\n";
            open(FHAND,">>$fname") || 
                die "Unable to open $fname for appending. ",
                "Aborting program.\n";
            close(FHAND);
            return;                    # successful exit status
        }    
        else {               #  $ans eq 'O'verwrote
            print STDERR "\n Overwriting $fname...\n";
        }    
    }    
    open(FHAND, ">$fname") || 
        die "Unable to open $fname for writing.  Aborting program.\n";
    close(FHAND);
}

sub open_for_read {
    my ($FHAND, $fname) = @_;

    open($$FHAND, "$fname") || 
        die "Unable to open $fname for reading.  Aborting program.\n";
}

sub open_for_write {
    my ($FHAND, $fname) = @_;

    open($$FHAND, ">$fname") || 
        die "Unable to open $fname for writing.  Aborting program.\n";
}

sub open_for_append {
    my ($FHAND, $fname) = @_;
    
    open ($$FHAND, ">>$fname") ||
        die "FATAL:  Unable to open output file ",
        &print_filename($fname), "\n\n";
}

# Function: tempname
# by SE, modification by TMJL
# Returns a unique temporary filename. 
#
# Normally puts temp files to /tmp. This directory can
# be overridden by an environment variable TMPDIR.
#

sub tempname {
    my ($temp_dir, $exten) = @_;
    my ($name);        
    
    $name = "$temp_dir/tscan$$"."$exten";
    return $name;                               
}

sub print_filename {
    my ($fname) = @_;
    if ($fname eq "-") {
        $fname = "Standard output";
    }
    return $fname;
}

sub rev_comp_seq {
    my ($seq) = @_;
    my ($seqlen) = length($seq);
    my ($i, $j, $rcseq);

    $rcseq = 'X' x $seqlen;        # pre-extending string for efficiency
    for ($i = ($seqlen - 1), $j = 0; $i > -1; $i--, $j++) {
        substr($rcseq, $j, 1) = $comp_map{(substr($seq, $i, 1))};
    }
    return $rcseq;
}

sub min {
    my ($a, $b) = @_;
    if ($a < $b) {
        return ($a); }
    else {
        return ($b);
    }
}

sub max {
    my ($a, $b) = @_;
    if ($a > $b) {
        return ($a);
    }
    else {
        return ($b);
    }
}

sub seg_overlap {
    my ($seg1_a, $seg1_b, $seg2_a, $seg2_b) = @_;

    if ((($seg1_a >= $seg2_a) && ($seg1_a <= $seg2_b)) ||
        (($seg1_b >= $seg2_a) && ($seg1_b <= $seg2_b)) ||
        (($seg2_a >= $seg1_a) && ($seg2_a <= $seg1_b)) ||
        (($seg2_b >= $seg1_a) && ($seg2_b <= $seg1_b)))  {
        return 1;
    }
    else {
        return 0;
    }
}

sub error_exit_status {
    my ($prog_name, $seq_name) = @_;

    if ($? != 0) {
        print STDERR "$prog_name could not complete successfully for $seq_name.\n",
            "Possible memory allocation problem or missing file. (Exit code=",$?,").\n\n";
        return 1;
    }
    else {
        return 0;
    }
}

sub trim
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

1;
