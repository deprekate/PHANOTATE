# tRNAscanSE/Tscan.pm
# This class contains parameters and functions for running tRNAscan used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2011 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::Tscan;

use strict;
use tRNAscanSE::Utils;


sub new {
    my $class = shift;
    my $self = {};

    initialize($self);

    bless ($self, $class);
    return $self;
}

sub DESTROY
{
    my $self = shift;
}

sub initialize
{
    my $self = shift;
    
    # set to non-zero if you do NOT want redundant, overlapping hits
    #  found by tRNAscan merged into one hit
    $self->{keep_tscan_repeats} = 0;
    
    $self->{tscan_params} = "-s";   # parameter set to be used for tRNAscan
                                    # default is "-s" strict params
                                    # default for prokaryotes should be relaxed
                                    # params "-r"
    
    $self->{tscan_version} = 1.4;   # version of tRNAscan used by tRNAscan-SE

    $self->{tscan_bin} = "trnascan-1.4";
    
    $self->{tscan_mask} = 1;        # Bit-wise masks for source of tRNA hits
}

sub keep_tscan_repeats
{
    my $self = shift;
    if (@_) { $self->{keep_tscan_repeats} = shift; }
    return $self->{keep_tscan_repeats};
}

sub tscan_params
{
    my $self = shift;
    if (@_) { $self->{tscan_params} = shift; }
    return $self->{tscan_params};
}

sub tscan_version
{
    my $self = shift;
    if (@_) { $self->{tscan_version} = shift; }
    return $self->{tscan_version};
}

sub tscan_bin
{
    my $self = shift;
    if (@_) { $self->{tscan_bin} = shift; }
    return $self->{tscan_bin};
}

sub tscan_mask
{
    my $self = shift;
    return $self->{tscan_mask};
}

sub set_bin {
    
    my $self = shift;
    my $bindir = shift;
    
    # choose correct name for version being run
    # only version 1.4 is provided with distribution

    if ($self->{tscan_version} == 1.4) {
        $self->{tscan_bin} = "trnascan-1.4";
    }
    elsif ($self->{tscan_version} == 1.39) {
        $self->{tscan_bin} = "trnascan-1.39";
    }
    elsif ($self->{tscan_version} == 2) {
        $self->{tscan_bin} = "TRNAscan";
    }
    elsif ($self->{tscan_version} == 1.3) {             
        $self->{tscan_bin} = "trnascan-1.3";
    }
    else {
        die "FATAL:  Illegal tRNAscan version.\n\n";
    }

    if ($^O =~ /^MSWin/) {
        $self->{tscan_bin} .= ".exe";
    }

    if (!(-x $self->{tscan_bin})) {
        $self->{tscan_bin} = $bindir.$self->{tscan_bin};
        if (!(-x $self->{tscan_bin})) {
            die "FATAL: Unable to find ".$self->{tscan_bin}." executable\n\n";
        }
    }
}

sub run_tRNAscan {
    
    my $self = shift;
    my ($tmp_fa, $tmp_raw, $start_index, $lib_dir, $seq_name) = @_;
    my $tscan_version = $self->{tscan_version};
    my $tscan_bin = $self->{tscan_bin};
    my $tscan_params = $self->{tscan_params};

    # version provided with distribution

    if ($tscan_version == 1.4) {
        # run default tRNAscan 1.4 using selected param set
        system ("$tscan_bin -i $start_index -c $tscan_params $tmp_fa > $tmp_raw");
        if (&error_exit_status("tRNAscan", $seq_name)) {
            return -1;
        }
    }
    
    # run tRNAscan without conservative ambiguous base pairing rules
    # not available in distribution version

    elsif ($tscan_version == 1.39) {
        system ("$tscan_bin -c $tscan_params $tmp_fa > $tmp_raw"); 
    }

    # run tRNAscan v2.0, not available in distribution version

    elsif ($tscan_version == 2) {
        system ("$tscan_bin -SEQ $tmp_fa -TEMPLATE SEtemplate -OUTPUT $tmp_raw > /dev/null");
        }

    # run original tRNAscan 1.3, not available in distribution version

    elsif ($tscan_version == 1.3) {             
        if (!(-r "./TPCsignal")) {
            system ("ln -s ".$lib_dir."TPCsignal TPCsignal");
        }
        if (!(-r "./Dsignal")) {
            system ("ln -s ".$lib_dir."Dsignal Dsignal");
        }
        system ("reformat -ld genbank $tmp_fa > tmp.gb");
        system ("$tscan_bin tmp.gb $tmp_raw > /dev/null");
        system ("rm tmp.gb");
    }
    else {
        die "FATAL:  Illegal tRNAscan version.\n\n";
    }
}


# Append tRNAscan verbose output to 
#   result file with header tag

sub append_verbfile {
    
    my $self = shift;
    my ($verb_file, $tmp_fa, $seq_name) = @_;

    &open_for_append(\*TSCANVERB, $verb_file);    
    print TSCANVERB "\n>>>> tRNA-Scan verbose output for <$seq_name>\n\n";
    close TSCANVERB;
    system ("cat tscan.verb.out >> $verb_file");
}

# extract trna hits from raw result file while weeding out repeated hits
# save non-redundant hits in "hit_list" array

sub process_tRNAscan_hits {
    
    my $self = shift;
    my $constants = shift;
    my $gc = shift;
    my $stats = shift;
    my $seq_name = shift;
    my $r_hit_list = shift;
    my $tmp_raw = $constants->tmp_raw();
    
    my ($istart, $iend, $from, $to, $intron, $trnact, $len, $score,
          $anticodon, $iso_type, $sense_strand, $pos, $i);

    $trnact = 0;               # trna count for this sequence
    $istart = 0; $iend = 0;    # intron bounds
    $from = 0; $to = 0;        # tRNA bounds
    $len = 0;                  # tRNA length
    $intron = 0;               # intron present? flag
    $anticodon = '';
    $iso_type = '';        
    $score = 0;
    
    # open trnascan raw output file for current seq
    
    open (TSCANRAW, "$tmp_raw")  ||
        die ("FATAL: Unable to open temp raw output file $tmp_raw\n\n");
    
    # parse one complete hit per call 
    while ($self->parse_tscan_hit($constants, $gc, \*TSCANRAW, \$from, \$to, \$sense_strand,
                                  \$istart, \$iend, \$intron, \$len, \$iso_type,
                                  \$anticodon, \$pos))  {
        
        if ($self->{keep_tscan_repeats} ||
            (!$self->merge_repeat_hit($stats, $r_hit_list, \$trnact, $from, $to,
                                      $sense_strand, $iso_type,$score)))

            # if NOT a repeat hit, put it on the hit list 
        {            
            # check to see if tscan 1.3 has incorrectly reported
            #  start/end index (happens occassionally) 
            
            if ((abs($to - $from) + 1) != $len) {
                if ($sense_strand) {
                    $to = $from + $len - 1; }
                else {
                    $to = $from - $len + 1; }
            }
            
            $i=0;
            while (($i < scalar(@$r_hit_list)) && ($r_hit_list->[$i]{position} < $pos)) {
                $i++;
            }
            
            # save non-redundant hit 
            splice(@$r_hit_list, $i ,0, {
                seqname => $seq_name, 
                start => $from, end => $to,
                type => $iso_type, acodon => $anticodon,
                istart => $istart, iend => $iend,
                sen_strand => $sense_strand,
                position => $pos, score => 0,
                source => $self->{tscan_mask},
            });   
            
            $trnact++;        
            $stats->increment_trnatotal();
            
        }         
        
    }        # while (&Parse_tscan_hit), more hits to process for cur seq    
}

sub parse_tscan_hit {

    my $self = shift;
    my $constants = shift;
    my $gc = shift;
    my ($TSCANRAW, $r_from, $r_to, $r_sense_strand,
          $r_istart, $r_iend, $r_intron, $r_len, $r_type, $r_anticodon, $r_pos) = @_;

    my ($trna_seq) = '';

    # clear intron info parsing each hit
    $$r_istart = 0;  $$r_iend = 0;  $$r_intron = 0;

    if ($self->{tscan_version} <= 1.4)  {

        while (<$TSCANRAW>) {
            if (/^start position=\s*(\d+)\s*end position=\s*(\d+)/o)
            {  
                $$r_from = $1;
                $$r_to = $2; 
                if ($$r_from < $$r_to) {
                    $$r_sense_strand = 1;
                    $$r_pos = $$r_from  }
                else {                
                    $$r_sense_strand = 0;
                    $$r_pos = $constants->REALLY_BIG_NUMBER() - $$r_from + 1;
                }
            }
                                
            elsif (/^potential tRNA sequence=\s(.+)\n/o)  {
                $trna_seq = $1;  $$r_len = length($trna_seq);
            }                        
            elsif (/^tRNA predict as a tRNA-\s*(\S+)\s*: anticodon (\S+)/o) {
                $$r_type = $1;
                $$r_anticodon = $2;
            }
            elsif (/^anticodon includes unknown bases/o) {
                $$r_type = $gc->undef_isotype();
                $$r_anticodon = $gc->undef_anticodon();
            }
            elsif (/^potential intron between positions\s*(\d+)\s*(\d+)/o) { 
                $$r_istart = $1;
                $$r_iend = $2; 
                $$r_intron = 1;
            }
            # flag for end of current tRNA hit info
            elsif (/^number of base pairing in the anticodon/o)  {
                return 1;
            } 
            elsif (/^number of predicted tRNA=(\d+)/o) {
                return 0;        # end of hits for this seq 
            }
        }
        return 0;                # reached end of raw hits file
    }                               

    else {
        die "FATAL: Illegal tRNAscan version selected.\n\n";
    }
}        

# check current hit for redundancy against all previous hits in hitlist
#
# if it IS a repeat, merge it with overlapping hit and return 1
# if it doesn't overlap with any hits, return 0

sub merge_repeat_hit  {

    my $self = shift;
    my $stats = shift;
    my ($r_hit_list, $r_trnact, $from, $to, $sense_strand, $iso_type, $score) = @_;
    my ($i);

    foreach $i (0..(scalar(@$r_hit_list) - 1)) {
        
        if ($sense_strand) {
            if (($r_hit_list->[$i]{sen_strand} == 1) &&
                (&seg_overlap($from, $to, $r_hit_list->[$i]{start},
                             $r_hit_list->[$i]{end}))) 
            {
                $r_hit_list->[$i]{start} = &min($from, $r_hit_list->[$i]{start});
                $r_hit_list->[$i]{end} = &max($to, $r_hit_list->[$i]{end});
                $r_hit_list->[$i]{source} = $r_hit_list->[$i]{source} | $self->{tscan_mask};
                $r_hit_list->[$i]{type} = $iso_type;
                $r_hit_list->[$i]{score} = $score;
    
                # check to see if extended endpoint overlaps
                #  i+1 hit's start boundary
                # if so, combine hit[i] and hit[i+1] into one
                #  hit and delete hit[i+1]
                if (($i != (scalar(@$r_hit_list) - 1)) && ($r_hit_list->[$i+1]{sen_strand})
                    && ($r_hit_list->[$i]{end} >= $r_hit_list->[$i+1]{start})) 
                {
                    $r_hit_list->[$i]{end} = &max($r_hit_list->[$i]{end}, $r_hit_list->[$i+1]{end});
                    $r_hit_list->[$i]{source} = $r_hit_list->[$i]{source} | $r_hit_list->[$i+1]{source};

                    splice(@$r_hit_list,$i+1,1);          # toss out overlapping hit 
                    $$r_trnact--;
                    $stats->decrement_trnatotal();
                }   
                return 1;                                 # exit loop immediately
            }
        }
        else         # else (antisense) strand 
        {                
            if (($r_hit_list->[$i]{sen_strand} == 0) &&
                (&seg_overlap($to, $from, $r_hit_list->[$i]{end}, $r_hit_list->[$i]{start}))) 
            {
                $r_hit_list->[$i]{start} = &max($from, $r_hit_list->[$i]{start});
                $r_hit_list->[$i]{end} = &min($to, $r_hit_list->[$i]{end});
                $r_hit_list->[$i]{source} = $r_hit_list->[$i]{source} | $self->{tscan_mask};
                $r_hit_list->[$i]{type} = $iso_type;
                $r_hit_list->[$i]{score} = $score;

                if (($i != (scalar(@$r_hit_list) - 1)) &&
                    ($r_hit_list->[$i]{end} <= $r_hit_list->[$i+1]{start}))
                {
                    $r_hit_list->[$i]{end} = &min($r_hit_list->[$i]{end}, $r_hit_list->[$i+1]{end});
                    $r_hit_list->[$i]{source} = $r_hit_list->[$i]{source} | $r_hit_list->[$i+1]{source};

                    splice(@$r_hit_list,$i+1,1);          # toss out overlapping hit 
                    $$r_trnact--;
                    $stats->decrement_trnatotal();
                }
                return 1;                                 # exit loop immediately
            }
        } # else (antisense) strand
        
    }  # for each (hit)                        

    return 0;                                             # current hit is not a repeat
}

1;