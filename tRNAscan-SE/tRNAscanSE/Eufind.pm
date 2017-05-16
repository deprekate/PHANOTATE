# tRNAscanSE/Eufind.pm
# This class contains parameters and functions for running eufindtRNA used in tRNAscan-SE.
#
# --------------------------------------------------------------
# This module is part of the tRNAscan-SE program.
# Copyright (C) 2011 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::Eufind;

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
    $self->{eufind_params} = "-r";      # relaxed params to be used with 
                                        # eufindtRNA program by default
                                        # this option selects tRNAs,  
                                        # not looking for poly T 
                                        # pol III termination signal

    $self->{eufind_intscore} = -32.10;  # Intermediate score cutoff for use
                                        # with eufindtRNA
#    $self->{eufind_Totscore} = -31.8;  # Total score cutoff for use
                                        # with eufindtRNA in non-relaxed mode
        
    $self->{eufind_bin} = "eufindtRNA";

    $self->{eufind_mask} = 2;           # Bit-wise masks for source of tRNA hits
}

sub eufind_params
{
    my $self = shift;
    if (@_) { $self->{eufind_params} = shift; }
    return $self->{eufind_params};
}

sub eufind_intscore
{
    my $self = shift;
    if (@_) { $self->{eufind_intscore} = shift; }
    return $self->{eufind_intscore};
}

sub eufind_bin
{
    my $self = shift;
    if (@_) { $self->{eufind_bin} = shift; }
    return $self->{eufind_bin};
}

sub eufind_mask
{
    my $self = shift;
    return $self->{eufind_mask};
}

sub set_bin {
    
    my $self = shift;
    my $bindir = shift;
    
    if ($^O =~ /^MSWin/) {
        $self->{eufind_bin} .= ".exe";
    }

    if (!(-x $self->{eufind_bin})) {
        $self->{eufind_bin} = $bindir.$self->{eufind_bin};
        if (!(-x $self->{eufind_bin})) {
            die "FATAL: Unable to find ".$self->{eufind_bin}." executable\n\n";
        }
    }
}

sub run_eufind {
    
    my $self = shift;
    my ($tmp_fa, $start_index, $max_int_len, $seq_name) = @_;
    my $eufind_bin = $self->{eufind_bin};
    my $eufind_intscore = $self->{eufind_intscore};
    my $eufind_params = $self->{eufind_params};

    # run default Eufind using selected param set
    my $eufind_output = 
        `$eufind_bin -i $start_index -F -I $eufind_intscore -l $max_int_len $eufind_params $tmp_fa`;
    if (&error_exit_status("EufindtRNA",$seq_name)) {
        $eufind_output = "";
    }
    return $eufind_output;
}

sub process_Eufind_hits {

    my $self = shift;
    my $constants = shift;
    my $stats = shift;
    my $r_hit_list = shift;
    my $eufind_output = shift;
    
    my ($istart, $iend, $from,$ to, $intron, $trnact, $len, $seq_name,
          $anticodon, $iso_type, $sense_strand, $score, $pos, $i, @eufind_lines);

    $trnact = 0;               # trna count for this sequence
    $istart = 0; $iend = 0;    # intron bounds
    $from = 0; $to = 0;        # tRNA bounds
    $len = 0;                  # tRNA length
    $intron = 0;               # intron present? flag
    $anticodon = '';
    $iso_type = '';        
    $score = 0.0;
    
    @eufind_lines = split(/\n/, $eufind_output);
    foreach (@eufind_lines) {
        if (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/o)
        {
            $seq_name = $1;
            $trnact = $2; 
            $from = $3;
            $to = $4;
            $iso_type = $5;
            $anticodon = $6;
            $score = $9;
            
            $istart = 0;
            $iend = 0;
            
            if ($from < $to)  {
                $len = $to - $from + 1;
                $pos = $from;                
                $sense_strand = 1;     # flag for forward or reverse strand
            }
            else  { 
                $len = $from - $to + 1;;
                $pos = $constants->REALLY_BIG_NUMBER() - $from + 1;
                $sense_strand = 0;
            }
            
            if ($from == $to) {
                print STDERR "Error reading EufindtRNA results: ",
                    "tRNA of length 0"; 
            }
            
            if (!$self->merge_repeat_hit($stats, $r_hit_list, \$trnact, $from, $to,
                                         $sense_strand, $iso_type, $score)) {
            
                # insert non-redundant hit in order
                # 'Merge_repeat_hits' depends on list being in order

                $i=0;
                while (($i < scalar(@$r_hit_list)) && ($r_hit_list->[$i]{position} < $pos)) {
                    $i++;
                }
                       
                splice(@$r_hit_list, $i, 0, {
                    seqname => $seq_name, 
                    start => $from, end => $to,
                    type => $iso_type, acodon => $anticodon,
                    istart => 0, iend => 0,
                    sen_strand => $sense_strand,
                    position => $pos, score => $score,
                    source => $self->{eufind_mask}
                });   
                
                $trnact++;        
                $stats->increment_trnatotal();
                
            }
        }
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
                (&seg_overlap($from, $to, $r_hit_list->[$i]{start}, $r_hit_list->[$i]{end}))) 
            {
                $r_hit_list->[$i]{start} = &min($from, $r_hit_list->[$i]{start});
                $r_hit_list->[$i]{end} = &max($to, $r_hit_list->[$i]{end});
                $r_hit_list->[$i]{source} = $r_hit_list->[$i]{source} | $self->{eufind_mask};
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
                (&seg_overlap($to,$from,$r_hit_list->[$i]{end}, $r_hit_list->[$i]{start}))) 
            {
                $r_hit_list->[$i]{start} = &max($from, $r_hit_list->[$i]{start});
                $r_hit_list->[$i]{end} = &min($to, $r_hit_list->[$i]{end});
                $r_hit_list->[$i]{source} = $r_hit_list->[$i]{source} | $self->{eufind_mask};
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
        
    } # for each (hit)                        

    return 0;                                             # current hit is not a repeat
}
1;
